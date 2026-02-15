"""Functions to convert ERA5 Land data files to DATM stream files.

This module contains functionality to control the interpretation of parameters
and reading/writing of ERA5 Land source files and DATM7 three-stream output
files. The conversion functionality is imported from the `convert_data` module.

Functions
---------
convert_era5_file
    Converts a single ERA5 Land GRIB file to a DATM7 threestream netCDF file,
    optionally using the last and first time steps ofadjacent files to properly
    compute cumulative variables.
convert_monthly_era5_files
    Converts multiple ERA5 Land GRIB files to DATM7 threestream netCDF files,
    assuming one file per calendar month, and handling adjacent files as needed
    for cumulative variables.
process_mask_and_nulls
    Detect and process masked non-null values and unmasked null values in the
    source ERA5 Land data, based on a provided mask file. Warns and/or raises
    errors for these values based on the provided parameters, and can optionally
    write a file with the locations of unmasked null values, and call
    `process_unmasked_nulls` to attempt to fill unmasked null values.
process_unmasked_nulls
    Attempt to fill unmasked null values in the source ERA5 Land data, given a
    provided mask.

Attributes
----------
LOCATION_DIM_ATTR_NAME : str
    The name of the Dataset / netCDF file attribute that holds the name of the
    flattened location dimension, for Datasets and files where the data is
    flattened from a multi-dimensional time/latitude/longitude grid into a
    single location dimension with latitude/longitude coordinates and times in
    one-dimensional coordinate variables. The name of the dimension may have to
    vary to avoid collisions with existing dimensions and variable names in the
    dataset, and is therefore stored in an attribute rather than being fixed.
"""
from collections import Counter
from collections.abc import (
    Callable,
    Hashable,
    Iterator,
    Mapping,
    Sequence,
)
import copy
import enum
import itertools
import logging
from pathlib import Path
import time
import typing as tp
from typing import NamedTuple
import warnings

from korsbakken_python_utils.containers.dataobject import UniformTypeDataObject
import numpy as np
import xarray as xr

from .convert_data import (
    era5land_from_linear_time,
    era5land_to_linear_time,
    make_datm_ds,
    round_coords,
)
from .datm_streams import Datm7Stream
from .dimensions import (
    Datm7Dim,
    ERA5_LINEARIZED_TIME_DIM,
    Era5LandDim,
    ERA5LandTimeLayout,
)
from .file_io import (
    open_era5land_grib,
    write_datm_nc,
)
from .file_name_parsing import (
    DuplicateFilesError,
    FilesAlreadyExistError,
    FilesNotFoundError,
    resolve_file_paths,
)
from .logger_registry import register_logger
from .masking import (
    ERA5_LAND_TO_DATM_MASK_VAR,
    Era5LandToDatmMaskFileDim,
    MaskedValuesHandling,
    UnmaskedNullsHandling,
    UnmaskedNullsProcessing,
)
from .types import YearMonth



LOCATION_DIM_ATTR_NAME: tp.Final[str] = "location_dim"

logger: logging.Logger = logging.getLogger(__name__)
register_logger(logger)


def convert_era5_file(
        *,
        source_file: Path,
        next_source_file: Path|None = None,
        previous_source_file: Path|None = None,
        output_streams: Sequence[Datm7Stream | str] = tuple(Datm7Stream),
        output_files: (
            Mapping[Datm7Stream, Path]
            | Callable[[Datm7Stream], Path]
            | None
        ) = None,
        keep_first_last_dates: bool = False,
        round_lat_to: float|None = None,
        round_lon_to: float|None = None,
        mask_file: Path|str|None = None,
        if_masked_values: MaskedValuesHandling = MaskedValuesHandling.RAISE,
        if_unmasked_nulls: UnmaskedNullsHandling = UnmaskedNullsHandling.WARN,
        unmasked_nulls_processing: UnmaskedNullsProcessing = \
            UnmaskedNullsProcessing.NONE,
        null_value_file: Path|str|None = None,
        eager: bool = True,
        disable_dask: bool = False,
) -> None:
    """Converts an ERA5 Land GRIB file to a DAMT7 threestream netCDF file.

    Parameters
    ----------
    source_file : Path
        Path to the ERA5 Land GRIB file for the current month or other time
        period.
    next_source_file : Path | None, optional
        Path to the ERA5 Land GRIB file for the subsequent time period. Some
        variables require data from the first time step of the next period,
        notably radiation and precipitation variables, which are cumulative in
        the ERA5 land data, but intensity-based in DATM7. If None, these
        variables will have missing values in the last time step of the output.
        Optional, by default None.
    previous_source_file : Path | None, optional
        Path to the ERA5 Land GRIB file for the previous time period. This file
        may be needed for the same reasons as `next_source_file`, namely that
        cumulative variables need to be differenced, which requires having the
        last time step of the previous period available. If None, the first time
        step of the output will have missing values for these variables.
    output_streams : Sequence[Datm7Stream | str], optional
        Sequence of Datm7Stream enum values or their string representations,
        indicating which output streams to generate. Optional, by default all
        streams.
    output_files : Mapping[Datm7Stream, Path] | Callable[[Datm7Stream], Path] | None, optional
        Output files for the requested DATM7 streams. Either as a mapping from
        Datm7Stream enums to Path objects, or as a callable that takes a
        Datm7Stream enum as its only argument and returns the Path for that
        stream. Mappings must contain entries for all requested streams, and
        callables must return valid Paths for all requested streams. Optional,
        None by default, in which case output files will be equal to the source
        file with an underscore and the stream ID appended, before a `.nc`
        extension.
    keep_first_last_dates : bool, optional
        Whether to keep the first and last dates from the source ERA5 Land
        dataset in the output DATM7 datasets. By default False. Normally, the
        one-month ERA5 Land files contain the last date of the previous month
        and the first timestep (midnight) of the first day in the next month,
        but with null values. These are normally discarded in the conversion in
        order to ensure that the output files only contain data for a single
        calendar month, and don\'t introduce overlaps in the file streams. Set
        this flag to keep these dates in the output.
    round_lat_to : float | None, optional
        If provided, round latitude values to the nearest multiple of this
        value. This can be useful to ensure that the latitude values in the
        output files are exactly the same as the expected latitude values in the
        DATM7 files, which are multiples of 0.1 degrees. The coordinates in many
        of the ERA5 Land files are very slightly off from multiples of 0.1
        degrees, and this can lead to issues when comparing or merging with
        other data. By default, no rounding is applied.
    round_lon_to : float | None, optional
        If provided, round longitude values to the nearest multiple of this
        value. See `round_lat_to` for details. By default, no rounding is
        applied.
    mask_file : Path | str, optional
        Path to a netCDF file containing a mask for the ERA5 Land grid points.
        See the documentation for the `mask_file` parameter in the
        `convert_monthly_era5_files` function for details.
    if_masked_values : MaskedValuesHandling, optional
        How to handle non-null values found in the masked areas. See the
        documentation for the `if_masked_values` parameter in the
        `convert_monthly_era5_files` function for details.
    if_unmasked_nulls : UnmaskedNullsHandling, optional
        How to handle null values found in the unmasked areas. See the
        documentation for the `if_unmasked_nulls` parameter in the
        `convert_monthly_era5_files` function for details.
    unmasked_nulls_processing : UnmaskedNullsProcessing, optional
        How to process null values found in the unmasked areas. See the
        documentation for the `unmasked_nulls_processing` parameter in the
        `convert_monthly_era5_files` function for details.
    null_value_file : Path | str, optional
        File in which to save the locations of unmasked null values. The output
        file will contain boolean variables with the same names as the ERA5 Land
        variables that were processed, and will have the value True wherever
        unmasked null values were found in the original data, and False
        elsewhere. This parameter is ignored if `mask_file` is not provided. It
        *is* used if `if_unmasked_nulls` is set to UnmaskedNullsHandling.IGNORE
        and if `unmasked_nulls_processing` is not UnmaskedNullsProcessing.NONE.
        Those parameters control logging/reporting and filling of unmasked null
        values, not whether to output a file with null value points.
    eager : bool, optional
        Whether to load the entire dataset into memory before processing. This
        can significantly speed up some operations, while setting it to False
        may lead require less memory usage but lead to some data being reloaded
        and even reprocesseed multiple times, which can slow down the
        processing. By default True. Set to False if you run into memory issues.
    disable_dask: bool = False, optional
        Whether to disable dask chunking when opening the ERA5 Land GRIB file.
        This may be useful if you have enough memory to hold the entire dataset
        in memory and don't want to install dask. But note that it may lead to
        less parallalism and higher memory usage. By default False.
    """
    if not disable_dask:
        from korsbakken_python_utils.dask_utils import (
            get_registered_progress_bars,
            set_global_progress_bar,
        )
        if len(get_registered_progress_bars()) == 0:
            set_global_progress_bar()
    if any(
        _round_to is not None and (
            not isinstance(_round_to, (int, float))
            or _round_to <= 0
        ) for _round_to in (round_lat_to, round_lon_to)
    ):
        raise ValueError(
            '`round_lat_to` and `round_lon_to` must be positive, non-zero '
            'numbers if provided.'
        )
    source_file = Path(source_file)
    output_streams = [Datm7Stream(_s) for _s in output_streams]
    output_files_mapping: dict[Datm7Stream, Path]
    if output_files is None:
        output_files_mapping = {
            _stream: source_file.with_name(
                source_file.stem + f'_{_stream.value}.nc'
            )
            for _stream in output_streams
        }
    elif callable(output_files):
        output_files_mapping = {
            _stream: output_files(_stream)
            for _stream in output_streams
        }
    else:
        output_files_mapping = dict(output_files)
    missing_output_streams = set(output_streams) - set(
        output_files_mapping.keys()
    )
    if missing_output_streams:
        raise ValueError(
            'The `output_files` parameter is missing entries for the following '
            f'requested streams: {missing_output_streams}'
        )
    if mask_file is not None:
        try:
            mask_file = Path(mask_file)
        except TypeError as _type_err:
            mask_file_str: str = str(mask_file)
            # Shorten the mask_file_str so that it won't flood the logs
            if len(mask_file_str) > 512:
                mask_file_str = (
                    mask_file_str[:256] + '... (truncated) ...' + mask_file_str[-256:]
                )
            error_msg: str = (
                'Received a value for `mask_file` that could not be converted '
                f'to a valid Path object: {mask_file_str}'
            )
            logger.error(
                msg=error_msg,
                exc_info=_type_err,
                extra={'mask_file_value': mask_file},
            )
            raise TypeError(error_msg) from _type_err
        if not mask_file.exists():
            raise FileNotFoundError(
                f'The specified mask file does not exist: {mask_file!s}'
            )
        if (
                (if_masked_values == MaskedValuesHandling.IGNORE)
                and (if_unmasked_nulls == UnmaskedNullsHandling.IGNORE)
                and (unmasked_nulls_processing == UnmaskedNullsProcessing.NONE)
                and (null_value_file is None)
        ):
            logger.warning(
                'A mask file was provided, but all options for handling masked '
                'non-null values and unmasked null values are set to ignore or '
                'none, and no null value file is specified. This means that the '
                'mask file will not actually be used for anything. Please check '
                'your parameters and make sure this is what you intended.'
            )

    logger.info(
        f'Opening ERA5 Land GRIB files:\n'
        f'  Main file: {source_file!s}\n'
        f'  Next timestep file: {next_source_file!s}\n'
        f'  Previous timestep file: {previous_source_file!s}\n'
    )
    source_ds: xr.Dataset = open_era5land_grib(
        file=source_file,
        next_file=next_source_file,
        previous_file=previous_source_file,
        chunks='auto' if not disable_dask else None,
    )
    logger.info('Finished opening ERA5 Land GRIB files.')
    latlon_rounding_map: dict[Hashable, float] = {
        _dim: _round_to
        for _dim, _round_to in zip(
            (Era5LandDim.LAT, Era5LandDim.LON),
            (round_lat_to, round_lon_to),
        )
        if _round_to is not None
    }
    if len(latlon_rounding_map) > 0:
        logger.info(
            'Rounding latitude and longitude values to nearest:\n'
            '    ' + ',  '.join(
                f'{_dim}: {_round_to}'
                for _dim, _round_to in latlon_rounding_map.items()
            )
            + 'The mask file will be rounded in the same way if provided.'
        )
        source_ds = round_coords(
            source_ds,
            round_to=latlon_rounding_map,
        )
    if eager:
        logger.info('Loading dataset into memory...')
        if not disable_dask:
            source_ds = source_ds.persist()
        else:
            source_ds = source_ds.load()
        logger.info('Finished loading dataset into memory.')
    else:
        logger.info(
            'Lazy loading requested, not loading datasets into memory yet.'
        )
    if mask_file is not None:
        mask_ds: xr.Dataset = xr.open_dataset(
            mask_file,
            chunks='auto' if not disable_dask else None,
        )
        if len(latlon_rounding_map) > 0:
            mask_ds = round_coords(
                mask_ds,
                round_to=latlon_rounding_map,
            )
        unmasked_null_ds: xr.Dataset | None
        masked_nonnull_ds: xr.Dataset | None
        source_ds = process_mask_and_nulls(
            source=source_ds,
            mask=mask_ds,
            if_masked_values=if_masked_values,
            if_unmasked_nulls=if_unmasked_nulls,
            unmasked_nulls_processing=unmasked_nulls_processing,
            null_value_file=null_value_file,
            source_files=(source_file,),
        ).filled_ds  # Should return a namedtuple with processed dataset (filled
                     # nulls) as one of the fields, and we only need to keep
                     # that here.

    logger.info('Converting and writing output DATM7 netCDF files...')
    for _target_stream in output_streams:
        logger.info(f'  Processing stream {_target_stream.value}...')
        _target_ds: xr.Dataset = make_datm_ds(
            source=source_ds,
            target_stream=_target_stream,
            eager=eager,
        )
        logger.info(
            f'    Writing output for stream {_target_stream.value} to '
            f'{output_files_mapping[_target_stream]}...'
        )
        if not keep_first_last_dates:
            time_arr: xr.DataArray = _target_ds[Datm7Dim.TIME]
            _target_ds = _target_ds.sel(
                {
                    Datm7Dim.TIME: (
                        time_arr.dt.date != time_arr.dt.date.min()
                    ) & (
                        time_arr.dt.date != time_arr.dt.date.max()
                    )
                }
            )
        write_datm_nc(
            _target_ds,
            output_file=output_files_mapping[_target_stream],
            stream=_target_stream,
        )
        logger.info(f'    Finished writing {_target_stream.value}.')
    logger.info('Finished converting ERA5 Land GRIB file.')
###END def convert_era5_file


class _YearMonthFilepathFunction(tp.Protocol):
    """Protocol for functions to make file paths from year and month."""
    def __call__(self, year: int, month: int) -> Path: ...
###END class _YearMonthFilepathFunction

class _YearMonthStreamFilepathFunction(tp.Protocol):
    """Protocol for functions to make file paths from year, month and stream."""
    def __call__(self, year: int, month: int, stream: Datm7Stream) -> Path: ...
###END class _YearMonthStreamFilenameFunction

def convert_monthly_era5_files(
        *,
        source_files: Sequence[Path|str] | _YearMonthFilepathFunction | str,
        output_files: (
            Sequence[Mapping[Datm7Stream, Path|str]]
            | _YearMonthStreamFilepathFunction
            | str
        ),
        start_year_month: tuple[int, int] | None = None,
        end_year_month: tuple[int, int] | None = None,
        next_source_file: Path|str|None = None,
        previous_source_file: Path|str|None = None,
        round_lat_to: float|None = None,
        round_lon_to: float|None = None,
        mask_file: Path|str|None = None,
        if_masked_values: MaskedValuesHandling = MaskedValuesHandling.RAISE,
        if_unmasked_nulls: UnmaskedNullsHandling = UnmaskedNullsHandling.WARN,
        unmasked_nulls_processing: UnmaskedNullsProcessing = UnmaskedNullsProcessing.NONE,
        null_value_files: _YearMonthFilepathFunction | str | None = None,
) -> None:
    """Converts multiple monthly ERA5 Land GRIB files to DATM7 threestream netCDF files.

    Parameters
    ----------
    source_files : Sequence[Path | str] | Callable | str
        Either a sequence of file paths to the ERA5 Land GRIB files to convert,
        or a function that takes year, month and stream as parameters with the
        names `year`, `month` and `stream`, in that order, and returns the
        corresponding file path. Can also be a format string with parameters
        `year`, `month` and `stream`, in which case the file paths will be
        generated by calling `.format(year=year, month=month, stream=stream)`.
        If a sequence is given (other than a scalar string), the elements can be
        either Path instances or strings, which will be connverted to Path.
        **NB!** A string parameter will *not* be interpreted as a single file path.
        Also note that absolute paths should be used. The basis for relative
        paths may not be predictable.
    output_files : Sequence[dict[Datm7Stream, Path|str]] | Callable | str
        Output files, either as a sequence of dictionaries mapping Datm7Stream
        (or their string equivalents) to file paths, or as a function that takes
        year, month and stream as parameters with the names `year`, `month` and
        `stream`, in that order, and returns the corresponding file path. Can
        also be a format string with parameters `year`, `month` and `stream`, in
        which case the file paths will be generated by calling
        `.format(year=year, month=month, stream=stream)`. If a sequence of
        dictionaries is given, it must have the same length as `source_files`.
    start_year_month : tuple[int, int], optional
        Tuple of (year, month) for the first month to convert, inclusive. If
        `source_files` and `output_files` are sequences, this parameter is
        ignored and can be left unspecified (None). If `source_files` and
        `output_files` are functions or format strings, this parameter must be
        specified, or a ValueError will be raised.
    end_year_month : tuple[int, int], optional
        Tuple of (year, month) for the last month to convert, inclusive. If
        `source_files` and `output_files` are sequences, this parameter is
        ignored and can be left unspecified (None). If `source_files` and
        `output_files` are functions or format strings, this parameter must be
        specified, or a ValueError will be raised.
    next_source_file : Path | str, optional
        Path to the ERA5 Land GRIB file for the month after `end_year_month`.
        This is optional, but some variables may have missing data in the last
        time step in the last output file if this is not provided.
    previous_source_file : Path | str, optional
        Path to the ERA5 Land GRIB file for the month before `start_year_month`.
        This is optional, but some variables may have missing data in the first
        time step in the first output file if this is not provided.
    round_lat_to : float, optional
        If provided, round latitude values to the nearest multiple of this
        value. This can be useful to ensure that the latitude values in the
        output files are exactly the same as the expected latitude values in the
        DATM7 files, which are multiples of 0.1 degrees. The coordinates in many
        of the ERA5 Land files are very slightly off from multiples of 0.1
        degrees, and this can lead to issues when comparing or merging with
        other data. By default, no rounding is applied.
    round_lon_to : float, optional
        If provided, round longitude values to the nearest multiple of this
        value. See `round_lat_to` for details. By default, no rounding is
        applied.
    mask_file : Path | str, optional
        Path to a netCDF file containing a mask for the ERA5 Land grid points.
        This is used to check for non-null values in the masked areas and null
        values in the unmasked areas, which may indicate issues with the source
        data. The mask file should have a boolean variable named `mask` with
        dimensions `latitude` and `longitude`, where True indicates unmasked
        points that should have valid data values, and False indicates masked
        points that should all have null values. The mask file must be on the
        same grid as the ERA5 Land files, and the latitude and longitude values
        must match those in the ERA5 Land files after rounding with the
        `round_lat_to` and `round_lon_to` parameters if specified. The mask file
        can be created using the `create_era5land_to_datm_mask_file.py` script
        based on an existing mask or a data file with the same grid as the ERA5
        Land files and representative values. **NB!** The mask file is *not*
        used to actually mask out values in the output files, and is not
        included as a mask variable in the output files. It is only used for
        checking the source data and reporting masked non-null and unmasked null
        values.
    if_masked_values : MaskedValuesHandling, optional
        How to handle non-null values found in the masked areas (where the mask
        is False). By default, MaskedValuesHandling.RAISE, which will raise a
        ValueError if any non-null values are found in the masked areas. If set
        to MaskedValuesHandling.WARN, a warning will be logged for each variable
        that has non-null values in the masked areas, but no error will be
        raised. If set to MaskedValuesHandling.IGNORE, any non-null values in
        the masked areas will be ignored and no warnings or errors will be
        raised. This parameter is ignored if `mask_file` is not provided.
    if_unmasked_nulls : UnmaskedNullsHandling, optional
        How to handle null values found in the unmasked areas (where the mask is
        True). By default, UnmaskedNullsHandling.WARN, which will log a warning
        for each variable that has null values in the unmasked areas, but will
        not raise an error. If set to UnmaskedNullsHandling.RAISE, a ValueError
        will be raised if any null values are found in the unmasked areas. If set
        to UnmaskedNullsHandling.IGNORE, any null values in the unmasked areas
        will be ignored and no warnings or errors will be raised. This parameter
        is ignored if `mask_file` is not provided.
    unmasked_nulls_processing : UnmaskedNullsProcessing, optional
        If any null values are found in the unmasked areas, this parameter
        controls whether to attempt to fill these values, and if so, how. By
        default UnmaskedNullsProcessing.NONE, which means that no filling will
        be attempted. If set to UnmaskedNullsProcessing.LINEAR, the null values
        in the unmasked areas will be filled using 1D linear interpolation along
        the time dimension. No other options are currently supported. This
        parameter is ignored if `mask_file` is not provided. It *is* used if
        `if_unmasked_nulls` is set to UnmaskedNullsHandling.IGNORE, since that
        option only controls whether to log and report the presence of unmasked
        nulls to the user.
    null_value_files : Callable | None, optional
        Files in which to save the locations of unmasked null values for each
        source file. This parameter must be provided as a function that takes
        year and month as input and returns a valid file name string, or a
        format string following the same formatting rules as `source_file`. The
        output files will contain boolean variables with the same names as the
        ERA5 Land variables that were processed, and will have the value True
        wherever unmasked null values were found in the original data, and False
        elsewhere. This parameter is ignored if `mask_file` is not provided. It
        *is* used if `if_unmasked_nulls` is set to UnmaskedNullsHandling.IGNORE
        and if `unmasked_nulls_processing` is not UnmaskedNullsProcessing.NONE.
        Those parameters control logging/reporting and filling of unmasked null
        values, respectively, not whether to output a file with null value
        points.
    """
    # TODO: The parsing of source_files and output_files should be refactored to
    # use the resolve_file_paths function.
    if isinstance(source_files, str):
        _source_files: str = copy.copy(source_files)
        def _source_files_func(year: int, month: int) -> Path:
            return Path(_source_files.format(year=year, month=month))
        source_files = _source_files_func
    if callable(source_files):
        source_files_func: _YearMonthFilepathFunction|None = source_files
        if not (
            isinstance(start_year_month, tuple)
            and isinstance(end_year_month, tuple)
            and len(start_year_month) == 2
            and len(end_year_month) == 2
        ):
            raise ValueError(
                'If `source_files` is a function or format string, '
                '`start_year_month` and `end_year_month` must be specified as '
                '(year, month) tuples.'
        )
        use_source_files: list[Path] = [
            source_files_func(year=_year, month=_month)
            for _year, _month in YearMonth.range(
                YearMonth(*start_year_month),
                YearMonth(*end_year_month),
            )
        ]
    elif isinstance(source_files, Sequence):
        use_source_files = [Path(_f) for _f in source_files]
        source_files_func = None
    else:
        raise TypeError(
            '`source_files` must be either a sequence of file paths, a '
            'function, or a format string.'
        )
    missing_source_files: list[Path] = [
        _f
        for _f in use_source_files
        if not _f.exists()
    ]
    if len(missing_source_files) > 0:
        raise FileNotFoundError(
            'The following source files do not exist:\n' + ',\n'.join(
                str(_f) for _f in missing_source_files
            )
        )
    logger.debug(
        f'Using source files: ' + ', '.join(str(_f) for _f in use_source_files)
    )
    if isinstance(output_files, str):
        _output_files: str = copy.copy(output_files)
        def _output_files_func(year: int, month: int, stream: Datm7Stream) -> Path:
            return Path(_output_files.format(year=year, month=month, stream=stream))
        output_files = _output_files_func
    if callable(output_files):
        if not (
            isinstance(source_files, str)
            or callable(source_files)
        ):
            raise ValueError(
                'If `output_files` is a function or format string, '
                '`source_files` must also be a function or format string.'
            )
        assert isinstance(start_year_month, tuple)
        assert isinstance(end_year_month, tuple)
        output_files_func: _YearMonthStreamFilepathFunction = output_files
        use_output_files: list[dict[Datm7Stream, Path]] = [
            {
                _stream: output_files_func(
                    year=_year,
                    month=_month,
                    stream=_stream,
                ).resolve()
                for _stream in Datm7Stream
            }
            for _year, _month in YearMonth.range(
                YearMonth(*start_year_month),
                YearMonth(*end_year_month),
            )
        ]
    elif isinstance(output_files, Sequence):
        use_output_files = [
            {
                Datm7Stream(_stream): Path(_path).resolve()
                for _stream, _path in _mapping.items()
            }
            for _mapping in output_files
        ]
        if len(use_output_files) != len(use_source_files):
            raise ValueError(
                'If `output_files` is a sequence, it must have the same length '
                'as `source_files`.'
            )
    else:
        raise TypeError(
            '`output_files` must be either a sequence of file paths, a '
            'function, or a format string.'
        )
    # Count occurrences of the same file in use_output_files to avoid duplicates
    output_file_counts: Counter[Path] = Counter(
        itertools.chain.from_iterable(
            _mapping.values() for _mapping in use_output_files
        )
    )
    duplicate_output_files: list[Path] = [
        _file
        for _file, _count in output_file_counts.items()
        if _count > 1
    ]
    if len(duplicate_output_files) > 0:
        raise ValueError(
            'The following output files are specified more than once:\n'
            + ',\n'.join(str(_f) for _f in duplicate_output_files)
        )
    # Check for existing output files
    existing_output_files: list[Path] = [
        _f for _f in output_file_counts.keys() if _f.exists()
    ]
    if len(existing_output_files) > 0:
        raise FileExistsError(
            'The following output files already exist, please move or delete '
            'them, or choose different names:\n' + ',\n'.join(
                str(_f) for _f in existing_output_files
            )
        )
    logger.debug(
        f'Using output files: ' + ', '.join(
            str(_f) for _mapping in use_output_files for _f in _mapping.values()
        )
    )
    if (
            previous_source_file is None
            and source_files_func is not None
            and start_year_month is not None
    ):
        previous_source_file = source_files_func(
            year=start_year_month[0] if start_year_month[1] > 1 else start_year_month[0] - 1,
            month=start_year_month[1] - 1 if start_year_month[1] > 1 else 12,
        )
        if not previous_source_file.exists():
            previous_source_file = None
    elif isinstance(previous_source_file, str):
        previous_source_file = Path(previous_source_file)
    logger.debug(
        f'Using previous source file: {previous_source_file!s}'
        if previous_source_file is not None
        else 'No previous source file specified.'
    )
    if (
            next_source_file is None
            and source_files_func is not None
            and end_year_month is not None
    ):
        next_source_file = source_files_func(
            year=end_year_month[0] if end_year_month[1] < 12 else end_year_month[0] + 1,
            month=end_year_month[1] + 1 if end_year_month[1] < 12 else 1,
        )
        if not next_source_file.exists():
            next_source_file = None
    elif isinstance(next_source_file, str):
        next_source_file = Path(next_source_file)
    logger.debug(
        f'Using next source file: {next_source_file!s}'
        if next_source_file is not None
        else 'No next source file specified.'
    )
    if null_value_files is not None:
        if start_year_month is None or end_year_month is None:
            null_files_field_values: dict[str, tuple[int, ...]] = {}
        else:
            year_months: Iterator[YearMonth] = YearMonth.range(
                YearMonth(*start_year_month),
                YearMonth(*end_year_month),
            )
            year_values: tuple[int, ...]
            month_values: tuple[int, ...]
            year_values, month_values = zip(*year_months)
            null_files_field_values: dict[str, tuple[int, ...]] = {
                'year': year_values,
                'month': month_values,
            }
        try:
            use_null_value_files: Sequence[Path|None] = resolve_file_paths(
                paths=null_value_files,
                field_values=null_files_field_values,
                check_not_exists=True,
                check_duplicates=True,
            )
        except TypeError as _type_err:
            error_msg: str = (
                'Encountered an error when trying to generate null '
                'value file paths from the provided `null_value_files` '
                'argument. It is likely that you provided a format or '
                'function that does not have "year" and "month" as '
                'parameters, but other issues can be the cause. Please '
                'check the stack trace for details.'
            )
            logger.error(
                msg=error_msg,
                exc_info=_type_err,
            )
            raise TypeError(error_msg) from _type_err
        except FilesAlreadyExistError as _exist_err:
            error_msg = (
                'Encountered an error when trying to generate null value '
                'file paths from the provided `null_value_files` argument. '
                'Some of the generated file paths already exist. Please '
                'move or delete these files, or choose a different format or '
                'function for `null_value_files`. The existing files are:\n'
                + ',\n'.join(str(_f) for _f in _exist_err.existing_files)
            )
            logger.error(
                msg=error_msg,
                exc_info=_exist_err,
                extra={
                    'existing_files': _exist_err.existing_files,
                },
            )
            raise FilesAlreadyExistError(
                error_msg,
                existing_files=_exist_err.existing_files,
            ) from _exist_err
        except DuplicateFilesError as _dup_err:
            error_msg = (
                'Encountered an error when trying to generate null value '
                'file paths from the provided `null_value_files` argument. '
                'Some of the generated file paths are duplicates. Please '
                'check your format or function for `null_value_files` to '
                'ensure that it generates unique file paths for each month. '
                'The duplicate files are:\n'
                + ',\n'.join(str(_f) for _f in _dup_err.duplicate_files)
            )
            logger.error(
                msg=error_msg,
                exc_info=_dup_err,
                extra={
                    'duplicate_files': _dup_err.duplicate_files,
                    'file_counts': _dup_err.file_counts,
                },
            )
            raise DuplicateFilesError(
                error_msg,
                duplicate_files=_dup_err.duplicate_files,
                file_counts=_dup_err.file_counts,
            ) from _dup_err
    else:
        use_null_value_files: Sequence[Path|None] = [None]*len(use_source_files)
    for (
            _previous_source,
            _source,
            _next_source,
            _output_mapping,
            _null_value_file,
    ) in zip(
            (previous_source_file, *use_source_files[:-1]),
            use_source_files,
            (*use_source_files[1:], next_source_file),
            use_output_files,
            use_null_value_files,
    ):
        # Ignore the mask_file if no options are set to use it, to avoid
        # generating a warning for every single monthly file conversion.
        if mask_file is not None:
            if (
                    (if_masked_values == MaskedValuesHandling.IGNORE)
                    and (if_unmasked_nulls == UnmaskedNullsHandling.IGNORE)
                    and (unmasked_nulls_processing == UnmaskedNullsProcessing.NONE)
                    and (_null_value_file is None)
            ):
                logger.warning(
                    'A mask file was provided, but all options for handling '
                    'masked non-null values and unmasked null values are set '
                    'to ignore or none, and no null value file is specified. '
                    'This means that the mask file will not actually be used '
                    'for anything. Please check your parameters and make sure '
                    'this is what you intended.'
                )
                use_mask_file: Path|str|None = None
            else:
                use_mask_file = mask_file
        else:
            use_mask_file = mask_file
        logger.info(
            f'Converting files:\n'
            f'  source_file: {_source!s}\n'
            f'  next_source_file: {_next_source!s}\n'
            f'  previous_source_file: {_previous_source!s}\n'
            f'  output_files: {_output_mapping!s}\n'
        )
        convert_era5_file(
            source_file=_source,
            next_source_file=_next_source,
            previous_source_file=_previous_source,
            output_streams=list(_output_mapping.keys()),
            output_files=_output_mapping,
            round_lat_to=round_lat_to,
            round_lon_to=round_lon_to,
            mask_file=use_mask_file,
            if_masked_values=if_masked_values,
            if_unmasked_nulls=if_unmasked_nulls,
            unmasked_nulls_processing=unmasked_nulls_processing,
            null_value_file=_null_value_file,
        )
###END def convert_monthly_era5_files


class ProcessMaskAndNullsResult(NamedTuple):
    """Named tuple for the result of process_mask_and_nulls.

    This class is a subclass of NamedTuple. The fields below are listed in the
    order they appear in the tuple.

    Fields
    ------
    filled_ds : xr.Dataset
        The dataset resulting from processing unmasked null values in the source
        Dataset passed to `process_mask_and_nulls`.
    unmasked_null_ds : xr.Dataset | None
        If `if_unmasked_nulls` is not UnmaskedNullsHandling.IGNORE or if
        `null_value_file` is not None, this will be a Dataset with boolean
        values indicating the locations of unmasked null values in the original
        source dataset, with the same variable names and dimensions as the
        original source dataset. If `if_unmasked_nulls` is
        UnmaskedNullsHandling.IGNORE, this will be None.
    masked_nonnull_ds : xr.Dataset | None
        If `if_masked_values` is not MaskedValuesHandling.IGNORE, this will be a
        Dataset with boolean values indicating the locations of non-null values
        in the masked areas in the original source dataset, with the same
        variable names and dimensions as the original source dataset. If
        `if_masked_values` is MaskedValuesHandling.IGNORE, this will be None.
    """
    filled_ds: xr.Dataset
    unmasked_null_ds: xr.Dataset | None
    masked_nonnull_ds: xr.Dataset | None
###END class ProcessMaskAndNullsResult


def process_mask_and_nulls(
        *,
        source: xr.Dataset,
        mask: xr.Dataset,
        if_masked_values: MaskedValuesHandling,
        if_unmasked_nulls: UnmaskedNullsHandling,
        unmasked_nulls_processing: UnmaskedNullsProcessing,
        null_value_file: Path|str|None,
        source_files: Sequence[Path|str] | None = None,
        source_time_layout: ERA5LandTimeLayout = ERA5LandTimeLayout.DATE_STEP,
) -> ProcessMaskAndNullsResult:
    """Checks for masked non-null and unmasked null values in the source dataset
    and optionally fills unmasked null values.

    This function performs the mask and null-value related checks and processing
    for the `convert_era5_file` and `convert_monthly_era5_files` functions.

    Parameters
    ----------
    source : xr.Dataset
        The source dataset to check and process. This should be the dataset
        opened from the ERA5 Land GRIB file, after any coordinate rounding has
        been applied.
    mask : xr.Dataset
        The mask dataset. See the `mask_file` parameter in `convert_era5_file`
        for details.
    if_masked_values : MaskedValuesHandling
        How to handle non-null values found in the masked areas. See the
        `if_masked_values` parameter in `convert_era5_file` for details.
    if_unmasked_nulls : UnmaskedNullsHandling
        How to handle null values found in the unmasked areas. See the
        `if_unmasked_nulls` parameter in `convert_era5_file` for details.
    unmasked_nulls_processing : UnmaskedNullsProcessing
        If any null values are found in the unmasked areas, this parameter
        controls whether to attempt to fill these values, and if so, how. See
        the `unmasked_nulls_processing` parameter in `convert_era5_file` for
        details.
    null_value_file : Path | str | None
        Files in which to save the locations of unmasked null values in the
        source file. See the `null_value_file` parameter in `convert_era5_file`
        for details.
    source_files : Sequence[Path|str] | None, optional
        If provided, this should be a sequence of file paths to the source files
        being processed. This is only used for logging and error messages, to
        provide more context to the user about which file is being processed and
        which files have issues with masked non-null or unmasked null values. If
        not provided, the log messages will not include file-specific context.
    source_time_layout : ERA5LandTimeLayout, optional
        The time layout of the source dataset. This is needed to correctly
        handle the time dimension when checking for masked non-null and unmasked
        null values, and when filling unmasked null values if requested. By
        default, ERA5LandTimeLayout.DATE_STEP, which means that the time
        dimension in the source dataset is expected to be in the format of a
        date dimension (e.g. 'time') and a step dimension (e.g. 'step') that
        indicates the time offset from the date. If set to
        ERA5LandTimeLayout.LINEAR, the time dimension is expected to be a single
        linearized time dimension (e.g. 'time') that combines the date and step
        information. If called from the `convert_era5_file` function, the file
        will usually be in `DATE_STEP` layout (the default) and will later be
        linearized in the `make_datm_ds` function.

    Returns
    -------
    ProcessMaskAndNullsResult
        A NamedTuple object with the results. See the documentation of
        `ProcessMaskAndNullsResult` for information about the fields.
    """
    coord_tolerance: float = 1e-15  # Tolerance for coordinate alignment between source and mask datasets, in coordinate units (degrees).
    return_unmasked_null_ds: xr.Dataset | None
    return_masked_nonnull_ds: xr.Dataset | None
    if_masked_values = MaskedValuesHandling(if_masked_values)
    if_unmasked_nulls = UnmaskedNullsHandling(if_unmasked_nulls)
    unmasked_nulls_processing = (
        UnmaskedNullsProcessing(unmasked_nulls_processing)
    )
    logger.info(
        msg = (
            'Processing mask and null values for source files: '
            + (
                ', '.join(str(_f) for _f in source_files)
                if source_files is not None
                else 'N/A (source file names not provided)'
            )
        ),
        extra={
            'source_files': source_files,
        }
    )
    # Shortcut the processing if the parameters are trivial
    if (
            if_masked_values == MaskedValuesHandling.IGNORE
            and if_unmasked_nulls == UnmaskedNullsHandling.IGNORE
            and unmasked_nulls_processing == UnmaskedNullsProcessing.NONE
            and null_value_file is None
    ):
        return ProcessMaskAndNullsResult(
            filled_ds=source,
            unmasked_null_ds=None,
            masked_nonnull_ds=None,
        )
    # Check that the mask dataset has the expected structure and dimensions
    if ERA5_LAND_TO_DATM_MASK_VAR not in mask.variables:
        raise KeyError(
            f'The mask dataset does not contain the expected mask variable '
            f'"{ERA5_LAND_TO_DATM_MASK_VAR}".'
        )
    if set(mask[ERA5_LAND_TO_DATM_MASK_VAR].dims) != {
            Era5LandToDatmMaskFileDim.LAT,
            Era5LandToDatmMaskFileDim.LON,
    }:
        raise ValueError(
            f'The mask variable "{ERA5_LAND_TO_DATM_MASK_VAR}" in the mask '
            f'dataset does not have the expected dimensions '
            f'"{Era5LandToDatmMaskFileDim.LAT}" and '
            f'"{Era5LandToDatmMaskFileDim.LON}". Found dimensions: '
            f'{mask[ERA5_LAND_TO_DATM_MASK_VAR].dims}.'
        )
    # Crop the mask to the source dataset. Raise an error if the result contains
    # any NaN values, which indicates that the source dataset region is not
    # fully contained within the mask region. At the same time, rename the mask
    # dimensions to match the source dataset.
    try:
        mask_aligned: xr.DataArray = mask[ERA5_LAND_TO_DATM_MASK_VAR].rename(
            {
                Era5LandToDatmMaskFileDim.LAT: Era5LandDim.LAT,
                Era5LandToDatmMaskFileDim.LON: Era5LandDim.LON,
            }
        ).sel(
            {
                Era5LandDim.LAT: source[Era5LandDim.LAT],
                Era5LandDim.LON: source[Era5LandDim.LON],
            },
            method='nearest',
            tolerance=coord_tolerance,
        ).reindex_like(source, tolerance=coord_tolerance, method='nearest')
    except KeyError as _key_err:
        error_msg = (
            'Could not align the mask dataset with the source dataset. Most'
            'likely, the source data is not fully contained within the region '
            'covered by the mask data. Please check the coordinates in the '
            'mask file against the source data file.'
        )
        logger.error(
            msg=error_msg,
            exc_info=_key_err,
            extra={
                'source_lat_range': (
                    float(source[Era5LandDim.LAT].min().item()),
                    float(source[Era5LandDim.LAT].max().item()),
                ),
                'source_lon_range': (
                    float(source[Era5LandDim.LON].min().item()),
                    float(source[Era5LandDim.LON].max().item()),
                ),
                'mask_lat_range': (
                    float(mask[Era5LandToDatmMaskFileDim.LAT].min().item()),
                    float(mask[Era5LandToDatmMaskFileDim.LAT].max().item()),
                ),
                'mask_lon_range': (
                    float(mask[Era5LandToDatmMaskFileDim.LON].min().item()),
                    float(mask[Era5LandToDatmMaskFileDim.LON].max().item()),
                ),
                'source_files': source_files,
            }
        )
        raise KeyError(error_msg) from _key_err
    try:
        mask_aligned = mask_aligned.astype(bool)
    except ValueError as _value_err:
        error_msg = (
            f'The mask variable "{ERA5_LAND_TO_DATM_MASK_VAR}" could not be '
            f'converted to boolean values. Please check the values in the mask '
            f'file and ensure that they can be interpreted as boolean (e.g. 0 '
            f'and 1, or True and False).'
        )
        logger.error(
            msg=error_msg,
            exc_info=_value_err,
            extra={
                'source_files': source_files,
            }
        )
        raise ValueError(error_msg) from _value_err

    # Check for non-null values in the masked areas if needed
    if if_masked_values != MaskedValuesHandling.IGNORE:
        logger.info('Checking for non-null values in masked areas...')
        masked_nonull_ds: xr.Dataset = (
            source
            .where(~mask_aligned)  # Set to null wherever unmasked
            .notnull()  # Set all non-null values (now only in masked areas) to True.
        )
        masked_nonnull_vars: list[str] = []
        for _var in masked_nonull_ds.data_vars:
            logger.debug(
                f'Checking variable "{_var}" for non-null values in masked '
                'areas...'
            )
            if masked_nonull_ds[_var].any().compute().item():
                masked_nonnull_vars.append(str(_var))
        if len(masked_nonnull_vars) > 0:
            masked_nonnull_msg: str = (
                f'Found non-null values in masked areas for variables: '
                + ', '.join(masked_nonnull_vars)
            )
            if if_masked_values == MaskedValuesHandling.RAISE:
                logger.error(
                    msg=masked_nonnull_msg,
                    extra={
                        'masked_nonnull_vars': masked_nonnull_vars,
                        'source_files': source_files,
                    },
                )
                raise ValueError(masked_nonnull_msg)
            elif if_masked_values == MaskedValuesHandling.WARN:
                logger.warning(
                    msg=masked_nonnull_msg,
                    extra={
                        'masked_nonnull_vars': masked_nonnull_vars,
                        'source_files': source_files,
                    },
                )
                warnings.warn(masked_nonnull_msg)
            else:
                error_msg: str = (
                    'Invalid value for `if_masked_values`: '
                    f'{if_masked_values!r}. This should have been caught by '
                    'the MaskedValuesHandling enum, so this indicates a bug.'
                )
                logger.error(
                    msg=error_msg,
                    extra={
                        'source_files': source_files,
                    },
                )
                raise RuntimeError(error_msg)
        else:
            logger.info('No non-null values found in masked areas.')
        return_masked_nonnull_ds = masked_nonull_ds
    else:
        logger.debug('Skipping check for non-null values in masked areas.')
        return_masked_nonnull_ds = None

    # Check for null values in the unmasked areas if needed.
    if (
            if_unmasked_nulls != UnmaskedNullsHandling.IGNORE
            or null_value_file is not None
    ):
        logger.info('Checking for null values in unmasked areas...')
        # We expect the umasked null values to be relatively sparse, so the
        # output below will be a flattened dataset of boolean values.
        # Flatten along a new dimension, named 'location', but with underscores
        # added as needed to make the name unique.
        source_existing_ids: set[Hashable] = {*source.dims, *source.variables}
        location_dim: str = 'location'
        variable_dim: str = 'variable'
        while location_dim in source_existing_ids:
            location_dim += '_'
        while variable_dim in source_existing_ids:
            variable_dim += '_'
        unmasked_null_ds: xr.Dataset = (
            source
            .isnull()  # Set all null values to True.
            .where(mask_aligned, other=False)  # Set to False wherever masked
            .stack({location_dim: tuple(source.dims)}, create_index=False)  # Flatten all dimensions into a single location dimensions, but wait to create a MultiIndex until we have pruned the result (the coordinates should still be kept)
            .to_dataarray(dim=variable_dim)  # Convert to DataArray with a variable dimension
            .pipe(lambda da: da.sel({location_dim: da.any(dim=variable_dim)}))  # Keep only locations where at least one variable is non-null (i.e., where the value of the data at this point is True for at least one variable).
            .to_dataset(dim=variable_dim)  # Convert back to Dataset with variables as data variables
            .assign_attrs(
                {
                    'description': (
                        'Boolean, flattened dataset of coordinates where '
                        'unmasked null values were found in the original '
                        'source. Each variable has the same name as a variable '
                        'in the original source dataset, and is True where that'
                        'variable had an unmasked null value, False '
                        'elsewhere. Points where no variables had unmasked '
                        'null values are not included. The latitude, longitude '
                        'and time dimensions have been flattened into a single '
                        'location dimension. The name of the location '
                        'dimension is stored in the '
                        f'`{LOCATION_DIM_ATTR_NAME}` attribute. The original '
                        'coordinates are stored in coordinate variables with '
                        'the same names as the dimensions in the source file.'
                    ),
                    'source_files': (
                        ', '.join(str(_f) for _f in source_files)
                        if source_files is not None
                        else 'N/A (source file names not provided)'
                    ),
                    LOCATION_DIM_ATTR_NAME: location_dim,
                }
            )
        ).compute()  # Compute here to avoid doing it multiple times in the checks below.
        unmasked_null_vars: list[str] = []
        for _var in unmasked_null_ds.data_vars:
            logger.debug(
                f'Checking variable "{_var}" for null values in unmasked '
                'areas...'
            )
            if unmasked_null_ds[_var].any().compute().item():
                unmasked_null_vars.append(str(_var))
        if len(unmasked_null_vars) > 0:
            if null_value_file is not None:
                logger.info(
                    msg=(
                        'Found unmasked null values for variables: '
                        + ', '.join(unmasked_null_vars)
                        + f'. Saving locations to file: {null_value_file!s}'
                    ),
                    extra={
                        'unmasked_null_vars': unmasked_null_vars,
                        'null_value_file': null_value_file,
                        'source_files': source_files,
                    },
                )
                _write_unmasked_nulls_file(
                    ds=unmasked_null_ds,
                    path=null_value_file,
                    reset_location_dim_index=location_dim,
                )
            return_unmasked_null_ds = unmasked_null_ds
            unmasked_null_msg: str = (
                f'Found null values in unmasked areas for variables: '
                + ', '.join(unmasked_null_vars)
            )
            if if_unmasked_nulls == UnmaskedNullsHandling.RAISE:
                logger.error(
                    msg=unmasked_null_msg,
                    extra={
                        'unmasked_null_vars': unmasked_null_vars,
                        'source_files': source_files,
                    },
                )
                raise ValueError(unmasked_null_msg)
            elif if_unmasked_nulls == UnmaskedNullsHandling.WARN:
                logger.warning(
                    msg=unmasked_null_msg,
                    extra={
                        'unmasked_null_vars': unmasked_null_vars,
                        'source_files': source_files,
                    },
                )
                warnings.warn(unmasked_null_msg)
            elif if_unmasked_nulls == UnmaskedNullsHandling.IGNORE:
                pass
            else:
                error_msg: str = (
                    'Invalid value for `if_unmasked_nulls`: '
                    f'{if_unmasked_nulls!r}. This should have been caught by '
                    'the UnmaskedNullsHandling enum, so this indicates a bug.'
                )
                logger.error(
                    msg=error_msg,
                    extra={
                        'source_files': source_files,
                    },
                )
                raise RuntimeError(error_msg)
        else:
            logger.info('No null values found in unmasked areas.')
            return_unmasked_null_ds = None
    else:
        logger.debug('Skipping explicit check for null values in unmasked areas.')
        return_unmasked_null_ds = None

    # Fill unmasked null values if requested.
    if unmasked_nulls_processing != UnmaskedNullsProcessing.NONE:
        logger.info(
            msg=(
                'Processing unmasked null values using option '
                f'{unmasked_nulls_processing!s}...'
            ),
            extra={
                'unmasked_nulls_processing': unmasked_nulls_processing,
                'source_files': source_files,
            }
        )
        source = process_unmasked_nulls(
            source=source,
            mask=mask_aligned,
            unmasked_null_ds=return_unmasked_null_ds,
            processing_method=unmasked_nulls_processing,
            time_layout=source_time_layout,
        )
        logger.info(
            msg='Finished processing unmasked null values.',
            extra={
                'source_files': source_files,
            },
        )
    else:
        logger.debug(
            msg='Skipping processing of unmasked null values.',
            extra={
                'source_files': source_files,
            }
        )

    return ProcessMaskAndNullsResult(
        filled_ds=source,
        unmasked_null_ds=return_unmasked_null_ds,
        masked_nonnull_ds=return_masked_nonnull_ds,
    )

### END def process_mask_and_nulls


def _write_unmasked_nulls_file(
        ds: xr.Dataset,
        path: Path|str,
        *,
        reset_location_dim_index: str|None = None,
) -> None:
    """Writes a dataset of unmasked null value locations to a NetCDF file.

    This is used to save the locations of unmasked null values in the source
    dataset, if any are found, for later inspection by the user. The dataset
    should be in the format produced by the `process_mask_and_nulls` function,
    with boolean values indicating the locations of unmasked null values for
    each variable, and with all dimensions flattened into a single location
    dimension.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to write. This should be in the format produced by the
        `process_mask_and_nulls` function for the `unmasked_null_ds` field.
    path : Path | str
        The file path where to save the dataset. This should be a .nc file.
    reset_location_dim_index : str | None, optional
        Specify that the location dimension index should be reset, and the name
        of that dimension. This is needed when the location dimension uses a
        MultiIndex (which it usually will), since xarray cannot save
        MultiIndexes to NetCDF files.

    Returns
    -------
    None
    """
    if reset_location_dim_index is not None:
        ds = ds.reset_index(reset_location_dim_index, drop=False)
    try:
        ds.to_netcdf(path)
        logger.info(f'Successfully saved unmasked null value locations to file: {path!s}')
    except Exception as _err:
        error_msg = (
            f'An error occurred while trying to save unmasked null value '
            f'locations to file: {path!s}.'
        )
        logger.error(
            msg=error_msg,
            exc_info=_err,
            extra={
                'path': path,
            }
        )
        raise _err
### END def _write_unmasked_nulls_file


def process_unmasked_nulls(
        *,
        source: xr.Dataset,
        mask: xr.DataArray,
        processing_method: UnmaskedNullsProcessing,
        unmasked_null_ds: xr.Dataset | None,
        preserve_masked_values: bool = False,
        time_layout: ERA5LandTimeLayout = ERA5LandTimeLayout.DATE_STEP,
        source_files: Sequence[Path|str] | None = None,
) -> xr.Dataset:
    """Fills unmasked null values in the source dataset.

    This function fills unmasked null values in the source dataset using the
    specified method. By default, non-null values in the masked areas become
    nulls in the returned Dataset, but this can be modified with the
    `preserve_masked_values` parameter.

    Parameters
    ----------
    source : xr.Dataset
        The source dataset to process. This should be the dataset opened from
        the ERA5 Land GRIB file, after any coordinate rounding has been applied,
        and after aligning with the mask dataset.
    mask : xr.DataArray
        The aligned mask DataArray, with boolean values indicating the masked
        and unmasked areas. This should have the same spatial dimensions and the
        same coordinate values in the same order as `source`. True for unmasked
        points (that should have be filled with non-null values), False for
        masked points (that should have null values or not be used).
    processing_method : UnmaskedNullsProcessing
        The method to use for filling unmasked null values. See the
        `unmasked_nulls_processing` parameter in `convert_era5_file` for
        details.
    unmasked_null_ds : xr.Dataset | None
        If available, this should be the dataset of unmasked null value
        locations produced by the `process_mask_and_nulls` function. It is not
        required, but can improve performance if it has been precomputed,
        especially in cases where most of the variables do not have any unmasked
        nulls.
    preserve_masked_values : bool, optional
        If True, non-null values in the masked areas in the source dataset will
        be preserved in the returned dataset. If False (the default), masked
        points in the returned Dataset will be set to null values.
    time_layout : ERA5LandTimeLayout
        The time layout of the source dataset (whether a `time` dimension with
        dates and a `step` dimension with intra-day offsets, or a single
        linearized `time` dimension). Optional, defaults to
        `ERA5LandTimeLayout.DATE_STEP`.
    source_files : Sequence[Path|str] | None, optional
        If provided, this should be a sequence of file paths to the source files
        being processed. This is only used for logging and error messages, to
        provide more context to the user about which file is being processed and
        which files have issues with unmasked null values. If not provided, the
        log messages will not include file-specific context.
    """
    time_layout = ERA5LandTimeLayout(time_layout)
    # Shortcut the process if the method is NONE
    if processing_method == UnmaskedNullsProcessing.NONE:
        if preserve_masked_values:
            return source
        else:
            return source.where(mask)
    # Shortcut the process if the caller reports no unmasked null values
    if unmasked_null_ds is not None and (
            not any(
                unmasked_null_ds[_var].any().compute().item()
                for _var in unmasked_null_ds.data_vars
            ) or not any (
                _size > 0 for _size in unmasked_null_ds.sizes.values()
            )
    ):
        logger.debug(
            msg=(
                'No unmasked null values found in the provided unmasked_null_ds, '
                'skipping filling process.'
                + (
                    f' Source files: {list(source_files)!s}'
                    if source_files is not None
                    else ' Source files: N/A (source file names not provided)'
                )
            ),
            extra={
                'source_files': source_files,
            },
        )

    # Confirm that the mask and source datasets have compatible dimensions
    if (
            (not set(mask.dims).issubset(set(source.dims)))
            or (mask.sizes != {_dim: source.sizes[_dim] for _dim in mask.dims})
    ):
        error_msg = (
            'The mask dataset does not have compatible dimensions with the '
            'source dataset. Please check that the mask dataset has the same '
            'spatial dimensions as the source dataset, with the same coordinate '
            'values in the same order. Mask dataset dimensions and sizes: '
            f'{mask.sizes!s}. Source dataset dimensions and sizes: '
            f'{source.sizes!s}.'
        )
        logger.error(
            msg=error_msg,
            extra={
                'mask_sizes': mask.sizes,
                'source_sizes': source.sizes,
                'source_files': source_files,
            },
        )
        raise ValueError(error_msg)

    # We only support `LINEAR` for now, so fail if we received anything
    # else.
    if processing_method != UnmaskedNullsProcessing.LINEAR:
        error_msg: str = (
            'Received unknown unmasked nulls processing method: '
            f'{processing_method!s}. This should not be possible at this point '
            'in the code, and indicates either a bug or that this function has '
            'not been updated to support newly added processing methods.'
        )
        logger.error(
            msg=error_msg,
            extra={
                'processing_method': processing_method,
                'source_files': source_files,
            },
        )
        raise RuntimeError(error_msg)

    # Fill unmasked null values using linear interpolation along the time
    # dimension.
    if time_layout == ERA5LandTimeLayout.DATE_STEP:
        source = era5land_to_linear_time(
            source=source,
            preserve_source_time_coord=True,
            preserve_source_time_component_coords=True,
        )
    logger.info(
        msg=(
            'Filling unmasked null values using linear interpolation along '
            'the linearized time dimension.'
            + (
                f' Source files: {list(source_files)!s}'
                if source_files is not None
                else ' Source files: N/A (source file names not provided)'
            )
        ),
        extra={
            'processing_method': processing_method,
            'time_layout': time_layout,
            'source_files': source_files,
        },
    )
    start_process_time: float = time.process_time()
    source_filled = source.where(mask).interpolate_na(
        dim=ERA5_LINEARIZED_TIME_DIM,
        method='linear',
    )
    if preserve_masked_values:
        source_filled = source_filled.combine_first(source)
    done_filling_time: float = time.process_time()
    filling_time_consumed: float = done_filling_time - start_process_time
    logger.info(
        msg=(
            'Finished filling unmasked null values in '
            f'{filling_time_consumed/1000.0:.3G} ms.'
        )
    )
    if time_layout == ERA5LandTimeLayout.DATE_STEP:
        source_filled = era5land_from_linear_time(
            source=source_filled,
            fast_unstack=True,
        )
        done_unstacking_time: float = time.process_time()
        unstacking_time_consumed: float = (
            done_unstacking_time - done_filling_time
        )
        logger.info(
            msg=(
                'Finished unstacking linearized time dimension in '
                f'{unstacking_time_consumed/1000.0:.3G} ms.'
            )
        )
    return source_filled
### END def process_unmasked_nulls
