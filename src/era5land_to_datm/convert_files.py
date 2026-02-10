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
"""
from collections import Counter
from collections.abc import (
    Callable,
    Iterator,
    Mapping,
    Sequence,
)
import copy
import enum
import itertools
import logging
from pathlib import Path
import typing as tp

from korsbakken_python_utils.containers.dataobject import UniformTypeDataObject
import numpy as np
import xarray as xr

from .convert_data import make_datm_ds
from .datm_streams import Datm7Stream
from .dimensions import (
    Datm7Dim,
    Era5LandDim,
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
    MaskedValuesHandling,
    UnmaskedNullsHandling,
    UnmaskedNullsProcessing,
)
from .types import YearMonth



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
    latlon_rounding_map: dict[Era5LandDim, float] = {
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
        )
        source_ds = round_coords(
            source_ds,
            rounding_map=latlon_rounding_map,
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
        Land files and representative values.
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
        be attempted. If set to UnmaskedNullsProcessing.FILL_NEAREST, the null
        values in the unmasked areas will be filled using nearest-neighbor
        values from the same variable in the same time step. No other options
        are currently supported. This parameter is ignored if `mask_file` is not
        provided. It *is* used if `if_unmasked_nulls` is set to
        UnmaskedNullsHandling.IGNORE, since that option only controls whether to
        log and report the presence of unmasked nulls to the user.
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
            mask_file=mask_file,
            if_masked_values=if_masked_values,
            if_unmasked_nulls=if_unmasked_nulls,
            unmasked_nulls_processing=unmasked_nulls_processing,
            null_value_file=_null_value_file,
        )
###END def convert_monthly_era5_files
