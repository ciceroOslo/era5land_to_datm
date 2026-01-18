"""Script to convert ERA5 Land GRIB files to netCDF files suitable for DATM7.

This script converts the ERA5 Land data files to the threestream format used by
the CRUJRA mode of DATM7. It does not use the more recent ERA5 mode, which did
not seem to be stable or consistently documented at the time of writing
(January 2026).

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
from collections.abc import (
    Callable,
    Mapping,
    Sequence
)
import enum
import logging
from pathlib import Path
import typing as tp

from korsbakken_python_utils.containers.dataobject import UniformTypeDataObject
import numpy as np
import xarray as xr

from era5land_to_datm.convert_data import make_datm_ds
from era5land_to_datm.datm_streams import Datm7Stream
from era5land_to_datm.file_io import (
    open_era5land_grib,
    write_datm_nc,
)



logger: logging.Logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)



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
        eager: bool = True,
        disable_dask: bool = False,
        logging_level: int = logging.INFO,
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
    root_logger: logging.Logger = logging.getLogger()
    if root_logger.level != logging_level:
        root_logger.setLevel(logging_level)
    global logger
    logger.setLevel(logging_level)
    if not disable_dask:
        from korsbakken_python_utils.dask_utils import (
            get_registered_progress_bars,
            set_global_progress_bar,
        )
        if len(get_registered_progress_bars()) == 0:
            set_global_progress_bar()
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
        f'  Main file: {source_file}\n'
        f'  Next timestep file: {next_source_file}\n'
        f'  Previous timestep file: {previous_source_file}\n'
    )
    source_ds: xr.Dataset = open_era5land_grib(
        file=source_file,
        next_file=next_source_file,
        previous_file=previous_source_file,
        chunks='auto' if not disable_dask else None,
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
        write_datm_nc(
            _target_ds,
            output_file=output_files_mapping[_target_stream],
            stream=_target_stream,
        )
        logger.info(f'    Finished writing {_target_stream.value}.')
    logger.info('Finished converting ERA5 Land GRIB file.')
###END def convert_era5_file


class _YearMonthStreamFilepathFunction(tp.Protocol):
    """Protocol for functions to make file paths from year, month and stream."""
    def __call__(self, year: int, month: int, stream: Datm7Stream) -> Path: ...
###END class _YearMonthStreamFilenameFunction

def convert_monthly_era5_files(
        *,
        source_files: Sequence[Path|str] | _YearMonthStreamFilepathFunction | str,
        output_files: Sequence[Path|str] | _YearMonthStreamFilepathFunction | str,
        start_year_month: tp.Tuple[int, int] | None = None,
        end_year_month: tp.Tuple[int, int] | None = None,
        next_source_file: Path|str|None = None,
        previous_source_file: Path|str|None = None,
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
    output_files : Sequence[Path | str] | Callable | str
        Output files, in the same format as `source_files`. If a sequence of
        string or Path instances is given, it must have the same length as
        `source_files`.
    start_year_month : Tuple[int, int]
        Tuple of (year, month) for the first month to convert, inclusive. If
        `source_files` and `output_files` are sequences, this parameter is
        ignored and can be left unspecified (None). If `source_files` and
        `output_files` are functions or format strings, this parameter must be
        specified, or a ValueError will be raised.
    end_year_month : Tuple[int, int] | None
        Tuple of (year, month) for the last month to convert, inclusive. If
        `source_files` and `output_files` are sequences, this parameter is
        ignored and can be left unspecified (None). If `source_files` and
        `output_files` are functions or format strings, this parameter must be
        specified, or a ValueError will be raised.
    next_source_file : Path | str | None, optional
        Path to the ERA5 Land GRIB file for the month after `end_year_month`.
        This is optional, but some variables may have missing data in the last
        time step in the last output file if this is not provided.
    previous_source_file : Path | str | None, optional
        Path to the ERA5 Land GRIB file for the month before `start_year_month`.
        This is optional, but some variables may have missing data in the first
        time step in the first output file if this is not provided.
    """
    if isinstance(source_files, str):
        source_files_func: Callable|None = source_files.format
    if callable(source_files):
        source_files_func = source_files
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
            source_files_func(year=_year, month=_month, stream=_stream)
            for _year in range(start_year_month[0], end_year_month[0] + 1)
            for _month in range(start_year_month[1], end_year_month[1] + 1)
            for _stream in Datm7Stream
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
    logger.debug(f'Using source files: {use_source_files}')
    if isinstance(output_files, str):
        output_files_func: Callable = output_files.format
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
        output_files_func = output_files
        use_output_files: list[Path] = [
            output_files_func(year=_year, month=_month, stream=_stream)
            for _year in range(start_year_month[0], end_year_month[0] + 1)
            for _month in range(start_year_month[1], end_year_month[1] + 1)
            for _stream in Datm7Stream
        ]
    elif isinstance(output_files, Sequence):
        use_output_files = [Path(_f) for _f in output_files]
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
    existing_output_files: list[Path] = [
        _f
        for _f in use_output_files
        if _f.exists()
    ]
    if len(existing_output_files) > 0:
        raise FileExistsError(
            'The following output files already exist, please move or delete '
            'them, or choose different names:\n' + ',\n'.join(
                str(_f) for _f in existing_output_files
            )
        )
    logger.debug(f'Using output files: {use_output_files}')
    if previous_source_file is None and source_files_func is not None:



if __name__ == '__main__':
    import argparse
 
    parser = argparse.ArgumentParser(
        description=(
            'Convert an ERA5 Land GRIB file to DATM7 threestream netCDF files.'
        ),
    )
    parser.add_argument(
        'source_file',
        type=Path,
        help='Path to the ERA5 Land GRIB file for the current time period.',
    )
    parser.add_argument(
        '--next-source-file',
        type=Path,
        default=None,
        help=(
            'Path to the ERA5 Land GRIB file for the subsequent time period. '
            'Needed to properly compute some variables in the last time step.'
        ),
    )
    parser.add_argument(
        '--previous-source-file',
        type=Path,
        default=None,
        help=(
            'Path to the ERA5 Land GRIB file for the previous time period. '
            'Needed to properly compute some variables in the first time step.'
        ),
    )
    parser.add_argument(
        '--output-streams',
        type=str,
        nargs='+',
        default=[_s.value for _s in Datm7Stream],
        help=(
            'List of DATM7 streams to generate. By default, all streams are '
            'generated.'
        ),
    )
    parser.add_argument(
        '--lazy',
        action='store_true',
        help=(
            'Lazy-load data with dask to reduce memory usage. May slow down '
            'the processing considerably due to some dating being read '
            'multiple times, and some processing steps having poor '
            'dask support.'
        ),
    )
    args = parser.parse_args()
    convert_era5_file(
        source_file=args.source_file,
        next_source_file=args.next_source_file,
        previous_source_file=args.previous_source_file,
        output_streams=args.output_streams,
        eager=not args.lazy,
    )
