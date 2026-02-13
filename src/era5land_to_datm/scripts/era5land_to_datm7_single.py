#!/usr/bin/env python3
"""Script to convert a single ERA5 Land GRIB file to netCDF suitable for DATM7.

The script converts a single main ERA5 Land GRIB file, with optional previous
and next files to ensure proper handling of cumulative variables at the edges of
the time range.

This script converts the ERA5 Land data files to the threestream format used by
the CRUJRA mode of DATM7. It does not use the more recent ERA5 mode, which did
not seem to be stable or consistently documented at the time of writing
(January 2026).
"""
import logging
from pathlib import Path

from era5land_to_datm.convert_files import convert_era5_file
from era5land_to_datm.datm_streams import Datm7Stream
from era5land_to_datm.logger_registry import (
    register_logger,
    set_logging_level,
)


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
    parser.add_argument(
        '--keep-first-last-dates',
        action='store_true',
        help=(
            'Keep the first and last dates from the source ERA5 Land file in '
            'the converted DATM files. Normally, the one-month ERA5 Land files '
            'contain the last date of the previous month and the first '
            'timestep (midnight) of the first day in the next month, but with '
            'null values. These are normally discarded in the conversion in '
            'order to ensure that the output files only contain data for a '
            'single calendar month, and don\'t introduce overlaps in the file '
            'streams. Set this flag to keep these dates in the output.'
        ),
    )
    parser.add_argument(
        '--log-level',
        type=str,
        default='INFO',
        help=(
            'Logging level for logging done by this script and its submodules '
            '(not third-party libraries). One of DEBUG, INFO, WARNING, ERROR, '
            'CRITICAL. Default is INFO.'
        ),
    )
    parser.add_argument(
        '--global-log-level',
        type=str,
        default='WARNING',
        help=(
            'Global logging level for third-party packages. One of DEBUG, '
            'INFO, WARNING, ERROR, CRITICAL. Default is WARNING.'
        ),
    )
    args = parser.parse_args()
    logging.basicConfig(level=args.global_log_level.upper())
    set_logging_level(getattr(logging, args.log_level.upper()))
    logger = logging.getLogger(__name__)
    register_logger(logger)
    logger.debug(f'Parsed arguments: {args}')
    logger.debug(
        f'Calling `convert_era5_file` with parameters:'
        f'\n  source_file: {args.source_file}'
        f'\n  next_source_file: {args.next_source_file}'
        f'\n  previous_source_file: {args.previous_source_file}'
        f'\n  output_streams: {args.output_streams}'
        f'\n  keep_first_last_dates: {args.keep_first_last_dates}'
        f'\n  eager: {not args.lazy}'
    )
    convert_era5_file(
        source_file=args.source_file,
        next_source_file=args.next_source_file,
        previous_source_file=args.previous_source_file,
        output_streams=args.output_streams,
        keep_first_last_dates=args.keep_first_last_dates,
        eager=not args.lazy,
    )
