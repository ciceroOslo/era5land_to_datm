"""Convert ERA5-Land monthly data to DATM7 format for multiple years and months.

This script converts a series of ERA5 Land GRIB files with data for one month
each. Note that the output is also in the form of monthly files, not full-year
files as is used in the standard data set for CRUJRA in CESM.

The script converts the data to the threestream format used by the CRUJRA mode
of DATM7, with three different sets of files per month. It does not use the more
recent ERA5 mode, which did not seem to be stable or consistently documented at
the time of writing (January 2026).
"""
import logging
from pathlib import Path

from era5land_to_datm.convert_files import convert_monthly_era5_files


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Convert a series of monthly ERA5 Land GRIB files to DATM7 '
            'threestream netCDF files.'
        ),
    )

    parser.add_argument(
        'source_files',
        type=Path,
        nargs='+',
        help='Paths to the ERA5 Land GRIB files to convert.',
    )
    parser.add_argument(
        '--output-files',
        type=str,
        nargs=1,
        required=True,
        help=(
            'Output file path pattern, with {year}, {month} and {stream} '
            'placeholders for the year, month and stream of each file. The '
            'string pattern will be processed using the `str.format()` Python '
            'method.'
        ),
    )
    parser.add_argument(
        '--start-year-month',
        type=int,
        nargs=2,
        required=True,
        help=(
            'Start year and month (inclusive) for the conversion. Must be '
            'provided as two integers that are valid year and month values, '
            'separated by a space. For example: 2000 01.\n'
            'Note that the script will attempt to read the files for one month '
            'before the specified start month, using the pattern provided in '
            '--output-files. If present, it will be usd to compute cumulative '
            'variables that need the last value of the previous month. If not '
            'present, those variables will have missing values in the first '
            'time step. Also note that if a file with the expected pattern for '
            'the previous month is present but is not in the right format or '
            'contains invalid data, the script may fail or produce incorrect '
            'results.'
        ),
    )
    parser.add_argument(
        '--end-year-month',
        type=int,
        nargs=2,
        required=True,
        help=(
            'End year and month (inclusive) for the conversion. in the same '
            'format as start-year-month. All monthly files from '
            '--start-year-month through --end-year-month must be present in '
            'the specified source files. Note that the script will attempt to '
            'read the files for one month after the specified end month, using '
            'the pattern provided in --output-files. If present, it will be '
            'usd to compute cumulative variables that need the first value of '
            'the next month. If not present, those variables will have missing '
            'values in the last time step. Also note that if a file with the '
            'expected pattern for the next month is present but is not in the '
            'right format or contains invalid data, the script may fail or '
            'produce incorrect results.'
        ),
    )
    parser.add_argument(
        '--log-level',
        type=str,
        default='INFO',
        help=(
            'Logging level. One of DEBUG, INFO, WARNING, ERROR, CRITICAL. '
            'Default is INFO.'
        ),
    )
    args = parser.parse_args()
    logging.basicConfig(level=args.log_level.upper())
    logger = logging.getLogger(__name__)
    logger.debug(f'Parsed arguments: {args}')
    logger.debug(
        f'Calling `convert_monthly_era5_files` with parameters:'
        f'\n  source_files: {args.source_files}'
        f'\n  output_files: {args.output_files}'
        f'\n  start_year_month: {args.start_year_month}'
        f'\n  end_year_month: {args.end_year_month}'
    )
    convert_monthly_era5_files(
        source_files=args.source_files,
        output_files=args.output_files[0],
        start_year_month=args.start_year_month,
        end_year_month=args.end_year_month,
    )
