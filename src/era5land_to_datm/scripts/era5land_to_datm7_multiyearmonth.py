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
from era5land_to_datm.logger_registry import (
    register_logger,
    set_logging_level,
)
from era5land_to_datm.masking import (
    MaskedValuesHandling,
    UnmaskedNullsHandling,
    UnmaskedNullsProcessing,
)



def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Convert a series of monthly ERA5 Land GRIB files to DATM7 '
            'threestream netCDF files.'
        ),
        epilog=(
            '**NB!** If you will enter the processed files into CIME XML '
            'database files for use with CESM, you should ensure that the file '
            'names contain year and month in the format `YYYY-MM`. See the '
            'documentation of the `--output-files` argument.'
        ),
    )
    parser.add_argument(
        '--source-dir',
        type=Path,
        required=True,
        help=(
            'Directory containing the source ERA5 Land GRIB files. If not an '
            'absolute path, it will be resolved relative to the current '
            'working directory. It is recommended to use absolute paths to '
            'avoid possible issues with what directory Python sets as the '
            'current working directory when running the script.'
        ),
    )
    parser.add_argument(
        '--source-files',
        type=str,
        required=True,
        help=(
            'Pattern for the file names of the source ERA5 Land GRIB files. '
            'The pattern must include {year} and {month} placeholders for the '
            'year and month of each file. The string pattern will be processed '
            'using the `str.format()` Python method. **NB!** If the years '
            'and/or months in the file names use leading zeros (e.g., month 01 '
            'for January), the placeholders must include the appropriate '
            'formatting specifiers, e.g., {month:02d} for two-digit months '
            'with leading zero.'
        ),
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        required=True,
        help=(
            'Directory where the output DATM7 netCDF files will be saved. If '
            'not an absolute path, it will be resolved relative to the '
            'current working directory. It is recommended to use absolute '
            'paths to avoid possible issues with what directory Python sets '
            'as the current working directory when running the script.'
        ),
    )
    parser.add_argument(
        '--output-files',
        type=str,
        required=True,
        help=(
            'Output file path pattern, with {year}, {month} and {stream} '
            'placeholders for the year, month and stream of each file. The '
            'string pattern will be processed using the `str.format()` Python '
            'method. The same comment about formatting of year and month '
            'applies as for --source-files.\n'
            '**NB!** The XML database files used by CIME under CESM require '
            'that the file names contain year and month in the format YYYY-MM '
            'if you want to parametrize the file paths by year and month '
            'rather than manually listing every single file. In this case, the '
            'pattern must include the year and month placeholders as '
            '{year:04d}-{month:02d}`.'
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
            '--source-files. If present, it will be usd to compute cumulative '
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
            'the pattern provided in --source-files. If present, it will be '
            'usd to compute cumulative variables that need the first value of '
            'the next month. If not present, those variables will have missing '
            'values in the last time step. Also note that if a file with the '
            'expected pattern for the next month is present but is not in the '
            'right format or contains invalid data, the script may fail or '
            'produce incorrect results.'
        ),
    )
    parser.add_argument(
        '--round-lat-to',
        type=float,
        default=0.1,
        help=(
            'Round latitude values to the nearest multiple of this value. This '
            'can be useful to avoid issues with very small floating point '
            'differences in latitude values between the source ERA5 Land files '
            'and the expected latitude values in the DATM7 files. The '
            'coordinates in many of the ERA5 Land files are very slightly off '
            'from multiples of 0.1 degrees, and are therefore by default '
            'rounded to the nearest multiple of 0.1 degrees to avoid issues '
            'when comparing or merging with other data. To not round at all, '
            'set this to 0.'
        ),
    )
    parser.add_argument(
        '--round-lon-to',
        type=float,
        default=0.1,
        help=(
            'Round longitude values to the nearest multiple of this value. '
            'Same as --round-lat-to but for longitude values. To not round at '
            'all, set this to 0.'
        ),
    )
    parser.add_argument(
        '--mask-file',
        type=Path,
        help=(
            'Path to a netCDF file containing a mask file. The file should be '
            'produced using the `create_era5land_to_datm_mask_file.py` script '
            'and follow the dimension and variable naming from there. Use the '
            'options `--if-masked-values`, `--if-unmasked-nulls`, '
            '`--unmaked_nulls_processing` and `--null-value-files` to specify '
            'respectively how to report and whether to fail if non-null values '
            'are found in the masked areas or nulls are found in the unmasked '
            'areas, whether and how to attempt to fill unmasked null values, '
            'and to specify file name patterns for netCDF files in which to '
            'put boolean values that are True whereever unmasked null values '
            'were found.'
        ),
    )
    parser.add_argument(
        '--if-masked-values',
        type=MaskedValuesHandling,
        default=MaskedValuesHandling.RAISE.value,
        choices=[_s.value for _s in MaskedValuesHandling],
        help=(
            'How to handle non-null values found in the masked-out area '
            'defined by the mask file provided with --mask-file. Options are:\n'
            f'- {MaskedValuesHandling.RAISE.value}:\n'
            '    Raise an error if any non-null values are found in the masked '
            'area.\n'
            f'- {MaskedValuesHandling.WARN.value}:\n'
            '    Log a warning if any non-null values are found in the masked '
            'area, but do not raise an error.\n'
            f'- {MaskedValuesHandling.IGNORE.value}:\n'
            '    Ignore any non-null values found in the masked area and do not '
            'raise an error or log a warning.'
        ),
    )
    parser.add_argument(
        '--if-unmasked-nulls',
        type=UnmaskedNullsHandling,
        default=UnmaskedNullsHandling.WARN.value,
        choices=[_s.value for _s in UnmaskedNullsHandling],
        help=(
            'How to handle null values found in the unmasked area defined by '
            'the mask file provided with --mask-file. Options are:\n'
            f'- {UnmaskedNullsHandling.RAISE.value}:\n'
            '    Raise an error if any null values are found in the unmasked '
            'area.\n'
            f'- {UnmaskedNullsHandling.WARN.value}:\n'
            '    Log a warning if any null values are found in the unmasked '
            'area, but do not raise an error.\n'
            f'- {UnmaskedNullsHandling.IGNORE.value}:\n'
            '    Ignore any null values found in the unmasked area and do not '
            'raise an error or log a warning.'
        ),
    )
    parser.add_argument(
        '--unmasked-nulls-processing',
        type=UnmaskedNullsProcessing,
        default=UnmaskedNullsProcessing.LINEAR.value,
        choices=[_s.value for _s in UnmaskedNullsProcessing],
        help=(
            'Whether and how to attempt to fill null values found in the '
            'unmasked area defined by the mask file provided with --mask-file. '
            'Options are:\n'
            f'- {UnmaskedNullsProcessing.LINEAR.value}:\n'
            '    Attempt to fill null values in the unmasked area by filling '
            'them with the nearest non-null value in the same variable. This '
            'is the default option.\n'
            f'- {UnmaskedNullsProcessing.NONE.value}:\n'
            '    Leave null values in the unmasked area as they are.'
        ),
    )
    parser.add_argument(
        '--null-value-files',
        type=str,
        help=(
            'Pattern for file names to save boolean netCDF files in which to '
            'report where unmasked null values were found. The pattern follows '
            'the same formatting rules as --source-files. The output files '
            'will contain boolean variables with the same names as the ERA5 '
            'Land variables that were processed, and will have the value True '
            'wherever unmasked null values were found in the original data, '
            'and False elsewhere. If not provided, or if the pattern is an '
            'empty string, these files will not be saved. This option is '
            'independent of the --if-unmasked-nulls option, but if that '
            f'option is set to {UnmaskedNullsHandling.RAISE.value}, the script '
            'will raise an error and not save any more files once the first '
            'unmasked null value is encountered. To get a report of all '
            'unmasked null values, that option must be set to either warn or '
            'ignore.'
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

    source_dir: Path = args.source_dir.resolve()
    output_dir: Path = args.output_dir.resolve()
    source_files_pattern: str = str(source_dir / args.source_files)
    output_files_pattern: str = str(output_dir / args.output_files)
    logger.debug(
        f'Calling `convert_monthly_era5_files` with parameters:'
        f'\n  source_files: {source_files_pattern}'
        f'\n  output_files: {output_files_pattern}'
        f'\n  start_year_month: {args.start_year_month}'
        f'\n  end_year_month: {args.end_year_month}'
    )
    convert_monthly_era5_files(
        source_files=source_files_pattern,
        output_files=output_files_pattern,
        start_year_month=tuple(args.start_year_month),
        end_year_month=tuple(args.end_year_month),
        round_lat_to=args.round_lat_to,
        round_lon_to=args.round_lon_to,
        mask_file=args.mask_file,
        if_masked_values=args.if_masked_values,
        if_unmasked_nulls=args.if_unmasked_nulls,
        unmasked_nulls_processing=args.unmasked_nulls_processing,
        null_value_files=args.null_value_files,
    )

###END def main

if __name__ == '__main__':
    main()
