"""Code to download ERA5 Land GRIB files from completed requests."""

# %%
# Imports
# %%
import logging
import os
import re
import time

import ecmwf.datastores as ecmwfds

from era5land_to_datm.download import (
    DownloadedFile,
    DownloadFilesResult,
    get_remotes,
    get_remote_yearmonth,
    get_remote_varset,
    make_grib_filename,
    remotes_dict_by_vars_and_yearmonth,
    retrieve_available_files,
)
from era5land_to_datm.remote import (
    CachedRemote,
)
from era5land_to_datm.types import (
    YearMonth,
)
from era5land_to_datm.variables import (
    VarSet,
)

# %%
# Dynamically import logging_common.py from the scripts/ directory (i.e., the
# same directory as this script file). We need to import the file with a path
# set using `__file__`, because the working directory may not be the same as the
# script directory when this script is run, and it will usually also not be in
# the Python module search path. The script is also not inside a package, so we
# cannote use `__name__` or a relative import (`__name__` will be `__main__`).
# %%
import sys
import pathlib
script_dir: pathlib.Path = pathlib.Path(__file__).parent.resolve()
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))
import logging_common
# logging_common = importlib.import_module('logging_common')

# %%
# Initialize logging
# %%
logger: logging.Logger = logging_common.initialize_logger()

# %%
# Get the Remote instances for the previously sent requests
# %%
remotes: dict[str, CachedRemote] = get_remotes(limit=300)
# %%
remotes_dict: dict[VarSet, dict[YearMonth, CachedRemote]] = (
    remotes_dict_by_vars_and_yearmonth(remotes.values())
)

# %%
# Check the status of the requests
# %%
statuses: dict[VarSet, dict[YearMonth, str]] = {
    _var_set: {
        _year_month: _remote.status
        for _year_month, _remote in _year_month_remotes.items()
    }
    for _var_set, _year_month_remotes in remotes_dict.items()
}
print('Request statuses:')

for _var_set, _year_month_statuses in statuses.items():
    print(f'VarSet: {_var_set}')
    for (_year, _month), _status in _year_month_statuses.items():
        print(f'  {_year:04d}-{_month:02d}: {_status}')

# %%
# Set the download directory
# %%
download_dir = pathlib.Path(
    '/cluster/shared/noresm/inputdata/cicero_mods'
    '/ERA5-Land_NorSink/original_GRIB_NorwayRect/'
)
# %%
# Change to the download directory, without attempting to create it.
# %%
os.chdir(download_dir)
print(f'Changed working directory to {pathlib.Path.cwd()}')

# %%
# Retrieve the available files for the completed requests
# %%
download_result: DownloadFilesResult[CachedRemote] = (
    retrieve_available_files(remotes.values())
)

# %%
# Print summary of downloaded files
# %%
print('Download summary:\n')
print('Downloaded the following files:')
for _remote, _filepath in download_result.downloaded_files:
    print(f'  {_filepath}')

print(
    'The following files were already present locally  and were not downloaded '
    'again:'
)
for _remote, _filepath in download_result.existing_files:
    print(f'  {_filepath}')

print(
    'The following files were already present locally  and were not downloaded '
    'again:'
)
for _remote, _filepath in download_result.existing_files:
    print(f'  {_filepath}')

print('\nThe following requests had no available files to download:')
_var_set: VarSet
_year_month: YearMonth
_year_month_dict: dict[YearMonth, CachedRemote]
for _var_set, _year_month_dict in remotes_dict_by_vars_and_yearmonth(
    download_result.no_files_remotes
).items():
    print(f'  VarSet: {_var_set}')
    for (_year, _month) in _year_month_dict.keys():
        print(f'    {_year:04d}-{_month:02d}')

print('\nThe following downloads failed:')
failed_downloads: dict[VarSet, dict[YearMonth, Exception]] = {}
for _remote, _error in download_result.failed_remotes:
    _var_set = get_remote_varset(_remote)
    _year_month = get_remote_yearmonth(_remote)
    failed_downloads.setdefault(_var_set, {})[_year_month] = _error
for _var_set, _year_month_errors in failed_downloads.items():
    print(f'  VarSet: {_var_set}')
    for (_year, _month), _error in _year_month_errors.items():
        print(f'    {_year:04d}-{_month:02d}: {repr(_error)}')

# %%
# Check and write all files that are available locally. **NB!** This code
# searches for files using a regex. It does not change when
# `make_grib_filename` changes, so you must change it manually if you change
# the code in `make_grib_filename` or provide a custom file name pattern.
# %%
# Create a regex pattern to extract year and month from filenames, then make a
# dict with YearMonth keys and a tuple of file size and file name as values.
year_month_filename_pattern: re.Pattern = re.compile(
    r'^.*_(?P<year>\d{4})_(?P<month>\d{2})\.grib$'
)
existing_files_info: dict[YearMonth, tuple[int, str]] = {}
for _path in download_dir.glob('*.grib'):
    _match = year_month_filename_pattern.match(_path.name)
    if _match is not None:
        _year: int = int(_match.group('year'))
        _month: int = int(_match.group('month'))
        existing_files_info[
            YearMonth(_year, _month)
        ] = (_path.stat().st_size, _path.name)
    else:
        logger.warning(
            f'File {_path} does not match expected filename pattern; '
            'skipping...'
        )
# %%
existing_file_sizes: dict[YearMonth, int] = {
    _year_month: _file_info[0]
    for _year_month, _file_info in sorted(existing_files_info.items())
}

# %%
# Delete the remotes for the files that were successfully downloaded.
#
# The first cell below deletes the downloaded instances, the second cell deletes
# the ones that were already present locally and not downloaded again.
#
# Follow the same pattern manually if you wish to delete any remotes that had
# errors associated with them.
#
# **NB!** Since this may not be wanted, the variable
# `clean_up_downloaded_remotes` must be set to True, or the deletions will not
# be carried out.
# %%
clean_up_downloaded_remotes: bool = False
deleted_remotes: list[tuple[YearMonth, DownloadedFile[CachedRemote]]] = []
logger.info('Deleting downloaded remotes...')
if clean_up_downloaded_remotes:
    for _remote, _filepath in download_result.downloaded_files:
        _year_month = get_remote_yearmonth(_remote)
        logger.info(
            f'Deleting Remote with id {_remote.request_id} for '
            f'{_year_month.year:04d}-{_year_month.month:02d} after '
            f'downloading file {_filepath}...'
        )
        _remote.delete()
        deleted_remotes.append(
            (_year_month, DownloadedFile(_remote, _filepath))
        )
        time.sleep(0.5)
# %%
logger.info('Deleting remotes for files already present locally...')
if clean_up_downloaded_remotes:
    for _remote, _filepath in download_result.existing_files:
        _year_month = get_remote_yearmonth(_remote)
        logger.info(
            f'Deleting Remote with id {_remote.request_id} for '
            f'{_year_month.year:04d}-{_year_month.month:02d} after '
            f'not downloading file {_filepath} (already present locally)...'
        )
        _remote.delete()
        deleted_remotes.append(
            (_year_month, DownloadedFile(_remote, _filepath))
        )
        time.sleep(0.5)
