"""Code to download ERA5 Land GRIB files from completed requests."""

# %%
# Imports
# %%
import logging

import ecmwf.datastores as ecmwfds

from era5land_to_datm.download import (
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

import logging_common

# %%
# Initialize logging
# %%
logger: logging.Logger = logging_common.initialize_logger()

# %%
# Get the Remote instances for the previously sent requests
# %%
remotes: dict[str, CachedRemote] = get_remotes()
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
# Retrieve the available files for the completed requests
# %%
download_result: DownloadFilesResult = (
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
_year_month_dict: dict[YearMonth, ecmwfds.Remote]
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
