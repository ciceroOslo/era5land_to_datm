"""Code to download ERA5 Land GRIB files from completed requests."""

# %%
# Imports
# %%
import logging

import ecmwf.datastores as ecmwfds

import era5land_to_datm as etd
from era5land_to_datm.download import (
    DownloadFilesResult,
    get_remotes,
    get_remote_yearmonth,
    get_remote_varset,
    make_grib_filename,
    remotes_dict_by_vars_and_yearmonth,
    retrieve_available_files,
)
from era5land_to_datm.types import (
    VarSet,
    YearMonth,
)

import logging_common

# %%
# Initialize logging
# %%
logger: logging.Logger = logging_common.initialize_logger()

# %%
# Get the Remote instances for the previously sent requests
# %%
remotes: dict[str, ecmwfds.Remote] = etd.download.get_remotes()
remotes_dict: dict[VarSet, dict[YearMonth, ecmwfds.Remote]] = (
    remotes_dict_by_vars_and_yearmonth(remotes)
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
for _remote, _filepath in download_result.downloaded_files.items():
    print('Downloaded the following files:')
    print(f'  {_filepath}')

print('\nThe following requests had no available files to download:')
for _var_set, _year_month in remotes_dict_by_vars_and_yearmonth(
    download_result.no_files_remotes
).items():
    print(f'  VarSet: {_var_set}')
    for (_year, _month) in _year_month.keys():
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
