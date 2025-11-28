"""Code to download ERA5 Land GRIB files from completed requests."""

# %%
# Imports
# %%
import ecmwf.datastores as ecmwfds

import era5land_to_datm as etd
from era5land_to_datm.download import (
    DownloadFilesResult,
    make_grib_filename,
    get_remotes,
    retrieve_available_files,
)
from era5land_to_datm.types import (
    YearMonth,
)

# %%
# Get the Remote instances for the previously sent requests
# %%
remotes: dict[YearMonth, ecmwfds.Remote] = etd.download.get_remotes()

# %%
# Check the status of the requests
# %%
statuses: dict[YearMonth, str] = {
    _year_month: _remote.status for _year_month, _remote in remotes.items()
}
print('Request statuses:')
for (_year, _month), _status in statuses.items():
    print(f'  {_year:04d}-{_month:02d}: {_status}')

# %%
# Retrieve the available files for the completed requests
# %%
download_result: DownloadFilesResult = retrieve_available_files(remotes)

# %%
# Print summary of downloaded files
# %%
print('Download summary:\n')
print('  Downloaded following files for following years/months:')
for (_year, _month), _path in download_result.downloaded_files.values():
    print(f'    {_year:04d}-{_month:02d}: {str(_path)}')
print('\n  No files available for following years/months:')
