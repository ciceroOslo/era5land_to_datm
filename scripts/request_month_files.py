"""Code to request ERA5 Land GRIB files per month for selected variables."""

# %%
# Imports
# %%

import ecmwf.datastores as ecmwfds

import era5land_to_datm as etd
from era5land_to_datm.download import (
    EcmwfDatastoreRequest,
    ERA5LAND_FOR_DATM_DATASET_ID,
    make_era5land_request,
    send_ecmwf_datastore_request,
)
from era5land_to_datm.types import (
    EcmwfDatasetId,
    YearMonth,
)
from era5land_to_datm.variables import (
    Era5LandVarMapping,
    Era5LandVar,
    VarSet,
    era5land_grib_varnames,
    era5land_request_varnames,
    era5_datm_vars,
)

# %% 
# Set the ranges of months and years to download
# %%
years_months: list[YearMonth] = [
    (_year, _month) for _year in range(2000, 2024 + 1)
    for _month in range(1, 12 + 1)
]

# %%
# Select variables to download. If not otherwise specified, use the defaults,
# which are the variables needed for DATM.
# %%
download_vars: VarSet = era5_datm_vars

# %%
# Get the request dictionary to use for downloading the data
# %%
request: EcmwfDatastoreRequest = make_era5land_request(
    dataset_id=ERA5LAND_FOR_DATM_DATASET_ID,
    years_months=years_months,
    variables=download_vars,
)

# %%
# Send the requests and store the Remote instances that are created
# %%
remotes: dict[YearMonth, ecmwfds.Remote] =  send_ecmwf_datastore_request(request)
