"""Code to request ERA5 Land GRIB files per month for selected variables."""

# %%
# Imports
# %%
import logging

import ecmwf.datastores as ecmwfds

import era5land_to_datm as etd
from era5land_to_datm.requests import (
    EcmwfDatastoreRequest,
    EcmfwRequestError,
    ERA5LAND_FOR_DATM_DATASET_ID,
    create_era5land_request,
    send_ecmwf_datastore_request,
)
from era5land_to_datm.regions import (
    NORWAY_RECT_0125_REGION,
)
from era5land_to_datm.types import (
    YearMonth,
)
from era5land_to_datm.variables import (
    VarSet,
    era5_datm_vars,
)

import logging_common

# %%
# Initialize logging
# %%
logger: logging.Logger = logging_common.initialize_logger()

# %% 
# Set the ranges of months and years to download
# %%
years_months: list[YearMonth] = [
    YearMonth(year=year, month=month) for year, month in [
        (2012, 10),
        (2011, 2),
        (2011, 1),
        (2010, 2),
        (2012, 9),
        (2010, 10),
        (2012, 1),
        (2013, 2),
        (2010, 6),
        (2010, 4),
        (2011, 8),
        (2012, 3),
        (2012, 2),
        (2012, 8),
        (2011, 4)
    ]
]

print(f'Requesting data for {len(years_months)} year-months:')
print(',\n'.join(f'{ym.year:04d}-{ym.month:02d}' for ym in years_months))

# %%
# Select variables to download. If not otherwise specified, use the defaults,
# which are the variables needed for DATM.
# %%
download_vars: VarSet = era5_datm_vars

# %%
# Get the request dictionary to use for downloading the data
# %%
requests: dict[YearMonth, EcmwfDatastoreRequest] = create_era5land_request(
    dataset_id=ERA5LAND_FOR_DATM_DATASET_ID,
    years_months=years_months,
    variables=download_vars,
    area=NORWAY_RECT_0125_REGION,
)

# %%
# Send the requests and store the Remote instances that are created
# %%
try:
    remotes: dict[VarSet, dict[YearMonth, ecmwfds.Remote]] =  (
        send_ecmwf_datastore_request(requests.values())
    )
except EcmfwRequestError as e:
    logger.error(f"Error sending ECMWF datastore request: {e}")
    request_error = e
    raise
