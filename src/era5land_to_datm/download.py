"""Code to make requests and download ERA5 Land variables and GRIB files.

Classes
-------
EcmfwDatastoreRequest
    Pydantic model representing a request dictionary to be submited through
    `ecmwf.datastores.Client.submit`.

Attributes
----------
ERA5LAND_FOR_DATM_DATASET_ID : EcmwfDatasetId
    The dataset ID for the ERA5-Land dataset to be used for conversion to DATM
    forcing data.

Functions
---------
make_era5land_request
    Create an EcmfwDatastoreRequest for downloading ERA5-Land data for the
    specified variables and time ranges.
send_ecmwf_datastore_request
    Send the given EcmfwDatastoreRequest to the ECMWF data store and return
    a mapping from YearMonth to the corresponding Remote instances.
"""
import typing as tp

import pydantic
import ecmwf.datastores as ecmwfds

from .types import (
    EcmwfDatasetId,
    Era5LandVar,
    YearMonth,
)
