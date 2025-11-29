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
import datetime
import typing as tp

import pydantic
import ecmwf.datastores as ecmwfds

from .types import (
    EcmwfDatasetId,
    YearMonth,
)
from .variables import (
    Era5LandVar,
    VarSet,
)



ERA5LAND_FOR_DATM_DATASET_ID: tp.Final[EcmwfDatasetId] = (
    EcmwfDatasetId.ERA5LANDHRLY
)

type DayInMonthInt = tp.Annotated[int, pydantic.Field(ge=1, le=31)]

def _validate_full_hour_time(
    time_value: datetime.time,
) -> datetime.time:
    """Validate that the given time represents a full hour (i.e., minutes
    and seconds are 0).

    Parameters
    ----------
    time_value : datetime.time
        The time to validate.

    Returns
    -------
    datetime.time
        The validated time.

    Raises
    ------
    ValueError
        If the given time does not represent a full hour.
    """
    if (
            time_value.minute != 0 
            or time_value.second != 0
            or time_value.microsecond != 0
    ):
        raise ValueError(
            'Time must represent a full hour (i.e., minutes, seconds, and '
            'microseconds are 0).'
        )
    return time_value
###END def _validate_full_hour_time

type FullHourTime = tp.Annotated[
    datetime.time,
    pydantic.Field(
        description=(
            'A time representing a full hour (i.e., minutes, seconds, and microseconds are 0).'
        )
    ),
    pydantic.AfterValidator(_validate_full_hour_time),
]

type LatitudeFloat = tp.Annotated[
    float,
    pydantic.Field(ge=-90.0, le=90.0),
]
type LongitudeFloat = tp.Annotated[
    float,
    pydantic.Field(ge=-360.0, le=360.0),
]


class EcmwfAreaTuple(tp.NamedTuple):
    """A namedtuple representing a geographical area for ECMWF data requests.

    Attributes
    ----------
    north : float
        The northern bound in degrees (-90.0 to 90.0).
    west : float
        The western bound in degrees (-360.0 to 360.0).
    south : float
        The southern bound in degrees (-90.0 to 90.0).
    east : float
        The eastern bound in degrees (-360.0 to 360.0).
    """
    north: LatitudeFloat
    west: LongitudeFloat
    south: LatitudeFloat
    east: LongitudeFloat
###END class EcmwfAreaTuple

class EcmwfDatastoreRequest(pydantic.BaseModel):
    """Pydantic model representing a request dictionary to be submited through
    `ecmwf.datastores.Client.submit`.

    The fields in this model have the data type that is logically appropriate
    for each field. Most will be converted to the appropriate string format
    when the request is sent. The dict representation to be sent to ECMWF can be
    obtained through the `.to_request_dict()` method.

    Note that the dataset ID is included as a separate field in this model, but
    is not included in the request dictionary itself, as it is specified
    separately when submitting the request.

    Fields
    ------
    dataset_id : EcmwfDatasetId
        The dataset ID to request data from.
    variable : VarSet
        Frozenset with the set of ERA5-Land variables to request.
    year : int
        The year to request data for. Will be converted to a YYYY string when
        the request is sent.
    month : int
        The month to request data for, as an integer from 1 to 12. Will be
        converted to a string when the request is sent.
    day : tuple[int, ...]
        Tuple of integers representing the days of the month to request data
        for, as integers from 1 to 31. Note that days from 29 through 31 may be
        included regardless of whether they exist in the given month.
    time : tuple[datetime.time, ...]
        Tuple of datetime.time instances representing the times of day to
        request data for. Note that for most datasets, only full hours are
        valid. Will be converted to 'HH:MM' strings when the request is sent.
    data_format : Literal['grib', 'netcdf']
        The data format to request. Note that ECMWF considers 'netcdf' to be
        experimental, and may lead to larger files and longer processing times.
        Optional, default is 'grib'.
    download_format : Literal['unarchived', 'zip']
        The download format to request. Optional, default is 'unarchived'.
    area : tuple[float, float, float, float] or None
        The geographical area to request data for, as a tuple of four floats
        representing the North, West, South, and East bounds in degrees. The
        latitudes (elements 0 and 2) must be between -90.0 and 90.0, and the
        longitudes (elements 1 and 3) must be between -360.0 and 360.0. South
        (element 2) must be less than North (element 0), and East (element 3)
        must be greater than West (element 1). Optional, set to None to request
        global data. Default is None.

    Methods
    -------
    to_request_dict
        Convert the EcmwfDatastoreRequest instance to a dictionary suitable
        for submission to `ecmwf.datastores.Client.submit`. In the resulting
        dict, most of the numeric fields will be converted to strings as
        expected by ECMWF.
    """

    dataset_id: EcmwfDatasetId
    variable: VarSet
    year: int
    month: int = pydantic.Field(ge=1, le=12)
    day: tuple[DayInMonthInt, ...]
    time: tuple[FullHourTime, ...]
    data_format: tp.Literal['grib', 'netcdf'] = 'grib'
    download_format: tp.Literal['unarchived', 'zip'] = 'unarchived'
    area: EcmwfAreaTuple | None = None

    def to_request_dict(self) -> dict[str, tp.Any]:
        """Convert the EcmwfDatastoreRequest instance to a dictionary suitable
        for submission to `ecmwf.datastores.Client.submit`. In the resulting
        dict, most of the numeric fields will be converted to strings as
        expected by ECMWF.

        Note that the request dictionary does not include the dataset ID, which
        is specified separately when submitting the request.

        Returns
        -------
        dict[str, tp.Any]
            The request dictionary.
        """
        request_dict: dict[str, tp.Any] = {
            'variable': sorted(
                _var.request_varname
                for _var in self.variable
            ),
            'year': f'{self.year:04d}',
            'month': f'{self.month:02d}',
            'day': [f'{_day:02d}' for _day in self.day],
            'time': [f'{_time.hour:02d}:00' for _time in self.time],
            'data_format': self.data_format,
            'download_format': self.download_format,
        }
        if self.area is not None:
            request_dict['area'] = [
                f'{self.area.north}',
                f'{self.area.west}',
                f'{self.area.south}',
                f'{self.area.east}',
            ]
        return request_dict
    ###END def EcmwfDatastoreRequest.to_request_dict
