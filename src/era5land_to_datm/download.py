"""Code to make requests and download ERA5 Land variables and GRIB files.

Classes
-------
EcmfwDatastoreRequest
    Pydantic model representing a request dictionary to be submited through
    `ecmwf.datastores.Client.submit`.
EcmfwRequestError
    Exception raised when there is an error sending a request to the ECMWF data
    store in the function `send_ecmwf_datastore_request`. This exception may be
    raised after some requests have been successfully sent, in which case the
    `successful_requests` attribute contains a dictionary of the
    `EcmwfDatastoreRequest instances, and the `remotes` attribute a dictionary
    of the corresponding Remote instances, in the same format as the return
    value of `send_ecmwf_datastore_request`.

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
    a nested dictionary from VarSet and then YearMonth to the corresponding
    Remote instance.
"""
import datetime
import logging
import time
import typing as tp

import pydantic
import ecmwf.datastores as ecmwfds

from .regions import (
    EcmwfArea,
)
from .types import (
    EcmwfDatasetId,
    YearMonth,
)
from .variables import (
    Era5LandVar,
    VarSet,
)



default_logger: logging.Logger = logging.getLogger(__name__)

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
    area : EcmwfAreaTuple | None
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
    area: EcmwfArea | None = None

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
                self.area.north,
                self.area.west,
                self.area.south,
                self.area.east,
            ]
        return request_dict
    ###END def EcmwfDatastoreRequest.to_request_dict


def create_era5land_request(
    dataset_id: EcmwfDatasetId,
    years_months: tp.Iterable[YearMonth],
    variables: VarSet,
    *,
    days: tp.Iterable[DayInMonthInt] = range(1, 31+1),
    times: tp.Iterable[FullHourTime] = (
        datetime.time(hour=_h) for _h in range(0, 23+1)
    ),
    data_format: tp.Literal['grib', 'netcdf'] = 'grib',
    download_format: tp.Literal['unarchived', 'zip'] = 'unarchived',
    area: EcmwfArea | None = None,
) -> dict[YearMonth, EcmwfDatastoreRequest]:
    """Create EcmwfDatastoreRequest instances for downloading ERA5-Land data
    for the specified variables and time ranges.

    Parameters
    ----------
    dataset_id : EcmwfDatasetId
        The dataset ID to request data from.
    years_months : Iterable[YearMonth]
        Iterable of YearMonth instances representing the year and month
        combinations to create requests for.
    variables : VarSet
        Frozenset with the set of ERA5-Land variables to request.
    days : Iterable[int], optional
        Iterable of integers representing the days of the month to request
        data for, as integers from 1 to 31. Note that days from 29 through 31
        may be included regardless of whether they exist in the given month.
        Default is all days from 1 to 31.
    times : Iterable[datetime.time], optional
        Iterable of datetime.time instances representing the times of day to
        request data for. Note that for most datasets, only full hours are
        valid. Default is all full hours from 00:00 to 23:00.
    data_format : Literal['grib', 'netcdf'], optional
        The data format to request. Note that ECMWF considers 'netcdf' to be
        experimental, and may lead to larger files and longer processing times.
        Default is 'grib'.
    download_format : Literal['unarchived', 'zip'], optional
        The download format to request. Default is 'unarchived'.
    area : EcmwfAreaTuple | None, optional
        The geographical area to request data for, specified as a tuple of
        (north, west, south, east) coordinates. If None, the entire dataset
        area is requested. Default is None.

    Returns
    -------
    dict[YearMonth, EcmwfDatastoreRequest]
        Dictionary mapping YearMonth instances to EcmwfDatastoreRequest instances
        representing the requests to be made.
    """
    requests: dict[YearMonth, EcmwfDatastoreRequest] = {}
    for year, month in years_months:
        requests[YearMonth(year=year, month=month)] = EcmwfDatastoreRequest(
            dataset_id=dataset_id,
            variable=variables,
            year=year,
            month=month,
            day=tuple(days),
            time=tuple(times),
            data_format=data_format,
            download_format=download_format,
            area=area,
        )
    return requests
###END def create_era5land_request


class EcmfwRequestError(Exception):
    """Exception raised when there is an error sending a request to the ECMWF
    data store in the function `send_ecmwf_datastore_request`. This exception
    may be raised after some requests have been successfully sent, in which case
    the `successful_requests` attribute contains a dictionary of the
    corresponding remote instances, in the same format as the return value of
    `send_ecmwf_datastore_request`.
    """

    def __init__(
        self,
        message: str,
        successful_requests: dict[VarSet, dict[YearMonth, EcmwfDatastoreRequest]],
        remotes: dict[VarSet, dict[YearMonth, ecmwfds.Remote]],
    ) -> None:
        """Initialize the EcmfwRequestError instance.

        Parameters
        ----------
        message : str
            The error message.
        successful_requests : dict[VarSet, dict[YearMonth, EcmwfDatastoreRequest]]
            Dictionary of successfully sent requests, mapping VarSet to a
            dictionary mapping YearMonth to the corresponding
            EcmwfDatastoreRequest instance.
        remotes : dict[VarSet, dict[YearMonth, ecmwfds.Remote]]
            Dictionary of Remote instances corresponding to the successfully
            sent requests, mapping VarSet to a dictionary mapping YearMonth to
            the corresponding Remote instance.
        """
        super().__init__(message)
        self.successful_requests = successful_requests
        self.remotes = remotes
    ###END def EcmfwRequestError.__init__

###END class EcmfwRequestError


def send_ecmwf_datastore_request(
    requests: tp.Iterable[EcmwfDatastoreRequest],
    *,
    delay: float = 3.0,
    logger: logging.Logger = default_logger,
) -> dict[VarSet, dict[YearMonth, ecmwfds.Remote]]:
    """Send the given EcmwfDatastoreRequest to the ECMWF data store and return
    a nested dictionary from VarSet and then YearMonth to the corresponding
    Remote instance.

    Parameters
    ----------
    requests : Iterable[EcmwfDatastoreRequest]
        Iterable of EcmwfDatastoreRequest instances to send.
    delay : float, optional
        Delay in seconds between sending each request. Default is 3.0 seconds.
    logger : logging.Logger, optional
        Logger to use for logging messages. Default is a module-level logger.

    Returns
    -------
    dict[VarSet, dict[YearMonth, ecmwfds.Remote]]
        Nested dictionary mapping VarSet to a dictionary mapping YearMonth to
        the corresponding Remote instance.

    Raises
    ------
    EcmfwRequestError
        If there is an error sending any of the requests. The
        `successful_requests` attribute contains a dictionary of the
        successfully sent requests, and the `remotes` attribute contains a
        dictionary of the corresponding Remote instances.
    """
    if not isinstance(delay, (int, float)):
        raise TypeError('`delay` must be a number representing seconds.')
    if delay < 0.0:
        raise ValueError('`delay` must be non-negative.')
    client: ecmwfds.Client = ecmwfds.Client()
    remotes: dict[VarSet, dict[YearMonth, ecmwfds.Remote]] = {}
    successful_requests: dict[VarSet, dict[YearMonth, EcmwfDatastoreRequest]] = {}
    try:
        for (_request_num, _request) in enumerate(requests):
            if _request_num > 0 and delay > 0.0:
                time.sleep(delay)
            logger.info(
                f'Submitting request for year={_request.year}, '
                f'month={_request.month:02d}...'
            )
            remote: ecmwfds.Remote = client.submit(
                collection_id=_request.dataset_id.collection_id,
                request=_request.to_request_dict(),
            )
            logger.info(
                f'Request submitted with ID {remote.request_id}, status '
                f'{remote.status}.')
            var_set: VarSet = _request.variable
            year_month: YearMonth = (
                YearMonth(year=_request.year, month=_request.month)
            )
            remotes.setdefault(var_set, dict())[year_month] = remote
            successful_requests.setdefault(var_set, dict())[year_month] = _request
    except Exception as e:
        raise EcmfwRequestError(
            'Error sending request to ECMWF data store.',
            successful_requests=successful_requests,
            remotes=remotes,
        ) from e
    return remotes
###END def send_ecmwf_datastore_request
