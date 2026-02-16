"""Custom types used in era5land_to_datm package.

Types
-----
YearMonth
    A namedtuple representing a year and month (e.g., (2020, 1) for January
    2020).
VarSet
    A frozenset of Era5LandVar instances representing a set of variables.

Enums
-----
EcmwfDatasetId
    Enumeration of ECMWF dataset IDs.
Era5LandVar
    Enumeration of ERA5-Land variables available for download.
"""
from collections.abc import Callable
import datetime
import enum
import typing as tp

import ecmwf.datastores as ecmwfds
import numpy as np



class YearMonth(tp.NamedTuple):
    """A namedtuple representing a year and month.

    Attributes
    ----------
    year : int
        The year (e.g., 2020).
    month : int
        The month (1-12).

    Class methods
    -------------
    range(start: YearMonth, end: YearMonth) -> Iterable[YearMonth]
        Generate a sequence of YearMonth instances from start to end, inclusive.
    """
    year: int
    month: int

    @classmethod
    def range(
            cls,
            start: tp.Self,
            end: tp.Self,
    ) -> tp.Iterator[tp.Self]:
        """Generate a sequence of YearMonth instances from start to end, inclusive.

        Parameters
        ----------
        start : YearMonth
            The starting YearMonth (inclusive).
        end : YearMonth
            The ending YearMonth (inclusive).

        Yields
        ------
        YearMonth
            An iterator of YearMonth instances from start to end, inclusive.
        """
        if (start.year, start.month) > (end.year, end.month):
            raise ValueError('start must be less than or equal to end')
        current_year: int = start.year
        current_month: int = start.month
        while (current_year, current_month) <= (end.year, end.month):
            yield cls(year=current_year, month=current_month)
            if current_month == 12:
                current_month = 1
                current_year += 1
            else:
                current_month += 1
###END class YearMonth


class RemoteErrorTuple[RemoteT: ecmwfds.Remote](tp.NamedTuple):
    """A namedtuple representing an error associated with a Remote instance.

    Attributes
    ----------
    remote : ecmwf.datastores.Remote
        The Remote instance associated with the error.
    error : Exception
        The exception representing the error.
    """
    remote: RemoteT
    error: Exception
###END class RemoteErrorTuple


class EcmwfDatasetId(enum.StrEnum):
    """Enumeration of ECMWF dataset IDs.

    Only select datasets are included. A full list can be found using
    `ecmwf.datastores.Client.get_collections()`. The attribute `collection_ids`
    of the returned Collections instance contains all dataset IDs. `Collection`
    instances with full data for each collection can then be obtained by calling
    `ecmwf.datastores.Client.get_collection(<dataset_id>)`, or by inspecing the
    contents of `.json['collections']` of the Collections instance.

    Properties
    ----------
    collection_id : str
        The collection ID corresponding to the dataset ID of a given member.
    collection_title : str
        The collection title corresponding to the dataset ID of a given member.
    """

    ERA5LANDHRLY = 'reanalysis-era5-land'

    @property
    def collection_id(self) -> str:
        """The collection ID corresponding to this dataset ID."""
        return self.value
    ###END def EcmwfDatasetId.collection_id

    @property
    def collection_title(self) -> str:
        """The collection title corresponding to this dataset ID."""
        return _ecmwf_collection_titles[self]
    ###END def EcmwfDatasetId.collection_title

###END class EcmwfDatasetId

_ecmwf_collection_titles: tp.Final[dict[EcmwfDatasetId, str]] = {
    EcmwfDatasetId.ERA5LANDHRLY: 'ERA5-Land hourly data from 1950 to present',
}


DATM7_DATAVAR_DTYPE: tp.Final[np.dtype] = np.dtype('float32')
DATM7_COORD_DTYPE: tp.Final[np.dtype] = np.dtype('float32')
DATM7_TIME_DTYPE: tp.Final[np.dtype] = np.dtype('float64')
DATM7_CALENDAR: tp.Final[str] = 'noleap'
DATM7_NC_FILE_FORMAT: tp.Final[str] = 'NETCDF4'

class NcTimeEncoding(tp.TypedDict):
    units: str | Callable[[np.datetime64], str]
    calendar: str
###END class NcTimeEncoding


def make_datm7_time_units(start_time: np.datetime64) -> str:
    """Make time units string for DATM7 time coordinate, given the start time.

    Parameters
    ----------
    start_time : np.datetime64
        The start time of the dataset, which will be used as the reference time
        in the returned time units string.

    Returns
    -------
    str
        A time units string for the DATM7 time coordinate, with the reference
        time set to the given start_time. The units will be days since the start
        of the day to which `start_time` belongs in UTC.
    """
    start_time_string: str = (
        start_time
        .astype(datetime.datetime)
        .strftime('%Y-%m-%d')
    )
    return f'days since {start_time_string}'
