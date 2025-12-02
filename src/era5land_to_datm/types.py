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
import enum
import typing as tp

import ecmwf.datastores as ecmwfds



class YearMonth(tp.NamedTuple):
    """A namedtuple representing a year and month.

    Attributes
    ----------
    year : int
        The year (e.g., 2020).
    month : int
        The month (1-12).
    """
    year: int
    month: int
###END class YearMonth


class RemoteErrorTuple(tp.NamedTuple):
    """A namedtuple representing an error associated with a Remote instance.

    Attributes
    ----------
    remote : ecmwf.datastores.Remote
        The Remote instance associated with the error.
    error : Exception
        The exception representing the error.
    """
    remote: ecmwfds.Remote
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
