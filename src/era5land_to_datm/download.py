"""Code to download ERA5 Land GRIB files from completed requests.

Classes
-------
DownloadFilesResult
    Represents the result of attempting to download files for multiple requests
    through the `retrieve_available_files` function. Provides information on
    which files were successfully downloaded, which were not available yet, and
    which downloads failed.

Protocols
---------
FilenameFormatCallable
    A protocol defining a callable that takes a VarSet and YearMonth and
    returns a string filename. Used for custom naming of downloaded GRIB files.

Functions
---------
get_remotes
    Retrieve Remote instances for previously sent requests, as a dictionary
    keyed by request id.
get_remote_varset
    Get the VarSet associated with a given Remote instance (i.e., the set of
    variables requested in that Remote).
get_remote_yearmonth
    Get the YearMonth associated with a given Remote instance (i.e., the year
    and month requested in that Remote).
make_grib_filename
    Create a standardized filename for the downloaded GRIB file based on the
    set of variables and year/month requested.
remotes_dict_by_vars_and_yearmonth
    Organize Remote instances into a nested dictionary keyed by VarSet and
    YearMonth, based on the `variable`, `year`, and `month` attributes of each
    Remote.
retrieve_available_files
    Attempt to download available files for multiple Remote instances,
    returning a DownloadFilesResult with details on successes, unavailable
    files, and failures.

Attributes
----------
GRIB_FILENAME_FORMAT_DEFAULT : str | FilenameFormatCallable
    The default format for naming downloaded GRIB files, which is used by the
    `make_grib_filename` function if no custom format is provided.
"""
import logging
import string
import typing as tp

import ecmwf.datastores as ecmwfds

from .types import (
    YearMonth,
)
from .variables import (
    VarSet,
)



logger: logging.Logger = logging.getLogger(__name__)


class FilenameFormatCallable(tp.Protocol):
    """A protocol defining a callable that takes a VarSet and YearMonth and
    returns a string filename. Used for custom naming of downloaded GRIB files.
    """

    def __call__(
            self,
            var_set: VarSet,
            year_month: YearMonth,
    ) -> str:
        """Generate a filename based on the provided VarSet and YearMonth.

        Parameters
        ----------
        var_set : VarSet
            The set of variables for which the file is being named.
        year_month : YearMonth
            The year and month for which the file is being named.

        Returns
        -------
        str
            The generated filename.
        """
        ...

    ###END def FilenameFormatCallable.__call__

###END class FilenameFormatCallable


GRIB_FILENAME_FORMAT_DEFAULT: str | FilenameFormatCallable = (
    'era5land_{variables}_{year:04d}_{month:02d}.grib'
)


def make_grib_filename(
        var_set: VarSet,
        year_month: YearMonth,
        filename_format: str | FilenameFormatCallable = (
            GRIB_FILENAME_FORMAT_DEFAULT
        ),
) -> str:
    """Create a standardized filename for the downloaded GRIB file based on the
    set of variables and year/month requested.

    Parameters
    ----------
    var_set : VarSet
        The set of variables for which the file is being named, as a set or
        frozenset.
    year_month : YearMonth
        The year and month for which the file is being named, as a named tuple
        with fields `year` and `month`.
    filename_format : str | FilenameFormatCallable, optional
        The format to use for the filename. If a string is provided, it should
        be a format string compatible with `str.format()`, using the fields
        `variables`, `year`, and `month`. These fields will be substituted by the
        following:
          * `variables`: The shortnames (variable shortnames used in the grib files)
            of each variable in alphanumeric order, joined by underscores.
          * `year`: The year as an integer (formatting will be applied as per
            the format string).
          * `month`: The month as an integer (formatting will be applied as per
            the format string).
        If a callable is provided, it should conform to the
        `FilenameFormatCallable` protocol, which receives the variables as a
        `VarSet` (frozenset) and the year and month as a `YearMonth` (named
        tuple), in parameters named `var_set` and `year_month`, respecitvely.
        Default is `GRIB_FILENAME_FORMAT_DEFAULT`, which is a format string of
        the form `'era5land_{variables}_{year:04d}_{month:02d}.grib'`.

    Returns
    -------
    str
        The generated filename based on the provided VarSet and YearMonth.
    """
    if isinstance(filename_format, str):
        vars_str: str = '_'.join(
            sorted(_var for _var in var_set)
        )
        filename: str = filename_format.format(
            variables=vars_str,
            year=year_month.year,
            month=year_month.month,
        )
    elif callable(filename_format):
        filename: str = filename_format(
            var_set=var_set,
            year_month=year_month,
        )
    else:
        raise TypeError(
            'filename_format must be a str or a callable conforming to the '
            'FilenameFormatCallable protocol.'
        )
    return filename
###END def make_grib_filename


def get_remotes(
        *,
        limit: int = 100,
) -> dict[str, ecmwfds.Remote]:
    """Retrieve Remote instances for previously sent requests.

    Parameters
    ----------
    limit : int, optional
        The maximum number of Remote instances to retrieve. Default is 100. A
        log message at log level DEBUG will state what limit was used, and will
        be followed by a message at level WARNING if there are more the
        number of instances hits this limit (which may mean there are more jobs
        on the server that were not retrieved).

    Returns
    -------
    dict[str, ecmwfds.Remote]
        A dictionary mapping request ids to their corresponding Remote
        instances.
    """
    client: ecmwfds.Client = ecmwfds.Client()
    logger.debug(f'Retrieving up to {limit} request IDs from server...')
    request_ids: list[str] = client.get_jobs(limit=limit).request_ids
    if len(request_ids) >= limit:
        logger.warning(
            f'Received {len(request_ids)} request IDs from server, the same '
            'number as the limit. There may be more requests on the server.'
        )
    else:
        logger.debug(
            f'Received {len(request_ids)} request IDs from server.'
        )
    remotes: dict[str, ecmwfds.Remote] = {}
    for request_id in request_ids:
        logger.debug(f'Retrieving Remote for request ID {request_id}...')
        remotes[request_id] = client.get_remote(request_id)
    return remotes
###END def get_remotes


def get_remote_varset(
        remote: ecmwfds.Remote,
) -> 'VarSet':
    """Get the VarSet associated with a given Remote instance.

    Parameters
    ----------
    remote : ecmwfds.Remote
        The Remote instance to get the VarSet for.

    Returns
    -------
    VarSet
        The VarSet corresponding to the variables requested in the Remote.
    """
    if not isinstance(remote, ecmwfds.Remote):
        raise TypeError(
            f'Expected ecmwfds.Remote, got {type(remote)}'
        )
    return VarSet(remote.request['variable'])
###END def get_remote_varset


def get_remote_yearmonth(
        remote: ecmwfds.Remote,
) -> YearMonth:
    """Get the YearMonth associated with a given Remote instance.

    Parameters
    ----------
    remote : ecmwfds.Remote
        The Remote instance to get the YearMonth for.

    Returns
    -------
    YearMonth
        The YearMonth corresponding to the year and month requested in the
        Remote.
    """
    if not isinstance(remote, ecmwfds.Remote):
        raise TypeError(
            f'Expected ecmwfds.Remote, got {type(remote)}'
        )
    year: int = int(remote.request['year'])
    month: int = int(remote.request['month'])
    return YearMonth(year=year, month=month)
###END def get_remote_yearmonth


def remotes_dict_by_vars_and_yearmonth(
        remotes: tp.Iterable[ecmwfds.Remote],
) -> dict[VarSet, dict[YearMonth, ecmwfds.Remote]]:
    """Organize Remote instances into a nested dictionary keyed by VarSet and
    YearMonth.

    Parameters
    ----------
    remotes : Iterable[ecmwfds.Remote]
        An iterable sequence of Remote instances to organize.

    Returns
    -------
    dict[VarSet, dict[YearMonth, ecmwfds.Remote]]
        A nested dictionary where the outer keys are VarSet instances, the inner
        keys are YearMonth instances, and the values are the corresponding
        Remote instances.
    """
    remotes_dict: dict[VarSet, dict[YearMonth, ecmwfds.Remote]] = {}
    for remote in remotes:
        var_set: VarSet = get_remote_varset(remote)
        year_month: YearMonth = get_remote_yearmonth(remote)
        remotes_dict.setdefault(var_set, {})[year_month] = remote
    return remotes_dict
###END def remotes_dict_by_vars_and_yearmonth
