"""Code to download ERA5 Land GRIB files from completed requests.

Classes
-------
DownloadFilesResult
    Represents the result of attempting to download files for multiple requests
    through the `retrieve_available_files` function. Provides information on
    which files were successfully downloaded, which were not available yet, and
    which downloads failed.

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
"""
import logging
import typing as tp

import ecmwf.datastores as ecmwfds

from .variables import (
    VarSet,
)



logger: logging.Logger = logging.getLogger(__name__)


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
