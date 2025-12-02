"""Code to download ERA5 Land GRIB files from completed requests.

Classes
-------
DownloadFilesResult
    Represents the result of attempting to download files for multiple requests
    through the `retrieve_available_files` function. Provides information on
    which files were successfully downloaded, which were not available yet, and
    which downloads failed.
DownloadedFile
    A namedtuple representing a successfully downloaded file, with the Remote
    instance and the local Path as elements (in that order), with field names
    `remote` and `path`, respectively.

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
get_dict_by_vars_and_yearmonth_values
    Convert a nested dictionary of objects keyed by VarSet and YearMonth into a
    list of the nested object values.
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
from dataclasses import dataclass
import logging
import os
from pathlib import Path
import string
import typing as tp

import ecmwf.datastores as ecmwfds

from .types import (
    RemoteErrorTuple,
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


class DownloadedFile(tp.NamedTuple):
    """A namedtuple representing a successfully downloaded file.

    Attributes
    ----------
    remote : ecmwf.datastores.Remote
        The Remote instance associated with the downloaded file.
    path : str
        The local file path where the downloaded file is saved.
    """
    remote: ecmwfds.Remote
    path: Path
###END class DownloadedFile


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
            sorted(_var.grib_varname for _var in var_set)
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
    request: dict[str, tp.Any] = remote.request
    year: int = int(request['year'])
    month: int = int(request['month'])
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


def get_dict_by_vars_and_yearmonth_values[_ObjType](
        nested_dict: dict[VarSet, dict[YearMonth, _ObjType]],
) -> list[_ObjType]:
    """Convert a nested dictionary of objects keyed by VarSet and YearMonth
    into a list of the nested object values.

    Parameters
    ----------
    nested_dict : dict[VarSet, dict[YearMonth, object]]
        A nested dictionary where the outer keys are VarSet instances, the inner
        keys are YearMonth instances, and the values are objects of any single
        type.

    Returns
    -------
    list[_ObjType]
        A list of all the objects contained in the nested
        dictionary.
    """
    return [
        _obj
        for _year_month_dict in nested_dict.values()
        for _obj in _year_month_dict.values()
    ]
###END def get_dict_by_vars_and_yearmonth_values


@dataclass(kw_only=True)
class DownloadFilesResult:
    """Represents the result of attempting to download files for multiple
    requests.

    Attributes
    ----------
    downloaded_files : list[DownloadedFile]
        A list of DownloadedFile instances (namedtuples) representing
        successfully downloaded files, each containing the Remote instance and
        the local Path of the downloaded file.
    existing_files : list[DownloadedFile]
        A list of DownloadedFile instances (namedtuples) representing files that
        were already present locally and therefore not downloaded again, each
        containing the Remote instance and the local Path of the existing file.
    no_files_remotes : list[ecmwf.datastore.Remote]
        A list of Remote instances for which no files were available to
        download.
    failed_remotes : list[RemoteErrorTuple]
        A list of RemoteErrorTuple instances (namedtuples with fields `remote`
        and `exception`) representing Remote instances for which the download
        failed, along with the associated exception.
    """
    downloaded_files: list[DownloadedFile]
    existing_files: list[DownloadedFile]
    no_files_remotes: list[ecmwfds.Remote]
    failed_remotes: list[RemoteErrorTuple]
###END class DownloadFilesResult


def retrieve_available_files(
        remotes: tp.Iterable[ecmwfds.Remote],
        *,
        download_dir: Path|str|None = None,
        filename_format: str | FilenameFormatCallable = (
            GRIB_FILENAME_FORMAT_DEFAULT
        ),
) -> DownloadFilesResult:
    """Attempt to download available files for multiple Remote instances.

    Parameters
    ----------
    remotes : Iterable[ecmwfds.Remote]
        An iterable sequence of Remote instances to attempt to download files
        for.
    download_dir : Path | str | None, optional
        The directory to save downloaded files to. If None, the current working
        directory is used. Default is None.
    filename_format : str | FilenameFormatCallable, optional
        The format to use for naming downloaded GRIB files. The filenames will
        be generated using the `make_grib_filename` function with this format.
        See the documentation of that function for details. The default is
        set by the module attribute `GRIB_FILENAME_FORMAT_DEFAULT`.

    Returns
    -------
    DownloadFilesResult
        An instance of DownloadFilesResult containing details on which files
        were successfully downloaded, which were not available yet, and which
        downloads failed.
    """
    downloaded_files: list[DownloadedFile] = []
    existing_files: list[DownloadedFile] = []
    no_files_remotes: list[ecmwfds.Remote] = []
    failed_remotes: list[RemoteErrorTuple] = []

    if download_dir is not None:
        download_dir_path: Path = Path(download_dir)
    else:
        download_dir_path: Path = Path.cwd()
    if not download_dir_path.exists():
        raise FileNotFoundError(
            f'Download directory {download_dir_path} does not exist.'
        )
    if not download_dir_path.is_dir():
        raise NotADirectoryError(
            f'Download directory {download_dir_path} is not a directory.'
        )
    if not os.access(download_dir_path, os.W_OK):
        raise PermissionError(
            f'Download directory {download_dir_path} is not writeable.'
        )

    for _remote in remotes:
        var_set: VarSet = get_remote_varset(_remote)
        year_month: YearMonth = get_remote_yearmonth(_remote)
        filename: str = make_grib_filename(
            var_set=var_set,
            year_month=year_month,
            filename_format=filename_format,
        )
        filepath: Path = download_dir_path / filename
        if filepath.exists():
            logger.debug(
                f'File {filepath} already exists locally; skipping download '
                f'for Remote with variables {var_set} for '
                f'{year_month.year:04d}-{year_month.month:02d} (delete it if '
                'you want to re-download).'
            )
            existing_files.append(
                DownloadedFile(
                    remote=_remote,
                    path=filepath,
                )
            )
            continue
        try:
            logger.debug(
                f'Checking status of remote with id {_remote.request_id} for '
                f'variables {var_set} in {year_month.year:04d}-'
                f'{year_month.month:02d}...'
            )
            status: str = _remote.status
            if status == 'failed':
                failed_remotes.append(
                    RemoteErrorTuple(
                        remote=_remote,
                        error=RuntimeError(
                            f'Remote with id {_remote.request_id} has status '
                            f'"failed".'
                        ),
                    )
                )
                logger.warning(
                    f'Remote with id {_remote.request_id} has status "failed"; '
                    'skipping download.'
                )
                continue
            if _remote.results_ready:
                logger.debug(
                    f'Results ready for remote with id {_remote.request_id}; '
                    'downloading...'
                )
                returned_path: str = _remote.download(target=str(filepath))
                if Path(returned_path).resolve() != filepath.resolve():
                    logger.error(
                        f'Downloaded file path {returned_path} does not match '
                        f'target path {filepath}; something unexpected '
                        'happened.'
                    )
                downloaded_files.append(
                    DownloadedFile(
                        remote=_remote,
                        path=Path(returned_path),
                    )
                )
            else:
                logger.info(
                    f'No files available yet for remote with id '
                    f'{_remote.request_id}; skipping download.'
                )
                no_files_remotes.append(_remote)
        except Exception as e:
            logger.error(
                f'Error occurred while attempting to download files for '
                f'remote with id {_remote.request_id}: {e}'
            )
            failed_remotes.append(
                RemoteErrorTuple(
                    remote=_remote,
                    error=e,
                )
            )
    return DownloadFilesResult(
        downloaded_files=downloaded_files,
        existing_files=existing_files,
        no_files_remotes=no_files_remotes,
        failed_remotes=failed_remotes,
    )
###END def retrieve_available_files
