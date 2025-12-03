"""Extensions to the ecmwf.datastores Remote class and related functionality.

Classes
-------
CachedClient
    A subclass of ecmwf.datastores.Client that caches Remote instances to
    avoid repeated API calls to access the same Remote.
CachedRemote
    A subclass of ecmwf.datastores.Remote that caches the `request` attribute,
    to avoid repeated API calls to access properties that usually do not change.

Functions
---------
get_client
    Get an ecmwf.datastores.Client instance for making requests. The function
    accepts the same parameters as the Client constructor, but caches the
    created Client instance for reuse on subsequent calls, so that only one
    Client instance is created for a given set of parameters.
"""
import functools
import inspect
import typing as tp

import attrs
import ecmwf.datastores as ecmwfds



class CachedRemote(ecmwfds.Remote):
    """A subclass of ecmwf.datastores.Remote that caches the `request` attribute,
    to avoid repeated API calls to access properties that usually do not change.

    Properties
    ----------
    request : dict[str, tp.Any]
        The request dictionary for this Remote instance. The request dictionary
        is cached after the first access to avoid repeated API calls to access
        properties that usually do not change.

    Methods
    -------
    clear_cached_request
        Clear the cached request dictionary, forcing a refresh on the next
        access to the `request` property.

    Class Methods
    -----------------------
    from_plain_remote
        Create a CachedRemote instance from a plain ecmwf.datastores.Remote
        instance.
    """

    __slots__ = ('_cached_request',)

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._cached_request: dict[str, tp.Any]|None = None
    ###END def CachedRemote.__init__

    @property
    def request(self) -> dict[str, tp.Any]:
        """The request dictionary for this Remote instance.

        The request dictionary is cached after the first access to avoid
        repeated API calls to access properties that usually do not change.
        """
        if self._cached_request is None:
            self._cached_request = super().request
        return self._cached_request
    ###END def CachedRemote.request

    def clear_cached_request(self) -> None:
        """Clear the cached request dictionary, forcing a refresh on the next
        access to the `request` property.
        """
        self._cached_request = None
    ###END def CachedRemote.clear_cached_request

    @classmethod
    def from_plain_remote(
            cls,
            remote: ecmwfds.Remote
    ) -> 'CachedRemote':
        """Create a CachedRemote instance from a plain ecmwf.datastores.Remote
        instance.

        NB! This method is initialized with the attributes of the plain Remote,
        and does not perform a deep copy. Any changes to mutable attributes of
        the plain Remote after calling this method will also affect the created
        CachedRemote instance, and vice-versa. Ideally, the caller should delete
        the plain Remote instance after calling this method or let it pass out
        of scope and be garbage-collected.

        Parameters
        ----------
        remote : ecmwf.datastores.Remote
            The plain Remote instance to convert.
        """
        remote_attrs: dict[str, tp.Any] = attrs.asdict(remote)
        cached_remote: tp.Self = cls(**remote_attrs)
        return cached_remote
    ###END def CachedRemote.from_plain_remote

###END class CachedRemote


class CachedClient(ecmwfds.Client):
    """A subclass of ecmwf.datastores.Client that caches Remote instances to
    avoid repeated API calls to access the same Remote.
    """

    __slots__ = ('_remote_cache',)

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._remote_cache: dict[str, CachedRemote] = {}
    ###END def CachedClient.__init__

    def get_remote(self, request_id: str) -> 'CachedRemote':
        """Get a CachedRemote instance for the given request ID.

        If a CachedRemote instance for the given request ID has already been
        created, it is returned from the cache. Otherwise, a new CachedRemote
        instance is created and stored in the cache before being returned.
        """
        if request_id not in self._remote_cache:
            self._remote_cache[request_id] = CachedRemote.from_plain_remote(
                super().get_remote(request_id)
            )
        return self._remote_cache[request_id]
    ###END def CachedClient.get_remote

###END class CachedClient


@functools.lru_cache(maxsize=None)
def get_client(
        url: str|None = None,
        key: str|None = None,
        verify: bool = True,
        timeout: float | tuple[float, float] = 60,
        progress: bool = True,
        cleanup: bool = True,
        sleep_max: float = 120,
        retry_after: float = 120,
        maximum_tries: int = 500,
) -> CachedClient:
    """Get a CachedClient instance for making requests.
    
    The function accepts the same parameters as the Client constructor, except
    for the non-hashable `request` parameter. Note that the defaults can be
    different.All parameters are passed directly to the Client constructor.
    Please see the documentation of `ecmwf.datastores.Client` for details.
    """
    local_values: dict[str, tp.Any] = locals()
    parameters: dict[str, tp.Any] = {
        _param.name: _param for _param in inspect.signature(
            ecmwfds.Client.__init__
        ).parameters.values()
    }
    init_args: dict[str, tp.Any] = {}
    for param_name in parameters:
        if param_name == 'request':
            continue
        init_args[param_name] = locals()[param_name]
    client: CachedClient = CachedClient(**init_args)
    return client
###END def get_client
