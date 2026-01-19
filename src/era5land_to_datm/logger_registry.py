"""Module to register loggers and set package-wide logging settings.

The module maintains a registry of loggers used in this package, in order to
coordinate log levels and other settings across Loggers in different modules.

The module defines a default logging level. This level will be applied to all
loggers that register if they have not already set an explicit level (i.e., the
Logger's `level` attribute is `logging.NOTSET`). If an explicit level has been
set at the time of registration, that level will not be overridden at
registration time, but will be overwritten later if the `set_logging_level`
function is called without specifying specific Loggers.

The default logging level is initially set to `logging.INFO`, but can be changed
using the `set_logging_level` function.

Functions
---------
register_logger:
    Register and add a logger to the logging system of this package.
get_logger:
    Get a logger registered in this package by its name.
get_loggers:
    Get all loggers registered in this package.
set_logging_level:
    Set the logging level for all loggers in this package, or optionally for a
    specific set of loggers by name.
get_current_logging_level:
    Get the current default logging level.
"""
from collections.abc import (
    Iterable,
)
import logging



_registered_loggers: dict[str, logging.Logger] = {}
_default_logging_level: int = logging.INFO


def register_logger(logger: logging.Logger) -> None:
    """Register and add a logger to the logging system of this package.

    If the logger has not set an explicit logging level (i.e., its `level`
    attribute is `logging.NOTSET`), the default logging level of this package
    will be applied to it. The current default logging level can be queried
    using the `get_current_logging_level` function.

    Parameters
    ----------
    logger : logging.Logger
        The logger to register.
    """
    global _registered_loggers, _default_logging_level
    if logger.level == logging.NOTSET:
        logger.setLevel(_default_logging_level)
    if (
            logger.name in _registered_loggers
    ) and (
            logger is not _registered_loggers[logger.name]
    ):
        raise ValueError(
            f'A different Logger instance with the name "{logger.name}" is '
            'already registered.'
        )
    _registered_loggers[logger.name] = logger
###END def register_logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger registered in this package by its name.

    Parameters
    ----------
    name : str
        The name of the logger.

    Returns
    -------
    logging.Logger
        The logger with the given name.

    Raises
    ------
    KeyError
        If no logger with the given name is registered.
    """
    return _registered_loggers[name]
###END def get_logger

def get_loggers() -> dict[str, logging.Logger]:
    """Get all loggers registered in this package.

    Returns
    -------
    dict[str, logging.Logger]
        A dictionary mapping logger names to logger instances. Note that a new
        dictionary is returned each time (but the instances it contains are the
        same). Modifying the returned dictionary will not affect the registry,
        but modifying the state of the Logger instances will affect the
        registered Loggers themselves.
    """
    return dict(_registered_loggers)
###END def get_loggers


def set_logging_level(
        level: int,
        *,
        logger_names: Iterable[str] | None = None,
) -> None:
    """Set either the default logging level for this package, or the logging
    level for a specific set of loggers by name.

    Parameters
    ----------
    level : int
        The logging level to set. Must be one of the logging levels defined in
        the `logging` module (e.g., `logging.DEBUG`, `logging.INFO`,
        `logging.WARNING`, etc., as integer values, or using the constants
        defined in the `logging` module).
    logger_names : Iterable[str] | None, optional
        The names of the loggers to set the logging level for. If None (the
        default), the default logging level for this package will be set to the
        given level, and the logging level of all registered loggers will be set
        to the given level.
    """
    global _registered_loggers, _default_logging_level
    if logger_names is None:
        _default_logging_level = level
        for _logger in _registered_loggers.values():
            _logger.setLevel(level)
    else:
        missing_names: set[str] = set(
            name for name in logger_names if name not in _registered_loggers
        )
        if len(missing_names) > 0:
            raise KeyError(
                f'No loggers found with the following names: '
                f'{", ".join(missing_names)}'
            )
        for _name in logger_names:
            _logger = _registered_loggers[_name]
            _logger.setLevel(level)
###END def set_logging_level


def get_current_logging_level() -> int:
    """Get the current default logging level for this package.

    Returns
    -------
    int
        The current default logging level.
    """
    global _default_logging_level
    return _default_logging_level
###END def get_current_logging_level
