"""Common logging functionality for ERA5-Land to DATM scripts.

Functions
---------
initialize_logger
    Initializes the root logger, by default to log both to the console and to a
    file in the current working directory.

Attributes
----------
default_log_filename_prefix : str
    The default prefix for log filenames created by initialize_logger. The
    prefix will be suffixed by an underscore, timestamp in the format
    YYYYMMDDTHHMMSS (UTC), and '.log' extension.
is_initialized : bool
    Whether the logger has been initialized.
"""
from datetime import (
    datetime,
    timezone,
)
import logging
import os
from pathlib import Path
import sys
import typing as tp

default_log_filename_prefix: str = 'era5land_to_datm_script'

is_initialized: bool = False

def initialize_logger(
    *,
    log_to_console: bool = True,
    log_to_file: bool = True,
    log_filename_prefix: str = default_log_filename_prefix,
    log_dir: Path|str|None = None,
    log_level: int = logging.INFO,
    override_existing: bool = False,
) -> logging.Logger:
    """Initializes the root logger.

    By default, logs to both the console (stdout) and to a file in the current
    working directory. The log filename is constructed by appending an
    underscore, timestamp in the format YYYYMMDDTHHMMSS (UTC), and '.log'
    extension to the given prefix.

    Parameters
    ----------
    log_to_console : bool, optional
        Whether to log to the console (stdout). Default is True.
    log_to_file : bool, optional
        Whether to log to a file. Default is True.
    log_filename_prefix : str, optional
        The prefix for the log filename. Default is
        'era5land_to_datm_script'.
    log_dir : Path | str | None, optional
        The directory to save the log file in. If None, uses the current working
        directory. Default is None.
    log_level : int, optional
        The logging level. Default is logging.INFO.
    override_existing : bool, optional
        If True, allows re-initialization of the logger. Default is False. If
        False, calling this function multiple times will have no effect after
        the first call, and the same logger instance will be returned.

    Returns
    -------
    logging.Logger
        The initialized root logger.
    """
    global is_initialized
    if is_initialized and not override_existing:
        return logging.getLogger()
    logger: logging.Logger = logging.getLogger()
    logger.setLevel(log_level)

    formatter = logging.Formatter(
        fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%dT%H:%M:%S%z'
    )

    log_path: Path|None = None

    if log_to_console:
        console_handler: logging.StreamHandler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    if log_to_file:
        timestamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
        log_filename = f'{log_filename_prefix}_{timestamp}.log'
        if log_dir is not None:
            log_path = Path(log_dir) / log_filename
        else:
            log_path = Path.cwd() / log_filename
        file_handler: logging.FileHandler = logging.FileHandler(log_path)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.info('Logger initialized.')
    if log_to_file:
        assert log_path is not None  # for mypy
        logger.info(f'Logging to file: {log_path.resolve()}')
    is_initialized = True
    return logger
###END def initialize_logger
