"""File input and output of ERA5 and DATM data files.

Functions
---------
open_era5land_grib
    Opens an ERA5 Land GRIB file, optionally with the next file to obtain values
    in the last time step when these are NaN. Returns data as an xarray Dataset.
write_datm_nc
    Writes a DATM netCDF file for a given stream from a given xarray Dataset.
make_datm7_encoding_dict
    Make a dictionary of encoding settings to use when writing the output DATM7.
"""
from collections.abc import (
    Callable,
    Hashable,
)
import logging
from pathlib import Path
import typing as tp

import numpy as np
import xarray as xr

from .dimensions import (
    Datm7Dim,
    Era5LandDim,
)
from .datm_streams import Datm7Stream
from .types import (
    DATM7_CALENDAR,
    DATM7_COORD_DTYPE,
    DATM7_DATAVAR_DTYPE,
    DATM7_NC_FILE_FORMAT,
    DATM7_TIME_DTYPE,
    make_datm7_time_units,
)



logger = logging.getLogger(__name__)


def open_era5land_grib(
        file: Path,
        *,
        next_file: Path|None = None,
        previous_file: Path|None = None,
        chunks: dict|int|tp.Literal['auto']|None = 'auto',
        date_dim: str = Era5LandDim.DATE,
        step_dim: str = Era5LandDim.STEP,
) -> xr.Dataset:
    """Opens an ERA5 Land GRIB file, optionally with the next file for
    cumulative variables.

    Parameters
    ----------
    file : Path
        Path to the ERA5 Land GRIB file.
    next_file : Path | None, optional
        Path to the next ERA5 Land GRIB file for cumulative variables. If given,
        it will be lazily opened, and the first time step will be extracted and
        concatenated to the main dataset along the time dimension. Optional, by
        default None.
    previous_file : Path | None, optional
        Path to the previous ERA5 Land GRIB file for cumulative variables. If
        given, it will be lazily opened, and the last date will be extracted and
        concatenated to the main dataset along the time dimension. Optional, by
        default None.
    chunks : dict | int | 'auto' | None, optional
        Chunking option for xarray when opening the dataset. Passed directly to
        `xarray.open_dataset()`. By default 'auto'. To disable chunking, set to
        False. In order to use xarray's lazy loading without dask, set to None.
    date_dim : str, optional
        Name of the time/date dimension in the ERA5 Land dataset. By default
        given by the string enum `Era5LandDim.DATE`.
    step_dim : str, optional
        Name of the "step" dimension (intra-date timediff relative to midnight)
        in the ERA5 Land dataset. By default given by the string enum
        `Era5LandDim.STEP`.
    Returns
    -------
    xr.Dataset
        The ERA5 Land dataset.

    Raises
    ------
    ValueError
        If the last date in the source file does not match the first date in the
        next file, if the last intra-date time step in the source file does not
        match the *last* intra-date time step in the next file, or if the next
        file contains variables that are not present in the source file.
    """
    era5_ds: xr.Dataset = xr.open_dataset(
        file,
        chunks=chunks,
    )
    if next_file is not None:
        next_ds: xr.Dataset = xr.open_dataset(
            next_file,
            chunks='auto',
        ).isel({date_dim: 0, step_dim: -1})
        source_last_date = era5_ds[date_dim].isel({date_dim: -1})
        next_first_date = next_ds[date_dim]
        if source_last_date != next_first_date:
            raise ValueError(
                f'The last date in the source file ({source_last_date}) does not '
                f'match the first date in the next file ({next_first_date}).'
            )
        if not set(era5_ds.variables).issuperset(set(next_ds.variables)):
            raise ValueError(
                'The next file contains variables that are not present in the '
                'source file.'
            )
        source_last_step = era5_ds[step_dim].isel({step_dim: -1})
        next_last_step = next_ds[step_dim]
        if source_last_step != next_last_step:
            raise ValueError(
                f'The last intra-date time step in the source file '
                f'({source_last_step}) does not match the last intra-date time '
                f'step in the next file ({next_last_step}).'
            )
        for _var in era5_ds.data_vars:
            era5_ds[_var].loc[
                {
                    date_dim: source_last_date,
                    step_dim: source_last_step,
                }
            ] = next_ds[_var]
    if previous_file is not None:
        previous_ds: xr.Dataset = xr.open_dataset(
            previous_file,
            chunks='auto',
        ).isel({date_dim: -1})
        source_first_date = era5_ds[date_dim].isel({date_dim: 0})
        previous_last_date = previous_ds[date_dim]
        if source_first_date != previous_last_date:
            raise ValueError(
                f'The first date in the source file ({source_first_date}) does not '
                f'match the last date in the previous file ({previous_last_date}).'
            )
        if not set(era5_ds.variables).issuperset(set(previous_ds.variables)):
            raise ValueError(
                'The previous file contains variables that are not present in the '
                'source file.'
            )
        source_step_coords = era5_ds[step_dim]
        previous_step_coords = previous_ds[step_dim]
        if not (source_step_coords == previous_step_coords).all():
            error_msg: str = (
                'The intra-date time steps in the source file do not match the '
                'intra-date time steps in the previous file. Source steps: '
                f'{source_step_coords.to_numpy()}, previous steps: '
                f'{previous_step_coords.to_numpy()}.'
            )
            logger.error(
                msg=error_msg,
                extra={
                    'source_step_coords': source_step_coords,
                    'previous_step_coords': previous_step_coords,
                },
            )
            raise ValueError(error_msg)
        for _var in era5_ds.data_vars:
            # We can't index on more than one dimension for a dask-backed array,
            # so we need to first create a combined array for the first date,
            # and then assign that to the first date of the dataset.
            first_date_array: xr.DataArray = (
                era5_ds[_var]
                .sel({date_dim: source_first_date})
                .combine_first(previous_ds[_var])
            )
            era5_ds[_var].loc[
                {
                    date_dim: source_first_date,
                }
            ] = first_date_array
    return era5_ds
###END def open_era5land_grib


def write_datm_nc(
        ds: xr.Dataset,
        *,
        output_file: Path,
        stream: Datm7Stream,
        clobber: bool = False,
        encoding_dict: dict[Hashable, dict[str, str|float]] | None = None,
) -> Path:
    """Writes a DATM netCDF file for a given stream from a given xarray Dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The xarray Dataset containing the data to write.
    output_file : Path
        The path to the output netCDF file.
    stream : Datm7Stream
        The DATM7 stream enum indicating which stream is being written. In the
        current implementation, this parameter is not used, but may be in the
        future. It should therefore be set to proper value, and is mandatory,
        even though it currently has no effect. If necessary, an empty string
        or a dummy value can be used, but this may not be future-proof.
    clobber : bool, optional
        Whether to overwrite an existing file at the output path. By default
        False, in which case a FileExistsError will be raised if the output file
        already exists.
    encoding_dict : dict[Hashable, dict[str, str|float]] | None, optional
        A dictionary of encoding settings to pass to
        `xarray.Dataset.to_netcdf()` when writing the output file. The keys
        should be variable or dimension names in the dataset, and the values
        should be dictionaries of encoding settings for those variables or
        dimensions. If None (the default), a default encoding dictionary will be
        generated using the `make_datm7_encoding_dict()` function, with default
        parameters. If a dictionary is provided, it will be used as is, and the
        `make_datm7_encoding_dict()` function will not be called.

    Returns
    -------
    Path
        The path to the written netCDF file.
    """
    output_file = Path(output_file)
    if output_file.exists() and not clobber:
        raise FileExistsError(
            f'The output file {output_file} already exists. Set the parameter '
            '`clobber=True` if you want to overwrite it.'
        )
    if encoding_dict is None:
        encoding_dict = make_datm7_encoding_dict(ds)
    if encoding_dict.get(Datm7Dim.TIME, {}).get('calendar') == 'noleap':
        # Calendar may have been set to noleap even if the source data contain
        # February 29, so drop that date to avoid errors when writing to netCDF.
        time_arr: xr.DataArray = ds[Datm7Dim.TIME]
        ds = ds.sel(
            {
                Datm7Dim.TIME: (
                    ~((time_arr.dt.month == 2) & (time_arr.dt.day == 29))
                )
            }
        )
    ds.to_netcdf(
        path=output_file,
        mode='w',
        format='NETCDF4',
        encoding=encoding_dict,
    )
    return output_file
###END def write_datm_nc


def make_datm7_encoding_dict(
        output_ds: xr.Dataset,
        *,
        time_dim: str = Datm7Dim.TIME,
        space_dims: tuple[str, str] = (Datm7Dim.LAT, Datm7Dim.LON),
        data_var_dtype: np.dtype = DATM7_DATAVAR_DTYPE,
        coord_dtype: np.dtype = DATM7_COORD_DTYPE,
        time_dtype: np.dtype = DATM7_TIME_DTYPE,
        time_units: str|Callable[[np.datetime64], str] = make_datm7_time_units,
        time_calendar: str = DATM7_CALENDAR,
        nc_file_format: str = DATM7_NC_FILE_FORMAT,
) -> dict[Hashable, dict[str, str|float]]:
    """Make a dictionary of encoding settings for writing output DAT7 netCDF.

    Parameters
    ----------
    output_ds : xr.Dataset
        The xarray Dataset that will be written to netCDF, for which the
        encoding dictionary will be made. **NB!** It is assumed that this
        Dataset is already sorted in ascending time order.
    time_dim : str, optional
        The name of the time dimension in the output dataset. By default given
        by the string enum `Datm7Dim.TIME`.
    space_dims : tuple[str, str], optional
        The names of the spatial dimensions (latitude and longitude) in the
        output dataset, in the order (lat_dim, lon_dim). By default given by the
        string enums `Datm7Dim.LAT` and `Datm7Dim.LON`.
    data_var_dtype : np.dtype, optional
        The data type to use for data variables in the output netCDF file. By
        default given by the constant `DATM7_DATAVAR_DTYPE`.
    coord_dtype : np.dtype, optional
        The data type to use for coordinate variables in the output netCDF file.
        By default given by the constant `DATM7_COORD_DTYPE`.
    time_units : str | Callable[[np.datetime64], str], optional
        The time units to use for the time coordinate in the output netCDF file.
        Can be a string or a callable that takes the start time as a
        np.datetime64 and returns a string. By default given by the function
        `make_datm7_time_units`, which generates a time units string with the
        reference time set to the start of the day of the given start time.
    time_calendar : str, optional
        The calendar to use for the time coordinate in the output netCDF file.
        By default given by the constant `DATM7_CALENDAR`.
    nc_file_format : str, optional
        The netCDF file format to use for the output file. By default given by
        the constant `DATM7_NC_FILE_FORMAT`.
    """
    DTYPE_ATTR: str = 'dtype'
    UNITS_ATTR: str = 'units'
    CALENDAR_ATTR: str = 'calendar'
    FILL_VALUE_ATTR: str = '_FillValue'
    FILL_VALUE_FLOAT: float = 1e36

    start_date: np.datetime64 = np.datetime64(
        (
            output_ds[time_dim]
            .isel({time_dim: 0})
            .compute()
            .item()
        ),
        'ns',
    ).astype('datetime64[s]')  # This weird step-wise conversion is needed to make sure we get something that can be converted to a Python datetime object, and which doesn't get turned into a big int
    if callable(time_units):
        time_units = time_units(start_date)
    time_encoding: dict[str, str|float] = {
        UNITS_ATTR: time_units,
        CALENDAR_ATTR: time_calendar,
        DTYPE_ATTR: str(time_dtype),
        FILL_VALUE_ATTR: FILL_VALUE_FLOAT,
    }
    data_var_encoding: dict[str, str|float] = {
        DTYPE_ATTR: str(data_var_dtype),
        FILL_VALUE_ATTR: FILL_VALUE_FLOAT,
    }
    coord_encoding: dict[str, str|float] = {
        DTYPE_ATTR: str(coord_dtype),
        FILL_VALUE_ATTR: FILL_VALUE_FLOAT,
    }

    encoding: dict[Hashable, dict[str, str|float]] = {
        time_dim: time_encoding,
    } | {
        _dim: coord_encoding for _dim in space_dims
    } | {
        var: data_var_encoding for var in output_ds.data_vars
    }

    return encoding
###END def make_datm7_encoding_dict
