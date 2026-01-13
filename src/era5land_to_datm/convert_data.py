"""Functions to convert ERA5 Land data to DATM data structures.

Functions
---------
make_datm_ds
    Creates a DATM xarray Dataset from an opened ERA5 Land xarray Dataset.
decumulate_era5land_var
    Differences cumulative ERA5 Land variables along the intra-day step
    dimension, to produce a DataArray where each value gives the cumulated value
    for the previous time step only. The value of each time step can then be
    averaged with the next time step and divided by 2 times the time step length
    to produce an average rate value.
"""
import xarray as xr

from .datm_streams import (
    Datm7Stream,
    datm7_stream_variables,
)
from .dimensions import (
    Era5LandDim,
    ERA5LandTimeLayout,
)
from .variables import (
    Datm7Attr,
    Datm7Coord,
    Datm7Var,
    Era5LandVar,
    datm7_required_era5_vars,
    era5_datm_vars,
    era5land_grib_varnames,
    era5land_grib_varnames_reverse,
    era5_cumulative_vars,
)



def make_datm_ds(
        source: xr.Dataset,
        *,
        target_stream: Datm7Stream,
        eager: bool = True,
        time_layout: ERA5LandTimeLayout | None = None,
) -> xr.Dataset:
    """Creates a DATM xarray Dataset from an opened ERA5 Land xarray Dataset.

    Parameters
    ----------
    source : xr.Dataset
        The opened ERA5 Land xarray Dataset. The dataset should have been opened
        using `open_era5land_grib`, with both the main source file and the
        "next" file if `target_stream` includes cumulative variables. If the
        dataset has been opened in some other way, the cumulative variable
        values for the final time step must have been added in some other way.
        At present, only datasets with the two-dimensional date+intradate
        timestep layout are supported. The dataset will be converted to linear
        time layout during the processing for compatilibility with DATM.
    target_stream : Datm7Stream
        The target DATM7 stream to create.
    eager: bool, optional
        Whether to load the entire dataset into memory before processing. This
        can significantly speed up some operations, while setting it to False
        may lead require less memory usage but lead to some data being reloaded
        and even reprocesseed multiple times, which can slow down the
        processing. By default True. Set to False if you run into memory issues.
    time_layout : ERA5LandTimeLayout | None, optional
        The time layout of the source ERA5 Land dataset. **Note:** Currently,
        only the two-dimensional date+intradate timestep layout is supported
        (`ERA5LandTimeLayout.DATE_STEP`). If None, that layout will be assumed,
        but in future versions, None may be used to auto-detect the layout.
        Passing `ERA5LandTimeLayout.LINEAR` will raise a NotImplementedError.
    """
    if time_layout is None:
        time_layout = ERA5LandTimeLayout.DATE_STEP
    if time_layout != ERA5LandTimeLayout.DATE_STEP:
        raise NotImplementedError(
            'Currently, only ERA5 Land datasets with the date+intradate '
            'timestep layout are supported.'
        )

    target_vars: frozenset[Datm7Var] = datm7_stream_variables[target_stream]
    required_era5_vars: frozenset[Era5LandVar] = frozenset(
        _var for _target_var in target_vars
        for _var in datm7_required_era5_vars[_target_var]
    )
    required_era5_girb_varnames: frozenset[str] = frozenset(
        era5land_grib_varnames[_var]
        for _var in required_era5_vars
    )
    missing_source_vars = required_era5_girb_varnames - set(
        source.data_vars.keys()
    )
    if missing_source_vars:
        raise ValueError(
            'The source ERA5 Land dataset is missing the following required '
            f'variables for target stream {target_stream}: '
            f'{missing_source_vars}'
        )

    target_ds: xr.Dataset = make_datm_base(
        source=source,
        target_stream=target_stream,
    )

    cumulative_required_vars: frozenset[Era5LandVar] = (
        required_era5_vars & era5_cumulative_vars
    )
    decumulated_source: xr.Dataset = source.copy(deep=False)
    for _cum_var in cumulative_required_vars:
        _cum_varname = era5land_grib_varnames[_cum_var]
        decumulated_source[_cum_varname] = decumulate_era5land_var(
            source[_cum_varname],
        )
    for _target_var in target_vars:
        target_ds[_target_var.value] = make_target_var(
            target_var=_target_var,
            source=decumulated_source,
        )
    target_ds = postprocess_converted_datm_ds(
        target_ds=target_ds,
        target_stream=target_stream,
        source=source,
        eager=eager,
    )
    return target_ds
###END def make_datm_ds


def decumulate_era5land_var(
        source: xr.DataArray,
        *,
        step_dim: str = Era5LandDim.STEP,
) -> xr.DataArray:
    """Differences cumulative ERA5 Land variables along the intra-day step
    dimension, to produce a DataArray where each value gives the cumulated value
    for the previous time step only. The value of each time step can then be
    averaged with the next time step and divided by 2 times the time step length
    to produce an average rate value.

    Parameters
    ----------
    source : xr.DataArray
        The source ERA5 Land cumulative variable DataArray. Must have a step
        dimension.
    step_dim : str, optional
        The name of the intra-day step dimension. By default `Era5LandDim.STEP`.
    """
    if step_dim not in source.dims:
        raise ValueError(
            f'The source DataArray must have a step dimension {step_dim}.'
        )
    decumulated: xr.DataArray = xr.concat(
        [
            source.isel({step_dim: 0}),
            source.diff(
                dim=step_dim,
                label='upper',
            ),
        ],
        dim=step_dim,
    )
    return decumulated
###END def decumulate_era5land_var
