"""Functions to convert ERA5 Land data to DATM data structures.

Functions
---------
make_datm_ds
    Creates a DATM xarray Dataset from an opened ERA5 Land xarray Dataset.
make_target_var
    Creates a target DATM7 variable DataArray from the source ERA5 Land
    Dataset.
make_datm_base
    Creates the base DATM xarray Dataset with coordinates, dimensions, and
    Dataset attributes set.
era5land_to_linear_time
    Converts an ERA5 Land Dataset with a two-dimensional date+intradate step
    time layout to a Dataset with a linearized one-dimensional time layout. No
    other changes are made to variable values, coordinates or attributes.
era5land_from_linear_time
    Converts an ERA5 Land Dataset with a linearized one-dimensional time layout
    back to a Dataset with a two-dimensional date+intradate step time layout.
postprocess_converted_datm_ds
    Postprocesses a converted DATM xarray Dataset after all variables have been
    created, to perform any final adjustments needed.
add_target_var_attrs
    A utilitity function that adds attributes to a target DATM7 variable after
    converting/computing the raw values.
set_target_dims_and_coords
    A utility function that sets the dimensions and coordinate variables for
    a target DATM7 variable after converting/computing the raw values (at which
    point it may still have dimensions and coordinates from the source ERA5
    Land Dataset).
compute_average_rate
    Computes average rate values from decumulated ERA5 Land variable values. The
    function assumes that the value at each point is the cumulated value for the
    previous time step only. The average rate is therefore computed by taking
    the average of the value at each point and the value at the next point along
    the time dimension, and dividing by 2 times the time step length.
comopute_specific_humidity
    Computes specific humidity from temperature, pressure, and dewpoint
    temperature DataArrays. This method is used to commpute the DATM7 QBOT
    variable from ERA5 Land T2M, SP, and D2M variables.
decumulate_era5land_var
    Differences cumulative ERA5 Land variables along the intra-day step
    dimension, to produce a DataArray where each value gives the cumulated value
    for the previous time step only. The value of each time step can then be
    averaged with the next time step and divided by 2 times the time step length
    to produce an average rate value.
round_coords
    Rounds coordinate values to the nearest value from a given mapping.

Attributes
----------
value_conversion_funcs : Mapping[Datm7Var, Callable[[xr.Dataset], xr.DataArray]]
    Mapping of DATM7 variables to functions that convert an ERA5 Land Dataset
    to the corresponding DATM7 variable DataArray.

Exception Classes
-----------------
UnexpectedEra5LandUnitsError
    Raised when one or more ERA5 Land variable has unexpected units.
"""
from collections.abc import (
    Callable,
    Hashable,
    Iterable,
    Mapping,
)
import datetime
import functools
import logging
import typing as tp

import numpy as np
import xarray as xr

from .datm_streams import (
    Datm7Stream,
    datm7_stream_variables,
)
from .dimensions import (
    Datm7Dim,
    ERA5_LINEARIZED_TIME_DIM,
    Era5LandDim,
    Era5LandLinearizedTimeDimId,
    ERA5LandTimeLayout,
    LinearizedTimeDimId,
    datm7_dim_order,
    era5_to_datm7_dim_map,
)
from .logger_registry import register_logger
from .meteorology import (
    specific_humidity_from_dewpoint_pressure,
)
from .variables import (
    Datm7Attr,
    Datm7Coord,
    Datm7Var,
    Era5LandCoord,
    Era5LandVar,
    datm7_required_era5_vars,
    datm7_var_attrs,
    era5_datm_vars,
    era5land_grib_varnames,
    era5land_grib_varnames_reverse,
    era5_cumulative_vars,
    era5_var_units,
)



logger = logging.getLogger(__name__)
register_logger(logger)


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
    logger.debug(
        f'Starting `make_datm_ds` for target stream {target_stream}.'
    )
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
    required_era5_girb_varnames: set[Hashable] = set(
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

    check_era5land_units(
        source=source,
        variables=required_era5_vars,
    )

    cumulative_required_vars: frozenset[Era5LandVar] = (
        required_era5_vars & era5_cumulative_vars
    )
    if len(cumulative_required_vars) > 0:
        logger.debug(
            f'Decumulating cumulative ERA5 Land variables for target stream '
            f'{target_stream}: {cumulative_required_vars}'
        )
    source_decumulated: xr.Dataset = source.copy(deep=False)
    for _cum_var in cumulative_required_vars:
        logger.debug(
            f'    Decumulating {era5land_grib_varnames[_cum_var]}...'
        )
        _cum_varname = era5land_grib_varnames[_cum_var]
        source_decumulated[_cum_varname] = decumulate_era5land_var(
            source[_cum_varname],
        )

    if time_layout == ERA5LandTimeLayout.DATE_STEP:
        logger.debug(
            'Converting ERA5 Land dataset to linear time layout...'
        )
        source_1d_time: xr.Dataset = era5land_to_linear_time(
            source=source_decumulated,
        )
    else:
        logger.debug(
            'ERA5 Land dataset already in linear time layout, no conversion '
            'of time dimension layout is needed.'
        )
        source_1d_time: xr.Dataset = source_decumulated
    del source_decumulated

    target_ds: xr.Dataset = make_datm_base(
        source=source_1d_time,
        target_stream=target_stream,
    )
    logger.debug(
        f'Created base DATM dataset for target stream {target_stream}.'
    )

    for _target_var in target_vars:
        logger.debug(
            f'Processing and adding target variable {_target_var}...'
        )
        # Merge the target variable into the target Dataset rather than simply
        # assigning it, to avoid overwriting index coordinate attributes.
        target_ds = target_ds.merge(
            {
                _target_var.value: make_target_var(
                    target_var=_target_var,
                    source=source_1d_time,
                ),
            },
            compat='override',
            join='left',
            combine_attrs='override',
        )
    del source_1d_time
    logger.debug(
        f'Converted all variables, starting postprocessing for target stream '
        f'{target_stream}...'
    )
    target_ds = postprocess_converted_datm_ds(
        target_ds=target_ds,
        target_stream=target_stream,
        source=source,
    )
    logger.debug(
        f'Finished processing target stream {target_stream}, returning from '
        f'`make_datm_ds`.'
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

    NB! The diff method that is applied in this function apperas to fail with
    lazy-loaded DataArrays based on grib files whose data has not been loaded
    yet. The data will in `source` will therefore be loaded and persisted after
    calling this function.

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
    # Load data if not backed by a dask array; Persist it but keep as a dask
    # array if is backed by a dask array.
    if source.chunks is not None:
        source = source.persist()
    else:
        source = source.load()
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


def era5land_to_linear_time(
        source: xr.Dataset,
        *,
        source_date_dim: str = Era5LandDim.DATE.value,
        source_step_dim: str = Era5LandDim.STEP.value,
        output_time_dim: str = ERA5_LINEARIZED_TIME_DIM,
        source_time_coord: str|None = Era5LandCoord.TIME_LINEAR.value,
        preserve_source_time_coord: bool = False,
        preserve_source_time_component_coords: bool = False,
        source_time_component_coords_rename: Mapping[str, str] | None = None,
) -> xr.Dataset:
    """Converts an ERA5 Land Dataset with a two-dimensional date+intradate step
    time layout to a Dataset with a linearized one-dimensional time layout. No
    other changes are made to variable values, coordinates or attributes.

    Parameters
    ----------
    source : xr.Dataset
        The source ERA5 Land Dataset with a two-dimensional date+intradate step
        time layout.
    source_date_dim : str, optional
        The name of the date dimension in the source Dataset. By default
        `Era5LandDim.DATE`.
    source_step_dim : str, optional
        The name of the intradate step dimension in the source Dataset. By
        default `Era5LandDim.STEP`.
    output_time_dim : str, optional
        The name to use for the linearized time dimension in the output Dataset.
        By default `ERA5_LINEARIZED_TIME_DIM`.
    source_time_coord : str|None, optional
        Name of a variable in the `source` that holds full time coordinate
        values for each point (assumed to be of dtype `numpy.datetime64`). This
        is typically the variable `valid_time` in ERA5 Land GRIB files (given
        by the enum `Era5LandCoord.TIME_LINEAR`). The variable must be
        two-dimensional, with dimensions equal to `source_date_dim` and
        `source_step_dim`. If None, the time coordinate variable will be
        computed from the indexes of the date and step dimensions, which will
        require extra computation. By default `Era5LandCoord.TIME_LINEAR`.
    preserve_source_time_coord : bool, optional
        Whether to preserve the source time coordinate variable in the output
        Dataset (in addition to setting a possibly renamed version of it as the
        time index coordinate). If False, the source time coordinate variable
        will have its values preserved in the new time index coordinate, but
        attributes may be dropped (since `xarray.Dataset.set_index` does not
        necessarily preserve attributes). If True, the source time coordinate
        variable will be kept in its original form with original attributes in
        the output Dataset, except for having its dimensions stacked as for the
        rest of the Dataset. By default False.
    preserve_source_time_component_coords : bool, optional
        Whether to preserve the source date and step dimension coordinates in
        the output Dataset. If False, these coordinates will be dropped, and
        only the new linearized time coordinate will be kept. If True, they will
        be kept in the output dataset, but will have their values copied to fill
        every point along the new linearized time dimension.
    source_time_component_coords_rename : Mapping[str, str] | None, optional
        If `preserve_source_time_component_coords` is True, this mappping can be
        used to specify a custom renaming for the original source date and step
        dimension coordinates in the output Dataset. The keys should be the
        original coordinate names and the values the names in the output. If
        None or not specified, they will be renamed to the coordinate names
        given in the `Era5LandLinearizedTimeDimId` enum. If specified, the
        caller is responsible for ensuring both that the keys are the correct
        original names and that the values/target names do not overlap with any
        dimension, coordinate or data variable names in the output. Set to an
        empty dict if you do not want any renaming, but note that this will
        most likely result in an error because of name collisions between the
        source date coordinate and the new time dimension name.

    Returns
    -------
    xr.Dataset
        A new xarray Dataset with the same variables as `source`, but with a
        linearized one-dimensional time dimension instead of the two-dimensional
        date+intradate step layout.
    """
    # We need a temporary time dimension name, since the desired dimension
    # name may collide with one of the source dimensions.
    temp_time_dim_name: str = f'{source_date_dim}_{source_step_dim}_stacked'
    if source_time_coord is None:
        temp_time_coord_name: str \
            = f'{source_date_dim}_{source_step_dim}_time_coord'
        source = source.assign_coords(
            {
                temp_time_coord_name: (
                    source[source_date_dim] + source[source_step_dim]
                )
            }
        )
    else:
        temp_time_coord_name = source_time_coord
    if source_time_component_coords_rename is None:
        source_time_component_coords_rename = {
            source_date_dim: Era5LandLinearizedTimeDimId.DATE.value,
            source_step_dim: Era5LandLinearizedTimeDimId.STEP.value,
        }
    # Define functions to set the index for the new time dimension and to
    # optionally drop the source date and step dimension coordinates, depending
    # on whether we should preserve variables for the original coordinates or
    # not.
    set_index_func: Callable[[xr.Dataset], xr.Dataset] = (
        (
            lambda ds: ds.set_index({temp_time_dim_name: temp_time_coord_name})
        ) if not preserve_source_time_coord else (
            lambda ds: ds.assign_coords(
                {temp_time_dim_name: ds[temp_time_coord_name]}
            )
        )
    )
    transform_time_component_coords_func: Callable[[xr.Dataset], xr.Dataset] = (
        (
            lambda ds: ds.rename(source_time_component_coords_rename)
        )  if preserve_source_time_component_coords else (
            lambda ds: ds.drop_vars(
                [source_date_dim, source_step_dim],
                errors='ignore',
            )
        )
    )
    # Pipe the source Dataset through each of the necessary transformations.
    output_ds: xr.Dataset = (
        source
        .stack(
            {temp_time_dim_name: (source_date_dim, source_step_dim)},
            create_index=False,
        )
        .pipe(set_index_func)
        .pipe(transform_time_component_coords_func)
        .rename({temp_time_dim_name: str(output_time_dim)})
    )
    return output_ds
###END def era5land_to_linear_time


def era5land_from_linear_time(
        source: xr.Dataset,
        *,
        source_time_dim: str = ERA5_LINEARIZED_TIME_DIM,
        target_date_dim: str = Era5LandDim.DATE.value,
        target_step_dim: str = Era5LandDim.STEP.value,
        fast_unstack: bool = False,
        source_date_coord: str|None = Era5LandLinearizedTimeDimId.DATE.value,
        source_step_coord: str|None = Era5LandLinearizedTimeDimId.STEP.value,
        target_time_coord_name: str | None = Era5LandCoord.TIME_LINEAR.value,
) -> xr.Dataset:
    """Converts an ERA5 Land Dataset with a linearized one-dimensional time layout
    back to a Dataset with a two-dimensional date+intradate step time layout. No
    other changes are made to variable values, coordinates or attributes.

    Parameters
    ----------
    source : xr.Dataset
        The source ERA5 Land Dataset with a linearized one-dimensional time
        layout. **NB!** The Dataset must not contain any variables, coordinates
        or dimension names that are equal to any of the names given in the
        other parameters of this function, *or* their defaults if not they are
        not provided explicitly, unless they are used as intended by those
        parameter values. Any name collisions will most likely result in an
        unspecified error, but could also lead to silent errors in the data.
    source_time_dim : str, optional
        The name of the linearized time dimension in the source Dataset. By
        default `ERA5_LINEARIZED_TIME_DIM`.
    target_date_dim : str, optional
        The name to use for the date dimension in the output Dataset. By default
        `Era5LandDim.DATE`.
    target_step_dim : str, optional
        The name to use for the intradate step dimension in the output Dataset.
        By default `Era5LandDim.STEP`.
    fast_unstack : bool, optional
        Whether to assume that the source Dataset has a sorted and complete time
        index, so unstack through direct reshaping rather than by using the
        `xarray.Dataset.unstack` method (which unstacks based on coordinate
        values, and is much slower for a large dataset). If True, the caller
        must ensure that `source` is sorted along the time dimension, that the
        first time point is 01:00 of the first date and the last time point is
        00:00 of the day after the last date, and that the spacing between time
        values is exactly 1 hour with no missing points. If not, the behavior is
        undefined (will probably crash, but can lead to silent errors).
    source_date_coord : str|None, optional
        If the source Dataset has a coordinate variable for the date component
        of the time dimension, the name of that variable. If None, it is assumed
        that there is no such variable. If specified, its values will be used to
        set the values of the date coordinate in the output Dataset. If not, it
        will be computed from the time index values in `source` (which may be
        slightly slower). By default `Era5LandLinearizedTimeDimId.DATE`. **NB!**
        If not specified, the computed values will be stored in a new coordinate
        variable with the default name, and any existing variable with the same
        name will be overwritten.
    source_step_coord : str|None, optional
        The name of the coordinate variable for the step component of the time
        dimension in the source Dataset, if present. Same comments apply as for
        `source_date_coord`. By default `Era5LandLinearizedTimeDimId.STEP`.
        **NB!** If not specified, the computed values will be stored in a new
        coordinate variable with the default name, and any existing variable
        with the same name will be overwritten.
    target_time_coord_name : str | None, optional
        The name to use for a 2D time coordinate variable in the output Dataset.
        If this variable is already present in the source Dataset, it will be
        used as the target variable without checking its values. If not present,
        the time index coordinate of `source` will simply be renamed to this
        name. If None, no 2D time coordinate variable will be set, and the
        time index coordinate of `source` will be dropped. By default
        `Era5LandCoord.TIME_LINEAR`.
    """
    if source_date_coord is None:
        source_date_coord = Era5LandLinearizedTimeDimId.DATE
        # Get the date from the time index coordinate. Since the last step for
        # each date is midnight and therefore belongs to the next date, we need
        # to subtract 1 second before flooring to date.
        source = source.assign_coords(
            {
                source_date_coord: (
                    source[source_time_dim] - np.timedelta64(1, 's')
                ).dt.floor('D')
            },
        )
    if source_step_coord is None:
        source_step_coord = Era5LandLinearizedTimeDimId.STEP
        # Get the step from the time index coordinate by taking the difference
        # between the time and the date
        source = source.assign_coords(
            {
                source_step_coord: (
                    source[source_time_dim] - source[source_date_coord]
                )
            },
        )
    if fast_unstack:
        target_ds: xr.Dataset = (
            source
            .coarsen({source_time_dim: 24}, boundary='exact')
            .construct(
                {
                    source_time_dim: (source_date_coord, source_step_coord),
                }
            )
        )
    else:
        target_ds: xr.Dataset = (
            source
            .set_index(
                {
                    source_time_dim: (source_date_coord, source_step_coord),
                },
            )
            .unstack(source_time_dim)
        )
    if (
            (target_time_coord_name is None)
            or (target_time_coord_name in target_ds.variables)
    ):
        target_ds = target_ds.drop_vars(source_time_dim)
    else:
        target_ds = target_ds.rename({source_time_dim: target_time_coord_name})
    target_ds = target_ds.rename(
        {
            source_date_coord: target_date_dim,
            source_step_coord: target_step_dim,
        }
    )
    return target_ds
###END def era5land_from_linear_time


def postprocess_converted_datm_ds(
        target_ds: xr.Dataset,
        *,
        target_stream: Datm7Stream,
        source: xr.Dataset,
        dim_order: Iterable[Datm7Dim] = datm7_dim_order,
) -> xr.Dataset:
    """Postprocesses a converted DATM xarray Dataset after all variables have
    been created, to perform any final adjustments needed.

    Currently the following postprocessing steps are performed in this function:
    - Reorder the dimensions to match the expected order in the DATM7 data
      (CESM may crash otherwise).
    - Sort coordinates in ascending order by and for each dimension, as expected
      by CESM.

    Data type conversion, float precision, fill values and time units are set in
    the `convert_files` module, since these attributes are specific to the
    output file, and not necessary if the output is used for further processing
    in Python.

    Parameters
    ----------
    target_ds : xr.Dataset
        The converted DATM xarray Dataset.
    target_stream : Datm7Stream
        The target DATM7 stream.
    source : xr.Dataset
        The source ERA5 Land Dataset.
    dim_order : Iterable[Datm7Dim], optional
        The desired order of dimensions in the output Dataset. By default given
        by the module-level attribute `datm7_dim_order`.

    Returns
    -------
    xr.Dataset
        The postprocessed DATM xarray Dataset. This will be a new Dataset
        object. The original `target_ds` will not be modified.
    """
    target_ds = target_ds.copy(deep=False)
    use_dim_order: list[str] = [str(_dim) for _dim in dim_order]
    if list(target_ds.dims) != use_dim_order:
        logger.debug(
            f'Reordering dimensions from {target_ds.dims} to {use_dim_order} '
            f'for target stream {target_stream}...'
        )
        target_ds = target_ds.transpose(*use_dim_order)
        for _var in target_ds.data_vars:
            target_ds[_var] = target_ds[_var].transpose(*use_dim_order)
    logger.debug(
        'Sorting coordinates in ascending order for target stream '
        f'{target_stream}...'
    )
    target_ds = target_ds.sortby(use_dim_order, ascending=True)
    return target_ds
###END def postprocess_converted_datm_ds


def make_target_var(
        target_var: Datm7Var,
        source: xr.Dataset,
) -> xr.DataArray:
    """Creates a target DATM7 variable DataArray from the source ERA5 Land
    Dataset.

    **NB!** This function normally assumes that cumulative ERA5 Land variables
    for radiation and precipitation have already been decumulated using
    `decumulate_era5land_var`. It uses the functions in the module-level
    dictionary `value_conversion_funcs` to perform the conversion, and the
    functions for radiation and precipitation variables in that dictionary
    assume that cumulative variables have been decumulated by a simple diff,
    meaning that the value at each time step gives the cumulated value for the
    previous time step only.

    Parameters
    ----------
    target_var : Datm7Var
        The target DATM7 variable to create.
    source : xr.Dataset
        The source ERA5 Land Dataset.
    """
    value_arr: xr.DataArray = value_conversion_funcs[target_var](source)
    new_arr = set_target_dims_and_coords(
        arr=value_arr,
        target_var=target_var,
        source=source,
    )
    new_arr: xr.DataArray = add_target_var_attrs(
        arr=new_arr,
        target_var=target_var,
        source=source,
    )
    return new_arr
###END def make_target_var


def compute_average_rate(
        source: xr.DataArray,
        *,
        time_step_dim: str = ERA5_LINEARIZED_TIME_DIM,
        time_step_seconds: float|int|xr.DataArray|None = None,
) -> xr.DataArray:
    """Computes average rate values from decumulated ERA5 Land variable values.

    The function assumes that the value at each point is the cumulated value for
    the previous time step only. The average rate is therefore computed by
    taking the average of the value at each point and the value at the next
    point along the time dimension, and dividing by 2 times the time step
    length.

    Parameters
    ----------
    source : xr.DataArray
        The source ERA5 Land variable DataArray.
    time_step_dim : str, optional
        The name of the time step dimension. By default given by the
        module-level attribute `ERA5_LINEARIZED_TIME_DIM`. If
        `time_step_seconds` is not specified, the dimnension must have an
        index (coordinate variable), and it must be of dtype `numpy.datetime64`.
    time_step_seconds : float|int|xr.DataArray|None, optional
        The time step length(s) in seconds. If None, the time step length will
        be inferred from the time coordinate values along the time step
        dimension. By default None.

    Returns
    -------
    xr.DataArray
        The average rate DataArray. The values will have the same units as
        `source` in the numerator, divided by seconds in the denominator. You
        will need to convert manually if you want other units than per second.
    """
    if time_step_seconds is None:
        time_coords: xr.DataArray = source[time_step_dim]
        time_step_seconds = time_coords.diff(
            dim=time_step_dim,
            label='upper',
        ).dt.total_seconds()
    shifted_source: xr.DataArray = source.shift({time_step_dim: -1})
    if isinstance(time_step_seconds, xr.DataArray):
        shifted_time_step_seconds: xr.DataArray|float|int = (
            time_step_seconds
            .shift({time_step_dim: -1})
        )
    else:
        shifted_time_step_seconds = time_step_seconds
    return (
        (source + shifted_source)
        / (time_step_seconds + shifted_time_step_seconds)
    )
###END def compute_average_rate

value_conversion_funcs: dict[Datm7Var, Callable[[xr.Dataset], xr.DataArray]] = {
    Datm7Var.TBOT: lambda source: \
        source[era5land_grib_varnames[Era5LandVar.T2M]],
    Datm7Var.PSRF: lambda source: \
        source[era5land_grib_varnames[Era5LandVar.SP]],
    Datm7Var.QBOT: lambda source: \
        specific_humidity_from_dewpoint_pressure(
            pressure=source[era5land_grib_varnames[Era5LandVar.SP]],
            dewpoint=source[era5land_grib_varnames[Era5LandVar.D2M]],
        ),
    Datm7Var.WIND: lambda source: \
        xr.ufuncs.sqrt(
            source[era5land_grib_varnames[Era5LandVar.U10]]**2 +
            source[era5land_grib_varnames[Era5LandVar.V10]]**2
        ),
    Datm7Var.FSDS: lambda source: compute_average_rate(
        source[era5land_grib_varnames[Era5LandVar.SSRD]],
    ),
    Datm7Var.PRECTmms: lambda source: compute_average_rate(
        source[era5land_grib_varnames[Era5LandVar.TP]],
    ) * 1000.0,  # Multiply by 1000, to convert from cumulative m to mm/sec rate.
    Datm7Var.FLDS: lambda source: compute_average_rate(
        source[era5land_grib_varnames[Era5LandVar.STRD]],
    ),
}


def add_target_var_attrs(
        arr: xr.DataArray,
        *,
        target_var: Datm7Var,
        source: xr.Dataset,
) -> xr.DataArray:
    """A utilitity function that adds attributes to a target DATM7 variable
    after converting/computing the raw values.

    Parameters
    ----------
    arr : xr.DataArray
        The target DATM7 variable DataArray to add attributes to.
    target_var : Datm7Var
        The target DATM7 variable id.
    source : xr.Dataset
        The source ERA5 Land Dataset. This parameter is not used in the current
        implememntation, but is included for possible future use cases where
        source-specific attributes may need to be added. If you need to use the
        function without a source Dataset, you can at present pass an empty
        Dataset.

    Returns
    -------
    xarray.DataArray
        The target DATM7 variable DataArray with attributes added. Note that a
        new DataArray is returned, the input array is not modified in place,
        but the underlying data is not copied.
    """
    return arr.drop_attrs().assign_attrs(
        **datm7_var_attrs[target_var],
    )
###END def add_target_var_attrs


def set_target_dims_and_coords(
        arr: xr.DataArray,
        *,
        target_var: Datm7Var,
        source: xr.Dataset,
        dim_map: Mapping[Era5LandDim|LinearizedTimeDimId, Datm7Dim] \
            = era5_to_datm7_dim_map,
) -> xr.DataArray:
    """Set the dimensions and coordinate variables for a target DATM7 variable
    after converting/computing the raw values.

    The function drops all existing coordinates other than the index (dimension)
    coordinates, and renames the dimensions and index coordinates to the names
    used in DATM7.

    Note that the function assumes that the source Dataset has already been
    converted to linear time using `era5land_to_linear_time`.

    Parameters
    ----------
    arr : xr.DataArray
        The target DATM7 variable DataArray to set dimensions and coordinates
        for.
    target_var : Datm7Var
        The target DATM7 variable id.
    source : xr.Dataset
        The source ERA5 Land Dataset. This parameter is not used in the current
        implememntation, but is included for possible future use cases where
        source-specific attributes may need to be added. If you need to use the
        function without a source Dataset, you can at present pass an empty
        Dataset.
    dim_map : Mapping[Era5LandDim|LinearizedTimeDimId, Datm7Dim], optional
        Mapping of ERA5 Land dimensions (or linearized time dimension) to
        corresponding DATM7 dimensions. By default the module-level
        `era5_to_datm7_dim_map` dictionary. Note that the time dimension must
        have been linearized already, so that a single time dimension in `arr`
        source maps to the DATM7 time dimension.

    Returns
    -------
    xarray.DataArray
        A new xarray.DataArray with the same values as `arr`, but with
        dimensions and coordinates appropriate for the target DATM7 variable.
    """
    rename_mapping: dict[Hashable, Hashable] = {
        str(_source_dim): str(_target_dim)
        for _source_dim, _target_dim in dim_map.items()
        if _source_dim in arr.dims
    }
    return arr.reset_coords(drop=True).rename(rename_mapping)
###END def set_target_dims_and_coords


class UnexpectedEra5LandUnitsError(Exception):
    """Raised when one or more ERA5 Land variable has unexpected units.

    Attributes
    ----------
    unexpected_vars : dict[Era5LandVar|Era5LandCoord, str]
        The ERA5 Land variables and coordinates that had unexpected units,
        mapped to the unepected units found.
    expected_units : dict[Era5LandVar|Era5LandCoord, str]
        Mapping from the same ERA5 Land variables as in `unexpected_vars`, to
        the expected units.

    Methods
    -------
    default_error_msg() -> str
        Creates a default error message listing the variables with unexpected
        units, the units found, and the expected units.
    """

    def __init__(
            self,
            *args,
            unexpected_vars: dict[Era5LandVar|Era5LandCoord, str],
            expected_units: dict[Era5LandVar|Era5LandCoord, str|None]|None = None,
            **kwargs,
    ) -> None:
        """
        Parameters
        ----------
        *args
            Positional arguments to pass to the base Exception class, which
            can include a custom error message as its first positional argument.
            If no positional arguments are given, a default error message will
            be created using the `default_error_msg` method.
        unexpected_vars : dict[Era5LandVar|Era5LandCoord, str]
            The ERA5 Land variables and coordinates that had unexpected units,
            mapped to the unepected units found.
        expected_units : dict[Era5LandVar|Era5LandCoord, str|None]|None, optional
            Mapping from the same ERA5 Land variables as in `unexpected_vars`,
            to the expected units. If None, the expected units will be looked up
            from the module-level `era5_var_units` dictionary. By default None.
        """
        self.unexpected_vars: dict[Era5LandVar|Era5LandCoord, str] \
            = unexpected_vars
        self.expected_units: dict[Era5LandVar|Era5LandCoord, str|None] \
            = expected_units if expected_units is not None else {
                _var: era5_var_units.get(_var, None)
                for _var in unexpected_vars.keys()
            }
        if len(args) == 0:
            args = (self.default_error_msg(),)
        super().__init__(*args, **kwargs)
    ###END def UnexpectedEra5LandUnitsError.__init__

    def default_error_msg(self) -> str:
        """Creates a default error message listing the variables with unexpected
        units, the units found, and the expected units.
        """
        msg_lines = [
            'The following ERA5 Land variables have unexpected units:'
        ]
        for _var, _found_units in self.unexpected_vars.items():
            _expected_units = self.expected_units[_var]
            msg_lines.append(
                f' - {_var.value}: found "{_found_units}", expected '
                f'"{_expected_units}"'
            )
        return '\n'.join(msg_lines)
    ###END def UnexpectedEra5LandUnitsError.default_error_msg

###END class UnexpectedEra5LandUnitsError


def check_era5land_units(
        source: xr.Dataset,
        *,
        variables: Iterable[Era5LandVar|Era5LandCoord],
        expected_units: dict[Era5LandVar|Era5LandCoord, str|None] | None \
            = None,
        raise_error: bool = True,
) -> dict[Era5LandVar|Era5LandCoord, str]:
    """Checks that the specified ERA5 Land variables and coordinates in the
    source Dataset have the expected units.

    Parameters
    ----------
    source : xr.Dataset
        The source ERA5 Land Dataset.
    variables : Iterable[Era5LandVar|Era5LandCoord]
        The ERA5 Land variables and coordinates to check.
    expected_units : dict[Era5LandVar|Era5LandCoord, str], optional
        Mapping from ERA5 Land variables and coordinates to their expected
        units. By default the module-level `era5_var_units` dictionary.
    raise_error : bool, optional
        Whether to raise an UnexpectedEra5LandUnitsError if any variables have
        unexpected units. If False, the function will simply return a mapping
        of the variables with unexpected units to the units found. By default
        True.
    """
    unexpected_vars: dict[Era5LandVar|Era5LandCoord, str] = {}
    if expected_units is None:
        expected_units = {
            _var: era5_var_units.get(_var, None)
            for _var in variables
        }
    for _var in variables:
        _varname = (
            era5land_grib_varnames[_var]
            if isinstance(_var, Era5LandVar) else _var.value
        )
        _found_units: str = source[_varname].attrs.get('units', '')
        _expected_units = expected_units.get(_var, '')
        if _found_units != _expected_units:
            unexpected_vars[_var] = _found_units
    if unexpected_vars and raise_error:
        raise UnexpectedEra5LandUnitsError(
            unexpected_vars=unexpected_vars,
            expected_units=expected_units,
        )
    return unexpected_vars
###END def check_era5land_units


_case_title_target_stream_defaults: dict[Datm7Stream, str] = {
    Datm7Stream.PREC: 'Precipitation',
    Datm7Stream.SOLR: 'Downward Shortwave Radiation',
    Datm7Stream.TPQWL: 'Temperature, Pressure, Winds, Humidity, and Downward ' \
        'Longwave Radiation',
}

def make_datm_base(
        *,
        source: xr.Dataset,
        target_stream: Datm7Stream,
        creation_date: str | datetime.date | None = None,
        case_title: str | None = None,
        other_attrs: dict[str, str] | None = None,
) -> xr.Dataset:
    """Creates the base DATM xarray Dataset with coordinates, dimensions, and
    Dataset attributes set.

    Parameters
    ----------
    source : xr.Dataset
        The source ERA5 Land Dataset.
    target_stream : Datm7Stream
        The target DATM7 stream to create.
    creation_date : str | datetime.date | None, optional
        The creation date to set in the DATM Dataset attributes. If a string,
        it should be in the format 'YYYY-MM-DD'. If a `datetime.date` object is
        given, it will be converted to a string in the same format. If None,
        the current date will be used. By default None.
    case_title : str | None, optional
        The case title to set in the DATM Dataset attributes. If None, an
        appropriate default will be generated based on the target stream.
        By default None.
    other_attrs : dict[str, str] | None, optional
        Additional Dataset attributes to set in the DATM Dataset. By default
        None.

    Returns
    -------
    xr.Dataset
    """
    target_stream = Datm7Stream(target_stream)

    creation_date_obj: datetime.date
    if creation_date is None:
        creation_date_obj = datetime.date.today()
    elif isinstance(creation_date, str):
        creation_date_obj = datetime.date.fromisoformat(creation_date)
    elif isinstance(creation_date, datetime.date):
        creation_date_obj = creation_date
    else:
        raise TypeError(
            'creation_date must be a str in format "YYYY-MM-DD", a '
            'datetime.date object, or None.'
        )

    if case_title is None:
        case_title = (
            'ERA5 Land 1-Hourly Atmospheric Forcing: '
            f'{_case_title_target_stream_defaults[target_stream]}.'
        )

    if other_attrs is None:
        other_attrs = {}

    target_ds = xr.Dataset(
        coords={
            # NB! The time dimension must be the outermost dimension in order
            # not to cause potentially silent crashes in the PIO library, and
            # possibly elsewhere in DATM (at least as of January 2026).
            str(Datm7Coord.TIME): xr.DataArray(
                data=source[ERA5_LINEARIZED_TIME_DIM].data,
                dims=(Datm7Dim.TIME.value,),
                attrs=datm7_var_attrs[Datm7Coord.TIME],
            ),
            str(Datm7Coord.LAT): xr.DataArray(
                data=source[Era5LandDim.LAT].data,
                dims=(Datm7Dim.LAT.value,),
                attrs=datm7_var_attrs[Datm7Coord.LAT],
            ),
            str(Datm7Coord.LON): xr.DataArray(
                data=source[Era5LandDim.LON].data,
                dims=(Datm7Dim.LON.value,),
                attrs=datm7_var_attrs[Datm7Coord.LON],
            ),
        },
        attrs={
            'case_title': case_title,
            'conventions': '',
            'creation_date': creation_date_obj.isoformat(),
            **other_attrs,
        }
    )

    # Create the 2-dimensions lat/lon variables LATIXY and LONGXY by
    # broadcasting (should be constant along the other dimension).
    for _var, _coord in (
            (Datm7Coord.LATIXY, Datm7Coord.LAT),
            (Datm7Coord.LONGXY, Datm7Coord.LON),
    ):
        target_ds[str(_var)] = target_ds[str(_coord)].broadcast_like(
            target_ds[[str(Datm7Coord.LAT), str(Datm7Coord.LON)]]
        )

    return target_ds
###END def make_datm_base


def round_coords[_XrObj: xr.Dataset | xr.DataArray](
        xr_obj: _XrObj,
        round_to: Mapping[Hashable, float] | float,
        *,
        use_numpy_digit_number_rounding: bool = False,
) -> _XrObj:
    """Round coordinate values in an xarray Dataset.

    NB! By default, this function rounds simply by dividing by the number being
    rounded to, rounding to the nearest integer with `numpy.round`, and then
    multiplying back. This means that for values close to the middle of a
    rounding interval, the results may not be as expected due to limited float
    precision. If you want to round to a specific number of digits with
    predictable outcomes instead, you can set `use_numpy_digit_number_rounding`
    to True. See parameter description below.

    Parameters
    ----------
    xr_obj : xarray.Dataset or xarray.DataArray
        The xarray object to round coordinates in. The function will return a
        new object with the same data variables and values, but with rounded
        coordinate values.
    round_to : Mapping[str, float] | float
        The number to round to. The coordinates will be rounded to the nearest
        multiple of tis number. If a mapping, the keys should be coordinate
        names and the values the rounding intervals to round those coordinates
        to. If a single float or int is given, it will be used as the rounding
        interval for all coordinates in `xr_obj.coords`. Note that the function
        will also round data variables if `round_to` is a mapping and contains
        keys that are names of data variables in `xr_obj`. No warning or error
        will be raised in this case. Data variables will not be rounded if
        `round_to` is a single number.
    use_numpy_digit_number_rounding : bool, optional
        Whether to round by using `numpy.round` with a specific number of digits
        instead of by dividing and multiplying. If True, the round_to will be
        interpreted as the number of digits after the decimal point to round to,
        rather than as the actual number to round to. See the documentation of
        `numpy.round` for more details on how the rounding works in this case.

    Returns
    -------
    xarray.Dataset or xarray.DataArray
        A new xarray object with the same data variables and values as `xr_obj`,
        but with rounded coordinate values.

    Raises
    ------
    KeyError
        If `round_to` is a mapping and contains keys that are not names of
        coordinates or data variables in `xr_obj`.
    TypeError
        If `use_numpy_digit_number_rounding` is True and any of the rounding
        values in `round_to` are not integers.
    """
    if isinstance(round_to, (float, int)):
        round_to_mapping: Mapping[Hashable, float] = {
            coord: float(round_to) for coord in xr_obj.coords
        }
    else:
        if any(
                _var_name not in xr_obj.variables
                for _var_name in round_to.keys()
        ):
            raise KeyError(
                'All keys in round_to mapping must be names of coordinates or '
                'data variables in xr_obj.'
            )
        round_to_mapping = round_to
    rounded_obj: _XrObj = xr_obj.copy(deep=False)
    if use_numpy_digit_number_rounding:
        if any(
                int(_rounding) != _rounding
                for _rounding in round_to_mapping.values()
        ):
            raise TypeError(
                'When use_numpy_digit_number_rounding is True, the values in '
                'round_to must be integers representing the number of digits to '
                'round to.'
            )
        rounding_func: Callable[[xr.DataArray, float], xr.DataArray] = (
            lambda _arr, _decimals: _arr.round(int(_decimals))
        )
    else:
        rounding_func: Callable[[xr.DataArray, float], xr.DataArray] = (
            lambda _arr, _rounding: (
                (_arr / _rounding).round() * _rounding
            )
        )
    for _var_name, _rounding in round_to_mapping.items():
        rounded_obj[_var_name] = rounding_func(
            rounded_obj[_var_name],
            _rounding,
        )
    return rounded_obj
###END def round_coords
