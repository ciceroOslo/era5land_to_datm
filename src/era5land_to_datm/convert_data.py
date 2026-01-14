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
decumulate_era5land_var
    Differences cumulative ERA5 Land variables along the intra-day step
    dimension, to produce a DataArray where each value gives the cumulated value
    for the previous time step only. The value of each time step can then be
    averaged with the next time step and divided by 2 times the time step length
    to produce an average rate value.

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
    Iterable,
)
import datetime
import functools

import xarray as xr

from .datm_streams import (
    Datm7Stream,
    datm7_stream_variables,
)
from .dimensions import (
    Datm7Dim,
    ERA5_LINEARIZED_TIME_DIM,
    Era5LandDim,
    ERA5LandTimeLayout,
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

    check_era5land_units(
        source=source,
        variables=required_era5_vars,
    )

    target_ds: xr.Dataset = make_datm_base(
        source=source,
        target_stream=target_stream,
    )

    cumulative_required_vars: frozenset[Era5LandVar] = (
        required_era5_vars & era5_cumulative_vars
    )
    source_decumulated: xr.Dataset = source.copy(deep=False)
    for _cum_var in cumulative_required_vars:
        _cum_varname = era5land_grib_varnames[_cum_var]
        source_decumulated[_cum_varname] = decumulate_era5land_var(
            source[_cum_varname],
        )

    if time_layout == ERA5LandTimeLayout.DATE_STEP:
        source_1d_time: xr.Dataset = era5land_to_linear_time(
            source=source_decumulated,
        )
    else:
        source_1d_time: xr.Dataset = source_decumulated
    del source_decumulated
    for _target_var in target_vars:
        target_ds[_target_var.value] = make_target_var(
            target_var=_target_var,
            source=source_1d_time,
        )
    del source_1d_time
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
    add_target_var_attrs(
        arr=value_arr,
        target_var=target_var,
        source=source,
    )
    set_target_dims_and_coords(
        arr=value_arr,
        target_var=target_var,
        source=source,
    )
    return value_arr
###END def make_target_var

value_conversion_funcs: dict[Datm7Var, Callable[[xr.Dataset], xr.DataArray]] = {
    Datm7Var.TBOT: lambda source: \
        source[era5land_grib_varnames[Era5LandVar.T2M]],
    Datm7Var.PSRF: lambda source: \
        source[era5land_grib_varnames[Era5LandVar.SP]],
    Datm7Var.QBOT: lambda source: \
        compute_specific_humidity(
            temperature=source[era5land_grib_varnames[Era5LandVar.T2M]],
            pressure=source[era5land_grib_varnames[Era5LandVar.SP]],
            dewpoint=source[era5land_grib_varnames[Era5LandVar.D2M]],
        ),
    Datm7Var.WIND: lambda source: \
        xr.ufuncs.sqrt(
            source[era5land_grib_varnames[Era5LandVar.U10]]**2 +
            source[era5land_grib_varnames[Era5LandVar.V10]]**2
        ),
} | {
    _target_var: lambda source: compute_average_rate(
        source[era5land_grib_varnames[_source_var]],
    ) for _target_var, _source_var in (
        (Datm7Var.FSDS, Era5LandVar.SSRD),
        (Datm7Var.PRECTmms, Era5LandVar.TP),
        (Datm7Var.FLDS, Era5LandVar.STRD),
    )
}


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
    if isinstance(creation_date, str):
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
            # TODO: **NB!** THE ERA5 TIME COORDINATES MUST BE CONVERTED TO `cftime.DatetimeNoLeap`!
            str(Datm7Coord.TIME): xr.DataArray(
                data=source[ERA5_LINEARIZED_TIME_DIM].data,
                dims=(Datm7Dim.TIME.value,),
                attrs=datm7_var_attrs[Datm7Coord.TIME],
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
    # broadcasting (should be cvonstant along the other dimension).
    for _var, _coord in (
            (Datm7Coord.LATIXY, Datm7Coord.LAT),
            (Datm7Coord.LONGXY, Datm7Coord.LON),
    ):
        target_ds[str(_var)] = target_ds[str(_coord)].broadcast_like(
            target_ds[[str(Datm7Coord.LAT), str(Datm7Coord.LON)]]
        )

    return target_ds
###END def make_datm_base
