"""Functions for metorological computations.

This module is currently only used for the function
`specific_humidity_from_dewpoint_pressure`, which computes specific humidity
from dewpoint temperature and pressure.

Functions
---------
specific_humidity_from_dewpoint_pressure
    Computes specific humidity from dewpoint and pressure DataArrays. The
    function simply calls the function
    `metpy.calc.specific_humidity_from_dewpoint` from the `metpy` package, and
    lightly adapts the DataArray returned from that function.
"""
import typing as tp

from metpy.calc import specific_humidity_from_dewpoint
import xarray as xr



def specific_humidity_from_dewpoint_pressure(
        *,
        pressure: xr.DataArray,
        dewpoint: xr.DataArray,
        **kwargs,
) -> xr.DataArray:
    """Compute specific humidity from dewpoint and pressure DataArrays.

    The function in its current implementation calls the function
    `metpy.calc.specific_humidity_from_dewpoint` from the `metpy` package. It
    requires that `presssure` and `dewpoint` DataArrays either have `units` as
    an attribute, or use the `pint` package to specify proper units.

    The return value is dimensionless (kg of water vapor per kg of moist air).

    Parameters
    ----------
    pressure : xr.DataArray
        The pressure DataArray. Must have the `units` attribute set or use the
        `pint` or `metpy` package to specify proper units.
    dewpoint : xr.DataArray
        The dewpoint DataArray. Must have the `units` attribute set or use the
        `pint` or `metpy` package to specify proper units.
    **kwargs
        Additional keyword arguments passed to the underlying calculation
        function, currently `metpy.calc.specific_humidity_from_dewpoint`. Note
        that this function by default calls that function with the `phase`
        parameter set to `'auto'` (as opposed to the default value of
        `'liquid'`). This ensures that the calculation tries to interpolate
        appropriately between formulas for humidity over liquid water and over
        ice. The calculation can be sped up (by roughly a factor of 2) by
        setting `phase` to either `'liquid'` or `'solid'` if the user knows that
        only one of those regimes is relevant.

    Returns
    -------
    xr.DataArray
        The specific humidity DataArray. The DataArray is equal to what is
        returned by the `metpy.calc.specific_humidity_from_dewpoint` function,
        but dequantified using the `.metpy.dequantify` accessor function to
        convert the underlying data array from a `pint` object to a regular
        `numpy.ndarray`.
    """
    kwargs = {'phase': 'auto'} | kwargs
    return (
        specific_humidity_from_dewpoint(
            pressure=pressure,
            dewpoint=dewpoint,
            **kwargs,
        )
        .metpy.dequantify()
    )
###END def specific_humidity_from_dewpoint_pressure
