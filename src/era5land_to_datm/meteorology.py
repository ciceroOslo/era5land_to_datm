"""Functions for metorological computations.

This module is currently only used for the function
`specific_humidity_from_dewpoint_pressure`, which computes specific humidity
from dewpoint temperature and pressure. DATM requires a data field for specific
humidity, but ERA5-Land provides dewpoint temperature instead.

Functions
---------
specific_humidity_from_dewpoint_pressure
    Computes specific humidity from dewpoint and pressure DataArrays. The
    function uses the method in Ambaum (2020) to compute the partial pressure of
    water vapour from the dewpoint, and then computes specific humidity from
    there using pressure data and the relative molecular weight of water and dry
    air.

Attributes
----------
T0 : float
    Temperature at the triple point of water in Kelvin. Used as reference
    temperature in the computation of specific humidity.
e0 : float
    Saturation vapour pressure at the triple point of water in Pascals. Used as
    reference vapour pressure in the computation of specific humidity.
Rv : float
    Specific gas constant for water vapour in J kg^-1 K^-1. Used in the
    computation of specific humidity. The value used is 462.52 J kg^-1 K^-1.
    Taken from Ambaum (2020).
dCp : float
    Difference in specific heat capacity of liquid water minus water vapour at
    constant pressure at the triple point of water, in J kg^-1 K^-1. Used in the
    computation of specific humidity. The value used is 212 J kg^-1 K^-1. Taken
    from Ambaum (2020).

References
----------
Ambaum, M. H. P. (2020). Accurate, simple equation for saturated vapour pressure
over water and ice. Quarterly Journal of the Royal Meteorological Society,
146(733), 4252-4258. https://doi.org/10.1002/qj.3899
"""