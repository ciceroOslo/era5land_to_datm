"""Custom types used in era5land_to_datm package.

Types
-----
YearMonth
    A namedtuple representing a year and month (e.g., (2020, 1) for January
    2020).
VarSet
    A frozenset of Era5LandVar instances representing a set of variables.

Enums
-----
Era5LandVar
    Enumeration of ERA5-Land variables available for download.
"""
import enum
import typing as tp



class YearMonth(tp.NamedTuple):
    """A namedtuple representing a year and month.

    Attributes
    ----------
    year : int
        The year (e.g., 2020).
    month : int
        The month (1-12).
    """
    year: int
    month: int
###END class YearMonth


class Era5LandVar(enum.StrEnum):
    """Enumeration of ERA5-Land variables available for download."""

    D2M = '2m_dewpoint_temperature'
    SP = 'surface_pressure'
    T2M = '2m_temperature'
    TP = 'total_precipitation'
    U10 = '10m_u_component_of_wind'
    V10 = '10m_v_component_of_wind'
    SSRD = 'surface_solar_radiation_downwards'
    STRD = 'surface_thermal_radiation_downwards'

###END class Era5LandVar


VarSet = frozenset[Era5LandVar]

