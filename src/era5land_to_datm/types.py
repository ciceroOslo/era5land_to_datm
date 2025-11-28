"""Custom types used in era5land_to_datm package.

Types
-----
YearMonth
    A namedtuple representing a year and month (e.g., (2020, 1) for January
    2020).
VarSet
    A frozenset of Era5LandVar instances representing a set of variables.
"""
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
