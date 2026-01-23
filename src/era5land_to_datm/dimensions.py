"""Definitions of dimension names used in ERA5-Land to DATM conversion.

Enums
-----
Era5Dim
    String enum for ERA5 Land dimension ids.
ERA5LandTimeLayout
    Enumeration of possibble time layouts in ERA5 land data (linear, or date
    plus intradate timediff step).
Datm7Dim
    String enum for DATM dimension ids.

Attributes
----------
ERA5_LINEARIZED_TIME_DIM: Final[str]
    The name to use for the time dimension in ERA5 Land Dataset instances after
    linearizing.
"""
import enum
import typing as tp



class Era5LandDim(enum.StrEnum):
    """String enum for ERA5 Land dimension ids."""

    DATE = 'time'
    STEP = 'step'
    LAT = 'latitude'
    LON = 'longitude'

###END class Era5Dim


class ERA5LandTimeLayout(enum.StrEnum):
    """Enumeration of possibble time layouts in ERA5 land data (linear, or date
    plus intradate timediff step).
    """

    LINEAR = 'linear'
    DATE_STEP = 'date_step'

###END class ERA5LandTimeLayout


class Datm7Dim(enum.StrEnum):
    """String enum for DATM dimension ids."""

    TIME = 'time'
    LAT = 'lat'
    LON = 'lon'

###END class Datm7Dim


LinearizedTimeDimId = tp.NewType('LinearizedTimeDimId', str)
ERA5_LINEARIZED_TIME_DIM: tp.Final[LinearizedTimeDimId] = LinearizedTimeDimId(
    Datm7Dim.TIME
)
