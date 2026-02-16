"""Definitions of dimension names used in ERA5-Land to DATM conversion.

Enums
-----
Era5LandDim
    String enum for ERA5 Land dimension ids.
Era5LandLinearizedTimeDimId
    String enum for the names to use for the time dimension in the ERA5 Land
    datasets after linearizing the time dimension, and for the coordinates that
    represent the values of the original date and (intra-date) step dimensions.
ERA5LandTimeLayout
    Enumeration of possibble time layouts in ERA5 land data (linear, or date
    plus intradate timediff step).
Datm7Dim
    String enum for DATM dimension ids.
Era5

Attributes
----------
ERA5_LINEARIZED_TIME_DIM: Final[str]
    The name to use for the time dimension in ERA5 Land Dataset instances after
    linearizing. This is equal to `Era5LandLinearizedTimeDimId.TIME`, and should
    usually be the same as `Datm7Dim.TIME` and be used for the time dimension in
    the converted DATM datasets.
datm7_dim_order: Final[tuple[Datm7Dim, Datm7Dim, Datm7Dim]]
    The order of dimensions to use in the converted DATM datasets.
era5_to_datm7_dim_map : dict[Era5LandDim, Datm7Dim]
    Mapping of ERA5 Land dimensions to corresponding DATM7 dimensions, after
    linearizing the ERA5 Land time dimension with
    `convert_data.era5land_to_linear_time`.
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


class Era5LandLinearizedTimeDimId(enum.StrEnum):
    """String enum for the names to use for the time dimension in the ERA5 Land
    datasets after linearizing the time dimension, and for the coordinates that
    represent the values of the original date and (intra-date) step dimensions.

    Members
    -------
    TIME
        The name to use for the time dimension in ERA5 Land Dataset instances
        after linearizing. This is equal to `ERA5_LINEARIZED_TIME_DIM`, and
        should usually be the same as `Datm7Dim.TIME` and be used for the time
        dimension in the converted DATM datasets.
    DATE
        The name to use for the coordinate representing the original date values.
    STEP
        The name to use for the coordinate representing the original intradate
        step values.
    """
    TIME = 'time'
    DATE = 'date'
    STEP = 'step'
###END class Era5LandLinearizedTimeDimId


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


datm7_dim_order: tp.Final[tuple[Datm7Dim, Datm7Dim, Datm7Dim]] = (
    Datm7Dim.TIME,
    Datm7Dim.LAT,
    Datm7Dim.LON,
)


LinearizedTimeDimId = tp.NewType('LinearizedTimeDimId', str)
ERA5_LINEARIZED_TIME_DIM: tp.Final[LinearizedTimeDimId] = LinearizedTimeDimId(
    Datm7Dim.TIME
)
try:
    assert ERA5_LINEARIZED_TIME_DIM == Era5LandLinearizedTimeDimId.TIME
except AssertionError as _ass_err:
    raise AssertionError(
        f'ERA5_LINEARIZED_TIME_DIM ({ERA5_LINEARIZED_TIME_DIM!r}) must be '
        'equal to Era5LandLinearizedTimeDimId.TIME '
        f'({Era5LandLinearizedTimeDimId.TIME!r})'
    ) from _ass_err


era5_to_datm7_dim_map: dict[Era5LandDim|LinearizedTimeDimId, Datm7Dim] = {
    Era5LandDim.LAT: Datm7Dim.LAT,
    Era5LandDim.LON: Datm7Dim.LON,
    ERA5_LINEARIZED_TIME_DIM: Datm7Dim.TIME,
}
