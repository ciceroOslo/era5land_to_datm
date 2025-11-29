"""Code for specifying regions for ECWMF data requests."""
import typing as tp

import pydantic

type LatitudeFloat = tp.Annotated[
    float,
    pydantic.Field(ge=-90.0, le=90.0),
]
type LongitudeFloat = tp.Annotated[
    float,
    pydantic.Field(ge=-360.0, le=360.0),
]


class EcmwfArea(pydantic.BaseModel):
    """Pydantic model representing a geographical area for ECMWF data requests.

    Attributes
    ----------
    north : LatitudeFloat
        The northern bound in degrees (-90.0 to 90.0).
    west : LongitudeFloat
        The western bound in degrees (-360.0 to 360.0).
    south : LatitudeFloat
        The southern bound in degrees (-90.0 to 90.0).
    east : LongitudeFloat
        The eastern bound in degrees (-360.0 to 360.0).

    Validators
    ----------
    validate_north_south
        Validate that the south bound is less than the north bound.
    validate_east_west
        Validate that the east bound is greater than the west bound.

    Methods
    -------
    to_tuple
        Convert the EcmwfArea instance to an EcmwfAreaTuple.
    """

    north: LatitudeFloat
    west: LongitudeFloat
    south: LatitudeFloat
    east: LongitudeFloat

    @pydantic.model_validator(mode='after')
    def validate_north_south(
        self,
    ) -> tp.Self:
        """Validate that the south bound is less than the north bound."""
        if self.south >= self.north:
            raise ValueError(
                'The south bound must be less than the north bound.'
            )
        return self
    ###END def EcmwfArea.validate_north_south

    @pydantic.model_validator(mode='after')
    def validate_east_west(
        self,
    ) -> tp.Self:
        """Validate that the east bound is greater than the west bound."""
        if self.east <= self.west:
            raise ValueError(
                'The east bound must be greater than the west bound.'
            )
        return self
    ###END def EcmwfArea.validate_east_west

    @pydantic.model_validator(mode='before')
    @classmethod
    def parse_from_tuple[_DataType](
        cls,
        data: tp.Sequence[float] | _DataType,
    ) -> dict[str, float] | _DataType:
        """Parse an EcmwfArea from a sequence of four floats representing the
        North, West, South, and East bounds.

        Parameters
        ----------
        data : Sequence[float] | _DataType
            The data to parse. If a sequence of four floats, it is interpreted
            as (north, west, south, east). Otherwise, it is returned unchanged.

        Returns
        -------
        dict[Literal['north', 'west', 'south', 'east'], float] | _DataType
            If the input data was a sequence of four floats, a dictionary with
            the corresponding keys and values is returned. Otherwise, the input
            data is returned unchanged.
        """
        if (
                isinstance(data, tuple)
                and hasattr(data, '_fields')
                and hasattr(data, '_asdict')
        ):
            return tp.cast(tp.NamedTuple, data)._asdict()
        if (
            isinstance(data, tp.Sequence)
            and len(data) == 4
            and all(isinstance(_v, float) for _v in data)
        ):
            return {
                'north': data[0],
                'west': data[1],
                'south': data[2],
                'east': data[3],
            }
        return tp.cast(_DataType, data)
    ###END def EcmwfArea.parse_from_tuple

    def to_tuple(
        self,
    ) -> 'EcmwfAreaTuple':
        """Convert the EcmwfArea instance to an EcmwfAreaTuple.

        Returns
        -------
        EcmwfAreaTuple
            The corresponding EcmwfAreaTuple.
        """
        return EcmwfAreaTuple(
            north=self.north,
            west=self.west,
            south=self.south,
            east=self.east,
        )
    ###END def EcmwfArea.to_tuple

####END class EcmwfArea

class EcmwfAreaTuple(tp.NamedTuple):
    """A namedtuple representing a geographical area for ECMWF data requests.

    Attributes
    ----------
    north : float
        The northern bound in degrees (-90.0 to 90.0).
    west : float
        The western bound in degrees (-360.0 to 360.0).
    south : float
        The southern bound in degrees (-90.0 to 90.0).
    east : float
        The eastern bound in degrees (-360.0 to 360.0).
    """
    north: LatitudeFloat
    west: LongitudeFloat
    south: LatitudeFloat
    east: LongitudeFloat
###END class EcmwfAreaTuple
