"""Definitions and names of variables used in ERA5-Land to DATM conversion.

Enums
-----
Era5LandVar
    Enumeration of ERA5-Land variables available for download.
Datm7Var
    Enumeration of DATM7 variables.
Datm7Coord
    Enumeration of DATM7 coordinate variables.
Datm7Attr
    Enumeration of DATM7 variable attributes.

Classes
-------
Era5LandVarMapping
    A mapping from Era5LandVar to values of generic type.

Attributes
----------
era5land_request_varnames : dict[Era5LandVar, str]
    Mapping from Era5LandVar to the variable names used in the `variable` field
    of ERA5-Land data requests.
era5land_grib_varnames : dict[Era5LandVar, str]
    Mapping from Era5LandVar to the variable names used in the GRIB files that
    are downloaded from EAR5-Land data requests.
era5land_request_varnames_reverse : dict[str, Era5LandVar]
    Reverse mapping from variable names used in the `variable` field of
    ERA5-Land data requests to Era5LandVar.
era5land_grib_varnames_reverse : dict[str, Era5LandVar]
    Reverse mapping from variable names used in the GRIB files that are
    downloaded from EAR5-Land data requests to Era5LandVar.
era5_datm_vars : frozenset[Era5LandVar]
    The set of ERA5-Land variables needed for DATM.
datm7_var_ids : dict[Datm7Var|Datm7Coord, str]
    Variable IDs used in DATM7 netCDF files for each variable (including
    coordinate variables).
datm7_var_longname : dict[Datm7Var|Datm7Coord, str]
    Long name attribute used in DATM7 netCDF files for each variable.
datm7_var_units : dict[Datm7Var|Datm7Coord, str]
    Units attribute used in DATM7 netCDF files for each variable. Note that the
    time coordinate (Datm7Coord.TIME) does not have a unit atribute at all, and
    is not included in this mapping.
datm7_var_mode : dict[Datm7Var, str]
    Mode attribute used in DATM7 netCDF files for each data variable. The
    possible values are `'time-dependent'` and `'time-invariant'`.
datm7_var_attrs : dict[Datm7Var|Datm7Coord, dict[str, str]]
    Attributes of each variable used in DATM7 netCDF files.
"""
from collections import UserDict
import enum
import typing as tp

import pydantic
import pydantic_core.core_schema



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

    @property
    def request_varname(self) -> str:
        """The variable name used in the `variable` field of ERA5-Land data
        requests.
        """
        return era5land_request_varnames[self]
    ###END def Era5LandVar.request_varname

    @property
    def grib_varname(self) -> str:
        """The variable name used in the GRIB files that are downloaded from
        EAR5-Land data requests.
        """
        return era5land_grib_varnames[self]
    ###END def Era5LandVar.grib_varname

    @classmethod
    def _missing_(cls, value: object) -> 'Era5LandVar':
        """Check whether the value matches a grib varname or request varname,
        and return the corresponding Era5LandVar if so.
        """
        if isinstance(value, str):
            if (
                _value := era5land_request_varnames_reverse.get(value)
            ) is not None:
                return _value
            if (
                _value := era5land_grib_varnames_reverse.get(value)
            ) is not None:
                return _value
        return super()._missing_(value)
    ###END def Era5LandVar._missing_

###END class Era5LandVar

class VarSet(frozenset[Era5LandVar]):
    """A frozenset of Era5LandVar instances representing a set of variables.

    The class constructor accepts any iterable of Era5LandVar enum members or
    strings that can be converted into Era5LandVar enum members, but will raise
    an error if any element cannot be converted.
    """

    def __new__(
            cls,
            iterable: tp.Iterable[Era5LandVar|str],
            *args,
            **kwargs
    ) -> 'VarSet':
        return super().__new__(
            cls,
            (Era5LandVar(_value) for _value in iterable)
        )
    ###END def VarSet.__new__

    # We override __new__, so pydantic needs a custom schema
    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        source_type: tp.Any,
        handler: pydantic.GetCoreSchemaHandler,
    ) -> pydantic_core.core_schema.CoreSchema:
        input_schema = pydantic_core.core_schema.frozenset_schema(
            items_schema=pydantic_core.core_schema.enum_schema(
                cls=Era5LandVar,
                members=list(Era5LandVar),
                sub_type='str',
            )
        )
        return pydantic_core.core_schema.no_info_after_validator_function(
            function=lambda v: VarSet(v),
            schema=input_schema,
        )
    ###END def VarSet.__get_pydantic_core_schema__

###END class VarSet


class Era5LandVarMapping[_T](UserDict[Era5LandVar, _T]):
    """A mapping from Era5LandVar to values of generic type.

    This class is a simple subclass of dict with keys of type Era5LandVar and
    values of generic type _T.

    Init Parameters
    ---------------
    frozen : bool, optional
        If True, makes the mapping immutable after creation (note that as with
        all Python objects, this is a shallow immutability implemented only
        by overriding __setitem__). If the values are hashable, this also allows
        the mapping itself to be hashable. Default is False.
    """
    __slots__ = ('__frozen__',)

    def __init__(
        self,
        initial_data: tp.Mapping[Era5LandVar, _T] = {},
        *,
        frozen: bool = False,
        **kwargs: dict[Era5LandVar, _T],
    ) -> None:
        super().__init__(initial_data, **kwargs)
        self.__frozen__: tp.Final[bool] = frozen
    ###END def Era5LandVarMapping.__init__

    def __setitem__(self, key: Era5LandVar, value: _T) -> None:
        # We need to check whether self.__frozen__ exists, since it may not have
        # been set yet during the super().__init__ call in self.__init__
        if getattr(self, '__frozen__', False):
            raise TypeError(
                'This Era5LandVarMapping is frozen and cannot be modified.'
            )
        super().__setitem__(key, value)
    ###END def Era5LandVarMapping.__setitem__

    def __hash__(self) -> int:
        if not self.__frozen__:
            raise TypeError(
                'Only frozen Era5LandVarMapping instances are hashable.'
            )
        return hash(frozenset(self.items()))
    ###END def Era5LandVarMapping.__hash__

###END class Era5LandVarMapping


era5land_request_varnames: Era5LandVarMapping[str] = Era5LandVarMapping(
    {
        Era5LandVar.D2M: '2m_dewpoint_temperature',
        Era5LandVar.SP: 'surface_pressure',
        Era5LandVar.T2M: '2m_temperature',
        Era5LandVar.TP: 'total_precipitation',
        Era5LandVar.U10: '10m_u_component_of_wind',
        Era5LandVar.V10: '10m_v_component_of_wind',
        Era5LandVar.SSRD: 'surface_solar_radiation_downwards',
        Era5LandVar.STRD: 'surface_thermal_radiation_downwards',
    },
    frozen=True,
)

era5land_request_varnames_reverse: dict[str, Era5LandVar] = {
    _value: _key for _key, _value in era5land_request_varnames.items()
}

era5land_grib_varnames: Era5LandVarMapping[str] = Era5LandVarMapping(
    {
        Era5LandVar.D2M: 'd2m',
        Era5LandVar.SP: 'sp',
        Era5LandVar.T2M: 't2m',
        Era5LandVar.TP: 'tp',
        Era5LandVar.U10: 'u10',
        Era5LandVar.V10: 'v10',
        Era5LandVar.SSRD: 'ssrd',
        Era5LandVar.STRD: 'strd',
    },
    frozen=True,
)

era5land_grib_varnames_reverse: dict[str, Era5LandVar] = {
    _value: _key for _key, _value in era5land_grib_varnames.items()
}

era5_datm_vars: tp.Final[VarSet] = VarSet(
    {
        Era5LandVar.D2M,
        Era5LandVar.SP,
        Era5LandVar.T2M,
        Era5LandVar.TP,
        Era5LandVar.U10,
        Era5LandVar.V10,
        Era5LandVar.SSRD,
        Era5LandVar.STRD,
    }
)


class Datm7Var(enum.StrEnum):
    """Enumeration of DATM7 data variables.

    The string value of each member is in principle equal to the variable id
    that is used in DATM7 netCDF files, but you should use the ones given in
    the attribute `datm_var_ids` to be sure. The enums used here will not
    change and cannot be changed if different variable ids need to be used for
    some reason, but the mapping in `datm_varname` can be updated if needed.

    This enum class only contains data variables. Coordinate variables are
    instead given in the `Datm7Coord` enum.
    """

    TBOT = 'TBOT'
    PSRF = 'PSRF'
    QBOT = 'QBOT'
    WIND = 'WIND'
    FLDS = 'FLDS'
    FSDS = 'FSDS'
    PRECTmms = 'PRECTmms'

###END class Datm7Var

class Datm7Coord(enum.StrEnum):
    """Enumeration of DATM7 coordinate variables."""

    LONGXY = 'LONGXY'
    LATIXY = 'LATIXY'
    TIME = 'time'
    LON = 'lon'
    LAT = 'lat'

###END class Datm7Coord

datm7_var_ids: dict[Datm7Var|Datm7Coord, str] = {
    _member: _member.value for _member in Datm7Var
}

datm7_var_longname: dict[Datm7Var|Datm7Coord, str] = {
    Datm7Coord.TIME: 'observation time',
    Datm7Coord.LON: 'longitude',
    Datm7Coord.LAT: 'latitude',
    Datm7Coord.LONGXY: 'longitude',
    Datm7Coord.LATIXY: 'latitude',
    Datm7Var.PRECTmms: 'PRECTmms total precipitation',
    Datm7Var.FSDS: 'incident shortwave radiation',
    Datm7Var.TBOT: 'temperature at the lowest atm level',
    Datm7Var.PSRF: 'surface pressure at the lowest atm level',
    Datm7Var.QBOT: 'specific humidity at the lowest atm level',
    Datm7Var.WIND: 'wind at the lowest atm level',
    Datm7Var.FLDS: 'incident longwave radiation'
}

datm7_var_units: dict[Datm7Var|Datm7Coord, str] = {
    Datm7Coord.LON: 'degrees_east',
    Datm7Coord.LAT: 'degrees_north',
    Datm7Coord.LONGXY: 'degrees_east',
    Datm7Coord.LATIXY: 'degrees_north',
    Datm7Var.PRECTmms: 'mm H2O / sec',
    Datm7Var.FSDS: 'W/m**2',
    Datm7Var.TBOT: 'K',
    Datm7Var.PSRF: 'Pa',
    Datm7Var.QBOT: 'kg/kg',
    Datm7Var.WIND: 'm/s',
    Datm7Var.FLDS: 'W/m**2'
}

datm7_var_mode: dict[Datm7Var|Datm7Coord, str] = {
    Datm7Coord.LON: 'time-invariant',
    Datm7Coord.LAT: 'time-invariant',
    Datm7Coord.LONGXY: 'time-invariant',
    Datm7Coord.LATIXY: 'time-invariant',
    Datm7Var.PRECTmms: 'time-dependent',
    Datm7Var.FSDS: 'time-dependent',
    Datm7Var.TBOT: 'time-dependent',
    Datm7Var.PSRF: 'time-dependent',
    Datm7Var.QBOT: 'time-dependent',
    Datm7Var.WIND: 'time-dependent',
    Datm7Var.FLDS: 'time-dependent',
}

class Datm7Attr(enum.StrEnum):
    """Enumeration of DATM7 variable attributes."""

    LONG_NAME = 'long_name'
    UNITS = 'units'
    MODE = 'mode'

###END class Datm7Attr

datm7_var_attrs: dict[Datm7Var|Datm7Coord, dict[Datm7Attr, str]] = {
    _var: {
        _attrname: _dict[_var]
        for _attrname, _dict in (
            (Datm7Attr.LONG_NAME, datm7_var_longname),
            (Datm7Attr.UNITS, datm7_var_units),
            (Datm7Attr.MODE, datm7_var_mode),
        ) if _var in _dict
    }
    for _var in tuple(Datm7Var) + tuple(Datm7Coord)
}
