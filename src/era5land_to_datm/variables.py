"""Definitions and names of variables used in ERA5-Land to DATM conversion.

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
era5_datm_vars : frozenset[Era5LandVar]
    The set of ERA5-Land variables needed for DATM.
"""
from collections import UserDict
import typing as tp

from .types import (
    Era5LandVar,
    VarSet,
)



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
        if self.__frozen__:
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

era5_datm_vars: tp.Final[VarSet] = frozenset(
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
