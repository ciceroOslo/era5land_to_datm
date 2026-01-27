"""Functions to convert encoding or datatypes for existing data stream files.

The functions in this module are used to write netCDF files with altered
encoding settings or datatypes for variables, based on existing datasets or
netCDF files.

Attributes
----------
TIME_DIM: str
    Default name of the time dimension used in the datasets.

Functions
---------
merge_encoding_dicts
    Merges two or more dictionaries with per-variable encoding settings.
make_float_encoding_dict
    Makes a dict that sets all floating-point variables to a specified datatype.
make_calendar_encoding_dict
    Makes a dict that sets the time variable to a specified calendar type.
make_time_units_encoding_dict
    Makes a dict that sets the time variable to specified time units.
make_time_units_string
    Convenience function to make a time units string of the type `units since
    reference date/time` for a given unit and reference date/time.

Enums
-----
FillValue
    Enum with default fill values for different data types.
DType
    Enum with different data types used in encoding conversion.
CalType
    Enum with different calendar types used in encoding conversion.
EncodingField
    Enum with the names of different fields in a netCDF variable encoding
    dictionary.
"""
from collections.abc import (
    Hashable,
    Iterable,
    Mapping,
    Sequence
)
import datetime
import enum
import itertools
import logging
from pathlib import Path
import typing as tp

import netCDF4 as nc4
import numpy as np
import xarray as xr



logger: tp.Final[logging.Logger] = logging.getLogger(__name__)

TIME_DIM: tp.Final[str] = 'time'


class FillValue(enum.Enum):
    """Enum with default fill values for different data types."""

    FLOAT32 = nc4.default_fillvals['f4']
    FLOAT64 = nc4.default_fillvals['f8']
    FLOAT = nc4.default_fillvals['f4']

    INT8 = nc4.default_fillvals['i1']
    INT16 = nc4.default_fillvals['i2']
    INT32 = nc4.default_fillvals['i4']
    INT64 = nc4.default_fillvals['i8']
    UINT8 = nc4.default_fillvals['u1']
    UINT16 = nc4.default_fillvals['u2']
    UINT32 = nc4.default_fillvals['u4']
    UINT64 = nc4.default_fillvals['u8']

###END class FillValue


class DType(enum.StrEnum):
    """Enum with different data types used in encoding conversion."""

    FLOAT32 = 'float32'
    FLOAT64 = 'float64'
    INT8 = 'int8'
    INT16 = 'int16'
    INT32 = 'int32'
    INT64 = 'int64'
    UINT8 = 'uint8'
    UINT16 = 'uint16'
    UINT32 = 'uint32'
    UINT64 = 'uint64'

    @property
    def np_type(self) -> type[np.generic]:
        """The corresponding numpy type for this DType enum member."""
        return np.dtype(self.value).type
    ###END def DType.np_dtype

###END class DType


class EncodingField(enum.StrEnum):
    """Enum with the names of different fields in a netCDF variable encoding dictionary."""

    DATATYPE = 'dtype'
    FILLVALUE = '_FillValue'
    CALENDAR = 'calendar'
    UNITS = 'units'

###END class EncodingField


class CalType(enum.StrEnum):
    """Enum with different calendar types used in encoding conversion.

    The valid values are those supported by the Python netCDF4 library. To avoid
    confusion and inconsistency with the latest CF conventions, `gregorian` is
    not supported; use `proleptic_gregorian` if you want Gregorian leap year
    rules for all years, or `standard` if you want the historical calendar that
    uses Gregorian rules from 1582 onwards and Julian rules and dates before
    that (with invalid dates from 1582-10-05 through 1582-10-14).
    """

    STANDARD = 'standard'
    JULIAN = 'julian'
    PROLEPTIC_GREGORIAN = 'proleptic_gregorian'
    NOLEAP = 'noleap'
    ALL_LEAP = 'all_leap'
    THREESIXTY_DAY = '360_day'

###END class CalType


VarEncoding: tp.TypeAlias = dict[EncodingField, str|int|float]
EncodingDict: tp.TypeAlias = dict[Hashable, VarEncoding]


def merge_encoding_dicts(
        dicts: Iterable[EncodingDict],
) -> EncodingDict:
    """Merges two or more dictionaries with per-variable encoding settings.

    The function assumes that each dict is a two-level nested dictionary where
    the outer keys are variable names and the inner keys are encoding attributes
    for the corresponding variable. It does not assume that the same variables
    or the same encoding attributes are present in all input dictionaries. If
    an attribute for a variable is present in multiple input dictionaries, the
    settings from later dictionaries will override those from earlier ones.

    Parameters
    ----------
    dicts : Iterable[EncodingDict]
        Iterable of nested dictionaries to merge. Each dictionary should have
        variable names as outer keys and encoding setting dictionaries as
        values.

    Returns
    -------
    EncodingDict
        Merged dictionary with combined encoding settings for each variable.
    """
    variables: set[Hashable] = set(
        itertools.chain.from_iterable(
            _enc_dict.keys() for _enc_dict in dicts
        )
    )
    merged_dict: EncodingDict = {
        _var_name: {} for _var_name in variables
    }

    for _enc_dict in dicts:
        for _var_name, _var_enc in _enc_dict.items():
            if _var_name not in merged_dict:
                merged_dict[_var_name] = {}
            merged_dict[_var_name].update(_var_enc)
    return merged_dict
###END def merge_encoding_dicts


def make_float_encoding_dict(
        ds: xr.Dataset,
        dtype: DType | tp.Literal[ 'float32', 'float64', ],
        *,
        include_dtype: bool = True,
        include_fillvalue: bool = True,
) -> EncodingDict:
    """Makes a dict that sets all floating-point variables to a specified
    datatype.

    This function also sets the fill value by default. It can be used to set
    only the dtype or only the fill value by setting the corresponding keyword
    arguments to False.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset from which to get the variable names and datatypes.
    dtype : DType | tp.Literal['float32', 'float64']
        Target datatype for all floating-point variables.
    include_dtype : bool, optional
        Whether to include the datatype setting in the output dictionary.
        Default is True.
    include_fillvalue : bool, optional
        Whether to include the fill value setting in the output dictionary.
        Default is True.

    Returns
    -------
    EncodingDict
        Dictionary with per-variable encoding settings for all floating-point
        variables in the input dataset, setting their datatype to the specified
        target datatype.
    """
    if isinstance(dtype, str):
        dtype = DType[dtype.upper()]

    inner_dict: dict[EncodingField, str|int|float] = {}
    if include_dtype:
        inner_dict[EncodingField.DATATYPE] = dtype
    if include_fillvalue:
        inner_dict[EncodingField.FILLVALUE] = dtype

    encoding_dict: EncodingDict = {}

    for _var_name, _var in ds.data_vars.items():
        if np.issubdtype(_var.dtype, np.floating):
            encoding_dict[_var_name] = inner_dict.copy()

    return encoding_dict

###END def make_float_encoding_dict


def make_calendar_encoding_dict(
        cal_type: CalType
) -> EncodingDict:
    """Makes a dict that sets the time variable to a specified calendar type.

    Parameters
    ----------
    cal_type : CalType
        Target calendar type for the time variable.

    Returns
    -------
    EncodingDict
        Dictionary with per-variable encoding settings for the time variable,
        setting its calendar type to the specified target calendar type.
    """
    return {
        TIME_DIM: {
            EncodingField.CALENDAR: cal_type.value
        }
    }
###END def make_calendar_encoding_dict


def make_time_units_encoding_dict(
        reference_time: datetime.datetime | str,
        *,
        time_unit: tp.Literal['days', 'hours'] = 'days',
        dtype: DType | tp.Literal['float32', 'float64'] = DType.FLOAT32,
) -> EncodingDict:
    """Makes a dict that sets the time variable to specified time units.

    Parameters
    ----------
    reference_time : datetime.datetime | str
        Reference date/time for the time units. If a string is provided, it
        should be in a format that can be parsed by `datetime.datetime.fromisoformat`.
    time_unit : tp.Literal['days', 'hours'], optional
        Time unit to use for the time variable. Default is 'days'.
    dtype : DType | tp.Literal['float32', 'float64'], optional
        Datatype to use for the time variable. Default is DType.FLOAT32.

    Returns
    -------
    EncodingDict
        Dictionary with per-variable encoding settings for the time variable,
        setting its time units to the specified target time units. The units
        string will be in the format `time_unit since reference_time`, where
        `reference_time` is formatted as `YYYY-MM-DD HH:MM:SS`.
    """
    units_str = make_time_units_string(
        time_unit=time_unit,
        reference_time=reference_time,
    )

    return {
        TIME_DIM: {
            EncodingField.UNITS: units_str,
            EncodingField.DATATYPE: str(dtype),
        }
    }
###END def make_time_units_encoding_dict

def make_time_units_string(
        *,
        time_unit: tp.Literal['days', 'hours'],
        reference_time: datetime.datetime | str,
) -> str:
    """Convenience function to make a time units string of the type `units since
    reference date/time` for a given unit and reference date/time.

    Parameters
    ----------
    time_unit : tp.Literal['days', 'hours']
        Time unit to use for the time variable.
    reference_time : datetime.datetime | str
        Reference date/time for the time units. If a string is provided, it
        should be in a format that can be parsed by
        `datetime.datetime.fromisoformat`.

    Returns
    -------
    str
        Time units string in the format `time_unit since reference_time`, where
        `reference_time` is formatted as `YYYY-MM-DD HH:MM:SS`.
    """
    if isinstance(reference_time, str):
        reference_time = datetime.datetime.fromisoformat(reference_time)
    return f'{time_unit} since {reference_time:%Y-%m-%d %H:%M:%S}'
###END def make_time_units_string
