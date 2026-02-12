"""Functionality related to masking of ERA5 Land data.

Enums
-----
MaskedValuesHandling
    Enum for what to do if data contains non-null values in a masked region
    (i.e., in a region that is masked out and not intended to be used).
UnmaskedNullsHandling
    Enum for what to do if data contains null values in an unmasked region
    (i.e., in a region that is not masked out and is intended to be used).
UnmaskedNullsProcessing
    Enum for how to process null values in an unmasked region, if any are found.
Era5LandToDatmMaskFileDim
    String enum for dimension ids the mask file used in the `convert_files` and
    `convert_data` modules.

Attributes
----------
ERA5_LAND_TO_DATM_MASK_VAR : Final[str]
    The variable name in the mask file used in the `convert_files` and
    `convert_data` modules that contains the mask.
"""
import enum
import typing as tp



class MaskedValuesHandling(enum.Enum):
    """Enum for what to do if data contains non-null values in a masked region
    (i.e., in a region that is masked out and not intended to be used).

    Members
    -------
    RAISE
        Raise a ValueError if any non-null values are found in the masked areas.
    WARN
        Log a warning if any non-null values are found in the masked areas, but
        do not raise an error.
    IGNORE
        Do not report any non-null values found in the masked areas.
    """
    RAISE = 'raise'
    WARN = 'warn'
    IGNORE = 'ignore'
###END class MaskedValuesHandling


class UnmaskedNullsHandling(enum.Enum):
    """Enum for what to do if data contains null values in an unmasked region
    (i.e., in a region that is not masked out and is intended to be used).

    Members
    -------
    RAISE
        Raise a ValueError if any null values are found in the unmasked areas.
    WARN
        Log a warning if any null values are found in the unmasked areas, but do
        not raise an error.
    IGNORE
        Do not report any null values found in the unmasked areas.
    """
    RAISE = 'raise'
    WARN = 'warn'
    IGNORE = 'ignore'
###END class UnmaskedNullsHandling


class UnmaskedNullsProcessing(enum.Enum):
    """Enum for how to process null values in an unmasked region.

    Members
    -------
    NONE
        Leave the unmasked null values as nulls. **NB!** Note that this option
        may still imply that masked values are converted to nulls, so may not
        return the original data set unchanged.
    LINEAR
        Fill the unmasked null values with 1D linear interpolation (usually
        along the time dimension).
    """
    NONE = 'none'
    LINEAR = 'linear'
###END class UnmaskedNullsProcessing


class Era5LandToDatmMaskFileDim(enum.StrEnum):
    """String enum for dimension ids the mask file used in the `convert_files`
    and `convert_data` modules.

    Members
    -------
    LAT
        Latitude dimension id.
    LON
        Longitude dimension id.
    """
    LAT = 'latitude'
    LON = 'longitude'
###END class Era5LandToDatmMaskFileDim

ERA5_LAND_TO_DATM_MASK_VAR: tp.Final[str] = 'mask'
