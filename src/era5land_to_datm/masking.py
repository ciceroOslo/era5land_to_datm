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

Functions
---------
make_spatial_null_mask
    Make a latitude/longitude grid mask that is False for grid cells that have
    only null values for all time values and all data variables.
"""
import enum
import typing as tp

import xarray as xr

from .convert_data import (
    era5land_to_linear_time,
)
from .dimensions import (
    Era5LandLinearizedTimeDimId,
    ERA5LandTimeLayout,
)



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


def make_spatial_null_mask(
        data: xr.Dataset,
        *,
        time_dim: str = Era5LandLinearizedTimeDimId.TIME.value,
        time_layout: ERA5LandTimeLayout = ERA5LandTimeLayout.LINEAR,
        mask_var_name: str = ERA5_LAND_TO_DATM_MASK_VAR,
) -> xr.Dataset:
    """Make a mask that is False where the data have only null values.
    
    The function returns a mask which is False for grid cells where all
    variables have only null values for all time values, and True otherwise.

    Parameters
    ----------
    data : xr.Dataset
        The dataset for which to make the mask. Should have a time dimension and
        one or more data variables.
    time_dim : str, optional
        The name of the time dimension in the dataset, if the time dimension is
        not in the 2D date+step format used in original ERA5 Land data files. It
        is ignored if `time_layout` is `ERA5LandTimeLayout.DATE_STEP`. By
        default, given by `Era5LandLinearizedTimeDimId.TIME`.
    time_layout : ERA5LandTimeLayout, optional
        The time layout of the dataset. If `ERA5LandTimeLayout.DATE_STEP`, the
        function will look use the function `era5land_to_linear_time` with
        default settings to first convert the dataset to a linear time layout.
        This option assumes that dimensions, coordinate names and time step
        layout are identical to those used in original ERA5 Land data files. By
        default `ERA5LandTimeLayout.LINEAR` (i.e., do not convert the time
        layout).
    mask_var_name : str, optional
        The name to use for the mask variable in the returned Dataset. By default
        given by `ERA5_LAND_TO_DATM_MASK_VAR`.

    Returns
    -------
    xr.Dataset
        A dataset containing a single variable with name given by
        `mask_var_name` or its default, and boolean values, which is False for
        grid cells where all variables have only null values for all time
        values, and True otherwise. The dataset will have the same spatial
        dimensions and coordinates as the input dataset.
    """

    if time_layout == ERA5LandTimeLayout.DATE_STEP:
        data = era5land_to_linear_time(
            data,
            output_time_dim=time_dim,
        )

    # Create a unique variable dimension name to turn the Dataset into a single
    # DataArray
    variable_dim: str = 'variable'
    while variable_dim in (
            set(data.dims)
            | set(data.variables)
            | {mask_var_name}
    ):
        variable_dim += '_'
    mask = (
        data
        .to_dataarray(dim=variable_dim)
        .notnull()
        .any(dim=[time_dim, variable_dim])
    ).to_dataset(name=mask_var_name)

    return mask

###END def make_spatial_null_mask
