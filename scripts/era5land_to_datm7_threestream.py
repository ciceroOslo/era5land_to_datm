"""Script to convert ERA5 Land GRIB files to netCDF files suitable for DATM7.

This script converts the ERA5 Land data files to the threestream format used by
the CRUJRA mode of DATM7. It does not use the more recent ERA5 mode, which did
not seem to be stable or consistently documented at the time of writing
(January 2026).
"""
import enum
from pathlib import Path
import typing as tp

from korsbakken_python_utils.containers.dataobject import UniformTypeDataObject
import numpy as np
import xarray as xr



class Era5Dim(enum.StrEnum):
    """String enum for ERA5 Land dimension ids."""
    
    DATE = 'time'
    STEP = 'step'

###END class Era5Dim

class Era5Var(enum.StrEnum):
    """String enum for ERA5 Land variable identifiers."""


def main(
        *,
        source_file: Path,
        next_source_file: Path|None = None,
        output_file: Path,
        eager: bool = True,
) -> None:
    """Converts an ERA5 Land GRIB file to a DAMT7 threestream netCDF file.
    
    Parameters
    ----------
    source_file : Path
        Path to the ERA5 Land GRIB file for the current month or other time
        period.
    next_source_file : Path | None, optional
        Path to the ERA5 Land GRIB file for the subsequent time period. Some
        variables require data from the first time step of the next period,
        notably radiation variables, which are cumulative in the ERA5 land data,
        but intensity-based in DATM7. If None, these variables will have missing
        values in the last time step of the output. Optional, by default None.
    output_file : Path
        Path to the output netCDF file.
    eager : bool, optional
        Whether to load the entire dataset into memory before processing. This
        can significantly speed up some operations, while setting it to False
        may lead require less memory usage but lead to some data being reloaded
        and even reprocesseed multiple times, which can slow down the
        processing. By default True. Set to False if you run into memory issues.
    """
    source_ds: xr.Dataset = open_era5land_grib(
        file=source_file,
        next_file=next_source_file,
        use_chunks=not eager,
    )



def open_era5land_grib(
        file: Path,
        *,
        next_file: Path|None = None,
        use_chunks: bool = False,
        date_dim: str = Era5Dim.DATE,
        step_dim: str = Era5Dim.STEP,
) -> xr.Dataset:
    """Opens an ERA5 Land GRIB file, optionally with the next file for
    cumulative variables.
    
    Parameters
    ----------
    file : Path
        Path to the ERA5 Land GRIB file.
    next_file : Path | None, optional
        Path to the next ERA5 Land GRIB file for cumulative variables. If given,
        it will be lazily opened, and the first time step will be extracted and
        concatenated to the main dataset along the time dimension.  Optional, by
        default None.
    use_chunks : bool, optional
        Whether to use chunking (with dask). If False, the returned dataset will
        use xarray's default lazy loading mechanism without chunking. By default
        False.
    time_dim : str, optional
        Name of the time dimension in the dataset. By default given by the
        string enum `Era5Var.TIME`.
    """
    chunk_option: str|int|None = 'auto' if use_chunks else None
    era5_ds: xr.Dataset = xr.open_dataset(
        file,
        chunks=chunk_option,
    )
    if next_file is not None:
        next_ds: xr.Dataset = xr.open_dataset(
            next_file,
            chunks='auto',
        ).isel({date_dim: 0, step_dim: -1})
        source_last_date = era5_ds[date_dim].isel({date_dim: -1}).item()
        next_first_date = next_ds[date_dim].item()
        if source_last_date != next_first_date:
            raise ValueError(
                f'The last date in the source file ({source_last_date}) does not '
                f'match the first date in the next file ({next_first_date}).'
            )
        if not set(era5_ds.variables).issuperset(set(next_ds.variables)):
            raise ValueError(
                'The next file contains variables that are not present in the '
                'source file.'
            )
