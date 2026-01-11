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
        Path to the ERA5 Land GRIB file for the current month or other timer
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
        eager=eager,
    )



def open_era5land_grib(
        file: Path,
        *,
        next_file: Path|None = None,
        eager: bool = True,
        time_dim: str = Era5Var.TIME,
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
    eager : bool, optional
        Whether to load the entire dataset into memory before returning. Note
        that this option only affects the main file (`file`). `next_file` is
        always opened lazily before extracting the first time step.
    time_dim : str, optional
        Name of the time dimension in the dataset. By default given by the
        string enum `Era5Var.TIME`.
    """
