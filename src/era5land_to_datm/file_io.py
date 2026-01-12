"""File input and output of ERA5 and DATM data files.

Functions
---------
open_era5land_grib
    Opens an ERA5 Land GRIB file, optionally with the next file to obtain values
    in the last time step when these are NaN. Returns data as an xarray Dataset.
"""
from pathlib import Path

import xarray as xr



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
    date_dim : str, optional
        Name of the time/date dimension in the ERA5 Land dataset. By default
        given by the string enum `Era5Dim.DATE`.
    step_dim : str, optional
        Name of the "step" dimension (intra-date timediff relative to midnight)
        in the ERA5 Land dataset. By default given by the string enum
        `Era5Dim.STEP`.

    Returns
    -------
    xr.Dataset
        The ERA5 Land dataset.

    Raises
    ------
    ValueError
        If the last date in the source file does not match the first date in the
        next file, if the last intra-date time step in the source file does not
        match the *last* intra-date time step in the next file, or if the next
        file contains variables that are not present in the source file.
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
        source_last_step = era5_ds[step_dim].isel({step_dim: -1}).item()
        next_last_step = next_ds[step_dim].item()
        if source_last_step != next_last_step:
            raise ValueError(
                f'The last intra-date time step in the source file '
                f'({source_last_step}) does not match the last intra-date time '
                f'step in the next file ({next_last_step}).'
            )
        for _var in era5_ds.data_vars:
            era5_ds[_var].loc[
                {
                    date_dim: source_last_date,
                    step_dim: source_last_step,
                }
            ] = next_ds[_var]
    return era5_ds
###END def open_era5land_grib