"""Script to compare ERA5 forcing files included in DATM to ERA5 land grib files
downloaded from ECMWF.
"""

# %%
# Imports
# %%
from pathlib import Path

import xarray as xr


# %%
# Source directories
# %%
root_input_dir: Path = Path(
    '/cluster/shared/noresm/inputdata/'
)

datm_era5_dir: Path = (
    root_input_dir
    / 'atm/datm7/ERA5/'
)

downloaded_era5_dir: Path = (
    root_input_dir
    / 'cicero_mods/ERA5-Land_NorSink/original_GRIB_NorwayRect/'
)

# %%
# Select the input files. We try to choose a common year and month to be able to
# compare directly.
# %%
data_year: int = 2019  # Only 2019 currently available in DATM7 data
data_month: int = 5  # May

datm_era5_file: Path = (
    datm_era5_dir
    / f'ERA5.TL639.{data_year:04d}.{data_month:02d}.200618.nc'
)
downloaded_era5_file: Path = (
    downloaded_era5_dir
    / f'era5land_d2m_sp_ssrd_strd_t2m_tp_u10_v10_{data_year:04d}_{data_month:02d}.grib'
)

# %%
# Open the datasets
# %%
datm_ds: xr.Dataset = xr.open_dataset(datm_era5_file)
download_ds: xr.Dataset = xr.open_dataset(downloaded_era5_file, engine='cfgrib')
