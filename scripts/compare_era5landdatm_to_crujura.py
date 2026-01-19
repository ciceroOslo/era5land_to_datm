"""Cell script to compare converted ERA5 Land data files to CRU JRA data.

This script is a cell script intended to be run in Visual Studio Code. If
needed, it can be exported from VS Code a Jupyter notebook and run in Jupyter or
compatible tools.
"""
 # %% [markdown]
# # Cell script to compare converted ERA5 Land data files to CRU JRA data.
#
# This script is a cell script intended to be run in Visual Studio Code. If
# needed, it can be exported from VS Code a Jupyter notebook and run in Jupyter or
# compatible tools.

# %% [markdown]
# ## Imports
# %%
from pathlib import Path
import typing as tp

import dask
from dask.diagnostics.progress import ProgressBar
import xarray as xr

from era5land_to_datm.dimensions import (
    Datm7Dim,
)
from era5land_to_datm.datm_streams import (
    Datm7Stream,
    datm7_stream_variables,
)
from era5land_to_datm.variables import (
    Datm7Var,
)

# %% [markdown]
# ## Initialize `dask` progress bar
# %%
pbar: ProgressBar = ProgressBar()
pbar.register()

# %% [markdown]
# ## Data files and paths
# 
# Set the paths to the directories that hold data files to be compared in the
# variables `era5datm_dir` and `crujra_dir` below, for converted ERA5 Land and
# CRU JRA data, respectively. Set a glob pattern or an expklicit list of files
# to include and compare from each data set in the variables `era5datm_files` and
# `crujra_files`, respectively.
#
# The CRU JRA data files will be downselected to the spatial and temporal range
# of the ERA5 Land data files for comparison. Since the ERA5 Land files are
# assumed to have higher resolution (both spatially and temporally), they will
# be interpolated to the coordinates of the CRU JRA data for comparison. Simple
# linear interpolation will be used given that the purpose is a direct
# comparison, not to produce a data set for model runs or that otherwise need to
# do conservative mapping or be as precise as possible.
# %%
era5datm_dir: Path = Path.home() / 'src/repos/norsink/temp_data/era5land/converted'
crujura_dir: Path = Path.home() / 'src/repos/norsink/temp_data/datm7'

era5datm_streamglobs: dict[Datm7Stream, str] = {
    _stream: f'clmforc.ERA5Land_NorwayRect0.1x0.1.{_stream.value}.2019_*.nc'
    for _stream in Datm7Stream
}

era5datm_files: dict[Datm7Stream, list[Path]] = {
    _stream: sorted(era5datm_dir.glob(era5datm_streamglobs[_stream]))
    for _stream in Datm7Stream
}
crujra_files: list[Path] = sorted(
    list(
        crujura_dir.glob('clmforc.CRUJRAv2.5_0.5x0.5.*.2019.nc')
    )
)

# %% [markdown]
# ## Open and merge data files
# 
# Open and merge the data files from each data set into single xarray Datasets.
# %% [markdown]
# ### Open ERA5 Land data files for each stream
# %%
drop_vars: list[str] = ['LATIXY', 'LONGXY']  # We don't need these, and they would trigger costly compatibility checks.
era5datm_mfdatasets: dict[Datm7Stream, xr.Dataset] = {}
for _stream, _files in era5datm_files.items():
    era5datm_mfdatasets[_stream] = xr.open_mfdataset(
        _files,
        drop_variables=drop_vars,
        chunks='auto',
        concat_dim=Datm7Dim.TIME.value,
        parallel=True,
        combine_attrs='override',
        data_vars='minimal',
        coords='minimal',
    )

# %% [markdown]
# ### Merge ERA5 Land streams into single Dataset
# %%
era5datm_ds: xr.Dataset = xr.merge(
    tuple(era5datm_mfdatasets.values()),
    join='outer',
    compat='identical',
    combine_attrs='override',
)

# %% [markdown]
# ### Obtain coordinate ranges of ERA5 Land data
# %%
era5datm_coord_ranges: dict[Datm7Dim, slice] = {
    Datm7Dim.LAT: slice(
        float(era5datm_ds[Datm7Dim.LAT.value].min()),
        float(era5datm_ds[Datm7Dim.LAT.value].max()),
    ),
    Datm7Dim.LON: slice(
        float(era5datm_ds[Datm7Dim.LON.value].min()),
        float(era5datm_ds[Datm7Dim.LON.value].max()),
    ),
    Datm7Dim.TIME: slice(
        float(era5datm_ds[Datm7Dim.TIME.value].min()),
        float(era5datm_ds[Datm7Dim.TIME.value].max()),
    ),
}

# ## Downselect and merge CRU JRA data files
# %%
crujra_ds: xr.Dataset = xr.merge(
    (
        xr.open_dataset(
            _file,
            chunks='auto',
        ).sel(
            era5datm_coord_ranges,
        ) for _file in crujra_files
    ),
    join='outer',
    combine_attrs='no_conflicts',
)
