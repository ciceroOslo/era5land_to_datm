"""Code to download ERA5 Land GRIB files from completed requests.

Classes
-------
DownloadFilesResult
    Represents the result of attempting to download files for multiple requests
    through the `retrieve_available_files` function. Provides information on
    which files were successfully downloaded, which were not available yet, and
    which downloads failed.

Functions
---------
get_remotes
    Retrieve Remote instances for previously sent requests.
get_remote_varset
    Get the VarSet associated with a given Remote instance (i.e., the set of
    variables requested in that Remote).
get_remote_yearmonth
    Get the YearMonth associated with a given Remote instance (i.e., the year
    and month requested in that Remote).
make_grib_filename
    Create a standardized filename for the downloaded GRIB file based on the
    set of variables and year/month requested.
remotes_dict_by_vars_and_yearmonth
    Organize Remote instances into a nested dictionary keyed by VarSet and
    YearMonth, based on the `variable`, `year`, and `month` attributes of each
    Remote.
retrieve_available_files
    Attempt to download available files for multiple Remote instances,
    returning a DownloadFilesResult with details on successes, unavailable
    files, and failures.
"""