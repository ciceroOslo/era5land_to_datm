"""Definitions of dimension names used in ERA5-Land to DATM conversion.

Enums
-----
Era5Dim
    String enum for ERA5 Land dimension ids.
"""
import enum



class Era5LandDim(enum.StrEnum):
    """String enum for ERA5 Land dimension ids."""

    DATE = 'time'
    STEP = 'step'

###END class Era5Dim
