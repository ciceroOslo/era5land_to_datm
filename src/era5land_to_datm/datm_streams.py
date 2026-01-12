"""Definitions and IDs related to DATM7 data streams.

Enums
-----
Datm7Stream
    Enumeration of DATM7 data streams.

Attributes
----------
datm7_stream_variables : dict[Datm7Stream, frozenset[Datm7Var]]
    Mapping from Datm7Stream to the set of variable names included in that
    stream.
"""
import enum
import typing as tp

from era5land_to_datm.variables import Datm7Var



class Datm7Stream(enum.StrEnum):
    """Enumeration of DATM7 data streams.
    This enum contains IDs for the three data streams used in the DATM7

    threestream mode, used e.g. by the CRUNCEP and CRUJRA modes.

    The values correspond to the stream ID used in data file names, as of
    January 2026.
    """

    SOLR = 'Solr'
    PREC = 'Prec'
    TPQWL = 'TPQWL'

###END class Datm7Stream

datm7_stream_variables: tp.Final[dict[Datm7Stream, frozenset[Datm7Var]]] = {
    Datm7Stream.SOLR: frozenset(
        {
            Datm7Var.FSDS,
        }
    ),
    Datm7Stream.PREC: frozenset(
        {
            Datm7Var.PRECTmms
        }
    ),
    Datm7Stream.TPQWL: frozenset(
        {
            Datm7Var.TBOT,
            Datm7Var.PSRF,
            Datm7Var.QBOT,
            Datm7Var.WIND,
            Datm7Var.FLDS,
        }
    ),
}
