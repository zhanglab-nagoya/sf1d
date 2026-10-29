"""1D Scheutjens-Fleer self-consistent field theory for polymer adsorption
on one wall or between two walls."""

from .core import SF1D, critical_chis_transfer_matrix, gap_scan, parse_sequence

__all__ = ["SF1D", "gap_scan", "parse_sequence", "critical_chis_transfer_matrix"]
__version__ = "0.2.0"
