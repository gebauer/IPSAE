"""
IPSAE (Interface Protein Structure Assessment Engine)
A library for analyzing protein-protein interactions using AlphaFold PAE scores.
"""

from .core.calculator import IPSAECalculator
from .core.parser import PAEParser, StructureParser
from .models.data_models import ChainChainScore, ResidueScore

__version__ = "0.1.8"
__all__ = [
    "IPSAECalculator",
    "PAEParser",
    "StructureParser",
    "ChainChainScore",
    "ResidueScore"
]
