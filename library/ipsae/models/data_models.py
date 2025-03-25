from dataclasses import dataclass, asdict
from typing import List, Optional, Tuple
import json
from pathlib import Path
import numpy as np

@dataclass
class ChainChainScore:
    """Score for a chain-chain interaction."""
    chain1: str
    chain2: str
    ipsae_score: float  # Changed from score to ipsae_score
    pdockq_score: float
    lis_score: float
    avg_pae: float
    avg_distance: float
    num_interactions: int
    valid_interactions: List[Tuple[int, int]]  # List of (residue1, residue2) pairs
    pae_values: List[float]  # PAE values for valid interactions
    distances: List[float]  # Distances for valid interactions

@dataclass
class ResidueScore:
    """Score for a single residue."""
    chain: str
    residue_number: int
    residue_name: str
    ipsae_score: float  # Changed from score to ipsae_score
    avg_pae: float
    avg_distance: float
    num_interactions: int
    valid_interactions: List[Tuple[int, int]]  # List of (residue1, residue2) pairs
    pae_values: List[float]  # PAE values for valid interactions
    distances: List[float]  # Distances for valid interactions

@dataclass
class StructureData:
    """Data for a protein structure."""
    chains: List[str]
    residues: List[Tuple[str, int, str]]  # List of (chain, residue_number, residue_name)
    coordinates: List[np.ndarray]  # List of CA coordinates
    pae_matrix: Optional[np.ndarray] = None  # Optional PAE matrix

@dataclass
class IPSAEResults:
    """Results from IPSAE calculation."""
    chain_chain_scores: List[ChainChainScore]
    residue_scores: List[ResidueScore]
    structure_data: StructureData
    
    def save(self, file_path: Path) -> None:
        """Save results to a JSON file."""
        def convert_numpy(obj):
            """Recursively convert numpy arrays to lists."""
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, (list, tuple)):
                return [convert_numpy(item) for item in obj]
            elif isinstance(obj, dict):
                return {key: convert_numpy(value) for key, value in obj.items()}
            return obj
        
        # Convert all numpy arrays to lists
        data = {
            'chain_chain_scores': [
                {
                    'chain1': score.chain1,
                    'chain2': score.chain2,
                    'ipsae_score': score.ipsae_score,
                    'pdockq_score': score.pdockq_score,
                    'lis_score': score.lis_score,
                    'avg_pae': score.avg_pae,
                    'avg_distance': score.avg_distance,
                    'num_interactions': score.num_interactions,
                    'valid_interactions': score.valid_interactions,
                    'pae_values': convert_numpy(score.pae_values),
                    'distances': convert_numpy(score.distances)
                }
                for score in self.chain_chain_scores
            ],
            'residue_scores': [
                {
                    'chain': score.chain,
                    'residue_number': score.residue_number,
                    'residue_name': score.residue_name,
                    'ipsae_score': score.ipsae_score,
                    'avg_pae': score.avg_pae,
                    'avg_distance': score.avg_distance,
                    'num_interactions': score.num_interactions,
                    'valid_interactions': score.valid_interactions,
                    'pae_values': convert_numpy(score.pae_values),
                    'distances': convert_numpy(score.distances)
                }
                for score in self.residue_scores
            ],
            'structure_data': {
                'chains': self.structure_data.chains,
                'residues': self.structure_data.residues,
                'coordinates': convert_numpy(self.structure_data.coordinates),
                'pae_matrix': convert_numpy(self.structure_data.pae_matrix) if self.structure_data.pae_matrix is not None else None
            }
        }
        
        # Save to JSON file
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
