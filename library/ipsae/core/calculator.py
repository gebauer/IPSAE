from pathlib import Path
from typing import Union, List, Dict, Tuple, Optional
import numpy as np
from scipy.spatial.distance import cdist

from .parser import PAEParser, StructureParser
from ..models.data_models import ChainChainScore, IPSAEResults, ResidueScore, StructureData

class IPSAECalculator:
    """Main class for calculating IPSAE scores."""
    
    def __init__(
        self,
        pae_file: Union[str, Path],
        structure_file: Union[str, Path],
        pae_cutoff: float = 30.0,
        distance_cutoff: float = 8.0
    ):
        """Initialize the calculator.
        
        Args:
            pae_file: Path to the PAE matrix file (JSON or NPZ)
            structure_file: Path to the structure file (PDB or CIF)
            pae_cutoff: PAE value cutoff for valid interactions
            distance_cutoff: Distance cutoff in Å for valid interactions
        """
        self.pae_file = Path(pae_file)
        self.structure_file = Path(structure_file)
        self.pae_cutoff = pae_cutoff
        self.distance_cutoff = distance_cutoff
        
        # Initialize data structures
        self.structure_data: StructureData = None
        self.pae_matrix: np.ndarray = None
        
        # Load data
        self.load_data()
    
    def load_data(self) -> None:
        """Load structure and PAE data from files."""
        # Load structure data
        if self.structure_file.suffix.lower() == '.cif':
            self.structure_data = StructureParser.parse_cif(self.structure_file)
        elif self.structure_file.suffix.lower() == '.pdb':
            self.structure_data = StructureParser.parse_pdb(self.structure_file)
        else:
            raise ValueError(f"Unsupported structure file format: {self.structure_file.suffix}")
        
        # Load PAE matrix
        # Handle gzipped files by checking the original extension
        pae_suffix = self.pae_file.suffix.lower()
        if pae_suffix == '.gz':
            # For gzipped files, check the original extension
            pae_suffix = '.' + self.pae_file.stem.split('.')[-1].lower()
        
        if pae_suffix == '.json':
            self.pae_matrix = PAEParser.parse_json(self.pae_file)
        elif pae_suffix == '.npz':
            self.pae_matrix = PAEParser.parse_npz(self.pae_file)
        else:
            raise ValueError(f"Unsupported PAE file format: {pae_suffix}")
        
        # Store PAE matrix in structure data
        self.structure_data.pae_matrix = self.pae_matrix
    
    def calculate_chain_chain_scores(self) -> List[ChainChainScore]:
        """Calculate chain-chain interaction scores."""
        scores = []
        
        # Get unique chains
        chains = sorted(set(self.structure_data.chains))
        
        # Calculate distances between all CA atoms
        coordinates = np.array(self.structure_data.coordinates)
        distances = cdist(coordinates, coordinates)
        
        for i, chain1 in enumerate(chains):
            for j, chain2 in enumerate(chains[i+1:], i+1):
                # Get residue indices for each chain
                chain1_indices = [k for k, (c, _, _) in enumerate(self.structure_data.residues) if c == chain1]
                chain2_indices = [k for k, (c, _, _) in enumerate(self.structure_data.residues) if c == chain2]
                
                # Get valid interactions
                valid_interactions = []
                pae_values = []
                dist_values = []
                
                for idx1 in chain1_indices:
                    for idx2 in chain2_indices:
                        pae = self.pae_matrix[idx1, idx2]
                        dist = distances[idx1, idx2]
                        
                        if pae <= self.pae_cutoff and dist <= self.distance_cutoff:
                            valid_interactions.append((idx1, idx2))
                            pae_values.append(pae)
                            dist_values.append(dist)
                
                if valid_interactions:
                    # Calculate scores
                    avg_pae = np.mean(pae_values)
                    avg_distance = np.mean(dist_values)
                    
                    # Calculate ipSAE score (example formula, adjust as needed)
                    ipsae_score = 1.0 - (avg_pae / self.pae_cutoff)
                    
                    # Calculate pDockQ score (example formula, adjust as needed)
                    pdockq_score = 1.0 / (1.0 + np.exp(-0.5 * (avg_pae - 20)))
                    
                    # Calculate LIS score (example formula, adjust as needed)
                    lis_score = len(valid_interactions) / (len(chain1_indices) * len(chain2_indices))
                    
                    scores.append(ChainChainScore(
                        chain1=chain1,
                        chain2=chain2,
                        ipsae_score=ipsae_score,
                        pdockq_score=pdockq_score,
                        lis_score=lis_score,
                        avg_pae=avg_pae,
                        avg_distance=avg_distance,
                        num_interactions=len(valid_interactions),
                        valid_interactions=valid_interactions,
                        pae_values=pae_values,
                        distances=dist_values
                    ))
        
        return scores
    
    def calculate_residue_scores(self) -> List[ResidueScore]:
        """Calculate residue-level scores."""
        scores = []
        
        # Calculate distances between all CA atoms
        coordinates = np.array(self.structure_data.coordinates)
        distances = cdist(coordinates, coordinates)
        
        for i, (chain, res_num, res_name) in enumerate(self.structure_data.residues):
            # Get valid interactions for this residue
            valid_interactions = []
            pae_values = []
            dist_values = []
            
            for j in range(len(coordinates)):
                if i != j:  # Skip self-interactions
                    pae = self.pae_matrix[i, j]
                    dist = distances[i, j]
                    
                    if pae <= self.pae_cutoff and dist <= self.distance_cutoff:
                        valid_interactions.append((i, j))
                        pae_values.append(pae)
                        dist_values.append(dist)
            
            if valid_interactions:
                # Calculate scores
                avg_pae = np.mean(pae_values)
                avg_distance = np.mean(dist_values)
                
                # Calculate ipSAE score (example formula, adjust as needed)
                ipsae_score = 1.0 - (avg_pae / self.pae_cutoff)
                
                scores.append(ResidueScore(
                    chain=chain,
                    residue_number=res_num,
                    residue_name=res_name,
                    ipsae_score=ipsae_score,
                    avg_pae=avg_pae,
                    avg_distance=avg_distance,
                    num_interactions=len(valid_interactions),
                    valid_interactions=valid_interactions,
                    pae_values=pae_values,
                    distances=dist_values
                ))
        
        return scores
    
    def calculate(self) -> IPSAEResults:
        """Perform all calculations and return results."""
        self.load_data()
        
        chain_chain_scores = self.calculate_chain_chain_scores()
        residue_scores = self.calculate_residue_scores()
        
        return IPSAEResults(
            chain_chain_scores=chain_chain_scores,
            residue_scores=residue_scores,
            structure_data=self.structure_data
        )
