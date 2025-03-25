import json
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Union
import logging

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

from ..models.data_models import StructureData

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

class StructureParser:
    """Parser for structure files (CIF and PDB)."""
    
    @staticmethod
    def parse_cif(file_path: Union[str, Path]) -> StructureData:
        """Parse a CIF file and return structure data."""
        parser = MMCIFParser()
        structure = parser.get_structure('structure', str(file_path))
        
        chains = []
        residues = []
        coordinates = []
        
        for model in structure:
            for chain in model:
                chains.append(chain.id)
                for residue in chain:
                    residues.append((chain.id, residue.id[1], residue.resname))
                    for atom in residue:
                        if atom.id == 'CA':  # Only store CA coordinates
                            coordinates.append(atom.get_coord())
        
        return StructureData(
            chains=chains,
            residues=residues,
            coordinates=coordinates
        )
    
    @staticmethod
    def parse_pdb(file_path: Union[str, Path]) -> StructureData:
        """Parse a PDB file and return structure data."""
        parser = PDBParser()
        structure = parser.get_structure('structure', str(file_path))
        
        chains = []
        residues = []
        coordinates = []
        
        for model in structure:
            for chain in model:
                chains.append(chain.id)
                for residue in chain:
                    residues.append((chain.id, residue.id[1], residue.resname))
                    for atom in residue:
                        if atom.id == 'CA':  # Only store CA coordinates
                            coordinates.append(atom.get_coord())
        
        return StructureData(
            chains=chains,
            residues=residues,
            coordinates=coordinates
        )

class PAEParser:
    """Parser for PAE (Predicted Aligned Error) files."""
    
    @staticmethod
    def parse_json(file_path: Union[str, Path]) -> np.ndarray:
        """Parse a PAE JSON file and return the PAE matrix."""
        file_path = Path(file_path)
        logger.debug(f"Parsing PAE JSON file: {file_path}")
        
        if file_path.suffix == '.gz':
            with gzip.open(file_path, 'rt') as f:
                data = json.load(f)
        else:
            with open(file_path, 'r') as f:
                data = json.load(f)
        
        # Log the structure of the JSON data
        logger.debug(f"JSON data keys: {data.keys()}")
        
        # Extract PAE matrix from JSON data
        # The PAE matrix might be under a different key in the JSON
        pae_matrix = None
        possible_keys = ['pae_matrix', 'predicted_aligned_error', 'pae']
        
        for key in possible_keys:
            if key in data:
                pae_matrix = np.array(data[key])
                logger.debug(f"Found PAE matrix under key '{key}' with shape {pae_matrix.shape}")
                break
        
        if pae_matrix is None:
            raise ValueError(f"No PAE matrix found in JSON file. Available keys: {list(data.keys())}")
        
        # Ensure the matrix is 2D
        if pae_matrix.ndim == 1:
            # If it's a 1D array, reshape it to a square matrix
            n = int(np.sqrt(len(pae_matrix)))
            pae_matrix = pae_matrix.reshape(n, n)
            logger.debug(f"Reshaped 1D array to {n}x{n} matrix")
        
        logger.debug(f"Final PAE matrix shape: {pae_matrix.shape}")
        return pae_matrix
    
    @staticmethod
    def parse_npz(file_path: Union[str, Path]) -> np.ndarray:
        """Parse a PAE NPZ file and return the PAE matrix."""
        file_path = Path(file_path)
        logger.debug(f"Parsing PAE NPZ file: {file_path}")
        
        if file_path.suffix == '.gz':
            with gzip.open(file_path, 'rb') as f:
                data = np.load(f)
        else:
            data = np.load(str(file_path))
        
        # Log available keys in the NPZ file
        logger.debug(f"NPZ file keys: {data.files}")
        
        # Extract PAE matrix from NPZ data
        pae_matrix = data['pae_matrix']
        logger.debug(f"Loaded PAE matrix with shape {pae_matrix.shape}")
        
        # Ensure the matrix is 2D
        if pae_matrix.ndim == 1:
            # If it's a 1D array, reshape it to a square matrix
            n = int(np.sqrt(len(pae_matrix)))
            pae_matrix = pae_matrix.reshape(n, n)
            logger.debug(f"Reshaped 1D array to {n}x{n} matrix")
        
        logger.debug(f"Final PAE matrix shape: {pae_matrix.shape}")
        return pae_matrix
