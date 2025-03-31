"""
File scanner module for discovering and matching structure and PAE files.
This module is part of the IPSAE project and is used to find matching pairs
of structure and PAE files in a directory.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class FileMatch:
    """Represents a matched pair of structure and PAE files."""
    structure_file: Path
    pae_file: Path
    model_name: str
    seed: str
    model_number: int
    is_valid: bool = True

class FileScanner:
    """Scans directories for matching structure and PAE files."""
    
    # Default file patterns
    STRUCTURE_PATTERNS = [
        "*_unrelaxed_alphafold2_multimer_v3_model_*_seed_*.pdb",
        "*_unrelaxed_alphafold2_multimer_v3_model_*.pdb",
        "*.pdb"
    ]
    
    PAE_PATTERNS = [
        "*_scores_alphafold2_multimer_v3_model_*_seed_*.json",
        "*_scores_alphafold2_multimer_v3_model_*_seed_*.json.gz",
        "*_scores_alphafold2_multimer_v3_model_*.json",
        "*_scores_alphafold2_multimer_v3_model_*.json.gz"
    ]
    
    def __init__(
        self,
        root_dir: str,
        structure_patterns: Optional[List[str]] = None,
        pae_patterns: Optional[List[str]] = None
    ):
        """
        Initialize the file scanner.
        
        Args:
            root_dir: Root directory to scan for files
            structure_patterns: Optional list of patterns for structure files
            pae_patterns: Optional list of patterns for PAE files
        """
        self.root_dir = Path(root_dir)
        self.structure_patterns = structure_patterns or self.STRUCTURE_PATTERNS
        self.pae_patterns = pae_patterns or self.PAE_PATTERNS
        
        # Initialize file collections
        self.structure_files: Dict[str, Path] = {}
        self.pae_files: Dict[str, Path] = {}
        
        # Scan for files
        self._scan_files()
    
    def _scan_files(self) -> None:
        """Scan the root directory for structure and PAE files."""
        logger.info(f"Scanning directory: {self.root_dir}")
        
        # Find structure files
        for pattern in self.structure_patterns:
            for file_path in self.root_dir.rglob(pattern):
                if self._is_valid_structure_file(file_path):
                    key = self._get_file_key(file_path)
                    if key in self.structure_files:
                        logger.warning(f"Duplicate structure file found: {file_path}")
                    else:
                        self.structure_files[key] = file_path
        
        # Find PAE files
        for pattern in self.pae_patterns:
            for file_path in self.root_dir.rglob(pattern):
                if self._is_valid_pae_file(file_path):
                    key = self._get_file_key(file_path)
                    if key in self.pae_files:
                        logger.warning(f"Duplicate PAE file found: {file_path}")
                    else:
                        self.pae_files[key] = file_path
        
        logger.info(f"Found {len(self.structure_files)} structure files and {len(self.pae_files)} PAE files")
    
    def _is_valid_structure_file(self, file_path: Path) -> bool:
        """
        Validate a structure file.
        
        Args:
            file_path: Path to the structure file
            
        Returns:
            bool: True if the file is valid, False otherwise
        """
        # TODO: Implement actual validation
        return True
    
    def _is_valid_pae_file(self, file_path: Path) -> bool:
        """
        Validate a PAE file.
        
        Args:
            file_path: Path to the PAE file
            
        Returns:
            bool: True if the file is valid, False otherwise
        """
        # TODO: Implement actual validation
        return True
    
    def _get_file_key(self, file_path: Path) -> str:
        """
        Extract a unique key from a file path for matching.
        
        Args:
            file_path: Path to the file
            
        Returns:
            str: A unique key for matching
        """
        # Remove file extension and any compression extensions
        name = file_path.stem
        if name.endswith('.gz'):
            name = name[:-3]
        
        # Extract the base name (everything before the first underscore)
        return name.split('_')[0]
    
    def find_matches(self) -> List[FileMatch]:
        """
        Find matching pairs of structure and PAE files.
        
        Returns:
            List[FileMatch]: List of matched file pairs
        """
        matches = []
        
        # Group files by their base key
        for key in set(self.structure_files.keys()) & set(self.pae_files.keys()):
            structure_file = self.structure_files[key]
            pae_file = self.pae_files[key]
            
            # Extract model information
            model_info = self._extract_model_info(structure_file.name)
            if model_info:
                model_name, seed, model_number = model_info
                matches.append(FileMatch(
                    structure_file=structure_file,
                    pae_file=pae_file,
                    model_name=model_name,
                    seed=seed,
                    model_number=model_number
                ))
        
        logger.info(f"Found {len(matches)} matching pairs")
        return matches
    
    def _extract_model_info(self, filename: str) -> Optional[tuple[str, str, int]]:
        """
        Extract model information from a filename.
        
        Args:
            filename: Name of the file
            
        Returns:
            Optional[tuple[str, str, int]]: Tuple of (model_name, seed, model_number) if found
        """
        try:
            # Example filename: RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000.pdb
            parts = filename.split('_')
            
            # Find model number
            model_idx = parts.index('model')
            model_number = int(parts[model_idx + 1])
            
            # Find seed
            seed_idx = parts.index('seed')
            seed = parts[seed_idx + 1]
            
            # Model name is everything before the first underscore
            model_name = parts[0]
            
            return model_name, seed, model_number
        except (ValueError, IndexError):
            return None 