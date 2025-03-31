#!/usr/bin/env python3
"""
Test script for the IPSAE file scanner.
This script demonstrates how to use the FileScanner class to find matching
structure and PAE files in a directory.
"""

import argparse
import logging
import sys
from pathlib import Path

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

from scanner import FileScanner

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def main():
    parser = argparse.ArgumentParser(description='Test the IPSAE file scanner')
    parser.add_argument('--root_dir', type=str, required=True,
                      help='Root directory to scan for files')
    parser.add_argument('--structure_patterns', type=str, nargs='+',
                      help='Optional list of patterns for structure files')
    parser.add_argument('--pae_patterns', type=str, nargs='+',
                      help='Optional list of patterns for PAE files')
    
    args = parser.parse_args()
    
    # Initialize scanner
    scanner = FileScanner(
        root_dir=args.root_dir,
        structure_patterns=args.structure_patterns,
        pae_patterns=args.pae_patterns
    )
    
    # Find matches
    matches = scanner.find_matches()
    
    # Print results
    print("\nFound File Matches:")
    print("-" * 50)
    for match in matches:
        print(f"\nModel: {match.model_name}")
        print(f"Seed: {match.seed}")
        print(f"Model Number: {match.model_number}")
        print(f"Structure File: {match.structure_file}")
        print(f"PAE File: {match.pae_file}")
        print("-" * 50)
    
    print(f"\nTotal matches found: {len(matches)}")

if __name__ == "__main__":
    main() 