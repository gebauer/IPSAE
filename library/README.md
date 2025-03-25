# IPSAE Library

A Python library for scoring interprotein interactions in AlphaFold2 and AlphaFold3 models.

## Features

- Calculate IPSAE scores for protein-protein interactions
- Support for AlphaFold2 and AlphaFold3 output formats
- Chain-chain and residue-level scoring
- Integration with MolStar for structure visualization

## Installation

```bash
pip install -e .
```

## Usage

Basic usage:

```python
from ipsae.core.calculator import IPSAECalculator

# Initialize calculator
calculator = IPSAECalculator(
    pae_file_path="path/to/pae.json",
    structure_file_path="path/to/structure.cif",
    pae_cutoff=15,
    dist_cutoff=15
)

# Calculate scores
results = calculator.calculate()

# Access results
chain_chain_scores = results.chain_chain_scores
residue_scores = results.residue_scores
structure_data = results.structure_data
```

## Development

To set up the development environment:

1. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

2. Install development dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Run tests:
   ```bash
   pytest
   ```

## Citation

If you use this software in your research, please cite:

Original IPSAE algorithm:
```
Dunbrack, R. (2025). IPSAE: Scoring function for interprotein interactions in AlphaFold2 and AlphaFold3.
bioRxiv 2025.02.10.637595
```

This library:
```
Gebauer, J. and Dunbrack, R. (2025). IPSAE-lib: A comprehensive Python library for protein interaction analysis.
[Publication details pending]
``` 