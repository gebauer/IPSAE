# IPSAE Analysis Scripts

Scripts for batch processing and visualization of IPSAE results.

## Features

- Batch processing of multiple structures
- Directory scanning and file matching
- HTML report generation with MolStar integration
- Interactive visualization of results

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Directory Scanning

```python
from analysis.scripts.scanner import DirectoryScanner

# Initialize scanner
scanner = DirectoryScanner(
    root_dir="path/to/structures",
    pae_cutoff=15,
    dist_cutoff=15
)

# Scan directory
results = scanner.scan()

# Generate HTML report
scanner.generate_report("output.html")
```

### Command Line Interface

```bash
python -m analysis.scripts.scanner --input-dir /path/to/structures --output report.html
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

## Dependencies

- ipsae (local library)
- numpy
- biopython
- molstar
- [other dependencies will be listed here] 