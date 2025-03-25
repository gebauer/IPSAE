#!/usr/bin/env python3
"""
Script to generate output files from IPSAE JSON results.
This script creates the following files:
1. *_15_15.txt - Chain-chain interaction scores
2. *_15_15_byres.txt - Residue-level scores
3. *_15_15.pml - PyMOL visualization script
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

def load_results(json_file: str) -> Dict:
    """Load IPSAE results from JSON file."""
    with open(json_file, 'r') as f:
        return json.load(f)

def write_chain_chain_scores(results: Dict, output_file: str) -> None:
    """Write chain-chain interaction scores to file."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("\nChn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n")
        
        # Write scores for each chain pair
        for score in results['chain_chain_scores']:
            # Calculate additional metrics
            n0res = score['num_interactions']
            n0chn = len([r for r in results['structure_data']['residues'] if r[0] == score['chain1']])
            n0dom = len([r for r in results['structure_data']['residues'] if r[0] == score['chain2']])
            d0res = 1.24 * (n0res - 15) ** (1/3) - 1.8
            d0chn = 1.24 * (n0chn - 15) ** (1/3) - 1.8
            d0dom = 1.24 * (n0dom - 15) ** (1/3) - 1.8
            
            # Write asymmetric scores
            f.write(f"{score['chain1']}    {score['chain2']}     15   15   asym  {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    0.530    {score['ipsae_score']:.6f}      {score['pdockq_score']:.4f}     {score['pdockq_score']:.4f}    {score['lis_score']:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            f.write(f"{score['chain2']}    {score['chain1']}     15   15   asym  {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    0.530    {score['ipsae_score']:.6f}      {score['pdockq_score']:.4f}     {score['pdockq_score']:.4f}    {score['lis_score']:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0dom}   {n0chn}   {n0dom}   {n0chn}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            # Write max scores
            f.write(f"{score['chain1']}    {score['chain2']}     15   15   max   {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    0.530    {score['ipsae_score']:.6f}      {score['pdockq_score']:.4f}     {score['pdockq_score']:.4f}    {score['lis_score']:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")

def write_residue_scores(results: Dict, output_file: str) -> None:
    """Write residue-level scores to file."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("i   AlignChn ScoredChain  AlignResNum  AlignResType  AlignRespLDDT      n0chn  n0dom  n0res    d0chn     d0dom     d0res   ipTM_pae  ipSAE_d0chn ipSAE_d0dom    ipSAE \n")
        
        # Write scores for each residue
        for i, score in enumerate(results['residue_scores'], 1):
            # Calculate additional metrics
            n0chn = len([r for r in results['structure_data']['residues'] if r[0] == score['chain']])
            n0dom = len([r for r in results['structure_data']['residues'] if r[0] != score['chain']])
            n0res = score['num_interactions']
            d0chn = 1.24 * (n0chn - 15) ** (1/3) - 1.8
            d0dom = 1.24 * (n0dom - 15) ** (1/3) - 1.8
            d0res = 1.24 * (n0res - 15) ** (1/3) - 1.8
            
            # Write residue score
            f.write(f"{i:3d}   {score['chain']:8s} {score['chain']:8s} {score['residue_number']:11d} {score['residue_name']:11s} {15.0:11.2f} {n0chn:11d} {n0dom:7d} {n0res:7d} {d0chn:9.3f} {d0dom:9.3f} {d0res:9.3f} {score['ipsae_score']:9.4f} {score['ipsae_score']:9.4f} {score['ipsae_score']:9.4f} {score['ipsae_score']:9.4f}\n")

def write_pymol_script(results: Dict, output_file: str) -> None:
    """Write PyMOL visualization script."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Chn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS      n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n")
        
        # Write scores and coloring commands for each chain pair
        for score in results['chain_chain_scores']:
            # Calculate additional metrics
            n0res = score['num_interactions']
            n0chn = len([r for r in results['structure_data']['residues'] if r[0] == score['chain1']])
            n0dom = len([r for r in results['structure_data']['residues'] if r[0] == score['chain2']])
            d0res = 1.24 * (n0res - 15) ** (1/3) - 1.8
            d0chn = 1.24 * (n0chn - 15) ** (1/3) - 1.8
            d0dom = 1.24 * (n0dom - 15) ** (1/3) - 1.8
            
            # Write asymmetric scores
            f.write(f"# {score['chain1']}    {score['chain2']}     15   15   asym  {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    0.530    {score['ipsae_score']:.6f}      {score['pdockq_score']:.4f}     {score['pdockq_score']:.4f}    {score['lis_score']:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            
            # Write coloring command
            f.write(f"alias color_{score['chain1']}_{score['chain2']}, color magenta, chain {score['chain1']} and resi {'+'.join(str(r[1]) for r in results['structure_data']['residues'] if r[0] == score['chain1'])}; color marine, chain {score['chain2']} and resi {'+'.join(str(r[1]) for r in results['structure_data']['residues'] if r[0] == score['chain2'])}\n\n")
            
            # Write reverse direction scores
            f.write(f"# {score['chain2']}    {score['chain1']}     15   15   asym  {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    0.530    {score['ipsae_score']:.6f}      {score['pdockq_score']:.4f}     {score['pdockq_score']:.4f}    {score['lis_score']:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0dom}   {n0chn}   {n0dom}   {n0chn}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            # Write max scores
            f.write(f"# {score['chain1']}    {score['chain2']}     15   15   max   {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    {score['ipsae_score']:.6f}    0.530    {score['ipsae_score']:.6f}      {score['pdockq_score']:.4f}     {score['pdockq_score']:.4f}    {score['lis_score']:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            
            # Write coloring command for reverse direction
            f.write(f"alias color_{score['chain2']}_{score['chain1']}, color marine, chain {score['chain2']} and resi {'+'.join(str(r[1]) for r in results['structure_data']['residues'] if r[0] == score['chain2'])}; color magenta, chain {score['chain1']} and resi {'+'.join(str(r[1]) for r in results['structure_data']['residues'] if r[0] == score['chain1'])}\n\n")

def main():
    parser = argparse.ArgumentParser(description='Generate output files from IPSAE JSON results')
    parser.add_argument('--json_file', type=str, required=True,
                      help='Path to the IPSAE results JSON file')
    parser.add_argument('--output_dir', type=str, default='output',
                      help='Directory for output files (default: output)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load results
    results = load_results(args.json_file)
    
    # Generate output files
    base_name = Path(args.json_file).stem.replace('_ipsae_results', '')
    
    # Write chain-chain scores
    write_chain_chain_scores(results, output_dir / f"{base_name}_15_15.txt")
    
    # Write residue-level scores
    write_residue_scores(results, output_dir / f"{base_name}_15_15_byres.txt")
    
    # Write PyMOL script
    write_pymol_script(results, output_dir / f"{base_name}_15_15.pml")
    
    print(f"Generated output files in {output_dir}")

if __name__ == "__main__":
    main() 