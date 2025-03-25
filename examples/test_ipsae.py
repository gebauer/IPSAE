#!/usr/bin/env python3
"""
Test script for the IPSAE library.
This script demonstrates how to use the library to calculate IPSAE scores
for protein-protein interactions using PAE matrices and structure files.
"""

import argparse
from pathlib import Path
from ipsae.core.calculator import IPSAECalculator

def write_chain_chain_scores(results, output_file: str) -> None:
    """Write chain-chain interaction scores to file."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("\nChn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS       n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n")
        
        # Write scores for each chain pair
        for score in results.chain_chain_scores:
            # Calculate additional metrics
            n0res = score.num_interactions
            n0chn = len([r for r in results.structure_data.residues if r[0] == score.chain1])
            n0dom = len([r for r in results.structure_data.residues if r[0] == score.chain2])
            d0res = 1.24 * (n0res - 15) ** (1/3) - 1.8
            d0chn = 1.24 * (n0chn - 15) ** (1/3) - 1.8
            d0dom = 1.24 * (n0dom - 15) ** (1/3) - 1.8
            
            # Write asymmetric scores
            f.write(f"{score.chain1}    {score.chain2}     15   15   asym  {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    0.530    {score.ipsae_score:.6f}      {score.pdockq_score:.4f}     {score.pdockq_score:.4f}    {score.lis_score:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            f.write(f"{score.chain2}    {score.chain1}     15   15   asym  {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    0.530    {score.ipsae_score:.6f}      {score.pdockq_score:.4f}     {score.pdockq_score:.4f}    {score.lis_score:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0dom}   {n0chn}   {n0dom}   {n0chn}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            # Write max scores
            f.write(f"{score.chain1}    {score.chain2}     15   15   max   {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    0.530    {score.ipsae_score:.6f}      {score.pdockq_score:.4f}     {score.pdockq_score:.4f}    {score.lis_score:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")

def write_residue_scores(results, output_file: str) -> None:
    """Write residue-level scores to file."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("i   AlignChn ScoredChain  AlignResNum  AlignResType  AlignRespLDDT      n0chn  n0dom  n0res    d0chn     d0dom     d0res   ipTM_pae  ipSAE_d0chn ipSAE_d0dom    ipSAE \n")
        
        # Write scores for each residue
        for i, score in enumerate(results.residue_scores, 1):
            # Calculate additional metrics
            n0chn = len([r for r in results.structure_data.residues if r[0] == score.chain])
            n0dom = len([r for r in results.structure_data.residues if r[0] != score.chain])
            n0res = score.num_interactions
            d0chn = 1.24 * (n0chn - 15) ** (1/3) - 1.8
            d0dom = 1.24 * (n0dom - 15) ** (1/3) - 1.8
            d0res = 1.24 * (n0res - 15) ** (1/3) - 1.8
            
            # Write residue score
            f.write(f"{i:3d}   {score.chain:8s} {score.chain:8s} {score.residue_number:11d} {score.residue_name:11s} {15.0:11.2f} {n0chn:11d} {n0dom:7d} {n0res:7d} {d0chn:9.3f} {d0dom:9.3f} {d0res:9.3f} {score.ipsae_score:9.4f} {score.ipsae_score:9.4f} {score.ipsae_score:9.4f} {score.ipsae_score:9.4f}\n")

def write_pymol_script(results, output_file: str) -> None:
    """Write PyMOL visualization script."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Chn1 Chn2  PAE Dist  Type   ipSAE    ipSAE_d0chn ipSAE_d0dom  ipTM_af  ipTM_d0chn     pDockQ     pDockQ2    LIS      n0res  n0chn  n0dom   d0res   d0chn   d0dom  nres1   nres2   dist1   dist2  Model\n")
        
        # Write scores and coloring commands for each chain pair
        for score in results.chain_chain_scores:
            # Calculate additional metrics
            n0res = score.num_interactions
            n0chn = len([r for r in results.structure_data.residues if r[0] == score.chain1])
            n0dom = len([r for r in results.structure_data.residues if r[0] == score.chain2])
            d0res = 1.24 * (n0res - 15) ** (1/3) - 1.8
            d0chn = 1.24 * (n0chn - 15) ** (1/3) - 1.8
            d0dom = 1.24 * (n0dom - 15) ** (1/3) - 1.8
            
            # Write asymmetric scores
            f.write(f"# {score.chain1}    {score.chain2}     15   15   asym  {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    0.530    {score.ipsae_score:.6f}      {score.pdockq_score:.4f}     {score.pdockq_score:.4f}    {score.lis_score:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            
            # Write coloring command
            f.write(f"alias color_{score.chain1}_{score.chain2}, color magenta, chain {score.chain1} and resi {'+'.join(str(r[1]) for r in results.structure_data.residues if r[0] == score.chain1)}; color marine, chain {score.chain2} and resi {'+'.join(str(r[1]) for r in results.structure_data.residues if r[0] == score.chain2)}\n\n")
            
            # Write reverse direction scores
            f.write(f"# {score.chain2}    {score.chain1}     15   15   asym  {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    0.530    {score.ipsae_score:.6f}      {score.pdockq_score:.4f}     {score.pdockq_score:.4f}    {score.lis_score:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0dom}   {n0chn}   {n0dom}   {n0chn}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            # Write max scores
            f.write(f"# {score.chain1}    {score.chain2}     15   15   max   {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    {score.ipsae_score:.6f}    0.530    {score.ipsae_score:.6f}      {score.pdockq_score:.4f}     {score.pdockq_score:.4f}    {score.lis_score:.4f}      {n0res}   {n0chn}   {n0dom}   {d0res:.2f}   {d0chn:.2f}   {d0dom:.2f}   {n0chn}   {n0dom}   {n0chn}   {n0dom}   RAF1_KSR1_MEK1_9f755_unrelaxed_alphafold2_multimer_v3_model_1_seed_000\n")
            
            # Write coloring command for reverse direction
            f.write(f"alias color_{score.chain2}_{score.chain1}, color marine, chain {score.chain2} and resi {'+'.join(str(r[1]) for r in results.structure_data.residues if r[0] == score.chain2)}; color magenta, chain {score.chain1} and resi {'+'.join(str(r[1]) for r in results.structure_data.residues if r[0] == score.chain1)}\n\n")

def main():
    parser = argparse.ArgumentParser(description='Calculate IPSAE scores for protein structures')
    parser.add_argument('--pae_file', type=str, required=True,
                      help='Path to the PAE JSON or NPZ file')
    parser.add_argument('--structure_file', type=str, required=True,
                      help='Path to the structure file (CIF or PDB)')
    parser.add_argument('--pae_cutoff', type=float, default=30.0,
                      help='PAE cutoff value (default: 30.0)')
    parser.add_argument('--distance_cutoff', type=float, default=8.0,
                      help='Distance cutoff value in Å (default: 8.0)')
    parser.add_argument('--output_dir', type=str, default='output',
                      help='Directory for output files (default: output)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize calculator
    calculator = IPSAECalculator(
        pae_file=args.pae_file,
        structure_file=args.structure_file,
        pae_cutoff=args.pae_cutoff,
        distance_cutoff=args.distance_cutoff
    )
    
    # Calculate IPSAE scores
    results = calculator.calculate()
    
    # Print chain-chain interaction scores
    print("\nChain-Chain Interaction Scores:")
    print("-" * 50)
    for score in results.chain_chain_scores:
        print(f"Chain {score.chain1} - Chain {score.chain2}:")
        print(f"  IPSAE Score: {score.ipsae_score:.2f}")
        print(f"  pDockQ Score: {score.pdockq_score:.2f}")
        print(f"  LIS Score: {score.lis_score:.2f}")
        print(f"  Number of Interactions: {score.num_interactions}")
        print(f"  Average PAE: {score.avg_pae:.2f}")
        print(f"  Average Distance: {score.avg_distance:.2f} Å")
        print("-" * 50)
    
    # Print residue-level scores
    print("\nResidue-Level Scores:")
    print("-" * 50)
    for score in results.residue_scores:
        print(f"Chain {score.chain} Residue {score.residue_number} ({score.residue_name}):")
        print(f"  IPSAE Score: {score.ipsae_score:.2f}")
        print(f"  Number of Interactions: {score.num_interactions}")
        print(f"  Average PAE: {score.avg_pae:.2f}")
        print(f"  Average Distance: {score.avg_distance:.2f} Å")
        print("-" * 50)
    
    # Save results to JSON
    results.save(output_dir / "ipsae_results.json")
    print(f"\nResults saved to {output_dir}/ipsae_results.json")
    
    # Generate output files
    base_name = Path(args.structure_file).stem
    
    # Write chain-chain scores
    write_chain_chain_scores(results, output_dir / f"{base_name}_15_15.txt")
    
    # Write residue-level scores
    write_residue_scores(results, output_dir / f"{base_name}_15_15_byres.txt")
    
    # Write PyMOL script
    write_pymol_script(results, output_dir / f"{base_name}_15_15.pml")
    
    print(f"\nGenerated output files in {output_dir}")

if __name__ == "__main__":
    main() 