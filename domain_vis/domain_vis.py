import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
import requests
import io
import py3Dmol
import numpy as np

def load_alphafold_structure(chain_id):
    url = f"https://alphafold.ebi.ac.uk/files/{chain_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure("protein", io.StringIO(response.text))
        return structure
    else:
        raise ValueError(f"Failed to download structure for {chain_id}")

def parse_boundaries(boundary_string):
    boundaries = []
    for segment in boundary_string.split('_'):
        boundaries.extend([tuple(map(int, pair.split('-'))) for pair in segment.split(',')])
    return boundaries

def load_domain_boundaries(csv_file, chain_id):
    df = pd.read_csv(csv_file)
    protein_row = df[df['chain_id'] == chain_id].iloc[0]
    
    hmm_boundaries = parse_boundaries(protein_row['hmm_boundaries'])
    foldseek_boundaries = parse_boundaries(protein_row['foldseek_boundaries'])
    
    return {
        'hmm_boundaries': hmm_boundaries,
        'foldseek_boundaries': foldseek_boundaries
    }

def assign_colors(structure, domain_boundaries):
    colors = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # Check if it's a standard amino acid
                    resnum = residue.id[1]
                    if any(start <= resnum <= end for start, end in domain_boundaries['hmm_boundaries']) and any(start <= resnum <= end for start, end in domain_boundaries['foldseek_boundaries']):
                        colors[resnum] = "#776bcd"  # Soft purple
                    elif any(start <= resnum <= end for start, end in domain_boundaries['hmm_boundaries']):
                        colors[resnum] = "#22a7f0"  # Soft blue
                    elif any(start <= resnum <= end for start, end in domain_boundaries['foldseek_boundaries']):
                        colors[resnum] = "#e14b31"  # Soft red
                    else:
                        colors[resnum] = "#E0E0E0"  # Light gray
    return colors

def visualize_structure(structure, colors):
    viewer = py3Dmol.view(width=800, height=600)
    pdb_io = io.StringIO()
    pdb_writer = PDBIO()
    pdb_writer.set_structure(structure)
    pdb_writer.save(pdb_io)
    pdb_string = pdb_io.getvalue()
    viewer.addModel(pdb_string, "pdb")
    viewer.setStyle({'cartoon': {'colorscheme': {'prop': 'resi', 'map': colors}}})
    viewer.zoomTo()
    
    print("Legend:")
    print("HMM only: Blue")
    print("TED only: Red")
    print("HMM and TED overlap: Purple")
    print("No domain: Grey")
    
    return viewer

def save_visualization(viewer, output_file):
    html = viewer.render()
    with open(output_file, "w") as f:
        f.write(html)

def vis_domain(df, index):
    row = df.iloc[index]
    chain_id = row['chain_id']
    print(f"Visualizing protein with chain ID: {chain_id}")

    structure = load_alphafold_structure(chain_id)
    
    domain_boundaries = {
        'hmm_boundaries': parse_boundaries(row['hmm_boundaries']),
        'foldseek_boundaries': parse_boundaries(row['foldseek_boundaries'])
    }
    
    hmm_evalue = row['hmm_evalue']
    print(f"HMM E-value: {hmm_evalue}")
    overlap_percentage = row['overlap_percentage']
    print(f"Overlap Percentage: {overlap_percentage}")

    colors = assign_colors(structure, domain_boundaries)
    
    viewer = visualize_structure(structure, colors)
    
    viewer.show()
