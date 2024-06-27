import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
import requests
import io
import py3Dmol
import numpy as np
import ipywidgets as widgets
from IPython.display import display, clear_output
import plotly.graph_objects as go

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


def visualize_domains_and_plot(df):
    output = widgets.Output()
    
    # Sorting options
    sort_options = ['Index', 'HMM E-value', 'Overlap Percentage']
    sort_dropdown = widgets.Dropdown(options=sort_options, description='Sort by:')
    
    # Slider
    slider = widgets.IntSlider(value=0, min=0, max=len(df) - 1, step=1, description='Index:')
    
    # Buttons
    prev_button = widgets.Button(description='Previous')
    next_button = widgets.Button(description='Next')
    
    # Create initial scatter plot
    fig = go.FigureWidget()
    fig.add_trace(go.Scatter(
        x=df['hmm_evalue'],
        y=df['overlap_percentage'],
        mode='markers',
        marker=dict(color='blue', size=8),
        text=df['chain_id'],
        name='All Proteins'
    ))
    fig.add_trace(go.Scatter(
        x=[],
        y=[],
        mode='markers',
        marker=dict(color='red', size=12, symbol='star'),
        name='Selected Protein'
    ))
    fig.update_layout(
        title='HMM E-value vs Overlap Percentage',
        xaxis_title='HMM E-value',
        yaxis_title='Overlap Percentage',
        xaxis_type='log',
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        ),
        width=600,  # Increase the width of the plot
        height=600  # Adjust the height if needed
    )
    
    def update_visualization(index):
        with output:
            clear_output(wait=True)
            vis_domain(df, index)
            print(f"Showing protein {index + 1} of {len(df)}")
        
        # Update scatter plot
        selected_point = df.iloc[index]
        with fig.batch_update():
            fig.data[1].x = [selected_point['hmm_evalue']]
            fig.data[1].y = [selected_point['overlap_percentage']]
    
    def on_value_change(change):
        update_visualization(change['new'])
    
    def on_prev_click(b):
        slider.value = max(0, slider.value - 1)
    
    def on_next_click(b):
        slider.value = min(len(df) - 1, slider.value + 1)
    
    def on_sort_change(change):
        nonlocal df
        if change['new'] == 'HMM E-value':
            df = df.sort_values('hmm_evalue')
        elif change['new'] == 'Overlap Percentage':
            df = df.sort_values('overlap_percentage', ascending=False)
        else:
            df = df.sort_index()
        
        slider.max = len(df) - 1
        slider.value = 0
        update_visualization(0)
    
    slider.observe(on_value_change, names='value')
    prev_button.on_click(on_prev_click)
    next_button.on_click(on_next_click)
    sort_dropdown.observe(on_sort_change, names='value')
    
    # Layout
    controls = widgets.HBox([prev_button, slider, next_button])
    top_controls = widgets.VBox([sort_dropdown, controls])
    
    # Create a horizontal layout for PyMOL widget and Plotly plot
    layout = widgets.HBox([output, fig])
    
    display(top_controls, layout)
    
    update_visualization(0)  # Display the first visualization
