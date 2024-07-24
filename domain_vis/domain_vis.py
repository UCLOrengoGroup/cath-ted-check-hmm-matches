import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
import requests
import io
import py3Dmol
import numpy as np
import ipywidgets as widgets
from IPython.display import display, clear_output, HTML
import plotly.graph_objects as go
from ipywidgets import Text, Button, ToggleButtons
import plotly.express as px

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
    
    legend_html = """
    <div style="font-family: Arial, sans-serif; font-size: 14px; margin-bottom: 10px;">
        <strong>Legend:</strong>
        <span style="color: #22a7f0; margin-right: 10px;">■ HMM only: Blue</span>
        <span style="color: #e14b31; margin-right: 10px;">■ TED only: Red</span>
        <span style="color: #776bcd; margin-right: 10px;">■ HMM and TED overlap: Purple</span>
        <span style="color: #A9A9A9;">■ No domain: Grey</span>
    """
    display(HTML(legend_html))
    
    return viewer

def save_visualization(viewer, output_file):
    html = viewer.render()
    with open(output_file, "w") as f:
        f.write(html)

def vis_domain(df, index):
    row = df.iloc[index]
    chain_id = row['chain_id']
    hmm_evalue = row['hmm_evalue']
    overlap_percentage = row['overlap_percentage']
    plddt = row['plddt']
    packing_density = row['packing_density']
    
    structure = load_alphafold_structure(chain_id)
    
    domain_boundaries = {
        'hmm_boundaries': parse_boundaries(row['hmm_boundaries']),
        'foldseek_boundaries': parse_boundaries(row['foldseek_boundaries'])
    }
    
    colors = assign_colors(structure, domain_boundaries)
    
    viewer = visualize_structure(structure, colors)
    
    return viewer, chain_id, hmm_evalue, overlap_percentage, plddt, packing_density

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
    
    # Search box and button
    search_box = Text(placeholder='Enter Uniprot ID or chain_id', description='Search:', continuous_update=False)
    search_button = Button(description='Go')
    
    # Add toggle button for coloring
    color_toggle = widgets.ToggleButtons(
        options=['pLDDT', 'Packing Density'],
        description='Color by:',
        disabled=False,
        button_style='',
    )

    # Add toggle button for graph type
    graph_toggle = ToggleButtons(
        options=['HMM E-value vs Overlap', 'Overlap vs pLDDT'],
        description='Graph:',
        disabled=False,
        button_style='',
    )

    # Create initial scatter plot
    fig = go.FigureWidget()

    def update_scatter_plot(color_by, graph_type, selected_index=None):
        with fig.batch_update():
            if graph_type == 'HMM E-value vs Overlap':
                x_data = df['hmm_evalue']
                y_data = df['overlap_percentage']
                x_title = 'HMM E-value'
                y_title = 'Overlap Percentage'
                x_type = 'log'
            else:  # 'Overlap vs pLDDT'
                x_data = df['overlap_percentage']
                y_data = df['plddt']
                x_title = 'Overlap Percentage'
                y_title = 'pLDDT'
                x_type = 'linear'

            if color_by == 'pLDDT':
                colorscale = px.colors.sequential.Viridis
                color_values = df['plddt']
                colorbar_title = 'pLDDT'
            else:
                colorscale = px.colors.sequential.Plasma
                color_values = df['packing_density']
                colorbar_title = 'Packing Density'

            fig.data[0].x = x_data
            fig.data[0].y = y_data
            fig.data[0].marker.color = color_values
            fig.data[0].marker.colorscale = colorscale
            fig.data[0].marker.colorbar.title = colorbar_title

            fig.update_layout(
                title=f'{x_title} vs {y_title}',
                xaxis_title=x_title,
                yaxis_title=y_title,
                xaxis_type=x_type
            )

            if selected_index is not None:
                selected_point = df.iloc[selected_index]
                fig.data[1].x = [x_data.iloc[selected_index]]
                fig.data[1].y = [y_data.iloc[selected_index]]

    # Initial plot setup
    fig.add_trace(go.Scatter(
        x=df['hmm_evalue'],
        y=df['overlap_percentage'],
        mode='markers',
        marker=dict(
            size=8,
            colorscale=px.colors.sequential.Viridis,
            color=df['plddt'],
            colorbar=dict(title='pLDDT'),
            showscale=True
        ),
        text=df['chain_id'],
        name='All Proteins'
    ))
    fig.add_trace(go.Scatter(
        x=[],
        y=[],
        mode='markers',
        marker=dict(color='black', size=12, symbol='star'),
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
        width=600,
        height=600
    )

    def on_color_toggle(change):
        update_scatter_plot(change['new'], graph_toggle.value, slider.value)

    def on_graph_toggle(change):
        update_scatter_plot(color_toggle.value, change['new'], slider.value)

    color_toggle.observe(on_color_toggle, names='value')
    graph_toggle.observe(on_graph_toggle, names='value')

    def update_visualization(index):
        with output:
            clear_output(wait=True)
            viewer, chain_id, hmm_evalue, overlap_percentage, plddt, packing_density = vis_domain(df, index)
            print(f"Visualizing protein with chain ID: {chain_id}")
            print(f"HMM E-value: {hmm_evalue}")
            print(f"Overlap Percentage: {overlap_percentage}")
            print(f"pLDDT: {plddt}")
            print(f"Packing Density: {packing_density}")
            print(f"Showing protein {index + 1} of {len(df)}")
            viewer.show()
        
        # Update scatter plot
        update_scatter_plot(color_toggle.value, graph_toggle.value, index)

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
    
    def on_search_click(b):
        search_term = search_box.value.strip()
        if search_term:
            # Convert Uniprot ID to chain_id format if necessary
            if not search_term.startswith('AF-'):
                search_term = f'AF-{search_term}-F1-model_v4'
            
            matches = df[df['chain_id'].str.contains(search_term, case=False, regex=False)]
            if not matches.empty:
                index = df.index.get_loc(matches.index[0])
                slider.value = index
            else:
                print(f"No matches found for '{search_term}'")
    
    slider.observe(on_value_change, names='value')
    prev_button.on_click(on_prev_click)
    next_button.on_click(on_next_click)
    sort_dropdown.observe(on_sort_change, names='value')
    search_button.on_click(on_search_click)
    
    # Layout
    controls = widgets.HBox([prev_button, slider, next_button])
    search_controls = widgets.HBox([search_box, search_button])
    left_controls = widgets.VBox([sort_dropdown, controls, search_controls])
    
    # Create a container for the plot and its controls
    plot_controls = widgets.VBox([color_toggle, graph_toggle, fig])
    
    # Create a container for the visualization
    vis_container = widgets.VBox([output])
    
    # Main layout
    main_layout = widgets.VBox([
        left_controls, 
        widgets.HBox([vis_container, plot_controls])
    ])
    
    display(main_layout)
    
    update_visualization(0)