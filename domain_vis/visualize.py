import pandas as pd
from IPython.display import display, clear_output, HTML
import plotly.graph_objects as go
from ipywidgets import Text, Button, ToggleButtons
import plotly.express as px
import ipywidgets as widgets

from .core import (
    normalize_uniprot_to_chain_id,
    parse_boundaries,
    fetch_alphafold_pdb,
    load_structure_from_pdb_string,
    residue_color_map,
    build_py3dmol_viewer,
)


def _load_structure(chain_id: str):
    pdb_text = fetch_alphafold_pdb(chain_id)
    return load_structure_from_pdb_string(pdb_text)


def _visualize_structure(structure, colors):
    viewer = build_py3dmol_viewer(structure, colors)

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


def _vis_domain(df: pd.DataFrame, index: int):
    row = df.iloc[index]
    chain_id = row['chain_id']
    hmm_evalue = row.get('hmm_evalue')
    overlap_percentage = row.get('overlap_percentage')
    plddt = row.get('plddt')
    packing_density = row.get('packing_density')

    structure = _load_structure(chain_id)

    hmm_boundaries = parse_boundaries(row.get('hmm_boundaries', ""))
    ted_boundaries = parse_boundaries(row.get('foldseek_boundaries', ""))

    colors = residue_color_map(structure, hmm_boundaries, ted_boundaries)

    viewer = _visualize_structure(structure, colors)

    return viewer, chain_id, hmm_evalue, overlap_percentage, plddt, packing_density


def visualize_domains_and_plot(df: pd.DataFrame):
    output = widgets.Output()

    sort_options = ['Index', 'HMM E-value', 'Overlap Percentage']
    sort_dropdown = widgets.Dropdown(options=sort_options, description='Sort by:')

    slider = widgets.IntSlider(value=0, min=0, max=max(len(df) - 1, 0), step=1, description='Index:')

    prev_button = widgets.Button(description='Previous')
    next_button = widgets.Button(description='Next')

    search_box = Text(placeholder='Enter Uniprot ID or chain_id', description='Search:', continuous_update=False)
    search_button = Button(description='Go')

    has_match_column = 'match' in df.columns

    if has_match_column:
        color_toggle = widgets.ToggleButtons(
            options=['Match/No Match', 'pLDDT', 'Packing Density'],
            description='Color by:',
            disabled=False,
            button_style='',
        )
    else:
        color_toggle = widgets.ToggleButtons(
            options=['pLDDT', 'Packing Density'],
            description='Color by:',
            disabled=False,
            button_style='',
        )

    graph_toggle = ToggleButtons(
        options=['HMM E-value vs Overlap', 'Overlap vs pLDDT'],
        description='Graph:',
        disabled=False,
        button_style='',
    )

    fig = go.FigureWidget()

    def update_scatter_plot(color_by, graph_type, selected_index=None):
        with fig.batch_update():
            if graph_type == 'HMM E-value vs Overlap':
                x_data = df['hmm_evalue']
                y_data = df['overlap_percentage']
                x_title = 'HMM E-value'
                y_title = 'Overlap Percentage'
                x_type = 'log'
            else:
                x_data = df['overlap_percentage']
                y_data = df['plddt']
                x_title = 'Overlap Percentage'
                y_title = 'pLDDT'
                x_type = 'linear'

            if has_match_column and color_by == 'Match/No Match':
                fig.data[0].marker.color = df['match'].map({'match': 'green', 'no match': 'red'})
                fig.data[0].marker.showscale = False
                fig.data[0].marker.colorscale = None
                fig.data[0].name = 'All Proteins'

                # Auxiliary traces for legend
                match_df = df[df['match'] == 'match']
                no_match_df = df[df['match'] == 'no match']
                fig.data[2].x = match_df['hmm_evalue']
                fig.data[2].y = match_df['overlap_percentage']
                fig.data[2].marker.color = 'green'
                fig.data[2].name = 'Match'

                fig.data[3].x = no_match_df['hmm_evalue']
                fig.data[3].y = no_match_df['overlap_percentage']
                fig.data[3].marker.color = 'red'
                fig.data[3].name = 'No Match'
            else:
                if color_by == 'pLDDT':
                    colorscale = px.colors.sequential.Viridis
                    color_values = df['plddt']
                    colorbar_title = 'pLDDT'
                else:
                    colorscale = px.colors.sequential.Plasma
                    color_values = df['packing_density']
                    colorbar_title = 'Packing Density'

                fig.data[0].marker.color = color_values
                fig.data[0].marker.colorscale = colorscale
                fig.data[0].marker.showscale = True
                fig.data[0].marker.colorbar.title = colorbar_title
                fig.data[0].name = 'All Proteins'

                # Hide aux traces
                if len(fig.data) > 3:
                    fig.data[2].x, fig.data[2].y = [], []
                    fig.data[3].x, fig.data[3].y = [], []

            fig.data[0].x = x_data
            fig.data[0].y = y_data

            fig.update_layout(
                title=f'{x_title} vs {y_title}',
                xaxis_title=x_title,
                yaxis_title=y_title,
                xaxis_type=x_type
            )

            if selected_index is not None:
                fig.data[1].x = [x_data.iloc[selected_index]]
                fig.data[1].y = [y_data.iloc[selected_index]]

    # Initial traces
    fig.add_trace(go.Scatter(
        x=df['hmm_evalue'],
        y=df['overlap_percentage'],
        mode='markers',
        marker=dict(
            size=8,
            colorscale=px.colors.sequential.Viridis,
            color=df['plddt'] if 'plddt' in df.columns else None,
            colorbar=dict(title='pLDDT') if 'plddt' in df.columns else None,
            showscale='plddt' in df.columns,
            opacity=0.7,
        ),
        text=df['chain_id'] if 'chain_id' in df.columns else None,
        name='All Proteins'
    ))

    fig.add_trace(go.Scatter(
        x=[],
        y=[],
        mode='markers',
        marker=dict(color='black', size=12, symbol='star'),
        name='Selected Protein'
    ))

    # Aux traces for match/no-match legend positions
    fig.add_trace(go.Scatter(x=[], y=[], mode='markers', marker=dict(color='green', size=8, opacity=0.5), name='Match'))
    fig.add_trace(go.Scatter(x=[], y=[], mode='markers', marker=dict(color='red', size=8, opacity=0.5), name='No Match'))

    fig.update_layout(
        title='HMM E-value vs Overlap Percentage',
        xaxis_title='HMM E-value',
        yaxis_title='Overlap Percentage',
        xaxis_type='log',
        legend=dict(
            yanchor="top",
            y=1.05,
            xanchor="left",
            x=0.01,
            itemsizing='constant',
            bgcolor='rgba(255, 255, 255, 0.5)',
            title='Legend',
            orientation='h'
        ),
        width=700,
        height=650,
        margin=dict(l=50, r=10, t=50, b=50)
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
            viewer, chain_id, hmm_evalue, overlap_percentage, plddt, packing_density = _vis_domain(df, index)
            print(f"Visualizing protein with chain ID: {chain_id}")
            if 'ted_id' in df.columns:
                print(f"TED ID: {df.iloc[index]['ted_id']}")
            if 'HMM-superfamily' in df.columns:
                print(f"HMM SF: {df.iloc[index]['HMM-superfamily']}")
            if 'TED-superfamily' in df.columns:
                print(f"TED SF: {df.iloc[index]['TED-superfamily']}")
            print(f"HMM E-value: {hmm_evalue}")
            print(f"Overlap Percentage: {overlap_percentage}")
            print(f"pLDDT: {plddt}")
            print(f"Packing Density: {packing_density}")
            print(f"Showing protein {index + 1} of {len(df)}")
            viewer.show()

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
            if not search_term.startswith('AF-'):
                search_term = normalize_uniprot_to_chain_id(search_term)

            matches = df[df['chain_id'].astype(str).str.contains(search_term, case=False, regex=False)]
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

    controls = widgets.HBox([prev_button, slider, next_button])
    search_controls = widgets.HBox([search_box, search_button])
    left_controls = widgets.VBox([sort_dropdown, controls, search_controls])

    plot_controls = widgets.VBox([color_toggle, graph_toggle, fig])

    vis_container = widgets.VBox([output])

    main_layout = widgets.VBox([
        left_controls,
        widgets.HBox([vis_container, plot_controls])
    ])

    display(main_layout)

    update_visualization(0)
