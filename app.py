import dash
from dash import html, dcc, Input, Output, State, dash_table
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import os
import pickle
from datetime import datetime
import dash_bio as dashbio

print("Starting app initialization...")

# Initialize the Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)

# Global variables
gwas_df = None
trait_data = {
    'todo_list': [],
    'skipped_list': [],
    'done_list': [],
    'traits': {}  # {trait: {'status': str, 'annotations': list, 'last_updated': str}}
}
genes_df = None

def load_gwas_results():
    """Load GWAS results from TSV file"""
    df = pd.read_csv('input/combined.tsv', sep='\t', dtype={'CHROM': str, 'TRAIT': str})
    return df

def load_gene_annotations():
    """Load gene annotations"""
    
    # Load gene annotations with specific columns
    annot_df = pd.read_csv('annotation/ZmB73.annot.tsv', sep='\t', 
                          dtype={'gene_id': str},
                          low_memory=False)
    annot_df.columns = ['gene_id', 'name', 'type', 'short_description', 'curator_summary', 
                       'computational_description', 'defline']
    
    # Load GFF3 file
    gff_df = pd.read_csv('annotation/ZmB73.gff3', sep='\t', comment='#', header=None,
                        dtype={0: str})
    gff_df.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 
                     'strand', 'phase', 'attributes']
    
    # Filter for genes and extract gene IDs
    gene_df = gff_df[gff_df['feature'] == 'gene'].copy()
    gene_df['gene_id'] = gene_df['attributes'].str.extract(r'Name=([^;]+)')
    
    # Merge annotations with coordinates
    gene_df = gene_df[['chrom', 'start', 'end', 'strand', 'gene_id']].merge(
        annot_df[['gene_id', 'name', 'short_description', 'curator_summary', 
                 'computational_description', 'defline']], 
        on='gene_id', 
        how='left'
    )
    
    return gene_df

def init_trait_data():
    """Initialize trait data structure"""
    global trait_data
    traits = gwas_df['TRAIT'].dropna().unique()
    
    # Sort traits by their minimum p-value
    trait_min_pvals = []
    for trait in traits:
        trait_data = gwas_df[gwas_df['TRAIT'] == trait]
        min_pval = trait_data['PVAL'].min()
        trait_min_pvals.append((trait, min_pval))
    
    # Sort traits by p-value (ascending)
    sorted_traits = [trait for trait, _ in sorted(trait_min_pvals, key=lambda x: x[1])]
    
    # Initialize trait data
    trait_data = {
        'todo_list': sorted_traits,
        'skipped_list': [],
        'done_list': [],
        'traits': {trait: {
            'status': 'Todo',
            'annotations': [],
            'last_updated': datetime.now().isoformat()
        } for trait in sorted_traits}
    }
    save_trait_data()

def save_trait_data():
    """Save trait data to pickle file"""
    with open('trait_data.pkl', 'wb') as f:
        pickle.dump(trait_data, f)

def load_trait_data():
    """Load trait data from pickle file"""
    global trait_data
    if os.path.exists('trait_data.pkl'):
        with open('trait_data.pkl', 'rb') as f:
            old_data = pickle.load(f)
            
            # Check if we need to migrate from old format
            if isinstance(old_data, dict) and 'traits' not in old_data:
                print("Migrating old data format to new format...")
                new_data = {
                    'todo_list': [],
                    'skipped_list': [],
                    'done_list': [],
                    'traits': {}
                }
                
                # Migrate trait status
                if isinstance(old_data, dict) and '_order' in old_data:
                    for status in ['Todo', 'Done', 'Skipped']:
                        new_data[f'{status.lower()}_list'] = old_data['_order'][status]
                        for trait in old_data['_order'][status]:
                            new_data['traits'][trait] = {
                                'status': status,
                                'annotations': [],
                                'last_updated': datetime.now().isoformat()
                            }
                
                # Migrate annotations if they exist
                if os.path.exists('annotations.pkl'):
                    with open('annotations.pkl', 'rb') as f:
                        annotations = pickle.load(f)
                        for trait, ann_list in annotations.items():
                            if trait in new_data['traits']:
                                new_data['traits'][trait]['annotations'] = ann_list
                
                trait_data = new_data
                save_trait_data()
            else:
                trait_data = old_data
    else:
        init_trait_data()

def get_trait_status(status_filter):
    """Get traits filtered by status"""
    try:
        print(f"\nGetting traits with status: {status_filter}")
        
        if status_filter == 'All':
            all_traits = []
            for status in ['todo_list', 'skipped_list', 'done_list']:
                all_traits.extend(trait_data[status])
            return all_traits
        
        # Get traits in their stored order
        traits = trait_data[f'{status_filter.lower()}_list']
        
        # For Done and Skipped lists, show most recent first
        if status_filter in ['Done', 'Skipped']:
            traits = traits[::-1]  # Reverse the list
        
        print(f"Found {len(traits)} traits with status '{status_filter}'")
        print(f"Total traits in database: {sum(len(trait_data[s]) for s in ['todo_list', 'skipped_list', 'done_list'])}")
        print(f"Status counts: Todo={len(trait_data['todo_list'])}, "
              f"Done={len(trait_data['done_list'])}, "
              f"Skipped={len(trait_data['skipped_list'])}")
        
        if len(traits) == 0:
            print("WARNING: No traits found with requested status!")
        
        print("First few filtered traits:", traits[:5])
        return traits
    except Exception as e:
        print(f"Error in get_trait_status: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())
        return []

def update_trait_status(trait, new_status):
    """Update trait status"""
    try:
        print(f"\nUpdating trait status: {trait} -> {new_status}")
        
        # Find current status
        old_status = trait_data['traits'][trait]['status']
        
        # Remove from old status list
        old_list = f'{old_status.lower()}_list'
        if trait in trait_data[old_list]:
            trait_data[old_list].remove(trait)
        
        # Add to new status list
        new_list = f'{new_status.lower()}_list'
        trait_data[new_list].append(trait)
        
        # Update trait data
        trait_data['traits'][trait]['status'] = new_status
        trait_data['traits'][trait]['last_updated'] = datetime.now().isoformat()
        
        print(f"Updated status for {trait} to {new_status}")
        print(f"New status counts: Todo={len(trait_data['todo_list'])}, "
              f"Done={len(trait_data['done_list'])}, "
              f"Skipped={len(trait_data['skipped_list'])}")
        save_trait_data()
    except Exception as e:
        print(f"Error updating trait status: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())

def save_annotation(trait, snp_info, selected_genes):
    """Save annotation"""
    try:
        # Get existing gene IDs for this trait
        existing_gene_ids = {ann['gene_id'] for ann in trait_data['traits'][trait]['annotations']}
        
        # Filter out duplicates
        new_annotations = []
        for gene in selected_genes:
            if gene['gene_id'] not in existing_gene_ids:
                annotation = {
                    'chromosome': snp_info['chrom'],
                    'position': snp_info['pos'],
                    'p_value': snp_info['pval'],
                    'gene_id': gene['gene_id'],
                    'gene_name': gene['name'],
                    'distance': gene['distance'],
                    'short_description': gene['short_description'],
                    'curator_summary': gene['curator_summary'],
                    'computational_description': gene['computational_description'],
                    'defline': gene['defline']
                }
                trait_data['traits'][trait]['annotations'].append(annotation)
                existing_gene_ids.add(gene['gene_id'])  # Add to set to prevent future duplicates
                new_annotations.append(annotation)
                print(f"Added annotation for gene: {gene['gene_id']}")
            else:
                print(f"Skipped duplicate gene: {gene['gene_id']}")

        save_trait_data()
        print(f"Saved {len(new_annotations)} new annotations for trait {trait}")
    except Exception as e:
        print(f"Error in save_annotation: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())

def get_annotations(trait):
    """Get annotations for a trait"""
    return trait_data['traits'][trait]['annotations']

def create_empty_manhattan_plot():
    """Create an empty Manhattan plot"""
    fig = go.Figure()
    fig.update_layout(
        title='No trait selected',
        xaxis_title='Position',
        yaxis_title='-log10(p-value)',
        showlegend=False
    )
    return fig

def create_manhattan_plot(selected_trait):
    """Create Manhattan plot for selected trait"""
    trait_data = gwas_df[gwas_df['TRAIT'] == selected_trait].copy()
    
    if len(trait_data) == 0:
        return create_empty_manhattan_plot()
    
    trait_data['neg_log_pval'] = -np.log10(trait_data['PVAL'])
    
    # Sort chromosomes
    chroms = sorted(trait_data['CHROM'].unique(), key=lambda x: int(x) if x.isdigit() else float('inf'))
    
    # Calculate cumulative positions for x-axis
    cumulative_pos = 0
    chrom_positions = {}
    for chrom in chroms:
        chrom_positions[chrom] = cumulative_pos
        chrom_data = trait_data[trait_data['CHROM'] == chrom]
        if len(chrom_data) > 0:
            cumulative_pos += chrom_data['POS'].max()
    
    fig = go.Figure()
    
    # Add genome-wide significance line
    fig.add_shape(
        type="line",
        x0=0,
        x1=cumulative_pos,
        y0=-np.log10(5e-8),
        y1=-np.log10(5e-8),
        line=dict(
            color="red",
            width=1,
            dash="dash",
        ),
        layer="below"
    )
    
    # Add suggestive significance line
    fig.add_shape(
        type="line",
        x0=0,
        x1=cumulative_pos,
        y0=-np.log10(1e-5),
        y1=-np.log10(1e-5),
        line=dict(
            color="grey",
            width=1,
            dash="dash",
        ),
        layer="below"
    )
    
    # Add alternating background colors for chromosomes
    for i, chrom in enumerate(chroms):
        if i % 2 == 0:
            fig.add_shape(
                type="rect",
                x0=chrom_positions[chrom],
                x1=chrom_positions[chrom] + trait_data[trait_data['CHROM'] == chrom]['POS'].max(),
                y0=0,
                y1=trait_data['neg_log_pval'].max() * 1.1,
                fillcolor="rgba(0,0,0,0.05)",
                line=dict(width=0),
                layer="below"
            )
    
    # Add scatter points for each chromosome
    for chrom in chroms:
        chrom_data = trait_data[trait_data['CHROM'] == chrom]
        if len(chrom_data) > 0:
            x_pos = chrom_data['POS'] + chrom_positions[chrom]
            
            fig.add_trace(go.Scatter(
                x=x_pos,
                y=chrom_data['neg_log_pval'],
                mode='markers',
                name=f'Chr{chrom}',
                marker=dict(
                    size=8,
                    color='rgb(31, 119, 180)' if int(chrom) % 2 == 0 else 'rgb(255, 127, 14)',
                    opacity=0.7
                ),
                customdata=np.column_stack((
                    chrom_data['CHROM'],
                    chrom_data['POS'],
                    chrom_data['PVAL']
                )),
                hovertemplate=(
                    "Chromosome: %{customdata[0]}<br>"
                    "Position: %{customdata[1]:,}<br>"
                    "P-value: %{customdata[2]:.2e}<br>"
                    "<extra></extra>"
                )
            ))
    
    # Add chromosome labels
    chrom_labels = []
    chrom_positions_list = []
    for chrom in chroms:
        chrom_data = trait_data[trait_data['CHROM'] == chrom]
        if len(chrom_data) > 0:
            mid_point = chrom_positions[chrom] + (chrom_data['POS'].max() / 2)
            chrom_labels.append(f'Chr{chrom}')
            chrom_positions_list.append(mid_point)
    
    fig.update_layout(
        xaxis=dict(
            title='Chromosome',
            ticktext=chrom_labels,
            tickvals=chrom_positions_list,
            showgrid=False
        ),
        yaxis=dict(
            title='-log10(p-value)',
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)'
        ),
        showlegend=False,
        clickmode='event+select',
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=300,
        margin=dict(l=60, r=30, t=20, b=50)
    )
    
    return fig

def create_nearby_genes_table(snp_info):
    """Create table of nearby genes"""
    if snp_info is None:
        return dash_table.DataTable(
            id='genes-table',
            columns=[
                {'name': 'Gene ID', 'id': 'gene_id'},
                {'name': 'Name', 'id': 'name'},
                {'name': 'Location', 'id': 'location'},
                {'name': 'Distance (bp)', 'id': 'distance'},
                {'name': 'Short Description', 'id': 'short_description'},
                {'name': 'Curator Summary', 'id': 'curator_summary'},
                {'name': 'Computational Description', 'id': 'computational_description'},
                {'name': 'Defline', 'id': 'defline'}
            ],
            data=[],
            row_selectable='multi',
            selected_rows=[],
            page_size=10,
            style_table={'overflowX': 'auto'},
            style_cell={
                'textAlign': 'left',
                'padding': '8px',
                'whiteSpace': 'normal',
                'height': 'auto',
                'fontSize': '0.9rem',
                'backgroundColor': '#ffffff',
                'color': '#2c3e50'
            },
            style_header={
                'backgroundColor': '#6495ED',
                'color': 'white',
                'fontWeight': 'bold',
                'fontSize': '0.9rem'
            },
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': '#f8f9fa'
                }
            ],
            persistence=True,
            persistence_type='local'
        )
    
    chrom = str(snp_info['CHROM'])
    pos = int(snp_info['POS'])
    
    # Find nearby genes
    nearby_genes = genes_df[
        (genes_df['chrom'] == chrom) & 
        (
            (genes_df['start'].between(pos-100000, pos+100000)) |
            (genes_df['end'].between(pos-100000, pos+100000)) |
            ((genes_df['start'] <= pos) & (genes_df['end'] >= pos))
        )
    ].copy()
    
    if len(nearby_genes) == 0:
        return dash_table.DataTable(
            id='genes-table',
            columns=[
                {'name': 'Gene ID', 'id': 'gene_id'},
                {'name': 'Name', 'id': 'name'},
                {'name': 'Location', 'id': 'location'},
                {'name': 'Distance (bp)', 'id': 'distance'},
                {'name': 'Short Description', 'id': 'short_description'},
                {'name': 'Curator Summary', 'id': 'curator_summary'},
                {'name': 'Computational Description', 'id': 'computational_description'},
                {'name': 'Defline', 'id': 'defline'}
            ],
            data=[],
            row_selectable='multi',
            selected_rows=[],
            page_size=10,
            style_table={'overflowX': 'auto'},
            style_cell={
                'textAlign': 'left',
                'padding': '8px',
                'whiteSpace': 'normal',
                'height': 'auto',
                'fontSize': '0.9rem',
                'backgroundColor': '#ffffff',
                'color': '#2c3e50'
            },
            style_header={
                'backgroundColor': '#6495ED',
                'color': 'white',
                'fontWeight': 'bold',
                'fontSize': '0.9rem'
            },
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': '#f8f9fa'
                }
            ],
            persistence=True,
            persistence_type='local'
        )
    
    # Calculate distances
    nearby_genes['distance'] = nearby_genes.apply(
        lambda row: min(abs(row['start'] - pos), abs(row['end'] - pos)), 
        axis=1
    )
    
    # Sort by distance
    nearby_genes = nearby_genes.sort_values('distance')
    
    # Create table data
    rows = []
    for _, gene in nearby_genes.iterrows():
        rows.append({
            'gene_id': gene['gene_id'],
            'name': gene['name'],
            'location': f"Chr{chrom} {gene['start']}-{gene['end']} ({gene['strand']})",
            'distance': gene['distance'],
            'short_description': gene['short_description'],
            'curator_summary': gene['curator_summary'],
            'computational_description': gene['computational_description'],
            'defline': gene['defline']
        })
    
    return dash_table.DataTable(
        id='genes-table',
        columns=[
            {'name': 'Gene ID', 'id': 'gene_id'},
            {'name': 'Name', 'id': 'name'},
            {'name': 'Location', 'id': 'location'},
            {'name': 'Distance (bp)', 'id': 'distance'},
            {'name': 'Short Description', 'id': 'short_description'},
            {'name': 'Curator Summary', 'id': 'curator_summary'},
            {'name': 'Computational Description', 'id': 'computational_description'},
            {'name': 'Defline', 'id': 'defline'}
        ],
        data=rows,
        row_selectable='multi',
        selected_rows=[],
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'padding': '8px',
            'whiteSpace': 'normal',
            'height': 'auto',
            'fontSize': '0.9rem',
            'backgroundColor': '#ffffff',
            'color': '#2c3e50'
        },
        style_header={
            'backgroundColor': '#6495ED',
            'color': 'white',
            'fontWeight': 'bold',
            'fontSize': '0.9rem'
        },
        style_data_conditional=[
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': '#f8f9fa'
            }
        ],
        persistence=True,
        persistence_type='local'
    )

def update_annotations_table(selected_trait):
    """Update the annotations table for the selected trait"""
    if not selected_trait or selected_trait == "No traits available":
        return "No trait selected"
    
    annotations_list = get_annotations(selected_trait)
    
    if not annotations_list:
        return "No annotations for this trait"
    
    return dash_table.DataTable(
        id='annotations-table',
        columns=[
            {'name': 'Chromosome', 'id': 'chromosome'},
            {'name': 'Position', 'id': 'position'},
            {'name': 'P-value', 'id': 'p_value'},
            {'name': 'Gene ID', 'id': 'gene_id'},
            {'name': 'Gene Name', 'id': 'gene_name'},
            {'name': 'Distance', 'id': 'distance'},
            {'name': 'Short Description', 'id': 'short_description'},
            {'name': 'Curator Summary', 'id': 'curator_summary'},
            {'name': 'Computational Description', 'id': 'computational_description'},
            {'name': 'Defline', 'id': 'defline'}
        ],
        data=annotations_list,
        row_selectable='multi',
        selected_rows=[],
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'padding': '8px',
            'whiteSpace': 'normal',
            'height': 'auto',
            'fontSize': '0.9rem',
            'backgroundColor': '#ffffff',
            'color': '#2c3e50'
        },
        style_header={
            'backgroundColor': '#6495ED',  # Cornflower blue
            'color': 'white',
            'fontWeight': 'bold',
            'fontSize': '0.9rem'
        },
        style_data_conditional=[
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': '#f8f9fa'
            }
        ]
    )

# Initialize data
print("Starting data initialization...")
try:
    gwas_df = load_gwas_results()
    print(f"Loaded {len(gwas_df)} GWAS results")
    print(f"Found {len(gwas_df['TRAIT'].unique())} unique traits")
    
    genes_df = load_gene_annotations()
    print(f"Loaded {len(genes_df)} gene annotations")
    
    load_trait_data()
except Exception as e:
    print(f"Error during initialization: {str(e)}")
    import traceback
    print(traceback.format_exc())
    gwas_df = pd.DataFrame(columns=['CHROM', 'POS', 'PVAL', 'TRAIT'])
    genes_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'strand', 'gene_id', 
                                   'name', 'short_description', 'curator_summary',
                                   'computational_description', 'defline'])

# App layout
app.layout = html.Div([
    dbc.Container([
        # Main title
        html.H1("SNPDASH", className="text-center my-4", style={
            'color': '#2c3e50',
            'fontWeight': 'bold',
            'fontSize': '2.5rem'
        }),
        
        # Toast container for notifications
        dbc.Toast(
            id="toast",
            header="Warning",
            is_open=False,
            dismissable=True,
            icon="warning",
            style={"position": "fixed", "top": 66, "right": 10, "width": 350, "zIndex": 1000}
        ),
        
        # Status filter and trait selector
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H5("Filter by Status", style={'color': '#34495e', 'fontSize': '1.1rem'}),
                    dcc.Dropdown(
                        id='status-filter',
                        options=[
                            {'label': 'All', 'value': 'All'},
                            {'label': 'Todo', 'value': 'Todo'},
                            {'label': 'Skipped', 'value': 'Skipped'},
                            {'label': 'Done', 'value': 'Done'}
                        ],
                        value='Todo',
                        clearable=False,
                        style={'fontSize': '0.9rem'}
                    ),
                    html.Hr(style={'margin': '10px 0'}),
                    html.Div(id='status-counter', style={
                        'fontSize': '0.9rem',
                        'color': '#666',
                        'textAlign': 'center'
                    })
                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
                })
            ], width=12, style={'margin': '0 5%', 'maxWidth': '90%'})
        ]),
        html.Br(),

        # Current trait display
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H5("Current Trait", style={'color': '#34495e', 'fontSize': '1.1rem'}),
                    html.Div(id='current-trait-display', style={
                        'fontSize': '1rem',
                        'color': '#2c3e50',
                        'padding': '10px',
                        'backgroundColor': '#f8f9fa',
                        'borderRadius': '5px'
                    })
                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
                })
            ], width=12, style={'margin': '0 5%', 'maxWidth': '90%'})
        ]),
        html.Br(),
        
        # Manhattan plot
        dbc.Row([
            dbc.Col([
                html.Div([
                    dcc.Graph(id='manhattan-plot', style={'height': '300px'})
                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
                })
            ], width=12, style={'margin': '0 5%', 'maxWidth': '90%'})
        ]),
        
        # Navigation buttons
        dbc.Row([
            dbc.Col([
                html.Div([
                    # First row of buttons
                    dbc.Row([
                        dbc.Col([
                            dbc.Button("Previous Trait", id="prev-trait-btn", 
                                      color="secondary", className="me-2 mb-2", style={'fontSize': '0.9rem', 'width': '150px'}),
                            dbc.Button("Next Trait", id="next-trait-btn", 
                                      color="success", className="me-2 mb-2", style={'fontSize': '0.9rem', 'width': '150px'}),
                            dbc.Button("Export to CSV", id="export-csv-btn",
                                      color="info", className="me-2 mb-2", style={'fontSize': '0.9rem', 'width': '150px'})
                        ], className="text-center"),
                    ]),
                    # Second row of buttons
                    dbc.Row([
                        dbc.Col([
                            dbc.Button("Done Trait", id="done-trait-btn",
                                      color="success", className="me-2 mb-2", style={'fontSize': '0.9rem', 'width': '150px'}),
                            dbc.Button("Skip Trait", id="skip-trait-btn",
                                      color="warning", className="me-2 mb-2", style={'fontSize': '0.9rem', 'width': '150px'}),
                            dbc.Button("Mark Undone", id="mark-undone-btn",
                                      color="info", className="me-2 mb-2", style={'fontSize': '0.9rem', 'width': '150px'})
                        ], className="text-center"),
                    ])
                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
                })
            ], width=12, className="text-center mt-3", style={'margin': '0 5%', 'maxWidth': '90%'})
        ]),
        html.Br(),

        # Tables section
        dbc.Row([
            # Nearby genes table
            dbc.Col([
                html.Div([
                    dbc.Row([
                        dbc.Col([
                            html.H5("Nearby Genes", style={'color': '#34495e', 'fontSize': '1.1rem', 'marginLeft': '15px', 'display': 'inline-block'}),
                            dbc.Button("Save Genes", id="save-genes-btn",
                                      color="primary", className="ms-2", 
                                      style={'fontSize': '0.7rem', 'width': '80px', 'height': '25px', 'display': 'inline-block'})
                        ], width=12)
                    ]),
                    html.Div(id='nearby-genes-table')
                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
                })
            ], width=12, style={'margin': '0 5%', 'maxWidth': '90%', 'marginBottom': '20px'}),
            
            # Trait annotations table
            dbc.Col([
                html.Div([
                    dbc.Row([
                        dbc.Col([
                            html.H5("Trait Annotations", style={'color': '#34495e', 'fontSize': '1.1rem', 'marginLeft': '15px', 'display': 'inline-block'}),
                            dbc.Button("Clear Annotations", id="clear-annotations-btn",
                                      color="danger", className="ms-2", 
                                      style={'fontSize': '0.7rem', 'width': '120px', 'height': '25px', 'display': 'inline-block'})
                        ], width=12)
                    ]),
                    html.Div(id='trait-annotations-table')
                ], style={
                    'padding': '15px',
                    'backgroundColor': '#ffffff',
                    'borderRadius': '8px',
                    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
                })
            ], width=12, style={'margin': '0 5%', 'maxWidth': '90%'})
        ]),
        
        # Store selected genes and current SNP
        dcc.Store(id='selected-genes-store', data={'selected_rows': [], 'data': []}),
        dcc.Store(id='current-snp-store')
    ], fluid=True, style={'backgroundColor': '#f8f9fa', 'padding': '20px'})
])

@app.callback(
    [Output('current-trait-display', 'children'),
     Output('manhattan-plot', 'figure'),
     Output('manhattan-plot', 'selectedData'),
     Output('nearby-genes-table', 'children'),
     Output('current-snp-store', 'data'),
     Output('trait-annotations-table', 'children'),
     Output('status-counter', 'children'),
     Output("toast", "is_open"),
     Output("toast", "children")],
    [Input('next-trait-btn', 'n_clicks'),
     Input('prev-trait-btn', 'n_clicks'),
     Input('skip-trait-btn', 'n_clicks'),
     Input('mark-undone-btn', 'n_clicks'),
     Input('done-trait-btn', 'n_clicks'),
     Input('save-genes-btn', 'n_clicks'),
     Input('clear-annotations-btn', 'n_clicks'),
     Input('export-csv-btn', 'n_clicks'),
     Input('status-filter', 'value'),
     Input('manhattan-plot', 'selectedData')],
    [State('current-trait-display', 'children'),
     State('nearby-genes-table', 'children'),
     State('current-snp-store', 'data')],
    prevent_initial_call=False
)
def unified_callback(n_clicks_next, n_clicks_prev, n_clicks_skip, n_clicks_undone, n_clicks_done,
                    n_clicks_save, n_clicks_clear, n_clicks_export, status_filter, selected_data,
                    current_trait, genes_table, current_snp):
    try:
        print("\nUnified callback triggered...")
        ctx = dash.callback_context
        if not ctx.triggered:
            # Initial load - return first trait
            traits = get_trait_status(status_filter)
            if not traits:
                return "No traits available", create_empty_manhattan_plot(), None, "No trait selected", None, "No trait selected", get_status_counter_text(), False, ""
            
            # Get the first trait and its best SNP
            first_trait = traits[0]
            trait_df = gwas_df[gwas_df['TRAIT'] == first_trait]
            best_snp = trait_df.loc[trait_df['PVAL'].idxmin()]
            
            # Create initial SNP selection
            initial_snp = {
                'chrom': best_snp['CHROM'],
                'pos': best_snp['POS'],
                'pval': best_snp['PVAL'],
                'trait': first_trait
            }
            
            # Create initial selected data for the plot
            initial_selected_data = {
                'points': [{
                    'customdata': [best_snp['CHROM'], best_snp['POS'], best_snp['PVAL']]
                }]
            }
            
            return (first_trait, 
                   create_manhattan_plot(first_trait), 
                   initial_selected_data,
                   create_nearby_genes_table({'CHROM': best_snp['CHROM'], 'POS': best_snp['POS'], 'PVAL': best_snp['PVAL']}),
                   initial_snp,
                   update_annotations_table(first_trait),
                   get_status_counter_text(),
                   False, "")

        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        print(f"Button triggered: {button_id}")

        # Extract selected rows and data from the table
        selected_rows = []
        genes_data = []
        if isinstance(genes_table, dict) and 'props' in genes_table:
            selected_rows = genes_table['props'].get('selected_rows', [])
            genes_data = genes_table['props'].get('data', [])
        print(f"Selected rows: {selected_rows}")
        print(f"Genes data: {genes_data}")

        # Handle status filter change
        if button_id == 'status-filter':
            traits = get_trait_status(status_filter)
            if not traits:
                return "No traits available", create_empty_manhattan_plot(), None, "No trait selected", None, "No trait selected", get_status_counter_text(), False, ""
            
            # Get the first trait and its best SNP
            first_trait = traits[0]
            trait_df = gwas_df[gwas_df['TRAIT'] == first_trait]
            best_snp = trait_df.loc[trait_df['PVAL'].idxmin()]
            
            # Create initial SNP selection
            initial_snp = {
                'chrom': best_snp['CHROM'],
                'pos': best_snp['POS'],
                'pval': best_snp['PVAL'],
                'trait': first_trait
            }
            
            # Create initial selected data for the plot
            initial_selected_data = {
                'points': [{
                    'customdata': [best_snp['CHROM'], best_snp['POS'], best_snp['PVAL']]
                }]
            }
            
            return (first_trait, 
                   create_manhattan_plot(first_trait), 
                   initial_selected_data,
                   create_nearby_genes_table({'CHROM': best_snp['CHROM'], 'POS': best_snp['POS'], 'PVAL': best_snp['PVAL']}),
                   initial_snp,
                   update_annotations_table(first_trait),
                   get_status_counter_text(),
                   False, "")

        # Handle Manhattan plot selection
        if button_id == 'manhattan-plot':
            if not current_trait or current_trait == "No traits available":
                return dash.no_update, dash.no_update, dash.no_update, "No trait selected", None, "No trait selected", dash.no_update, False, ""
            
            if not selected_data or not selected_data['points']:
                trait_df = gwas_df[gwas_df['TRAIT'] == current_trait]
                if len(trait_df) == 0:
                    return dash.no_update, dash.no_update, dash.no_update, "No data available for this trait", None, dash.no_update, dash.no_update, False, ""
                
                best_snp = trait_df.loc[trait_df['PVAL'].idxmin()]
                point = {
                    'customdata': [best_snp['CHROM'], best_snp['POS'], best_snp['PVAL']]
                }
            else:
                point = selected_data['points'][0]
            
            chrom = point['customdata'][0]
            pos = point['customdata'][1]
            pval = point['customdata'][2]
            
            current_snp = {
                'chrom': chrom,
                'pos': pos,
                'pval': pval,
                'trait': current_trait
            }
            
            genes_table = create_nearby_genes_table({
                'CHROM': chrom,
                'POS': pos,
                'PVAL': pval
            })
            
            return dash.no_update, dash.no_update, dash.no_update, genes_table, current_snp, update_annotations_table(current_trait), dash.no_update, False, ""

        # Handle trait navigation and status updates
        traits = get_trait_status(status_filter)
        if not traits:
            return "No traits available", create_empty_manhattan_plot(), None, "No trait selected", None, "No trait selected", get_status_counter_text(), False, ""

        if not current_trait or current_trait not in traits:
            return traits[0], create_manhattan_plot(traits[0]), None, "No trait selected", None, "No trait selected", get_status_counter_text(), False, ""

        current_index = traits.index(current_trait)
        new_trait = current_trait
        show_toast = False
        toast_message = ""

        # Handle different button actions
        if button_id == 'next-trait-btn':
            next_index = (current_index + 1) % len(traits)
            new_trait = traits[next_index]
        elif button_id == 'prev-trait-btn':
            prev_index = (current_index - 1) % len(traits)
            new_trait = traits[prev_index]
        elif button_id == 'skip-trait-btn':
            update_trait_status(current_trait, 'Skipped')
            traits = get_trait_status(status_filter)
            new_trait = traits[0] if traits else "No traits available"
        elif button_id == 'mark-undone-btn':
            if current_trait in trait_data['traits']:
                trait_data['traits'][current_trait]['annotations'] = []
                save_trait_data()
            update_trait_status(current_trait, 'Todo')
            traits = get_trait_status(status_filter)
            new_trait = traits[0] if traits else "No traits available"
        elif button_id == 'done-trait-btn':
            print(f"Done button clicked. Selected rows: {selected_rows}, Genes data: {genes_data}")
            if not selected_rows or not genes_data or not current_snp:
                show_toast = True
                toast_message = "Please select at least one gene before marking the trait as done."
                return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, show_toast, toast_message
            
            selected_genes = []
            for row_idx in selected_rows:
                if row_idx < len(genes_data):
                    gene = genes_data[row_idx]
                    selected_genes.append({
                        'gene_id': gene['gene_id'],
                        'name': gene['name'],
                        'distance': gene['distance'],
                        'short_description': gene['short_description'],
                        'curator_summary': gene['curator_summary'],
                        'computational_description': gene['computational_description'],
                        'defline': gene['defline']
                    })
            
            if selected_genes:
                save_annotation(current_trait, current_snp, selected_genes)
                save_trait_data()  # Make sure to save after adding annotations
            
            update_trait_status(current_trait, 'Done')
            traits = get_trait_status(status_filter)
            new_trait = traits[0] if traits else "No traits available"
        elif button_id == 'save-genes-btn':
            print(f"Save button clicked. Selected rows: {selected_rows}, Genes data: {genes_data}")
            if not selected_rows or not genes_data or not current_snp:
                show_toast = True
                toast_message = "Please select at least one gene to save."
                return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, show_toast, toast_message
            
            selected_genes = []
            for row_idx in selected_rows:
                if row_idx < len(genes_data):
                    gene = genes_data[row_idx]
                    selected_genes.append({
                        'gene_id': gene['gene_id'],
                        'name': gene['name'],
                        'distance': gene['distance'],
                        'short_description': gene['short_description'],
                        'curator_summary': gene['curator_summary'],
                        'computational_description': gene['computational_description'],
                        'defline': gene['defline']
                    })
            
            if selected_genes:
                save_annotation(current_trait, current_snp, selected_genes)
                save_trait_data()  # Make sure to save after adding annotations
                show_toast = False
                toast_message = f"Saved {len(selected_genes)} genes for {current_trait}"
            
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, update_annotations_table(current_trait), dash.no_update, show_toast, toast_message
        elif button_id == 'clear-annotations-btn':
            if current_trait in trait_data['traits']:
                trait_data['traits'][current_trait]['annotations'] = []
                save_trait_data()
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, update_annotations_table(current_trait), dash.no_update, False, ""
        elif button_id == 'export-csv-btn':
            rows = []
            for trait, trait_info in trait_data['traits'].items():
                for ann in trait_info['annotations']:
                    row = ann.copy()
                    row['trait'] = trait
                    rows.append(row)
            
            if rows:
                df = pd.DataFrame(rows)
                cols = ['trait'] + [col for col in df.columns if col != 'trait']
                df = df[cols]
                df.to_csv('annotations.csv', index=False)
                show_toast = True
                toast_message = f"Exported {len(rows)} annotations to annotations.csv"
            
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, show_toast, toast_message

        # Update all outputs based on the new trait
        if new_trait == "No traits available":
            return new_trait, create_empty_manhattan_plot(), None, "No trait selected", None, "No trait selected", get_status_counter_text(), show_toast, toast_message
        
        # Get the best SNP for the new trait
        trait_df = gwas_df[gwas_df['TRAIT'] == new_trait]
        best_snp = trait_df.loc[trait_df['PVAL'].idxmin()]
        
        # Create SNP selection
        new_snp = {
            'chrom': best_snp['CHROM'],
            'pos': best_snp['POS'],
            'pval': best_snp['PVAL'],
            'trait': new_trait
        }
        
        # Create selected data for the plot
        selected_data = {
            'points': [{
                'customdata': [best_snp['CHROM'], best_snp['POS'], best_snp['PVAL']]
            }]
        }
        
        return (new_trait, 
               create_manhattan_plot(new_trait), 
               selected_data,
               create_nearby_genes_table({'CHROM': best_snp['CHROM'], 'POS': best_snp['POS'], 'PVAL': best_snp['PVAL']}),
               new_snp,
               update_annotations_table(new_trait),
               get_status_counter_text(),
               show_toast, 
               toast_message)

    except Exception as e:
        print(f"Error in unified_callback: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())
        return "Error occurred", create_empty_manhattan_plot(), None, "Error occurred", None, "Error occurred", get_status_counter_text(), True, f"An error occurred: {str(e)}"

def get_status_counter_text():
    """Helper function to get the status counter text"""
    todo_count = len(trait_data['todo_list'])
    done_count = len(trait_data['done_list'])
    skipped_count = len(trait_data['skipped_list'])
    total_count = todo_count + done_count + skipped_count
    
    return f"Total: {total_count} | Todo: {todo_count} | Done: {done_count} | Skipped: {skipped_count}"

if __name__ == '__main__':
    app.run_server(debug=True) 