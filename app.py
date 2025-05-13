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
trait_status = {}  # {trait: status}
annotations = {}   # {trait: [list of annotations]}
genes_df = None
trait_last_updated = {}  # {trait: timestamp}

def load_gwas_results():
    """Load GWAS results from TSV file"""
    print("Loading GWAS results...")
    df = pd.read_csv('input/combined.tsv', sep='\t', dtype={'CHROM': str, 'TRAIT': str})
    return df

def load_gene_annotations():
    """Load gene annotations"""
    print("Loading gene annotations...")
    
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

def init_trait_status():
    """Initialize trait status"""
    global trait_status
    traits = gwas_df['TRAIT'].dropna().unique()
    
    # Sort traits by their minimum p-value
    trait_min_pvals = []
    for trait in traits:
        trait_data = gwas_df[gwas_df['TRAIT'] == trait]
        min_pval = trait_data['PVAL'].min()
        trait_min_pvals.append((trait, min_pval))
    
    # Sort traits by p-value (ascending)
    sorted_traits = [trait for trait, _ in sorted(trait_min_pvals, key=lambda x: x[1])]
    
    trait_status = {
        'Todo': {trait: 'Todo' for trait in sorted_traits},
        'Done': {},
        'Skipped': {},
        '_order': {
            'Todo': sorted_traits,
            'Done': [],
            'Skipped': []
        }
    }
    save_trait_status()

def save_trait_status():
    """Save trait status to pickle file"""
    with open('trait_status.pkl', 'wb') as f:
        pickle.dump(trait_status, f)

def load_trait_status():
    """Load trait status from pickle file"""
    global trait_status
    if os.path.exists('trait_status.pkl'):
        with open('trait_status.pkl', 'rb') as f:
            old_status = pickle.load(f)
            
            # Check if we need to migrate from old format
            if isinstance(old_status, dict) and '_order' not in old_status:
                print("Migrating old trait_status format to new format...")
                new_status = {
                    'Todo': {},
                    'Done': {},
                    'Skipped': {},
                    '_order': {
                        'Todo': [],
                        'Done': [],
                        'Skipped': []
                    }
                }
                
                # Migrate traits to new format
                for trait, status in old_status.items():
                    if trait != '_trait_order_list':  # Skip old order list if it exists
                        new_status[status][trait] = status
                        new_status['_order'][status].append(trait)
                
                trait_status = new_status
                save_trait_status()
            else:
                trait_status = old_status
    else:
        init_trait_status()

def save_annotations():
    """Save annotations to pickle file"""
    with open('annotations.pkl', 'wb') as f:
        pickle.dump(annotations, f)

def load_annotations():
    """Load annotations from pickle file"""
    global annotations
    if os.path.exists('annotations.pkl'):
        with open('annotations.pkl', 'rb') as f:
            annotations = pickle.load(f)
    else:
        annotations = {}

def get_trait_status(status_filter):
    """Get traits filtered by status"""
    try:
        print(f"\nGetting traits with status: {status_filter}")
        
        if status_filter == 'All':
            all_traits = []
            for status in ['Todo', 'Done', 'Skipped']:
                all_traits.extend(trait_status['_order'][status])
            return all_traits
        
        # Get traits in their stored order
        traits = trait_status['_order'][status_filter]
        
        # For Done and Skipped lists, show most recent first
        if status_filter in ['Done', 'Skipped']:
            traits = traits[::-1]  # Reverse the list
        
        print(f"Found {len(traits)} traits with status '{status_filter}'")
        print(f"Total traits in database: {sum(len(trait_status[s]) for s in ['Todo', 'Done', 'Skipped'])}")
        print(f"Status counts: Todo={len(trait_status['Todo'])}, "
              f"Done={len(trait_status['Done'])}, "
              f"Skipped={len(trait_status['Skipped'])}")
        
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
        old_status = None
        for status in ['Todo', 'Done', 'Skipped']:
            if trait in trait_status[status]:
                old_status = status
                break
        
        if old_status:
            # Remove from old status
            del trait_status[old_status][trait]
            trait_status['_order'][old_status].remove(trait)
            
            # Add to new status
            trait_status[new_status][trait] = new_status
            trait_status['_order'][new_status].append(trait)
            
            print(f"Updated status for {trait} to {new_status}")
            print(f"New status counts: Todo={len(trait_status['Todo'])}, "
                  f"Done={len(trait_status['Done'])}, "
                  f"Skipped={len(trait_status['Skipped'])}")
            save_trait_status()
    except Exception as e:
        print(f"Error updating trait status: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())

def save_annotation(trait, snp_info, selected_genes):
    """Save annotation"""
    try:
        print(f"\nSaving annotation for trait: {trait}")
        print(f"SNP info: {snp_info}")
        print(f"Number of genes to save: {len(selected_genes)}")
        
        if trait not in annotations:
            annotations[trait] = []
        
        # Get existing gene IDs for this trait
        existing_gene_ids = {ann['gene_id'] for ann in annotations[trait]}
        
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
                annotations[trait].append(annotation)
                existing_gene_ids.add(gene['gene_id'])  # Add to set to prevent future duplicates
                new_annotations.append(annotation)
                print(f"Added annotation for gene: {gene['gene_id']}")
            else:
                print(f"Skipped duplicate gene: {gene['gene_id']}")

        save_annotations()
        print(f"Saved {len(new_annotations)} new annotations for trait {trait}")
    except Exception as e:
        print(f"Error in save_annotation: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())

def get_annotations(trait):
    """Get annotations for a trait"""
    return annotations.get(trait, [])

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
        return "No SNP selected"
    
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
        return "No genes found within 100kb"
    
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
    
    load_trait_status()
    load_annotations()
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
        dcc.Store(id='selected-genes-store'),
        dcc.Store(id='current-snp-store')
    ], fluid=True, style={'backgroundColor': '#f8f9fa', 'padding': '20px'})
])

@app.callback(
    [Output('current-trait-display', 'children'),
     Output('manhattan-plot', 'figure'),
     Output('manhattan-plot', 'selectedData')],
    Input('status-filter', 'value'),
    prevent_initial_call=False
)
def initialize_trait(status_filter):
    """Initialize the first trait when the app loads"""
    print("\nInitializing first trait...")
    traits = get_trait_status(status_filter)
    if not traits:
        print("No traits available for initialization")
        return "No traits available", create_empty_manhattan_plot(), None
    print(f"Initializing with first trait: {traits[0]}")
    return traits[0], create_manhattan_plot(traits[0]), None

@app.callback(
    [Output('current-trait-display', 'children', allow_duplicate=True),
     Output('manhattan-plot', 'figure', allow_duplicate=True),
     Output('manhattan-plot', 'selectedData', allow_duplicate=True)],
    [Input('next-trait-btn', 'n_clicks'),
     Input('prev-trait-btn', 'n_clicks'),
     Input('skip-trait-btn', 'n_clicks'),
     Input('mark-undone-btn', 'n_clicks'),
     Input('done-trait-btn', 'n_clicks')],
    [State('current-trait-display', 'children'),
     State('status-filter', 'value'),
     State('genes-table', 'selected_rows'),
     State('genes-table', 'data'),
     State('current-snp-store', 'data')],
    prevent_initial_call=True
)
def update_trait_and_plot(n_clicks_next, n_clicks_prev, n_clicks_skip, n_clicks_undone, n_clicks_done, 
                         current_trait, status_filter, selected_rows, genes_data, current_snp):
    try:
        print(f"\nUpdating trait and plot...")
        print(f"Status filter: {status_filter}")
        print(f"Current trait: {current_trait}")
        
        ctx = dash.callback_context
        if not ctx.triggered:
            return dash.no_update, dash.no_update, dash.no_update
        
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        print(f"Button triggered: {button_id}")
        
        # Handle Done button with gene selection
        if button_id == 'done-trait-btn':
            if not selected_rows or not genes_data or not current_snp:
                print("Cannot mark as done: No genes selected or no gene data available")
                return dash.no_update, dash.no_update, dash.no_update
            
            # Save selected genes first
            selected_genes = []
            for row_idx in selected_rows:
                if row_idx < len(genes_data):
                    gene = genes_data[row_idx]
                    print(f"Processing gene: {gene['gene_id']}")
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
                print(f"Saving {len(selected_genes)} genes for trait {current_trait}")
                save_annotation(current_trait, current_snp, selected_genes)
                print("Genes saved successfully")
                
                # Now mark as done and move to next trait
                print(f"Marking trait as done: {current_trait}")
                update_trait_status(current_trait, 'Done')
                traits = get_trait_status(status_filter)
                if len(traits) == 0:
                    print("No traits available after marking done")
                    return "No traits available", create_empty_manhattan_plot(), None
                new_trait = traits[0]
                print(f"Moving to next trait after marking done: {new_trait}")
                return new_trait, create_manhattan_plot(new_trait), None
            else:
                print("No genes to save")
                return dash.no_update, dash.no_update, dash.no_update
        
        traits = get_trait_status(status_filter)
        print(f"Found {len(traits)} traits for status {status_filter}")
        
        if len(traits) == 0:
            print("No traits available")
            return "No traits available", create_empty_manhattan_plot(), None
        
        if not current_trait or current_trait not in traits:
            print(f"Current trait {current_trait} not found in traits list, returning first trait")
            return traits[0], create_manhattan_plot(traits[0]), None
        
        current_index = traits.index(current_trait)
        print(f"Current index: {current_index}")
        
        if button_id == 'next-trait-btn':
            next_index = (current_index + 1) % len(traits)
            new_trait = traits[next_index]
            print(f"Moving to next trait: {new_trait}")
        elif button_id == 'prev-trait-btn':
            prev_index = (current_index - 1) % len(traits)
            new_trait = traits[prev_index]
            print(f"Moving to previous trait: {new_trait}")
        elif button_id == 'skip-trait-btn':
            print(f"Skipping trait: {current_trait}")
            update_trait_status(current_trait, 'Skipped')
            traits = get_trait_status(status_filter)
            if len(traits) == 0:
                print("No traits available after skipping")
                return "No traits available", create_empty_manhattan_plot(), None
            new_trait = traits[0]
            print(f"Moving to next trait after skip: {new_trait}")
        elif button_id == 'mark-undone-btn':
            print(f"Marking trait as undone: {current_trait}")
            # Clear annotations for this trait
            if current_trait in annotations:
                print(f"Clearing {len(annotations[current_trait])} annotations for trait {current_trait}")
                annotations[current_trait] = []
                save_annotations()
            update_trait_status(current_trait, 'Todo')
            traits = get_trait_status(status_filter)
            if len(traits) == 0:
                print("No traits available after marking undone")
                return "No traits available", create_empty_manhattan_plot(), None
            new_trait = traits[0]
            print(f"Moving to next trait after marking undone: {new_trait}")
        else:
            print(f"Unknown button: {button_id}")
            return dash.no_update, dash.no_update, dash.no_update
        
        return new_trait, create_manhattan_plot(new_trait), None
    except Exception as e:
        print(f"Error in update_trait_and_plot: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())
        return "Error occurred", create_empty_manhattan_plot(), None

@app.callback(
    [Output('nearby-genes-table', 'children'),
     Output('current-snp-store', 'data')],
    [Input('manhattan-plot', 'selectedData'),
     Input('current-trait-display', 'children')]
)
def update_nearby_genes(selected_data, selected_trait):
    if not selected_trait or selected_trait == "No traits available":
        return "No trait selected", None
    
    if not selected_data or not selected_data['points']:
        trait_data = gwas_df[gwas_df['TRAIT'] == selected_trait]
        if len(trait_data) == 0:
            return "No data available for this trait", None
        
        best_snp = trait_data.loc[trait_data['PVAL'].idxmin()]
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
        'trait': selected_trait
    }
    
    genes_table = create_nearby_genes_table({
        'CHROM': chrom,
        'POS': pos,
        'PVAL': pval
    })
    
    return genes_table, current_snp

@app.callback(
    Output('trait-annotations-table', 'children'),
    [Input('current-trait-display', 'children')]
)
def update_annotations_table(selected_trait):
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

@app.callback(
    Output('trait-annotations-table', 'children', allow_duplicate=True),
    Input('clear-annotations-btn', 'n_clicks'),
    State('current-trait-display', 'children'),
    prevent_initial_call=True
)
def clear_annotations(n_clicks, selected_trait):
    if not n_clicks or not selected_trait:
        raise dash.exceptions.PreventUpdate
    
    if selected_trait in annotations:
        print(f"Clearing all annotations for trait {selected_trait}")
        annotations[selected_trait] = []
        save_annotations()
    
    return update_annotations_table(selected_trait)

@app.callback(
    Output('export-csv-btn', 'n_clicks'),
    Input('export-csv-btn', 'n_clicks'),
    prevent_initial_call=True
)
def handle_export(n_clicks):
    if not n_clicks:
        return dash.no_update
    
    # Convert annotations to DataFrame
    rows = []
    for trait, trait_annotations in annotations.items():
        for ann in trait_annotations:
            row = ann.copy()
            row['trait'] = trait
            rows.append(row)
    
    if not rows:
        return dash.no_update
    
    df = pd.DataFrame(rows)
    # Reorder columns to put trait first
    cols = ['trait'] + [col for col in df.columns if col != 'trait']
    df = df[cols]
    csv_file = 'annotations.csv'
    df.to_csv(csv_file, index=False)
    return dash.no_update

@app.callback(
    Output('trait-annotations-table', 'children', allow_duplicate=True),
    Input('save-genes-btn', 'n_clicks'),
    [State('genes-table', 'selected_rows'),
     State('genes-table', 'data'),
     State('current-snp-store', 'data'),
     State('current-trait-display', 'children')],
    prevent_initial_call=True
)
def save_selected_genes(n_clicks_save, selected_rows, genes_data, current_snp, selected_trait):
    try:
        print("\nAttempting to save selected genes...")
        print(f"Save button clicks: {n_clicks_save}")
        print(f"Current SNP: {current_snp}")
        print(f"Selected trait: {selected_trait}")
        print(f"Selected rows: {selected_rows}")
        print(f"Genes data available: {genes_data is not None}")
        
        if not n_clicks_save or not current_snp or not selected_trait:
            print("Missing required data for saving")
            raise dash.exceptions.PreventUpdate
        
        if not selected_rows or not genes_data:
            print("No genes selected or no gene data available")
            return update_annotations_table(selected_trait)
        
        selected_genes = []
        for row_idx in selected_rows:
            if row_idx < len(genes_data):
                gene = genes_data[row_idx]
                print(f"Processing gene: {gene['gene_id']}")
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
            print(f"Saving {len(selected_genes)} genes for trait {selected_trait}")
            save_annotation(selected_trait, current_snp, selected_genes)
            print("Genes saved successfully")
        else:
            print("No genes to save")
        
        return update_annotations_table(selected_trait)
    except Exception as e:
        print(f"Error in save_selected_genes: {str(e)}")
        print("Full traceback:")
        import traceback
        print(traceback.format_exc())
        return update_annotations_table(selected_trait)

@app.callback(
    Output("toast", "is_open"),
    Output("toast", "children"),
    [Input('done-trait-btn', 'n_clicks')],
    [State('current-trait-display', 'children'),
     State('genes-table', 'selected_rows'),
     State('genes-table', 'data')],
    prevent_initial_call=True
)
def show_done_warning(n_clicks, selected_trait, selected_rows, genes_data):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    if not selected_rows or not genes_data:
        return True, "Please select at least one gene before marking the trait as done."
    
    return False, ""

@app.callback(
    Output('status-counter', 'children'),
    Input('status-filter', 'value')
)
def update_status_counter(status_filter):
    todo_count = len(trait_status['Todo'])
    done_count = len(trait_status['Done'])
    skipped_count = len(trait_status['Skipped'])
    total_count = todo_count + done_count + skipped_count
    
    return f"Total: {total_count} | Todo: {todo_count} | Done: {done_count} | Skipped: {skipped_count}"

if __name__ == '__main__':
    app.run_server(debug=True) 