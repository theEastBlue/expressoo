from dash import Dash, dcc, html, Input, Output, callback
import plotly.express as px
import pandas as pd
import json
import time
from pathlib import Path
import numpy as np
import plotly.graph_objects as go
import dash_auth

# ---------------------------
# Load metadata
# ---------------------------
metadata = pd.read_csv("meta.csv")
group_vars = ["Sex", "long.COVID"]
cell_types = metadata['predicted.celltype.l1'].unique()
module_scores = ["cd8cyto_rna1", "ifnk_rna1", "tsem_rna1", "isr_rna1", "ifns_rna1", "cd8e_rna1", "tcellse_rna1"]

# ---------------------------
# Load UMAP JSON
# ---------------------------
JSON_PATH = Path("umap_plot_cell_type.json")
SAMPLE_N = 50000
MARKER_SIZE_DEFAULT = 2
MARKER_OPACITY_DEFAULT = 0.6

with open(JSON_PATH) as f:
    umap_raw = json.load(f)

rows = []
for tr in umap_raw["data"]:
    name = tr.get("name") or tr.get("legendgroup") or "unknown"
    xs = tr.get("x", [])
    ys = tr.get("y", [])
    texts = tr.get("text", [None]*len(xs))
    for x, y, t in zip(xs, ys, texts):
        rows.append({"x": x, "y": y, "text": t, "celltype": name})

df_umap = pd.DataFrame(rows)
celltype_order = sorted(df_umap['celltype'].unique().tolist()) if not df_umap.empty else []

# ---------------------------
# Load 2nd UMAP (Mapping Score-based)
# ---------------------------
JSON_PATH_2 = Path("umap_plot_mapping_score.json")
with open(JSON_PATH_2) as f:
    umap_data2 = json.load(f)

main_trace = umap_data2["data"][0]
df_mapping = pd.DataFrame({
    "x": main_trace["x"],
    "y": main_trace["y"],
    "color": main_trace["marker"]["color"],
    "text": main_trace["text"]
})

def extract_score(text):
    import re
    match = re.search(r"mapping\.score:\s*([0-9.]+)", text)
    return float(match.group(1)) if match else None

df_mapping["score"] = df_mapping["text"].apply(extract_score)
score_min = max(0, np.floor(df_mapping["score"].min() * 4) / 4)
score_max = min(1.0, np.ceil(df_mapping["score"].max() * 4) / 4)
score_step = 0.05
score_marks = {str(round(s, 2)): f"{s:.2f}" for s in np.arange(score_min, score_max + 0.001, 0.25)}


VALID_USERNAME_PASSWORD_PAIRS = {
    'hello': 'world'
}

# ---------------------------
# Dash App
# ---------------------------
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server
app.title = "EXPresso"
auth = dash_auth.BasicAuth(
    app,
    VALID_USERNAME_PASSWORD_PAIRS
)


# ---------------------------
# Figure builder (WebGL)
# ---------------------------
def make_umap_fig(dataframe, marker_size=MARKER_SIZE_DEFAULT, opacity=MARKER_OPACITY_DEFAULT, showlegend=True):
    fig = px.scatter(
        dataframe,
        x="x",
        y="y",
        color="celltype",
        hover_data=["text"],
        title=f"UMAP Cell Types - {len(dataframe):,} points",
        render_mode="webgl",
        category_orders={"celltype": celltype_order}
    )
    fig.update_traces(marker=dict(size=marker_size, opacity=opacity))
    fig.update_layout(
        hovermode="closest",
        clickmode="event+select",
        legend_title_text="Cell Type",
        showlegend=showlegend,
        margin=dict(l=40, r=10, t=50, b=40)
    )
    return fig

# ---------------------------
# Layout
# ---------------------------
styles = {'pre': {'border': 'thin lightgrey solid', 'overflowX': 'scroll', 'maxHeight': '250px'}}

app.layout = html.Div([
    html.H2("☕️ EXPresso: Expression, X"),
    
    dcc.Tabs([
        # ---------------- Tab 1: Metadata Analysis ----------------
        dcc.Tab(label='Metadata Analysis', children=[
                    
        html.H4("Filter by Cell Type Hierarchy"),
        html.P("Hierarchy Level:"),
        dcc.Dropdown(
            id='hierarchy-level',
            options=[
                {'label': 'L1', 'value': 'predicted.celltype.l1'},
                {'label': 'L2', 'value': 'predicted.celltype.l2'},
                {'label': 'L3', 'value': 'predicted.celltype.l3'}
            ],
            value='predicted.celltype.l1'
        ),
        html.P("Select Cell Types:"),
        dcc.Dropdown(id='hierarchy-filter', multi=True),
        html.H4("Stacked Barplot"),
        html.P("Group by:"),
        dcc.Dropdown(id='stacked-group', options=[{'label': g, 'value': g} for g in group_vars], value='Sex'),
        dcc.Graph(id="stacked-bar"),
        
        # add heatmap

        html.H4("Boxplot"), #normalise this in the future for more meaningful comparisons
        html.P("Group by:"),
        dcc.Dropdown(id='box-group', options=[{'label': g, 'value': g} for g in group_vars], value='Sex'),
        html.P("Cell type:"),
        dcc.Dropdown(id='box-cell', options=[{'label': c, 'value': c} for c in cell_types], value='CD8 T'),
        dcc.Graph(id="box-plot"),

        html.H4("Violin plot (Gene Expression)"),
        html.P("Group by:"),
        dcc.Dropdown(id='violin-group', options=[{'label': g, 'value': g} for g in group_vars], value='Sex'),
        html.P("Module score / Gene:"), 
        dcc.Dropdown(id='violin-score', options=[{'label': s, 'value': s} for s in module_scores], value='ifnk_rna1'),
        html.P("Normalization:"),
        dcc.Dropdown(
            id="violin-norm-method",
            options=[
                {"label": "Log(1+x)", "value": "log1p"},
                {"label": "Raw counts", "value": "raw"}
            ],
            value="log1p"
        ),
        dcc.Graph(id="violin-plot"),

        html.H4("UMAP"),
        html.P("Color by:"),
        dcc.Dropdown(id='umap-color', options=[
            {'label': 'Cell type', 'value': 'predicted.celltype.l1'},
            {'label': 'Long COVID', 'value': 'long.COVID'}
        ], value='predicted.celltype.l1'),
        dcc.Graph(id="umap-plot")
        ]),

        # ---------------- Tab 2: Pool UMAP ----------------
        dcc.Tab(label='Pool UMAP', children=[
            html.H3("UMAP 1: Cell Types"),
            html.Div([
                html.Label("Dataset mode"),
                dcc.RadioItems(
                    id="mode-radio",
                    options=[{"label": f"Sampled ({min(SAMPLE_N, len(df_umap)):,})", "value": "sample"},
                            {"label": "Full", "value": "full"}],
                    value="sample",
                    inline=True
                ),
                html.Label("Marker size"),
                dcc.Slider(id='marker-size', min=1, max=6, step=1, value=MARKER_SIZE_DEFAULT,
                        marks={i: str(i) for i in range(1, 7)}),
                html.Label("Opacity"),
                dcc.Slider(id='marker-opacity', min=0.2, max=1.0, step=0.1, value=MARKER_OPACITY_DEFAULT,
                        marks={0.2: '0.2', 0.5: '0.5', 0.8: '0.8', 1.0: '1.0'}),
                html.Label("Legend"),
                dcc.Checklist(id='toggle-legend', options=[{"label": "Show legend", "value": "on"}], value=["on"]),
                dcc.Graph(id='umap-graph'),
                html.Div(id='perf-metrics')
            ]),

            html.Hr(),

            html.H3("UMAP 2: Mapping Score"),
            html.Div([
                html.Label("Marker Size"),
                dcc.Slider(id='marker-size2', min=1, max=6, step=1, value=3,
                        marks={i: str(i) for i in range(1, 7)}),
                html.Label("Opacity"),
                dcc.Slider(id='marker-opacity2', min=0.2, max=1.0, step=0.1, value=0.8,
                        marks={0.2: '0.2', 0.5: '0.5', 0.8: '0.8', 1.0: '1.0'}),
                html.Label("Mapping Score Range"),
                dcc.RangeSlider(
                    id='score-range',
                    min=score_min,
                    max=score_max,
                    step=score_step,
                    value=[score_min, score_max],
                    marks=score_marks,
                    tooltip={"placement": "bottom", "always_visible": False}
                ),
                dcc.Graph(id='umap-graph2'),
                html.Div(id='perf-metrics2')
            ])
        ])

    ])
])

# ---------------------------
# Callbacks: Metadata Analysis
# ---------------------------
@app.callback(
    Output('hierarchy-filter', 'options'),
    Output('hierarchy-filter', 'value'),
    Input('hierarchy-level', 'value')
)
def update_hierarchy_options(level):
    unique_types = metadata[level].dropna().unique().tolist()
    unique_types.sort()
    # Select all by default
    return [{'label': t, 'value': t} for t in unique_types], unique_types

@app.callback(
    Output("stacked-bar", "figure"),
    Input("stacked-group", "value"),
    Input("hierarchy-level", "value"),
    Input("hierarchy-filter", "value")
)
def update_stacked(group, level, selected_types):
    df_filtered = metadata[metadata[level].isin(selected_types)]
    counts = df_filtered.groupby([group, level]).size().reset_index(name="count")
    fig = px.bar(counts, x=group, y="count", color=level,
                 title=f"Cell type proportions grouped by {group}")
    return fig

@app.callback(
    Output("box-plot", "figure"),
    Input("box-group", "value"),
    Input("hierarchy-level", "value"),
    Input("hierarchy-filter", "value")
)
def update_box(group, level, selected_types):
    df_filtered = metadata[metadata[level].isin(selected_types)]
    fig = px.box(df_filtered, x=group, y="nCount_RNA", points="all", color=level,
                 title=f"nCount_RNA by {group} for selected cell types")
    return fig


@app.callback(
    Output("violin-plot", "figure"),
    Input("violin-group", "value"),
    Input("violin-score", "value"),
    Input("violin-norm-method", "value"),  # Keep normalization input
    Input("hierarchy-level", "value"),
    Input("hierarchy-filter", "value")
)
def update_violin(group, score, norm_method, level, selected_types):
    # Filter metadata for selected hierarchy
    df_filtered = metadata[metadata[level].isin(selected_types)]
    
    # Compute normalized gene expression for the selected score
    df_expr = calculate_gene_expression(df_filtered, genes=[score], method=norm_method)
    expr_col = score + "_expr"

    # Make violin plot
    fig = px.violin(df_expr, x=group, y=expr_col, color=level, box=True,
                    title=f"{norm_method} normalized {score} expression across selected cell types by {group}")
    return fig


@app.callback(
    Output("umap-plot", "figure"),
    Input("umap-color", "value"),
    Input("hierarchy-level", "value"),
    Input("hierarchy-filter", "value")
)
def update_umap(color, level, selected_types):
    if 'UMAP_1' not in metadata.columns or 'UMAP_2' not in metadata.columns:
        return px.scatter(title="UMAP coordinates not found in metadata")
    df_filtered = metadata[metadata[level].isin(selected_types)]
    fig = px.scatter(df_filtered, x='UMAP_1', y='UMAP_2', color=color,
                     facet_col='Sex' if 'Sex' in df_filtered.columns else None,
                     title=f"UMAP split by Sex, colored by {color}")
    return fig



# ---------------------------
# Callbacks: Pool UMAP
# ---------------------------
@app.callback(
    Output('umap-graph', 'figure'),
    Output('perf-metrics', 'children'),
    Input('mode-radio', 'value'),
    Input('marker-size', 'value'),
    Input('marker-opacity', 'value'),
    Input('toggle-legend', 'value')
)
def update_umap_viewer(mode, marker_size, opacity, legend_value):
    if mode == 'sample' and len(df_umap) > SAMPLE_N:
        dfx = df_umap.sample(n=SAMPLE_N, random_state=42)
        mode_label = f"Sampled ({SAMPLE_N:,})"
    else:
        dfx = df_umap
        mode_label = f"Full ({len(df_umap):,})"
    
    # Ensure that legend_value is always a list
    if not isinstance(legend_value, list):
        legend_value = []

    # Show legend if selected
    showlegend = ("on" in legend_value)  # Check if 'on' is in the list

    # Create the figure with updated marker size and opacity
    fig = make_umap_fig(dfx, marker_size=marker_size, opacity=opacity, showlegend=showlegend)
    
    metrics = [
        f"Points: {len(df_umap):,}",
        f"Cell types: {len(celltype_order):,}",
        f"Mode: {mode_label}"
    ]
    metrics_div = html.Ul([html.Li(m) for m in metrics])
    return fig, metrics_div

@app.callback(Output('hover-data', 'children'), Input('umap-graph', 'hoverData'))
def display_hover_data(hoverData):
    return json.dumps(hoverData, indent=2)

@app.callback(Output('click-data', 'children'), Input('umap-graph', 'clickData'))
def display_click_data(clickData):
    return json.dumps(clickData, indent=2)

@app.callback(Output('selected-data', 'children'), Input('umap-graph', 'selectedData'))
def display_selected_data(selectedData):
    return json.dumps(selectedData, indent=2)

@app.callback(Output('relayout-data', 'children'), Input('umap-graph', 'relayoutData'))
def display_relayout_data(relayoutData):
    return json.dumps(relayoutData, indent=2)


@app.callback(
    Output('umap-graph2', 'figure'),
    Output('perf-metrics2', 'children'),
    Input('marker-size2', 'value'),
    Input('marker-opacity2', 'value'),
    Input('score-range', 'value')
)
def update_umap2(marker_size, opacity, score_range):
    _t0 = time.time()

    # Ensure score_range is a valid tuple (low, high)
    if not isinstance(score_range, list) or len(score_range) != 2:
        score_range = [score_min, score_max]  # Default range if invalid

    low, high = score_range
    filtered_df = df_mapping[(df_mapping["score"] >= low) & (df_mapping["score"] <= high)]

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=filtered_df["x"],
        y=filtered_df["y"],
        mode='markers',
        marker=dict(
            color=filtered_df["color"],
            size=marker_size,  # Set marker size from slider
            opacity=opacity,  # Set opacity from slider
            line=dict(width=0)
        ),
        text=filtered_df["text"],
        hoverinfo="text",
        showlegend=False
    ))
    fig.update_layout(
        title=f"UMAP by Mapping Score (Filtered: {len(filtered_df):,} points)",
        xaxis_title="UMAP_1",
        yaxis_title="UMAP_2",
        hovermode="closest",
        margin=dict(l=40, r=10, t=50, b=40)
    )

    # Performance metrics
    return fig, f"Filtered points: {len(filtered_df):,} - Build time: {time.time() - _t0:.3f}s"


# ---------------------------
# Gene expression calculation
# ---------------------------
def calculate_gene_expression(df, genes=None, method="log1p"):
    if genes is None:
        genes = module_scores

    expr_df = df.copy()
    for g in genes:
        if g in expr_df.columns:
            if method == "log1p":
                # Clip values below -1 to avoid invalid log1p computation, or use np.maximum
                cleaned_vals = np.maximum(expr_df[g], -0.9999) 
                expr_df[g + "_expr"] = np.log1p(cleaned_vals)
            else:
                expr_df[g + "_expr"] = expr_df[g]
    return expr_df

# ---------------------------
if __name__ == "__main__":
    app.run(debug=True)
