import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from py_analysis.config.gen_var import RESULTS_DIR

KRESCFACTOR_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
ids = ["601", "601RTA", "601MF", "601L", "5S"]
HARD_CONS = False

# Load data for all KRESCFACTOR values
data_dict = {}
for kresc in KRESCFACTOR_values:
    if HARD_CONS:
        filename = RESULTS_DIR / f"nbfiles/nucbreathfe/breathstatefe_K{kresc}_hc.csv"
    else:
        filename = RESULTS_DIR / f"nbfiles/nucbreathfe/breathstatefe_K{kresc}_sc.csv"
    df = pd.read_csv(filename)
    data_dict[kresc] = df

# Function to compute df_exp_norm for a given E_ads
def compute_df_exp_norm(df, E_ads):
    df["Adsorp_F"] = -E_ads * (14 - (df["left_open"] + df["right_open"]))
    df["dF_total"] = df["dF"] + df["Adsorp_F"]
    df["df_exp"] = df["dF_total"].apply(lambda x: np.exp(-x))
    df["df_exp_norm"] = df.groupby("id")["df_exp"].transform(lambda x: x / x.sum())
    return df


# Add heatmaps for initial E_ads with dynamic zmin and zmax
for i, kresc in enumerate(KRESCFACTOR_values):

    fig = make_subplots(
    rows=1, cols=len(ids),
    subplot_titles=[f"ID={id_val}" for id_val in ids],
    shared_yaxes=True,
    horizontal_spacing=0.05
    )


    # Initial E_ads value
    initial_E_ads = 15



    df = compute_df_exp_norm(data_dict[kresc].copy(), initial_E_ads)
    for j, id_val in enumerate(ids):
        df_id = df[df['id'] == id_val]
        data_vals = df_id['df_exp_norm']
        # Apply log_scale logic: zmin is min non-zero, zmax is max
        nz = data_vals[data_vals > 0]
        vmin_ = nz.min() if not nz.empty else 1e-10  # Fallback for all-zero case
        vmax_ = data_vals.max()

        pivot_df = df_id.pivot(index='right_open', columns='left_open', values='df_exp_norm')
        heatmap = go.Heatmap(
            z=pivot_df.values,
            x=pivot_df.columns,
            y=pivot_df.index,
            colorscale='plasma_r',
            zmin=vmin_,
            zmax=vmax_,
            showscale=(j == 0)  # Show colorbar only for the first heatmap
        )
        fig.add_trace(heatmap, row=1, col=j+1)

# Add a slider for each row (KRESCFACTOR) with dynamic zmin and zmax updates
sliders = []
for i, kresc in enumerate(KRESCFACTOR_values):
    slider_steps = []
    for e in range(1, 31):
        # Compute updated data for this KRESCFACTOR
        df = compute_df_exp_norm(data_dict[kresc].copy(), e)
        z_values = []
        zmin_values = []
        zmax_values = []
        for id_val in ids:
            df_id = df[df['id'] == id_val]
            data_vals = df_id['df_exp_norm']
            # Dynamic zmin and zmax for each heatmap
            nz = data_vals[data_vals > 0]
            vmin_ = nz.min() if not nz.empty else 1e-10
            vmax_ = data_vals.max()
            pivot_df = df_id.pivot(index='right_open', columns='left_open', values='df_exp_norm')
            z_values.append(pivot_df.values)
            zmin_values.append(vmin_)
            zmax_values.append(vmax_)
        # Define which traces to update (all columns in this row)
        trace_indices = [i * len(ids) + j for j in range(len(ids))]
        step = {
            "method": "update",
            "label": str(e),
            "args": [
                {
                    "z": [z_values[j] for j in range(len(ids))],
                    "zmin": [zmin_values[j] for j in range(len(ids))],
                    "zmax": [zmax_values[j] for j in range(len(ids))]
                },
                {"frame": {"redraw": True}},
                trace_indices  # Update only this row's traces
            ]
        }
        slider_steps.append(step)
    
    sliders.append({
        "steps": slider_steps,
        "currentvalue": {"prefix": f"E_ads (K={kresc}): "},
        "pad": {"t": 20},
        "len": 0.5,
        "x": 0.25,
        "y": -0.1 - i * 0.1,  # Position below each row
        "name": f"slider_{kresc}"
    })

fig.update_layout(sliders=sliders)
fig.show()