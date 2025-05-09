import numpy as np
import plotly.graph_objects as go
from scipy.spatial import Delaunay

def plot_3d_filled_volume(
    df,
    id_col='id',
    x_col='left_open',
    y_col='right_open',
    z_col='dF_total',
    x_label=None,
    y_label=None,
    z_label=None,
    width=1000,
    height=700,
    colorscale='Plasma',
    grid_spacing=1,
    grid_color='black',
    grid_width=2,
    contour_z=False,
    lighting=dict(ambient=0.6, diffuse=0.5, roughness=1.0, specular=0.1),
    lightposition=dict(x=100, y=200, z=0),
    title="Free‐Energy Volume Plot"
):
    """
    Generate an interactive 3D filled volume plot with:
      - a closed mesh volume (Mesh3d) under the surface
      - flat shading + custom lighting
      - black grid‐lines on the floor and vertical walls
      - colored z‐contours projected to the floor
      - a slider to toggle different groups (by id)
    
    Parameters
    ----------
    df : pd.DataFrame
      Data containing at least these columns: [id_col, x_col, y_col, z_col].
    id_col : str
      Column name for grouping and slider.
    x_col : str
      Column name for the x‐axis.
    y_col : str
      Column name for the y‐axis.
    z_col : str
      Column name for the z‐axis (value).
    width, height : int
      Figure size in pixels.
    colorscale : str
      Plotly colorscale for the surface.
    grid_spacing : float
      Interval for black grid lines on x and y.
    grid_color : str
      Color of the grid lines.
    grid_width : float
      Width of the grid lines.
    contour_z : bool
      Whether to show colored z‐contours.
    lighting : dict
      Lighting settings for Mesh3d facets.
    lightposition : dict
      Light source position.
    title : str
      Figure title.

    Returns
    -------
    fig : go.Figure
      A Plotly Figure object ready to show().
    """
    ids = np.sort(df[id_col].unique())
    data = []

     # Axis titles
    x_title = x_label if x_label is not None else x_col
    y_title = y_label if y_label is not None else y_col
    z_title = z_label if z_label is not None else z_col

    for idx, id_val in enumerate(ids):
        df_id = df[df[id_col] == id_val]
        pivot = df_id.pivot(index=x_col, columns=y_col, values=z_col)
        X, Y = np.meshgrid(pivot.columns.values, pivot.index.values)
        Z = pivot.values 

        # triangulate for Mesh3d
        pts = np.column_stack([X.flatten(), Y.flatten()])
        tri = Delaunay(pts)
        tris = tri.simplices

        # build top + bottom vertices
        top_vs = np.column_stack([pts,     Z.flatten()])
        bot_vs = np.column_stack([pts, np.zeros_like(Z.flatten())])
        verts = np.vstack([top_vs, bot_vs])
        n_top = top_vs.shape[0]

        # faces: top surface
        i_top, j_top, k_top = tris.T
        # bottom surface (reversed)
        i_bot = tris[:, 0] + n_top
        j_bot = tris[:, 2] + n_top
        k_bot = tris[:, 1] + n_top

        # side walls along convex hull
        hull = tri.convex_hull
        i_side, j_side, k_side = [], [], []
        for a, b in hull:
            # two triangles per edge
            i_side += [a,     a+n_top]
            j_side += [b,     b]
            k_side += [a+n_top, b+n_top]

            i_side += [a+n_top, b+n_top]
            j_side += [b,      a+n_top]
            k_side += [b+n_top, a+n_top]

        # combine all faces
        i = np.hstack([i_top,     i_bot,     i_side])
        j = np.hstack([j_top,     j_bot,     j_side])
        k = np.hstack([k_top,     k_bot,     k_side])

        # the closed volume mesh
        mesh = go.Mesh3d(
            x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
            i=i, j=j, k=k,
            intensity=verts[:, 2],
            colorscale=colorscale,
            flatshading=True,
            lighting=lighting,
            lightposition=lightposition,
            visible=(idx == 0),
            name=f"{id_col}={id_val}"
        )
        data.append(mesh)

        # overlay for grid + contours
        data.append(
            go.Surface(
                x=X, y=Y, z=Z,
                showscale=False,
                opacity=0,
                contours=dict(
                    x=dict(show=True, start=X.min(), end=X.max(), size=grid_spacing,
                           color=grid_color, width=grid_width),
                    y=dict(show=True, start=Y.min(), end=Y.max(), size=grid_spacing,
                           color=grid_color, width=grid_width),
                    z=dict(show=contour_z, usecolormap=True,
                           highlightcolor="white", project_z=True)
                ),
                visible=(idx == 0),
                name=f"{id_col}={id_val}-grid"
            )
        )

    # slider toggles pairs of traces
    steps = []
    for i, id_val in enumerate(ids):
        vis = [False] * len(data)
        vis[2*i] = True
        vis[2*i+1] = True
        steps.append(dict(method="update", args=[{"visible": vis}], label=str(id_val)))

    sliders = [dict(active=0, pad={"t": 50}, steps=steps)]
    layout = go.Layout(
        title=title,
        sliders=sliders,
        scene=dict(
            xaxis_title=x_title,
            yaxis_title=y_title,
            zaxis_title=z_title
        ),
        width=width, height=height,
        margin=dict(l=0, r=0, t=50, b=0),
        showlegend=False
    )

    fig = go.Figure(data=data, layout=layout)
    return fig

