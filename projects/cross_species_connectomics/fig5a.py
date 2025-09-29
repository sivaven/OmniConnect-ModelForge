################################################################################
# File: fig5a.py
# Project: Cross-Species Connectomics
# Author: Siva Venkadesh 
# Date: 2025-09-27
# Description: Generates Figure 5a. Produces directed connectome delta plots to illustrate
#              species-specific efficiency shifts at the network level.
################################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

_dir = "data/"
# Load final dataset
df_final = pd.read_csv(_dir+"full_path_efficiency_table.csv")

short_map = {
    "FRO_Prefrontal": "Pref",
    "FRO_Precentral": "PreC",
    "FRO_Premotor": "PreM",
    "TEM_Hippocampus": "Hipp",
    "TEM_Amygdala": "Amyg",
    "THL_Thalamus": "Thal",
    "THL_Hypothalamus": "Hypo",
    "BG_CaudoPutamen": "CPu",
    "BG_Accumbens": "NAc",
    "BG_Pallidum": "Pall",
    "CIN_Anterior": "CingA",
    "CIN_Posterior": "CingP",
    "INS_Anterior": "InsA",
    "INS_Posterior": "InsP",
    "OLF_Anterior": "OlfA",
    "OLF_Piriform": "Pir",
    "TEM_Medial": "TemM",
    "TEM_Inferior": "TemI",
    "TEM_Superior": "TemS",
    "OCC_Medial": "OccM",
    "OCC_Lateral": "OccL"
}

def shorten(name):
    """Auto-added docstring. See manuscript methods for details."""
    return short_map.get(name, name) 

# --- pick one focal region+role per row ---
row_specs = [
    ("INS_Anterior (Source)", [("human","mouse"), ("human","marmoset"), ("human","rhesus")]),
    ("TEM_Superior (Source)", [("human","mouse"), ("human","marmoset"), ("human","rhesus")]),
    ("INS_Anterior (Target)", [("human","mouse"), ("human","marmoset"), ("human","rhesus")]),
    ("TEM_Inferior (Source)", [("rhesus","mouse"), ("human","mouse"), ("human","rhesus")]),
    ("OLF_Anterior (Source)", [("marmoset","mouse"), ("rhesus","marmoset"), ("human","marmoset")]),
]


# Build per-row dataframes from existing regions_of_interest
roi_map = {
    "INS_Anterior (Target)": df_final[df_final['source__target'].str.endswith("INS_Anterior")],
    "INS_Anterior (Source)": df_final[df_final['source__target'].str.startswith("INS_Anterior")],
    "TEM_Superior (Source)": df_final[df_final['source__target'].str.startswith("TEM_Superior")],
    "TEM_Inferior (Source)": df_final[df_final['source__target'].str.startswith("TEM_Inferior")],
    "OLF_Anterior (Source)": df_final[df_final['source__target'].str.startswith("OLF_Anterior")],
    "OLF_Piriform (Target)": df_final[df_final['source__target'].str.endswith("OLF_Piriform")]
}


def adjust_circular_layout(layout_dict, focal_region, offset=5):
    """Auto-added docstring. See manuscript methods for details."""
    adjusted_layout = {}
    for node, (x, y) in layout_dict.items():
        if node == focal_region:
            adjusted_layout[node] = (x, y)
        else:
            norm = (x**2 + y**2)**0.5
            scale = 1 + offset if norm != 0 else 1
            adjusted_layout[node] = (x * scale, y * scale)
    return adjusted_layout

def draw_delta_connectome_circular_adjusted(df_region, sp1, sp2, ax, region_name, circular_pos, top_n=10):
    """Auto-added docstring. See manuscript methods for details."""
    G = nx.DiGraph()
    df_region = df_region.copy()  
    df_region.loc[:, "Source"] = df_region["source__target"].apply(lambda x: x.split("__")[0])
    df_region.loc[:, "Target"] = df_region["source__target"].apply(lambda x: x.split("__")[1])
    
    is_target_mode = "Target" in region_name
    focal_region = region_name.split(" ")[0]
    node_deltas = {}
    edge_deltas = []

    for _, row in df_region.iterrows():
        src, tgt = row["Source"], row["Target"]
        val1, val2 = row[sp1], row[sp2]
        if val2 == 0:
            continue
        delta = np.log2(val1 / val2)
        partner = src if is_target_mode else tgt
        if partner not in node_deltas:
            node_deltas[partner] = []
        node_deltas[partner].append(delta)
        edge = (partner, focal_region) if is_target_mode else (focal_region, partner)
        edge_deltas.append((edge[0], edge[1], delta))

    avg_deltas_all = {k: np.mean(v) for k, v in node_deltas.items()}
    top_nodes = sorted(avg_deltas_all.items(), key=lambda x: abs(x[1]), reverse=True)[:top_n]
    avg_deltas = dict(top_nodes)

    for node, delta in avg_deltas.items():
        G.add_node(node, size= 200+ 200 * abs(delta), delta=delta)
    G.add_node(focal_region, size=200, delta=0)

    for src, tgt, delta in edge_deltas:
        if (src == focal_region and tgt in avg_deltas) or (tgt == focal_region and src in avg_deltas):
            G.add_edge(src, tgt, weight=abs(delta), color='darkcyan' if delta > 0 else 'indianred')

    adjusted_pos = adjust_circular_layout(circular_pos, focal_region)
    node_sizes = [G.nodes[n]["size"] for n in G.nodes]
    font_sizes = [8 + 6 * abs(G.nodes[n]["delta"]) if G.nodes[n]["delta"] != 0 else 10 for n in G.nodes]
    node_colors = []
    edgecolors = []
    for n in G.nodes:
        delta = G.nodes[n]["delta"]
        if n == focal_region:
            node_colors.append('white')
            edgecolors.append('black')
        else:
            node_colors.append('darkcyan' if delta > 0 else 'indianred')
            edgecolors.append('none')

    nx.draw_networkx_nodes(G, adjusted_pos, nodelist=G.nodes, node_size=node_sizes,
                           node_color=node_colors, edgecolors=edgecolors,
                           linewidths=1.5, ax=ax)
    for node in G.nodes:
        x, y = adjusted_pos[node]
        #ax.text(x, y, node, fontsize=font_sizes[list(G.nodes).index(node)], ha='center', va='center')
        ax.text(x, y, shorten(node), fontsize=16, ha='center', va='center')

    edge_colors = [G[u][v]['color'] for u, v in G.edges()]
    edge_widths = [1 + 1 * G[u][v]['weight'] for u, v in G.edges()]
    nx.draw_networkx_edges(G, adjusted_pos, width=edge_widths, edge_color=edge_colors, ax=ax,
                           arrows=True, arrowstyle='-|>', arrowsize=20, connectionstyle='arc3,rad=0.1')
    ax.set_title(f"{region_name}: {sp1} vs {sp2}", fontsize=10)
    ax.axis('off')
    
#    
#    
# --- plot 5 rows x 3 columns ---
nrows, ncols = 5, 3
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14, 24.5))



def main():
    all_nodes = set()
    for region_name, df_region in roi_map.items():
        df_region = df_region.copy()  
        df_region.loc[:, "Source"] = df_region["source__target"].apply(lambda x: x.split("__")[0])
        df_region.loc[:, "Target"] = df_region["source__target"].apply(lambda x: x.split("__")[1])
        is_target_mode = "Target" in region_name
        focal_region = region_name.split(" ")[0]
        for _, row in df_region.iterrows():
            src, tgt = row["Source"], row["Target"]
            partner = src if is_target_mode else tgt
            all_nodes.add(partner)
            all_nodes.add(focal_region)
    dummy_G = nx.DiGraph()
    dummy_G.add_nodes_from(all_nodes)
    circular_pos = nx.circular_layout(dummy_G)

    fig.subplots_adjust(hspace=0.8, wspace=0.4)
    for r, (region_name, comps) in enumerate(row_specs):
        df_region = roi_map[region_name].copy()
        for c in range(ncols):
            ax = axes[r, c]
            if c < len(comps):
                sp1, sp2 = comps[c]
                draw_delta_connectome_circular_adjusted(df_region.copy(), sp1, sp2, ax, region_name, circular_pos, top_n=5)
            else:
                ax.axis('off')
    plt.show()


if __name__ == '__main__':
    main()
