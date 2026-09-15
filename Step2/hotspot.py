#ref
##https://github.com/JingtaoLab/STO-analysis/blob/STO-analysis/02-regulon/step2_hotspot.py
#https://github.com/JingtaoLab/STO-analysis/blob/STO-analysis/02-regulon/step3_module_plot.py

import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import os
import pickle

counts_file = 'TP*_counts.txt'
pos_file = 'TP*_pos.txt'

pos = pd.read_csv(pos_file, index_col=0, sep=" ", header=None)
counts = pd.read_csv(counts_file, index_col=0, sep=" ")
counts = counts.loc[:, pos.index]
barcodes = pos.index.values
num_umi = counts.sum(axis=0)

gene_counts = (counts > 0).sum(axis=1)
valid_genes = gene_counts >= 50
counts = counts.loc[valid_genes]

hs = hotspot.Hotspot(counts, model='danb', latent=pos, umi_counts=num_umi)
hs.create_knn_graph(weighted_graph=True, n_neighbors=30)
hs_results = hs.compute_autocorrelations(jobs=4)
hs.results.to_csv(path_or_buf=os.path.join(output_data_path, f"{tp_prefix}_gene_results.csv"))# hs_genes = hs_results.query("FDR < 0.05").sort_values("Z", ascending=False).head(1000).index #https://github.com/search?q=hs_genes+%3D+hs_results.head%28+language%3APython+&type=code
hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False).head(1000).index
lcz = hs.compute_local_correlations(hs_genes, jobs=10)
modules = hs.create_modules(min_gene_threshold=50, core_only=False, fdr_threshold=0.05)  
modules.to_csv(path_or_buf=os.path.join(output_data_path, f"{tp_prefix}_cluster.csv"))
fig = plt.figure(figsize=(6, 6))
fig = hs.plot_local_correlations()
plt.savefig(os.path.join(output_data_path, f"{tp_prefix}_cor.png"), dpi=100, bbox_inches='tight')
plt.close(fig)
subdir = ['module_scores', 'module_genes', 'module_3d_json', 'module_2d_plot']
for i in subdir:
    os.makedirs(os.path.join(output_data_path, i), exist_ok=True)

module_scores = hs.calculate_module_scores()
max_module = hs.modules.max() if (hs.modules is not None and not hs.modules.empty) else 0
for module in range(1, max_module + 1):
    scores = hs.module_scores[module]
    plotly_df = pd.DataFrame({'x': hs.latent.iloc[:, 0],
            'y': hs.latent.iloc[:, 1],
            'z': hs.latent.iloc[:, 2],
            'scores': scores,
            'cell_id': adata.obs_names
        })
    score_path = os.path.join(output_data_path, "module_scores", f"module_{module}_scores.csv")
    plotly_df.to_csv(score_path)
    plotly_df['x'] = plotly_df['x'].round(1)
    plotly_df['y'] = plotly_df['y'].round(1)
    plotly_df['z'] = plotly_df['z'].round(1)
    plotly_df['scores'] = plotly_df['scores'].round(3)

    json_path = os.path.join(output_data_path, "module_3d_json", f"module_{module}_3d.json")
    plotly_df.to_json(json_path, orient='split', index=False, force_ascii=False)

        # 2D Plot
    if 'ggplot' in globals():
        plotly_df['slice_num'] = plotly_df['z'].rank(method='dense').astype(int) - 1
        p = ggplot(aes(x="x", y="y"), plotly_df) + geom_point(aes(color='scores'), size=1) + facet_wrap('slice_num', ncol=5)+ coord_fixed(ratio=1) + theme_void() + theme(strip_text=element_text(size=8))

        ggsave(p, os.path.join(output_data_path, f"module_2d_plot/module_{module}_2d.png"), width=15, height=15, verbose=False)
        ggsave(p, os.path.join(output_data_path, f"module_2d_plot/module_{module}_2d.pdf"), width=15, height=15, verbose=False)
        results = hs.results.join(hs.modules)
        results = results.loc[results.Module == module]
        results = results.sort_values('Z', ascending=False)
        results.to_csv(os.path.join(output_data_path, "module_genes", f"module_{module}_genes.csv"))