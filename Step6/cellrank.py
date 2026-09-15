# ref: https://github.com/theislab/moscot-framework_reproducibility/blob/main/notebooks/space/spatiotemporal/3_Full_embryo_CellRank_analysis/ZP_2023-04-20_spatiotemporal_fullembryo-cellrank.ipynb
#https://github.com/theislab/moscot-framework_reproducibility/blob/main/notebooks/space/spatiotemporal/4_Brain_analysis/1_CellRank_analysis/ZP_2023-04-20_spatiotemporal_brain-cellrank.ipynb
import os
import re
import glob

import numpy as np
import pandas as pd
import scipy.sparse as sp

import scanpy as sc
import anndata as ad
import cellrank as cr
from cellrank.estimators import GPCCA

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_OUTPUT_DIR = "../CellRank"

RUN_DATASETS = []
USE_LAYER = "raw_counts"
N_PCS = 50
N_NEIGHBORS = 30
N_TOP_GENES = 2000
N_DIFFMAP = 15
N_STATES = 6
PSEUDOTIME_WEIGHT = 0.8
TOP_DRIVERS = 1000
MAX_CELLS_PER_TP = 80000
MAX_CELLS_TOTAL = 150000
SEED = 42

np.random.seed(SEED)


def tp_number(path):
    return int(re.search(r"TP(\d+)", os.path.basename(path)).group(1))


def subsample(adata, max_cells):
    if adata.n_obs <= max_cells:
        return adata
    per_tp = max_cells // adata.obs["Sequence_TP"].nunique()
    picks = []
    for tp in adata.obs["Sequence_TP"].unique():
        pos = np.where(adata.obs["Sequence_TP"].values == tp)[0]
        picks.append(np.random.choice(pos, min(len(pos), per_tp), replace=False))
    return adata[np.sort(np.concatenate(picks))].copy()


def load_run(files, subset_key=None, subset_val=None):
    parts = []
    for path in files:
        lazy = sc.read_h5ad(path, backed="r")
        if subset_key is None:
            part = lazy.to_memory()
        else:
            mask = (lazy.obs[subset_key].astype(str) == str(subset_val)).values
            part = lazy[mask].to_memory()
        lazy.file.close()
        tp = os.path.basename(path).split("_")[0]
        part.obs_names = [f"{tp}_{name}" for name in part.obs_names]
        part.obs["Sequence_TP"] = tp
        part.obs["TP_Numeric"] = tp_number(path)
        parts.append(subsample(part, MAX_CELLS_PER_TP))
    merged = ad.concat(parts, join="outer", fill_value=0)
    return subsample(merged, MAX_CELLS_TOTAL)


def preprocess(adata):
    adata.X = (adata.layers[USE_LAYER] if USE_LAYER in adata.layers else adata.X).copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=N_TOP_GENES, subset=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.diffmap(adata, n_comps=N_DIFFMAP)
    early = np.where(adata.obs["TP_Numeric"].values == adata.obs["TP_Numeric"].min())[0]
    adata.uns["iroot"] = int(early[np.argmin(adata.obsm["X_diffmap"][early, 0])])
    sc.tl.dpt(adata)
    return adata


def build_kernel(adata):
    pseudotime = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime").compute_transition_matrix()
    connectivity = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    return PSEUDOTIME_WEIGHT * pseudotime + (1 - PSEUDOTIME_WEIGHT) * connectivity


def anno_keys(adata):
    return [k for k in POTENTIAL_KEYS if k in adata.obs.columns and 2 <= adata.obs[k].nunique() <= 50]


def save_figure(fig_dir, name):
    plt.savefig(os.path.join(fig_dir, f"{name}.pdf"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(fig_dir, f"{name}.png"), dpi=300, bbox_inches="tight")
    plt.close("all")


def export_drivers(drivers, path, dataset, run_name):
    rows = []
    for lineage in [c[:-5] for c in drivers.columns if c.endswith("_corr")]:
        top = drivers[f"{lineage}_corr"].dropna().sort_values(ascending=False).head(TOP_DRIVERS)
        for gene, score in top.items():
            row = {"Gene": gene, "Score": round(float(score), 4), "Lineage": lineage,
                   "Dataset": dataset, "TP": run_name}
            for stat in ["pval", "qval"]:
                if f"{lineage}_{stat}" in drivers.columns:
                    row[stat.capitalize()] = float(drivers.loc[gene, f"{lineage}_{stat}"])
            rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)


def run_analysis(dataset, files, run_name, subset_key=None, subset_val=None):
    out_dir = os.path.join(BASE_OUTPUT_DIR, dataset, run_name)
    csv_dir, fig_dir, h5ad_dir = (os.path.join(out_dir, d) for d in ["CSV", "Figures", "H5AD"])
    for d in [csv_dir, fig_dir, h5ad_dir]:
        os.makedirs(d, exist_ok=True)
    print(f"[{dataset}] {run_name} ({len(files)} timepoints)")

    adata = preprocess(load_run(files, subset_key, subset_val))
    kernel = build_kernel(adata)

    for key in anno_keys(adata):
        adata.obs[key] = adata.obs[key].astype("category")
        g = GPCCA(kernel)
        g.compute_schur()
        g.plot_spectrum(real_only=True)
        save_figure(fig_dir, f"{run_name}_{key}_spectrum")

        g.compute_macrostates(n_states=N_STATES, cluster_key=key)
        g.plot_macrostates(which="all", legend_loc="right", dpi=300, s=50, show=False)
        save_figure(fig_dir, f"{run_name}_{key}_macrostates")

        g.predict_terminal_states()
        g.plot_macrostates(which="terminal", legend_loc="right", dpi=300, s=50, show=False)
        save_figure(fig_dir, f"{run_name}_{key}_terminal_states")

        g.compute_fate_probabilities()
        g.plot_fate_probabilities()
        save_figure(fig_dir, f"{run_name}_{key}_fate_probabilities")

        cr.pl.circular_projection(adata, keys=key, legend_loc="right", figsize=(20, 20), dpi=200, s=5)
        save_figure(fig_dir, f"{run_name}_{key}_circular")
        export_drivers(g.compute_lineage_drivers(cluster_key=key),
                       os.path.join(csv_dir, f"{run_name}_{key}_drivers.csv"), dataset, run_name)

def common_values(files, key):
    values = None
    for path in files:
        lazy = sc.read_h5ad(path, backed="r")
        found = set(lazy.obs[key].astype(str))
        lazy.file.close()
        values = found if values is None else values & found
    return sorted(v for v in values if v.lower() not in ["nan", "unknown", ""])


def run_subsets(dataset, files):
    for key in POTENTIAL_KEYS:
        lazy = sc.read_h5ad(files[0], backed="r")
        present = key in lazy.obs.columns
        lazy.file.close()
        if not present:
            continue
        for value in common_values(files, key):
            clean = re.sub(r"[^0-9A-Za-z]+", "_", value).strip("_")
            run_analysis(dataset, files, f"Full_trajectory_{key}_{clean}",
                         subset_key=key, subset_val=value)


def main():
    for dataset in RUN_DATASETS:
        files = sorted(glob.glob(os.path.join(BASE_INPUT_DIR, dataset, ".h5ad")), key=tp_number)
        run_analysis(dataset, files, "Full_trajectory")
        run_subsets(dataset, files)


if __name__ == "__main__":
    main()
