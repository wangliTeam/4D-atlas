# ============================================================================
#gene regulatory network (GRN) inference
# Ref
#https://github.com/gpenglab/SpinalCordInjury/blob/v1.0.0/script/Co-expression/RM/pyscienic.py
import os
import glob
import pickle
import pandas as pd
import numpy as np
import re
import time
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.export import export2loom, add_scenic_metadata
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from distributed import Client, LocalCluster
import hotspot # Hotspot
import seaborn as sns
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpb
from sklearn import mixture
import plotly.express as px
import plotly as py
from typing import Optional, Mapping

ncores = 8
nthreads = 8
RANDOM_SEED = 42
MIN_MOD_THRESHOLD = 5 
motif_file_map = {
    "dmel": "motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl",
    "hg38": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    "mm": "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
}
ranking_file_map = {
    "dmel": "dm6_v10_clust.genes_vs_motifs.rankings.feather",
    "hg38": "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    "mm": "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
}
tf_file_map = {
    "dmel": "allTFs_dmel.txt",
    "hg38": "allTFs_hg38.txt",
    "mm": "allTFs_mm.txt"
}

def _derive_threshold(auc_mtx: pd.DataFrame, regulon_name: str) -> float:
    assert auc_mtx is not None and not auc_mtx.empty
    assert regulon_name in auc_mtx.columns
    data = auc_mtx[regulon_name].values.reshape(-1, 1)
    gmm = mixture.GaussianMixture(n_components=2, covariance_type='full').fit(data)
    avgs = gmm.means_
    stds = np.sqrt(gmm.covariances_.reshape(-1, 1))

    idx = np.argmax(avgs)
    threshold = max(avgs[idx] - 2 * stds[idx], 0)
    idx = np.argmin(avgs)
    lower_bound = avgs[idx] + 2 * stds[idx]

    return max(lower_bound, threshold)

def binarize(auc_mtx: pd.DataFrame, threshold_overides:Optional[Mapping[str,float]]=None) -> (pd.DataFrame, pd.Series):
    def derive_thresholds(auc_mtx):
        return pd.Series(index=auc_mtx.columns, data=[_derive_threshold(auc_mtx, name) for name in auc_mtx.columns])
    thresholds = derive_thresholds(auc_mtx)
    if threshold_overides is not None:
        thresholds[list(threshold_overides.keys())] = list(threshold_overides.values())
    return (auc_mtx > thresholds).astype(int), thresholds

def plot_binarization(auc_mtx: pd.DataFrame, regulon_name: str, bins: int=200, threshold=None, ax=None) -> None:
    if ax is None:
        ax=plt.gca()
    auc_mtx[regulon_name].hist(bins=bins,ax=ax)
    if threshold is None:
        threshold = _derive_threshold(auc_mtx, regulon_name)

    ylim = ax.get_ylim()
    ax.plot([threshold]*2, ylim, 'r:')
    ax.set_ylim(ylim)
    ax.set_xlabel('AUC')
    ax.set_ylabel('#')
    ax.set_title(regulon_name)

if __name__ == '__main__':
    RESOURCES_FOLDER = "../Motif_to_TF_annotation_databases/"
    DATABASE_FOLDER = "../Ranked_whole_genome_databases/"
    MM_TF_FOLDER = "../Transcription_factor_gene list/"
    input_data_path = "../Data"
    output_base_path = '../3DGRN'
    h5ad_files = sorted(glob.glob(os.path.join(input_data_path, '*_alignment.h5ad')))
    for in_h5ad in h5ad_files:
        base_name = os.path.basename(in_h5ad)
        timepont_label = base_name.replace('_alignment.h5ad', '')
        output_GRN_path = os.path.join(output_base_path, timepont_label)
        os.makedirs(output_GRN_path, exist_ok=True)
        os.chdir(output_GRN_path)
        REGULONS_FNAME = "regulons.p"
        MOTIFS_FNAME = "motifs.csv"
        REGULONS_DF_FNAME = "regulons.csv"
        AUCMTX_FNAME = "auc_mtx.csv"
        LOOM_FILE = "pyscenic.auc.loom"
        HS_PICKLE_FNAME = f"{timepont_label}.hs"
        DATABASES_GLOB = os.path.join(DATABASE_FOLDER, ranking_file_map[species])
        MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, motif_file_map[species])
        MM_TFS_FNAME = os.path.join(MM_TF_FOLDER, tf_file_map[species])
        adata = sc.read_h5ad(in_h5ad)
        exp = pd.DataFrame(adata.layers[counts_key])
        exp.index = adata.obs.index.tolist()
        gene_name = [i.replace("\"", "") for i in adata.var.index.tolist()]
        exp.columns = gene_name
        x_coords = adata.obsm['spatial_3D'][:, 0].tolist()
        y_coords = adata.obsm['spatial_3D'][:, 1].tolist()
        z_coords = adata.obsm['spatial_3D'][:, 2].tolist()
        coords=pd.DataFrame({'physical_x':x_coords,'physical_y':y_coords,'physical_z':z_coords})
        coords.index=adata.obs.index.tolist()
        slice_ID = adata.obs['slice_ID'].tolist()
        anno_data = adata.obs[cluster_key].tolist()
        tissue_annoData=pd.DataFrame({'annotation':anno_data,'slice_ID':slice_ID})
        tissue_annoData.index=adata.obs.index.tolist()
        exp.T.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_expdata.csv"))
        coords.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_phycoords.csv"))
        tissue_annoData.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_tissue_annoData.csv"))
        ex_matrix = pd.read_csv(os.path.join(output_GRN_path, f"{timepont_label}_expdata.csv"), index_col=0).T
        tf_names_full = load_tf_names(MM_TFS_FNAME)
        tf_names = [i for i in tf_names_full if i in ex_matrix.columns.values]
        db_fnames = glob.glob(DATABASES_GLOB)
        def name(fname):
            return os.path.splitext(os.path.basename(fname))[0]
        dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
        cluster = LocalCluster(n_workers=ncores, threads_per_worker=nthreads)
        adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
        adjacencies.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_adjacencies.csv"))
        modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
        modules = np.array(modules)
        np.save(os.path.join(output_GRN_path, f"{timepont_label}_modules.npy"), modules)
        with ProgressBar():
            modules = list(modules)
            df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
        regulons = df2regulons(df)
        df.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_{MOTIFS_FNAME}"))
        with open(REGULONS_FNAME, "wb") as f:
            pickle.dump(regulons, f)
        with open(REGULONS_FNAME, "rb") as f:
            regulons = pickle.load(f)
        name = []; tfs = []; targets = []; score = []; motif = []
        for i in regulons:
            name.append(i.name)
            tfs.append(i.transcription_factor)
            targets.append(','.join(i.genes))
            score.append(i.score)
            ct = list(i.context)
            if 'png' in ct[0]:
                motif.append(ct[0].split('.')[0])
            elif 'png' in ct[1]:
                motif.append(ct[1].split('.')[0])
            else:
                motif.append('')

        regulons_df = pd.DataFrame(data={'name': name, 'tfs': tfs, 'score': score, 'targets': targets, 'motif': motif})
        regulons_df.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_{REGULONS_DF_FNAME}"), index=False)
        auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
        auc_mtx.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_{AUCMTX_FNAME}"))
        auc_bin, thresholds = binarize(auc_mtx)
        auc_bin.to_csv(os.path.join(output_GRN_path, f"{timepont_label}_binarizedbyPyscenic.auc_mtx.csv"))
        plot_df = coords.join(auc_mtx)
        pltdata = tissue_annoData.join(plot_df)
        tf_columns = [col for col in pltdata.columns if col not in ['annotation', 'slice_ID', 'physical_x', 'physical_y', 'physical_z']]
        grouped = pltdata.groupby('annotation')[tf_columns].mean()
        fig, ax = plt.subplots(figsize=(max(12, len(grouped.index)*0.5), max(6, len(grouped.columns)*0.2)))
        sns.heatmap(grouped.T, cmap='Purples', linewidths=0.5, linecolor='gray', ax=ax)
        plt.title(f"{timepont_label} Average TF Activity")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        fig.savefig(os.path.join(output_GRN_path, f"{timepont_label}_TF_activity_heatmap.png"), dpi=300)
        plt.close(fig)
        long_df = grouped.reset_index().melt(id_vars='annotation', var_name='TF', value_name='mean_activity')
        long_df['size'] = (long_df['mean_activity'] - long_df['mean_activity'].min()) / (long_df['mean_activity'].max() - long_df['mean_activity'].min()) * 300 + 20
        fig, ax = plt.subplots(figsize=(max(12, len(grouped.index)*0.5), max(8, len(grouped.columns)*0.25)))
        sns.scatterplot(data=long_df, x='annotation', y='TF',
                        size='size', sizes=(20, 300),
                        hue='mean_activity', palette='Purples', ax=ax)
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)
        plt.title(f"{timepont_label} TF Activity per Annotation (Bubble Plot)")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        fig.savefig(os.path.join(output_GRN_path, f"{timepont_label}_TF_activity_bubbleplot.png"), dpi=300)
        plt.close(fig)

