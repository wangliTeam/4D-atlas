# Ref:https://decoupler.readthedocs.io/
import os
import sys
import json
import warnings
import argparse
import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
import scipy.sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gc
import logging
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Ignore warnings
warnings.filterwarnings("ignore")

adata = sc.read_h5ad(h5ad_path)
progeny = dc.op.progeny(organism="human")
progeny
dc.mt.ulm(data=adata, net=progeny)
score = dc.pp.get_obsm(adata=adata, key="score_ulm")
score
sc.pl.matrixplot(
    adata=score,
    var_names=score.var_names,
    groupby="celltype",
    dendrogram=True,
    standard_scale="var",
    colorbar_title="Z-scaled scores",
    cmap="RdBu_r",
)
hallmark = dc.op.hallmark(organism="human")
hallmark
dc.mt.ulm(data=adata, net=hallmark)
score = dc.pp.get_obsm(adata=adata, key="score_ulm")
score
df = dc.tl.rankby_group(adata=score, groupby="celltype", reference="rest", method="t-test_overestim_var")
df = df[df["stat"] > 0]
df
n_markers = 3
source_markers = (
    df.groupby("group")
    .head(n_markers)
    .drop_duplicates("name")
    .groupby("group")["name"]
    .apply(lambda x: list(x))
    .to_dict()
)
source_markers
sc.pl.matrixplot(
    adata=score,
    var_names=source_markers,
    groupby="celltype",
    dendrogram=True,
    standard_scale="var",
    colorbar_title="Z-scaled scores",
    cmap="RdBu_r",
    swap_axes=True,
)
