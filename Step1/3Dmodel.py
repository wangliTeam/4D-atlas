# ============================================================================
# Ref
# https://github.com/Flysta3D/multi_omics_atlas/blob/main/Stereo-seq_Basic_analysis/03.Alignment.py
#https://spateo-release.readthedocs.io/en/latest/tutorials/notebooks/3_alignment/5.2%203D%20reconstruction%20on%20Drosophila%20data%20with%20pairwise%20alignment.html

import os
import re
import math
import glob
import gc
import warnings
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import spateo as st
import scipy.sparse
import matplotlib.pyplot as plt
from typing import Union, List, Optional
from anndata import AnnData

# Default slice thickness
DEFAULT_THICKNESS = 10 
# # Custom uniform spacing
# st.tl.assign_z_coordinates(slices, z_spacing=10.0)
# # Use known tissue section thickness
# st.tl.assign_z_coordinates(slices, tissue_thickness=15.0)  # 15 µm sections
# # Handle missing sections or variable spacing
# st.tl.assign_z_coordinates(slices, z_spacing=[10.0, 30.0, 10.0])  # Larger gap for missing section

slices = [slice1, slice2, slice3, slice4]
spatial_key = 'spatial'
cluster_key = 'tissue'
align_key = "align_spatial"
st.pl.slices_2d(
    slices = slices,
    label_key = cluster_key,
    spatial_key = spatial_key,
    height=2,
    center_coordinate=True,
    show_legend=True,
    legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, 0) ,'ncol': 5, 'borderaxespad': -6, 'frameon': False},
    palette=palette,
    save_show_or_return='save',
    save_kwargs={"path": output_3d_model_dir, "prefix": f"{tp_label}_{current_key}_slices_2d_pre", "ext": "png", "dpi": 300, "bbox_inches": "tight"}
   
)

for i in range(1, len(slices)):
    st.tl.morpho_align(sliceA=slices[i], sliceB=slices[i-1], 
                       spatial_key=spatial_key, key_added=align_key)

st.tl.assign_z_coordinates(slices, spatial_key=align_key, 
                           tissue_thickness=DEFAULT_THICKNESS)
st.pl.slices_2d(
    slices=slices, 
    label_key=cluster_key, 
    spatial_key=align_key, 
    height=2,
    show_legend=True, 
    return_palette=True, 
    save_show_or_return='save',
    save_kwargs={"path": output_3d_model_dir, "prefix": f"{tp_label}_{current_key}_slices_2d", "ext": "png", "dpi": 300, "bbox_inches": "tight"}
    )
st.pl.overlay_slices_2d(
                slices=aligned_slices, spatial_key=align_key, 
                ncols=optimal_ncols,
                height=2, overlay_type='both', save_show_or_return='save',
                save_kwargs={"path": output_3d_model_dir, "prefix": f"{tp_label}_{current_key}_overlay", "ext": "png", "dpi": 300, "bbox_inches": "tight"}
            )

st.pl.three_d_multi_plot.multi_models(*slices, spatial_key=align_key)
