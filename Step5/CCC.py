# Ref
# https://spateo-release.readthedocs.io/en/latest/tutorials/notebooks/6_cci/1_cell-cell_communication_LR_coexpression_analysis.html

import os
import glob
import re
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import gc
import warnings
aligned_adata = sc.read_h5ad(file_path)
slice_ids = sorted(list(set(aligned_adata.obs['slice_ID'].values)), key=lambda x: int(x.replace('slice', '')))
print(slice_ids)
for slice_id in slice_ids:
    first_slice = aligned_adata[aligned_adata.obs['slice_ID'] == slice_id].copy()
    first_slice.uns["__type"] = "UMI"
    _, first_slice = st.tl.neighbors(first_slice, basis='spatial', spatial_key='align_spatial', n_neighbors=10)

    labels = first_slice.obs[cluster_key].unique()
    palette = sns.color_palette("tab10", len(labels))
    color_map = dict(zip(labels, palette.as_hex()))
    first_slice.uns['color_key'] = color_map
    st.plotting.plot_connections(
        first_slice,
        cat_key=cluster_key,
        save_show_or_return='save',
        save_kwargs={
            # "path": tp_output_dir,
            # "prefix": f"2D_{tp_label}_{slice_id}_connection",
            "prefix": os.path.join(tp_output_dir, f"2D_{tp_label}_{slice_id}_connection"),
            "dpi": 300,
            "ext": "png",
            "transparent": False,
            "close": True,
            "verbose": True
        },
        title_str=f"2D_{tp_label}_{slice_id}",
        title_fontsize=8,
        label_fontsize=6,
        figsize=(4, 4)
    )

    counts = first_slice.obs[cluster_key].value_counts()
    a = counts[counts >= 30].index.tolist()
    df = pd.DataFrame({
        "celltype_sender": np.repeat(a, len(a)),
        "celltype_receiver": list(a) * len(a),
    })
    df = df[df['celltype_sender'] != df['celltype_receiver']]
    df["celltype_pair"] = df["celltype_sender"].str.cat(df["celltype_receiver"], sep="-")
    df = df.reset_index(drop=True)

    res = {}
    dropped = []

    for idx, i in enumerate(df['celltype_pair']):
        s, r = i.split(sep='-')
        st.tl.prepare_cci_cellpair_adata(first_slice, sender_group=s, receiver_group=r, group=cluster_key, all_cell_pair=True)
        st.pl.space(
                first_slice,
                color=['spec'],
                pointsize=0.2,
                color_key={'other': '#D3D3D3', s: first_slice.uns['color_key'][s], r: first_slice.uns['color_key'][r]},
                show_legend='upper left',
                figsize=(4, 4),
                save_show_or_return='save',
                save_kwargs={
                    "prefix": os.path.join(tp_output_dir, f"2D_{tp_label}_{slice_id}_{s}_{r}_pair"),
                    "ext": "png"
                }
            )

        result = st.tl.find_cci_two_group(
                first_slice,
                path=db_dir,
                species=species,
                group=cluster_key,
                sender_group=s,
                receiver_group=r,
                filter_lr='outer',
                min_pairs=0,
                min_pairs_ratio=0,
                top=20,
            )
        if result is not None:
                res[i] = result
        else:
                dropped.append(idx)

        result = pd.DataFrame(columns=res[df['celltype_pair'][1]]['lr_pair'].columns)
        for l in df.index:
            if l not in dropped:
                res[df['celltype_pair'][l]]['lr_pair'] = res[df['celltype_pair'][l]
                                                             ]['lr_pair'].sort_values('lr_co_exp_ratio', ascending=False)[0:3]
                result = pd.concat([result, res[df['celltype_pair'][l]]
                                    ['lr_pair']], axis=0, join='outer')

        df_result = result.loc[result['lr_co_exp_num'] > 5]
        df_result.drop_duplicates(
            subset=['lr_pair', 'sr_pair', ], keep='first', inplace=True)
        df_result.to_csv(os.path.join(tp_output_dir, f"{tp_label}_{slice_id}_df_results.csv"), index=False)
