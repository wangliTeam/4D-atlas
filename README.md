# 4D-Atlas
A cross-species spatiotemporal (4D) transcriptomics pipeline.

**Step 1 — 3D reconstruction**: pairwise morphological alignment of tissue slices and z-coordinate assignment with spateo to build a 3D model.

**Step 2 — Spatial modules**: Hotspot detection of spatially variable genes and their co-expression modules.

**Step 3 — Pathway activity**: decoupler inference of PROGENy signalling pathway and MSigDB Hallmark gene-set activity.

**Step 4 — Gene regulatory network**: pySCENIC regulon inference and AUCell activity scoring.

**Step 5 — Cell–cell communication**: spateo ligand–receptor co-expression analysis across slices.

**Step 6 — Trajectory and key regulators**: CellRank GPCCA fate probabilities and lineage drivers; TOME state transitions and key transcription factors.

**Step 7 — Temporal clustering**: Mfuzz soft clustering of time-resolved expression profiles.
