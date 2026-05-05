# Materials and Methods

## Single-Cell RNA Sequencing Data Processing and Analysis

### Data Preprocessing

Single-cell RNA sequencing data were processed using Scanpy v1.9+ (Wolf et al., 2018) in Python 3.8+. Raw count matrices in either 10X Genomics HDF5 or Matrix Market (MTX) format were imported for individual samples (C3L-03405, C3L-03968, C3N-01334, C3N-02190, C3N-02784, and C3N-03188). Gene identifiers were converted to gene symbols, and duplicate gene names were resolved by making variable names unique.

### Quality Control

Quality control metrics were calculated using `sc.pp.calculate_qc_metrics()` to assess the number of genes expressed per cell, total counts per cell, and the percentage of mitochondrial gene expression. Cells were filtered based on the following criteria:
- Minimum of 200 genes expressed per cell
- Maximum of 10,000 genes expressed per cell
- Minimum of 1,000 UMI counts per cell
- Maximum of 10,000 total UMI counts per cell
- Maximum of 10% mitochondrial gene content

Genes expressed in fewer than 3 cells were excluded from downstream analysis. Mitochondrial genes were identified by the "MT-" prefix in gene symbols.

### Sample Integration

After individual sample preprocessing, filtered AnnData objects were merged using `sc.concat()` with an outer join strategy, preserving sample identity through a "sample_id" label. The merged dataset underwent an additional round of quality control using the same filtering criteria to ensure consistency across samples.

### Normalization and Highly Variable Gene Selection

Library sizes were normalized to 10,000 counts per cell using `sc.pp.normalize_total()`, followed by log-transformation (natural logarithm) using `sc.pp.log1p()`. Highly variable genes (HVGs) were identified using the Seurat flavor method (Stuart et al., 2019) implemented in Scanpy, selecting the top 3,000 most variable genes across samples using batch-aware HVG detection with the "sample" batch key to account for sample-specific variation.

### Data Scaling and Regression

Technical variation was mitigated by regressing out the effects of total counts per cell and the percentage of mitochondrial gene expression using `sc.pp.regress_out()`. Cell cycle effects were visualized using cell cycle scores calculated with `sc.tl.score_genes_cell_cycle()` based on S phase and G2/M phase marker genes from Regev lab (Tirosh et al., 2016), but cell cycle effects were not regressed out in the final merged dataset analysis to preserve biological variation.

### Dimensionality Reduction and Clustering

Principal component analysis (PCA) was performed on the scaled and log-normalized expression matrix of highly variable genes using `sc.tl.pca()` with the ARPACK SVD solver, computing 30 principal components. A nearest neighbor graph was constructed using `sc.pp.neighbors()` with 20 neighbors based on the first 30 principal components. Uniform Manifold Approximation and Projection (UMAP) was computed using `sc.tl.umap()` for visualization of the data in two dimensions.

Unsupervised clustering was performed using the Leiden algorithm (Traag et al., 2019) with a resolution parameter of 0.35, implemented via `sc.tl.leiden()`.

### Reference Mapping with scANVI

Cell type annotation was performed using supervised transfer learning with scANVI (Xu et al., 2021), implemented in scvi-tools v1.1.2+ (Gayoso et al., 2022). The GBmap core reference atlas, a comprehensive glioblastoma single-cell reference containing multiple annotated cell types across different studies, was used as the reference dataset.

#### Reference Model Training

The reference dataset was preprocessed to ensure raw counts were used for model training. A single-cell variational inference (scVI) model was trained on the reference atlas using the following hyperparameters:
- Batch key: "author" (to account for study-specific batch effects)
- Layer normalization: "both"
- Batch normalization: "none"
- Covariate encoding: True
- Dropout rate: 0.2
- Number of hidden layers: 2

The trained scVI model was then used to initialize a scANVI model using `scvi.model.SCANVI.from_scvi_model()` with the "annotation_level_3" cell type labels from the reference and "Unknown" as the unlabeled category. The scANVI reference model was trained for 200 epochs with early stopping enabled, checking validation metrics every 25 epochs.

#### Query Dataset Mapping

The preprocessed query dataset (merged samples) was prepared for reference mapping using `scvi.model.SCANVI.prepare_query_anndata()`, which aligned gene sets between query and reference datasets and filtered to retain only shared genes. A query-specific scANVI model was initialized and trained for 150 epochs with the following settings:
- Learning rate: 1×10⁻⁴
- Weight decay: 0.0
- Maximum epochs: 150

Cell type predictions were obtained using `scanvi_query.predict()`, and latent representations were extracted using `scanvi_query.get_latent_representation()`. The scANVI predictions were stored in the AnnData object as "predictions_scanvi" in the `.obs` metadata.

### Computational Environment

All analyses were performed using Python 3.8+ with the following key packages: scanpy (v1.9+), scvi-tools (v1.1.2+), anndata, numpy, pandas, and matplotlib. Parallel processing was enabled using 6 CPU cores for individual sample processing. Reference mapping was accelerated using GPU computing with PyTorch CUDA support. Random seeds were set (seed=0) for reproducibility of scVI/scANVI model training.

### Data Availability

Processed AnnData objects were saved at key processing stages:
- Individual samples after quality control: `{sample}_filtered.h5ad`
- Merged samples after quality control: `samples_merged_filtered_mg_mc_gc_mt_tc.h5ad`
- scANVI-annotated dataset: `All_samples_scANVI_annotated.h5ad`
- Final processed and clustered dataset: `All_samples_pp_clustered.h5ad`

## References

Gayoso, A., Lopez, R., Xing, G., Boyeau, P., Valiollah Pour Amiri, V., Hong, J., Wu, K., Jayasuriya, M., Mehlman, E., Langevin, M., Liu, Y., Samaran, J., Misrachi, G., Nazaret, A., Clivio, O., Xu, C., Ashuach, T., Gabitto, M., Lotfollahi, M., ... Yosef, N. (2022). A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology*, 40(2), 163-166.

Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W.M., Hao, Y., Stoeckius, M., Smibert, P., & Satija, R. (2019). Comprehensive Integration of Single-Cell Data. *Cell*, 177(7), 1888-1902.

Tirosh, I., Izar, B., Prakadan, S.M., Wadsworth, M.H., Treacy, D., Trombetta, J.J., Rotem, A., Rodman, C., Lian, C., Murphy, G., Fallahi-Sichani, M., Dutton-Regester, K., Lin, J.R., Cohen, O., Shah, P., Lu, D., Genshaft, A.S., Hughes, T.K., Ziegler, C.G., ... Garraway, L.A. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science*, 352(6282), 189-196.

Traag, V.A., Waltman, L., & van Eck, N.J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9(1), 5233.

Wolf, F.A., Angerer, P., & Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19(1), 15.

Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M.I., & Yosef, N. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. *Molecular Systems Biology*, 17(1), e9620.
