# scCrossTalk
Infer cell-cell communications based on CellTalkDB





# Install
```
# download the source package of scCrossTalk-1.0.tar.gz and install it
# ensure the right directory for scCrossTalk-1.0.tar.gz
install.packages(pkgs = 'scCrossTalk-1.0.tar.gz',repos = NULL, type = "source")
```
or
```
# install devtools and install scCrossTalk
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCrossTalk')
```

# Usage
`library(scCrossTalk)`
### Find highly expressed ligand-receptor pairs
Find highly expressed ligands and receptors between pairwise clusters using Z score for a `Seurat` object (>= 3.0.0) after log1p normalization, cluster analysis and tSNE or Umap dimensionality reduction
```
clu_pairs <- FindPairs(object = mouse_kidney_203_Seurat,
                       species = "Mouse",
                       use_LRdb = "LRdb",
                       revise_gene = T,
                       cell_min_pct = 0.25,
                       p_value = 0.05)
```

