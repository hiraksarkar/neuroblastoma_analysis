## Data descriptions 
The analysis contains data from different samples. Briefly the dataset contains
### Mouse NB (New)
#### Tumor
The tumor samples - `/home/meisl/Workplace/neuroblastoma/ThMYCN_210610/conos/conos/all.tumor_conos.rds`
```
'het1_tumor''het2_tumor''homo1_tumor''homo2_tumor''homo3_tumor'
```
#### BM (should contain one more sample)
The BM samples - `/home/meisl/Workplace/neuroblastoma/ThMYCN_210610/conos/conos/all.bm_conos.rds`
```
'het1_bm''het2_bm''homo1_bm''homo2_bm''wt1_bm''wt2_bm''wt3_bm'
```
### Human NB (old)

The viable cells 
NB01-26, NB34, NB36 `/home/meisl/Workplace/neuroblastoma/Figures/data/scon.new.conos.rds`

Annotation `/home/meisl/Workplace/neuroblastoma/Figures/data/newano.rds`

Nuc-Seq `/home/meisl/Workplace/neuroblastoma/Figures/data/nucseq`

### Public data
#### Mouse adrenal samples 
```
'adrenal_70' = '/home/hsarkar/Projects/ThMYCN/PMID33833454/GSM5067113_10x70.raw_feature_bc_matrix.h5',
'adrenal_71' = '/home/hsarkar/Projects/ThMYCN/PMID33833454/GSM5067114_10x71.raw_feature_bc_matrix.h5',
'adrenal_72' = '/home/hsarkar/Projects/ThMYCN/PMID33833454/GSM5067115_10x72.raw_feature_bc_matrix.h5'
```

### Notebooks
- [Initial Seurat Analysis](https://github.com/hiraksarkar/neuroblastoma_analysis/blob/master/notebook/Seurat_Analysis.ipynb)
- [Cell-cycle Analysis](https://github.com/hiraksarkar/neuroblastoma_analysis/blob/master/notebook/Cell_Cycle_Regression.ipynb)
- [Velocity Analysis](https://github.com/hiraksarkar/neuroblastoma_analysis/blob/master/notebook/ThMYCN_velocity.ipynb)
