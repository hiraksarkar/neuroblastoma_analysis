suppressPackageStartupMessages({
    library(tidyr)
    library(patchwork)
    library(Seurat)
    library(conos)  
    library(dplyr)
    library(ggplot2)
    library(infercnv)
    library(logger)
})

log_info('Reading count matrix...')

count_matrix = readRDS(
    glue::glue(
        '/home/hsarkar/Projects/neuroblastoma_analysis/',
        'infercnv_files/',
        'mouse_nb_adrenal_8_samples_mat.rds'
    )
)

log_info('Count matrix read...')

log_info('Creating infercnv object...')
# create the object
infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix=count_matrix,
    annotations_file= glue::glue(
            '/home/hsarkar/Projects/neuroblastoma_analysis/',
            'infercnv_files/',
            'mouse_nb_adrenal_8_samples_seurat_cc_regressed_annot.txt'
        ),
    delim="\t",
    gene_order_file = glue::glue(
            '/home/hsarkar/Projects/neuroblastoma_analysis/',
            'infercnv_files/mycn_gene_pos_mod.txt'
    ),
    ref_group_names = c(
        'B cells',
        'T cells',
        'Myeloid',
        'Neutrophils',
        'DC'
    )
) 

log_info('Running infercnv...')
infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff=0.5, 
    out_dir=glue::glue(
        '/home/hsarkar/Projects/neuroblastoma_analysis/',
        'infercnv_files/out_Jan_2/'
    ),
    cluster_by_groups=TRUE, 
    denoise=TRUE,
    num_threads = 30,
    HMM=TRUE,
    output_format = "pdf"
)


log_info('Storing infercnv object...')
saveRDS(
    infercnv_obj,
    glue::glue(
        '/home/hsarkar/Projects/neuroblastoma_analysis/',
        'infercnv_files/out_Jan_2/',
        'mouse_nb_adrenal_8_samples_seurat_infercnv.rds'
    )
)