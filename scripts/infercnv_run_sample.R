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

samples = c('het2_tumor','homo1_tumor','homo2_tumor','homo3_tumor')
for(sample_name in samples){
    count_matrix = readRDS(
        glue::glue(
                    '/home/hsarkar/Projects/neuroblastoma_analysis/',
                    'infercnv_files/',
                    'mouse_nb_{sample_name}_count_matrix.rmd'
        )
    )

    infercnv_obj = CreateInfercnvObject(
        raw_counts_matrix=count_matrix,
        annotations_file= glue::glue(
                '/home/hsarkar/Projects/neuroblastoma_analysis/',
                'infercnv_files/',
                'mouse_nb_{sample_name}_annot.txt'
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

    dir.create(file.path(
    glue::glue(
        '/home/hsarkar/Projects/neuroblastoma_analysis/',
        'infercnv_files/out_Jan_4/{sample_name}'
    )
    ), recursive = TRUE)

    infercnv_obj = infercnv::run(
        infercnv_obj,
        cutoff=0.5, 
        out_dir=glue::glue(
            '/home/hsarkar/Projects/neuroblastoma_analysis/',
            'infercnv_files/out_Jan_4/{sample_name}/'
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
            'infercnv_files/out_Jan_4/',
            '{sample_name}.rds'
        )
    )

}
