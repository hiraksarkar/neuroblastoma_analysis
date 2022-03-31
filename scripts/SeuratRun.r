library(Seurat)
library(SeuratData)
library(patchwork)
library(conos)
library(dplyr)

## read the 5 samples 
###---------------------------------------------
Conos2Seurat = function(
    con, 
    prefix=NULL, 
    nfeatures=2000,
    remove.cc.genes = FALSE 
){
    raw.count.objs = lapply(con$samples, function(x){
        raw.count = x$misc$rawCounts
        mito.genes = grep(pattern = "^mt-", x = colnames(raw.count), value = TRUE, ignore.case = TRUE)
        rpl.genes = grep(pattern = "^rpl", x = colnames(raw.count), value = TRUE, ignore.case = TRUE)
        rps.genes = grep(pattern = "^rps", x = colnames(raw.count), value = TRUE, ignore.case = TRUE)
        genes.removed = union(mito.genes, rpl.genes)
        genes.removed = union(genes.removed, rps.genes)
        
        raw.count = raw.count[, setdiff(colnames(raw.count), genes.removed) ]
        colnames(raw.count) = toupper(colnames(raw.count))
        t(raw.count)

    }
    )

    seurat.obj.list = lapply(names(raw.count.objs),function(sample.name)  CreateSeuratObject(
        counts = raw.count.objs[[sample.name]], 
        #project = paste("mnb5",sample.name, sep="_"), 
        min.cells = 3, 
        min.features = 500
    ))

    names(seurat.obj.list) = names(raw.count.objs)

    if (remove.cc.genes == TRUE){
        s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes

        s.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  s.genes )
        g2m.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  g2m.genes )
        cc.genes = union( s.genes.selected, g2m.genes.selected )
        seurat.obj.list = lapply(X = seurat.obj.list, FUN = function(x){
            x = x[setdiff( rownames(x), cc.genes),]
        })
    }

    seurat.obj.list <- lapply(X = seurat.obj.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
    })

    seurat.obj.list
}

CombineSeurat = function(obj.list.1, obj.list.2){
    n.sample = length(obj.list.1)
    combined.list = merge(
        obj.list.1[[1]],
        y = c(obj.list.1[2:n.sample], obj.list.2),
        merge.data = T
    )
    combined.list
}

IntegrateSeuratObjCC = function(
    seurat.obj.list, 
    nfeatures=2000,
    cc = TRUE,
    remove.cc.genes = FALSE
){
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes

    s.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  s.genes )
    g2m.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  g2m.genes )

    print('getting features ...')
    features <- SelectIntegrationFeatures(
        object.list = seurat.obj.list, 
        nfeatures = nfeatures
    )

    print('find anchors (slow) ...')
    tumor.anchors <- FindIntegrationAnchors(
        object.list = seurat.obj.list, 
        anchor.features = features
    )

    print('integrate ...')
    tumor.combined <- IntegrateData(
        anchorset = tumor.anchors
    )

    DefaultAssay(tumor.combined) <- "integrated"
    if (cc == TRUE){
        print('cell-cycle removal ...')
        tumor.combined <-  CellCycleScoring(
            tumor.combined, 
            s.features = s.genes.selected, 
            g2m.features = g2m.genes.selected, 
            set.ident = TRUE
        )

        tumor.combined <-  ScaleData(
            tumor.combined, 
            vars.to.regress = c("S.Score", "G2M.Score"),
            features = VariableFeatures(tumor.combined)
        )
    }

    tumor.combined = RunPCA(tumor.combined, npcs=30)
    tumor.combined = RunUMAP(tumor.combined, reduction="pca",dims=1:30)
    tumor.combined = FindNeighbors(tumor.combined, reduction="pca", dims=1:30)
    tumor.combined = FindClusters(tumor.combined, resolution = 0.5)

    return(tumor.combined)

}

IntegrateSeuratObjCCFast = function(
    seurat.obj.list, 
    nfeatures=2000,
    remove.cc.genes = FALSE
){
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes

    s.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  s.genes )
    g2m.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  g2m.genes )

    print('getting features ...')
    features <- SelectIntegrationFeatures(
        object.list = seurat.obj.list, 
        nfeatures = nfeatures
    )

    seurat.obj.list = lapply(X = seurat.obj.list, FUN = function(x){
        x = ScaleData(x, features = features, verbose = FALSE)
        x = RunPCA(x, features = features, verbose = FALSE)
    })

    print('find anchors (fast) ...')
    tumor.anchors <- FindIntegrationAnchors(
        object.list = seurat.obj.list,
        reference = c(1, 3),
        reduction = "rpca",
        dims = 1:50
    )

    print('integrate ...')
    tumor.combined <- IntegrateData(
        anchorset = tumor.anchors,
        dims = 1:50
    )

    DefaultAssay(tumor.combined) <- "integrated"
    print('cell-cycle removal ...')
    tumor.combined.cc <-  CellCycleScoring(
        tumor.combined, 
        s.features = s.genes.selected, 
        g2m.features = g2m.genes.selected, 
        set.ident = TRUE
    )

    tumor.combined.cc <-  ScaleData(
        tumor.combined.cc, 
        vars.to.regress = c("S.Score", "G2M.Score"),
        features = VariableFeatures(tumor.combined.cc)
    )

    return(list(
        cc = tumor.combined.cc,
        wocc = tumor.combined
        )
    )

    # tumor.combined.cc = RunPCA(tumor.combined.cc)
    # tumor.combined.cc = RunUMAP(tumor.combined.cc, reduction="pca",dims=1:50)
    # tumor.combined.cc = FindNeighbors(tumor.combined.cc, reduction="pca", dims=1:50)
    # tumor.combined.cc = FindClusters(tumor.combined.cc, resolution = 0.5)
    
    # tumor.combined = ScaleData(tumor.combined, verbose=FALSE)
    # tumor.combined = RunPCA(tumor.combined)
    # tumor.combined = RunUMAP(tumor.combined, reduction="pca",dims=1:50)
    # tumor.combined = FindNeighbors(tumor.combined, reduction="pca", dims=1:50)
    # tumor.combined = FindClusters(tumor.combined, resolution = 0.5)

    # return(list(
    #     cc = tumor.combined.cc,
    #     wocc = tumor.combined
    #     )
    # )

}

IntegrateSeuratObjFast = function(
    seurat.obj.list, 
    nfeatures=2000,
    remove.cc.genes = FALSE
){
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes

    s.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  s.genes )
    g2m.genes.selected = intersect(seurat.obj.list[[1]] %>% rownames,  g2m.genes )

    print('getting features ...')
    features <- SelectIntegrationFeatures(
        object.list = seurat.obj.list, 
        nfeatures = nfeatures
    )

    seurat.obj.list = lapply(X = seurat.obj.list, FUN = function(x){
        x = ScaleData(x, features = features, verbose = FALSE)
        x = RunPCA(x, features = features, verbose = FALSE)
    })

    print('find anchors (fast) ...')
    tumor.anchors <- FindIntegrationAnchors(
        object.list = seurat.obj.list,
        reference = c(1, 3),
        reduction = "rpca",
        dims = 1:50
    )

    print('integrate ...')
    tumor.combined <- IntegrateData(
        anchorset = tumor.anchors,
        dims = 1:50
    )

    DefaultAssay(tumor.combined) <- "integrated"

    tumor.combined = ScaleData(tumor.combined, verbose=FALSE)
    tumor.combined = RunPCA(tumor.combined)
    tumor.combined = RunUMAP(tumor.combined, reduction="pca",dims=1:50)
    tumor.combined = FindNeighbors(tumor.combined, reduction="pca", dims=1:50)
    tumor.combined = FindClusters(tumor.combined, resolution = 0.5)

    return(tumor.combined)

    # return(list(
    #     cc = tumor.combined.cc,
    #     wocc = tumor.combined
    #     )
    # )

}

covert_seurat_to_pagoda_app = function(
    SR,
    fname,
    app.title,
    assay_name = "RNA",
    n.cores=10
){
    library(pagoda2)
    library(igraph)
    
    p2 = basicP2proc(
        SR@assays$RNA@counts[ VariableFeatures(SR) ,], 
        n.cores = n.cores,
        get.largevis = FALSE,
        get.tsne = FALSE,
        n.cores = 10,
        min.cells.per.gene = 0,
        min.transcripts.per.cell = 0
    )

    go.env <- p2.generate.human.go(p2)
    p2$clusters$PCA$seurat_cluster = as.factor(SR@meta.data$seurat_cluster)
    names(p2$clusters$PCA$seurat_cluster) = rownames(SR@meta.data)
    p2$embeddings$PCA$tSNE = as.matrix(SR@reductions$umap@cell.embeddings)
    p2$clusters$PCA$timepoint = as.factor(SR@meta.data$orig.ident)
    names(p2$clusters$PCA$timepoint) = rownames(SR@meta.data)
    
    hdea <- p2$getHierarchicalDiffExpressionAspects(
        type='PCA',
        clusterName='seurat_cluster',
        z.threshold=3, 
        n.cores = 20
    )

    extraWebMetadata = NULL
    metadata.forweb <- list();
    metadata.forweb$timepoint <- p2.metadata.from.factor(
        p2$clusters$PCA$timepoint,
        displayname='timepoint'
    )
    metadata.forweb$leiden <- p2.metadata.from.factor(
        p2$clusters$PCA$seurat_cluster,
        displayname='seurat_cluster'
    )
    metadata.forweb$multilevel <- p2.metadata.from.factor(
        p2$clusters$PCA$multilevel,
        displayname='multilevel'
    )
    metadata.forweb <- c(metadata.forweb, extraWebMetadata)
    genesets <- hierDiffToGenesets(hdea)
    appmetadata = list(apptitle=app.title)

    p2w = make.p2.app(p2, 
        additionalMetadata = metadata.forweb, 
        geneSets = genesets, 
        dendrogramCellGroups = p2$clusters$PCA$seurat_cluster, 
        show.clusters=F, 
        appmetadata = appmetadata
    )

    p2w$serializeToStaticFast(binary.filename = fname)
    return(p2)
}

covert_seurat_to_pagoda_app2 = function(
    SR,
    fname,
    app.title,
    assay_name = "RNA",
    n.cores = 20
){
    library(pagoda2)
    library(igraph)
    
    p2 = basicP2proc(
        SR@assays$RNA@counts, 
        n.cores = n.cores,
        get.largevis = FALSE,
        get.tsne = FALSE
    )

    go.env <- p2.generate.human.go(p2)
    p2$clusters$PCA$seurat_cluster = as.factor(SR@meta.data$seurat_cluster)
    names(p2$clusters$PCA$seurat_cluster) = rownames(SR@meta.data)
    p2$embeddings$PCA$tSNE = as.matrix(SR@reductions$umap@cell.embeddings)
    p2$clusters$PCA$timepoint = as.factor(SR@meta.data$orig.ident)
    names(p2$clusters$PCA$timepoint) = rownames(SR@meta.data)
    
    hdea <- p2$getHierarchicalDiffExpressionAspects(
        type='PCA',
        clusterName='seurat_cluster',
        z.threshold=3, 
        n.cores = 20
    )

    extraWebMetadata = NULL
    metadata.forweb <- list();
    metadata.forweb$timepoint <- p2.metadata.from.factor(
        p2$clusters$PCA$timepoint,
        displayname='timepoint'
    )
    metadata.forweb$leiden <- p2.metadata.from.factor(
        p2$clusters$PCA$seurat_cluster,
        displayname='seurat_cluster'
    )
    metadata.forweb$multilevel <- p2.metadata.from.factor(
        p2$clusters$PCA$multilevel,
        displayname='multilevel'
    )
    metadata.forweb <- c(metadata.forweb, extraWebMetadata)
    genesets <- hierDiffToGenesets(hdea)
    appmetadata = list(apptitle=app.title)

    p2w = make.p2.app(p2, 
        additionalMetadata = metadata.forweb, 
        geneSets = genesets, 
        dendrogramCellGroups = p2$clusters$PCA$seurat_cluster, 
        show.clusters=F, 
        appmetadata = appmetadata
    )

    p2w$serializeToStaticFast(binary.filename = fname)
    return(p2)
}

store_seurat = function(obj,out_dir){
    obj_meta_data = merge(
        obj@meta.data,
        obj@reductions$umap@cell.embeddings %>% data.frame,
        by = "row.names"
    )
    dir.create(out_dir, recursive = TRUE)
    write.csv(
        obj_meta_data,
        glue::glue(
            out_dir,
            '/meta_data.csv'
        ),
        row.names = FALSE
    )
    write.csv(
        obj@reductions$pca@cell.embeddings,
        file = glue::glue(
            out_dir,
            '/pca.csv'
        ),
        quote = F,
        row.names = F
    )
    counts_matrix <- GetAssayData(obj, assay='RNA', slot='counts')
    writeMM(
        counts_matrix,
        file = glue::glue(
            out_dir,
            '/counts.mtx'
        )
    )
    write.table(
        data.frame('gene' = rownames(counts_matrix)),
        file = glue::glue(
            out_dir,
            '/genes.csv'
        ),
        quote = F,
        row.names = F,
        col.names = F
    )
    
}

store_count_matrix = function(count_obj,out_dir){

    dir.create(out_dir, recursive = TRUE)
    
    write.csv(
        colnames(count_obj),
        file = glue::glue(
            out_dir,
            '/cell_names.csv'
        ),
        quote = F,
        row.names = F
    )
    
    writeMM(
        count_obj,
        file = glue::glue(
            out_dir,
            '/counts.mtx'
        )
    )
    write.table(
        data.frame('gene' = rownames(count_obj)),
        file = glue::glue(
            out_dir,
            '/genes.csv'
        ),
        quote = F,
        row.names = F,
        col.names = F
    )
    
}

load_seurat = function(out_dir){
    meta_data = read.csv(
        glue::glue(
            out_dir,
            '/meta_data.csv'
        ),
        row.names = 1
    )

    pca_data = read.csv(
        glue::glue(
            out_dir,
            '/pca.csv'
        ),
        row.names = 1
    )
    
    count_mat = readMM(
        glue::glue(
            out_dir,
            '/counts.mtx'
        )
    )
    
    gene_names = read.csv(
        glue::glue(
            out_dir,
            '/genes.csv'
        ),
        header = F
    )$V1
    
    colnames(count_mat) = rownames(meta_data)
    rownames(count_mat) = gene_names
    
    combined.removed = CreateSeuratObject(
        counts = count_mat,
        meta.data = meta_data
    )
    
    pca_data = as.matrix(pca_data)
    rownames(pca_data) = rownames(meta_data)
    
    combined.removed <- ScaleData(combined.removed, verbose = FALSE)
    combined.removed <- RunPCA(combined.removed, npcs = 30, verbose = FALSE)
    combined.removed <- RunUMAP(combined.removed, reduction = "pca", dims = 1:30)
    
    combined.removed@reductions$pca@cell.embeddings = pca_data
    combined.removed@reductions$umap@cell.embeddings = 
        as.matrix(meta_data[,c('UMAP_1','UMAP_2')])
    
    return(combined.removed)
    
}

save_seurat_as_pagoda=function(SR, fout, app.title = 'adrenal_sr_10K'){
    library(pagoda2)
    library(igraph)

    p2 <- basicP2proc(SR@assays$RNA@counts, n.cores = 20)

    go.env <- p2.generate.human.go(p2)


    #p2$clusters$PCA$seurat_cluster = as.factor(SR@meta.data$seurat_cluster)
    p2$clusters$PCA$seurat_cluster = as.factor(SR@meta.data$integrated_snn_res.0.5)
    
    names(p2$clusters$PCA$seurat_cluster) = rownames(SR@meta.data)

    p2$embeddings$PCA$tSNE = as.matrix(SR@reductions$umap@cell.embeddings)
        #p2$embeddings$PCA = as.matrix(p2$embeddings$PCA@cell.embeddings)

    p2$clusters$PCA$timepoint = as.factor(SR@meta.data$orig.ident)
    names(p2$clusters$PCA$timepoint) = rownames(SR@meta.data)

    p2$embeddings$PCA$tSNE = as.matrix(SR@reductions$umap@cell.embeddings)
        #p2$embeddings$PCA = as.matrix(p2$embeddings$PCA@cell.embeddings)


    n.cores=20

    cat('Calculating hdea...\n')
    hdea <- p2$getHierarchicalDiffExpressionAspects(
        type='PCA',
        clusterName='multilevel',
        z.threshold=3, 
        n.cores = n.cores)

    extraWebMetadata = NULL
    
    metadata.forweb <- list();
    metadata.forweb$timepoint <- p2.metadata.from.factor(
        p2$clusters$PCA$timepoint,
        displayname='timepoint'
    )
    metadata.forweb$leiden <- p2.metadata.from.factor(
        p2$clusters$PCA$seurat_cluster,
        displayname='seurat_cluster'
    )
    metadata.forweb$multilevel <- p2.metadata.from.factor(
        p2$clusters$PCA$multilevel,
        displayname='multilevel'
    )
    metadata.forweb <- c(metadata.forweb, extraWebMetadata)
    genesets <- hierDiffToGenesets(hdea)
    appmetadata = list(apptitle=app.title)
    cat('Making KNN graph...\n')
    #p2$makeGeneKnnGraph(n.cores=n.cores)
    p2w = make.p2.app(
        p2, 
        additionalMetadata = metadata.forweb, 
        geneSets = genesets, 
        dendrogramCellGroups = p2$clusters$PCA$multilevel, 
        show.clusters=F, 
        appmetadata = appmetadata)

    p2w$serializeToStaticFast(binary.filename = fout)    
}