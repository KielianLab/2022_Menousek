## Made by Christopher M. Horn, MS
## Kielian Lab data
## Created: 2022-01-31
## Updated: 2023-03-28

## Notes: Cells were infected with dsRed expressing S. aureus, sacrificed, and sorted into dsRed positive/negative immune cells prior to scRNA-seq;
##        p == dsRed positive, m == dsRed negative




# Setting up environment -----
## Load in packages
packages <- c(
  'tidyverse',
  'Seurat',
  'patchwork',
  'SingleR',
  'SingleCellExperiment',
  'sctransform',
  'MAST',
  'plotly',
  'ggsci',
  'readxl',
  'fgsea',
  'data.table',
  'clusterExperiment',
  'pheatmap',
  'mgcv',
  'slingshot',
  'msigdbr',
  'celldex',
  'dittoSeq',
  'here',
  'reshape2'
)

invisible(lapply(packages, library, character.only = T))
rm(packages)

## Set options
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(12345)



# Initializing objects -----
## Setting up Seurat objects
## Load in the data
brain_p.data <- Read10X(data.dir = here('raw data', 'brain_p'))
brain_m.data <- Read10X(data.dir = here('raw data', 'brain_m'))
galea_p.data <- Read10X(data.dir = here('raw data', 'galea_p'))
galea_m.data <- Read10X(data.dir = here('raw data', 'galea_m'))

## Initialize Seurat object w/raw data
brain.p <- CreateSeuratObject(counts = brain_p.data, project = 'brain.p', min.cells = 3, min.features = 200)
brain.m <- CreateSeuratObject(counts = brain_m.data, project = 'brain.m', min.cells = 3, min.features = 200)
galea.p <- CreateSeuratObject(counts = galea_p.data, project = 'galea.p', min.cells = 3, min.features = 200)
galea.m <- CreateSeuratObject(counts = galea_m.data, project = 'galea.m', min.cells = 3, min.features = 200)

## Pre-processing and QC
brain.p[['percent.mt']] <- PercentageFeatureSet(brain.p, pattern = '^mt-')
brain.p <- subset(brain.p, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
brain.m[['percent.mt']] <- PercentageFeatureSet(brain.m, pattern = '^mt-')
brain.m <- subset(brain.m, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
galea.p[['percent.mt']] <- PercentageFeatureSet(galea.p, pattern = '^mt-')
galea.p <- subset(galea.p, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
galea.m[['percent.mt']] <- PercentageFeatureSet(galea.m, pattern = '^mt-')
galea.m <- subset(galea.m, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Create timing metadata
sample.id_1 <- rep('brain.p', length(brain.p@meta.data$orig.ident))
sample.id_2 <- rep('brain.m', length(brain.m@meta.data$orig.ident))
sample.id_3 <- rep('galea.p', length(galea.p@meta.data$orig.ident))
sample.id_4 <- rep('galea.m', length(galea.m@meta.data$orig.ident))
names(sample.id_1) <- rownames(brain.p@meta.data)
names(sample.id_2) <- rownames(brain.m@meta.data)
names(sample.id_3) <- rownames(galea.p@meta.data)
names(sample.id_4) <- rownames(galea.m@meta.data)
brain.p <- AddMetaData(brain.p, sample.id_1, col.name = 'Sample_origin')
brain.m <- AddMetaData(brain.m, sample.id_2, col.name = 'Sample_origin')
galea.p <- AddMetaData(galea.p, sample.id_3, col.name = 'Sample_origin')
galea.m <- AddMetaData(galea.m, sample.id_4, col.name = 'Sample_origin')

## Write individual object metadata to file
write.csv(brain.p@meta.data, here('output', 'brain.p_metadata.csv'))
write.csv(brain.m@meta.data, here('output', 'brain.m_metadata.csv'))
write.csv(galea.p@meta.data, here('output', 'galea.p_metadata.csv'))
write.csv(galea.m@meta.data, here('output', 'galea.m_metadata.csv'))

## Write QC metrics
brain_p.QC_mets <- dim(brain.p)
brain_m.QC_mets <- dim(brain.m)
galea_p.QC_mets <- dim(galea.p)
galea_m.QC_mets <- dim(galea.m)

write.csv(brain_p.QC_mets, here('output', 'QC', 'brain_p.QC_mets.csv'))
write.csv(brain_m.QC_mets, here('output', 'QC', 'brain_m.QC_mets.csv'))
write.csv(galea_p.QC_mets, here('output', 'QC', 'galea_p.QC_mets.csv'))
write.csv(galea_m.QC_mets, here('output', 'QC', 'galea_m.QC_mets.csv'))

brain_p.QC_mets.plot  <- VlnPlot(brain.p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
brain_m.QC_mets.plot <- VlnPlot(brain.m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
galea_p.QC_mets.plot  <- VlnPlot(galea.p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
galea_m.QC_mets.plot <- VlnPlot(galea.m, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave('brain_p.QC_mets.plot.png', plot = brain_p.QC_mets.plot, device = 'png', path = here('output', 'QC'))
ggsave('brain_m.QC_mets.plot.png', plot = brain_m.QC_mets.plot, device = 'png', path = here('output', 'QC'))
ggsave('galea_p.QC_mets.plot.png', plot = galea_p.QC_mets.plot, device = 'png', path = here('output', 'QC'))
ggsave('galea_m.QC_mets.plot.png', plot = galea_m.QC_mets.plot, device = 'png', path = here('output', 'QC'))

## Remove temp objects
rm(
  brain_p.data,
  brain_m.data,
  galea_p.data,
  galea_m.data,
  sample.id_1,
  sample.id_2,
  sample.id_3,
  sample.id_4,
  brain_p.QC_mets,
  brain_p.QC_mets.plot,
  brain_m.QC_mets,
  brain_m.QC_mets.plot,
  galea_p.QC_mets,
  galea_p.QC_mets.plot,
  galea_m.QC_mets,
  galea_m.QC_mets.plot
)

gc()



# Annotating cell types -----
## Annotate the cells w/SingleR
brain.p <- NormalizeData(brain.p)
brain.m <- NormalizeData(brain.m)
galea.p <- NormalizeData(galea.p)
galea.m <- NormalizeData(galea.m)

ref.se <- ImmGenData() # snapshot date: 2020-10-27

brain.p_sce <- as.SingleCellExperiment(brain.p)
brain.m_sce <- as.SingleCellExperiment(brain.m)
galea.p_sce <- as.SingleCellExperiment(galea.p)
galea.m_sce <- as.SingleCellExperiment(galea.m)
commonGenes.1 <- intersect(rownames(brain.p_sce), rownames(ref.se))
commonGenes.2 <- intersect(rownames(brain.m_sce), rownames(ref.se))
commonGenes.3 <- intersect(rownames(galea.p_sce), rownames(ref.se))
commonGenes.4 <- intersect(rownames(galea.m_sce), rownames(ref.se))
ref.se_1 <- ref.se[commonGenes.1,]
ref.se_2 <- ref.se[commonGenes.2,]
ref.se_3 <- ref.se[commonGenes.3,]
ref.se_4 <- ref.se[commonGenes.4,]
brain.p_sce <- brain.p_sce[commonGenes.1,]
brain.m_sce <- brain.m_sce[commonGenes.2,]
galea.p_sce <- galea.p_sce[commonGenes.3,]
galea.m_sce <- galea.m_sce[commonGenes.4,]
pred.brain_p <- SingleR(test = brain.p_sce, ref = ref.se_1, labels = ref.se_1$label.main)
pred.brain_m <- SingleR(test = brain.m_sce, ref = ref.se_2, labels = ref.se_2$label.main)
pred.galea_p <- SingleR(test = galea.p_sce, ref = ref.se_3, labels = ref.se_3$label.main)
pred.galea_m <- SingleR(test = galea.m_sce, ref = ref.se_4, labels = ref.se_4$label.main)
brain.p[['celltype']] <- pred.brain_p$pruned.labels
brain.m[['celltype']] <- pred.brain_m$pruned.labels
galea.p[['celltype']] <- pred.galea_p$pruned.labels
galea.m[['celltype']] <- pred.galea_m$pruned.labels

## Remove temp objects
rm(
  ref.se,
  brain.p_sce,
  brain.m_sce,
  galea.p_sce,
  galea.m_sce,
  commonGenes.1,
  commonGenes.2,
  commonGenes.3,
  commonGenes.4,
  ref.se_1,
  ref.se_2,
  ref.se_3,
  ref.se_4,
  pred.brain_p,
  pred.brain_m,
  pred.galea_p,
  pred.galea_m
)

gc()



# Object integration -----
## Integrating all cells from all days
brain.list <- c(brain.p, brain.m)
names(brain.list) <- c('brain.p', 'brain.m')
for (i in 1:length(brain.list))
  {brain.list[[i]] <- SCTransform(brain.list[[i]], verbose = T)}
brain.features <- SelectIntegrationFeatures(object.list = brain.list, nfeatures = 3000)
brain.list <- PrepSCTIntegration(object.list = brain.list, anchor.features = brain.features, verbose = T)
brain.anchors <- FindIntegrationAnchors(object.list = brain.list, normalization.method = 'SCT', anchor.features = brain.features, verbose = T)
brain.integrated <- IntegrateData(anchorset = brain.anchors, normalization.method = 'SCT', verbose = T)
brain.integrated <- RunPCA(brain.integrated, verbose = T)
brain.integrated <- RunUMAP(brain.integrated, dims = 1:30)
brain.integrated <- FindNeighbors(brain.integrated, dims = 1:30)
brain.integrated <- FindClusters(brain.integrated, resolution = 0.5)
saveRDS(brain.integrated, here('output', 'Joe_brain.integrated.rds'))

## Remove temp objects
rm(
  brain.list,
  brain.features,
  brain.anchors,
  brain.p,
  brain.m,
  i
)

gc()

galea.list <- c(galea.p, galea.m)
names(galea.list) <- c('galea.p', 'galea.m')
for (i in 1:length(galea.list))
  {galea.list[[i]] <- SCTransform(galea.list[[i]], verbose = T)}
galea.features <- SelectIntegrationFeatures(object.list = galea.list, nfeatures = 3000)
galea.list <- PrepSCTIntegration(object.list = galea.list, anchor.features = galea.features, verbose = T)
galea.anchors <- FindIntegrationAnchors(object.list = galea.list, normalization.method = 'SCT', anchor.features = galea.features, verbose = T)
galea.integrated <- IntegrateData(anchorset = galea.anchors, normalization.method = 'SCT', verbose = T)
galea.integrated <- RunPCA(galea.integrated, verbose = T)
galea.integrated <- RunUMAP(galea.integrated, dims = 1:30)
galea.integrated <- FindNeighbors(galea.integrated, dims = 1:30)
galea.integrated <- FindClusters(galea.integrated, resolution = 0.5)
saveRDS(galea.integrated, here('output', 'Joe_galea.integrated.rds'))

## Remove temp objects
rm(
  galea.list,
  galea.features,
  galea.anchors,
  galea.p,
  galea.m,
  i
)

gc()

## Output celltype composition of each cluster
brain_clust.comp <- brain.integrated@meta.data %>%
  select(seurat_clusters, celltype) %>%
  table() %>%
  as_tibble() %>%
  filter(n > 0) %>%
  arrange(seurat_clusters)

brain_split_clust.comp <- brain.integrated@meta.data %>%
  select(seurat_clusters, Sample_origin, celltype) %>%
  table() %>%
  as_tibble() %>%
  filter(n > 0) %>%
  arrange(seurat_clusters)

galea_clust.comp <- galea.integrated@meta.data %>%
  select(seurat_clusters, celltype) %>%
  table() %>%
  as_tibble() %>%
  filter(n > 0) %>%
  arrange(seurat_clusters)

galea_split_clust.comp <- galea.integrated@meta.data %>%
  select(seurat_clusters, Sample_origin, celltype) %>%
  table() %>%
  as_tibble() %>%
  filter(n > 0) %>%
  arrange(seurat_clusters)

write_csv(brain_clust.comp, here('output', 'brain', 'brain_clust_comp.csv'))
write_csv(brain_split_clust.comp, here('output', 'brain', 'brain_split_clust_comp.csv'))
write_csv(galea_clust.comp, here('output', 'galea', 'galea_clust_comp.csv'))
write_csv(galea_split_clust.comp, here('output', 'galea', 'galea_split_clust_comp.csv'))

rm(
  brain_clust.comp,
  brain_split_clust.comp,
  galea_clust.comp,
  galea_split_clust.comp
)

gc()



# Renaming/plotting -----
## Renaming clusters
brain_new.cluster.ids <- c(
  'Mono/Mac 1',
  'Microglia 1',
  'Microglia 2',
  'Granulocytes 1',
  'Microglia 3',
  'Granulocytes 2',
  'Mono/Mac/DC',
  'ILC/NKT/T/Tgd',
  'Mono/Mac 2',
  'ILC/NK',
  'B Cells',
  'Granulocytes 3',
  'Monocytes',
  'Basophils'
)

names(brain_new.cluster.ids) <- levels(brain.integrated)
brain.integrated <- RenameIdents(brain.integrated, brain_new.cluster.ids)

galea_new.cluster.ids <- c(
  'Granulocytes 1',
  'Granulocytes 2',
  'Granulocytes 3',
  'Granulocytes 4',
  'Granulocytes 5',
  'Granulocytes 6',
  'Granulocytes 7',
  'C7',
  'Granulocytes 8'
)

names(galea_new.cluster.ids) <- levels(galea.integrated)
galea.integrated <- RenameIdents(galea.integrated, galea_new.cluster.ids)

rm(
  brain_new.cluster.ids,
  galea_new.cluster.ids
)

gc()

## Save seurat objects
saveRDS(brain.integrated, here('output', 'Joe_brain.integrated.rds'))
saveRDS(galea.integrated, here('output', 'Joe_galea.integrated.rds'))

## Plotting w/new labels
brain_pca.origin <- DimPlot(brain.integrated, group.by = 'Sample_origin', reduction = 'pca', cols = c('slateblue1', 'tomato')) +
  theme_classic() +
  ggtitle('Split by Sample origin') +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

brain_pca.clus <- DimPlot(brain.integrated, reduction = 'pca', label = T, label.box = T, repel = T) +
  theme_classic() +
  NoLegend() +
  labs(x = 'PC 1', y = 'PC 2') + theme(axis.text = element_blank(), axis.ticks = element_blank())

brain_umap.origin <- DimPlot(brain.integrated, group.by = 'Sample_origin', cols = c('slateblue1', 'tomato')) +
  theme_classic() +
  ggtitle('Split by Sample origin') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

brain_umap.clus <- DimPlot(brain.integrated, label = T, label.box = T, repel = T) +
  theme_classic() +
  NoLegend() +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

galea_pca.origin <- DimPlot(galea.integrated, group.by = 'Sample_origin', reduction = 'pca', cols = c('slateblue1', 'tomato')) +
  theme_classic() +
  ggtitle('Split by Sample origin') +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

galea_pca.clus <- DimPlot(galea.integrated, reduction = 'pca', label = T, label.box = T, repel = T) +
  theme_classic() +
  NoLegend() +
  labs(x = 'PC 1', y = 'PC 2') + theme(axis.text = element_blank(), axis.ticks = element_blank())

galea_umap.origin <- DimPlot(galea.integrated, group.by = 'Sample_origin', cols = c('slateblue1', 'tomato')) +
  theme_classic() +
  ggtitle('Split by Sample origin') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

galea_umap.clus <- DimPlot(galea.integrated, label = T, label.box = T, repel = T) +
  theme_classic() +
  NoLegend() +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave('brain_pca.origin.png', plot = brain_pca.origin, device = 'png', path = here('output', 'brain'))
ggsave('brain_pca.cluster.png', plot = brain_pca.clus, device = 'png', path = here('output', 'brain'))
ggsave('brain_umap.origin.png', plot = brain_umap.origin, device = 'png', path = here('output', 'brain'))
ggsave('brain_umap.cluster.png', plot = brain_umap.clus, device = 'png', path = here('output', 'brain'))
ggsave('galea_pca.origin.png', plot = galea_pca.origin, device = 'png', path = here('output', 'galea'))
ggsave('galea_pca.cluster.png', plot = galea_pca.clus, device = 'png', path = here('output', 'galea'))
ggsave('galea_umap.origin.png', plot = galea_umap.origin, device = 'png', path = here('output', 'galea'))
ggsave('galea_umap.cluster.png', plot = galea_umap.clus, device = 'png', path = here('output', 'galea'))

rm(
  brain_pca.origin,
  brain_pca.clus,
  brain_umap.origin,
  brain_umap.clus,
  galea_pca.origin,
  galea_pca.clus,
  galea_umap.origin,
  galea_umap.clus
)

gc()



## Output cluster number to name cheat sheet
brain.num2name <- data.frame('Cluster_name' = levels(brain.integrated), 'Cluster_num' = levels(brain.integrated@meta.data$seurat_clusters))
galea.num2name <- data.frame('Cluster_name' = levels(galea.integrated), 'Cluster_num' = levels(galea.integrated@meta.data$seurat_clusters))

write_csv(brain.num2name, here('output', 'brain', 'brain_num 2 name cheat sheet.csv'))
write_csv(galea.num2name, here('output', 'galea', 'galea_num 2 name cheat sheet.csv'))

## Write the metadata to a file
write.csv(brain.integrated@meta.data, here('output', 'brain', 'brain.int_metadata.csv'))
write.csv(galea.integrated@meta.data, here('output', 'galea', 'galea.int_metadata.csv'))

rm(
  brain.num2name,
  galea.num2name
)

gc()

## Output m vs p for each cluster
brain.mp <- table(Idents(brain.integrated), brain.integrated$Sample_origin)
galea.mp <- table(Idents(galea.integrated), galea.integrated$Sample_origin)

write.csv(brain.mp, here('output', 'brain_mp.csv'))
write.csv(galea.mp, here('output', 'galea_mp.csv'))

rm(
  brain.mp,
  galea.mp
)

gc()

## Subsetting out
brain.microglia <- subset(brain.integrated, idents = c(
  'Microglia 1',
  'Microglia 2',
  'Microglia 3'
))

brain.granulocyte <- subset(brain.integrated, idents = c(
  'Granulocytes 1',
  'Granulocytes 2',
  'Granulocytes 3'
))

brain.monomac <- subset(brain.integrated, idents = c(
  'Mono/Mac 1',
  'Mono/Mac 2'
))




# Cluster-level DE -----
## Gene DE Analysis
DefaultAssay(brain.integrated) <- 'RNA'
brain.integrated <- NormalizeData(brain.integrated, verbose = T)
brain.integrated <- ScaleData(brain.integrated, verbose = T)

brain.markers <- FindAllMarkers(brain.integrated, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
write.csv(brain.markers, here('output', 'brain', 'DE', 'between cluster', 'full DE_brain.csv'))

for (i in seq_along(levels(brain.markers$cluster))) {
  if(i != length(brain.markers$cluster)) {
    DE <- brain.markers %>% filter(brain.markers$cluster %in% c(levels(brain.markers$cluster)[i]))
    write.csv(DE, here('output', 'brain', 'DE', 'between cluster', paste0('DE_brain_', i-1, '.csv')))
  }
}

brain.integrated$cell_origin <- paste(Idents(brain.integrated), brain.integrated$Sample_origin, sep = '_')
brain.integrated$cell <- Idents(brain.integrated)
Idents(brain.integrated) <- 'cell_origin'

DefaultAssay(galea.integrated) <- 'RNA'
galea.integrated <- NormalizeData(galea.integrated, verbose = T)
galea.integrated <- ScaleData(galea.integrated, verbose = T)

galea.markers <- FindAllMarkers(galea.integrated, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
write.csv(galea.markers, here('output', 'galea', 'DE', 'between cluster', 'full DE_galea.csv'))

for (i in seq_along(levels(galea.markers$cluster))) {
  if(i != length(galea.markers$cluster)) {
    DE <- galea.markers %>% filter(galea.markers$cluster %in% c(levels(galea.markers$cluster)[i]))
    write.csv(DE, here('output', 'galea', 'DE', 'between cluster', paste0('DE_galea_', i-1, '.csv')))
  }
}

galea.integrated$cell_origin <- paste(Idents(galea.integrated), galea.integrated$Sample_origin, sep = '_')
galea.integrated$cell <- Idents(galea.integrated)
Idents(galea.integrated) <- 'cell_origin'

# brain.m v. brain.p
brain.de_00 <- FindMarkers(brain.integrated, ident.2 = 'Mono/Mac 1_brain.m', ident.1 = 'Mono/Mac 1_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_01 <- FindMarkers(brain.integrated, ident.2 = 'Microglia 1_brain.m', ident.1 = 'Microglia 1_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_02 <- FindMarkers(brain.integrated, ident.2 = 'Microglia 2_brain.m', ident.1 = 'Microglia 2_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_03 <- FindMarkers(brain.integrated, ident.2 = 'Granulocytes 1_brain.m', ident.1 = 'Granulocytes 1_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_04 <- FindMarkers(brain.integrated, ident.2 = 'Microglia 3_brain.m', ident.1 = 'Microglia 3_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_05 <- FindMarkers(brain.integrated, ident.2 = 'Granulocytes 2_brain.m', ident.1 = 'Granulocytes 2_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_06 <- FindMarkers(brain.integrated, ident.2 = 'Mono/Mac/DC_brain.m', ident.1 = 'Mono/Mac/DC_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_07 <- FindMarkers(brain.integrated, ident.2 = 'ILC/NKT/T/Tgd_brain.m', ident.1 = 'ILC/NKT/T/Tgd_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_08 <- FindMarkers(brain.integrated, ident.2 = 'Mono/Mac 2_brain.m', ident.1 = 'Mono/Mac 2_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_09 <- FindMarkers(brain.integrated, ident.2 = 'ILC/NK_brain.m', ident.1 = 'ILC/NK_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_10 <- FindMarkers(brain.integrated, ident.2 = 'B Cells_brain.m', ident.1 = 'B Cells_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_11 <- FindMarkers(brain.integrated, ident.2 = 'Granulocytes 3_brain.m', ident.1 = 'Granulocytes 3_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_12 <- FindMarkers(brain.integrated, ident.2 = 'Monocytes_brain.m', ident.1 = 'Monocytes_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
brain.de_13 <- FindMarkers(brain.integrated, ident.2 = 'Basophils_brain.m', ident.1 = 'Basophils_brain.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)

brain_de.list <- list(
  brain.de_00,
  brain.de_01,
  brain.de_02,
  brain.de_03,
  brain.de_04,
  brain.de_05,
  brain.de_06,
  brain.de_07,
  brain.de_08,
  brain.de_10,
  brain.de_11,
  brain.de_13
)

rm(
  brain.de_00,
  brain.de_01,
  brain.de_02,
  brain.de_03,
  brain.de_04,
  brain.de_05,
  brain.de_06,
  brain.de_07,
  brain.de_08,
  brain.de_10,
  brain.de_11,
  brain.de_13
)

gc()

# galea.m v. galea.p
galea.de_00 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 1_galea.m', ident.1 = 'Granulocytes 1_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_01 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 2_galea.m', ident.1 = 'Granulocytes 2_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_02 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 3_galea.m', ident.1 = 'Granulocytes 3_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_03 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 4_galea.m', ident.1 = 'Granulocytes 4_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_04 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 5_galea.m', ident.1 = 'Granulocytes 5_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_05 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 6_galea.m', ident.1 = 'Granulocytes 6_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_06 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 7_galea.m', ident.1 = 'Granulocytes 7_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_07 <- FindMarkers(galea.integrated, ident.2 = 'C7_galea.m', ident.1 = 'C7_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)
galea.de_08 <- FindMarkers(galea.integrated, ident.2 = 'Granulocytes 8_galea.m', ident.1 = 'Granulocytes 8_galea.p', test.use = 'MAST', min.pct = 0, logfc.threshold = 0)

galea_de.list <- list(
  galea.de_00,
  galea.de_01,
  galea.de_02,
  galea.de_03,
  galea.de_04,
  galea.de_05,
  galea.de_06,
  galea.de_07,
  galea.de_08
)

rm(
  galea.de_00,
  galea.de_01,
  galea.de_02,
  galea.de_03,
  galea.de_04,
  galea.de_05,
  galea.de_06,
  galea.de_07,
  galea.de_08
)

gc()


for (i in 1:length(brain_de.list)) {
  brain_de.list[[i]] %>%
    arrange(desc(avg_log2FC)) %>% 
    write.csv(here('output', 'brain', 'DE', paste0('brain.de_', i-1, '.csv')))
}

for (i in 1:length(galea_de.list)) {
  galea_de.list[[i]] %>%
    arrange(desc(avg_log2FC)) %>% 
    write.csv(here('output', 'galea', 'DE', paste0('galea.de_', i-1, '.csv')))
}


rm(i)

gc()



## fgsea
GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

geneSet_list <- list(GO.set, HM.set, CP.set)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in seq_along(levels(brain.markers$cluster))) {
      glist <- brain.markers %>% filter(brain.markers$cluster %in% c(levels(brain.markers$cluster)[ii]))
      stats <- glist$avg_log2FC
      names(stats) <- toupper(glist$gene)
      stats <- sort(stats, decreasing = T)
      eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
      eaRes <- arrange(eaRes, desc(NES))
      fwrite(eaRes, file = here('output', 'brain', 'GSEA', 'between cluster', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_brain eaRes.', ii-1, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
}

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in seq_along(levels(galea.markers$cluster))) {
    glist <- galea.markers %>% filter(galea.markers$cluster %in% c(levels(galea.markers$cluster)[ii]))
    stats <- glist$avg_log2FC
    names(stats) <- toupper(glist$gene)
    stats <- sort(stats, decreasing = T)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(eaRes, file = here('output', 'galea', 'GSEA', 'between cluster', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_galea eaRes.', ii-1, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
}

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:length(brain_de.list)) {
    glist <- brain_de.list[[ii]]
    stats <- glist$avg_log2FC
    names(stats) <- toupper(rownames(glist))
    stats <- sort(stats, decreasing = T)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- filter(eaRes, padj <= 0.25)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(eaRes, file = here('output', 'brain', 'GSEA', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes.', ii-1, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
}

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:length(galea_de.list)) {
    glist <- galea_de.list[[ii]]
    stats <- glist$avg_log2FC
    names(stats) <- toupper(rownames(glist))
    stats <- sort(stats, decreasing = T)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- filter(eaRes, padj <= 0.25)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(eaRes, file = here('output', 'galea', 'GSEA', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes.', ii-1, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
}


rm(
  geneSet_list,
  m_list, 
  stats,
  glist, 
  eaRes,
  brain_de.list,
  galea_de.list,
  DE,
  brain.markers,
  galea.markers,
  GO.set,
  HM.set,
  CP.set,
  i,
  ii
)

gc()



DefaultAssay(brain.microglia) <- 'RNA'
brain.microglia <- NormalizeData(brain.microglia, verbose = T)
brain.microglia <- ScaleData(brain.microglia, verbose = T)

DefaultAssay(brain.granulocyte) <- 'RNA'
brain.granulocyte <- NormalizeData(brain.granulocyte, verbose = T)
brain.granulocyte <- ScaleData(brain.granulocyte, verbose = T)

DefaultAssay(brain.monomac) <- 'RNA'
brain.monomac <- NormalizeData(brain.monomac, verbose = T)
brain.monomac <- ScaleData(brain.monomac, verbose = T)

brain_microglia.markers <- FindAllMarkers(brain.microglia, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
brain_granulocyte.markers <- FindAllMarkers(brain.granulocyte, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
brain_monomac.markers <- FindAllMarkers(brain.monomac, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')

brain_microglia.filtered <- filter(brain_microglia.markers, brain_microglia.markers$p_val_adj <= 0.05)
write.csv(brain_microglia.filtered, here('output', 'brain', 'DE', 'microglia', 'microglia_full DE.csv'))

brain_granulocyte.filtered <- filter(brain_granulocyte.markers, brain_granulocyte.markers$p_val_adj <= 0.05)
write.csv(brain_granulocyte.filtered, here('output', 'brain', 'DE', 'granulocyte', 'granulocyte_full DE.csv'))

brain_monomac.filtered <- filter(brain_monomac.markers, brain_monomac.markers$p_val_adj <= 0.05)
write.csv(brain_monomac.filtered, here('output', 'brain', 'DE', 'monomac', 'monomac_full DE.csv'))

## fgsea

GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

geneSet_list <- list(GO.set, HM.set, CP.set)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in seq_along(levels(brain_microglia.markers$cluster))) {
    if(ii != length(brain_microglia.markers$cluster)) {
    glist <- brain_microglia.markers %>% filter(brain_microglia.markers$cluster %in% c(levels(brain_microglia.markers$cluster)[ii]))
    stats <- glist$avg_log2FC
    names(stats) <- toupper(glist$gene)
    stats <- sort(stats, decreasing = T)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- filter(eaRes, padj <= 0.25)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(eaRes, file = here('output', 'brain', 'GSEA', 'microglia', paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
  }
}

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in seq_along(levels(brain_granulocyte.markers$cluster))) {
    if(ii != length(brain_granulocyte.markers$cluster)) {
      glist <- brain_granulocyte.markers %>% filter(brain_granulocyte.markers$cluster %in% c(levels(brain_granulocyte.markers$cluster)[ii]))
      stats <- glist$avg_log2FC
      names(stats) <- toupper(glist$gene)
      stats <- sort(stats, decreasing = T)
      eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
      eaRes <- filter(eaRes, padj <= 0.25)
      eaRes <- arrange(eaRes, desc(NES))
      fwrite(eaRes, file = here('output', 'brain', 'GSEA', 'granulocyte', paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii, '.tsv')), sep="\t", sep2=c("", " ", ""))
    }
  }
}

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in seq_along(levels(brain_monomac.markers$cluster))) {
    if(ii != length(brain_monomac.markers$cluster)) {
      glist <- brain_monomac.markers %>% filter(brain_monomac.markers$cluster %in% c(levels(brain_monomac.markers$cluster)[ii]))
      stats <- glist$avg_log2FC
      names(stats) <- toupper(glist$gene)
      stats <- sort(stats, decreasing = T)
      eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
      eaRes <- filter(eaRes, padj <= 0.25)
      eaRes <- arrange(eaRes, desc(NES))
      fwrite(eaRes, file = here('output', 'brain', 'GSEA', 'monomac', paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii, '.tsv')), sep="\t", sep2=c("", " ", ""))
    }
  }
}

rm(
  brain.microglia,
  brain.granulocyte,
  brain.monomac,
  brain_microglia.markers,
  brain_granulocyte.markers,
  brain_monomac.markers,
  brain_microglia.filtered,
  brain_granulocyte.filtered,
  brain_monomac.filtered,
  GO.set,
  CP.set,
  HM.set,
  geneSet_list,
  m_list,
  glist,
  stats,
  eaRes,
  i,
  ii
)

gc()


# Brain & galea integration ------
cranio.list <- c(brain.integrated, galea.integrated)
names(cranio.list) <- c('brain', 'galea')
for (i in 1:length(cranio.list))
{cranio.list[[i]] <- SCTransform(cranio.list[[i]], verbose = T)}
cranio.features <- SelectIntegrationFeatures(object.list = cranio.list, nfeatures = 3000)
cranio.list <- PrepSCTIntegration(object.list = cranio.list, anchor.features = cranio.features, verbose = T)
cranio.anchors <- FindIntegrationAnchors(object.list = cranio.list, normalization.method = 'SCT', anchor.features = cranio.features, verbose = T)
cranio.integrated <- IntegrateData(anchorset = cranio.anchors, normalization.method = 'SCT', verbose = T)
cranio.integrated <- RunPCA(cranio.integrated, verbose = T)
cranio.integrated <- RunUMAP(cranio.integrated, dims = 1:30)
cranio.integrated <- FindNeighbors(cranio.integrated, dims = 1:30)
cranio.integrated <- FindClusters(cranio.integrated, resolution = 0.5)
saveRDS(cranio.integrated, here('output', 'Joe_cranio.integrated.rds'))

## Remove temp objects
rm(
  cranio.list,
  cranio.features,
  cranio.anchors,
  brain.integrated,
  galea.integrated,
  i
)

gc()

clust.comp <- prop.table(table(Idents(cranio.integrated), cranio.integrated$celltype), margin = 1)

write.csv(clust.comp, here('output', 'agg', 'clust_comp.csv'))

rm(clust.comp)

gc()

cranio_new.cluster.ids <- c(
  'Granulocytes 1',
  'Granulocytes 2',
  'Granulocytes/Microglia 1',
  'Granulocytes 3',
  'Microglia',
  'Mono/Mac',
  'Granulocytes/Monocytes',
  'Granulocytes/Microglia 2',
  'T/NK Cells',
  'Granulocytes/Microglia 3',
  'Mono/DC',
  'Granulocytes/Macrophages',
  'B Cells'
)

names(cranio_new.cluster.ids) <- levels(cranio.integrated)
cranio.integrated <- RenameIdents(cranio.integrated, cranio_new.cluster.ids)

rm(cranio_new.cluster.ids)

gc()

## Save seurat objects
saveRDS(cranio.integrated, here('output', 'Joe_cranio.integrated.rds'))

## Plotting w/new labels
cranio_pca.origin <- DimPlot(cranio.integrated, group.by = 'Sample_origin', reduction = 'pca', cols = c('#BC3C2999', '#0072B599', '#E1872799', '#20854E99'), pt.size = 2) +
  theme_classic() +
  ggtitle('Split by sample origin') +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

cranio_pca.clus <- DimPlot(cranio.integrated, reduction = 'pca', label = T, label.box = T, repel = T, pt.size = 2) +
  theme_classic() +
  NoLegend() +
  labs(x = 'PC 1', y = 'PC 2') + theme(axis.text = element_blank(), axis.ticks = element_blank())

cranio_umap.origin <- DimPlot(cranio.integrated, group.by = 'Sample_origin', cols = c('#BC3C2999', '#0072B599', '#E1872799', '#20854E99'), pt.size = 2) +
  theme_classic() +
  ggtitle('Split by sample origin') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

cranio_umap.origin2 <- DimPlot(cranio.integrated, group.by = 'Sample_origin', cols = c('tomato', 'slateblue1', 'tomato', 'slateblue1'), pt.size = 2) +
  theme_classic() +
  ggtitle('Split by sample origin') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

cranio_umap.clus <- DimPlot(cranio.integrated, label = T, label.box = T, repel = T, pt.size = 2) +
  theme_classic() +
  NoLegend() +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave('cranio_pca.origin.png', plot = cranio_pca.origin, device = 'png', path = here('output', 'agg'))
ggsave('cranio_pca.cluster.png', plot = cranio_pca.clus, device = 'png', path = here('output', 'agg'))
ggsave('cranio_umap.origin.png', plot = cranio_umap.origin, device = 'png', path = here('output', 'agg'))
ggsave('cranio_umap.origin2.png', plot = cranio_umap.origin2, device = 'png', path = here('output', 'agg'))
ggsave('cranio_umap.cluster.png', plot = cranio_umap.clus, device = 'png', path = here('output', 'agg'))

rm(
  cranio_pca.origin,
  cranio_pca.clus,
  cranio_umap.origin,
  cranio_umap.origin2,
  cranio_umap.clus
)

gc()

## Output cluster number to name cheat sheet
cranio.num2name <- data.frame('Cluster_name' = levels(cranio.integrated), 'Cluster_num' = levels(cranio.integrated@meta.data$seurat_clusters))

write_csv(cranio.num2name, here('output', 'agg', 'cranio_num 2 name cheat sheet.csv'))

## Write the metadata to a file
write.csv(cranio.integrated@meta.data, here('output', 'agg', 'cranio.int_metadata.csv'))

rm(cranio.num2name)

gc()

## Output m vs p for each cluster
cranio.origin <- table(Idents(cranio.integrated), cranio.integrated$Sample_origin)

write.csv(cranio.origin, here('output', 'cranio.origin.csv'))

rm(cranio.origin)

gc()



# Complex Heatmap ------
## Average single cell data & pull out normalized counts
## Counts are normalized by dividing the counts for a given feature by the total counts per cell, multiplying by a scale factor (default == 10,000), and then taking the natural log using log1p()
brain.avg <- AverageExpression(brain.integrated, return.seurat = T)
avg.norm_counts <- brain.avg@assays$RNA@data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')

galea.avg <- AverageExpression(galea.integrated, return.seurat = T)
avg.norm_counts <- galea.avg@assays$RNA@data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')

brain_gran.avg <- AverageExpression(brain.granulocyte, return.seurat = T)
avg.norm_counts <- brain_gran.avg@assays$RNA@data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')

brain_micro.avg <- AverageExpression(brain.microglia, return.seurat = T)
avg.norm_counts <- brain_micro.avg@assays$RNA@data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')

brain_mm.avg <- AverageExpression(brain.monomac, return.seurat = T)
avg.norm_counts <- brain_mm.avg@assays$RNA@data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')

## Grabbing pathway sets
# GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
# CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

## Grabbing gene lists for desired pathways
pathways <- HM.set %>%
   filter(gs_name == 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY') # Set this to be whatever pathway you'd like

# gene_list <- pathways$human_gene_symbol
# Brain
gene_list <- c('S100A9',
               'S100A8',
               'CYBB',
               'NCF4',
               'RAC2',
               'SEM1',
               'HIGD1A',
               'ELOB', 
               'HIF1A',
               'UBB',
               'VEGFA',
               'EGLN3',
               'RBX1',
               'UBE2D3',
               'EGLN1',
               'ELOC',
               'LIMD1',
               'UBA52',
               'HMOX1',
               'SEM1',
               'UBB',
               'RBX1',
               'CSNK2B',
               'CYBB',
               'PRDX6',
               'GPX1',
               'PRDX5',
               'PRDX1',
               'ATOX1',
               'NCF2',
               'NCF4',
               'P4HB',
               'PRDX2',
               'GSR',
               'CYBA',
               'SOD1',
               'SOD2')

# Galea
gene_list <- c('CYBB',
               'CYBA',
               'NOS2',
               'NCF4',
               'TCIRG1',
               'RAC2',
               'HVCN1',
               'PRDX1',
               'GPX1',
               'SOD2',
               'ATOX1',
               'TXN1',
               'TXNRD1',
               'CYCS',
               'NCF4',
               'SOD1',
               'GSR',
               'HMOX1',
               'SEM1',
               'UBB',
               'RBX1',
               'CSNK2B',
               'RBX1',
               'H13')

gene_list <- str_to_title(gene_list)

# MDSC v PMN genes
gene_list <- c('Lcn2',
               'Slpi',
               'Wfdc17',
               'Prdx5',
               'Retnlg',
               'Junb',
               'Ctsd',
               'Nfkbia',
               'Clec4e',
               'Acod1',
               'Cd274',
               'Pla2g7',
               'Fpr2',
               'Cxcr2',
               'Camp',
               'Anxa1',
               'Ptgs2',
               'Ngp',
               'Mmp8',
               'Ltf',
               'Ly6g',
               'Cd177')

# Brian granulocyte genes
gene_list <- c('Lcn2',
               'Slpi',
               'Wfdc17',
               'Prdx5',
               'Retnlg',
               'Junb',
               'Ctsd',
               'Nfkbia',
               'Clec4e',
               'Acod1',
               'Cd274',
               'Pla2g7',
               'Fpr2',
               'Cxcr2',
               'Camp',
               'Anxa1',
               'Ptgs2',
               'Ngp',
               'Mmp8',
               'Ltf',
               'Ly6g',
               'Cd177',
               'S100a9',
               'Mmp9',
               'Csf3r',
               'Eef1a1',
               'Lrg1',
               'Igfbp6',
               'Rpl10a',
               'Gadd45b',
               'Actb',
               'Fth1',
               'Cxcl2',
               'Sod2',
               'Esd',
               'Cclr2',
               'Isg15',
               'Rsad2',
               'Ifit1',
               'Isg20',
               'Trim30a',
               'Gbp2',
               'Rtp4',
               'Parp14',
               'Osal2',
               'Slfn1')

# Brain microglia
gene_list <- c('P2ry12',
               'Selplg',
               'Siglech',
               'Tmem119',
               'Tgfbr1',
               'Jun',
               'Spp1',
               'Cybb',
               'Apoe',
               'Saa3',
               'Il10',
               'Mki67',
               'Cdk1',
               'Cdk4',
               'Ccnb2',
               'Hist1h2ab',
               'Hist1h1b',
               'Cxcl2',
               'S100a9',
               'Cx3cr1',
               'Gpr34',
               'Cst3',
               'Sparc',
               'Rhob',
               'Arhgap5',
               'Hexb',
               'Malat1',
               'Stmn1',
               'Btg1',
               'Ppia',
               'Tuba1b',
               'Cd63',
               'Mt1',
               'C1qa',
               'Ckb',
               'Ly86',
               'Cd81',
               'Olfml3',
               'C1qc',
               'C1qb',
               'Trem2')

# Brain mono/mac
gene_list <- c('Spp1',
               'Ccl5',
               'Ifitm6',
               'Apoe',
               'Saa3',
               'H2-Ab1',
               'Arg1',
               'Mif',
               'Ms4a7',
               'Mafb',
               'Axl',
               'Acod1',
               'Cybb',
               'Lyz2',
               'Fn1',
               'Vim',
               'Tgfbi',
               'F13a1',
               'Ccr2',
               'S100a4',
               'Ahnak',
               'Chil3',
               'Apoc2')

## Filtering average Seurat object for genes in pathway
heatmap.data <- avg.norm_counts %>%
  filter(gene %in% gene_list)
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]

## Create vectors for column/row annotations
# Brain
celltype_anno <- c('Monocyte',
                   'Microglia',
                   'Microglia',
                   'Granulocyte',
                   'Microglia',
                   'Granulocyte',
                   'Monocyte',
                   'Monocyte',
                   'Granulocyte',
                   'Monocyte')

# Galea
celltype_anno <- c('Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'N/A',
                   'Granulocyte')

# Brain
gene_anno <- c('Cellular Response\nto Hypoxia',
               'General\nROS',
               'General\nROS',
               'GTPases Activate\nNADPH Oxidases',
               'GTPases Activate\nNADPH Oxidases',
               'Cellular Response\nto Hypoxia',
               'General\nROS',
               'Cellular Response\nto Hypoxia',
               'General\nROS',
               'Cellular Response\nto Hypoxia',
               'Hmox1\nSignature',
               'General\nROS',
               'General\nROS',
               'Cellular Response\nto Hypoxia',
               'General\nROS',
               'Cellular Response\nto Hypoxia',
               'Cellular Response\nto Hypoxia',
               'General\nROS',
               'Hmox1\nSignature',
               'General\nROS',
               'Cellular Response\nto Hypoxia',
               'Cellular Response\nto Hypoxia',
               'GTPases Activate\nNADPH Oxidases',
               'GTPases Activate\nNADPH Oxidases',
               'Hmox1\nSignature',
               'General\nROS',
               'General\nROS',
               'Cellular Response\nto Hypoxia',
               'Hmox1\nSignature',
               'Cellular Response\nto Hypoxia',
               'General\nROS',
               'General\nROS') 

# Galea
gene_anno <- c('Cytoprotection\nby Hmox1',
               'Detoxification\nof ROS',
               'Detoxification\nof ROS',
               'ROS & RNS\nProduction',
               'Cytoprotection\nby Hmox1',
               'Detoxification\nof ROS',
               'Detoxification\nof ROS',
               'Cytoprotection\nby Hmox1',
               'ROS & RNS\nProduction',
               'Detoxification\nof ROS',
               'Detoxification\nof ROS',
               'Detoxification\nof ROS',
               'Cytoprotection\nby Hmox1',
               'ROS & RNS\nProduction',
               'Detoxification\nof ROS',
               'ROS & RNS\nProduction',
               'Cytoprotection\nby Hmox1',
               'Detoxification\nof ROS',
               'Detoxification\nof ROS',
               'Cytoprotection\nby Hmox1',
               'ROS & RNS\nProduction',
               'ROS & RNS\nProduction') 

# MDSC v PMN genes
gene_anno <- c('PMN genes',
               'MDSC genes',
               'MDSC genes',
               'MDSC genes',
               'MDSC genes',
               'PMN genes',
               'MDSC genes',
               'MDSC genes',
               'PMN genes',
               'PMN genes',
               'PMN genes',
               'PMN genes',
               'MDSC genes',
               'MDSC genes',
               'MDSC genes',
               'PMN genes',
               'MDSC genes',
               'MDSC genes',
               'MDSC genes',
               'MDSC genes',
               'MDSC genes',
               'MDSC genes')

which(heatmap.data == max(heatmap.data), arr.ind = TRUE) # Use this to find the max value within the matrix so that you can set your upper bound
which(heatmap.data == min(heatmap.data), arr.ind = TRUE) # Use this to find the min value within the matrix so that you can set your lower bound

col_fun = circlize::colorRamp2(c(0, 7), c('white', 'red')) # Check to make sure that the upper bound of this range isn't cutting off the expression of some genes in the matrix
ComplexHeatmap::Heatmap(heatmap.data,
                        name = 'Average\nNormalized\nCounts',
                        col = col_fun,
                        rect_gp = grid::gpar(col = 'black', lwd = 2),
                        cluster_columns = T,
                        cluster_column_slices = F,
                        column_gap = grid::unit(5, 'mm'),
                        column_split = celltype_anno,
                        cluster_rows = T,
                        cluster_row_slices = F,
                        row_split = gene_anno,
                        row_gap = grid::unit(5, 'mm'),
                        show_parent_dend_line = F,
                        heatmap_legend_param = list(border = 'black'),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid::grid.text(sprintf('%.1f', heatmap.data[i, j]), x, y, gp = grid::gpar(fontsize = 10))
                        })

## Pathway heatmap for brain
pathway_list <- c('HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                  'HALLMARK_HYPOXIA',
                  'HALLMARK_APOPTOSIS',
                  'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
                  'HALLAMRK_INFLAMMATORY_RESPONSE',
                  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                  'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                  'HALLMARK_GLYCOLYSIS',
                  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                  'HALLMARK_PEROXISOME',
                  'HALLAMRK_IL6_JAK_STAT3_SIGNALING',
                  'REACTOME_CYTOPROTECTION_BY_HMOX1',
                  'REACTOME_REGULATION_OF_HMOX1_EXPRESSION_AND_ACTIVITY',
                  'REACTOME_RHO_GTPASES_ACTIVATE_NADPH_OXIDASES',
                  'REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES',
                  'WP_OXIDATIVE_STRESS')

## Pathway heatmap for galea
pathway_list <- c('HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                  'HALLMARK_TGF_BETA_SIGNALING',
                  'HALLMARK_INFLAMMATORY_RESPONSE',
                  'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
                  'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                  'HALLAMRK_IL6_JAK_STAT3_SIGNALING',
                  'HALLMARK_GLYCOLYSIS',
                  'HALLMARK_HYPOXIA',
                  'HALLMARK_P53_PATHWAY',
                  'REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',
                  'REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES',
                  'REACTOME_CYTOPROTECTION_BY_HMOX1',
                  'REACTOME_REGULATION_OF_HMOX1_EXPRESSION_AND_ACTIVITY',
                  'REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA',
                  'WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SUYSTEM_IN_MITOCHONDRIA')

# Withi-cluster brain
pathway_list <- c('HALLMARK_HYPOXIA',
                  'HALLAMRK_TNFA_SIGNALING_VIA_NFKB',
                  'WP_NRF2_PATHWAY',
                  'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX',
                  'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION',
                  'KEGG_LYSOSOME',
                  'REACTOME_HSF1_DEPENDENT_TRANSACTIVATION',
                  'WP_WNT_SIGNALING',
                  'WP_RAS_SIGNALING',
                  'REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS',
                  'HALLMARK_IL6_JAK_STAT3_SIGNALING',
                  'REACTOME_NEUTROPHIL_DEGRANULATION',
                  'REACTOME_COMPLEMENT_CASCADE',
                  'REACTOME_NUCLEAR_EVENTS_KINASE_AND_TRANSCRIPTION_FACTOR_ACTIVATION',
                  'REACTOME_NGF_STIMULATED_TRANSCRIPTION',
                  'REACTOME_AMINO_ACIDS_REGULATE_MTORC1',
                  'REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',
                  'REACTOME_IRON_UPTAKE_AND_TRANSPORT',
                  'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                  'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING',
                  'NABA_ECM_GLYCOPROTEINS',
                  'NABA_CORE_MATRISOME',
                  'REACTOME_INTERACTION_BETWEEN_L1_AND_ANKYRINS',
                  'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                  'REACTOME_TRANSLATION',
                  'KEGG_RIBOSOME',
                  'REACTOME_CELLULAR_RESPONSE_TO_STARVATION')

# Withi-cluster galea
pathway_list <- c('WP_MIRNAS_INVOLVEMENT_IN_THE_IMMUNE_RESPONSE_IN_SEPSIS',
                  'WP_NUCLEAR_RECEPTORS_METAPATHWAY',
                  'KEGG_LYSOSOME',
                  'REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',
                  'REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS',
                  'KEGG_OXIDATIVE_PHOSPHORYLATION',
                  'NABA_SECRETED_FACTORS',
                  'NABA_MATRISOME_ASSOCIATED',
                  'REACTOME_TRANSFERRIN_ENDOCYTOSIS_AND_RECYCLING',
                  'REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',
                  'REACTOME_AMINO_ACIDS_REGULATE_MTORC1',
                  'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                  'REACTOME_IRON_UPTAKE_AND_TRANSPORT',
                  'REACTOME_GLUTATHIONE_CONJUGATION',
                  'REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS',
                  'REACTOME_GPCR_LIGAND_BINDING',
                  'HALLMARK_ANGIOGENESIS',
                  'REACTOME_INTERLEUKIN_17_SIGNALING',
                  'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                  'HALLMARK_INFLAMMATORY_RESPONSE',
                  'HALLMARK_GLYCOLYSIS',
                  'REACTOME_SELENOAMINO_ACID_METABOLISM',
                  'PID_AVB3_INTEGRIN_PATHWAY',
                  'KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION',
                  'HALLMARK_HYPOXIA',
                  'KEGG_RIBOSOME',
                  'REACTOME_CELLULAR_RESPONSE_TO_STARVATION')

## Filtering pathway data for pathways in pathway_list
heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- dcast(heatmap.data, pathway ~ cluster_name, value.var = 'NES', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'pathway')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- dcast(heatmap2.data, pathway ~ cluster_name, value.var = 'padj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'pathway')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

## Create vectors for column/row annotations for brain; Check to make sure annotation order matches how the genes appear in the matrix
celltype_anno <- c('Lymphocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Lymphocyte',
                   'Lymphocyte',
                   'Microglia',
                   'Microglia',
                   'Microglia',
                   'Monocyte',
                   'Monocyte',
                   'Monocyte',
                   'Monocyte')

## Create vectors for column/row annotations for brain; Check to make sure annotation order matches how the genes appear in the matrix
celltype_anno <- c('N/A',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte')


which(heatmap.data == max(heatmap.data), arr.ind = TRUE) # Use this to find the max value within the matrix so that you can set your upper bound
which(heatmap.data == min(heatmap.data), arr.ind = TRUE) # Use this to find the min value within the matrix so that you can set your lower bound

col_fun = circlize::colorRamp2(c(-4, 0, 4), c('blue', 'white', 'yellow')) # Check to make sure that the upper bound of this range isn't cutting off the expression of some genes in the matrix
ht <- ComplexHeatmap::Heatmap(heatmap.data,
                        name = 'NES',
                        col = col_fun,
                        rect_gp = grid::gpar(col = 'black', lwd = 2),
                        cluster_columns = T,
                        cluster_column_slices = F,
                        column_gap = grid::unit(5, 'mm'),
                        column_names_gp = grid::gpar(fontsize = 8),
                        cluster_rows = F,
                        cluster_row_slices = F,
                        row_gap = grid::unit(5, 'mm'),
                        show_parent_dend_line = F,
                        heatmap_legend_param = list(border = 'black'),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if(heatmap2.data[i, j] < 0.0001) {
                            grid::grid.text('****', x, y)
                          } else if(heatmap2.data[i, j] < 0.001) {
                            grid::grid.text('***', x, y)
                          } else if(heatmap2.data[i, j] < 0.01) {
                            grid::grid.text('**', x, y) 
                          } else if(heatmap2.data[i, j] < 0.05) {
                            grid::grid.text('*', x, y)
                          }
                        })

ComplexHeatmap::draw(ht, padding = unit(c(40, 5, 2, 5), 'mm')) # bottom, left, top, right paddings
