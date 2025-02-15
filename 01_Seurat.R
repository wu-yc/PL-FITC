my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                         '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                         '#968175')
                         
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)

obj = NormalizeData(obj)
obj <- FindVariableFeatures(obj, verbose = F, nfeatures = 3000)
obj <- ScaleData(obj)
obj <- FindVariableFeatures(obj, verbose = F, nfeatures = 3000)
obj <- RunPCA(obj, verbose = FALSE)
obj <- obj %>% RunUMAP(dims = 1:50) %>% FindNeighbors(dims = 1:50)
obj <- obj %>% FindClusters(resolution = 0.1)

DimPlot(obj, reduction = "umap", label = TRUE, pt.size = .5, cols = my36colors, group.by = c("RNA_snn_res.0.1"), raster = F) & theme(aspect.ratio = 1)

quantile(obj$Label)
hist(log(obj$Label + 1))

VlnPlot(obj, features = c("nCount_RNA"), pt.size = 0, log = T) & geom_boxplot(outlier.shape = NA) & scale_fill_manual(values = my36colors)
VlnPlot(obj, features = c("nFeature_RNA"), pt.size = 0, log = T) & geom_boxplot(outlier.shape = NA) & scale_fill_manual(values = my36colors)

table(obj$RNA_snn_res.0.1)

library(SingleR)
load("/www/data/reference/singler_hpca.rda")
ref.se <- sce

obj_random = obj
input.count = obj_random@assays$RNA$counts
anno.cell.main = SingleR(test = input.count, ref = ref.se, labels = ref.se$label.main)
obj_random$pred_singler <- as.character(anno.cell.main$labels)
table(as.character(anno.cell.main$labels))

pred_singler.stat <- data.frame(table(obj_random$pred_singler))
pred_singler.stat <- subset(pred_singler.stat, pred_singler.stat$Freq > 20)
obj_random$pred_singler[!(obj_random$pred_singler %in% pred_singler.stat$Var1)] <- "other"

table(obj_random$pred_singler)
tmp <- data.frame(table(obj_random$pred_singler))

obj_random$Label_log = log(obj_random$Label + 1)
quantile(obj_random$Label_log, 0.01)
obj_random$Label_log[obj_random$Label_log < quantile(obj_random$Label_log, 0.01)] = quantile(obj_random$Label_log, 0.01)

DimPlot(obj_random, reduction = "tsne", label = F, pt.size = .5, group.by = c("pred_singler"), ncol = 2, cols = my36colors) +
  FeaturePlot(obj_random, reduction = "tsne", features = "Label_log", pt.size = .5, ncol = 1) & theme(aspect.ratio = 1)

Idents(obj_random) = obj_random$pred_singler
VlnPlot(obj_random, features = c("Label_log"), pt.size = 0, log = T) & geom_boxplot(outlier.shape = NA) & scale_fill_manual(values = my36colors)

DotPlot(obj_random, features = "Label_log", cols = c("lightgrey", "darkgreen"), group.by = "pred_singler")

hist(obj_random$Label_log)

library(dittoSeq)
dittoRidgePlot(obj_random, "Label_log", group.by = "pred_singler")

obj_random_T = subset(obj_random, pred_singler == "T_cells")
hist(obj_random_T$Label_log)
sum(obj_random_T$Label_log > 4) / ncol(obj_random_T)
sum(obj_random_T$Label_log > 4.5) / ncol(obj_random_T)

load(file = "/www/data/reference/ningzhangnature2022/reference10downsample50000.rda")

table(reference10$celltype_l1)
table(reference10$celltype_l2)

reference10 = subset(reference10, celltype_l1 != "Fb")
reference10 = subset(reference10, celltype_l1 != "Mast")
reference10 = subset(reference10, celltype_l1 != "Mono")
reference10 = subset(reference10, celltype_l1 != "Mu")
reference10 = subset(reference10, celltype_l1 != "Neu")

input.celltype = unique((reference10$celltype_l2))[unique((reference10$celltype_l2)) %like% "STMN"]
for (j in 1:length(input.celltype)) {
  reference10 = subset(reference10, celltype_l2 != input.celltype[j])
}

reference10 <- reference10 %>% RunUMAP(reduction = "harmony", dims = 1:50, return.model = TRUE)

DimPlot(reference10, reduction = "umap", group.by = "celltype_l1", label = TRUE, label.size = 3, repel = F)

quantile(reference10$nCount_RNA, 0.9)
quantile(reference10$nFeature_RNA)
reference10 = subset(reference10, nCount_RNA < 20000)

library(harmony)
reference10 <- NormalizeData(reference10)
reference10 <- FindVariableFeatures(reference10)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
reference10 <- CellCycleScoring(reference10, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
reference10 <- ScaleData(reference10, vars.to.regress = c("S.Score", "G2M.Score"))
reference10 <- RunPCA(reference10, verbose = FALSE)
reference10 <- reference10 %>% RunHarmony("orig.ident", plot_convergence = TRUE)
reference10 <- reference10 %>% RunUMAP(reduction = "harmony", dims = 1:20, return.model = TRUE) %>% FindNeighbors(reduction = "harmony", dims = 1:20)
DimPlot(reference10, reduction = "umap", label = TRUE, pt.size = .05, group.by = "celltype_l1", raster = F)

reference10 <- reference10 %>% FindClusters(resolution = 1)
DimPlot(reference10, reduction = "umap", label = TRUE, pt.size = .5, cols = my36colors, group.by = c("RNA_snn_res.1"), raster = F) & theme(aspect.ratio = 1)

library(scProgram)
FeatureMatrix = GetFeatures(obj = reference10, group.by = "RNA_snn_res.1", genenumber = 50, pct_exp = 0.1, mode = "fast")

obj = RunALRA(obj)
obj_tmp = obj
obj_tmp@assays$RNA@layers$data = obj_tmp@assays$alra@data

t1 = Sys.time()
anchors <- FindTransferAnchors(
  reference = reference10,
  query = obj_tmp,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30,
  k.filter = NA,
  k.anchor = 5,
  k.score = 30,
  max.features = 500,
  n.trees = 50
)
print(Sys.time() - t1)

obj2 <- MapQuery(
  anchorset = anchors,
  query = obj_tmp,
  reference = reference10,
  refdata = list(celltype_l1 = "celltype_l1", celltype_l2 = "celltype_l2"),
  reference.reduction = "pca",
  reduction.model = "umap"
)

DimPlot(reference10, reduction = "umap", group.by = "celltype_l1", label = TRUE, label.size = 3, repel = F) +
  DimPlot(obj2, reduction = "ref.umap", group.by = "predicted.celltype_l1", label = TRUE, label.size = 3, repel = F) & theme(aspect.ratio = 1)

DimPlot(obj2, reduction = "ref.umap", label = F, pt.size = .5, group.by = c("orig.ident"), ncol = 1, cols = my36colors[16:36]) & theme(aspect.ratio = 1)