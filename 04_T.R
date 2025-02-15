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

load(file = "/www/data/YXW/treateddata/obj2.rda")

obj_T = subset(obj2, predicted.celltype_l1 == "T")
table(obj_T$Label_type)

input.mat = obj_T@assays$RNA@layers$counts
row.names(input.mat) = row.names(obj2)
colnames(input.mat) = colnames(obj_T)
obj_T = CreateSeuratObject(counts = input.mat, meta.data = obj_T@meta.data, min.cells = 100)

library(SeuratWrappers)
obj_T = RunALRA(obj_T)

obj_T <- FindVariableFeatures(obj_T, verbose = F, nfeatures = 500)
obj_T = NormalizeData(obj_T)
obj_T <- FindVariableFeatures(obj_T, verbose = F, nfeatures = 3000)
obj_T <- ScaleData(obj_T)
obj_T <- FindVariableFeatures(obj_T, verbose = F, nfeatures = 3000)
obj_T <- RunPCA(obj_T, verbose = FALSE)
obj_T <- obj_T %>% RunUMAP(dims = 1:20) %>% FindNeighbors(dims = 1:20)
obj_T <- obj_T %>% RunTSNE(dims = 1:10)
obj_T <- obj_T %>% FindClusters(resolution = 0.5)
obj_T <- obj_T %>% RunUMAP(dims = 1:20)

obj_T$alra_snn_res.0.5 = obj_T$RNA_snn_res.0.5
obj_T$predicted.celltype_l2[!obj_T$predicted.celltype_l2 %like% "T_"] = "other"

DimPlot(obj_T, reduction = "umap", label = F, pt.size = .5, cols = my36colors, group.by = c("alra_snn_res.0.5", "Label_type", "predicted.celltype_l2"), raster = F) & theme(aspect.ratio = 1)
DimPlot(obj_T, reduction = "tsne", label = TRUE, pt.size = .5, cols = my36colors, group.by = c("alra_snn_res.0.5", "Label_type"), raster = F) & theme(aspect.ratio = 1)

dittoBarPlot(obj_T, "Label_type", group.by = "alra_snn_res.0.5", scale = "percent") + scale_fill_manual(values = my36colors)

difftable = FindAllMarkers(obj_T, only.pos = T, logfc.threshold = 0.5, return.thresh = 0.05)
top10 <- difftable %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

GetSignature = function(obj_random, gmtfile = "/www/data/signature/APC.gmt", usecores = 16){
  t1 = Sys.time()
  library(GSEABase)
  library(AUCell)
  countexp2 = obj_random@assays$alra@data
  cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores = usecores, plotStats = F, verbose = F)
  geneSets <- getGmt(gmtfile)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, verbose = F)
  signature.matrix <- t(data.frame(getAUC(cells_AUC)))
  print(dim(signature.matrix))
  print(Sys.time() - t1)
  return(signature.matrix)
}

obj_random = obj_T
signature.matrix.metab = GetSignature(obj_random, gmtfile = "/www/data/signature/tcell_basic_signature.gmt")
row.names(signature.matrix.metab) = colnames(obj_random)
obj_random2 = CreateSeuratObject(t(signature.matrix.metab), meta.data = obj_random@meta.data)

obj_random@meta.data = cbind(obj_random@meta.data, signature.matrix.metab)

FeaturePlot(obj_random, reduction = "umap", features = colnames(signature.matrix.metab), pt.size = .5, ncol = 4, cols = c("grey95", "red")) & theme(aspect.ratio = 1)