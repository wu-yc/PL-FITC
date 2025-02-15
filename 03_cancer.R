my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                         '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                         '#968175')
                         
load(file = "/www/data/YXW/treateddata/obj2.rda")

obj_cancer = subset(obj2, predicted.celltype_l1 == "Tumor")
table(obj_cancer$Label_type)

input.mat = obj_cancer@assays$RNA@layers$counts
row.names(input.mat) = row.names(obj2)
colnames(input.mat) = colnames(obj_cancer)
obj_cancer = CreateSeuratObject(counts = input.mat, meta.data = obj_cancer@meta.data, min.cells = 50)

obj_cancer <- FindVariableFeatures(obj_cancer, verbose = F, nfeatures = 500)
obj_cancer = NormalizeData(obj_cancer)
obj_cancer <- FindVariableFeatures(obj_cancer, verbose = F, nfeatures = 3000)
obj_cancer <- ScaleData(obj_cancer)
obj_cancer <- FindVariableFeatures(obj_cancer, verbose = F, nfeatures = 3000)
obj_cancer <- RunPCA(obj_cancer, verbose = FALSE)
obj_cancer <- obj_cancer %>% RunUMAP(dims = 1:20) %>% FindNeighbors(dims = 1:20)
obj_cancer <- obj_cancer %>% RunTSNE(dims = 1:10)
obj_cancer <- obj_cancer %>% FindClusters(resolution = 0.5)
obj_cancer <- obj_cancer %>% RunUMAP(dims = 1:20)

obj_cancer$alra_snn_res.0.5 = obj_cancer$RNA_snn_res.0.5

DimPlot(obj_cancer, reduction = "umap", label = F, pt.size = .5, cols = my36colors, group.by = c("alra_snn_res.0.5", "Label_type"), raster = F) & theme(aspect.ratio = 1)
DimPlot(obj_cancer, reduction = "tsne", label = TRUE, pt.size = .5, cols = my36colors, group.by = c("alra_snn_res.0.5", "Label_type"), raster = F) & theme(aspect.ratio = 1)

dittoBarPlot(obj_cancer, "Label_type", group.by = "alra_snn_res.0.5", scale = "percent") + scale_fill_manual(values = my36colors)

library(Seurat)
library(SeuratDisk)
system("rm -f /www/data/YXW/treateddata/sctour/PanB.h5Seurat")
system("rm -f /www/data/YXW/treateddata/sctour/PanB.h5ad")

SaveH5Seurat(obj_cancer, filename = "/www/data/YXW/treateddata/sctour/PanB.h5Seurat")
Convert("/www/data/YXW/treateddata/sctour/PanB.h5Seurat", dest = "h5ad")

mat_count = obj_cancer@assays$alra@data
dim(mat_count)
quantile(rowMeans(mat_count))

mat_count = subset(as.matrix(mat_count), rowMeans(mat_count) > 0.5)
quantile(colMeans(mat_count))
quantile((mat_count))

write.table(mat_count, "/www/data/YXW/treateddata/sctour/PanB_cancer_X_vargene.txt", col.names = T, row.names = T, sep="\t", quote = F)

meta = obj_cancer@meta.data[,c("Label_type", "alra_snn_res.0.1")]
write.csv(meta, "/www/data/YXW/treateddata/sctour/PanB_cancer_meta.csv", row.names = T, quote = F)

library(scProgram)
FeatureMatrix = GetFeatures(obj = obj_cancer, assay = "alra", group.by = "alra_snn_res.0.5", genenumber = 100, pct_exp = 0.1, mode = "fast")
difftable = FindAllMarkers(obj_cancer, only.pos = T, logfc.threshold = 0.5, return.thresh = 0.05, slot = "data", assay = "alra")
top10 <- difftable %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

avgexp = AggregateExpression(obj_cancer, group.by = "alra_snn_res.0.5", assay = "alra", slot = "data")$alra
input.features2 = unique(top10$gene)
avgexp = avgexp[input.features2,]
pheatmap(avgexp, color = colorRampPalette(c("white", "white", "white", "#52A85F"))(100), show_rownames = F, cluster_rows = F, cluster_cols = F, scale = "row")

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

obj_random = obj_cancer
signature.matrix.metab = GetSignature(obj_random, gmtfile = "/www/data/signature/h.all.v7.2.symbols.c2.cp.kegg.v7.2.symbols.gmt")
row.names(signature.matrix.metab) = colnames(obj_random)
obj_random2 = CreateSeuratObject(t(signature.matrix.metab), meta.data = obj_random@meta.data)
Idents(obj_random2) = obj_random2$Label_type

obj_random2 = NormalizeData(obj_random2)
diff.metab.celltype <- FindMarkers(obj_random2, ident.1 = "Positive", ident.2 = "Negative", logfc.threshold = -1, return.thresh = 1, assay = "RNA", slot = "counts")
diff.metab.celltype$gene = row.names(diff.metab.celltype)

VlnPlot(obj_random2, features = c(row.names(obj_random2)), pt.size = 0, log = F, ncol = 4) & geom_boxplot(outlier.shape = NA) & scale_fill_manual(values = my36colors) & theme(aspect.ratio = 1)

library(dittoSeq)
dittoHeatmap(obj_random2, row.names(obj_random2), assay = "RNA", slot = "counts", scale = "row", scaled.to.max = F, show_colnames = FALSE, show_rownames = T, order.by = "Label_type", annot.by = c("Label_type"))