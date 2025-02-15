load(file = "/www/data/YXW/treateddata/obj2.rda")

obj_macro = subset(obj2, predicted.celltype_l1 == "Mph")
table(obj_macro$Label_type)

input.mat = obj_macro@assays$RNA@layers$counts
row.names(input.mat) = row.names(obj2)
colnames(input.mat) = colnames(obj_macro)
obj_macro = CreateSeuratObject(counts = input.mat, meta.data = obj_macro@meta.data, min.cells = 50)

obj_macro <- FindVariableFeatures(obj_macro, verbose = F, nfeatures = 500)
obj_macro = NormalizeData(obj_macro)
obj_macro <- FindVariableFeatures(obj_macro, verbose = F, nfeatures = 500)
obj_macro <- ScaleData(obj_macro)
obj_macro <- FindVariableFeatures(obj_macro, verbose = F, nfeatures = 500)
obj_macro <- RunPCA(obj_macro, verbose = FALSE)
obj_macro <- obj_macro %>% RunUMAP(dims = 1:5) %>% FindNeighbors(dims = 1:5)
obj_macro <- obj_macro %>% FindClusters(resolution = 0.1)

obj_macro$alra_snn_res.0.1 = obj_macro$RNA_snn_res.0.1

DimPlot(obj_macro, reduction = "umap", label = TRUE, pt.size = .5, cols = my36colors, group.by = c("alra_snn_res.0.1", "Label_type"), raster = F) & theme(aspect.ratio = 1)

Idents(obj_macro) = obj_macro$Label_type
diff.metab.celltype <- FindMarkers(obj_macro, ident.1 = "Positive", ident.2 = "Negative", logfc.threshold = -1, return.thresh = 1)
diff.metab.celltype$gene = row.names(diff.metab.celltype)

GetSignature = function(obj_random, gmtfile = "/www/data/signature/APC.gmt", usecores = 16) {
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

obj_random = obj_macro
signature.matrix.metab = GetSignature(obj_random, gmtfile = "/www/data/signature/macrophage_basic_signature.gmt")
row.names(signature.matrix.metab) = colnames(obj_random)
obj_random2 = CreateSeuratObject(t(signature.matrix.metab), meta.data = obj_random@meta.data)
Idents(obj_random2) = obj_random2$Label_type

wilcox.test(signature.matrix.metab[,1] ~ obj_random2$Label_type)
wilcox.test(signature.matrix.metab[,2] ~ obj_random2$Label_type)
wilcox.test(signature.matrix.metab[,3] ~ obj_random2$Label_type)
wilcox.test(signature.matrix.metab[,4] ~ obj_random2$Label_type)
wilcox.test(signature.matrix.metab[,5] ~ obj_random2$Label_type)
wilcox.test(signature.matrix.metab[,6] ~ obj_random2$Label_type)
wilcox.test(signature.matrix.metab[,7] ~ obj_random2$Label_type)

VlnPlot(obj_random2, features = c(colnames(signature.matrix.metab)), pt.size = 0, log = F, ncol = 7) & geom_boxplot(outlier.shape = NA) & scale_fill_manual(values = my36colors) & theme(aspect.ratio = 2)

obj_random2 = NormalizeData(obj_random2)

library(dittoSeq)
dittoHeatmap(obj_random2, colnames(signature.matrix.metab), assay = "RNA", slot = "counts", scale = "row", scaled.to.max = F, show_colnames = FALSE, show_rownames = T, annot.by = c("Label_type"))

immunomodulator = read.table("/home/wuyingcheng/data/reference/immunomodulator.txt", header = F, sep = "\t")
table(immunomodulator$V2)
immunomodulator = subset(immunomodulator, immunomodulator$V2 %in% c("chemokine", "receptor"))

obj_macro = subset(obj2, predicted.celltype_l1 == "Mph")
table(obj_macro$Label_type)
obj_macro <- FindVariableFeatures(obj_macro, verbose = F, nfeatures = 10000)
obj_macro = NormalizeData(obj_macro)

Idents(obj_macro) = obj_macro$Label_type
diff.metab.celltype <- FindMarkers(obj_macro, ident.1 = "Positive", ident.2 = "Negative", logfc.threshold = -1, return.thresh = 1)
diff.metab.celltype$gene = row.names(diff.metab.celltype)

diff.metab.celltype = subset(diff.metab.celltype, gene %in% immunomodulator$V3)

diff.metab.celltype$logp = -log(diff.metab.celltype$p_val)

diff.metab.celltype$sig = "ns"
diff.metab.celltype$sig[diff.metab.celltype$avg_log2FC > 1.5 & diff.metab.celltype$p_val < .01] = "up"
table(diff.metab.celltype$sig)

diff.metab.celltype$size = abs(diff.metab.celltype$avg_log2FC)

diff.metab.celltype_label = subset(diff.metab.celltype, sig == "up")

quantile(diff.metab.celltype$avg_log2FC)
diff.metab.celltype$avg_log2FC[diff.metab.celltype$avg_log2FC < -3] = -3

ggplot(data = diff.metab.celltype, aes(x = avg_log2FC, y = logp, color = sig)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = diff.metab.celltype_label, mapping = aes(label = gene), show.legend = F, max.overlaps = 100000, size = 3) + 
  theme_bw() +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("grey50", my36colors[2])) +
  theme(legend.position = "right", legend.box = "right", axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlim(-4, 4) +
  NULL

table(diff.metab.celltype$V2)