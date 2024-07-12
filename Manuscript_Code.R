#MERGE PROJECT PROCESSED FILES
library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(foreach)
library(doParallel)
library(stringr)



tmp = list.files(pattern = 'barcodes.tsv.gz', recursive = T)
files = str_remove(tmp,basename(tmp))



alldata<- lapply(files,
                      function(sample){
                        print(sample)
                        Read10X(sample)
                      })

out = strsplit(tmp, '/', fixed =T)

samples = lapply(out, function(sample){
  sample[1]
})


samples = as.matrix(samples)
names(out) = samples
names(alldata) = samples
#combine
combined_SeuratObject <- lapply(samples,
                                function(sample){
                                  print(sample)
                                  CreateSeuratObject(alldata[[sample]], project = sample)
                                })
merge_SeuratObject = merge(combined_SeuratObject[[1]],combined_SeuratObject[-c(1)])

merge_SeuratObject[["percent.mt"]] <- PercentageFeatureSet(merge_SeuratObject, pattern = "^MT-")

merge_SeuratObject = subset(merge_SeuratObject, subset = percent.mt <=25 & nCount_RNA < 25000 & nFeature_RNA < 6000 & nCount_RNA > 100)


merge.list <- SplitObject(merge_SeuratObject, split.by = 'orig.ident')

merge = lapply(X = merge.list, FUN = function(x){
  x = subset(x,subset = nCount_RNA > 100)
  counts = x@assays$RNA@counts
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- decontX(sce)
  x <- CreateSeuratObject(round(decontXcounts(sce)))
  
})

merge.list <- lapply(X = merge.list, FUN = function(x){
  x = subset(x, subset = nCount_RNA > 100 & percent.mt <= 25)
})

#Doublet Detection Samplewise
merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

merge.list <- lapply(X = merge.list, FUN = function(x) {
  print(unique(x$orig.ident))
  x <- ScaleData(x, features = VariableFeatures(x), verbose = T)
  x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
})

merge.list <- lapply(X = merge.list, FUN = function(x) {
suppressMessages(require(DoubletFinder))
nExp <- round(ncol(x) * 0.05)  # expect 4% doublets
x <- doubletFinder_v3(x, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:20)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(x@meta.data)[grepl("DF.classification", colnames(x@meta.data))]

table(x@meta.data[,DF.name])
x = x[, x@meta.data[,DF.name] == "Singlet"]
})


merge <- merge(
  x = merge.list[[1]],
  y = merge.list[c(-1)]
)

#Merge all Objects that are currently in a List format into an actual SO format


#Integrate:

#Add metadata State
#Normalize
merge <- NormalizeData(object = merge, normalization.method = "LogNormalize", scale.factor = 10000)

#Percent Mitochondrial Genes
#QC Metric used to remove cells with overabundant Mitochondrial genes, typically associated with nuclear wash out during sequencing
merge[["percent.mt"]] <- PercentageFeatureSet(merge, pattern = "^MT-")

#Plot the nFeatures/counts/% Mito to get general idea about the quality of your data
pdf('QC_VlnPlot.pdf', height = 20, width = 20)
VlnPlot(merge, features = c("nFeature_RNA", "nCounts_RNA", "percent.mt"), ncol = 3, pt.size = .0)
dev.off()

#FIND VARIABLE GENES
merge<- FindVariableFeatures(object = merge, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merge_SeuratObject<- CellCycleScoring(merge_SeuratObject, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
merge_SeuratObject$CC.Difference <- merge_SeuratObject$S.Score - merge_SeuratObject$G2M.Score

save(merge,file="merge.RData")

#THIS STEP MAY TAKE A VERY LONG TIME
#Scale Data
merge<- ScaleData(object = merge, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"), features = VariableFeatures(merge))


#Data Visualization
#Run PCA and Determine Dimensions for 90% Variance
merge <- RunPCA(object = merge, features = VariableFeatures(object = merge))
merge <- RunHarmony(merge, "orig.ident", plot_convergence = F)
stdev <- merge@reductions$harmony@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)
#Find Neighbors + Find CLusters (without integration

merge <- FindNeighbors(merge, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")
merge <- FindClusters(merge, verbose = TRUE, resolution = 0.4)
merge <- RunUMAP(merge, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)


#Run UMAP and get unlabeled cluster UMAP and violin plot (without integration)
DimPlot(object = merge, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
merge[["UMAP_Clusters"]] <- Idents(object = merge)

#UMAP - Look how object looks pre-integration
merge = AddMetaData(merge, meta.data = meta)
DimPlot(object = merge, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(object = merge, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'GSE.SRA..Study.')
DimPlot(object = merge, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'State')


#Integrate 

#Integrate:
Idents(object = merge) <- 'GSE.SRA..Study.'
merge.list <- SplitObject(merge, split.by = 'GSE.SRA..Study.')

merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

merge.features <- SelectIntegrationFeatures(object.list = merge.list)
merge.list <- lapply(X = merge.list, FUN = function(x) {
  x <- ScaleData(x, features = merge.features, verbose = FALSE)
  x <- RunPCA(x, features = merge.features, verbose = FALSE)
})

merge.anchors <- FindIntegrationAnchors(object.list = merge.list, anchor.features = merge.features, reduction = "rpca")
#READ THE OUTPUT OF THIS FUNCTION
#Make sure the smallest number in Found (number) anchors is greater 100 and if not use the lowest anchor value in the k.weight parameter 

merge.integrated <- IntegrateData(anchorset = merge.anchors, k.weight = 100 )

DefaultAssay(merge.integrated) <- "integrated"
merge.integrated<- FindVariableFeatures(object = merge.integrated, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

merge<- ScaleData(object = merge, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"), features = VariableFeatures(merge))
merge.integrated <- RunPCA(merge.integrated, npcs = 30, verbose = FALSE)


#Calculate 90% Variance:
merge <- RunHarmony(merge, "orig.ident", plot_convergence = F)

st_dev <- merge.integrated@reductions$harmony@stdev
var <- st_dev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
sum(var[1:PCNum])/ sum(var)

#table(merge.integrated@meta.data$ID,merge.integrated@meta.data$DiseaseState)

merge.integrated <- FindNeighbors(object = merge.integrated, dims = 1:PCNum)
merge.integrated <- FindClusters(object = merge.integrated, resolution = 1.8)
merge.integrated <- RunUMAP(merge.integrated, reduction = "pca", dims = 1:PCNum, verbose = F)
DimPlot(merge.integrated, reduction = "umap", label = TRUE, split.by = 'Cell.Type',raster = F) + NoLegend()


DefaultAssay(merge.integrated) <- "RNA"

out =  FindAllMarkers(merge.integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
markers = out %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  
dir.create("./UMAP_plots/")
for (i in 1:length(gene_list)){
  file_png_name<-paste("./UMAP_plots/UMAP_",gene_list[i],".pdf",sep="")
  fplot<-FeaturePlot(object = merge.integrated, features = gene_list[i])
  pdf(file_png_name, width=8.5, height=8)
  print(fplot)
  dev.off()
}



###Cell Labeling
##Fibroblast
genes = c('PDPN','DCN','LUM','LRRC15','SRFP2','COL1A1','COL1A2','COL3A1','PI16','PDGFRA','FN1','VIM')
pdf('Fibroblast_Markers.pdf', height = 8, width = 8)
FeaturePlot(merge.integrated, features = genes, pt.size = 0)
dev.off()

#Perivascular cells
genes = c('ACTA2','PDGFRB','RGS5')
pdf('PV_Markers.pdf', height = 8, width = 8)
VlnPlot(combined, features = genes, pt.size = 0)
dev.off()

##Proliferating cells
pdf('Proliferating_Markers.pdf', height = 4, width = 4)
genes = c('MKI67','STMN1','TOP2A')
dev.off()


top3 = out %>% group_by(cluster) %>% top_n(avg_log2FC)
pdf('DotPlot_top3.pdf',height = 8, width = 24)
DotPlot(merge, features =unique(top3$gene),  dot.scale = 8, cols = c("blue", "red"))  + coord_flip()
dev.off()


library(ggplot2)
Ductal = readRDS('Ductal_Object.rds')
tmp = table(Ductal$seurat_clusters, Ductal$Name)

# create a dataset
cluster <- rep(rownames(tmp),each =ncol(tmp))
condition <- rep(colnames(tmp) , nrow(tmp))
value <- c(t(tmp))
data <- data.frame(cluster,condition,value)
data$cluster = as.factor(data$cluster)
colnames(data) = c('Name','Cluster', 'Proportion')
# Grouped
p1 <- ggplot(data, aes(fill= Cluster, y=Proportion, x=Name)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() +
  scale_fill_manual("legend", values = c("A" = "tomato", "B" = "darkred", "C" = "maroon", "D" = "dodgerblue3","Acinar" = "darkblue","E" = "springgreen"))
pdf('Ductal_Barplot.pdf', height = 4, width = 4)
print(p1)
dev.off()


Mye = readRDS('Myeloid_Object.rds')
tmp = table(Mye$Clusters, Mye$Name)

# create a dataset
cluster <- rep(rownames(tmp),each =ncol(tmp))
condition <- rep(colnames(tmp) , nrow(tmp))
value <- c(t(tmp))
data <- data.frame(cluster,condition,value)
data$cluster = as.factor(data$cluster)
colnames(data) = c('Name','Cluster', 'Proportion')
p1 = ggplot(data, aes(fill= Cluster, y=Proportion, x=Name)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() +
  scale_fill_manual("legend",values = c("Neutrophil" = "red","Macrophage" =  "#00FFFF","Monocyte"= "orange3","Dendritic" ="yellow4"))
pdf('Myeloid_Barplot.pdf', height = 4, width = 4)
print(p1)
dev.off()

#hold = data.frame(Name = rep(c('Name1','Name2'), each = 4), Cluster = rep(unique(data$Cluster),  2), Proportion = rep(0,4))
#hold1 = data.frame(Name = rep(c('Name3','Name4'), each = 4), Cluster = rep(unique(data$Cluster),  2), Proportion = rep(0,4))
#hold2 = data.frame(Name = rep(c('Name5','Name6'), each = 4), Cluster = rep(unique(data$Cluster),  2), Proportion = rep(0,4))

#data = rbind(data[1:88,],hold,data[89:256,],hold1, data[257:1456,],hold2,data[1457:nrow(data),])


Ductal = readRDS('Ductal_Object.rds')
tmp = table(Ductal$seurat_clusters, Ductal$TreatmentStatus)

# create a dataset
cluster <- rep(rownames(tmp),each =ncol(tmp))
condition <- rep(colnames(tmp) , nrow(tmp))
value <- c(t(tmp))
data <- data.frame(cluster,condition,value)
data$cluster = as.factor(data$cluster)
colnames(data) = c('Name','Cluster', 'Proportion')
p1 = ggplot(data, aes(fill= Cluster, y=Proportion, x=Treatment)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() +
scale_fill_manual(values = c("1" = "darkorange", "2" = "darkgreen", "3" = "darkblue", "4" = "darkred","5" = "deeppink","6" = "navajowhite4","7" = "springgreen","8" = "darkmagenta", "9" = "orangered"))
pdf('Ductal_Treatment.pdf',height = 6, width =6)
print(p1)
dev.off()

### Running inferCNV
library(infercnv)

genes = read.table('./Genes.txt', row.names = 1)
rownames(genes) = make.unique(genes[,1])
genes = genes[,-c(1)]
merge = subset(merge, subset = Clusters %in% c('DUCTAL','TNK','MYELOID','ENDOTHELIAL'))
options(scipen = 100)
for(i in unique(merge$Name)){
  if(file.exists(i)) next()
  tmp = subset(merge, subset = Name %in% i)
  
  counts = tmp@assays$RNA@counts
  label = tmp@meta.data[,c('Clusters')]
  label = as.data.frame(label)
  rownames(label) = colnames(tmp)
  out = names(which(table(tmp$Clusters) > 2))
  label1 = label[label$label %in% out,]
  label1 = as.character(label1)
  label1 = as.data.frame(label1)
  rownames(label1) = rownames(label)[label$label %in% out]
  out = out[!out %in% c('DUCTAL','Prolif/Ductal')]
  #label1$label1[label1$label1 == 'Prolif/Ductal'] = 'Ductal'
  counts = counts[,colnames(counts) %in% rownames(label1)]
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                    annotations_file=label1,
                                    gene_order_file=genes,
                                    ref_group_names= out) 

  dir.create(i)
  outdir = paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Analysis_V4/Migraine/inferCNV/',i,'/')
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=outdir, 
                             denoise=TRUE,no_prelim_plot = T,plot_probabilities=F,
                             HMM=TRUE, tumor_subcluster_partition_method = "random_trees",cluster_by_groups = F, # cluster or not
                             analysis_mode = "subclusters",
                             HMM_type = "i6",
                             BayesMaxPNormal = 0,
                             num_threads = 15 )
}

out = NULL
tmp = list.files(getwd(), pattern = '17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.infercnv_obj', recursive = T)
tmp1 = list.files(getwd(), pattern = '17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.pred_cnv_regions.dat', recursive = T)
for(i in 1:length(tmp)){
  expr = readRDS(tmp[i])
  regions = read.table(tmp1[i],header = T)
  hold = NULL
  for(j in 1:nrow(regions)){
    out1 = expr@expr.data[(rownames(expr@gene_order)[(expr@gene_order$chr %in% regions$chr[j] & expr@gene_order$start >= regions$start[j]) & (expr@gene_order$chr %in% regions$chr[j] & expr@gene_order$stop <= regions$end[j])]),]
  #hold = rbind(hold, colMeans((out-1)^2))
    hold = rbind(hold, colMeans((out1)))
    
  }
  hold[hold>5] = 5
  hold = (hold - 3)/2
  out = c(out, colMeans((hold)^2))
}



##Signature Scoring
library(survival)
library(survminer)
markers = c('TumorMarkers_NoMet.csv',
            'CoarseFibroblast_OnlyTumor.csv',
            'CoarseMyeloidMarkersNOTumor.csv')

            cox = NULL
scores = NULL
for(k in markers){
  out = read.csv(k)
  out = out[which(out$gene %in% rownames(data)),]
  out = out[!grepl('^MT',out$gene) & !grepl('^RPS',out$gene)& !grepl('^RPL',out$gene),]
  #out1 = out
  out1 = out %>% group_by(cluster) %>% top_n(25, avg_log2FC)
#top50 = out %>% group_by(cluster) %>% top_n(50, avg_log2FC)
#gene_set = list(top50$gene[top50$cluster == unique(top50$cluster)[1]])
gene_set = list(out1$gene[out1$cluster == unique(out1$cluster)[1] & out1$p_val_adj < 0.05])

for(i in 2:length(unique(out1$cluster))){
  gene_set[[i]] = out1$gene[out1$cluster == unique(out1$cluster)[i] & out1$p_val_adj < 0.05]}
names(gene_set) = unique(out1$cluster)
gene_sets = gene_set
system.time(assign('res', ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)))
mat = (res - rowMeans(res))/(rowSds(as.matrix(res)))[row(res)]
res = mat
scores = rbind(scores, res)
for(i in c(1:nrow(res))){
  if (quantile(res[i,],0.67) == 0) 
    clin$bin = ifelse(res[i,] > 0,1,0)
  else
    clin$bin = ifelse(res[i,] >= quantile(res[i,],0.67),1,0)

  km_fit_add <- survfit(Surv(time=clin$days_to_last_followup, event=clin$ vital_status) ~ bin, data = clin)
  plot = ggsurvplot(km_fit_add,conf.int=F,pval = T, surv.size=1.5,legend.labs = c('Low Score','High Score'),palette = c('blue','red'))+labs(x="Time (Months)")
  pdf(paste0('Score_67th_',rownames(res)[i],'_TCGA_SSGSEA_AllFDR25.pdf'), height = 6, width = 6)
  print(plot)
  dev.off()
  
  
}
sur = NULL
for(i in 1:nrow(res)){
  tmp <- summary(coxph( Surv(time=clin$days_to_last_followup, event=clin$vital_status) ~res[i,]))$coefficients
  sur = rbind(sur,tmp)
}
rownames(sur) = rownames(res)
cox = rbind(cox,sur)

}
write.csv(cox, file = 'SurvivalScores.csv')
write.csv(scores, file = 'SignatureScores.csv')

###Ductal Reactome Dot Plot
DuctalMarkers = read.csv('TumorMarkers_NoMet.csv')
for(i in 1:length(unique(out$cluster))){
  genes = out[out$cluster %in% unique(out$cluster)[i],]
  genes1 = genes[genes$p_val_adj < 0.05,]
  genes1 = genes1[!grepl('^MT',genes1$gene) & !grepl('^RPS',genes1$gene)& !grepl('^RPL',genes1$gene),]
  geneID = mapIds(org.Hs.eg.db, genes1$gene, 'ENTREZID', 'SYMBOL')
  genes_out[[i]]= enrichPathway(geneID)}

for(i in c(1:8)){
  genes[[i]] = peakAnnoList[[i]]@gene
}


compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 8, title = "Reactome Pathway Enrichment Analysis") + coord_flip() + theme(axis.text.x = element_text(angle = 90, size = 8))



###Forest Plot
df1 = read.csv('SurvivalScores')
ggforestplot::forestplot(
  df = df1[df1$X %in% c('Neutrophil','Macrophage','Dendritic','Monocyte'),],name = X,
  estimate = `coef`,
  pvalue = `Pr...z..`,
  psignif = 0.05,se = `se.coef.`,
  xlab = "",
  title = "",ci = 0.95,
  colour = Study, logodds = T
) + coord_flip() + theme_classic() + theme(axis.text.x = element_text(angle = 90))



##Spatial Analysis
#Processed ST data
ST = readRDS('Spatial_Objects.rds')
lapply(X = ST, FUN = function(x) {
  print(as.character(unique(x$Sample_ID)))
  print(dim(x))
  x <- NormalizeData(x, assay = "Spatial")
  x@meta.data$Decon_topics <-
    paste(x@meta.data$Location,
          x@meta.data$seurat_clusters,
          sep = "_")
  expr_values = as.matrix(x@assays$Spatial@data) #ST expr log
  nolog_expr = 2 ^ (expr_values) - 1 #ST expr nolog
  meta_data <-
    x@meta.data[, c("nCount_Spatial", "Decon_topics", "Location")]
  
  #Signature score
  for (cluster in names(clustermarkers_list)) {
    cluster_markers = clustermarkers_list[[cluster]][1:25]
    cluster_score <-
      apply(x@assays$Spatial@data[rownames(x@assays$Spatial@data) %in% cluster_markers, ], 2, mean)
    meta_data <- cbind(meta_data, cluster_score)
  }
  colnames(meta_data) <-
    c("nCount_Spatial",
      "Decon_topics",
      "Location",
      names(clustermarkers_list))
  
  #filter st and sc feature
  intersect_gene = intersect(rownames(sig_exp), rownames(nolog_expr))
  filter_sig = sig_exp[intersect_gene, ]
  filter_expr = nolog_expr[intersect_gene, ]
  filter_log_expr = expr_values[intersect_gene, ]
  
  #enrichment analysis
  enrich_matrix <-
    get_enrich_matrix(filter_sig = filter_sig,
                      clustermarkers_list = clustermarkers_list)
  enrich_result <-
    enrich_analysis(filter_log_expr = filter_log_expr,
                    enrich_matrix = enrich_matrix)
  
  dwls_results <-  spot_proportion_initial(
    enrich_matrix = enrich_matrix,
    enrich_result = enrich_result,
    filter_expr = filter_expr,
    filter_sig = filter_sig,
    clustermarkers_list = clustermarkers_list,
    meta_data = meta_data,
    malignant_cluster = "Myeloid",
    tissue_cluster = NULL,
    stromal_cluster = NULL
  )
  
  enrich_spot_proportion <- dwls_results
  binary_matrix <- ifelse(enrich_spot_proportion >= 0.01, 1, 0)
  
  spot_proportion <- spot_deconvolution(
    expr = filter_expr,
    meta_data = meta_data,
    ct_exp = filter_sig,
    enrich_matrix = enrich_matrix,
    binary_matrix = binary_matrix
  )
  
  
  DeconData <- as.data.frame(t(spot_proportion)) %>% tibble::rownames_to_column(., var = "cell_ID")
  
  write.csv(DeconData, file = paste0('Deconvolution_',as.character(unique(x$Sample_ID)),'_Spatial.csv'))
  gc()
})

plot_col <- colnames(DeconData)[2:ncol(DeconData)]

DeconPieplot(
  DeconData = DeconData,
  TumorST = TumorST,
  plot_col = plot_col,
  img_path = img_path,
  pie_scale = 0.4,
  scatterpie_alpha = 0.8,
  border_color = "grey"
)




tmp = list.files( 'lowres', recursive = T)

out = strsplit(tmp, '/', fixed =T)

samples = lapply(out, function(sample){
  sample[1]
})

img = lapply(out, function(sample){
  sample[2]
})


for(i in 8:length(samples)){
Duc = read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Analysis_V4/Migraine/Ductal_Subcluster/Deconvolution_',samples[[i]],'_New.csv'), row.names = 1)
img_path = paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/phs002371.v3.p1/',out[[i]])
tmp = readRDS(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/',samples[[i]],'_decon_FB_Coarse_thres_0.4.rds'))
Fib = tmp@meta.data[,31:38]
Fib$cell_ID = rownames(Fib)
tmp = anti_join(Duc,Fib, by = 'cell_ID')
tmp1 = anti_join(Fib,Duc, by = 'cell_ID')
DeconData = merge(tmp, tmp1, by = 'cell_ID', all = T)
DeconData[is.na(DeconData)] = 0
plot_col <- colnames(DeconData)[2:ncol(DeconData)]

TumorST =merge[[which(names(merge) %in% samples[[i]])]]

Dplot = DeconPieplot(
  DeconData = DeconData,
  TumorST = TumorST,
  plot_col = plot_col,
  img_path = img_path,
  pie_scale = 0.4,
  scatterpie_alpha = 0.8,
  border_color = "grey"
) +
  scale_fill_manual(values = c("X0" = "darkorange", "X1" = "darkgreen", "X2" = "darkblue", "X3" = "darkred","X4" = "deeppink","X5" = "navajowhite4","X7" = "springgreen","X8" = "darkmagenta", "X9" = "orangered","myFb_MMP11" = "cornsilk","myFb_metabolism" =  "deepskyblue","myFb_CXCL"= "forestgreen", "iFb_CXCL" ="plum", "iFb_C7" = "gray0",  "myFb_VMP1" = "pink2","myFb_CXCL10" = "darkorchid4", "Vascular_Fb" = "blue"))

png(paste0('DuctalSpatial_',samples[[i]],'.png'), height = 8, width = 8, units = 'in', res = 600)
print(Dplot)
dev.off()
}
out = list.files('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/Spatial_Realign/', 'lowres', recursive = T)
out = out[!grepl('fork',out)]

for(i in c(1:8)){
  Duc = read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Analysis_V4/Migraine/Ductal_Subcluster/Deconvolution_',names(merge)[i],'_New.csv'), row.names = 1)
  img_path = paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/Spatial_Realign/',out[[i]])
  tmp = readRDS(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/',names(merge)[i],'_decon_FB_Coarse_thres_0.4.rds'))
  Fib = tmp@meta.data[,31:38]
  Fib$cell_ID = rownames(Fib)
  tmp = anti_join(Duc,Fib, by = 'cell_ID')
  tmp1 = anti_join(Fib,Duc, by = 'cell_ID')
  DeconData = merge(tmp, tmp1, by = 'cell_ID', all = T)
  DeconData[is.na(DeconData)] = 0
  plot_col <- colnames(DeconData)[2:ncol(DeconData)]
  
  TumorST =merge[[which(names(merge) %in% names(merge)[i])]]
  
  Dplot = DeconPieplot(
    DeconData = DeconData,
    TumorST = TumorST,
    plot_col = plot_col,
    img_path = img_path,
    pie_scale = 0.4,
    scatterpie_alpha = 0.8,
    border_color = "grey"
  ) +
    scale_fill_manual(values = c("X0" = "darkorange", "X1" = "darkgreen", "X2" = "darkblue", "X3" = "darkred","X4" = "deeppink","X5" = "navajowhite4","X7" = "springgreen","X8" = "darkmagenta", "X9" = "orangered","myFb_MMP11" = "cornsilk","myFb_metabolism" =  "deepskyblue","myFb_CXCL"= "forestgreen", "iFb_CXCL" ="plum", "iFb_C7" = "gray0",  "myFb_VMP1" = "pink2","myFb_CXCL10" = "darkorchid4", "Vascular_Fb" = "blue"))
  
  png(paste0('DuctalSpatial_',names(merge)[i],'.png'), height = 8, width = 8, units = 'in', res = 600)
  print(Dplot)
  dev.off()
}


coords = list.files('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/phs002371.v3.p1', 'tissue_positions', recursive = T)


for(i in c(1:5,10:length(samples))){
  Duc = read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Analysis_V4/Migraine/Ductal_Subcluster/Deconvolution_',samples[[i]],'_New.csv'), row.names = 1)
  img_path = paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/phs002371.v3.p1/',samples[[i]],'/',img[[i]])
  tmp = readRDS(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/',samples[[i]],'_decon_FB_Coarse_thres_0.4.rds'))
  Fib = tmp@meta.data[,31:38]
  Fib$cell_ID = rownames(Fib)
  tmp = anti_join(Duc,Fib, by = 'cell_ID')
  tmp1 = anti_join(Fib,Duc, by = 'cell_ID')
  DeconData = merge(tmp, tmp1, by = 'cell_ID', all = T)
  DeconData[is.na(DeconData)] = 0
  #plot_col <- colnames(DeconData)[2:ncol(DeconData)]
  
  #TumorST =merge[[which(names(merge) %in% samples[[i]])]]
  cell_data =read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/phs002371.v3.p1/',coords[grepl(samples[[i]],coords)]),header = F)
  DeconData = DeconData[DeconData$cell_ID %in% paste0(as.name(samples[[i]]),':',cell_data$V1),]
  cell_data = cell_data[match(DeconData$cell_ID,paste0(as.name(samples[[i]]),':',cell_data$V1)),]
  DeconData$cell_ID = gsub(paste0(samples[[i]],':'),'',DeconData$cell_ID)
  
   #DeconData = DeconData[,-c(1)]
  #CellType = colnames(DeconData)[DeconData > 0]
  names = NULL
  CellType = NULL
  for(k in 1:nrow(DeconData)){
  names = c(names, rep(DeconData$cell_ID[k],rowSums(DeconData[k,-c(1)] >0)))
  CellType = c(CellType,colnames(DeconData)[which(DeconData[k,] > 0)][-c(1)])}
  
  #CellType = CellType[!grepl('cell_ID',CellType)]
  #cells = DeconData$cell_ID[]
  cell_data = as.data.frame(cell_data)
  cell_data = cell_data[match(names, cell_data$V1),]
  for(k in 1:10){
    cell_data$V5[duplicated(cell_data$V5)] = cell_data$V5[duplicated(cell_data$V5)] + 1
  }
  cell_data$CellType = CellType
  cell_data$Study = rep(samples[[i]],nrow(cell_data))
  
  all = rbind(all, cell_data)}

coords = list.files('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/Spatial_Realign/', 'tissue_positions', recursive = T)
coords = coords[!grepl('fork', coords)]
for(i in c(1:8)){
  Duc = read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Analysis_V4/Migraine/Ductal_Subcluster/Deconvolution_',names(merge)[i],'_New.csv'), row.names = 1)
  img_path = paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/Spatial_Realign/',out[[i]])
  tmp = readRDS(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/',names(merge)[i],'_decon_FB_Coarse_thres_0.4.rds'))
  Fib = tmp@meta.data[,31:38]
  Fib$cell_ID = rownames(Fib)
  tmp = anti_join(Duc,Fib, by = 'cell_ID')
  tmp1 = anti_join(Fib,Duc, by = 'cell_ID')
  DeconData = merge(tmp, tmp1, by = 'cell_ID', all = T)
  DeconData[is.na(DeconData)] = 0
  #plot_col <- colnames(DeconData)[2:ncol(DeconData)]
  
  #TumorST =merge[[which(names(merge) %in% samples[[i]])]]
  cell_data =read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/Spatial_Realign/',coords[i]),header = F)
  DeconData = DeconData[DeconData$cell_ID %in% paste0(as.name(names(merge)[i]),':',cell_data$V1),]
  cell_data = cell_data[match(DeconData$cell_ID,paste0(as.name(names(merge)[i]),':',cell_data$V1)),]
  DeconData$cell_ID = gsub(paste0(names(merge)[i],':'),'',DeconData$cell_ID)
  
  #DeconData = DeconData[,-c(1)]
  #CellType = colnames(DeconData)[DeconData > 0]
  names = NULL
  CellType = NULL
  for(k in 1:nrow(DeconData)){
    names = c(names, rep(DeconData$cell_ID[k],rowSums(DeconData[k,-c(1)] >0)))
    CellType = c(CellType,colnames(DeconData)[which(DeconData[k,] > 0)][-c(1)])}
  
  #CellType = CellType[!grepl('cell_ID',CellType)]
  #cells = DeconData$cell_ID[]
  cell_data = as.data.frame(cell_data)
  cell_data = cell_data[match(names, cell_data$V1),]
  for(k in 1:10){
    cell_data$V5[duplicated(cell_data$V5)] = cell_data$V5[duplicated(cell_data$V5)] + 1
  }
  cell_data$CellType = CellType
  cell_data$Study = rep(names(merge)[i],nrow(cell_data))
  
  all = rbind(all, cell_data)}
n = nrow(all)
y <- matrix(nrow = n, ncol = n)
spe1 <- SpatialExperiment(
  assay = y,
  colData =all, spatialCoordsNames = c('V5','V6'))


library(imcRtools)

spe1 = buildSpatialGraph(
  spe1,
  'Study',
  type =  "knn",
  k = 12,
  directed = T,
  max_dist =200 ,
  threshold = NULL,
  coords = c("V5", "V6")
)

spe1 <- buildSpatialGraph(spe1, img_id = "Study",
                                 type = "knn",
                                 threshold = 200, coords = c("V5", "V6")

out <- testInteractions(spe1,
                        group_by = "Study",
                        label = "CellType",
                        method = "histocat",
                        colPairName = "knn_interaction_graph",
                        BPPARAM = BiocParallel::SerialParam(RNGseed = 123))


DeconData[is.na(DeconData)] = 0
out = Barplot(DeconData, plot_col)
out = cbind(out, rep(names[i],nrow(out)))
combined = rbind(combined, out )}


Barplot <- function(DeconData = DeconData,
                         plot_col = plot_col) {
  metadata <- DeconData
  #metadata[, "Location"] <- TumorST@meta.data$Location[match(metadata$cell_ID, rownames(TumorST@meta.data))]
  metadata[,"Location"] <- 'nMAL'
  Location <- do.call(rbind, lapply(unique(metadata$Location), function(x) {
    metadata_split <- metadata[metadata$Location == x, ] %>% as.data.frame()
    sub <- colSums(metadata_split[, plot_col], na.rm = T) %>%
      data.frame() %>%
      tibble::rownames_to_column() %>%
      set_colnames(., c("Types", "Sum")) %>%
      dplyr::mutate(Per = 100 * Sum / sum(Sum))
  }))
  #Location$Group <- rep(unique(metadata$Location), each = length(plot_col))
  Location$Types <- factor(Location$Types, levels = plot_col)
  return(Location)}
  
  colnames(combined)[4] = 'Sample'
ggplot(combined, aes(x = Sample, y = Per, fill = Types)) +
    geom_bar(stat = "identity") +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("X0" = "darkorange", "X1" = "darkgreen", "X2" = "darkblue", "X3" = "darkred","X4" = "deeppink","X5" = "navajowhite4","X7" = "springgreen","X8" = "darkmagenta", "X9" = "orangered","Neutrophil" = "red","Macrophage" =  "#00FFFF","Monocyte"= "orange3","Dendritic" ="yellow4","myFb_MMP11" = "cornsilk","myFb_metabolism" =  "deepskyblue","myFb_CXCL"= "forestgreen", "iFb_CXCL" ="plum", "iFb_C7" = "gray0",  "myFb_VMP1" = "pink2","myFb_CXCL10" = "darkorchid4", "Vascular_Fb" = "blue"))

  return(Barplot)




cell_groupings <- read.tree(file = paste0(cnv_outdir,"/infercnv.17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt",sep = ""))
infercnv.label <- dendextend::cutree(cell_groupings,k=8)
infercnv.label <- as.data.frame(infercnv.label)
cnv_table <- readRDS(paste0(cnv_outdir,'17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.infercnv_obj'))
cnv_score_table <- as.matrix(cnv_table@expr.data)

cnv_score_tableA <- abs(cnv_score_table-3)

cell_scores_CNV <- as.data.frame(colSums(cnv_score_tableA)) %>% set_colnames(.,c('cnv_score'))
rownames(cell_scores_CNV) <- gsub("\\.",'-',rownames(cell_scores_CNV))
cell_scores_CNV = cell_scores_CNV[match(names(infercnv.label), rownames(cell_scores_CNV)),]
cell_scores_CNVA <- data.frame('CNVLabel' = infercnv.label ,"cnv_score" = cell_scores_CNV)


###Spatial Plots
tmp = list.files('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/phs002371.v3.p1', 'lowres', recursive = T)

out = strsplit(tmp, '/', fixed =T)

samples = lapply(out, function(sample){
  sample[1]
})

img = lapply(out, function(sample){
  sample[2]
})


for(i in c(1:5,8:length(samples))){
  Duc = read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Analysis_V4/Migraine/Ductal_Subcluster/Deconvolution_',samples[[i]],'_New.csv'), row.names = 1)
  img_path = paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/phs002371.v3.p1/',samples[[i]],'/',img[[i]])
  Fib = read.csv(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/Neutrophils/Deconvolution_',samples[[i]],'_Myeloid.csv'),row.names = 1)
  tmp = anti_join(Duc,Fib, by = 'cell_ID')
  tmp1 = anti_join(Fib,Duc, by = 'cell_ID')
  tmp2 = readRDS(paste0('/share/studies/Pancreatic_Cancer_Center/Pancreatic_Cancer_SC_Database/Spatial/',samples[[i]],'_decon_FB_Coarse_thres_0.4.rds'))
  Fib = tmp2@meta.data[,31:38]
  Fib$cell_ID = rownames(Fib)
  tmp3 = merge(tmp, tmp1, by = 'cell_ID', all = T)
  tmp = anti_join(tmp3,Fib, by = 'cell_ID')
  tmp1 = anti_join(Fib,tmp3, by = 'cell_ID')
  DeconData = merge(tmp, tmp1, by = 'cell_ID', all = T)
  DeconData[is.na(DeconData)] = 0
  plot_col <- colnames(DeconData)[2:ncol(DeconData)]
  
  TumorST =merge[[which(names(merge) %in% samples[[i]])]]
  
  Dplot = DeconPieplot(
    DeconData = DeconData,
    TumorST = TumorST,
    plot_col = plot_col,
    img_path = img_path,
    pie_scale = 0.4,
    scatterpie_alpha = 0.8,
    border_color = "grey"
  ) +
    scale_fill_manual(values = c("X0" = "darkorange", "X1" = "darkgreen", "X2" = "darkblue", "X3" = "darkred","X4" = "deeppink","X5" = "navajowhite4","X7" = "springgreen","X8" = "darkmagenta", "X9" = "orangered","Neutrophil" = "red","Macrophage" =  "#00FFFF","Monocyte"= "orange3","Dendritic" ="yellow4","myFb_MMP11" = "cornsilk","myFb_metabolism" =  "deepskyblue","myFb_CXCL"= "forestgreen", "iFb_CXCL" ="plum", "iFb_C7" = "gray0",  "myFb_VMP1" = "pink2","myFb_CXCL10" = "darkorchid4", "Vascular_Fb" = "blue"))
  
  
  #cols = c("myFb_MMP11" = "cornsilk","myFb_metabolism" =  "deepskyblue","myFb_CXCL"= "forestgreen", "iFb_CXCL" ="plum", "iFb_C7" = "gray0",  "myFb_VMP1" = "pink2","myFb_CXCL10" = "darkorchid4", "Vascular_Fb" = "blue")
  
  png(paste0('DuctalSpatial_',samples[[i]],'.png'), height = 8, width = 8, units = 'in', res = 600)
  print(Dplot)
  dev.off()
}

###Correlation Plots
library(corrplot)
scores = read.csv('SignatureScores.csv',row.names = 1)
cort = cor(t(scores))
pdf('CorrelationPlot.pdf', height = 4, width = 4)
corrplot(cort,tl.cex = 0.5,method = 'color',order = 'hclust', col =  rev(COL2('RdBu')), addrect = 3, rect.col = 'blue')
dev.off()