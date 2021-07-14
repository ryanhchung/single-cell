Integrated single-cell RNA sequencing analyses suggest developmental paths of cancer-associated fibroblasts with gene expression dynamics
================
Authors : Heechul(Ryan) Chung, M.S. (heechulrchung@gmail.com) & Chang Ohk Sung, M.D., Ph.D. (co.sung@amc.seoul.kr)

Asan Center for Cancer Genome Discovery In Collaboration with Dana-Farber Cancer Institute, Asan Medical Center

University of Ulsan College of Medicine

``` r
setwd("/Users/ryanmachine/")
library(tidyverse)
library(Matrix)
library(Seurat)
```


### If you don't have a merged raw matrix, then follow this procedure first to merge all 10X runs
### Or, If you already have a merged raw matrix, then go to the next part and load your merged raw matrix
The main paths to use your 10x runs should be the parent path of your file.
For example, my path looks like this : "ovarymatrix/ovary1-T/filtered_feature_bc_matrix/~.gz" (barcode.tsv.gz, features.tsv.gz, matrix.mtx.gz)
Therefore, My parent path is **ovarymatrix**
``` r
ovarypath <- "ovarymatrix"

#Subfolder1 is just used to obtain Subfolder2
subfolder1 <- list.dirs(ovarypath, recursive =TRUE)[-1]

#Subfolder2 -> "ovarymatrix/ovary1-T/filtered_feature_bc_matrix", "ovarymatrix/ovary1-N/filtered_feature_bc_matrix"...
subfolder2 <- list.dirs(subfolder1, recursive=FALSE)

#Generate Seuart object for each sample named ovary_1, ovary_2, ovary_3...
i <- 1
for (file in subfolder2){
  ov_seurat <- Read10X(data.dir = paste0(file))
  assign(paste0("ovary_", i),CreateSeuratObject(counts = ov_seurat, project = "ovary"))
  i <- i+1
}

#Merge all 10x runs
mat_ovarian <- merge(ovary_1, y=c(ovary_2, ovary_3, ovary_4, ovary_5, ovary_6, ovary_7, ovary_8, ovary_9, ovary_10, ovary_11, ovary_12, ovary_13, 
                               ovary_14, ovary_15, ovary_16, ovary_17, ovary_18, ovary_19, ovary_20, ovary_21, ovary_22, ovary_23, ovary_24,
                               ovary_25, ovary_26, ovary_27, ovary_28, ovary_29, ovary_30, ovary_31, ovary_32, ovary_33),
                    add.cell.ids = c("ovary01-N", "ovary01-T", "ovary02-N", "ovary02-T", "ovary03-N", "ovary03-T",
                                     "ovary04-N", "ovary04-T", "ovary05-N", "ovary05-T", "ovary06-N", "ovary06-T", "ovary07-N", "ovary07-T", "ovary08-N",
                                     "ovary08-T", "ovary09-N", "ovary09-T", "ovary10-N", "ovary10-T", "ovary11-T", "ovary14-T", "ovary15-T", "ovary16-T",
                                     "ovary17-T", "ovary18-T", "ovary19-T", "ovary20-T", "ovary21-T", "ovary22-T", "ovary23-T", "ovary24-T", "ovary25-T"),
                    project = "ovary")
```
### Load scRNA-seq merged matrix

``` r
ovariananno <- read.csv("OvarianCancer/OvC_counts/OvC_metadata.csv.gz")
ovarian_fib_tumor <- subset(ovariananno, ovariananno$CellFromTumor==TRUE)
ovarian_fib_normal <- subset(ovariananno, ovariananno$CellFromTumor==FALSE)

feature_ovarian <- paste0("OvarianCancer/OvC_counts/genes.tsv")
barcodes_ovarian <- paste0("OvarianCancer/OvC_counts/barcodes.tsv")
ovarian_mat<- paste0("OvarianCancer/OvC_counts/matrix.mtx")
mat_ovarian <- readMM(file = ovarian_mat)

#Skip below two lines if you don't need to regress out cellcycle.
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ovarian_feature.names = read.delim(feature_ovarian, 
                                  header = FALSE,
                                  stringsAsFactors = FALSE)
ovarian_barcode.names = read.delim(barcodes_ovarian, 
                                  header = FALSE,
                                  stringsAsFactors = FALSE)
colnames(mat_ovarian) = ovarian_barcode.names$V1
rownames(mat_ovarian) = ovarian_feature.names$V2

mat_ovarian <- mat_ovarian[, !duplicated(colnames(mat_ovarian))] ###Remove duplicated colnames if exist###
```
Divide by tumor and normal cells
``` r
mat_ovarian_fib_tumor <- mat_ovarian[,colnames(mat_ovarian) %in% ovarian_fib_tumor$Cell]
mat_ovarian_fib_normal <- mat_ovarian[,colnames(mat_ovarian) %in% ovarian_fib_normal$Cell]
```

### Pre-processing of scRNA-seq raw data
Making tumor fibroblast data - Repeat the same thing with normal fibroblast data too!
``` r
mat_seurat <- CreateSeuratObject(counts = mat_ovarian_fib_tumor, project = "Pan_cancer_col", min.cells = 3, min.features = 200)
#If you don't need percent.mt, then just skip it
mat_seurat[["percent.mt"]] <- PercentageFeatureSet(mat_seurat, pattern = "^MT-")
#Subsetting might be changed based on your standards.
mat_seurat <- subset(mat_seurat, subset = nCount_RNA  > 401 & nFeature_RNA > 201 & nFeature_RNA < 6000 & percent.mt < 25)
mat_seurat <- NormalizeData(mat_seurat)
mat_seurat <- FindVariableFeatures(mat_seurat, selection.method = "vst", mean.cutoff = c(0.125, 3), dispersion.cutoff = c(0.5, Inf))

all.genes <- rownames(mat_seurat)
mat_seurat <- ScaleData(mat_seurat,features = all.genes) 

#Regressing out cellcycle might not be essential to your experiment. 
#If you don't need it, just skip below two lines and perform RunPCA with mat_seurat
mat_cellcycle <- CellCycleScoring(mat_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mat_regressout <- ScaleData(mat_cellcycle, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "orig.ident", "nCount_RNA"), features = all.genes)

mat_pca <- RunPCA(mat_regressout, features = VariableFeatures(object = mat_regressout))

#All numbers below four lines might be changed
ElbowPlot(mat_pca, 30) 
mat_nei <- FindNeighbors(mat_pca, dims = 1:13)
mat_clu <- FindClusters(mat_nei, resolution = c(0.2, 2))
mat.tsne <- RunTSNE(mat_clu, dims = 1:13, method = 'FIt-SNE')

DimPlot(mat.tsne, reduction = 'tsne') + labs(title = 'Ovarian_Fib_Tumor')
```

### SingleR & Marker-gene based fibroblast selection 
It might be changed if you want other subtypes such as mast cells, T-cells, etc.
``` r
#Fibrobalst marker genes -> COL1A1, DCN, BGN
plot1 <- FeaturePlot(mat.tsne, features = "COL1A1")
col1 <- filter(plot1$data, plot1$data$COL1A1 > 0)

plot2 <- FeaturePlot(mat.tsne, features = "DCN")
dcn <- filter(plot2$data, plot2$data$DCN > 0)

plot3 <- FeaturePlot(mat.tsne, features = "BGN")
bgn <- filter(plot3$data, plot3$data$BGN > 0)

ovarian_tfib <- mat_ovarian_fib_tumor[, colnames(mat_ovarian_fib_tumor) %in% c(rownames(col1), rownames(dcn), rownames(bgn))]

library(SingleR)
library(celldex)

hpca.se <- HumanPrimaryCellAtlasData()
data_tumor <- NormalizeData(mat_ovarian_fib_tumor)

singler_ovarian_tumor <- SingleR(test = data_tumor, ref = hpca.se, assay.type.test=1,
                                  labels = hpca.se$label.main)
singler_ovarian_tfib <- rownames(subset(singler_ovarian_tumor, labels == "Fibroblasts"))

overlap_tumor <- intersect(colnames(ovarian_tfib), singler_ovarian_tfib)

library(VennDiagram)
venn.plot <- draw.pairwise.venn(length(colnames(ovarian_tfib)), length(singler_ovarian_tfib), 
                                length(overlap_tumor), c("Original", "SingleR"))

grid.draw(venn.plot)
require(gridExtra)
grid.arrange(gTree(children=venn.plot), top="Ovary_Tumor_Fibroblasts")
```


### Unsupervised k-means clustering of fibroblasts followed by a validation with PCA analysis
``` r
overlap_ovarian_tfib <- as.data.frame(data_tumor[,colnames(data_tumor) %in% overlap_tumor])
#Make overlap_ovarian_nfib with the same process!

#wssplot -> to get optimized k-mer value for clustering
wssplot <- function(tocluster, nc=15, seed=1234){
  wss <- (nrow(tocluster)-1)*sum(apply(tocluster,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(tocluster, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}
wssplot(t(overlap_ovarian_tfib))

#kmeans clustering with k from wssplot
ovary_cluster <- kmeans(t(overlap_ovarian_tfib), 4)

#PCA analysis
#remove all rows with 0s after transpose
ovary_tfib_without_zeros <- t(overlap_ovarian_tfib)[, colSums(t(overlap_ovarian_tfib) != 0) > 0]
pca <- prcomp(ovary_tfib_without_zeros, center = T, scale. = T)

pcaplot <- as.data.frame(pca$x)
temp <- data.frame(matrix(0, ncol = 3, nrow = nrow(pcaplot)))
rownames(temp) <- rownames(pcaplot)
temp$X1 <- pcaplot$PC1
temp$X2 <- pcaplot$PC2
temp$X3 <- ovary_cluster$cluster

#Clusters combined from kmeans and PCA
plot(temp$X1, temp$X2, col = unlist(temp$X3) , main="Fibroblast_Clusters (Ovary)",
     ylab="PC2",
     xlab="PC1",
)
```

### Expression of CAF related genes
``` r
caf.gene <- c("TNC","LOX","LOXL1","TGFB1","VEGFA","PDGFRA","VCAM1","FAP","CAV1","VIM","S100A4")
ovary.caf.maker = NULL
for (i in 1:length(caf.gene)){
  aa = overlap_ovarian_tfib[caf.gene[i],]
  aa = data.frame(aa)
  aa = t(aa)
  aa = data.frame(aa)
  aa$gene = caf.gene[i]
  aa$group = ovary_cluster$cluster
  names(aa) = c("value","gene","group")
  ovary.caf.maker = rbind(ovary.caf.maker, aa)
  print (i)
}

View(ovary.caf.maker)
ovary.caf.maker$x1 = factor(ovary.caf.maker$gene, 
                            levels = c("TNC","LOX","LOXL1","TGFB1","VEGFA","PDGFRA","VCAM1","FAP","CAV1","VIM","S100A4"))
boxplot(value ~ group + x1, data = ovary.caf.maker,
        at = c(1:4,6:9,11:14,16:19,21:24,26:29,31:34,36:39,41:44,46:49,51:54), col = c("lightseagreen","purple"), 
        frame = F, outpch=NA, ylim=c(0,5), main = "ovary.caf", las=2 )

#PRRX1 expression across tumor clusters
boxplot(as.numeric(overlap_ovarian_tfib["PRRX1",]) ~ ovary_cluster$cluster)
```

### Spearman Correlation between tumor clusters and normal clusters
``` r
spearman <- matrix(nrow = length(colnames(overlap_ovarian_tfib)), ncol = length(colnames(overlap_ovarian_nfib)))
for(i in 1:length(colnames(overlap_ovarian_tfib))){
  for (j in 1:length(colnames(overlap_ovarian_nfib))){
    spearman[i, j] <- cor.test(overlap_ovarian_tfib[,i], overlap_ovarian_nfib[,j], method = 'spearman')$estimate
  }
}
rownames(spearman) <- colnames(overlap_ovarian_tfib)
colnames(spearman) <- colnames(overlap_ovarian_nfib)

library(Rfast)
#spearman <- as.matrix(spearman)
r1 <- rowMaxs(spearman, value = TRUE)
r2 <- cbind(1:nrow(spearman), max.col(spearman, 'first'))

col1 <- as.data.frame(rownames(spearman)) #tumor
col2 <- c()
for(i in 1:length(rownames(spearman))){
  col2[i] <- colnames(spearman)[r2[i, 2]]
}
col2 <- as.data.frame(col2)
colsum <- cbind(col1, col2)
colnames(colsum) <- c("Tumor", "Normal")

tumorcluster <- c()
normalcluster <- c()

# ****caf1, t1, caf2, t4 -> gene expression data divided from overlap_ovarian_tfib****#
#processed based on the cluster information. It can be changed by your own cluster.
for(i in 1:length(rownames(spearman))){
  if(colsum[i, 1] %in% colnames(ovary_caf1)){
    tumorcluster[i] <- 1
  }
  else if(colsum[i, 1] %in% colnames(ovary_t2)){
    tumorcluster[i] <- 2
  }
  else if(colsum[i, 1] %in% colnames(ovary_caf2)){
    tumorcluster[i] <- 3
  }
  else if(colsum[i, 1] %in% colnames(ovary_t4)){
    tumorcluster[i] <- 4
  }
}
```
n1, n2, n3, n4 -> gene expression data divided from overlap_ovarian_nfib base on the cluster information. It can be changed by your own cluster.

```r
for(i in 1:length(rownames(spearman))){
  if(colsum[i, 2] %in% colnames(ovary_n1)){
    normalcluster[i] <- 1
  }
  else if(colsum[i, 2] %in% colnames(ovary_n2)){
    normalcluster[i] <- 2
  }
  else if(colsum[i, 2] %in% colnames(ovary_n3)){
    normalcluster[i] <- 3
  }
  else if(colsum[i, 2] %in% colnames(ovary_n4)){
    normalcluster[i] <- 4
  }
}

tumorcluster <- as.data.frame(tumorcluster)
normalcluster <- as.data.frame(normalcluster)
colsum <- cbind(colsum, tumorcluster, normalcluster)

####GGalluvial -> To visualize correlation results between CAF and corresponding precursor normal resident fibroblast.###
library(ggalluvial)
library(ggplot2)

#Remember! It might be changed based on your clustered data!
colsum$tumorcluster[colsum$tumorcluster == 2] = "Cluster_2"
colsum$tumorcluster[colsum$tumorcluster == 1] = "CAF_1"
colsum$tumorcluster[colsum$tumorcluster == 3] = "CAF_2"
colsum$tumorcluster[colsum$tumorcluster == 4] = "Cluster_4"
colsum$normalcluster[colsum$normalcluster == 4] = "Cluster_4"
colsum$normalcluster[colsum$normalcluster == 3] = "Cluster_3"
colsum$normalcluster[colsum$normalcluster == 2] = "Cluster_2"
colsum$normalcluster[colsum$normalcluster == 1] = "Cluster_1"


p <- ggplot(colsum, 
            aes(axis1 = tumorcluster, axis2 = normalcluster))+
  geom_alluvium(aes(fill = tumorcluster), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Cluster_Tumor", "Cluster_Normal"), expand = c(.05, .05)) +
  ggtitle("Correlations between Tumor and Normal clusters of ovarian samples")

a = p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent", color = "white"),
              panel.background = element_rect(fill = "transparent", color = "white")) 
```

From the above, I determined normal cluster 2 as a precursor normal resident fibroblast **(tr-MSCF)**. Now mark it as **ovary_pnrf**.

pnrf cluster might be different with your data

### DEG extraction between CAF and corresponding tr-MSCF. If you have more than one CAF, then perform twice!
``` r
library(limma)

ovary_label = c()
ovary_deg_data <- cbind(ovary_pnrf, ovary_CAF1) #ovary_CAF1, ovary_CAF2...
for (i in 1:ncol(ovary_deg_data)){
  if(i <= ncol(ovary_pnrf)){
    ovary_label[i] <- "PNRF"
  }
  else{
    ovary_label[i] <- "CAF1"
  }
}
ovary_label <- as.matrix(ovary_label)
design <- model.matrix(~0 + ovary_label)
#Before setting colnames, double-check 'design' to validate the column order! It might be opposite in some cases.
colnames(design) <- c("PNRF", "CAF1")

fit <- lmFit(ovary_deg_data,design)
cont <- makeContrasts(CAF1-PNRF,levels=design)
fit.cont <- contrasts.fit(fit,cont)
fit.cont <- eBayes(fit.cont)
res <- topTable(fit.cont,number=Inf)

#Volcanoplot
library(EnhancedVolcano)

vol <- EnhancedVolcano(res, 
                       lab = rownames(res),
                       x = 'logFC',
                       y = 'adj.P.Val',
                       FCcutoff = 0.5849625007211562, #log2(1.5) = 0.5849625007211562, Change if you use different FC.
                       title = "Ovary PNRF vs CAF1")

voldata <- vol$data %>% filter(Sig == "FC_P") ###DEG lists.
```

### Identification of three distinct CAF subtypes -> myCAF, iCAF, apCAF
``` r
apcaf <- c("COL1A1", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "SLPI")
icaf <- c("IL6", "PDGFRA", "CXCL12", "CFD", "DPT", "LMNA", "AGTR1", "HAS1", "CXCL1", "CXCL2", "CCL2") #IL8
mycaf <- c("ACTA2", "TAGLN", "MMP11", "MYL9", "HOPX", "POSTN", "TPM1", "TPM2")

#Three distince CAF genes expression from CAF1, CAF2
#You might use boxplot instead of vioplot if you want
library(vioplot)
ovary_caf1_ap <- ovary_caf1[apcaf,]
ovary_caf1_ap_mean <- as.data.frame(apply(ovary_caf1_ap,2,mean))
ovary_caf2_ap <- ovary_caf2[apcaf,]
ovary_caf2_ap_mean <- as.data.frame(apply(ovary_caf2_ap, 2, mean))

vioplot(unlist(ovary_caf1_ap_mean), unlist(ovary_caf2_ap_mean), 
        col = c("royalblue4", "yellow"),  names=c("1", "2"), ylab = "Signature",
        main="apCAF")

ovary_caf1_i <- ovary_caf1[icaf,]
ovary_caf1_i_mean <- as.data.frame(apply(ovary_caf1_i,2,mean))
ovary_caf2_i <- ovary_caf2[icaf,]
ovary_caf2_i_mean <- as.data.frame(apply(ovary_caf2_i, 2, mean))

vioplot(unlist(ovary_caf1_i_mean), unlist(ovary_caf2_i_mean), 
        col = c("royalblue4", "yellow"),  names=c("1", "2"),
        main="iCAF", ylab="Signature")

ovary_caf1_my <- ovary_caf1[mycaf,]
ovary_caf1_my_mean <- as.data.frame(apply(ovary_caf1_my,2,mean))
ovary_caf2_my <- ovary_caf2[mycaf,]
ovary_caf2_my_mean <- as.data.frame(apply(ovary_caf2_my, 2, mean))

vioplot(unlist(ovary_caf1_my_mean), unlist(ovary_caf2_my_mean),
        col = c("royalblue4", "yellow"),  names=c("1", "2"),
        main="myCAF", ylab="Signature")
```

### Trajectory analysis with monocle R package
``` r
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

biocLite("monocle")
source("http://bioconductor.org/biocLite.R")
biocLite()
library(monocle)

ov_pnrfcaf <- cbind(ov_pnrf, ov_caf1, ov_caf2)

ov_cluster <- data.frame(matrix(NA, nrow = ncol(ov_pnrfcaf), ncol = 4))
colnames(ov_cluster) <- c("cluster", "origin", "day", "cell_type")
rownames(ov_cluster) <- colnames(ov_pnrfcaf)

ov_cluster$origin[rownames(ov_cluster) %in% colnames(ov_caf1)] <- "tumor"
ov_cluster$origin[rownames(ov_cluster) %in% colnames(ov_caf2)] <- "tumor"
ov_cluster$origin[rownames(ov_cluster) %in% colnames(ov_pnrf)] <- "normal"

ov_cluster$cluster[rownames(ov_cluster) %in% colnames(ov_caf1)] <- 1
ov_cluster$cell_type[rownames(ov_cluster) %in% colnames(ov_caf1)] <- "CAF1"

ov_cluster$cluster[rownames(ov_cluster) %in% colnames(ov_caf2)] <- 2
ov_cluster$cell_type[rownames(ov_cluster) %in% colnames(ov_caf2)] <- "CAF2"

ov_cluster$cluster[rownames(ov_cluster) %in% colnames(ov_pnrf)] <- 3
ov_cluster$cell_type[rownames(ov_cluster) %in% colnames(ov_pnrf)] <- "PNRF"

ov_cluster$day[ov_cluster$origin == "tumor"] <- 1
ov_cluster$day[ov_cluster$origin == "normal"] <- 0

#ov_pnrfcaf
HSMM_expr_matrix <- as.matrix(ov_pnrfcaf)
HSMM_sample_sheet <- ov_cluster
#HSMM_gene_annotation <- read.delim("gene_metadata_ovary.txt", row.names = 1)
test <- data.frame(matrix(0, nrow = nrow(ov_pnrfcaf), ncol = 2))
rownames(test) <- rownames(ov_pnrfcaf)
colnames(test) <- c("biotype", "gene_short_name")
test$biotype <- "coding"
test$gene_short_name <- rownames(test)
HSMM_gene_annotation <- test

rownames(HSMM_gene_annotation) <- rownames(ov_pnrfcaf)
HSMM_gene_annotation$gene_short_name <- rownames(ov_pnrfcaf)

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd, featureData = fd)


names(HSMM_expr_matrix) == row.names(HSMM_sample_sheet)

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~origin")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM_myo <- reduceDimension(HSMM, max_components = 2,
                            method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)
plot_cell_trajectory(HSMM_myo, color_by = "cell_type")
```

### Disease and functional enrichment analysis
``` r
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
``` 
Use ov_upgenes and ov_downgenes obtained from DEG analysis. Upgenes and downgenes are divided base on logFC values.
``` r
ov_upgenes <- voldata %>% filter(logFC > 0)
ov_downgenes <- voldata %>% filter(logFC > 0)
ov_genes <- rbind(ov_upgenes, ov_downgenes)

common <- subset(ov_genes, select = logFC)

t <- rownames(common)
et <- bitr(t, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")

edo <- enrichDGN(et$ENTREZID) #DisGeNet
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
up <- setNames(unlist(common, use.names=F),rownames(common))

cnetplot(edox, foldChange=up, circular = TRUE, colorEdge = TRUE)
```

### GSVA and clinical significance of CAF signature genes
```r
library(GSEABase)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(GSVAdata)
library(dplyr)

setwd("/Users/ryanmachine/")
#You should download your own TCGA datasets from https://gdac.broadinstitute.org, clinical data from https://www.cbioportal.org/
tcga_ovary <- read.delim("Downloads/single_cell_material/TCGA_data/ovary/ovary.rna.dat.txt", row.names = 1)
upgenes <- rownames(ov_upgenes)
downgenes <- rownames(ov_downgenes)

#########GSVA with significant gene sets (Upregulated genes in CAF, downregulated genes in CAF obtained by DEG analysis)
upgenes <- list(upgenes)
downgenes <- list(downgenes)
gsva_up=gsva(as.matrix(tcga_ovary), upgenes, method="gsva", kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
gsva_df_up=data.frame(t(gsva_up), stringsAsFactors = F, check.names = F)

gsva_down=gsva(as.matrix(tcga_ovary), downgenes, method="gsva", kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
gsva_df_down=data.frame(t(gsva_down), stringsAsFactors = F, check.names = F)

ovary_gsva <- cbind(gsva_df_up, gsva_df_down)
colnames(ovary_gsva) <- c("Up", "Down")

#Use your own clinical data: I already pre-processed overall survival data of ovary cancer patient's data which is linked to TCGA data
ovary_overall <- read.csv("Downloads/single_cell_material/TCGA_data/ovary/ov_survival_overall.csv", row.names =1)
ovary_gsva_ov <- ovary_gsva[row.names(ovary_overall),]

ovary_ov <- data.frame(matrix(0, ncol = 5, nrow = nrow(ovary_gsva_ov)))
ovary_ov$X1 <- ovary_overall$Overall.Survival..Months.
ovary_ov$X2 <- ovary_overall$Overall.Survival.Status
ovary_ov$X3 <- ovary_gsva_ov$Up
ovary_ov$X4 <- ovary_gsva_ov$Down
ovary_ov$X5 <- 0
colnames(ovary_ov) <- c("Overall.Survival..Months.", "Overall.Survival.Status", "Up", "Down", "CAF")
rownames(ovary_ov) <- rownames(ovary_overall)
```



### Survival analysis
```r
library(survival)
library(survminer)
```
Calculate the median of GSVA score of upregulated genes to separate samples between two groups.
```r
upmed <- median(ovary_ov$Up)
g1 <- ovary_ov %>% filter(Up > upmed)
g2 <- ovary_ov %>% filter(!(Up > upmed))
ovary_ov$CAF <- 0
ovary_ov$CAF[rownames(ovary_ov) %in% rownames(g1)] <- 1
ovary_ov$CAF[rownames(ovary_ov) %in% rownames(g2)] <- 2

ovary_ov$CAF <- factor(ovary_ov$CAF)

fit2 <- survfit(Surv(Overall.Survival..Months., 
                     as.numeric(Overall.Survival.Status))~CAF, data = ovary_ov)
ggsurvplot(fit2, legend.title="overall_Survival",
           conf.int = T, pval = T)
```

### Validation 1 - Slingshot trajectory analysis
To perform slingshot, you should use **raw count matrix**!
``` r
library(slingShot)
crc2 <- read.delim("Downloads/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt", row.names = 1)

crc2_mscf <- read.csv("GSVA_Fibroblasts/colon/SMC_nfib_cluster2.csv", row.names = 1)
crc2_rf <- read.csv("GSVA_Fibroblasts/SMC_CRC/Normal/SMC_nfib_cluster1.csv", row.names = 1)
crc2_noncaf <- read.csv("GSVA_Fibroblasts/SMC_CRC/Tumor/SMC_tfib_cluster2.csv", row.names = 1)
crc2_caf <- read.csv("GSVA_Fibroblasts/colon/SMC_tfib_cluster1.csv", row.names = 1)

crc2 <- c(colnames(crc2_mscf), colnames(crc2_rf), colnames(crc2_noncaf), colnames(crc2_caf))

sce <- SingleCellExperiment(assays = List(counts = as.matrix(crc2)))
rm(crc2)

geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

sce <- sce[,crc2]

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

library(mclust, quietly = TRUE)

cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

colData(sce)$kmeans[colnames(sce) %in% colnames(crc2_mscf)] <- 1
colData(sce)$kmeans[colnames(sce) %in% colnames(crc2_rf)] <- 2
colData(sce)$kmeans[colnames(sce) %in% colnames(crc2_noncaf)] <- 3
colData(sce)$kmeans[colnames(sce) %in% colnames(crc2_caf)] <- 4
cl2 <- colData(sce)$kmeans
names(cl2) <- rownames(colData(sce))

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce_shot <- slingshot(sce, clusterLabels = "kmeans", reducedDim = "PCA")

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_shot$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce_shot)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_shot), lwd=2, col='black')

plot(reducedDims(sce_shot)$PCA, col = brewer.pal(9,'Set1')[sce_shot$kmeans], pch=16, asp = 1)
lines(SlingshotDataSet(sce_shot), lwd=2, type = 'lineages', col = 'black')

lin1 <- getLineages(rd1, cl2, start.clus = '1', end.clus = '4')
plot(rd1, col = c("yellow", "gray75", "darkgreen", "mediumpurple")[cl2], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black', show.constraints = TRUE)

temp <- rd1
temp <- as.data.frame(temp)
temp$pseudotime <- sce_shot$slingPseudotime_1
temp$cluster <- sce_shot$kmeans

#setwd("/Users/ryanmachine/Dropbox/Single_cell/Manuscript/submission/CTM_letter_revision1/slingshot/")
#write.csv(temp, "crc2_data.csv")
#setwd("/Users/ryanmachine/")
a <- ggplot(temp, aes(x = pseudotime, y = PC1)) + geom_point(colour = c("yellow", "gray75", "darkgreen", "mediumpurple")[temp$cluster])
a = a + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent", color = "white"),
              panel.background = element_rect(fill = "transparent", color = "white")) 
```

### Validation 2 - RNA velocity
``` r
library(SeuratDisk)
library(SeuratWrappers)
```
You need **loom** format files generated from **CellRanger** and **velocyto** using **FASTQ** files to perform RNA velocity
``` r
setwd("//home/miruware/data/ryanchung/jupyter_codes/velocity/")
filepath <- "//home/miruware/data/ryanchung/jupyter_codes/velocity/CRC/"
files <- list.files(filepath, pattern=".loom")

mscf <- read.csv("pancrc_nfib_cluster1.csv", row.names = 1)
caf1 <- read.csv("crc1_caf_subgroup1.csv", row.names = 1)
caf2 <- read.csv("crc1_caf_subgroup2.csv", row.names = 1)

caf <- c(colnames(mscf), colnames(caf1), colnames(caf2))

setwd("//home/miruware/data/ryanchung/jupyter_codes/velocity/CRC/")
i<- 1
for (file in files){
  test <- ReadVelocity(file = file)
  assign(paste0("CRC_", i), as.Seurat(test))
  i <- i+1
}

bm <- merge(x = CRC_1, y = c(CRC_2, CRC_3, CRC_4, CRC_5, CRC_6, CRC_7, CRC_8, CRC_9, CRC_10, CRC_11, CRC_12, CRC_13, 
                             CRC_14,CRC_15, CRC_16, CRC_17, CRC_18, CRC_19, CRC_20, CRC_21, CRC_22, CRC_23, CRC_24, 
                             CRC_25, CRC_26, CRC_27), merge.data = TRUE)

temp <- colnames(bm)
temp2 <- gsub(temp, pattern = "hg19:", replace = "")
temp2 <- gsub(temp2, pattern = "x", replace = "")
bm <- RenameCells(bm, new.names = temp2)


bm[['RNA']] <- bm[['spliced']]
bm[['percent.mt']] <- PercentageFeatureSet(bm, pattern = "^MT-")
bm_msc <- subset(bm, subset = nCount_spliced > 401 & nFeature_spliced > 201 & nFeature_spliced < 6000 & percent.mt < 25)
#bm_msc <- bm_msc[,caf] <- If you only want to investigate the velocity of mscf, caf1, caf2 
bm_msc <- SCTransform(object = bm_msc)
bm_msc <- RunPCA(object = bm_msc, verbose = FALSE)
bm_msc <- RunUMAP(bm_msc, dims = 1:20)
bm_msc <- FindNeighbors(object = bm_msc, dims = 1:20)
bm_msc <- FindClusters(object = bm_msc)


temp <- bm_msc
temp$seurat_clusters <- as.character(temp$seurat_clusters)
temp$seurat_clusters[colnames(temp) %in% colnames(mscf)] <- "tr-MSCF"
temp$seurat_clusters[colnames(temp) %in% colnames(caf1)] <- "iCAF"
temp$seurat_clusters[colnames(temp) %in% colnames(caf2)] <- "myCAF"


setwd("//home/miruware/data/ryanchung/jupyter_codes/velocity/")
DefaultAssay(temp) <- "RNA"
SaveH5Seurat(temp, filename = "Pan.h5Seurat")
Convert("Pan.h5Seurat", dest = "h5ad")
```
### Next step using python -> scVelo.ipynb
