# Load R packages/libraries
suppressMessages(require(Seurat))
suppressMessages(require(scater))
suppressMessages(require(Matrix))
library(scran)
suppressMessages(require(tidyverse))
suppressMessages(require(igraph))
suppressMessages(require(clustree))
library(org.Dr.eg.db)
library("pheatmap")
library("RColorBrewer")
library(gridExtra)

source("vlnplot.R")
source("multiplot.R")

setwd("MN")

# Create Seurat objects, store metadata and QC measures####

# Read raw counts and metadata. Islet plates are 174 and 358
counts1 <- read.table("358.354.353.merge.rpkmforgenes_counts.csv", header=TRUE, sep=',', row.names='gene')
counts2 <- read.table("174.merge.rpkmforgenes_counts.csv", header=TRUE, sep=',', row.names='gene')


dim(counts1)
#[1] 32612  1152
dim(counts2)
#[1] 32612   384

# Add gene names to IDs
genes <- rownames(counts1)
table(genes == rownames(counts2))
#  TRUE 
#32612 
symb <- mapIds(org.Dr.eg.db, keys=genes, keytype="ENSEMBL", column="SYMBOL")
m <- match(rownames(counts1),names(symb)) 
newnames <- apply(cbind(as.vector(symb)[m],rownames(counts1)),1,paste,collapse=":")
rownames(counts1) <- newnames

m <- match(rownames(counts2),names(symb)) # To make sure that gene ids match
newnames <- apply(cbind(as.vector(symb)[m],rownames(counts2)),1,paste,collapse=":")
rownames(counts2) <- newnames

# Read metadata
info1 <- read.table("358.354.353.sample_info.csv", header=TRUE, sep=',', row.names='SM')
info2 <- read.table("174.sample_info.csv", header=TRUE, sep=',', row.names='SM')

# Create Seurat objects, take only genes expressed in more than 3 cells
seu_c1 <- CreateSeuratObject(counts1, project = "SS2_18_c", meta.data = info1, min.cells =3)
seu_c2 <- CreateSeuratObject(counts2, project = "SS2_18_c", meta.data = info2, min.cells =3)

# Read and store mapping info
star1 <- read.table("358.354.353.multiqc_star.txt", header=TRUE, sep='\t', row.names='Sample')
star2 <- read.table("174.multiqc_star.txt", header=TRUE, sep='\t', row.names='Sample')
seu_c1 <- AddMetaData(seu_c1, star1$total_reads, col.name = "total_reads")
seu_c1 <- AddMetaData(seu_c1, star1$uniquely_mapped, col.name = "uniquely_mapped")
seu_c1 <- AddMetaData(seu_c1, star1$uniquely_mapped_percent, col.name = "uniquely_mapped_percent")
seu_c1 <- AddMetaData(seu_c1, star1$multimapped_percent, col.name = "multimapped_percent")
seu_c2 <- AddMetaData(seu_c2, star2$total_reads, col.name = "total_reads")
seu_c2 <- AddMetaData(seu_c2, star2$uniquely_mapped, col.name = "uniquely_mapped")
seu_c2 <- AddMetaData(seu_c2, star2$uniquely_mapped_percent, col.name = "uniquely_mapped_percent")
seu_c2 <- AddMetaData(seu_c2, star2$multimapped_percent, col.name = "multimapped_percent")

# Select plate 358 ==> Remove plates 353 and 354 (V2a)
selected <- WhichCells(seu_c1, expression = plate == "p358" )
seu_c1 <- subset(seu_c1, cells = selected)

# Merged plates#### -This section can be skipped for reproducing paper analysis
# This section demonstrates why data needed integration.
seu_c <- merge(seu_c1, seu_c2)

# Quality control, filtering and normalizing
# Calculate proportion of mitochondrial genes
mt.genes <- read.table("mt_genes.txt")
mt.genes <- lapply(mt.genes, as.character)$V1   #Convert data frame to string array...
m <- match(mt.genes,names(symb)) # To make sure that gene ids match
mt.genes <- apply(cbind(as.vector(symb)[m],mt.genes),1,paste,collapse=":")
C<-GetAssayData(object = seu_c, slot = "counts")
percent.mt <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
seu_c <- AddMetaData(seu_c, percent.mt, col.name = "percent.mito")

# Calculate proportion of spikein
spikein <- rownames(seu_c)[grep(":ERCC-",rownames(seu_c))]
percent.spike <- colSums(C[spikein,])/Matrix::colSums(C)*100
seu_c <- AddMetaData(seu_c, percent.spike, col.name = "percent.spike")

# Calculate proportion of ribosomal proteins
rib.genes <- read.table("rpls_genes.txt")
rib.genes <- lapply(rib.genes, as.character)$V1   #Convert data frame to string array...
m <- match(rib.genes,names(symb)) # To make sure that gene ids match
rib.genes <- apply(cbind(as.vector(symb)[m],rib.genes),1,paste,collapse=":")
percent.rib <- colSums(C[rib.genes,])/Matrix::colSums(C)*100
seu_c <- AddMetaData(seu_c, percent.rib, col.name = "percent.ribo")

# Calculate proportion of GFP (pEGFP-N1)
gfp <- rownames(seu_c)[grep(":pEGFP-N1",rownames(seu_c))]
percent.gfp <- C[gfp,]/Matrix::colSums(C)*100
seu_c <- AddMetaData(seu_c, percent.gfp, col.name = "percent.gfp")

# Plot QC and filtering

# Number of reads per cell and number of features / detected genes
VlnPlot(object = seu_c, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")

# Scatter plot
FeatureScatter(seu_c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")

# Percent mitochondrial and ribosomal and spikein and pgfp
VlnPlot(seu_c, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
VlnPlot(seu_c, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)

# Filter out cells with high (>20%) spike-in
seu_c.filt <- subset(seu_c, subset = percent.spike < 20 )
Idents(seu_c) <- 'plate'
Idents(seu_c.filt) <- 'plate'

table(Idents(seu_c))
#p358 p174 
#384  384 

table(Idents(seu_c.filt))
#p358 p174 
#267  383 

# Plot the proportion of mitochondrial counts against some of the other QC metrics.
p1 <- FeatureScatter(seu_c.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
p2 <- FeatureScatter(seu_c.filt, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by="plate")
p3 <- FeatureScatter(seu_c.filt, feature1 = "nFeature_RNA", feature2 = "percent.mito", group.by="plate")
multiplot(p1, p2, p3, cols = 1)

# Mitochondrial filtering
# Those with higher % than that are not likely to contain much biological signal
seu_c.filt <- subset(seu_c.filt,  subset = percent.mito < 30)
table(Idents(seu_c.filt))
# p358 p174 
#247  367

#plot after filtering
VlnPlot(seu_c.filt, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
VlnPlot(seu_c.filt, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)

# Number of reads per cell and number of features / detected genes
VlnPlot(object = seu_c.filt, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")

# Remove cells with very high gene detection (possibly doublets)
seu_c.filt <- subset(seu_c.filt, subset = nFeature_RNA < 10000)
table(Idents(seu_c.filt))
#p358 p174 
#245  367 

# Filter the cells with low gene detection (low quality libraries) 
seu_c.filt <- subset(seu_c.filt, subset = nFeature_RNA > 1000)
table(Idents(seu_c.filt))
# p358 p174 
#244  362 

#plot after filtering
VlnPlot(seu_c.filt, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
VlnPlot(seu_c.filt, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)

# Scatter and violin plots of nCounts and nFeature
FeatureScatter(seu_c.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")

# Regress out number of detected genes, number of reads, percent_ribo
#NOTE: typos here. this code is actually only regressing nFeature
seu_c.filt <- ScaleData(seu_c.filt, vars.to.regress = c("nFeature_RNA", "nCounts_RNA", "percent_ribo"))

# Clustering

# Normalize
scale.factor_c <- mean(colSums(GetAssayData(seu_c.filt, slot = "counts")))
scale.factor_c
#[1] 435291
seu_c.filt <- NormalizeData(object = seu_c.filt, normalization.method = "LogNormalize", scale.factor = scale.factor_c)

# Feature selection
seu_c.filt <- FindVariableFeatures(seu_c.filt, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(seu_c.filt), 10)
# [1] "tph2:ENSDARG00000057239"    "ddc:ENSDARG00000016494"     "hspa1b:ENSDARG00000056210" 
#[4] "urp2:ENSDARG00000067781"    "crhbp:ENSDARG00000024831"   "her15.2:ENSDARG00000054560"
#[7] "otpa:ENSDARG00000014201"    "socs3a:ENSDARG00000025428"  "hspb1:ENSDARG00000041065"  
#[10] "calb2a:ENSDARG00000041062" 

vplot_c <- VariableFeaturePlot(seu_c.filt)
top10_c <- head(VariableFeatures(seu_c.filt), 10)
labels <- c("tph2","ddc","hspa1b","urp2","crhbp","her15.2","otpa","socs3a","hspb1","calb2a")
LabelPoints(plot = vplot_c, points = top10_c, repel = TRUE, xnudge=0, ynudge=0, label=labels)

#clustering
# PCA on variable genes
# Default 50 PCs
seu_c.filt <- RunPCA(seu_c.filt, features = VariableFeatures(seu_c.filt), do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# Plot PCAs
PCAPlot(object = seu_c.filt, dims=c(1,2),group.by="plate")


ElbowPlot(seu_c.filt)
seu_c.filt <- JackStraw(seu_c.filt)
seu_c.filt <- ScoreJackStraw(seu_c.filt, dims = 1:20)
pcs <- which(JS(object = seu_c.filt[['pca']], slot = 'overall')[, 2] < 1e-3)
pcs
# [1]  1  2  3  4  5  6  7  8 12 13 15 16 17 18 19

seu_c.filt <- FindNeighbors(seu_c.filt, reduction = "pca", dims = 1:10)
seu_c.filt <- FindClusters(seu_c.filt, resolution = 0.6)
table(Idents(seu_c.filt))

# 10 dim 0.6 res:
#  0   1   2   3   4   5   6 
#186  82  75  74  70  62  57 

# Check clusters over resolutions
res <- c(0.4,0.5,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4)
tmp=seu_c.filt
tmp <- FindClusters(tmp, resolution = res, print.output = 0)
clustree(tmp)

# Run UMAP and tSNE
seu_c.filt <- RunUMAP(seu_c.filt, dims = 1:10)
seu_c.filt <- RunTSNE(seu_c.filt, dims = 1:10)
DimPlot(seu_c.filt, reduction = "umap", group.by = "plate")
DimPlot(seu_c.filt, reduction = "umap")

#here you can see that clusters are segregated by plate. Not good, probable batch effect. Integrate plates!


# Integrate plates####

#normalize and find var features on each plate then integrate
scale.factor_c1 <- mean(colSums(GetAssayData(seu_c1, slot = "counts")))
scale.factor_c1
#[1] 299428.4
scale.factor_c2 <- mean(colSums(GetAssayData(seu_c2, slot = "counts")))
scale.factor_c2
#[1] 431300.9

seu_c1 <- NormalizeData(seu_c1, normalization.method = "LogNormalize", scale.factor = scale.factor_c2, verbose = FALSE)
seu_c2 <- NormalizeData(seu_c2, normalization.method = "LogNormalize", scale.factor = scale.factor_c2, verbose = FALSE)

#find variable features
seu_c1 <- FindVariableFeatures(seu_c1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
head(VariableFeatures(seu_c1), 10)
#[1] "tph2:ENSDARG00000057239"     "urp2:ENSDARG00000067781"     "ddc:ENSDARG00000016494"      "tubb5:ENSDARG00000037997"    "nfixb:ENSDARG00000061836"   
#[6] "adcyap1b:ENSDARG00000027740" "otpa:ENSDARG00000014201"     "hspb1:ENSDARG00000041065"    "fev:ENSDARG00000009242"      "calb2a:ENSDARG00000041062"  

seu_c2 <- FindVariableFeatures(seu_c2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
head(VariableFeatures(seu_c2), 10)
#[1] "tph2:ENSDARG00000057239"    "ddc:ENSDARG00000016494"     "mbpa:ENSDARG00000036186"    "crhbp:ENSDARG00000024831"   "hspa1b:ENSDARG00000056210" 
#[6] "urp2:ENSDARG00000067781"    "calb2a:ENSDARG00000041062"  "her15.2:ENSDARG00000054560" "nppc:ENSDARG00000068126"    "tac1:ENSDARG00000014490"   

#integrate datasets
all.genes<-row.names(counts1) 
list<-list(seu_c1, seu_c2)
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:30)
# Computing 2000 integration features
# Scaling features for provided objects
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=00s  
# Finding all pairwise anchors
# |                                                  | 0 % ~calculating  Running CCA
# Merging objects
# Finding neighborhoods
# Finding anchors
# Found 1135 anchors
# Filtering anchors
# Retained 1085 anchors
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
p174.p358_i <- IntegrateData(anchorset = anchors, dims = 1:30, features.to.integrate = all.genes)
# Merging dataset 2 into 1
# Extracting anchors for merged samples
# Finding integration vectors
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Integrating data
head(VariableFeatures(p174.p358_i), 10)
#[1] "tph2:ENSDARG00000057239"    "ddc:ENSDARG00000016494"     "urp2:ENSDARG00000067781"    "calb2a:ENSDARG00000041062"  "hspa1b:ENSDARG00000056210" 
#[6] "fev:ENSDARG00000009242"     "socs3a:ENSDARG00000025428"  "tac1:ENSDARG00000014490"    "nppc:ENSDARG00000068126"    "cited4b:ENSDARG00000101009"

DefaultAssay(p174.p358_i) <- "integrated"

# Quality control, filtering and normalizing####
# Calculate proportion of mitochondrial genes
setwd("MN/QC")
mt.genes <- read.table("mt_genes.txt")
mt.genes <- lapply(mt.genes, as.character)$V1   #Convert data frame to string array...
m <- match(mt.genes,names(symb)) # To make sure that gene ids match
mt.genes <- apply(cbind(as.vector(symb)[m],mt.genes),1,paste,collapse=":")
C<-GetAssayData(object = p174.p358_i, assay = "RNA", slot = "counts")
percent.mt <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
p174.p358_i <- AddMetaData(p174.p358_i, percent.mt, col.name = "percent.mito")

# Calculate proportion of spikein
spikein <- rownames(p174.p358_i)[grep(":ERCC-",rownames(p174.p358_i))]
percent.spike <- colSums(C[spikein,])/Matrix::colSums(C)*100
p174.p358_i <- AddMetaData(p174.p358_i, percent.spike, col.name = "percent.spike")

# Calculate proportion of ribosomal proteins
rib.genes <- read.table("rpls_genes.txt")
rib.genes <- lapply(rib.genes, as.character)$V1   #Convert data frame to string array...
m <- match(rib.genes,names(symb)) # To make sure that gene ids match
rib.genes <- apply(cbind(as.vector(symb)[m],rib.genes),1,paste,collapse=":")
percent.rib <- colSums(C[rib.genes,])/Matrix::colSums(C)*100
p174.p358_i <- AddMetaData(p174.p358_i, percent.rib, col.name = "percent.ribo")

# Calculate proportion of GFP (pEGFP-N1)
gfp <- rownames(p174.p358_i)[grep(":pEGFP-N1",rownames(p174.p358_i))]
percent.gfp <- C[gfp,]/Matrix::colSums(C)*100
p174.p358_i <- AddMetaData(p174.p358_i, percent.gfp, col.name = "percent.gfp")

# Plot QC and filtering

# Number of reads per cell and number of features / detected genes
png('1.VlnPlot.nCount-nFeature_byPlate.png', width=480, height=250)
VlnPlot(object = p174.p358_i, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")
dev.off()

pdf('1.VlnPlot.nCount-nFeature_byPlate.pdf')
hist(p174.p358_i@meta.data[["nFeature_RNA"]])
dev.off()

# Scatter plot
png('2.Scatter.nCount-nFeature.rpkm.png', width=480, height=250)
FeatureScatter(p174.p358_i, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
dev.off()

# Percent mitochondrial and ribosomal and spikein and pgfp
png('3.VlnPlot_pctmito-pctribo.png')
VlnPlot(p174.p358_i, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
dev.off()

png('3b.VlnPlot_pctspike-pctgfp.png')
VlnPlot(p174.p358_i, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)
dev.off()


# Filter out cells with high (>20%) spike-in
p174.p358_i.filt <- subset(p174.p358_i, subset = percent.spike < 20)

Idents(p174.p358_i) <- 'plate'
Idents(p174.p358_i.filt) <- 'plate'
# How many cells are left?
table(Idents(p174.p358_i))
#p358 p174 
#384  384 

table(Idents(p174.p358_i.filt))
#p358 p174 
#268  383 

# Plot the proportion of mitochondrial counts against some of the other QC metrics.
png('3b.Scatter.mito_count.mito_feature.p174.358.png')
p1 <- FeatureScatter(p174.p358_i, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
p2 <- FeatureScatter(p174.p358_i, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by="plate")
p3 <- FeatureScatter(p174.p358_i, feature1 = "nFeature_RNA", feature2 = "percent.mito", group.by="plate")
grid.arrange(p1, p2, p3, ncol = 1)
dev.off()

# Mitochondrial filtering
p174.p358_i.filt <- subset(p174.p358_i.filt,  subset = percent.mito < 30)
table(Idents(p174.p358_i.filt))
# p358 p174 
#247  367

# Percent mitochondrial and ribosomal and spikein and pgfp
png('4.VlnPlot_pctmito-pctribo.cellfilt.png')
VlnPlot(p174.p358_i.filt, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
dev.off()

png('4b.VlnPlot_pctspike-pctgfp.counts.cellfilt.png')
VlnPlot(p174.p358_i.filt, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)
dev.off()

#cells with 0  mito and 0 ribo disappeared=probably the high spike-ins, empty slots

# Scatter plot
png('5.Scatter.count_feature_by plate.png')
p1 <- FeatureScatter(p174.p358_i.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
p2 <- FeatureScatter(p174.p358_i.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
grid.arrange(p1, p2, ncol = 1)
dev.off()

# Number of reads per cell and number of features / detected genes
png('5.VlnPlot.count_feature_by plate.png')
VlnPlot(object = p174.p358_i.filt, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")
dev.off()

# Remove cells with very high gene detection (possibly doublets)
p174.p358_i.filt <- subset(p174.p358_i.filt, subset = nFeature_RNA < 10000)
length(selected)
table(Idents(p174.p358_i.filt))
#p358 p174 
#245  367 

# Filter the cells with low gene detection (low quality libraries) 
p174.p358_i.filt <- subset(p174.p358_i.filt, subset = nFeature_RNA > 1000)
length(selected)
table(Idents(p174.p358_i.filt))
# p358 p174 
#244  362 

# Plots after all cell filtering
# Percent mitochondrial and ribosomal and spikein and pgfp
# after all cell filtering
png('6.VlnPlot_pctmito-pctribo.cellfilt.png')
VlnPlot(p174.p358_i.filt, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
dev.off()

png('6b.VlnPlot_pctspike-pctgfp.cellfilt.png')
VlnPlot(p174.p358_i.filt, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)
dev.off()

# Scatter and violin plots of nCounts and nFeature
png('7.Scatter_nCount-nFeature.byPlate.cellfilt.png', width=700, height=600)
FeatureScatter(p174.p358_i.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
dev.off()

png('7b.VlnPlot_nCount-nFeature.byPlate.cellfilt.png', width=700, height=600)
VlnPlot(object = p174.p358_i.filt, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")
dev.off()


# Clustering####
setwd("MN")
# Feature selection
seu_c.filt2 <- FindVariableFeatures(seu_c.filt2, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(p174.p358_i.filt), 10)
#  [1] "tph2:ENSDARG00000057239"    "ddc:ENSDARG00000016494"     "urp2:ENSDARG00000067781"    "calb2a:ENSDARG00000041062"  "hspa1b:ENSDARG00000056210" 
#[6] "fev:ENSDARG00000009242"     "socs3a:ENSDARG00000025428"  "tac1:ENSDARG00000014490"    "nppc:ENSDARG00000068126"    "cited4b:ENSDARG00000101009"

# Regress out number of detected genes, number of reads, percent_ribo 
# still with the typos, actually regressing only bt nfeature
p174.p358_i.filt <- ScaleData(p174.p358_i.filt, vars.to.regress = c("nFeature_RNA", "nCounts_RNA", "percent_ribo"))


# PCA on variable genes
# Default 50 PCs
p174.p358_i.filt <- RunPCA(p174.p358_i.filt, features = VariableFeatures(p174.p358_i.filt), do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# PC_ 1 
# Positive:  pkd1l2a:ENSDARG00000105344, urp1:ENSDARG00000093493, myo3b:ENSDARG00000006892, pkd2l1:ENSDARG00000022503, c2cd4a:ENSDARG00000061416, tal2:ENSDARG00000042041, urp2:ENSDARG00000067781, slc6a1b:ENSDARG00000039647, espn:ENSDARG00000076414, gad1b:ENSDARG00000027419 
# gad1a:ENSDARG00000093411, si:ch1073-70f20.1:ENSDARG00000074930, ntn1b:ENSDARG00000022531, foxa:ENSDARG00000087094, gad2:ENSDARG00000015537, tal1:ENSDARG00000019930, esm1:ENSDARG00000104370, nppc:ENSDARG00000068126, NA:ENSDARG00000104929, scinlb:ENSDARG00000058348 
# gpx3:ENSDARG00000043342, ldlrad2:ENSDARG00000090297, bag3:ENSDARG00000039486, rhpn2:ENSDARG00000014577, pard3ba:ENSDARG00000092671, rfx4:ENSDARG00000026395, fosb:ENSDARG00000055751, si:dkey-30j16.3:ENSDARG00000037587, eps8l2:ENSDARG00000058108, fat1b:ENSDARG00000019063 
# Negative:  cd99l2:ENSDARG00000056722, chata:ENSDARG00000015854, neurod1:ENSDARG00000019566, prph:ENSDARG00000028306, nfixb:ENSDARG00000061836, slc18a3a:ENSDARG00000006356, nhlh2:ENSDARG00000025495, ebf3a:ENSDARG00000100244, ndufa4:ENSDARG00000056108, alcama:ENSDARG00000026531 
# mex3b:ENSDARG00000058369, adcyap1b:ENSDARG00000027740, nfixa:ENSDARG00000043226, gap43:ENSDARG00000099744, slc5a7a:ENSDARG00000074860, mafbb:ENSDARG00000070542, hapln1a:ENSDARG00000089769, sox11a:ENSDARG00000077811, si:ch211-232b12.5:ENSDARG00000087186, zgc:65851:ENSDARG00000012281 
# scn4ba:ENSDARG00000099031, tubb5:ENSDARG00000037997, glrba:ENSDARG00000052782, lhx4:ENSDARG00000039458, si:dkey-7j14.5:ENSDARG00000097528, hunk:ENSDARG00000091317, nat8l:ENSDARG00000077256, LOC100331987:ENSDARG00000100594, nxph1:ENSDARG00000033447, fsta:ENSDARG00000052846 
# PC_ 2 
# Positive:  fev:ENSDARG00000009242, tph2:ENSDARG00000057239, ddc:ENSDARG00000016494, slc6a4a:ENSDARG00000061165, ucn3l:ENSDARG00000087241, pth3r:ENSDARG00000018418, ntsr1:ENSDARG00000077577, lmx1bb:ENSDARG00000068365, calb2a:ENSDARG00000041062, prkg1b:ENSDARG00000031702 
# gch1:ENSDARG00000070453, npr1b:ENSDARG00000018750, mt2:ENSDARG00000041623, rxfp3.3b:ENSDARG00000059348, htr1d:ENSDARG00000054124, pcsk1:ENSDARG00000002600, prdx1:ENSDARG00000058734, gstp1:ENSDARG00000104068, cygb1:ENSDARG00000099371, dpp6b:ENSDARG00000024744 
# plcxd3:ENSDARG00000054794, si:dkeyp-72h1.1:ENSDARG00000095347, cpne5a:ENSDARG00000061466, NA:ENSDARG00000114226, kcnip3b:ENSDARG00000017880, qdpra:ENSDARG00000040190, zfhx3:ENSDARG00000103057, si:dkey-27j5.5:ENSDARG00000026383, rimbp2:ENSDARG00000001154, col14a1a:ENSDARG00000005762 
# Negative:  isl1:ENSDARG00000004023, vim:ENSDARG00000010008, si:dkey-56m19.5:ENSDARG00000068432, mafbb:ENSDARG00000070542, neurod1:ENSDARG00000019566, NA:pEGFP-N1, nfixa:ENSDARG00000043226, nhlh2:ENSDARG00000025495, chata:ENSDARG00000015854, myo1b:ENSDARG00000024694 
# fgfr4:ENSDARG00000069105, ebf3a:ENSDARG00000100244, atp1a3b:ENSDARG00000104139, hapln1a:ENSDARG00000089769, nfixb:ENSDARG00000061836, pkd2l1:ENSDARG00000022503, sall1a:ENSDARG00000074319, plxna4:ENSDARG00000019328, urp2:ENSDARG00000067781, ldlrad2:ENSDARG00000090297 
# alcama:ENSDARG00000026531, si:ch1073-70f20.1:ENSDARG00000074930, ccdc88c:ENSDARG00000053713, slc18a3a:ENSDARG00000006356, pkd1l2a:ENSDARG00000105344, necab2:ENSDARG00000056745, foxa:ENSDARG00000087094, btbd6b:ENSDARG00000032369, scinlb:ENSDARG00000058348, esm1:ENSDARG00000104370 
# PC_ 3 
# Positive:  ndufa4:ENSDARG00000056108, esrrga:ENSDARG00000004861, bmp16:ENSDARG00000103679, tac1:ENSDARG00000014490, nefma:ENSDARG00000021351, rspo2:ENSDARG00000079570, scn4ba:ENSDARG00000099031, aclya:ENSDARG00000099079, neflb:ENSDARG00000012426, slc5a7a:ENSDARG00000074860 
# pcdh19:ENSDARG00000034344, glrba:ENSDARG00000052782, LOC100331987:ENSDARG00000100594, pvalb6:ENSDARG00000009311, grin1b:ENSDARG00000025728, slc18a3a:ENSDARG00000006356, fsta:ENSDARG00000052846, scn4bb:ENSDARG00000060319, chata:ENSDARG00000015854, cbx7a:ENSDARG00000038025 
# prph:ENSDARG00000028306, zgc:65851:ENSDARG00000012281, pdyn:ENSDARG00000087798, npas1:ENSDARG00000015876, si:dkey-33c12.3:ENSDARG00000057881, atp2a3:ENSDARG00000060978, chka:ENSDARG00000041078, ak5:ENSDARG00000012555, ptprna:ENSDARG00000058646, atp1a3b:ENSDARG00000104139 
# Negative:  hapln1a:ENSDARG00000089769, gja1b:ENSDARG00000041799, notch1a:ENSDARG00000103554, her4.5:ENSDARG00000056729, selenop:ENSDARG00000093549, cldn7a:ENSDARG00000036376, hspb15:ENSDARG00000078411, nfixb:ENSDARG00000061836, her4.2:ENSDARG00000094426, sulf2a:ENSDARG00000018423 
# perp:ENSDARG00000063572, zfp36l1a:ENSDARG00000016154, ebf3a:ENSDARG00000100244, neurod4:ENSDARG00000003469, sox19a:ENSDARG00000010770, mgst1.1:ENSDARG00000032618, ncam1a:ENSDARG00000056181, si:ch1073-303k11.2:ENSDARG00000088247, dla:ENSDARG00000010791, her15.2:ENSDARG00000054560 
# nhlh2:ENSDARG00000025495, stm:ENSDARG00000035694, tubb5:ENSDARG00000037997, tmod4:ENSDARG00000020890, her15.1:ENSDARG00000054562, si:ch211-251b21.1:ENSDARG00000007275, ebf2:ENSDARG00000042525, abtb2b:ENSDARG00000062000, myo10l3:ENSDARG00000074143, dlb:ENSDARG00000004232 
# PC_ 4 
# Positive:  zfp36l1a:ENSDARG00000016154, selenop:ENSDARG00000093549, si:ch1073-303k11.2:ENSDARG00000088247, gja1b:ENSDARG00000041799, perp:ENSDARG00000063572, npas1:ENSDARG00000015876, mgst1.1:ENSDARG00000032618, cldn7a:ENSDARG00000036376, si:ch211-251b21.1:ENSDARG00000007275, atp1b4:ENSDARG00000053262 
# hspb15:ENSDARG00000078411, esrrga:ENSDARG00000004861, swap70b:ENSDARG00000057286, myo1f:ENSDARG00000078734, atp2a3:ENSDARG00000060978, bmp16:ENSDARG00000103679, cnn2:ENSDARG00000035858, snap23.1:ENSDARG00000012874, slc5a7a:ENSDARG00000074860, LOC100334443:ENSDARG00000040503 
# si:ch211-66e2.5:ENSDARG00000091579, NA:ENSDARG00000054352, tac1:ENSDARG00000014490, rgs5a:ENSDARG00000002644, ptgdsb.1:ENSDARG00000027088, id1:ENSDARG00000040764, slc7a2:ENSDARG00000037097, cyfip1:ENSDARG00000044345, LOC100331987:ENSDARG00000100594, si:ch211-183d21.1:ENSDARG00000092499 
# Negative:  adcyap1b:ENSDARG00000027740, ebf3a:ENSDARG00000100244, necab2:ENSDARG00000056745, neurod1:ENSDARG00000019566, nhlh2:ENSDARG00000025495, nfixa:ENSDARG00000043226, insm1a:ENSDARG00000091756, nanos1:ENSDARG00000109337, plxna2:ENSDARG00000060372, nxph1:ENSDARG00000033447 
# mafbb:ENSDARG00000070542, NA:pEGFP-N1, elavl3:ENSDARG00000014420, olfm1b:ENSDARG00000014053, rassf1:ENSDARG00000004840, stap2a:ENSDARG00000092810, ebf2:ENSDARG00000042525, nfixb:ENSDARG00000061836, tuba8l3:ENSDARG00000070155, olig4:ENSDARG00000052610 
# ncam1a:ENSDARG00000056181, znf804a:ENSDARG00000027079, lrch4:ENSDARG00000068258, c1qtnf4:ENSDARG00000024299, si:dkeyp-34c12.1:ENSDARG00000071083, draxin:ENSDARG00000058256, gria2b:ENSDARG00000052765, tubb5:ENSDARG00000037997, hoxb8a:ENSDARG00000056027, myo1b:ENSDARG00000024694 
# PC_ 5 
# Positive:  islr2:ENSDARG00000051875, bmp16:ENSDARG00000103679, hunk:ENSDARG00000091317, olig4:ENSDARG00000052610, neurod4:ENSDARG00000003469, pdyn:ENSDARG00000087798, rbpjb:ENSDARG00000052091, sox4a:ENSDARG00000004588, dnajb1b:ENSDARG00000041394, lima1a:ENSDARG00000101441 
# ebf2:ENSDARG00000042525, usp43b:ENSDARG00000088072, itga3a:ENSDARG00000037917, lgi1b:ENSDARG00000058421, gab1:ENSDARG00000037018, snap25b:ENSDARG00000058117, atp2a3:ENSDARG00000060978, myo1b:ENSDARG00000024694, dcc:ENSDARG00000104282, gpc2:ENSDARG00000104217 
# hapln1a:ENSDARG00000089769, dhrsx:ENSDARG00000079280, syt1a:ENSDARG00000030614, slc18a3a:ENSDARG00000006356, dlc:ENSDARG00000002336, npas1:ENSDARG00000015876, shtn3:ENSDARG00000056519, insm1a:ENSDARG00000091756, pnoca:ENSDARG00000025024, esrrga:ENSDARG00000004861 
# Negative:  rbfox3a:ENSDARG00000010083, trarg1a:ENSDARG00000090481, calb2b:ENSDARG00000036344, c1qtnf4:ENSDARG00000024299, ppp1r13ba:ENSDARG00000004377, gabra5:ENSDARG00000070730, necab2:ENSDARG00000056745, chrna2b:ENSDARG00000057025, cbln2b:ENSDARG00000077151, eef1a2:ENSDARG00000006838 
# pcdh7a:ENSDARG00000078898, kcna1a:ENSDARG00000062942, nptx1l:ENSDARG00000074671, slit1b:ENSDARG00000099446, nefmb:ENSDARG00000043697, calca:ENSDARG00000056590, chgb:ENSDARG00000076500, igfbp5b:ENSDARG00000025348, mast3b:ENSDARG00000086505, zgc:73226:ENSDARG00000023759 
# LOC101886435:ENSDARG00000104767, hpca:ENSDARG00000018397, cbln4:ENSDARG00000061240, slc17a6b:ENSDARG00000041150, cntn3a.1:ENSDARG00000062880, rnd3b:ENSDARG00000007396, gria2b:ENSDARG00000052765, hoxb1b:ENSDARG00000054033, megf11:ENSDARG00000062686, si:ch211-207d6.2:ENSDARG00000031658

# Plot PCAs
png('12.PCAplot.PC1vsPC2.vst.plate.counts.png')
PCAPlot(object = p174.p358_i.filt, dims=c(1,2),group.by="plate")
dev.off()

#with the integration the two plates are better mixed

ElbowPlot(p174.p358_i.filt)
# I would say 10 or 15
p174.p358_i.filt <- JackStraw(p174.p358_i.filt)
p174.p358_i.filt <- ScoreJackStraw(p174.p358_i.filt, dims = 1:20)
pcs <- which(JS(object = p174.p358_i.filt[['pca']], slot = 'overall')[, 2] < 1e-3)
pcs
#  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 15 17


###
# Cluster with Seurat
# Seurat should be fine for >500 cells
# UMAP and tSNE just for visualization

# DECIDE HOW MANY PCs for clustering (from elbow plot AND JACKSTRAW)
p174.p358_i.filt <- FindNeighbors(p174.p358_i.filt, reduction = "pca", dims = 1:10)
p174.p358_i.filt <- FindClusters(p174.p358_i.filt, resolution = 0.9)
table(Idents(p174.p358_i.filt))

#10 dim res 0.9      
# 0   1   2   3   4   5   6   7 
# 114 103  75  74  72  62  53  53 

# Run UMAP and tSNE
p174.p358_i.filt <- RunUMAP(p174.p358_i.filt, dims = 1:10)
p174.p358_i.filt <- RunTSNE(p174.p358_i.filt, dims = 1:10)

p1 <- DimPlot(p174.p358_i.filt, reduction = "umap", group.by = "plate")
p2 <- DimPlot(p174.p358_i.filt, reduction = "umap")
p3 <- DimPlot(p174.p358_i.filt, reduction = "tsne", group.by = "plate")
p4 <- DimPlot(p174.p358_i.filt, reduction = "tsne")
grid.arrange (p1, p2, p3,p4, ncol = 2)
#with the integration the two plates are better mixed

pdf('13.TSNE_Seurat.clusters.p358p174.pdf')
DimPlot(p174.p358_i.filt, reduction = "tsne")
dev.off()

# Check clusters over resolutions
res <- c(0.4,0.5,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4)
tmp=p174.p358_i.filt
tmp <- FindClusters(tmp, resolution = res, print.output = 0)
png('14.clustree.seurat.clusters.p174.3582.png')
clustree(tmp)
dev.off()
# 

# Differential expression

p174.p358_i.filt.markers <- FindAllMarkers(p174.p358_i.filt, assay = "RNA", test.use= "MAST", min.pct=0.25, logfc.threshold=0.3) 
p174.p358_i.filt.markers.pos <- FindAllMarkers(p174.p358_i.filt, assay = "RNA", test.use = "MAST", min.pct=0.25, logfc.threshold=0.3, only.pos = TRUE)

# Write list to file
write.table(p174.p358_i.filt.markers %>% group_by(cluster), file="New.allMarkerGenes.new.p174.358.txt",quote=FALSE,row.names=FALSE)
write.table(p174.p358_i.filt.markers.pos %>% group_by(cluster), file="New.allMarkerGenes.onlyPos.new.p174.358.txt",quote=FALSE,row.names=FALSE)

# Heatmap top10 marker genes (sorted by padj) for each cluster
# use the integrated assay for clustering and other analysis but always use the raw data (RNA) when plotting gene expression

DefaultAssay(p174.p358_i.filt) <- "RNA"
feat = p174.p358_i.filt.markers %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('15.heatmap.markerGenes.top10padj.p358.174.png', width=920, height=1200)
DoHeatmap(p174.p358_i.filt, features = feat$gene,  slot="data")
dev.off()

feat.pos = p174.p358_i.filt.markers.pos %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('16.heatmap.markerGenes.top10padj.onlyPos.p358.174.png', width=920, height=1200)
DoHeatmap(p174.p358_i.filt, features = feat.pos$gene,  slot="data")
dev.off()

#reordering clusters for Ext Data Fig. 2
my_levels.tfac <- c(0, 4, 2, 6, 1, 7, 3, 5)
p174.p358_i.filt@active.ident <- factor(x = p174.p358_i.filt@active.ident, levels = my_levels.tfac)

pdf('dotplot.islet.MNmarkers.p358.174.pdf')
DotPlot(p174.p358_i.filt, features = mn.markers,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()

# Working on MN only#####

# subsetting MN clusters 0,2,4,6 that have MN markers
MN.only <- subset(x = p174.p358_i.filt, subset = seurat_clusters == 0 | seurat_clusters == 2 | seurat_clusters ==4 | seurat_clusters == 6 )
DefaultAssay(MN.only) <- "integrated"
MN.only <- FindVariableFeatures(MN.only, selection.method = "vst", nfeatures = 2000)
MN.only <- ScaleData(MN.only, vars.to.regress = c("nFeature_RNA", "nCounts_RNA", "percent_ribo"))

#PCA 
MN.only <- RunPCA(MN.only, features = VariableFeatures(MN.only), verbose = FALSE)
PCAPlot(object = MN.only, dims=c(1,2), group.by="seurat_clusters")
ElbowPlot(MN.only)
MN.only <- JackStraw(MN.only)
MN.only <- ScoreJackStraw(MN.only, dims = 1:20)
JackStrawPlot(MN.only, dims = 1:20)
pcs <- which(JS(object = MN.only[['pca']], slot = 'overall')[, 2] < 1e-3)
pcs
#[1]  1  2  3  4  5  6  7 20

#reclustering
MN.only <- FindNeighbors(MN.only, reduction = "pca", dims = 1:9)
MN.only <- FindClusters(MN.only, resolution = 0.8)
MN.only <- RunUMAP(MN.only, dims = 1:9)
MN.only <- RunTSNE(MN.only, dims = 1:9)
table(Idents(MN.only))
#9pcs 0.8
# 0  1  2  3  4 
#98 93 51 47 27 

#clustree
DefaultAssay(MN.only)<-"integrated"
tmp=MN.only
tmp <- FindClusters(tmp, resolution = res, print.output = 0)
pdf('17.clustree.seurat.reclustering.MN.p174.3582.pdf')
clustree(tmp)
dev.off()

p1 <- DimPlot(MN.only, reduction = "tsne", group.by = "plate", pt.size = 2)
p2 <- DimPlot(MN.only, reduction = "tsne", pt.size = 2)
p3 <- DimPlot(MN.only, reduction = "umap", group.by = "plate", pt.size = 2)
p4 <- DimPlot(MN.only, reduction = "umap", pt.size = 2)


#p3 <- DimPlot(tmp, reduction = "tsne", group.by = "pctspike")
png('18.TSNE_UMAP_Seurat.MNonlyclusters.png', width=480, height=800)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

pdf('19.TSNE_Seurat.MNonlyclusters.pdf')
DimPlot(MN.only, reduction = "tsne")
dev.off()

# Differential expression analysis
DefaultAssay(MN.only) <- "RNA"
MN.only.markers <- FindAllMarkers(MN.only, test.use= "MAST", min.pct=0.25, logfc.threshold=0.3) 
MN.only.markers.pos <- FindAllMarkers(MN.only, test.use = "MAST", min.pct=0.25, logfc.threshold=0.3, only.pos = TRUE)


# Write list to file
write.table(MN.only.markers %>% group_by(cluster), file="MN.allMarkerGenes.txt",quote=FALSE,row.names=FALSE)
write.table(MN.only.markers.pos %>% group_by(cluster), file="MN.allMarkerGenes.onlyPos.txt",quote=FALSE,row.names=FALSE)

# reorder clusters
my_levels <- c(0, 2, 4, 3, 1)
MN.only@active.ident <- factor(x = MN.only@active.ident, levels = my_levels)

# Heatmap top10 marker genes (sorted by padj) for each cluster
DefaultAssay(MN.only) <- "RNA"
feat = MN.only.markers %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('24.heatmap.MN.markerGenes.top10padj.png', width=920, height=1200)
DoHeatmap(MN.only, features = feat$gene,  slot="data")
dev.off()

feat.pos = MN.only.markers.pos %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('20.heatmap.MN.markerGenes.top10padj.onlyPos.png', width=920, height=1200)
DoHeatmap(MN.only, features = feat.pos$gene,  slot="data")
dev.off()

#Save object for integration with V2a
saveRDS(MN.only, "MN.RDS")

#Plots for figures####
setwd("MN/Figures")
DefaultAssay(MN.only) <- "RNA"
#Fig.1a
pdf('UMAP_Seurat.MNonlyclusters.pdf')
DimPlot(MN.only, reduction = "umap")+xlim(-6,6)+ylim(-6,6)
dev.off()

#Fig.1b
Fig_MN_cluster_markers<-read.table("MN markers per cluster for figure.txt")

pdf('20.dotplot.MN.markerGenes.top10padj.onlyPos.pdf')
DotPlot(MN.only, features = Fig_MN_cluster_markers,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()

#Fig.1c
Fig.1c<-read.table("Fig.1c.txt")
Fig.1c<- lapply(Fig.1c, as.character)$V1 
pdf("Fig.1c.pdf")
DotPlot(MN.only, features = Fig.1c, cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()

#Fig1e-h
Fig1e.h<-read.table("Fig1e-h.txt")
Fig1e.h<- lapply(Fig1e.h, as.character)$V1 

pdf("Fig1e-h.pdf", width = 5, height = 13)
StackedVlnPlot(MN.only, features = Fig1e.h)
dev.off()


#Ext data Fig2a,b
pdf('13.UMAP_Seurat.clusters.p358p174.pdf')
DimPlot(p174.p358_i.filt, reduction = "umap")
dev.off()

# Ext Data Fig. 2c-g
pdf('neuronalmarkers.p174.p358.for figure.pdf', width =  6, height = 6)
StackedVlnPlot(p174.p358_i.filt, features = c("elavl3:ENSDARG00000014420", "NA:ENSDARG00000055052", "elavl4:ENSDARG00000045639"))
dev.off()

pdf('MNmarkers.p174.p358.for figure.pdf')
StackedVlnPlot(p174.p358_i.filt, features = c("mnx1:ENSDARG00000035984", "slc18a3a:ENSDARG00000006356", "chata:ENSDARG00000015854"))
dev.off()

nonMNclusters<-read.table("nonMNclusters.markers.txt")
nonMNclusters <- lapply(nonMNclusters, as.character)$V1 

pdf('nonMNmarkers.pdf', width =  6, height = 18)
StackedVlnPlot(p174.p358_i.filt, features = nonMNclusters)
dev.off()

#reordering clusters for Ext Data Fig. 2
my_levels.tfac <- c(0, 4, 2, 6, 1, 7, 3, 5)
p174.p358_i.filt@active.ident <- factor(x = p174.p358_i.filt@active.ident, levels = my_levels.tfac)

islet_cluster_markers<-read.table("islet markers per cluster.txt")

pdf('dotplot.islet.markerGenes.p358.174.pdf')
DotPlot(p174.p358_i.filt, features = islet_cluster_markers,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()

#Fig. 2a,b,g,h and Ext Data Fig.3a
DefaultAssay(MN.only) <- "RNA"

pdf('esrrga featureplot.pdf')
FeaturePlot(MN.only, features = "esrrga:ENSDARG00000004861",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-6,6)+ylim(-6,6)
dev.off()

pdf('hoxb13 featureplot.pdf')
FeaturePlot(MN.only, features = "hoxb13a:ENSDARG00000056015",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-6,6)+ylim(-6,6)
dev.off()

pdf('grin1b featureplot.pdf')
FeaturePlot(MN.only, features = "grin1b:ENSDARG00000025728",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-6,6)+ylim(-6,6)
dev.off()

pdf('pvalb6 featureplot.pdf')
FeaturePlot(MN.only, features = "pvalb6:ENSDARG00000009311",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-6,6)+ylim(-6,6)
dev.off()

pdf('chrna2b featureplot.pdf')
FeaturePlot(MN.only, features = "chrna2b:ENSDARG00000057025",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-6,6)+ylim(-6,6)
dev.off()

pdf('neurod1 featureplot.pdf')
FeaturePlot(MN.only, features = "neurod1:ENSDARG00000019566",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-6,6)+ylim(-6,6)
dev.off()

#Fig.3b
pdf("esrrga Violin plot.pdf")
VlnPlot(MN.only, features = "esrrga:ENSDARG00000004861")
dev.off()


#Fig.6b
BestSF<-read.table("Slow_Fast_combined_MNandV2a.txt")
BestSF<- lapply(BestSF, as.character)$V1 
DoHeatmap(MN.only, features = BestSF, slot="data")
pdf('BestSF_MN.pdf')
DotPlot(MN.only, features = BestSF,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()

setwd("MN")


#Extract counts for Fig2m analysis outside R

probes.counts<- FetchData(MN.only, vars = probes, slot = 'counts')
head(probes.counts)
cluster<- FetchData(MN.only, vars = 'seurat_clusters')
cluster<-as.character(cluster$seurat_clusters)
probes.counts$cluster<- cluster
write.csv(probes.counts, "probes.counts.MNs.csv")

#clustering by TFs####
setwd("MN/TF")
#for this analysis we excluded the immature cluster
MN.noRibo <- subset(x = MN.only, subset = seurat_clusters == 0 | seurat_clusters == 2 | seurat_clusters ==3 | seurat_clusters == 4 )
tfactors<-read.table("TF MN.txt")
tfactors <- lapply(tfactors, as.character)$V1 

DefaultAssay(MN.noRibo)
DefaultAssay(MN.noRibo) <- "integrated"
MN.noRibo <- ScaleData(MN.noRibo, features = tfactors, vars.to.regress = c("nFeature_RNA", "nCounts_RNA", "percent_ribo"))

# PCA on transcription factor genes only
MN.TF <- RunPCA(MN.TF, features = tfactors, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
# PC_ 1 
# Positive:  jun:ENSDARG00000043531, sox11a:ENSDARG00000077811, sox4a:ENSDARG00000004588, bhlhe41:ENSDARG00000041691, zfhx3:ENSDARG00000103057, junbb:ENSDARG00000104773, mafgb:ENSDARG00000100097, fosab:ENSDARG00000031683, zfhx4:ENSDARG00000075542, si:ch73-386h18.1:ENSDARG00000073944 
# npas3:ENSDARG00000109820, pbx1b:ENSDARG00000101131, lbx1a:ENSDARG00000018321, pbx3b:ENSDARG00000013615, junba:ENSDARG00000074378, meis1b:ENSDARG00000012078, atf3:ENSDARG00000007823, egr2a:ENSDARG00000044098, cebpb:ENSDARG00000042725, dachd:ENSDARG00000069440 
# mxi1:ENSDARG00000040884, tp53:ENSDARG00000035559, maff:ENSDARG00000028957, klf11b:ENSDARG00000013794, pknox1.1:ENSDARG00000018765, klf6a:ENSDARG00000029072, LOC110439853:ENSDARG00000112811, si:dkey-14k9.3:ENSDARG00000045595, hsf2:ENSDARG00000053097, smad1:ENSDARG00000027199 
# Negative:  nfixb:ENSDARG00000061836, neurod1:ENSDARG00000019566, sox2:ENSDARG00000070913, nfixa:ENSDARG00000043226, zbtb18:ENSDARG00000028228, neurod6b:ENSDARG00000020794, irx3a:ENSDARG00000101076, irx5a:ENSDARG00000034043, mafbb:ENSDARG00000070542, nkx2.2b:ENSDARG00000101549 
# hoxb6b:ENSDARG00000026513, stat5a:ENSDARG00000019392, mnx1:ENSDARG00000035984, neurod2:ENSDARG00000016854, neurod6a:ENSDARG00000040008, isl2a:ENSDARG00000003971, mafa:ENSDARG00000015890, scrt1b:ENSDARG00000040214, hoxb7a:ENSDARG00000056030, sp5a:ENSDARG00000076571 
# lhx4:ENSDARG00000039458, znf740a:ENSDARG00000070939, irx1a:ENSDARG00000101831, si:ch211-161m3.6:ENSDARG00000091393, arid3a:ENSDARG00000070843, znf1040:ENSDARG00000098986, zgc:113295:ENSDARG00000074768, si:ch211-232b12.5:ENSDARG00000087186, plag1:ENSDARG00000051926, znf1147:ENSDARG00000100652 
# PC_ 2 
# Positive:  zfhx4:ENSDARG00000075542, npas4b:ENSDARG00000087753, irf5:ENSDARG00000045681, lhx4:ENSDARG00000039458, sall3b:ENSDARG00000057586, fosl2:ENSDARG00000040623, lbx1a:ENSDARG00000018321, sall1a:ENSDARG00000074319, myog:ENSDARG00000009438, epas1a:ENSDARG00000008697 
# stat2:ENSDARG00000031647, junba:ENSDARG00000074378, rxrba:ENSDARG00000078954, bcl11aa:ENSDARG00000061352, irf2:ENSDARG00000040465, spic:ENSDARG00000012435, klf11b:ENSDARG00000013794, nkx6.1:ENSDARG00000022569, id3:ENSDARG00000054823, npas4a:ENSDARG00000055752 
# sp8b:ENSDARG00000056666, her13:ENSDARG00000007097, atf3:ENSDARG00000007823, zbtb40:ENSDARG00000091762, raraa:ENSDARG00000056783, scxa:ENSDARG00000035695, klf4:ENSDARG00000079922, lrrfip1b:ENSDARG00000095170, zgc:153116:ENSDARG00000077205, znfl2a:ENSDARG00000008333 
# Negative:  nr2f1a:ENSDARG00000052695, nr2f2:ENSDARG00000040926, tshz3b:ENSDARG00000103361, meis2a:ENSDARG00000098240, pbx1b:ENSDARG00000101131, mafga:ENSDARG00000018109, mnx1:ENSDARG00000035984, mxd4:ENSDARG00000054031, rbpjb:ENSDARG00000052091, bhlhe23:ENSDARG00000037588 
# si:cabz01069013.3:ENSDARG00000088367, im:7141269:ENSDARG00000101756, grhl2a:ENSDARG00000058342, smad1:ENSDARG00000027199, isl2b:ENSDARG00000053499, tshz1:ENSDARG00000005026, scrt2:ENSDARG00000056175, nfkb1:ENSDARG00000105261, nkx6.2:ENSDARG00000104735, rfx5:ENSDARG00000063258 
# si:ch211-222k6.1:ENSDARG00000071581, bhlhe41:ENSDARG00000041691, sox12:ENSDARG00000025847, nkx1.2lb:ENSDARG00000099427, erg:ENSDARG00000077304, prox3:ENSDARG00000088810, pknox1.1:ENSDARG00000018765, etv1:ENSDARG00000101959, onecutl:ENSDARG00000040253, isl2a:ENSDARG00000003971 
# PC_ 3 
# Positive:  foxb1a:ENSDARG00000089042, id2b:ENSDARG00000029544, nfia:ENSDARG00000062420, camta1a:ENSDARG00000077428, foxb1b:ENSDARG00000053650, tcf20:ENSDARG00000078348, prox3:ENSDARG00000088810, brpf3a:ENSDARG00000070513, tshz3b:ENSDARG00000103361, mxi1:ENSDARG00000040884 
# relb:ENSDARG00000086173, max:ENSDARG00000024844, zbtb38:ENSDARG00000100607, znf292a:ENSDARG00000016763, hnrnpa0l:ENSDARG00000036161, brd1b:ENSDARG00000051798, fosab:ENSDARG00000031683, dachd:ENSDARG00000069440, mxd4:ENSDARG00000054031, erg:ENSDARG00000077304 
# mafba:ENSDARG00000017121, si:ch73-386h18.1:ENSDARG00000073944, pbx3a:ENSDARG00000089262, pknox1.1:ENSDARG00000018765, zbtb34:ENSDARG00000013492, si:dkey-14k9.3:ENSDARG00000045595, fosb:ENSDARG00000055751, nfil3-6:ENSDARG00000087188, stat2:ENSDARG00000031647, mef2b:ENSDARG00000093170 
# Negative:  tbr1a:ENSDARG00000076259, runx3:ENSDARG00000052826, neurod4:ENSDARG00000003469, cbfb:ENSDARG00000040917, mafk:ENSDARG00000100947, nr1d4a:ENSDARG00000031161, nfat5b:ENSDARG00000102556, lmx1al:ENSDARG00000077915, rfx4:ENSDARG00000026395, sp7:ENSDARG00000019516 
# zgc:101562:ENSDARG00000040179, zgc:171901:ENSDARG00000074365, mixl1:ENSDARG00000069252, znf1157:ENSDARG00000094736, LOC100334054:ENSDARG00000100883, znf1040:ENSDARG00000098986, si:ch73-221f6.1:ENSDARG00000104904, tp63:ENSDARG00000044356, zic1:ENSDARG00000015567, LOC563292:ENSDARG00000098424 
# zgc:112083:ENSDARG00000037405, si:ch211-110e21.4:ENSDARG00000076807, rfx6:ENSDARG00000041702, epas1b:ENSDARG00000057671, onecut2:ENSDARG00000090387, glis2a:ENSDARG00000078388, znf281b:ENSDARG00000035910, si:rp71-1g18.1:ENSDARG00000053792, atf2:ENSDARG00000023903, elk1:ENSDARG00000078066 
# PC_ 4 
# Positive:  foxb1a:ENSDARG00000089042, foxb1b:ENSDARG00000053650, klf6b:ENSDARG00000038561, LOC563292:ENSDARG00000098424, LOC101882267:ENSDARG00000099172, si:ch211-232b12.5:ENSDARG00000087186, zbtb34:ENSDARG00000013492, tshz1:ENSDARG00000005026, nhlh2:ENSDARG00000025495, runx3:ENSDARG00000052826 
# mafk:ENSDARG00000100947, hey2:ENSDARG00000013441, mafbb:ENSDARG00000070542, lmx1al:ENSDARG00000077915, stat2:ENSDARG00000031647, neurod4:ENSDARG00000003469, klf12b:ENSDARG00000032197, mef2ca:ENSDARG00000029764, sall1a:ENSDARG00000074319, irx4b:ENSDARG00000036051 
# sp6:ENSDARG00000099880, mllt10:ENSDARG00000045401, znf438:ENSDARG00000061311, rxrgb:ENSDARG00000004697, mkxa:ENSDARG00000075450, sp7:ENSDARG00000019516, LOC100334054:ENSDARG00000100883, zgc:101562:ENSDARG00000040179, si:ch73-138e16.2:ENSDARG00000104894, nfat5b:ENSDARG00000102556 
# Negative:  meis2a:ENSDARG00000098240, isl2b:ENSDARG00000053499, nfkb2:ENSDARG00000038687, onecut1:ENSDARG00000007982, camta1b:ENSDARG00000007824, bhlhe23:ENSDARG00000037588, znf770:ENSDARG00000070786, pbx1b:ENSDARG00000101131, junbb:ENSDARG00000104773, dbpb:ENSDARG00000057652 
# znf1001:ENSDARG00000103777, fosab:ENSDARG00000031683, smad4a:ENSDARG00000075226, im:7141269:ENSDARG00000101756, znf740a:ENSDARG00000070939, en2a:ENSDARG00000026599, zgc:175107:ENSDARG00000075417, LOC100534733:ENSDARG00000102411, si:ch211-166g5.4:ENSDARG00000063555, si:ch211-161m3.4:ENSDARG00000103375 
# nkx2.5:ENSDARG00000018004, fosl1b:ENSDARG00000079373, myclb:ENSDARG00000034956, mafb:ENSDARG00000076520, zeb2a:ENSDARG00000062338, nkx6.2:ENSDARG00000104735, LOC100535338:ENSDARG00000112895, znf1155:ENSDARG00000100192, si:ch1073-127d16.1:ENSDARG00000087674, nfic:ENSDARG00000043210 
# PC_ 5 
# Positive:  foxb1a:ENSDARG00000089042, si:ch211-232b12.5:ENSDARG00000087186, cebpa:ENSDARG00000036074, cebpb:ENSDARG00000042725, mafbb:ENSDARG00000070542, znf1184:ENSDARG00000101694, egr3:ENSDARG00000089156, camta1a:ENSDARG00000077428, cebpg:ENSDARG00000036073, arnt:ENSDARG00000021855 
# si:ch211-79k12.2:ENSDARG00000090679, foxl2a:ENSDARG00000042180, id2b:ENSDARG00000029544, mafga:ENSDARG00000018109, sox19a:ENSDARG00000010770, atf2:ENSDARG00000023903, tshz3a:ENSDARG00000087394, znf1001:ENSDARG00000103777, tshz1:ENSDARG00000005026, foxd3:ENSDARG00000021032 
# si:dkey-77f5.4:ENSDARG00000097808, stat5b:ENSDARG00000055588, foxb1b:ENSDARG00000053650, si:ch73-221f6.1:ENSDARG00000104904, znf142:ENSDARG00000061373, znf982:ENSDARG00000088507, zgc:173720:ENSDARG00000100479, erg:ENSDARG00000077304, znf362b:ENSDARG00000062231, znf516:ENSDARG00000070809 
# Negative:  nkx6.2:ENSDARG00000104735, si:cabz01069013.3:ENSDARG00000088367, hey2:ENSDARG00000013441, nfkb1:ENSDARG00000105261, zeb2a:ENSDARG00000062338, myt1b:ENSDARG00000102879, prdm5:ENSDARG00000006288, onecut1:ENSDARG00000007982, elf2b:ENSDARG00000079626, tbx16l:ENSDARG00000006939 
# en2a:ENSDARG00000026599, epas1b:ENSDARG00000057671, sall1b:ENSDARG00000075891, zgc:162936:ENSDARG00000069957, znf1059:ENSDARG00000074009, zgc:66443:ENSDARG00000057707, hnf1ba:ENSDARG00000006615, stat2:ENSDARG00000031647, si:ch211-253b8.2:ENSDARG00000089230, bhlhe23:ENSDARG00000037588 
# znf451:ENSDARG00000076939, nfe2l2a:ENSDARG00000042824, nr2f1b:ENSDARG00000017168, egr2a:ENSDARG00000044098, mntb:ENSDARG00000073988, tbx2b:ENSDARG00000006120, znf507:ENSDARG00000052164, stat4:ENSDARG00000028731, irf5:ENSDARG00000045681, tshz3b:ENSDARG00000103361 
# Warning messages:
#   1: In PrepDR(object = object, features = features, verbose = verbose) :
#   The following 342 features requested have not been scaled (running reduction without them): NA:ENSDARG00000003803, sp5b:ENSDARG00000005846, spdef:ENSDARG00000029930, NA:ENSDARG00000039935, her15.2:ENSDARG00000054560, si:ch211-216l23.1:ENSDARG00000059707, tbx3b:ENSDARG00000061509, LOC571757:ENSDARG00000062621, arid3c:ENSDARG00000067729, znf618:ENSDARG00000071469, zc3h6:ENSDARG00000075165, LOC571865:ENSDARG00000079972, znf521:ENSDARG00000086825, LOC797319:ENSDARG00000086831, NA:ENSDARG00000088039, LOC100537716:ENSDARG00000088051, NA:ENSDARG00000090304, NA:ENSDARG00000090740, NA:ENSDARG00000091679, LOC100151333:ENSDARG00000092765, si:ch73-120g24.5:ENSDARG00000095030, LOC100535113:ENSDARG00000095605, NA:ENSDARG00000096107, NA:ENSDARG00000098122, NA:ENSDARG00000098248, NA:ENSDARG00000098811, NA:ENSDARG00000099025, znf1070:ENSDARG00000099048, NA:ENSDARG00000099070, NA:ENSDARG00000099670, NA:ENSDARG00000099890, NA:ENSDARG00000099945, NA:ENSDARG00000100262, NA:ENSDARG00000100329, NA:ENSDARG000 [... truncated]
# 2: In PrepDR(object = object, features = features, verbose = verbose) :
#   The following 4 features requested have zero variance (running reduction without them): nkx1.2la:ENSDARG00000006350, cdx4:ENSDARG00000036292, twist3:ENSDARG00000019646, her2:ENSDARG00000038205


png('Tfac.PCAplot.PC1vsPC2.vst.plate.counts.png')
PCAPlot(object = MN.TF, dims=c(1,2),group.by="plate")
dev.off()

ElbowPlot(MN.TF)

MN.TF <- JackStraw(MN.TF)
MN.TF <- ScoreJackStraw(MN.TF, dims = 1:20)
pcs <- which(JS(object = MN.TF[['pca']], slot = 'overall')[, 2] < 1e-3)
pcs
# [1] 1  2  4  5 12

###
# Cluster with Seurat

# DECIDE HOW MANY PCs for clustering (from elbow plot AND JACKSTRAW)
MN.TF <- FindNeighbors(MN.TF, reduction = "pca", dims = 1:12)
MN.TF <- FindClusters(MN.TF, resolution = 0.5)
table(Idents(MN.TF))

#pcs 12 res 0.5
#  0   1   2 
#140  50  33

# Check clusters over resolutions
res <- c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4)
tmp=MN.TF
tmp <- FindClusters(tmp, resolution = res, print.output = 0)
png('Tfac.clustree.seurat.clusters.TFnoribo.png')
clustree(tmp)
dev.off()
# clusters are stable

# Run UMAP and tSNE
MN.TF <- RunUMAP(MN.TF, dims = 1:12)
MN.TF <- RunTSNE(MN.TF, dims = 1:12)
p1 <- plotUMAP(object = MN.noRibo.TF2, colour_by="plate")
p2 <- DimPlot(MN.TF, reduction = "umap", group.by = "plate")
p3 <- DimPlot(MN.TF, reduction = "umap")
p4 <- DimPlot(MN.TF, reduction = "tsne", group.by = "plate")
p5 <- DimPlot(MN.TF, reduction = "tsne")
multiplot(p1, p2, p3, ncol = 1)
multiplot(p1, p2, p3,p4,p5, ncol = 2)

p1 <- DimPlot(MN.TF, reduction = "tsne", group.by = "plate", pt=2)
p2 <- DimPlot(MN.TF, reduction = "tsne", pt=2)
png('TfacNoribo.TSNE_Seurat.clusters.png', width=480, height=800)
multiplot(p1, p2, ncol = 1)
dev.off()


#plot for Ext Data Fig.5a
setwd("MN/Figures")
MN.only <- AddMetaData(MN.only, MN.TF@meta.data["seurat_clusters"], col.name = "tfac_clusters")
pdf('TFac.UMAP_MN.pdf')
DimPlot(MN.only, reduction = "umap", group.by = "tfac_clusters")+xlim(-6,6)+ylim(-6,6)
dev.off()
setwd("MN")

###
# Differential expression 
DefaultAssay(MN.TF) <- "RNA"
MN.TF.markers <- FindAllMarkers(MN.TF, test.use= "MAST", min.pct=0.25, logfc.threshold=0.3) 
MN.TF.markers.pos <- FindAllMarkers(MN.TF, test.use = "MAST", min.pct=0.25, logfc.threshold=0.3, only.pos = TRUE)

# Write list to file
write.table(MN.TF.markers %>% group_by(cluster), file="TfacNoribo.allMarkerGenes.txt",quote=FALSE,row.names=FALSE)
write.table(MN.TF.markers.pos %>% group_by(cluster), file="TfacNoribo.allMarkerGenes.onlyPos.txt",quote=FALSE,row.names=FALSE)

#reordering clusters for figure
my_levels.tfac <- c(0, 2, 1)
MN.TF@active.ident <- factor(x = MN.TF@active.ident, levels = my_levels.tfac)

# Heatmap top10 marker genes (sorted by padj) for each cluster
DefaultAssay(MN.TF) <- "RNA"
feat = MN.TF.markers %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('Tfac.heatmap.markerGenes.top10padj.png', width=920, height=1200)
DoHeatmap(MN.TF, features = feat$gene,  slot="data")
dev.off()

feat.pos = MN.TF.markers.pos %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('Tfac.heatmap.markerGenes.top10padj.onlyPos.png', width=920, height=1200)
DoHeatmap(MN.TF, features = feat.pos$gene,  slot="data")
dev.off()

setwd("MN/Figures")
#dotplot top10 TF DE for Ext Data Fig.5b
top10tfactors<-read.table("MN.TF_DE_first10.txt")
top10tfactors <- lapply(top10tfactors, as.character)$V1 


pdf('Tfac.tfactors.top10padj.onlyPos.pdf')
DotPlot(MN.only.tfac, features = top10tfactors,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()