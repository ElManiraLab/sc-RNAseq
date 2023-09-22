# Load R packages/libraries and functions
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
library('ape')

source("multiplot.R")
source("vlnplot.R")


setwd("V2a")

# Object creation, QC and filtering####
# Create Seurat objects, store metadata and QC measures

# Read raw counts and metadata (the file includes three plates, only 353 and 354 are V2a)
counts <- read.table("C358.354.353.merge.rpkmforgenes_counts.csv", header=TRUE, sep=',', row.names='gene')
dim(counts)
#[1] 32612  1152


# Add gene names to IDs
genes <- rownames(counts)
table(genes == rownames(counts))
#  TRUE 
#32612 
symb <- mapIds(org.Dr.eg.db, keys=genes, keytype="ENSEMBL", column="SYMBOL")
#'select()' returned 1:many mapping between keys and columns
m <- match(rownames(counts),names(symb)) 
newnames <- apply(cbind(as.vector(symb)[m],rownames(counts)),1,paste,collapse=":")
rownames(counts) <- newnames

# Read metadata
info <- read.table("358.354.353.sample_info.csv", header=TRUE, sep=',', row.names='SM')

# Create Seurat objects, take only genes expressed in more than 3 cells
seu_c <- CreateSeuratObject(counts, project = "SS2_18_c", meta.data = info, min.cells =3)
#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

# Read and store mapping info
star <- read.table("358.354.353.multiqc_star.txt", header=TRUE, sep='\t', row.names='Sample')
seu_c <- AddMetaData(seu_c, star$total_reads, col.name = "total_reads")
seu_c<- AddMetaData(seu_c, star$uniquely_mapped, col.name = "uniquely_mapped")
seu_c <- AddMetaData(seu_c, star$uniquely_mapped_percent, col.name = "uniquely_mapped_percent")
seu_c <- AddMetaData(seu_c, star$multimapped_percent, col.name = "multimapped_percent")

# Exclude plate 358 (MN)
selected <- WhichCells(seu_c, expression = plate == "p358", invert = TRUE )
seu_c <- subset(seu_c, cells = selected)


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

setwd("V2a/QC")
png('1.VlnPlot.nCount-nFeature_byPlate.png', width=480, height=250)
VlnPlot(object = seu_c, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")
dev.off()

# Scatter plot
png('2.Scatter.nCount-nFeature.png', width=480, height=250)
FeatureScatter(seu_c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
dev.off()
 

# Percent mitochondrial and ribosomal and spikein and pgfp
png('3.VlnPlot_pctmito-pctribo-pctspike-pctgfp.counts.png')
p1 <- VlnPlot(seu_c, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
p2 <- VlnPlot(seu_c, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)
multiplot(p1, p2, cols = 1)
dev.off()


# Filter out cells with high (>25%) spike-in
selected <- names(which(seu_c$percent.spike < 25))
length(selected)
#[1] 754

# Subset to keep only cells with spikein < 25%
seu_c.filt <- subset(seu_c, cells = selected)
Idents(seu_c) <- 'plate'
Idents(seu_c.filt) <- 'plate'
# How many cells are left?
table(Idents(seu_c))
#p353 p354 
#384  384 
table(Idents(seu_c.filt))
#p353 p354 
#376  378  

# Plot the proportion of mitochondrial counts against some of the other QC metrics.
png('3b.Scatter.mito_count.mito_feature.p174.358.png')
p1 <- FeatureScatter(seu_c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
p2 <- FeatureScatter(seu_c, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by="plate")
p3 <- FeatureScatter(seu_c, feature1 = "nFeature_RNA", feature2 = "percent.mito", group.by="plate")
multiplot(p1, p2, p3, cols = 1)
dev.off()

# Mitochondrial filtering
selected <- WhichCells(seu_c.filt, expression = percent.mito < 30)
length(selected)
#[1] 691
seu_c.filt <- subset(seu_c.filt, cells = selected)
table(Idents(seu_c.filt))
# p353 p354 
#339  352

# Percent mitochondrial and ribosomal and spikein and pgfp
png('4.VlnPlot_pctmito-pctribo-pctspike-pctgfp.counts.cellfilt.png')
p1 <- VlnPlot(seu_c.filt, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
p2 <- VlnPlot(seu_c.filt, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)
multiplot(p1, p2, cols = 1)
dev.off()

# Scatter plot
png('5.Scatter.count_feature_by plate.png')
FeatureScatter(seu_c.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
dev.off()

# Number of reads per cell and number of features / detected genes
png('5.VlnPlot.count_feature_by plate.png')
VlnPlot(object = seu_c.filt, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")
dev.off()

# Remove cells with very high gene detection (possibly doublets)
selected <- WhichCells(seu_c.filt, expression = nFeature_RNA < 9000)
length(selected)
#[1] 684
seu_c.filt <- subset(seu_c.filt, cells = selected)
table(Idents(seu_c.filt))
#p353 p354 
#339  345 

# Filter the cells with low gene detection (low quality libraries) 
selected <- WhichCells(seu_c.filt, expression = nFeature_RNA > 2100)
length(selected)
#593
seu_c.filt <- subset(seu_c.filt, cells = selected)
table(Idents(seu_c.filt))
#p353 p354 
# 282  311  

# Plots after all cell filtering
png('6.VlnPlot_pctmito-pctribo-pctspike-pctgfp.cellfiltered.png')
p1 <- VlnPlot(seu_c.filt, features = c("percent.mito", "percent.ribo"), group.by="plate", pt.size = 0.1)
p2 <- VlnPlot(seu_c.filt, features = c("percent.spike", "percent.gfp"), group.by="plate", pt.size = 0.1)
multiplot(p1, p2, cols = 1)
dev.off()

# Scatter and violin plots of nCounts and nFeature
png('7.Scatter.VlnPlot_nCount-nFeature.byPlate.cellfiltered.png', width=700, height=600)
p1 <- FeatureScatter(seu_c.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="plate")
p2 <- VlnPlot(object = seu_c.filt, features = c("nCount_RNA", "nFeature_RNA"), ncol = 2, pt.size = 0.1, group.by="plate")
multiplot(p1, p2, cols = 1)
dev.off()

setwd("V2a")

# Cluster with Seurat####
# Normalize
scale.factor_c <- mean(colSums(GetAssayData(seu_c.filt, slot = "counts")))
scale.factor_c
#[1] 339780.7
seu_c.filt <- NormalizeData(object = seu_c.filt, normalization.method = "LogNormalize", scale.factor = scale.factor_c)

# Feature selection
seu_c.filt <- FindVariableFeatures(seu_c.filt, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(seu_c.filt), 10)
# [1] "chgb:ENSDARG00000076500"    "tac1:ENSDARG00000014490"    "plp1b:ENSDARG00000011929"  
#[4] "atp1a3b:ENSDARG00000104139" "NA:ENSDARG00000104418"          
#[7] "angptl3:ENSDARG00000044365" "hoxb13a:ENSDARG00000056015" "fabp7a:ENSDARG00000007697" 
#[10] "tubb5:ENSDARG00000037997"  "mbpa:ENSDARG00000036186"

vplot_c <- VariableFeaturePlot(seu_c.filt)
top10_c <- head(VariableFeatures(seu_c.filt), 10)
labels <- c("chgb","tac1","plp1b","atp1a3b","NA:ENSDARG00000104418","angptl3","hoxb13a","fabp7a","tubb5","mbpa")
png('11.Variable_features.vst.counts.png')
LabelPoints(plot = vplot_c, points = top10_c, repel = TRUE, xnudge=0, ynudge=0, label=labels)
dev.off()

# Regress out number of detected genes, number of reads, percent_ribo
seu_c.filt <- ScaleData(seu_c.filt, vars.to.regress = c("nFeature_RNA", "nCounts_RNA", "percent_ribo"))

# PCA on variable genes
seu_c.filt <- RunPCA(seu_c.filt, features = VariableFeatures(seu_c.filt), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
# PC_ 1 
# Positive:  gja1b:ENSDARG00000041799, si:ch211-251b21.1:ENSDARG00000007275, gpr37l1b:ENSDARG00000056774, zfp36l1b:ENSDARG00000021443, prom1b:ENSDARG00000034007, si:ch73-173h19.3:ENSDARG00000104210, slc6a11b:ENSDARG00000087981, slc7a10b:ENSDARG00000051730, msrb2:ENSDARG00000018459, her6:ENSDARG00000006514 
# sc:d189:ENSDARG00000102858, slc4a4a:ENSDARG00000013730, entpd2b:ENSDARG00000044795, slc25a48:ENSDARG00000021250, cyp3c1:ENSDARG00000015575, s100b:ENSDARG00000057598, itgb5:ENSDARG00000012942, sox19a:ENSDARG00000010770, swap70b:ENSDARG00000057286, hspb15:ENSDARG00000078411 
# slc3a2a:ENSDARG00000036427, entpd1:ENSDARG00000045066, s100v2:ENSDARG00000070702, selenop:ENSDARG00000093549, atp1b4:ENSDARG00000053262, ca9:ENSDARG00000102300, fgfbp3:ENSDARG00000040162, si:ch211-212k18.5:ENSDARG00000103316, acsl2:ENSDARG00000078399, cyp2aa6:ENSDARG00000103590 
# Negative:  atp1b1b:ENSDARG00000076833, ache:ENSDARG00000031796, gapdhs:ENSDARG00000039914, eno1a:ENSDARG00000022456, pcdh17:ENSDARG00000027041, nefma:ENSDARG00000021351, islr2:ENSDARG00000051875, scg2b:ENSDARG00000038574, zgc:65851:ENSDARG00000012281, olfm3a:ENSDARG00000071493 
# pcbp3:ENSDARG00000054378, si:dkey-183i3.9:ENSDARG00000104326, egr4:ENSDARG00000077799, nr5a2:ENSDARG00000100940, si:dkey-35i13.1:ENSDARG00000097648, oxr1b:ENSDARG00000063310, ndnf:ENSDARG00000062936, spon1b:ENSDARG00000023694, esrrga:ENSDARG00000004861, shox2:ENSDARG00000075713 
# pvalb6:ENSDARG00000009311, meis3:ENSDARG00000002795, shisa9a:ENSDARG00000045145, angptl3:ENSDARG00000044365, kif26ab:ENSDARG00000015016, prdm1a:ENSDARG00000002445, b3glcta:ENSDARG00000073917, slc6a1a:ENSDARG00000045944, si:ch211-222l21.1:ENSDARG00000076532, lrrtm1:ENSDARG00000052713 
# PC_ 2 
# Positive:  slc6a11b:ENSDARG00000087981, gpr37l1b:ENSDARG00000056774, s100v2:ENSDARG00000070702, sox19a:ENSDARG00000010770, slc7a10b:ENSDARG00000051730, si:ch211-251b21.1:ENSDARG00000007275, acsl2:ENSDARG00000078399, msrb2:ENSDARG00000018459, sc:d189:ENSDARG00000102858, ca8:ENSDARG00000039098 
# si:ch211-180a12.2:ENSDARG00000087857, hspb15:ENSDARG00000078411, entpd2b:ENSDARG00000044795, cblc:ENSDARG00000078217, zgc:101744:ENSDARG00000038694, gbx1:ENSDARG00000071418, f3a:ENSDARG00000099124, si:ch73-173h19.3:ENSDARG00000104210, slc25a38a:ENSDARG00000059805, mgst1.1:ENSDARG00000032618 
# s1pr1:ENSDARG00000042690, atic:ENSDARG00000016706, cnn3a:ENSDARG00000099510, tmem176:ENSDARG00000078659, rdh8a:ENSDARG00000028048, cyp1b1:ENSDARG00000068934, aqp1a.1:ENSDARG00000023713, slc25a48:ENSDARG00000021250, ppap2d:ENSDARG00000069940, prom1b:ENSDARG00000034007 
# Negative:  NA:ENSDARG00000096508, ccl34b.1:ENSDARG00000093608, ccl35.1:ENSDARG00000103466, apoc1:ENSDARG00000092170, si:busm1-266f07.2:ENSDARG00000031745, ncf1:ENSDARG00000033735, havcr2:ENSDARG00000091692, cxcr3.2:ENSDARG00000041041, fcer1gl:ENSDARG00000104077, si:ch211-226m16.2:ENSDARG00000036785 
# lgals9l1:ENSDARG00000025903, si:dkey-27i16.2:ENSDARG00000098293, si:dkey-206d17.12:ENSDARG00000101034, havcr1:ENSDARG00000040178, ms4a17a.10:ENSDARG00000095695, si:ch73-86n18.1:ENSDARG00000056379, arpc1b:ENSDARG00000027063, cd74b:ENSDARG00000036628, nckap1l:ENSDARG00000075748, ch25h:ENSDARG00000045190 
# mpeg1.1:ENSDARG00000055290, tspan36:ENSDARG00000024540, ctss2.2:ENSDARG00000013771, cd74a:ENSDARG00000009087, mhc2dab:ENSDARG00000079105, cd83:ENSDARG00000079553, spi1b:ENSDARG00000000767, vsir:ENSDARG00000068784, cmklr1:ENSDARG00000090890, lcp1:ENSDARG00000023188 
# PC_ 3 
# Positive:  pcdh17:ENSDARG00000027041, sulf2b:ENSDARG00000013838, gapdhs:ENSDARG00000039914, zfhx4:ENSDARG00000075542, nr5a2:ENSDARG00000100940, si:dkey-183i3.9:ENSDARG00000104326, olfm3a:ENSDARG00000071493, meis3:ENSDARG00000002795, zfhx3:ENSDARG00000103057, si:ch73-386h18.1:ENSDARG00000073944 
# ca8:ENSDARG00000039098, shox2:ENSDARG00000075713, sox19a:ENSDARG00000010770, eno1a:ENSDARG00000022456, s100v2:ENSDARG00000070702, zgc:101744:ENSDARG00000038694, pon2:ENSDARG00000016856, atp1b1b:ENSDARG00000076833, spon1b:ENSDARG00000023694, pcbp3:ENSDARG00000054378 
# slc25a38a:ENSDARG00000059805, f3a:ENSDARG00000099124, si:ch211-191a24.4:ENSDARG00000037528, eno1b:ENSDARG00000013750, ier2b:ENSDARG00000086881, dbx1a:ENSDARG00000086393, nefma:ENSDARG00000021351, islr2:ENSDARG00000051875, zgc:112437:ENSDARG00000009215, atic:ENSDARG00000016706 
# Negative:  nfixb:ENSDARG00000061836, slc18a3a:ENSDARG00000006356, chata:ENSDARG00000015854, nfixa:ENSDARG00000043226, neurod6a:ENSDARG00000040008, neurod1:ENSDARG00000019566, zeb2a:ENSDARG00000062338, slc5a7a:ENSDARG00000074860, rp1l1b:ENSDARG00000088377, slc6a5:ENSDARG00000067964 
# nrp2a:ENSDARG00000096546, hoxc3a:ENSDARG00000070339, aff2:ENSDARG00000052242, gad2:ENSDARG00000015537, ebf3a:ENSDARG00000100244, zfpm2a:ENSDARG00000040123, oprd1b:ENSDARG00000037159, tnc:ENSDARG00000021948, prkg1b:ENSDARG00000031702, postnb:ENSDARG00000104267 
# bhlhe22:ENSDARG00000058039, igsf9b:ENSDARG00000010408, mcama:ENSDARG00000089643, nrp1a:ENSDARG00000102153, zbtb16a:ENSDARG00000007184, cdh23:ENSDARG00000007561, neurod6b:ENSDARG00000020794, NA:ENSDARG00000075022, arhgap5:ENSDARG00000061294, slc32a1:ENSDARG00000059775 
# PC_ 4 
# Positive:  zfhx4:ENSDARG00000075542, zfhx3:ENSDARG00000103057, si:ch211-222l21.1:ENSDARG00000076532, si:ch73-386h18.1:ENSDARG00000073944, pcdh19:ENSDARG00000034344, nr5a2:ENSDARG00000100940, sulf2b:ENSDARG00000013838, pcdh17:ENSDARG00000027041, nxph1:ENSDARG00000033447, hapln1a:ENSDARG00000089769 
# prdm1a:ENSDARG00000002445, spon1b:ENSDARG00000023694, tubb5:ENSDARG00000037997, b3glcta:ENSDARG00000073917, dlb:ENSDARG00000004232, si:dkey-35i13.1:ENSDARG00000097648, olfm3a:ENSDARG00000071493, sema5a:ENSDARG00000058821, shox2:ENSDARG00000075713, hhip:ENSDARG00000060397 
# si:dkey-183i3.9:ENSDARG00000104326, gpc2:ENSDARG00000104217, ebf2:ENSDARG00000042525, fsta:ENSDARG00000052846, NA:ENSDARG00000013312, igfn1.3:ENSDARG00000097613, kif21b:ENSDARG00000009733, htr1aa:ENSDARG00000093745, rasd4:ENSDARG00000099238, sall1a:ENSDARG00000074319 
# Negative:  slc5a7a:ENSDARG00000074860, chata:ENSDARG00000015854, calb2a:ENSDARG00000041062, sema3aa:ENSDARG00000019235, aclya:ENSDARG00000099079, slc18a3a:ENSDARG00000006356, nfixb:ENSDARG00000061836, ndnfl:ENSDARG00000076462, megf11:ENSDARG00000062686, rgs5b:ENSDARG00000017860 
# neurod6a:ENSDARG00000040008, gbx2:ENSDARG00000002933, fgfr2:ENSDARG00000058115, chgb:ENSDARG00000076500, rxfp3.2b:ENSDARG00000061846, sema3ab:ENSDARG00000042210, eno1a:ENSDARG00000022456, gadd45ba:ENSDARG00000027744, neurod1:ENSDARG00000019566, junba:ENSDARG00000074378 
# oprd1b:ENSDARG00000037159, hoxa10b:ENSDARG00000031337, grm8b:ENSDARG00000076508, ptprz1b:ENSDARG00000020871, igsf21b:ENSDARG00000056084, slc27a1a:ENSDARG00000006240, si:dkey-178e17.1:ENSDARG00000080000, ier2a:ENSDARG00000099195, gapdhs:ENSDARG00000039914, cbln2b:ENSDARG00000077151 
# PC_ 5 
# Positive:  tubb5:ENSDARG00000037997, nfixa:ENSDARG00000043226, hapln1a:ENSDARG00000089769, ebf2:ENSDARG00000042525, dlb:ENSDARG00000004232, cd99l2:ENSDARG00000056722, neurod4:ENSDARG00000003469, hoxb8a:ENSDARG00000056027, myo1b:ENSDARG00000024694, her4.5:ENSDARG00000056729 
# nhlh2:ENSDARG00000025495, nfixb:ENSDARG00000061836, her15.2:ENSDARG00000054560, cdk6:ENSDARG00000070228, olig4:ENSDARG00000052610, her15.1:ENSDARG00000054562, spegb:ENSDARG00000009567, her4.2:ENSDARG00000094426, rfx4:ENSDARG00000026395, adcyap1b:ENSDARG00000027740 
# stap2a:ENSDARG00000092810, si:ch211-194k22.8:ENSDARG00000096711, phc2b:ENSDARG00000013224, si:dkey-56m19.5:ENSDARG00000068432, crip2:ENSDARG00000070670, igsf9bb:ENSDARG00000069467, zbtb16a:ENSDARG00000007184, flrt3:ENSDARG00000076895, zgc:92242:ENSDARG00000029443, chd7:ENSDARG00000075211 
# Negative:  fosab:ENSDARG00000031683, eno1a:ENSDARG00000022456, atp1b1b:ENSDARG00000076833, igfn1.3:ENSDARG00000097613, ier2b:ENSDARG00000086881, egr4:ENSDARG00000077799, gadd45ba:ENSDARG00000027744, rxfp2b:ENSDARG00000019660, trpc2b:ENSDARG00000003344, egr1:ENSDARG00000037421 
# si:dkey-185m8.2:ENSDARG00000099026, gapdhs:ENSDARG00000039914, gpat3:ENSDARG00000016048, ier2a:ENSDARG00000099195, nlrp3:ENSDARG00000078620, fasn2:ENSDARG00000096762, muc5f:ENSDARG00000104914, junbb:ENSDARG00000104773, shisa8b:ENSDARG00000088204, ache:ENSDARG00000031796 
# gjb1a:ENSDARG00000035553, NA:ENSDARG00000075022, si:dkey-39a18.1:ENSDARG00000045669, mertka:ENSDARG00000074695, igsf9b:ENSDARG00000010408, ndnfl:ENSDARG00000076462, nes:ENSDARG00000088805, mmrn2a:ENSDARG00000076135, si:dkey-204a24.11:ENSDARG00000058685, il10:ENSDARG00000078147 
# Plot PCAs
png('12.PCAplot.PC1vsPC2.vst.plate.counts.png')
PCAPlot(object = seu_c.filt, dims=c(1,2),group.by="plate")
dev.off()

ElbowPlot(seu_c.filt)
#Jackstraw calculates which ones are the significant PCs.
seu_c.filt <- JackStraw(seu_c.filt)
seu_c.filt <- ScoreJackStraw(seu_c.filt, dims = 1:20)
pcs <- which(JS(object = seu_c.filt[['pca']], slot = 'overall')[, 2] < 1e-3)
pcs
# [1]  1  2  3  4  5 10 12 13 14 18

seu_c.filt <- FindNeighbors(seu_c.filt, reduction = "pca", dims = 1:9)
seu_c.filt <- FindClusters(seu_c.filt, resolution = 0.5)
table(Idents(seu_c.filt))
# 0   1   2   3   4 
# 221 137 125  77  33

# Check clusters over resolutions
res <- c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4)
tmp=seu_c.filt
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:9)
tmp <- FindClusters(tmp, resolution = res, print.output = 0)
pdf('19.Clustree.V2a.excribo.pdf', height = 10, width = 10)
clustree(tmp)
dev.off()

# Run UMAP and tSNE
seu_c.filt <- RunUMAP(seu_c.filt, dims = 1:9)
seu_c.filt <- RunTSNE(seu_c.filt, dims = 1:9)

pdf('18.TSNE_Seurat.clusters.pdf')
DimPlot(seu_c.filt, reduction = "tsne")
dev.off()

pdf('18.UMAP_Seurat.clusters.pdf')
DimPlot(seu_c.filt, reduction = "umap")
dev.off()

p1 <- DimPlot(seu_c.filt, reduction = "umap", group.by = "plate")
p2 <- DimPlot(seu_c.filt, reduction = "umap")
p3 <- DimPlot(seu_c.filt, reduction = "tsne", group.by = "plate")
p4 <- DimPlot(seu_c.filt, reduction = "tsne")
multiplot(p1, p2, p3, p4, cols = 2)

p1 <- DimPlot(seu_c.filt, reduction = "tsne", group.by = "ident", pt=2)
p2 <- DimPlot(seu_c.filt, reduction = "tsne", pt=2)
png('18b.TSNE_Seurat.clusters.p353.354.counts.png', width=480, height=800)
multiplot(p1, p2, cols = 1)
dev.off()

pdf('18b.TSNE_Seurat.clusters.p353.354.counts.pdf', width=10, height=10)
DimPlot(seu_c.filt, reduction = "tsne", pt=2)
dev.off()

# Differential expression analysis and glia exclusion ####
V2a_incl.glia.markers <- FindAllMarkers(seu_c.filt, test.use= "MAST", min.pct=0.25, logfc.threshold=0.3) 
V2a_incl.glia.markers.pos <- FindAllMarkers(seu_c.filt, test.use = "MAST", min.pct=0.25, logfc.threshold=0.3, only.pos = TRUE)

# Write list to file -> This is the differential gene expression analysis used in the paper
write.table(V2a_incl.glia.markers %>% group_by(cluster), file="AllMarkerGenes.V2a.incl.glia.txt",quote=FALSE,row.names=FALSE)
write.table(V2a_incl.glia.markers.pos %>% group_by(cluster), file="AllMarkerGenes.onlyPos.V2a.incl.glia.txt",quote=FALSE,row.names=FALSE)

# Heatmap top10 marker genes (sorted by padj) for each cluster
feat = V2a_incl.glia.markers %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('20.heatmap.markerGenes.top10padj.counts.p353.354.png', width=920, height=1200)
DoHeatmap(seu_c.filt, features = feat$gene,  slot="data")
dev.off()

feat.pos = V2a_incl.glia.markers.pos %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('21.heatmap.markerGenes.top10padj.onlyPos.counts.p353.354.png', width=920, height=1200)
DoHeatmap(seu_c.filt, features = feat.pos$gene,  slot="data")
dev.off()


#plot glia genes to exclude cluster 7
glia.markers<-read.table("glia.txt")
glia.markers <- lapply(glia.markers, as.character)$V1 
DoHeatmap(seu_c.filt, features = glia.markers, slot="data")
pdf('glia.markers.pdf', width=10, height=10)
DotPlot(seu_c.filt, features = glia.markers, dot.scale = 7) + RotatedAxis() + coord_flip()
dev.off()

#exclude cluster 4 glia
V2a <- subset(x = seu_c.filt, subset = seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters ==3)

#get all gene names in dataset (for gene ontology and other things):
all.genes <-rownames(V2a@assays$RNA)
write(all.genes, file="all.genes.V2a.txt")

#rename and order clusters 
V2a <- RenameIdents(V2a, '1' = 'Mito', '2' = 'Shox2', '0' = 'VaChTa', '3' = 'Ribo')

p1 <- DimPlot(V2a, reduction = "umap", group.by = "plate", pt=2)
p2 <- DimPlot(V2a, reduction = "umap", pt=2)
png('18c.UMAP_V2a.only.png', width=480, height=800)
multiplot(p1, p2, cols = 1)
dev.off()

pdf('18c.UMAP_V2a.only.pdf', width=10, height=10)
DimPlot(V2a, reduction = "umap", pt=2)
dev.off()

#Extract counts of RNAscope probe genes for analysis outside R (Figure 5p)
probes <- c("tubb5:ENSDARG00000037997", "esrrga:ENSDARG00000004861", "pvalb6:ENSDARG00000009311", "shox2:ENSDARG00000075713","neurod1:ENSDARG00000019566", "slc18a3a:ENSDARG00000006356")
validation.counts<- FetchData(V2a, vars = probes, slot = 'counts')
head(validation.counts)
cluster<- FetchData(V2a, vars = 'seurat_clusters')
cluster<-as.character(cluster$seurat_clusters)
validation.counts$cluster<- cluster
write.csv(validation.counts, "validation.counts.V2a.csv")


#save object as rds for integration with MNs
saveRDS(V2a, "V2a.rds")

#extract cluster ID for harmony
write.table(V2a_4clusters@active.ident, file = "V2a_clusterID.txt")


#Cluster by TF####

setwd("V2a/TF")

tfs<-read.table("TF V2a.txt")
tfs <- lapply(tfs, as.character)$V1 

V2a.tfac <- V2a
#we decided to exclude the immature V2a from this analysis
V2a.tfac <- subset(x = V2a.tfac, subset = seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 2)
V2a.tfac <- ScaleData(V2a.tfac, features = tfs, vars.to.regress = c("nFeature_RNA", "nCounts_RNA", "percent_ribo"))

# PCA on transcription factor genes only
V2a.tfac <- RunPCA(V2a.tfac, features = tfs, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
# PC_ 1 
# Positive:  nfixb:ENSDARG00000061836, neurod6a:ENSDARG00000040008, nfia:ENSDARG00000062420, neurod1:ENSDARG00000019566, vsx2:ENSDARG00000005574, aff2:ENSDARG00000052242, zeb2a:ENSDARG00000062338, nfixa:ENSDARG00000043226, hoxb3a:ENSDARG00000029263, gbx2:ENSDARG00000002933 
# hoxc3a:ENSDARG00000070339, neurod6b:ENSDARG00000020794, hoxc9a:ENSDARG00000092809, zfpm2a:ENSDARG00000040123, bhlhe22:ENSDARG00000058039, yaf2:ENSDARG00000015966, hoxb7a:ENSDARG00000056030, hoxb8a:ENSDARG00000056027, camta1a:ENSDARG00000077428, prdm8b:ENSDARG00000054683 
# twist1a:ENSDARG00000030402, hoxa10b:ENSDARG00000031337, adcyap1b:ENSDARG00000027740, dachc:ENSDARG00000003142, nhlh2:ENSDARG00000025495, prox3:ENSDARG00000088810, hoxb9a:ENSDARG00000056023, ebf3a:ENSDARG00000100244, hoxc8a:ENSDARG00000070346, hoxc6a:ENSDARG00000070343 
# Negative:  zfhx4:ENSDARG00000075542, zfhx3:ENSDARG00000103057, nr5a2:ENSDARG00000100940, si:ch73-386h18.1:ENSDARG00000073944, prdm1a:ENSDARG00000002445, tox:ENSDARG00000032317, shox2:ENSDARG00000075713, meis2a:ENSDARG00000098240, onecut1:ENSDARG00000007982, pou3f1:ENSDARG00000009823 
# nkx6.2:ENSDARG00000104735, nr5a1a:ENSDARG00000103176, pou2f2a:ENSDARG00000019658, esrrb:ENSDARG00000100847, meis3:ENSDARG00000002795, id4:ENSDARG00000045131, mkxb:ENSDARG00000114182, pcbp3:ENSDARG00000054378, tox2:ENSDARG00000100055, bhlhe41:ENSDARG00000041691 
# ahrrb:ENSDARG00000052618, meis2b:ENSDARG00000077840, znf503:ENSDARG00000018492, scxa:ENSDARG00000035695, nr2f1a:ENSDARG00000052695, pcbp4:ENSDARG00000024276, hoxb13a:ENSDARG00000056015, myca:ENSDARG00000045695, nr0b1:ENSDARG00000056541, jupb:ENSDARG00000059067 
# PC_ 2 
# Positive:  zbtb16a:ENSDARG00000007184, nfixa:ENSDARG00000043226, nr5a2:ENSDARG00000100940, prdm1a:ENSDARG00000002445, zfhx3:ENSDARG00000103057, hoxc3a:ENSDARG00000070339, hoxc6a:ENSDARG00000070343, onecut1:ENSDARG00000007982, hoxc1a:ENSDARG00000070337, hoxc8a:ENSDARG00000070346 
# zfhx4:ENSDARG00000075542, sox5:ENSDARG00000011582, nkx6.2:ENSDARG00000104735, nr5a1a:ENSDARG00000103176, hoxb8a:ENSDARG00000056027, klf6a:ENSDARG00000029072, jupa:ENSDARG00000070787, id4:ENSDARG00000045131, hoxb3a:ENSDARG00000029263, irx2a:ENSDARG00000001785 
# gli2a:ENSDARG00000025641, stat4:ENSDARG00000028731, ftr54:ENSDARG00000043850, nfat5b:ENSDARG00000102556, mecom:ENSDARG00000060808, si:ch211-155e24.3:ENSDARG00000105127, hoxc4a:ENSDARG00000070338, hoxc5a:ENSDARG00000070340, zeb2a:ENSDARG00000062338, pou6f2:ENSDARG00000086362 
# Negative:  dacha:ENSDARG00000010132, esrrga:ENSDARG00000004861, pbx3b:ENSDARG00000013615, tefa:ENSDARG00000039117, cited4a:ENSDARG00000035990, jdp2b:ENSDARG00000020133, sp8a:ENSDARG00000011870, tob1a:ENSDARG00000032619, bhlhe40:ENSDARG00000004060, fosb:ENSDARG00000055751 
# dachd:ENSDARG00000069440, meis1b:ENSDARG00000012078, rps6ka1:ENSDARG00000033437, chchd2:ENSDARG00000059304, vsx2:ENSDARG00000005574, junbb:ENSDARG00000104773, cbx7a:ENSDARG00000038025, jund:ENSDARG00000067850, meis1a:ENSDARG00000002937, nfil3-6:ENSDARG00000087188 
# egr4:ENSDARG00000077799, fosab:ENSDARG00000031683, hnrnpk:ENSDARG00000018914, calcoco1b:ENSDARG00000016391, pknox2:ENSDARG00000055349, dachc:ENSDARG00000003142, klf13:ENSDARG00000061368, sox11b:ENSDARG00000095743, lhx4:ENSDARG00000039458, ppargc1b:ENSDARG00000101569 
# PC_ 3 
# Positive:  hoxc11a:ENSDARG00000070351, hoxc12a:ENSDARG00000070352, hoxc12b:ENSDARG00000103133, hoxa11b:ENSDARG00000007009, hoxa10b:ENSDARG00000031337, hoxd11a:ENSDARG00000059267, hoxd10a:ENSDARG00000057859, hoxc10a:ENSDARG00000070348, hoxa9b:ENSDARG00000056819, hoxd12a:ENSDARG00000059263 
# hoxc11b:ENSDARG00000102631, hoxc13a:ENSDARG00000070353, hoxa11a:ENSDARG00000104162, nkx6.2:ENSDARG00000104735, hoxd13a:ENSDARG00000059256, insm1a:ENSDARG00000091756, zeb2a:ENSDARG00000062338, rybpb:ENSDARG00000053459, yaf2:ENSDARG00000015966, sox21b:ENSDARG00000008540 
# onecut1:ENSDARG00000007982, znf536:ENSDARG00000103648, hoxa9a:ENSDARG00000105013, zfpm2a:ENSDARG00000040123, pou2f1b:ENSDARG00000007996, hoxd9a:ENSDARG00000059274, rorb:ENSDARG00000033498, pou2f2a:ENSDARG00000019658, chd7:ENSDARG00000075211, pou3f3b:ENSDARG00000095896 
# Negative:  hoxb8a:ENSDARG00000056027, pbx3b:ENSDARG00000013615, hoxc8a:ENSDARG00000070346, hoxb6b:ENSDARG00000026513, meis1b:ENSDARG00000012078, sp8a:ENSDARG00000011870, meis2b:ENSDARG00000077840, npas1:ENSDARG00000015876, si:dkey-68o6.5:ENSDARG00000099247, adcyap1b:ENSDARG00000027740 
# npas3:ENSDARG00000109820, myca:ENSDARG00000045695, hoxc6a:ENSDARG00000070343, znf703:ENSDARG00000035563, meis2a:ENSDARG00000098240, nrip1b:ENSDARG00000068894, mrtfab:ENSDARG00000076229, hoxb5b:ENSDARG00000054030, tle3a:ENSDARG00000002787, vgll2a:ENSDARG00000041706 
# lmx1al:ENSDARG00000077915, mef2aa:ENSDARG00000031756, pax3b:ENSDARG00000028348, klf7b:ENSDARG00000043821, emx3:ENSDARG00000020417, mecom:ENSDARG00000060808, fosab:ENSDARG00000031683, atf5a:ENSDARG00000068096, hoxb8b:ENSDARG00000054025, six7:ENSDARG00000070107 
# PC_ 4 
# Positive:  gli2a:ENSDARG00000025641, stat4:ENSDARG00000028731, sox5:ENSDARG00000011582, hoxc12b:ENSDARG00000103133, hoxc12a:ENSDARG00000070352, six7:ENSDARG00000070107, etv4:ENSDARG00000018303, ftr54:ENSDARG00000043850, znf366:ENSDARG00000040116, emx3:ENSDARG00000020417 
# lmx1al:ENSDARG00000077915, ftr85:ENSDARG00000034707, vgll2a:ENSDARG00000041706, si:dkey-68o6.5:ENSDARG00000099247, hoxc11a:ENSDARG00000070351, nfat5b:ENSDARG00000102556, dmrt1:ENSDARG00000007349, hoxc11b:ENSDARG00000102631, pax3b:ENSDARG00000028348, esrrgb:ENSDARG00000011696 
# six1a:ENSDARG00000039304, dlx5a:ENSDARG00000042296, fhl3a:ENSDARG00000034643, ftr97:ENSDARG00000055436, hoxc13a:ENSDARG00000070353, foxc1b:ENSDARG00000055398, foxp3b:ENSDARG00000078279, hoxc13b:ENSDARG00000113877, lef1:ENSDARG00000031894, mecom:ENSDARG00000060808 
# Negative:  hoxc9a:ENSDARG00000092809, hoxc8a:ENSDARG00000070346, hoxb8a:ENSDARG00000056027, hoxc3a:ENSDARG00000070339, hoxc6a:ENSDARG00000070343, zbtb16a:ENSDARG00000007184, hoxb7a:ENSDARG00000056030, hoxb3a:ENSDARG00000029263, hoxc1a:ENSDARG00000070337, pou3f1:ENSDARG00000009823 
# nfixa:ENSDARG00000043226, twist1a:ENSDARG00000030402, ahrrb:ENSDARG00000052618, smarcc2:ENSDARG00000077946, aff2:ENSDARG00000052242, sp8b:ENSDARG00000056666, neurod1:ENSDARG00000019566, si:dkeyp-72h1.1:ENSDARG00000095347, hoxb6b:ENSDARG00000026513, hoxb6a:ENSDARG00000010630 
# si:dkey-89b17.4:ENSDARG00000075545, jmjd1cb:ENSDARG00000079939, hoxc4a:ENSDARG00000070338, drap1:ENSDARG00000041203, shox2:ENSDARG00000075713, hoxb10a:ENSDARG00000011579, pou3f2b:ENSDARG00000076262, si:ch73-138e16.4:ENSDARG00000102363, klf6a:ENSDARG00000029072, adcyap1b:ENSDARG00000027740 
# PC_ 5 
# Positive:  hoxc13a:ENSDARG00000070353, hoxc13b:ENSDARG00000113877, hoxb13a:ENSDARG00000056015, hoxa13b:ENSDARG00000036254, hmgb3a:ENSDARG00000056725, ftr66:ENSDARG00000005327, camta1b:ENSDARG00000007824, si:ch211-232b12.5:ENSDARG00000087186, foxi3b:ENSDARG00000009550, znf703:ENSDARG00000035563 
# myt1la:ENSDARG00000008209, lmx1a:ENSDARG00000020354, notch1b:ENSDARG00000052094, dnmt3ab:ENSDARG00000015566, si:cabz01054394.7:ENSDARG00000089715, smarca4a:ENSDARG00000077226, tshz2:ENSDARG00000079201, epc1b:ENSDARG00000060054, brca2:ENSDARG00000079015, atf5a:ENSDARG00000068096 
# foxo3b:ENSDARG00000042904, irf1a:ENSDARG00000043492, mrgbp:ENSDARG00000028894, sox11a:ENSDARG00000077811, tbx19:ENSDARG00000079187, nfil3-5:ENSDARG00000094965, tcf7l2:ENSDARG00000004415, hoxc12a:ENSDARG00000070352, maml1:ENSDARG00000076466, nr1d2a:ENSDARG00000003820 
# Negative:  npas4a:ENSDARG00000055752, junba:ENSDARG00000074378, hoxd10a:ENSDARG00000057859, egr4:ENSDARG00000077799, zeb2a:ENSDARG00000062338, pou3f1:ENSDARG00000009823, hoxd9a:ENSDARG00000059274, hoxc10a:ENSDARG00000070348, hoxc9a:ENSDARG00000092809, znf503:ENSDARG00000018492 
# cbx7a:ENSDARG00000038025, hoxc3a:ENSDARG00000070339, nr4a1:ENSDARG00000000796, egr1:ENSDARG00000037421, fosaa:ENSDARG00000040135, xbp1:ENSDARG00000035622, tob1a:ENSDARG00000032619, cbx7b:ENSDARG00000087181, zgpat:ENSDARG00000027403, rbm8a:ENSDARG00000016516 
# hoxa10b:ENSDARG00000031337, junbb:ENSDARG00000104773, ftr38:ENSDARG00000095678, maff:ENSDARG00000028957, tle2a:ENSDARG00000042484, nr5a2:ENSDARG00000100940, pcgf5a:ENSDARG00000071007, fosab:ENSDARG00000031683, phf19:ENSDARG00000078050, stat4:ENSDARG00000028731 

# Plot PCAs
png('Tfac.PCAplot.PC1vsPC2.vst.plate.counts.png')
PCAPlot(object = V2a.tfac, dims=c(1,2),group.by="plate")
dev.off()

ElbowPlot(V2a.tfac)

V2a.tfac <- JackStraw(V2a.tfac)
V2a.tfac <- ScoreJackStraw(V2a.tfac, dims = 1:20)
pcs <- which(JS(object = V2a.tfac[['pca']], slot = 'overall')[, 2] < 1e-3)
pcs
# [1]  1  2  3  4  7  9 

# Cluster with Seurat
# DECIDE HOW MANY PCs for clustering (from elbow plot AND JACKSTRAW)
V2a.tfac <- FindNeighbors(V2a.tfac, reduction = "pca", dims = pcs)
V2a.tfac <- FindClusters(V2a.tfac, resolution = 0.5)
table(Idents(V2a.tfac))

#pcs 0.5
# 0   1   2   3 
# 148 124 123  88  

# Check clusters over resolutions
res <- c(0.4,0.5,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4)
tmp=V2a.tfac
tmp <- FindClusters(tmp, resolution = res, print.output = 0)
png('Tfac.clustree.v2a.png')
clustree(tmp)
dev.off()

# Run UMAP and tSNE
V2a.tfac <- RunUMAP(V2a.tfac, dims = pcs)
V2a.tfac <- RunTSNE(V2a.tfac, dims = pcs)
p1 <- DimPlot(V2a.tfac, reduction = "umap", group.by = "plate")
p2 <- DimPlot(V2a.tfac, reduction = "umap")
p3 <- DimPlot(V2a.tfac, reduction = "tsne", group.by = "plate")
p4 <- DimPlot(V2a.tfac, reduction = "tsne")
multiplot(p1, p2, p3,p4, cols = 2)

#plot new tSNE with old clusters mn
p1 <- DimPlot(V2a.tfac, reduction = "umap", group.by = "plate", pt=2)
p2 <- DimPlot(V2a.tfac, reduction = "umap", pt=2)

png('Tfac.umap_V2a.png', width=480, height=800)
multiplot(p1, p2, cols = 1)
dev.off()

#plot on normal V2a umap
V2a <- AddMetaData(V2a, V2a.tfac@meta.data["seurat_clusters"], col.name = "tfac_clusters")
pdf('TFac.UMAP_V2a4clusters.pdf')
DimPlot(V2a, reduction = "umap", group.by = "tfac_clusters")+xlim(-7.5,5)+ylim(-6,6)
dev.off()

# Differential expression analisys
V2a.tfac.markers <- FindAllMarkers(V2a.tfac, test.use= "MAST", min.pct=0.25, logfc.threshold=0.3) 
V2a.tfac.markers.pos <- FindAllMarkers(V2a.tfac, test.use = "MAST", min.pct=0.25, logfc.threshold=0.3, only.pos = TRUE)


# Write list to file
write.table(V2a.tfac.markers %>% group_by(cluster), file="Tfac.allMarkerGenes.txt",quote=FALSE,row.names=FALSE)
write.table(V2a.tfac.markers.pos %>% group_by(cluster), file="Tfac.allMarkerGenes.onlyPos.txt",quote=FALSE,row.names=FALSE)
#the list of differentially expressed TFs has been made starting from this and filtering with the list of all tfs (TF NEW.txt)

# Heatmap top10 marker genes (sorted by padj) for each cluster
feat = V2a.tfac.markers %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('Tfac.heatmap.markerGenes.top10padj.png', width=920, height=1200)
DoHeatmap(V2a.tfac, features = feat$gene,  slot="data")
dev.off()

feat.pos = V2a.tfac.markers.pos %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj) %>% arrange(cluster,p_val_adj)
png('Tfac.heatmap.markerGenes.top10padj.onlyPos.png', width=920, height=1200)
DoHeatmap(V2a.tfac, features = feat.pos$gene,  slot="data")
dev.off()


#plots for figures####

setwd("V2a/figures")

#markers per cluster ig 4c
Fig_V2a_cluster_markers <- read.table("V2a.markers.per.cluster.txt")
Fig_V2a_cluster_markers <- lapply(Fig_V2a_cluster_markers, as.character)$V1 

pdf('Dotplot.V2a.markers per clusters.pdf',  width=8, height=9)
DotPlot(V2a, features = Fig_V2a_cluster_markers,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()

#UMAP Fig 4b
pdf('UMAP_V2a.only.pdf', width=10, height=10)
DimPlot(V2a, reduction = "umap", pt=2)
dev.off()

#interesting genes and markers Fig 4e-g
MarkersFig4 <- c("tubb5:ENSDARG00000037997", "tuba1a:ENSDARG00000001889", "nme2b.1:ENSDARG00000103791", "cnp:ENSDARG00000070822", 
                "slc17a6a:ENSDARG00000001127", "slc17a6b:ENSDARG00000041150",
                "esrrga:ENSDARG00000004861", "shox2:ENSDARG00000075713", "slc18a3a:ENSDARG00000006356", "chata:ENSDARG00000015854")


pdf("MarkersFig4 Violin plots.pdf", width = 7, height = 10)
StackedVlnPlot(V2a, features = MarkersFig4)
dev.off()

#validation genes Fig.5a, d, g
pdf('esrrga featureplot.pdf')
FeaturePlot(V2a, features = "esrrga:ENSDARG00000004861",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-7.5,5)+ylim(-6,6)
dev.off()

pdf('shox2 featureplot.pdf')
FeaturePlot(V2a, features = "shox2:ENSDARG00000075713",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-7.5,5)+ylim(-6,6)
dev.off()

pdf('vachta featureplot.pdf')
FeaturePlot(V2a, features = "slc18a3a:ENSDARG00000006356",  reduction = "umap", pt=2, cols = c("#E0E0E0", "#0066A6"))+xlim(-7.5,5)+ylim(-6,6)
dev.off()

#V2a trascription factor clusters Ext Data Fig.6a
pdf('TFac.UMAP_V2a.pdf')
DimPlot(V2a, reduction = "umap", group.by = "tfac_clusters")+xlim(-7.5,5)+ylim(-6,6)
dev.off()

#dotplot DE TFs Ext Data Fig.6b
top10tfactors<-read.table("V2a.TF_DE_first10.txt")
top10tfactors <- lapply(top10tfactors, as.character)$V1 

pdf('Tfac.tfactors.top10padj.onlyPos.pdf')
DotPlot(V2a.tfac, features = top10tfactors,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()
#cluster 0 is ribo, to be removed

#plot selected best slow/fast for MN AND V2a Fig. 6b
BestSF<-read.table("Slow_Fast_combined_MNandV2a.txt")
BestSF<- lapply(BestSF, as.character)$V1
DoHeatmap(V2a, features = BestSF, slot="data")
pdf('BestSF_V2a.pdf')
DotPlot(V2a, features = BestSF,  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()
