
#Harmony intergation ####
# Load R packages/libraries

#install.packages("harmony")
library(Seurat)
library(cowplot)
library(harmony)
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

setwd("Integration MN_V2a")

#extract count tables from analyzed objects to know which cells to keep from the raw data
#also ocreate table with neuron identity to add as MetaData

MN <- readRDS("MN.rds")
V2aexc6 <- readRDS("V2a.rds")

MNcounts <- MN@assays$RNA@counts
V2acounts <- V2aexc6@assays$RNA@counts

MNcells<-MNcounts@Dimnames[[2]]
V2acells<-V2acounts@Dimnames[[2]]

CellsToKeep <- c(MNcells, V2acells)
length(CellsToKeep)
#[1] 876

write.table(MN@active.ident, file = "MN_clusterID.txt")
write.table(V2a@active.ident, file = "V2a_clusterID.txt")
#manually merge the two files into "metadata.txt", also add cell type and column name and remove "


#start from raw counts
counts1 <- read.table("C:/Users/irene/Desktop/New analysis/358.354.353.merge.rpkmforgenes_counts.csv", header=TRUE, sep=',', row.names='gene')
counts2 <- read.table("C:/Users/irene/Desktop/New analysis/174.merge.rpkmforgenes_counts.csv", header=TRUE, sep=',', row.names='gene')
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

#create object and add MetaData 
#Saving before Harmony to compare later
MN_V2a <- CreateSeuratObject(counts = cbind(counts1, counts2), project = "HarmonyInteg", min.cells = 3) 

#keep only the cells I have in final objects from prev analysis
length(colnames(MN_V2a))
#[1] 1536

MN_V2a <- MN_V2a[,!colnames(MN_V2a) %in% CellsToKeep, invert = TRUE]

length(colnames(MN_V2a))
#[1] 876

metadata <- read.table("Metadata.txt", header=TRUE, sep='\t', row.names='Sample')
MN_V2a <- AddMetaData(MN_V2a, metadata$orig.cluster, col.name = "orig.cluster")
MN_V2a <- AddMetaData(MN_V2a, metadata$celltype, col.name = "celltype")

#normalize etc.
scale.factor <- mean(colSums(GetAssayData(MN_V2a, slot = "counts")))
scale.factor
#[1] 394901.5
MN_V2a <- NormalizeData(object = MN_V2a, normalization.method = "LogNormalize", scale.factor = scale.factor)

# Feature selection
MN_V2a <- FindVariableFeatures(MN_V2a, selection.method = "vst", nfeatures = 2000)

head(VariableFeatures(MN_V2a), 10)
# [1] "chgb:ENSDARG00000076500"    "her15.2:ENSDARG00000054560" "bmp16:ENSDARG00000103679"   "tubb5:ENSDARG00000037997"  
# [5] "foxp4:ENSDARG00000076120"   "tac1:ENSDARG00000014490"    "hoxb13a:ENSDARG00000056015" "NA:ENSDARG00000104418"     
# [9] "angptl3:ENSDARG00000044365" "tesca:ENSDARG00000028346"  
MN_V2a <- ScaleData(MN_V2a)
MN_V2a <- RunPCA(MN_V2a, features = VariableFeatures(MN_V2a), do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# PC_ 1 
# Positive:  isl1:ENSDARG00000004023, tubb5:ENSDARG00000037997, mafbb:ENSDARG00000070542, vim:ENSDARG00000010008, si:dkey-56m19.5:ENSDARG00000068432, slit3:ENSDARG00000034268, myo10l3:ENSDARG00000074143, rassf1:ENSDARG00000004840, si:ch211-222l21.1:ENSDARG00000076532, tuba8l3:ENSDARG00000070155 
# tuba1a:ENSDARG00000001889, lima1a:ENSDARG00000101441, klf7b:ENSDARG00000043821, cnp:ENSDARG00000070822, shtn3:ENSDARG00000056519, NA:ENSDARG00000099464, prph:ENSDARG00000028306, myo1b:ENSDARG00000024694, alcama:ENSDARG00000026531, dlb:ENSDARG00000004232 
# isl2a:ENSDARG00000003971, si:ch211-232b12.5:ENSDARG00000087186, ppp1r14ba:ENSDARG00000044541, lmo4b:ENSDARG00000054749, eef1a1l1:ENSDARG00000020850, mnx1:ENSDARG00000035984, ebf2:ENSDARG00000042525, chd7:ENSDARG00000075211, stard13b:ENSDARG00000098954, olig4:ENSDARG00000052610 
# Negative:  vsx2:ENSDARG00000005574, arhgap5:ENSDARG00000061294, NA:ERCC-00002, NA:ERCC-00074, ryr3:ENSDARG00000071331, NA:ERCC-00130, npas4a:ENSDARG00000055752, gbx2:ENSDARG00000002933, egr4:ENSDARG00000077799, sema3aa:ENSDARG00000019235 
# chgb:ENSDARG00000076500, bhlhe22:ENSDARG00000058039, gadd45ba:ENSDARG00000027744, mbpa:ENSDARG00000036186, fxyd7:ENSDARG00000097110, NA:ERCC-00108, ier2a:ENSDARG00000099195, calb2a:ENSDARG00000041062, NA:ERCC-00009, atp1b1b:ENSDARG00000076833 
# grm8b:ENSDARG00000076508, egr1:ENSDARG00000037421, sema3ab:ENSDARG00000042210, si:ch211-199f5.1:ENSDARG00000102185, aff2:ENSDARG00000052242, ptprua:ENSDARG00000092638, slc6a1a:ENSDARG00000045944, oprd1b:ENSDARG00000037159, sema3fb:ENSDARG00000055373, histh1l1:ENSDARG00000035519 
# PC_ 2 
# Positive:  nefma:ENSDARG00000021351, pcbp3:ENSDARG00000054378, bmp16:ENSDARG00000103679, atp1b1b:ENSDARG00000076833, tac1:ENSDARG00000014490, pcdh17:ENSDARG00000027041, si:dkey-35i13.1:ENSDARG00000097648, scg2b:ENSDARG00000038574, zgc:65851:ENSDARG00000012281, atp2a3:ENSDARG00000060978 
# neflb:ENSDARG00000012426, tns1b:ENSDARG00000020845, pcdh19:ENSDARG00000034344, pvalb6:ENSDARG00000009311, sncga:ENSDARG00000034423, atp1a3b:ENSDARG00000104139, rspo2:ENSDARG00000079570, npas1:ENSDARG00000015876, agrn:ENSDARG00000079388, rgs3a:ENSDARG00000099746 
# anxa13l:ENSDARG00000013613, inab:ENSDARG00000053248, LOC100331987:ENSDARG00000100594, prph:ENSDARG00000028306, itga3a:ENSDARG00000037917, pcdh11:ENSDARG00000098652, fsta:ENSDARG00000052846, pdyn:ENSDARG00000087798, ttc27:ENSDARG00000007918, zfhx3:ENSDARG00000103057 
# Negative:  nfixb:ENSDARG00000061836, neurod6a:ENSDARG00000040008, nfixa:ENSDARG00000043226, neurod1:ENSDARG00000019566, nhlh2:ENSDARG00000025495, NA:pEGFP-N1, hoxb8a:ENSDARG00000056027, neurod6b:ENSDARG00000020794, zeb2a:ENSDARG00000062338, neurod4:ENSDARG00000003469 
# mcama:ENSDARG00000089643, nrp2a:ENSDARG00000096546, zfpm2a:ENSDARG00000040123, ebf2:ENSDARG00000042525, ebf3a:ENSDARG00000100244, hoxc3a:ENSDARG00000070339, bhlhe22:ENSDARG00000058039, dlb:ENSDARG00000004232, vsx2:ENSDARG00000005574, oprd1b:ENSDARG00000037159 
# dla:ENSDARG00000010791, cox4i2:ENSDARG00000022509, olig4:ENSDARG00000052610, insm1a:ENSDARG00000091756, hapln1a:ENSDARG00000089769, sulf2a:ENSDARG00000018423, myo1b:ENSDARG00000024694, hoxc9a:ENSDARG00000092809, gbx2:ENSDARG00000002933, ntsr1:ENSDARG00000077577 
# PC_ 3 
# Positive:  chata:ENSDARG00000015854, slc18a3a:ENSDARG00000006356, slc5a7a:ENSDARG00000074860, rgs5b:ENSDARG00000017860, ret:ENSDARG00000055305, tac1:ENSDARG00000014490, prph:ENSDARG00000028306, aff2:ENSDARG00000052242, esrrga:ENSDARG00000004861, pdyn:ENSDARG00000087798 
# sncga:ENSDARG00000034423, ddit3:ENSDARG00000059836, bmp16:ENSDARG00000103679, itga3a:ENSDARG00000037917, nrp1a:ENSDARG00000102153, sox11a:ENSDARG00000077811, zeb2a:ENSDARG00000062338, igsf21b:ENSDARG00000056084, spega:ENSDARG00000039256, oprd1b:ENSDARG00000037159 
# atp1b1a:ENSDARG00000013144, rspo2:ENSDARG00000079570, hunk:ENSDARG00000091317, postnb:ENSDARG00000104267, npas1:ENSDARG00000015876, prkg1b:ENSDARG00000031702, sema3ab:ENSDARG00000042210, dock1:ENSDARG00000099093, npr3:ENSDARG00000035253, kif26bb:ENSDARG00000024575 
# Negative:  nr5a2:ENSDARG00000100940, sulf2b:ENSDARG00000013838, zfhx4:ENSDARG00000075542, nxph1:ENSDARG00000033447, spon1b:ENSDARG00000023694, olfm3a:ENSDARG00000071493, shox2:ENSDARG00000075713, meis3:ENSDARG00000002795, prdm1a:ENSDARG00000002445, anos1a:ENSDARG00000012896 
# zfhx3:ENSDARG00000103057, epha4b:ENSDARG00000011600, ndnf:ENSDARG00000062936, b3glcta:ENSDARG00000073917, hapln4:ENSDARG00000018542, NA:ENSDARG00000013312, pcdh17:ENSDARG00000027041, kif26ab:ENSDARG00000015016, hhip:ENSDARG00000060397, lhx4:ENSDARG00000039458 
# LOC567790:ENSDARG00000094282, angptl3:ENSDARG00000044365, grm2a:ENSDARG00000004150, shdb:ENSDARG00000062109, btbd6b:ENSDARG00000032369, drd3:ENSDARG00000032131, slc6a1a:ENSDARG00000045944, htr1aa:ENSDARG00000093745, rdh10a:ENSDARG00000058730, calb1:ENSDARG00000031598 
# PC_ 4 
# Positive:  necab2:ENSDARG00000056745, calca:ENSDARG00000056590, rnd3b:ENSDARG00000007396, calb1:ENSDARG00000031598, bmp1b:ENSDARG00000028053, pcdh10a:ENSDARG00000099729, frzb:ENSDARG00000018383, calb2b:ENSDARG00000036344, chrna2b:ENSDARG00000057025, mecom:ENSDARG00000060808 
# nefmb:ENSDARG00000043697, sox2:ENSDARG00000070913, tnfaip8l3:ENSDARG00000088709, kcna1a:ENSDARG00000062942, ncs1a:ENSDARG00000055229, irx1a:ENSDARG00000101831, megf11:ENSDARG00000062686, NA:ENSDARG00000117377, atp1a3b:ENSDARG00000104139, zgc:194990:ENSDARG00000069407 
# cbln2b:ENSDARG00000077151, isl2a:ENSDARG00000003971, chga:ENSDARG00000008829, mnx1:ENSDARG00000035984, lhx4:ENSDARG00000039458, chka:ENSDARG00000041078, chgb:ENSDARG00000076500, gadd45gb.1:ENSDARG00000016725, irx3a:ENSDARG00000101076, cdh6:ENSDARG00000014522 
# Negative:  pdyn:ENSDARG00000087798, itga3a:ENSDARG00000037917, pnoca:ENSDARG00000025024, hunk:ENSDARG00000091317, bmp16:ENSDARG00000103679, dock1:ENSDARG00000099093, glra4a:ENSDARG00000006865, her4.5:ENSDARG00000056729, her15.1:ENSDARG00000054562, atp1b1a:ENSDARG00000013144 
# rasd4:ENSDARG00000099238, cox4i1l:ENSDARG00000012388, her15.2:ENSDARG00000054560, her4.2:ENSDARG00000094426, npas1:ENSDARG00000015876, olig4:ENSDARG00000052610, pcdh19:ENSDARG00000034344, her4.1:ENSDARG00000056732, notch1a:ENSDARG00000103554, arrdc2:ENSDARG00000020761 
# sox11a:ENSDARG00000077811, ret:ENSDARG00000055305, neurod4:ENSDARG00000003469, stmn1b:ENSDARG00000033655, ebf2:ENSDARG00000042525, hsp90aa1.2:ENSDARG00000024746, atp2a3:ENSDARG00000060978, cox4i2:ENSDARG00000022509, kalrna:ENSDARG00000104119, si:dkeyp-110a12.4:ENSDARG00000087784 
# PC_ 5 
# Positive:  stmn1b:ENSDARG00000033655, hoxb13a:ENSDARG00000056015, zfhx4:ENSDARG00000075542, stmn4l:ENSDARG00000043932, hsp70l:ENSDARG00000055723, nr5a2:ENSDARG00000100940, ppp1r14ba:ENSDARG00000044541, zfhx3:ENSDARG00000103057, arhgdig:ENSDARG00000004034, vasnb:ENSDARG00000102565 
# mecom:ENSDARG00000060808, sox4a:ENSDARG00000004588, gap43:ENSDARG00000099744, fgfr3:ENSDARG00000004782, zbtb16a:ENSDARG00000007184, fsta:ENSDARG00000052846, tubb5:ENSDARG00000037997, hoxa13a:ENSDARG00000100312, kif26ab:ENSDARG00000015016, cd99l2:ENSDARG00000056722 
# dclk1b:ENSDARG00000104664, spegb:ENSDARG00000009567, plppr3a:ENSDARG00000010144, hunk:ENSDARG00000091317, nrgna:ENSDARG00000039626, sox11a:ENSDARG00000077811, prph:ENSDARG00000028306, hapln1a:ENSDARG00000089769, NA:ENSDARG00000092467, pappaa:ENSDARG00000101090 
# Negative:  her15.1:ENSDARG00000054562, her4.2:ENSDARG00000094426, her15.2:ENSDARG00000054560, her4.5:ENSDARG00000056729, her4.1:ENSDARG00000056732, notch1a:ENSDARG00000103554, ndnfl:ENSDARG00000076462, scn4bb:ENSDARG00000060319, nefmb:ENSDARG00000043697, tmprss4b:ENSDARG00000012860 
# her9:ENSDARG00000056438, tppp3:ENSDARG00000030463, atp1a3b:ENSDARG00000104139, eef1a2:ENSDARG00000006838, islr2:ENSDARG00000051875, megf11:ENSDARG00000062686, sema3aa:ENSDARG00000019235, abhd6a:ENSDARG00000060756, chgb:ENSDARG00000076500, zgc:65851:ENSDARG00000012281 
# her13:ENSDARG00000007097, angptl3:ENSDARG00000044365, rxfp3.2b:ENSDARG00000061846, her6:ENSDARG00000006514, thbs2b:ENSDARG00000073810, calb2a:ENSDARG00000041062, pvalb6:ENSDARG00000009311, sox19a:ENSDARG00000010770, zbtb4:ENSDARG00000105255, sncga:ENSDARG00000034423 

#check differences between datasets
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = MN_V2a, reduction = "pca", pt.size = 1, group.by = "celltype")
p2 <- VlnPlot(object = MN_V2a, features = "PC_1", group.by = "celltype", pt.size = .1)
plot_grid(p1,p2)
p1 <- DimPlot(object = MN_V2a, reduction = "pca", pt.size = 1, group.by = "orig.cluster")
p2 <- VlnPlot(object = MN_V2a, features = "PC_1", group.by = "orig.cluster", pt.size = .1)
plot_grid(p1,p2)


#run harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
Harmony_MN_V2a <- MN_V2a %>% 
  RunHarmony("celltype", plot_convergence = TRUE)
# Harmony converged after 2 iterations
# Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details on syntax validity
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = Harmony_MN_V2a, reduction = "harmony", pt.size = 1, group.by = "celltype")
p2 <- VlnPlot(object = Harmony_MN_V2a, features = "harmony_1", group.by = "celltype", pt.size = .1)
plot_grid(p1,p2)


#Clustering####
Harmony_MN_V2a <- Harmony_MN_V2a %>% 
  RunUMAP(reduction = "harmony", dims = pcs) %>% 
  RunTSNE(reduction = "harmony", dims = pcs) %>%
  FindNeighbors(reduction = "harmony", dims = pcs) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

ElbowPlot(Harmony_MN_V2a)
Harmony_MN_V2a <- JackStraw(Harmony_MN_V2a)
Harmony_MN_V2a <- ScoreJackStraw(Harmony_MN_V2a, dims = 1:20)
pcs <- which(JS(object = Harmony_MN_V2a[['pca']], slot = 'overall')[, 2] < 1e-3)
pcs
#[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 17 18

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "celltype", pt.size = 2, split.by = 'celltype')
DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "orig.cluster", pt.size = 2, split.by = 'celltype')
DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "celltype", pt.size = 3)
DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "orig.cluster", pt.size = 3)

DimPlot(Harmony_MN_V2a, reduction = "tsne", group.by = "celltype", pt.size = 3)
DimPlot(Harmony_MN_V2a, reduction = "tsne", group.by = "orig.cluster", pt.size = 3)

options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(Harmony_MN_V2a, reduction = "umap", label = TRUE, pt.size = 3)

png("harmony_MN2plates_V2a_UMAP_clusters_celltype.png", width=920, height=500)
p1 <- DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "celltype", pt.size = 3)
p2 <- DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "orig.cluster", pt.size = 3)
plot_grid(p1,p2)
dev.off()

#Figure 6a
setwd("Integration MN_V2a/Figure")
pdf("harmony_MN2plates_V2a_UMAP_clusters_celltype.pdf", width=12, height=6)
p1 <- DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "celltype", pt.size = 2)
p2 <- DimPlot(Harmony_MN_V2a, reduction = "umap", group.by = "orig.cluster", pt.size = 2)
plot_grid(p1,p2)
dev.off()
setwd("Integration MN_V2a")

saveRDS(Harmony_MN_V2a, "C:/Users/Irene/Desktop/New analysis/Integration MN + V2a/Harmony/Integrated_MN_V2a.rds")

#UMAP of merged dataset without harmony for comparison####
MN_V2a<-RunUMAP(MN_V2a, reduction = "pca", dims = 1:20)

png("UMAP_ MN2_V2a_merged_nonharmonized_clusters_celltype.png", width=920, height=500)
p1 <- DimPlot(MN_V2a, reduction = "umap", group.by = "celltype", pt.size = 3)
p2 <- DimPlot(MN_V2a, reduction = "umap", group.by = "orig.cluster", pt.size = 3)
plot_grid(p1,p2)
dev.off()


