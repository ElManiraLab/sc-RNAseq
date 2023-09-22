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

#Read data from GSE232801 and reproduce figures####
#Citation "Kelly JJ, Wen H, Brehm P. Single cell RNA-seq analysis of spinal locomotor circuitry in larval zebrafish. bioRxiv 2023 Sep 11. PMID: 37333232"

spine <- readRDS("GSE232801_Integrated_spine.RDS")
MNall <- readRDS("GSE232801_Motor_Neuron.RDS")
neurons <- readRDS("GSE232801_Neurons.RDS")

my_cols <- c('1'='#9bddb1','2'='#8cdaec','3'='#3cb464','4'='#d48c84','5'='#a89a49',
             '6'='#d6cfa2','7'='#b45248','8'='#167288','9'='#643c6a','10'='#836394')

setwd("Larva_analysis")
pdf("all MN umap.pdf")
DimPlot(MNall,reduction="umap", cols=my_cols)
dev.off()

pdf("all MN tsne.pdf")
DimPlot(MNall,reduction="tsne", cols=my_cols)
dev.off()

pdf("Reproduce Fig.5c_allMN integrated.pdf")
DimPlot(MNall,reduction="tsne",label=FALSE)+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position = "none")
dev.off()

pdf("Reproduce Fig.5c_allMN integrated_2.pdf")
FeaturePlot(MNall,reduction = "tsne",label = FALSE,features = "nr2f1a")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position = "none")
dev.off()

pdf("Reproduce Fig.5c_allMN integrated_3.pdf")
FeaturePlot(MNall,reduction = "tsne",label = FALSE,features = "chga")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position = "none")
dev.off()

my_levels <- c('7', '8', '5', '10', '6', '2', '3', '1', '4', '9')
MNall@active.ident <- factor(x = MNall@active.ident, levels = my_levels)
my_cols <- c('1'='#9bddb1','2'='#8cdaec','3'='#3cb464','4'='#d48c84','5'='#a89a49',
             '6'='#d6cfa2','7'='#b45248','8'='#167288','9'='#643c6a','10'='#836394')

#MN_Plot markers and adult genes####
FeaturePlot(neurons, features = c("mnx1"),  reduction = "umap", pt=0.3, cols = c("#E0E0E0", "#0066A6"))
FeaturePlot(neurons, features = c("foxp1b"),  reduction = "umap", pt=0.3, cols = c("#E0E0E0", "#0066A6"))

pdf("mnx1_umap.pdf", width = 10, height = 7)
FeaturePlot(neurons, features = c("mnx1"),  reduction = "umap", pt=0.3, cols = c("#E0E0E0", "#0066A6"))
dev.off()

pdf("mnx1_tsne.pdf", width = 10, height = 7)
FeaturePlot(neurons, features = c("mnx1"),  reduction = "tsne", pt=0.3, cols = c("#E0E0E0", "#0066A6"))
dev.off()

MNs<-subset(neurons,idents = c(26,4,1,2,7,20))

#markers per cluster adult 
DotPlot(MNs, features = c("itga3a", "bmp16", "ret","pcdh11", "oxr1a", "grik2", "grin1b", "esrrga", "tmprss4b", "scn4bb", "pvalb6", "atp1a3b", "apoa", "glrba", "stmn3", "chrna2b", "calb1", "fndc4a", "pdch10a", "sox2", "necab2", "ramp1", "hoxb13", "zfhx4", "hoxc13b", "esrrb", "nhlh2", "olig4", "ebf2", "tuba8l5", "elavl3", "rassf1", "nkx6.1"),  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)

#TF per cluster adult
DotPlot(MNs, features = c("rest", "gtf2i", "npas3", "npas1", "esrrga","bhlhe41", "trim13", "pcbp3", "nr2f2", "nupr1a", "bhlhe40", "ebf1b", "hoxb3", "tbpl2", "mecom", "neurod1", "pknox2", "nfixb", "isl2a", "mnx1", "lhx4", "ebf3a", "irx5a", "zfhx4", "hoxb13a", "hoxc13b", "mnx2b", "hoxc13a", "hoxa13a", "esrrb", "zfhx3", "sp8b"),  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
DotPlot(MNs, features = c("itga3a", "bmp16", "ret","pcdh11", "oxr1a", "grik2", "grin1b", "esrrga", "tmprss4b", "scn4bb", "pvalb6", "atp1a3b", "apoa", "glrba", "stmn3", "chrna2b", "calb1", "fndc4a", "pdch10a", "sox2", "necab2", "ramp1", "hoxb13", "zfhx4", "hoxc13b", "esrrb", "nhlh2", "olig4", "ebf2", "tuba8l5", "elavl3", "rassf1", "nkx6.1"),  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)

my_levels_2 <- c('26','4','1','2','7','20')
MNs@active.ident <- factor(x = MNs@active.ident, levels = my_levels_2)
my_cols_MN <- c('1'='#9bddb1','2'='#d6cfa2','7'='#d48c84','4'='#167288','20'='#836394',
            '26'='#b45248')

pdf("larva_markers_featureplot.pdf", width = 10, height =5)
FeaturePlot(MNs, features = c("calca", "chrna2b", "pcdh9", "esrrga", "foxp1b"),  reduction = "umap", pt=0.3, cols = c("#E0E0E0", "#0066A6"))
dev.off()

#Ext Data Fig.8a-d
setwd("Larva_analysis/Figures")

pdf("larva_markers_violin plot.pdf", width = 10, height = 7)
StackedVlnPlot(MNs, features = c( "pcdh9","calca", "foxp1b"), cols=my_cols_MN)
dev.off()

pdf("adultmarkers_violin plot.pdf", width = 10, height = 30)
StackedVlnPlot(MNs, features = c("esrrga","grin1b", "pvalb6", "neurod1","chrna2b", "ebf3a"), cols=my_cols_MN)
dev.off()

pdf("MNs.umap.pdf")
DimPlot(MNs, reduction = "umap", cols=my_cols_MN)
dev.off()
setwd("Larva_analysis")

#V2a_plots####

pdf("all neurons umap.pdf")
plot<-DimPlot(neurons,reduction="umap")+ NoLegend()
LabelClusters(plot = plot, id = "ident")
dev.off()

pdf("vsx2_umap.pdf", width = 10, height = 7)
FeaturePlot(neurons, features = c("vsx2"),  reduction = "umap", pt=0.3, cols = c("#E0E0E0", "#0066A6"))
dev.off()

v2a<-subset(neurons,idents = c(19,29,23))#v2a

#TF per cluster adult
pdf('V2aDotplot.adult.markers.pdf',  width=7, height=6)
DotPlot(v2a, features = c("dacha", "oxr1b", "isl2a", "esrrga", "pbx3b", "pcp4a", "tmem88b", "tns1b", "angptl3", "zfhx3", "zfhx4", "nr5a2", "olfm3a", "ntng1a", "syt13", "pcdh17", "spon1b", "chata", "slc18a3a", "fgfr2", "aclya", "igsf21b", "sema3ab", "prkg1b"),  cols = c("white", "#0066A6")) + coord_flip() +scale_x_discrete(limits=rev)
dev.off()

pdf("V2a_TF.pdf", width = 10, height = 20)
StackedVlnPlot(v2a, features = c("dacha", "pbx3b", "meis1b", "esrrga", "sp8a", "zfhx3", "zfhx4", "nr5a2", "meis2a", "shox2", "nfixb", "neurod1", "ebf1a", "neurod6a"))
dev.off()

#Ext Data Fig.8e-h
setwd("Larva_analysis/Figures")
my_cols <- c('1'='#9bddb1','2'='#8cdaec','3'='#3cb464','4'='#d48c84','5'='#a89a49',
             '6'='#d6cfa2','7'='#b45248','8'='#167288','9'='#643c6a','10'='#836394')

my_levels_v2a <- c('23', '19', '29')
v2a@active.ident <- factor(x = v2a@active.ident, levels = my_levels_v2a)

my_cols_v2a <- c('23'='#b45248','29'='#167288','19'='#a89a49')

pdf("V2a.umap.pdf")
DimPlot(v2a, reduction = "umap", cols=my_cols_v2a)
dev.off()

pdf("V2a_all markers.pdf", width = 10, height = 15)
StackedVlnPlot(v2a, features = c("esrrga", "sp8a", "shox2", "zfhx3", "slc18a3a", "ebf1a"), cols=my_cols_v2a)
dev.off()
