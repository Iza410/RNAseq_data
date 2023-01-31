setwd('~/Documents/UniBern/6_2022/5_RNA-seq/results')

library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)


# read file
data<-read.table('featureCounts.txt', header=TRUE, row.names = 1)

#file preparation
##remove first 5 columns
data<-data[, -c(1:5)]
##rename columns- leave only name of sample
colnames(data)<- sub("X.data.users.ibiedron.rnaseq.results.sorted.", "", colnames(data))
colnames(data)<- sub("_sort.bam", "", colnames(data))
## remove the last character (= replicate number), create sample names
sampleGroup<-sub(".$", "", colnames(data)) 

########################################
########################################annotation file
########################################
#library(AnnotationDbi)
#library(org.Hs.eg.db)
#library(EnsDb.Hsapiens.v86)
#library(dplyr)
#library(biomaRt)

#returns a list of BioMart databases to which biomaRt can connect to
#biolist <- as.data.frame(listMarts())

#select a BioMart database and dataset to use
#ensembl = useMart("ensembl")

#which datasets are available in the selected BioMart by using the function listDatasets()
#esemblist <- as.data.frame(listDatasets(ensembl))

#To select a dataset we can update the Mart object using the function useDataset()
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

#The listFilters() function shows you all available filters in the selected dataset.
#filters = listFilters(ensembl)

#Attributes define the values we are interested in to retrieve. For example we want 
#to retrieve the gene symbols or chromosomal coordinates. The listAttributes() 
#function displays all available attributes in the selected dataset.
#attributes = listAttributes(ensembl)

#The getBM() function is the main query function in biomaRt. It has four main arguments:
#attributes: is a vector of attributes that one wants to retrieve (= the output of the query).
#filters: is a vector of filters that one wil use as input to the query.
#values: a vector of values for the filters. In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument (see examples below).
#mart: is an object of class Mart, which is created by the useMart() function.

#t2g<-getBM(attributes=c('ensembl_gene_id', 'ensembl_peptide_id', 'chromosome_name','start_position','end_position'), mart = ensembl)

# Note, some IDs get >1 mapping, others remain empty.
#select <- AnnotationDbi::select
#ann <- select(org.Hs.eg.db, keys = as.character(row.names(data)), keytype = 'ENSEMBL', columns = c("ENTREZID", "SYMBOL", "GENENAME", "ALIAS", "MAP"))

# Add chrom info:
#ann <- merge(x=ann, y=t2g, by.x="ENSEMBL", by.y="ensembl_gene_id", all.x=TRUE)

# Make a dataframe where gene from S190910.mat is listed only once. If there are multiple "SYMBOL", take the first. If "SYMBOL" is NA, replace with ENS id:
#data.ann <- ann[!duplicated(ann["ENSEMBL"]),]
#data.ann <- data.ann[c("ENSEMBL", "SYMBOL", "MAP", "chromosome_name")]
#data.ann$SYMBOL <- ifelse(is.na(data.ann$SYMBOL), data.ann$ENSEMBL, data.ann$SYMBOL)

# set ENSEMBL as rownames
#rownames(data.ann) <- data.ann[,1]
#data.ann <- data.ann[-1]
#write.csv(data.ann, "anno.csv", row.names=TRUE)
anno <- read.csv("anno_f.csv", row.names = 1, sep=";")

#library(tidyverse)
#sum(duplicated(anno$SYMBOL)) #107 duplicated gene names
########################################
########################################
########################################
#gene names
anno$ens <- row.names(anno)
data$ens <-row.names(data)
df_sym <- merge(data, anno, by = 'ens')
rownames(df_sym) <- df_sym$SYMBOL
df <- df_sym[ , !names(df_sym) %in% c("ens","MAP", "chromosome_name", "SYMBOL")]
####################################



colData<-data.frame(condition=factor(sampleGroup))

#create deseq object
ds_obj<-DESeqDataSetFromMatrix(df, colData, formula(~condition))

#run DESeq2
dds<-DESeq(ds_obj, betaPrior = TRUE)
res <-results(dds)

#hist(res$pvalue, breaks=50, col="grey")

# Plot dispersion
#png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
#dev.off()

#visualize differentially expressed genes
plotMA(dds)


# apply a regularized log transformation, ignoring information about experimental groups
rld <-rlog(dds, blind=TRUE)
#hist(assay(rld))

#PCA using the most variably expressed genes
png("pca.png", width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 1200)
plotPCA(rld, intgroup=c("condition")) + 
  theme_linedraw() +
  scale_colour_brewer(palette = "Set2") 
dev.off()

#Hierarchical Clustering Heatmap
library("pheatmap")
expmatrix <- SummarizedExperiment::assay(dds)
sampleDists <- dist(t(expmatrix))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$condition)
colnames(sampleDistMatrix) <- NULL
colors <- grDevices::colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)

png("hchmap.png", width     = 10,
    height    = 10,
    units     = "in",
    res       = 1200)
pheatmap::pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,col=colors, fontsize = 15)
dev.off()

############################heatmap
#interesting genes, heatmap based on 50 most highly expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
hmcol <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
heatmap(assay(rld)[select,], col = hmcol, trace="none", margin=c(10,6),
        labCol=colnames(dds), cexRow = 1/log10(length(select)), Colv = NA)






#######################################################
#create all pairs to compare expression
variants <- unique(sampleGroup)
pairs <- combn(variants, 2, FUN = NULL, simplify = TRUE)

for (n in 1:ncol(pairs)) {
  assign(paste0('res_', pairs[1,n], '_', pairs[2,n]),
         results(dds, contrast=  c("condition", pairs[1,n], pairs[2,n]) ,alpha = 0.05))
  
}

##################################HER2 vs Normal_1
#Extract the DE test results for at least one pairwise contrast
summary(res_HER2_Normal) #up 9477, down 5000
sum(res_HER2_Normal$padj < 0.05, na.rm=TRUE) #14477

#Mean expression against log-fold change
ME_HER2_Normal <- plotMA(res_HER2_Normal, alpha=0.05)
show(ME_HER2_Normal)
##########################volcano plot
#convert res to a dataframe and remove all genes that don't have an adjusted p-value
res_volcano_HNOR <- as.data.frame(res_HER2_Normal)
res_volcano_HNOR <- res_volcano_HNOR[!is.na(res_volcano_HNOR$padj),]

res_volcano_HNOR <- res_volcano_HNOR %>%
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "up",
                               log2FoldChange <= -2 & padj < 0.05 ~ "down",
                               TRUE~ "ns" ))
res_volcano_HNOR %>%
  count(Expression)

res_volcano_HNOR$gene=row.names(res_volcano_HNOR)

sig_il_genes <- res_volcano_HNOR %>%
  filter(gene %in% c("RN7SL3", "PPP1R15A", "RN7SL1",
                     "ENSG00000263934", "SCARNA7", "ENSG00000265735",
                     "RMRP", "RN7SL4P", "RNU5A-1",
                     "MAFF", "RNU4-2", "SNORA73B", "ATF3"))



cols <- c("up" = "firebrick3", "down" = "dodgerblue3", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 2) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
#ggplot
volcano_HER2_Normal <- ggplot(data=res_volcano_HNOR, aes(x=log2FoldChange, 
                                                         y=-log10(padj), 
                                                         fill = Expression,    
                                                         size = Expression,
                                                         alpha = Expression)) + 

  geom_point(aes(colour = 'black'), 
             alpha = 0.5, 
             shape = 21,
             size = 2) +


  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  geom_label_repel(data = sig_il_genes, # Add labels last to appear as the top layer  
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1,
                   show_guide=FALSE,
                   alpha = 0.75) +
  scale_colour_manual(values = cols) +
  xlim(-15,15) +
  labs(title = "Gene expression changes in HER2 versus Normal",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_linedraw() 

volcano_HER2_Normal

###EnhancedVolcano
volcano_en_HER2_Normal <- EnhancedVolcano(res_volcano_HNOR,
                lab = rownames(res_volcano_HNOR),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'HER2 versus Normal',
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 2,
                labSize = 3.0,
                colAlpha = 1,
                col=c("black","black", "black", "red3"),
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5)
volcano_en_HER2_Normal
###################################################################
#################################HER2 vs NonTNBC_2
#Extract the DE test results for at least one pairwise contrast
summary(res_HER2_NonTNBC) #uo 2615, down 2633
sum(res_HER2_NonTNBC$padj < 0.05, na.rm=TRUE) #5248

#Mean expression against log-fold change
ME_HER2_NonTNBC <- plotMA(res_HER2_NonTNBC, alpha=0.05)
show(ME_HER2_NonTNBC)
##########################volcano plot
#convert res to a dataframe and remove all genes that don't have an adjusted p-value
res_volcano_HNon <- as.data.frame(res_HER2_NonTNBC)
res_volcano_HNon <- res_volcano_HNon[!is.na(res_volcano_HNon$padj),]
#ggplot
res_volcano_HNon <- res_volcano_HNon %>%
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "up",
                                log2FoldChange <= -2 & padj < 0.05 ~ "down",
                                TRUE~ "ns" ))
res_volcano_HNon %>%
  count(Expression)

res_volcano_HNon$gene=row.names(res_volcano_HNon)

sig_il_genes_HNon <- res_volcano_HNon %>%
  filter(gene %in% c("AGR3", "ESR1", "CYP2B7P",
                     "SNORA53", "NKAIN1", "TMEM178B",
                     "GRB7", "SEL1L3", "ST6GALNAC2",
                     "RIMS4", "KCNK15", "DSCAM-AS1", 
                     "GFRA1", "HIF1A-AS3", "SKP2"))

#ggplot
volcano_HNon <- ggplot(data=res_volcano_HNon, aes(x=log2FoldChange, 
                                                         y=-log10(padj), 
                                                         fill = Expression,    
                                                         size = Expression,
                                                         alpha = Expression)) + 
  geom_point(aes(colour = 'black'), 
             alpha = 0.5, 
             shape = 21,
             size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  geom_label_repel(data = sig_il_genes_HNon, # Add labels last to appear as the top layer  
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1,
                   show_guide=FALSE,
                   alpha = 0.75) +
  scale_colour_manual(values = cols) +
  xlim(-15,15) +
  labs(title = "Gene expression changes in \nHER2 versus Non-TNBC",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_linedraw() 
png("her_non.png", width     = 4,
    height    = 4,
    units     = "in",
    res       = 1200)
volcano_HNon
dev.off()


###EnhancedVolcano
volcano_en_HER2_NonTNBC <- EnhancedVolcano(res_volcano_HNon,
                                          lab = rownames(res_volcano_HNon),
                                          x = 'log2FoldChange',
                                          y = 'pvalue')
volcano_en_HER2_NonTNBC

###################################################################
#################################HER2 vs TNBC_3
#Extract the DE test results for at least one pairwise contrast
summary(res_HER2_TNBC) #uo 2051, down 3636
sum(res_HER2_TNBC$padj < 0.05, na.rm=TRUE) #5687

#Mean expression against log-fold change
ME_HER2_TNBC <- plotMA(res_HER2_TNBC, alpha=0.05)
show(ME_HER2_TNBC)
##########################volcano plot
#convert res to a dataframe and remove all genes that don't have an adjusted p-value
res_volcano_HTn <- as.data.frame(res_HER2_TNBC)
res_volcano_HTn <- res_volcano_HTn[!is.na(res_volcano_HTn$padj),]
#ggplot
res_volcano_HTn <- res_volcano_HTn %>%
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "up",
                                log2FoldChange <= -2 & padj < 0.05 ~ "down",
                                TRUE~ "ns" ))
res_volcano_HTn %>%
  count(Expression)

res_volcano_HTn$gene=row.names(res_volcano_HTn)

sig_il_genes_HTNBC <- res_volcano_HTn %>%
  filter(gene %in% c("SNORA53", "GRB7", "SNORA22",
                     "CDKN2A", "ERBB2", "RNU5A-1",
                     "SNORA65", "SNORA80E", "VTRNA1-1",
                     "MT-TM", "SNORA37", "SNORA22C", "SNORA74B",
                     "FBN3","PPP1R14A", "LOC124902611"))

#ggplot
volcano_HTNBC <- ggplot(data=res_volcano_HTn, aes(x=log2FoldChange, 
                                                  y=-log10(padj), 
                                                  fill = Expression,    
                                                  size = Expression,
                                                  alpha = Expression)) + 
  geom_point(aes(colour = 'black'), 
             alpha = 0.5, 
             shape = 21,
             size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2,2), linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  geom_label_repel(data = sig_il_genes_HTNBC, # Add labels last to appear as the top layer  
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1,
                   show_guide=FALSE,
                   alpha = 0.75) +
  scale_colour_manual(values = cols) +
  xlim(-15,15) +
  labs(title = "Gene expression changes in \nHER2 versus TNBC",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_linedraw() 

png("her_tnbc.png", width     = 4,
    height    = 4,
    units     = "in",
    res       = 1200)
volcano_HTNBC
dev.off()
###EnhancedVolcano
volcano_en_HER2_TNBC <- EnhancedVolcano(res_volcano_HTn,
                                           lab = rownames(res_volcano_HTn),
                                           x = 'log2FoldChange',
                                           y = 'pvalue')
volcano_en_HER2_TNBC

###################################################################
#################################NonTNBC vs Normal_4
#Extract the DE test results for at least one pairwise contrast
summary(res_NonTNBC_Normal) #up 8723, down 4064
sum(res_NonTNBC_Normal$padj < 0.05, na.rm=TRUE) #12787

#Mean expression against log-fold change
ME_NonTNBC_Normal <- plotMA(res_NonTNBC_Normal, alpha=0.05)
show(ME_NonTNBC_Normal)
##########################volcano plot
#convert res to a dataframe and remove all genes that don't have an adjusted p-value
res_volcano_NonN <- as.data.frame(res_NonTNBC_Normal)
res_volcano_NonN <- res_volcano_NonN[!is.na(res_volcano_NonN$padj),]
#ggplot
res_volcano_NonN <- res_volcano_NonN %>%
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "up",
                                log2FoldChange <= -2 & padj < 0.05 ~ "down",
                                TRUE~ "ns" ))
res_volcano_NonN %>%
  count(Expression)

res_volcano_NonN$gene=row.names(res_volcano_NonN)

sig_il_genes_NonN <- res_volcano_NonN %>%
  filter(gene %in% c("RN7SL3", "RN7SL4P", "SCARNA7",
                     "RN7SL5P", "RN7SL1", "SNORD3A",
                     "PPP1R15A", "RMRP", "MAFF",
                     "UBC", "SNORA73B", "ATF3", "RASD1"))

#ggplot
volcano_NonTNBC_Normal <- ggplot(data=res_volcano_NonN, aes(x=log2FoldChange, 
                                                  y=-log10(padj), 
                                                  fill = Expression,    
                                                  size = Expression,
                                                  alpha = Expression)) + 
  geom_point(aes(colour = 'black'), 
             alpha = 0.5, 
             shape = 21,
             size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2,2), linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  geom_label_repel(data = sig_il_genes_NonN, # Add labels last to appear as the top layer  
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1,
                   show_guide=FALSE,
                   alpha = 0.75) +
  scale_colour_manual(values = cols) +
  xlim(-15,15) +
  labs(title = "Gene expression changes in NonTNBC versus Normal",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_linedraw() 

volcano_NonTNBC_Normal



###EnhancedVolcano
volcano_en_NonTNBC_Normal <- EnhancedVolcano(res_volcano_NonN,
                                        lab = rownames(res_volcano_NonN),
                                        x = 'log2FoldChange',
                                        y = 'pvalue')
volcano_en_NonTNBC_Normal

###################################################################
#################################NonTNBC vs TNBC_5
#Extract the DE test results for at least one pairwise contrast
summary(res_NonTNBC_TNBC) #up 934, down 812
sum(res_NonTNBC_TNBC$padj < 0.05, na.rm=TRUE) #1746

#Mean expression against log-fold change
ME_NonTNBC_TNBC <- plotMA(res_NonTNBC_TNBC, alpha=0.05)
show(ME_NonTNBC_TNBC)
##########################volcano plot
#convert res to a dataframe and remove all genes that don't have an adjusted p-value
res_volcano_NonT <- as.data.frame(res_NonTNBC_TNBC)
res_volcano_NonT <- res_volcano_NonT[!is.na(res_volcano_NonT$padj),]
#ggplot
res_volcano_NonT <- res_volcano_NonT %>%
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "up",
                                log2FoldChange <= -2 & padj < 0.05 ~ "down",
                                TRUE~ "ns" ))
res_volcano_NonT %>%
  count(Expression)

res_volcano_NonT$gene=row.names(res_volcano_NonT)

sig_il_genes_NonT <- res_volcano_NonT %>%
  filter(gene %in% c("ESR1", "AGR3", "ELF5",
                     "AGR2", "TFF1", "ROPN1B",
                     "FSIP1", "ACADSB", "FOXA1",
                     "ACTG2", "KCNK15", "CYP2B7P", "PRRT3"))

#ggplot
volcano_NonT <- ggplot(data=res_volcano_NonT, aes(x=log2FoldChange, 
                                                            y=-log10(padj), 
                                                            fill = Expression,    
                                                            size = Expression,
                                                            alpha = Expression)) + 
  geom_point(aes(colour = 'black'), 
             alpha = 0.5, 
             shape = 21,
             size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  geom_label_repel(data = sig_il_genes_NonT, # Add labels last to appear as the top layer  
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1,
                   show_guide=FALSE,
                   alpha = 0.75) +
  scale_colour_manual(values = cols) +
  xlim(-15,15) +
  labs(title = "Gene expression changes in NonTNBC versus TNBC",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_linedraw() 

volcano_NonT

###EnhancedVolcano
volcano_en_NonTNBC_TNBC <- EnhancedVolcano(res_volcano_NonT,
                                             lab = rownames(res_volcano_NonT),
                                             x = 'log2FoldChange',
                                             y = 'pvalue')
volcano_en_NonTNBC_TNBC

###################################################################
#################################Normal vs TNBC_6
#Extract the DE test results for at least one pairwise contrast
summary(res_Normal_TNBC) #up 4281, down 11944
sum(res_Normal_TNBC$padj < 0.05, na.rm=TRUE) #16225

#Mean expression against log-fold change
ME_Normal_TNBC <- plotMA(res_Normal_TNBC, alpha=0.05)
show(ME_Normal_TNBC)
##########################volcano plot
#convert res to a dataframe and remove all genes that don't have an adjusted p-value
res_volcano_NT <- as.data.frame(res_Normal_TNBC)
res_volcano_NT <- res_volcano_NT[!is.na(res_volcano_NT$padj),]
#ggplot
res_volcano_NT <- res_volcano_NT %>%
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "up",
                                log2FoldChange <= -2 & padj < 0.05 ~ "down",
                                TRUE~ "ns" ))
res_volcano_NT %>%
  count(Expression)

res_volcano_NT$gene=row.names(res_volcano_NT)

sig_il_genes_NT <- res_volcano_NT %>%
  filter(gene %in% c("RN7SL5P", "RN7SL3", "SCARNA7",
                     "RN7SL4P", "PPP1R15A", "RN7SL1",
                     "UBC", "RMRP", "ATF3",
                     "IL6", "ENSG00000281181", "ENSG00000280800", "MAFF"))

#ggplot
volcano_NT <- ggplot(data=res_volcano_NT, aes(x=log2FoldChange, 
                                                  y=-log10(padj), 
                                                  fill = Expression,    
                                                  size = Expression,
                                                  alpha = Expression)) + 
  geom_point(aes(colour = 'black'), 
             alpha = 0.5, 
             shape = 21,
             size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2,2), linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  geom_label_repel(data = sig_il_genes_NT, # Add labels last to appear as the top layer  
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1,
                   show_guide=FALSE,
                   alpha = 0.75) +
  scale_colour_manual(values = cols) +
  xlim(-15,15) +
  labs(title = "Gene expression changes in NonTNBC versus TNBC",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)") +
  theme_linedraw() 

volcano_NT

###EnhancedVolcano
volcano_en_Normal_TNBC <- EnhancedVolcano(res_volcano_NT,
                                           lab = rownames(res_volcano_NT),
                                           x = 'log2FoldChange',
                                           y = 'pvalue')
volcano_en_Normal_TNBC




library("gplots")
#interesting genes, heatmap based on 35 most highly expressed genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:35]
hmcol <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
bcol <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(100)
heatmap(assay(rld)[select,], col = hmcol, trace="none", margin=c(10,6),
        labCol=colnames(dds), cexRow = 1/log10(length(select)), Colv = NA)

png("heat_exp1.png", width     = 15,
    height    = 13,
    units     = "in",
    res       = 1200)
heatmap.2(assay(rld)[select,], col = bcol, trace="none", margin=c(10,8), 
           Rowv=TRUE, labCol=colnames(dds),  
          Colv = TRUE,keysize=0.75, key.par = list(cex=0.75))
dev.off()

####genes
library(ggpubr)
cdk <- plotCounts(dds, gene="CDKN2A", intgroup="condition", returnData = TRUE )
p1 <- ggplot(cdk, aes(x=condition, y=count, color=condition)) + 
  theme_linedraw() + geom_point(size=3) +
  scale_colour_brewer(palette = "Set2") 

snora <- plotCounts(dds, gene="SNORA53", intgroup="condition", returnData = TRUE )
p2 <- ggplot(snora, aes(x=condition, y=count, color=condition)) + 
  theme_linedraw() + geom_point(size=3) +
  scale_colour_brewer(palette = "Set2") 

esr <- plotCounts(dds, gene="ESR1", intgroup="condition", returnData = TRUE )
p3 <- ggplot(esr, aes(x=condition, y=count, color=condition)) + 
  theme_linedraw() + geom_point(size=3) +
  scale_colour_brewer(palette = "Set2") 

png("genes.png", width     = 9,
    height    = 4,
    units     = "in",
    res       = 1200)
ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()
