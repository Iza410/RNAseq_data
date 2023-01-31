setwd('~/Documents/UniBern/6_2022/5_RNA-seq/results')
library(ggpubr)
library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(RColorBrewer)
library(ggplot2)


data<-read.table('featureCounts.txt', header=TRUE, row.names = 1)
data<-data[, -c(1:5)]
colnames(data)<- sub("X.data.users.ibiedron.rnaseq.results.sorted.", "", colnames(data))
colnames(data)<- sub("_sort.bam", "", colnames(data))
sampleGroup<-sub(".$", "", colnames(data)) 
colData<-data.frame(condition=factor(sampleGroup))
#create deseq object
ds_obj<-DESeqDataSetFromMatrix(data, colData, formula(~condition))
#run DESeq2
dds<-DESeq(ds_obj, betaPrior = TRUE)
res <-results(dds)
variants <- unique(sampleGroup)
pairs <- combn(variants, 2, FUN = NULL, simplify = TRUE)

for (n in 1:ncol(pairs)) {
  assign(paste0('res_', pairs[1,n], '_', pairs[2,n]),
         results(dds, contrast=  c("condition", pairs[1,n], pairs[2,n]) ,alpha = 0.05))
  
}

# SET THE DESIRED ORGANISM 
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#her2- nontnbc
df_HN <- as.data.frame(res_HER2_NonTNBC)

#her2-tnbc
df_HT <- as.data.frame(res_HER2_TNBC)
##################################################GO
#H-TNBC
all_genes_HT <- as.character(rownames(res_HER2_TNBC))
signif_res_HT <- res_HER2_TNBC[res$padj < 0.05 & !is.na(res_HER2_TNBC$padj), ]
signif_genes_HT <- as.character(rownames(signif_res_HT))

ego_HT <- enrichGO(gene          = signif_genes,
                universe = all_genes_HT,
                keyType = "ENSEMBL",
                OrgDb         = organism,
                ont           = "BP" ,#"BP", #biological process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego_HT)
p1 <- dotplot(ego_HT, showCategory=20) + ggtitle("HER2 vs TNBC") 
  
barplot(ego_HT, showCategory=30) 

#H-NonTNBC
# Extract significant results
all_genes_HN <- as.character(rownames(res_HER2_NonTNBC))
signif_res_HN <- res_HER2_NonTNBC[res$padj < 0.05 & !is.na(res_HER2_NonTNBC$padj), ]
signif_genes_HN <- as.character(rownames(signif_res_HN))

ego_HN <- enrichGO(gene          = signif_genes_HN,
                universe = all_genes_HN,
                keyType = "ENSEMBL",
                OrgDb         = organism,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
p2 <- dotplot(ego, showCategory=20) + ggtitle("HER2 vs Non-TNBC") 

png("go.png", width     = 10,
    height    = 9,
    units     = "in",
    res       = 1200)
ggarrange(p1, p2, labels = c("A", "B"), ncol=2, nrow=1)
dev.off()
###########################################################################
####GSEA
##her2-nonTNBC
# log2 fold change 
original_gene_list_HN <- df_HN$log2FoldChange

# name the vector
names(original_gene_list_HN) <- row.names(df_HN)

# omit any NA values 
gene_list_HN<-na.omit(original_gene_list_HN)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_HN = sort(gene_list_HN, decreasing = TRUE)


#GO Gene Set Enrichment Analysis
gse_HN <- gseGO(geneList=gene_list_HN, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

p3 <- dotplot(gse_HN, showCategory=10, split=".sign", font.size = 10) + facet_grid(.~.sign) + ggtitle("HER2 vs Non-TNBC")

gseaplot(gse_HN, by = "all", title = gse_HN$Description[1], geneSetID = 1)

##her2-tnbc
original_gene_list_HT <- df_HT$log2FoldChange

# name the vector
names(original_gene_list_HT) <- row.names(df_HT)

# omit any NA values 
gene_list_HT<-na.omit(original_gene_list_HT)

# sort the list in decreasing order (required for clusterProfiler)
gene_list_HT = sort(gene_list_HT, decreasing = TRUE)


#GO Gene Set Enrichment Analysis
gse_HT <- gseGO(geneList=gene_list_HT, 
                ont ="ALL", 
                keyType = "ENSEMBL", 
                nPerm = 10000, 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none")

p4 <- dotplot(gse_HT, showCategory=10, split=".sign", font.size = 10) + facet_grid(.~.sign) + ggtitle("HER2 vs TNBC")

gseaplot(gse_HT, by = "all", title = gse_HN$Description[1], geneSetID = 1)

png("gsea.png", width     = 10,
    height    = 9,
    units     = "in",
    res       = 1200)
ggarrange(p4, p3, labels = c("A", "B"), ncol=2, nrow=1)
dev.off()


#####GSE KEGG
#her2-nontnbc
ids<-bitr(names(gene_list_HN), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df_HN[row.names(df_HN) %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
p5 <- dotplot(kk2, showCategory = 10, title = "HER2 vs non-TNBC" , split=".sign",font.size = 10) + facet_grid(.~.sign)

#her2-tnbc
ids<-bitr(names(gene_list_HT), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_idsHT = ids[!duplicated(ids[c("ENSEMBL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
dfHT = df_HT[row.names(df_HT) %in% dedup_idsHT$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
dfHT$Y = dedup_idsHT$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_listHT <- dfHT$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_listHT) <- dfHT$Y

# omit any NA values 
kegg_gene_listHT<-na.omit(kegg_gene_listHT)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_listHT = sort(kegg_gene_listHT, decreasing = TRUE)

kkHT <- gseKEGG(geneList     = kegg_gene_listHT,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
p6 <-dotplot(kkHT, showCategory = 10, title = "HER2 vs TNBC" , split=".sign", font.size = 10) + facet_grid(.~.sign)

png("kegg.png", width     = 10,
    height    = 9,
    units     = "in",
    res       = 1200)
ggarrange(p6, p5, labels = c("A", "B"), ncol=2, nrow=1)
dev.off()
