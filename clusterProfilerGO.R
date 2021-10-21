setwd("C:\\Users\\Javan\\Desktop\\NelsonSoares\\candidaProject\\DifferentialsPx")
#https://www.learnbyexample.org/r-bar-plot-ggplot2/#coloring-a-bar-graph
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html ; updated analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(pathview)
library(proteus)
library(org.Cauris.eg.db)
library(org.Mm.eg.db)
keytypes(org.Cauris.eg.db)

data <- read.csv("SA01-SC01.csv",header = T,sep = ',')

colnames(data)

data <- dplyr::select(data, X,EffectSize,pValue) ; dim(data)

data = subset(data,EffectSize >= 1.5 | pValue < 0.05 ) ;dim(data)#| EffectSize <= -1);dim(data)

#write.csv(data,"sa01-sc01.csv")

gene <- data$X# extract Gene names

# this translates the HUGO gene symbols to ENTREZID
gene.df <- bitr(gene, fromType = "UNIPROT", toType = "ENTREZID",OrgDb = org.Cauris.eg.db) ; dim(gene.df)

keytypes(org.Cauris.eg.db)

# Make a geneList for some future functions
geneList <- gene.df$ENTREZID
names(geneList) <- as.character(gene.df$UNIPROT)
geneList <- sort(geneList, decreasing = TRUE)

# gene enrichment analysis cnplots are commented out as they look crazy with a large number of proteins
## BP
ego_BP2 <- enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    readable = TRUE,
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

head(ego_BP2,10) #check the first 10 entries from ego_BP

df = as.data.frame(ego_BP2)

write.csv(df,"LTBI_B1vsPPD_GO_BP.csv")

ego2 <- simplify(ego_BP2) ; dim(ego2) # remove redundant GO terms first 

dotplot(ego2, showCategory=24,x="count",font.size = 9,title=" ")

#Barplot
barplot(ego2, 
        drop = TRUE, 
        showCategory = 20, 
        title = " ",
        font.size = 9,
        x="count")

#Treeplot,
edox2 <- pairwise_termsim(ego_BP2)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "median",nCluster=7)
aplot::plot_list(list(p1, p2)) + plot_annotation(tag_levels='A')

cnetplot(ego2, categorySize="pvalue", foldChange=geneList)


###############################GSEA########################

# reading in data from deseq2
df = read.csv("SA01-SB01.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$EffectSize

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#Run gsea
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "GENENAME", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Sc.sgd.db, 
             pAdjustMethod = "none")

keytypes(org.Sc.sgd.db)

#========================gprofiler2=======================
# install from CRAN
#https://f1000researchdata.s3.amazonaws.com/manuscripts/29655/97ec495b-9096-4146-b91b-4f17bb2e7ebd_24956_-_hedi_peterson_v2.pdf?doi=10.12688/f1000research.24956.2&numberOfBrowsableCollections=29&numberOfBrowsableInstitutionalCollections=4&numberOfBrowsableGateways=29
#install.packages("gprofiler2")
# load the package
library(gprofiler2)

#Load the data
data <- read.csv("C:\\Users\\Javan\\Desktop\\NelsonSoares\\cancerProject\\DifferentialPx\\PC-PD.csv",header = T,sep = ',')

colnames(data)

#subset the data
data <- dplyr::select(data, X,EffectSize,pValue) ; dim(data)

#Filter the most important proteins
data = subset(data,EffectSize >= 1.5 | pValue < 0.05 ) ;dim(data)#| EffectSize <= -1);dim(data)

names(data)[1] = "Gene_symbol"
names(data)[2] = "logFC"
names(data)[3] = "FDR_p"

#query the most the data
gostres = gost(query = c(data$X),
               organism = "hsapiens"
               )
gostplot(gostres, interactive = TRUE)

View(gostres)

p1 = gostplot(gostres, interactive = F)
publish_gostplot(p1, highlight_terms = c("GO:0072594", "GO:0051641","GO:0045047","GO:0006612",
                                         "GO:0051649","GO:0033365","GO:0003723"))

#==================================pathfindR=================================
install.packages("pathfindR")
library(pathfindR)

data <- read.csv("SA01-SB01.csv",header = T,sep = ',')

colnames(data)

data <- dplyr::select(data, X,EffectSize,pValue) ; dim(data)

data = subset(data,EffectSize >= 1.5 | pValue < 0.05 ) ;dim(data)

names(data)[1] = "Gene_symbol"
names(data)[2] = "logFC"
names(data)[3] = "FDR_p"

write.csv(data,"SA01_SB01.csv")

data = read.csv("SA01_SB01.csv",header = T,sep = ',')

output_df <- run_pathfindR(data)

























