setwd("C:\\Users\\Javan\\Desktop\\NelsonSoares\\candidaProject\\DifferentialsPx/")
## Prepare GO annotation
library(mgsa)
library(enrichplot)
library(clusterProfiler)

# read GO annotation file
# Note that I manually downloaded the *.goa file, and is in working directory.
GAF <- readGAF(filename="candidauris.goa")

# extract relevant info.
# unfortunately could not find accessor functions for all required
# info, thus sometimes had to utilize object slots directly (with @)
mapping.index <-  GAF@itemName2ItemIndex
ID.annotations <- itemAnnotations(GAF)

GO.sets <- GAF@sets
GO.annotation <- setAnnotations(GAF)

# create a 2-column data frame with GOID and ID index
# after little further processing, this will be used as input for clusterProfiler
GO.df <- data.frame("GOID" = rep(names(GO.sets), sapply(GO.sets, length)),
                    "ID.index" = unlist(GO.sets),  row.names = NULL)


# do some processing for objects GO.annotation and GO.df
# in both remove category 'all',
# and to GO.df also add column with Uniprot ids

# GO.annotation
GO.annotation <- GO.annotation[GO.annotation[,"term"] != "all", ]
GO.annotation[,"GOID"] <- rownames(GO.annotation)

# GO.df
GO.df <- GO.df[GO.df[,"GOID"] != "all", ]
GO.df[,"UNIPROTKB"] <- names(mapping.index [GO.df[,"ID.index"] ])


## objects are now ready for use with clusterProfiler!


## Functional analysis with clusterProfiler.
#All proteins
df1 <- read.csv("SE01-SD01.csv",header = T,sep = ','); dim(df1)

df2 = subset(df1,EffectSize >= 1.5 | pValue < 0.05);dim(df2)

# Perform GO ORA analysis
# Note the use of arguments TERM2GENE and TERM2NAME. Their column
# order is important!
# Also: to show proof-of-principle no cutoff on significance
# is applied! i.e Cutoff=1

res.GO.ora <- enricher(
  gene=df2$X,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe=ID.annotations,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 1,
  TERM2GENE = GO.df[ ,c("GOID","UNIPROTKB")],
  TERM2NAME = GO.annotation[ ,c("GOID", "term")]
)

# check
as.data.frame(res.GO.ora)[1:15,]

#Generate plots
barplot(res.GO.ora, showCategory=30,x="count",font.size = 9)

dotplot(res.GO.ora, showCategory=25,x="count",font.size = 9)


# Also perform KEGG-based analyses.
# note: only 10841 (~52%) uniprot ids have been annotated by KEGG!
uni2kegg <- bitr_kegg( df1, fromType="uniprot",
                       toType='kegg', organism='Cau')

# KEGG ORA.
kk <- enrichKEGG(
  gene =  proteins.random,
  organism = "sly",
  keyType = "uniprot",
  pvalueCutoff = 1,
  pAdjustMethod = "none",
  universe = proteins.all,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 1
)

# check
as.data.frame(kk)[1:15,]

# object kk can be used to generate same plots as above!


# GSEA with KEGG pathways.
# for GSEA a per-gene metric is required that is used for ranking.
# assume dataset consists of 17.5k uniprot ids.
# make up fold changes (= ranking metric) ranging from -15 to 15.

protList <- runif(n = 17500, min = -15, max = 15)
names(protList) <- base::sample( proteins.all, 17500)
# for GSEA, sort on FC
protList <- sort(protList, decreasing = TRUE)

# KEGG GSEA
kk2 <- gseKEGG(geneList = protList,
               organism     = "sly",
               keyType      = "uniprot",
               eps          = 0.0,
               minGSSize    = 10,
               maxGSSize    = 500,
               pAdjustMethod = "none",
               pvalueCutoff = 0.9,
               verbose      = FALSE)

# check
as.data.frame(kk2)[1:15,]

# again, object kk2 can be used to generate same plots as above!

# Lastly, compare multiple clusters using KEGG.
# First generate list with multiple clusters
clusterData <- list(
  "Cluster 1" = base::sample(proteins.all, 19500),
  "Cluster 2" = base::sample(proteins.all, 18900),
  "Cluster 3" = base::sample(proteins.all, 22000),
  "Cluster 4" = base::sample(proteins.all, 20500) )

# run compareCluster
xx <- compareCluster(clusterData, fun="enrichKEGG",
                     organism="sly",  keyType = "uniprot",
                     pAdjustMethod = "none", pvalueCutoff=1.0)

# plot
emapplot( pairwise_termsim(xx) )