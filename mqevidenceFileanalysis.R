library(proteus)
evidenceFile <- "C:\\Users\\Javan\\Desktop\\NelsonSoares\\fourthProject\\evidence.txt"
metadataFile <- "C:\\Users\\Javan\\Desktop\\NelsonSoares\\fourthProject\\met.txt"

evi <- readEvidenceFile(evidenceFile)
meta <- read.delim(metadataFile, header=TRUE, sep="\t")
pepdat <- makePeptideTable(evi, meta)
prodat <- makeProteinTable(pepdat)
prodat.med <- normalizeData(prodat)
res <- limmaDE(prodat.med,conditions = c("Tamoxifen_5uM","DMSO"))

plotVolcano_live(prodat.med, res)

plotIntensities(prodat.med, id='P37802', log=TRUE)

# For LFQ analysis
library(DEP)
library(MSnbase)
run_app("LFQ")

#================Boxplot=======================
library(dplyr)
setwd("C:\\Users\\Javan\\Desktop\\NelsonSoares\\fourthProject\\DEPs")

dmso_tamo = read.csv("dmso_tamoxifen.csv",header = T,sep = ",")

df = select(dmso_tamo,GeneName,EffectSize)

boxplot(df$EffectSize)

barplot(df$EffectSize,names.arg = df$GeneName,col = rainbow(6),
        main = "DMSO_Tamoxifen down and upregulated proteins",
        ylim = c(-3,3),horiz = F)
