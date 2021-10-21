setwd("C:\\Users\\Javan\\Desktop\\UFS_collaboration\\SA_provinces\\SNPs") #set working directory
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(BioVenn)
library(cowplot)
library(viridis)
library(hrbrthemes)
library(grid)      # combining plots
library(gridExtra) # combining plots
library(ggpubr)    # combining plots
library(patchwork) # combining plots
#Load the first wave data
firstwave = read.csv("results.csv",header = T, sep = ',')
firstwave

firstwave = select(firstwave,Lineage) #grab the Lineage column

df = table(firstwave$Lineage) #create the table with frequency of observations

df = as.data.frame(df) #convert the table into a dataframe

#df = df[-c(8), ] #delete the column with none variants

names(df)[1] = "Variant_type"
names(df)[2] = "Infected_individuals"

#write.csv(df,"wv_3.csv")

#write.csv(df,"first_wave.csv") #Write the summarized variant types

# Barplot
ggplot(df, aes(x=Variant_type, y=Infected_individuals,fill=Variant_type)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle(" ")+
  xlab("Variant types") +
  ylab("Number of infected individuals") 

#Getting variant trend
vaNtdata = read.csv("./wv_1.csv",header = T,sep = ',')

# Stacked 
ggplot(vaNtdata, aes(fill=Wave, y=Infected_individuals, x=Variant_type)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_cowplot() +
  ylab("Number of infected individuals") +
  facet_wrap(~Wave) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text = element_text(size = 8))

#Spike proteins mutation
f1varname = read.csv("GT_wv3.csv",header = T, sep = ',')
f1varname = na.omit(f1varname)

#Droping mutations on the other proteins keeping spike protein
f1varname = f1varname[!grepl("ORF", f1varname$protein),]
f1varname = f1varname[!grepl("NSP",f1varname$protein),]
f1varname = f1varname[!grepl("E",f1varname$protein),]
f1varname = f1varname[!grepl("N",f1varname$protein),]
f1varname = f1varname[!grepl("M",f1varname$protein),]

f1varname = select(f1varname,varname) #grab the varname column

df = table(f1varname$varname) #create the table with frequency of observations

df = as.data.frame(df) #convert the table into a dataframe

df = subset(df,Freq >= 80) ;dim(df)

names(df)[1] = "Amino.acid.mutations"
names(df)[2] = "Frequency"

#write.csv(df,"GT_prov_w3.csv") #Write the summarized variant types

# Barplot
ggplot(df, aes(x=Amino.acid.mutations, y=Frequency,fill=Amino.acid.mutations)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle("North West wave 3 SARS-COV-2 protein mutation")+
  xlab("Viral spike protein mutations") +
  ylab("Mutations frequency")+
  theme(axis.text = element_text(size = 8))

#============Provincial spike mutations==============
ggplot(vaNtdata, aes(Frequency, Amino.acid.mutations, fill = Wave)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size = 8)) +
  ylab("Spike protein mutations") +
  xlab("Variant frequency") +
  ggtitle(" ") +
  facet_wrap(~Wave) +
  ggtitle(" ")

















