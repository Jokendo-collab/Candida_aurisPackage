setwd("C:\\Users\\Javan\\Desktop\\UFS_collaboration\\SA_provinces\\lineages") #set working directory
library(dplyr)
library(ggplot2)
library(cowplot)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
#Load the first wave data
firstwave = read.csv("western_cape_wave3.csv",header = T, sep = ',')
firstwave
colnames(firstwave)

firstwave = select(firstwave,Lineage) #grab the Lineage column

df = table(firstwave$Lineage) #create the table with frequency of observations

df = as.data.frame(df) #convert the table into a dataframe

View(df)
df = df[-c(28), ]

names(df)[1] = "Variant_type"
names(df)[2] = "Infected_individuals"

df = subset(df,Infected_individuals >= 2 ) ;dim(df)

#write.csv(df,"wcp3.csv") #Write the summarized variant types

# Barplot
ggplot(df, aes(x=Variant_type, y=Infected_individuals,fill=Variant_type)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle(" ")+
  xlab("Variant types") +
  ylab("Number of infected individuals") 

#Density plot
df = read.csv("C:\\Users\\Javan\\Desktop\\UFS_collaboration\\SA_provinces\\SNPs\\provincialSpikeMutations.csv",header = T,sep = ',')
# Change density plot line colors by groups
ggplot(data=df, aes(x=Frequency, group=Province, fill=Province)) +
  geom_density(adjust=1.5) +
  theme_cowplot() +
  facet_wrap(~Province) 

ggplot(data=df, aes(x=Frequency, group=Province, fill=Province)) +
  geom_density(adjust=1.5,position = "fill") +
  theme_cowplot() +
  xlab("Spike protein mutation distribution") +
  ggtitle("SARS-COV-2 spike protein mutations in South Africa")

#Second wave dataframe
secondwave = read.csv("secondwave_November2020-March2021.csv",header = T, sep = ',')
secondwave

secondwave = select(secondwave,Lineage) #grab the Lineage column

df2 = table(secondwave$Lineage) #create the table with frequency of observations

df2 = as.data.frame(df2) #convert the table into a dataframe

df2 = df2[-c(5), ]

names(df2)[1] = "Variant_type"
names(df2)[2] = "Infected_individuals"

#write.csv(df2, "second_wave.csv")
# Barplot
ggplot(df2, aes(x=Variant_type, y=Infected_individuals,fill=Variant_type)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle("Second wave variant types in Free state")+
  xlab("Variant types") +
  ylab("Number of infected individuals") 

setwd("C:\\Users\\Javan\\Desktop\\UFS_collaboration")
#Third wave
thirdwave = read.csv("thirdwave_April 2021_September2021.csv",header = T, sep = ',')
thirdwave

thirdwave = select(thirdwave,Lineage) #grab the Lineage column

df3 = table(thirdwave$Lineage) #create the table with frequency of observations

df3 = as.data.frame(df3) #convert the table into a dataframe

df3 = df3[-c(12), ]

names(df3)[1] = "Variant_type"
names(df3)[2] = "Infected_individuals"
#write.csv(df3,"third_wave.csv")

# Barplot
ggplot(df3, aes(x=Variant_type, y=Infected_individuals,fill=Variant_type)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle("Third wave variant types in Free state")+
  xlab("Variant types") +
  ylab("Number of infected individuals") 

#Getting variant trend in South Africa
vaNtdata = read.csv("GT_province.csv",header = T,sep = ',')
vaNtdata = subset(vaNtdata,Frequency >= 50) ;dim(vaNtdata)

#variants classification by waves
ggplot(vaNtdata, aes(Frequency, Amino.acid.mutations, color = Wave)) + 
  geom_point(size=3.5) +
  theme_cowplot() +
  ggtitle(" ")+
  xlab("SARS-COV-2 variants frequency") +
  ylab("Spike protein mutations") +
  theme(axis.text = element_text(size = 8)) +
  facet_wrap(~Wave)

#####################################
ggplot(vaNtdata, aes(Frequency, Amino.acid.mutations, fill = Wave)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size = 7)) +
  ylab("Spike protein mutations") +
  xlab("Variant frequency") +
  ggtitle(" ") +
  facet_wrap(~Wave)

#samples which were analyzed from nine provinces
dd = read.csv("sampleNumber.csv",header = T,sep = ",")
ggplot(dd, aes(Analyzed.samples, Provincial.waves, fill = Provincial.waves)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size = 10)) +
  theme_cowplot()

#What variants are common and different
library(BioVenn)
df = read.csv("western_cape_province.csv",header = T,sep = ',')

First_wave = filter(df, Wave=="First")
First_wave = First_wave$Variant_type
Second_wave = filter(df, Wave=="Second")
Second_wave = Second_wave$Variant_type
Third_wave = filter(df, Wave=="Third")
Third_wave = Third_wave$Variant_type

biovenn <- draw.venn(First_wave, Second_wave, Third_wave, 
                     subtitle=" ", nrtype="abs",
                     title = "Western Cape",
                     xtitle = "First wave",
                     ytitle = "Second wave",
                     ztitle = "Third wave")

#Variant name and variant class analysis
setwd("C:\\Users\\Javan\\Desktop\\UFS_collaboration\\variantYpes")
f1varname = read.csv("firstWaveVariantTypes.csv",header = T, sep = ',')
f1varname = na.omit(f1varname)

f1varname = select(f1varname,varname) #grab the varname column

df = table(f1varname$varname) #create the table with frequency of observations

df = as.data.frame(df) #convert the table into a dataframe

df = subset(df,Freq >= 15 ) ;dim(df)

names(df)[1] = "Variant_name"
names(df)[2] = "Frequency"

#write.csv(df,"first_wave.csv") #Write the summarized variant types

# Barplot
ggplot(df, aes(x=Variant_name, y=Frequency,fill=Variant_name)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle(" ")+
  xlab("Viral protein mutations") +
  ylab("Mutations frequency")+
  theme(axis.text = element_text(size = 8)) 

#Second wave amino acids substitution
secondwave = read.csv("secondwavevarianttypes.csv",header = T, sep = ',')
f2varname = na.omit(secondwave)

f2varname = select(secondwave,varname) #grab the varname column

df2 = table(f2varname$varname) #create the table with frequency of observations

df2 = as.data.frame(df2) #convert the table into a dataframe

df2 = subset(df2,Freq >= 15 ) ;dim(df)

names(df2)[1] = "Variant_name"
names(df2)[2] = "Frequency"

#write.csv(df,"first_wave.csv") #Write the summarized variant types

# Barplot
 ggplot(df2, aes(x=Variant_name, y=Frequency,fill=Variant_name)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle(" ")+
  xlab("Viral protein mutations") +
  ylab("Mutations frequency")+
  theme(axis.text = element_text(size = 8)) 

#third wave
thirdwave = read.csv("thirdwavevariantTypes.csv",header = T, sep = ',')
f3varname = na.omit(thirdwave)

f3varname = select(thirdwave,varname) #grab the varname column

df3 = table(f3varname$varname) #create the table with frequency of observations

df3 = as.data.frame(df3) #convert the table into a dataframe

df3 = subset(df3,Freq >= 200 ) ;dim(df)

names(df3)[1] = "Variant_name"
names(df3)[2] = "Frequency"

#write.csv(df,"first_wave.csv") #Write the summarized variant types

# Barplot
ggplot(df3, aes(x=Variant_name, y=Frequency,fill=Variant_name)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle(" ")+
  xlab("Viral protein mutations") +
  ylab("Mutations frequency")+
  theme(axis.text = element_text(size = 8)) 

plot_grid(p1, p2,p3, labels = c('A', 'B',"C"))

grid.arrange(p1, p2,p3, nrow = 1, ncol = 3)


setwd("C:\\Users\\Javan\\Desktop\\UFS_collaboration\\variantYpes")
#Spike proteins mutations from UFS data
f1varname = read.csv("thirdwavevariantTypes.csv",header = T, sep = ',')
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

df = subset(df,Freq >= 10) ;dim(df)

names(df)[1] = "Amino.acid.mutations"
names(df)[2] = "Frequency"

# Barplot
ggplot(df, aes(x=Amino.acid.mutations, y=Frequency,fill=Amino.acid.mutations)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle(" ")+
  xlab("Viral spike protein mutations") +
  ylab("Mutation frequency")+
  theme(axis.text = element_text(size = 8))
