#Processing and initial differential expression analysis


#Libraries
library(readr)
library(tximport)
library(BiocManager)
library(edgeR)
library(limma)
library(DESeq2)
library(ctc)
library(Biobase)
library(gplots)
library(ape)
library(argparse)
library(pheatmap)
library(viridis)
library(tidyverse)
library(EnhancedVolcano)

setwd("C:/Users/teres/OneDrive/Documents/R/Data copy/data/Metatranscriptomics")
#####Initial Steps#####
###Preparing counts matrix
samples<-read_tsv("samples.txt") #samples.txt file is read in to get the sample ID column
files<-file.path(getwd(),"SalmonFiles",samples$ID,"quant.sf.gz") #create file path for all quant.sf files, and give them the sample IDs
names(files) <- (samples$ID)

#Read in file matching the transcript/isoform level information to gene level information
tx2gene<-read_tsv("SalmonFiles/Trinity.fasta.gene_trans_map", #this is bc we want to have gene level differential expression analysis
                  col_names = c("Gene","Transcript"))%>%
  relocate("Gene",.after="Transcript")

#The different lengths of isoforms can create bias in the data
txi.salmon.corrected <- tximport(files, type = "salmon", 
                                 tx2gene = tx2gene, #uses the gene/transcript mapping to give gene level results
                                 countsFromAbundance="lengthScaledTPM") #can correct the isoform length bias here
#result is suitable for DE with limna or edgeR, but not with DESeq2

#With everything below, 1 transcript = 1 gene. Isoforms have been accounted for
transcript_counts<- round(txi.salmon.corrected$counts)%>% #Round counts to nearest whole number
  as.data.frame()%>%
  mutate(transcript_id= rownames(txi.salmon.corrected$counts), .before= "BP_OG1") #New column, "gene", placing the transcript ID as the first column 
transcript_counts%>%  
  write_tsv("transcript_counts.csv") #Save file for future reference

###CreateDGE (Differential Gene Expression) object

rownames(transcript_counts)<- transcript_counts$transcript_id #make the transcript IDs the rownames
transcript_counts<- transcript_counts%>%select(-transcript_id) #remove the transcript ID column

#read in expdesign.csv
expdesign<-read_csv("expdesign.csv")%>% #experimental Design = sample metadata indicating which group samples are in
  as.data.frame() 

new_order = sort(colnames(transcript_counts))#sort cols of transcript counts to be in the same order as the rows in expdesign
transcript_counts <- transcript_counts[, new_order] 

rownames(expdesign)<-colnames(transcript_counts) #set the col names of transcript_counts (which are the sample names), to be the row_names of exp_design
view(expdesign)
view(transcript_counts) #make sure these both look correct before moving on

#Set data groups. Groups here consist of 2 factors: the spring that the sample came from, and the condition of flow
#flow state for stable springs is always labelled Constant, while for geysing springs, it can be On or Off
datagroups<- factor(paste(expdesign$Spring,expdesign$Flow, sep="."))
cbind(expdesign,Group=datagroups) #there is now a new column, where variables have the format SpringID.FlowState

d <- DGEList(counts=transcript_counts,group=factor(datagroups)) #create DGE object
d.full <- d #create copy of original DGE object 
dim(d.full) #352459 transcripts and 20 samples

#Filter out lowly expressed transcripts because these can skew analysis
#Common cutoff is to keep transcripts at least 10 counts per million (CPM)
apply(d$counts, 2, sum) 
keep <- rowSums(cpm(d)>10) >= 2
d <- d[keep,]
dim(d) #left with 23775 genes

#After filtering, reset the library size information
d$samples$lib.size <- colSums(d$counts)
d$samples

##Normalization:
#calcNormFactors() function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes
#default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
d <- calcNormFactors(d)
d
#adding variables to the object
d$samples<-d$samples%>% 
  mutate(flow= factor(paste(expdesign$Flow)))%>%mutate(spring=factor(paste(expdesign$Spring)))

#data exploration with an MDSplot
plotMDS(d, method="logFC", col = as.numeric(d$samples$group))

##Estimating the Dispersion
#The common dispersion
#Estimating the common dispersion gives an idea of overall variability across the dataset.
d1 <- estimateCommonDisp(d, verbose=T) 
d1 <- estimateTagwiseDisp(d1)
# BCV = biological coefficient of variation estimate by sqrt(dispersion)
plotBCV(d1)

## GLM will be used to fit dispersion better; prior to differential expression
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)

d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto") 
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)


##### Differential expression contrasts#####
#the glm identified model of dispersion is used to find differentially expressed genes
fit.glm <- glmFit(d2, design.mat)

#Contrasts of interest
glm.contrasts<-makeContrasts(levels=design.mat,
  GSvSB=(JJ.On+FC.On+JJ.Off+FC.Off)/4-(BP.Constant+PR.Constant+RC1.Constant+RC2.Constant)/4, #geysing vs stable
  JJvFC=(JJ.On+JJ.Off)/2 - (FC.On+FC.Off)/2,
  JJ.OnvOff=JJ.Off-JJ.On,
  FCvBP=(FC.On+FC.Off)/2-BP.Constant, #Bison Pool is the closest physically to the 2 geysing springs
  JJvBP=(JJ.On+JJ.Off)/2-BP.Constant,
  RCAvBP=(RC1.Constant-BP.Constant)
  )

#GEYSING VS STABLE
GSvSB.glm <- glmLRT(fit.glm, contrast=glm.contrasts[,"GSvSB"])
groupTags <- topTags(GSvSB.glm, n = Inf)
plotMD(GSvSB.glm) #volcano plot
abline(h=c(-2, 2), col="blue")
GSvSB.glmDE <- decideTestsDGE(GSvSB.glm, adjust.method = "fdr", p.value=0.1)
summary(GSvSB.glmDE) 
#create object of all the transcripts, along with results of DE of stable vs geysing
allSvG<-groupTags[["table"]]

allSvG<-allSvG%>%   
  mutate(gene=rownames(allSvG))

protein_info<-read.csv("bio_meta.csv") #bio_meta csv has the UniProt IDs for the transcripts

u_protein_info<- protein_info%>% 
  select(-X)%>%
  group_by(gene) %>% 
  filter(max_pid == max(max_pid))

allSvG<-left_join(allSvG, u_protein_info, by="gene")
write_csv(allSvG, "allGenes_with_SvG_pvals.csv") #write the file for visualizing results 


#JJvFC=(JJ.On+JJ.Off)/2 - (FC.On+FC.Off)/2,#####
JJvFC.glm <- glmLRT(fit.glm, contrast=glm.contrasts[,"JJvFC"])
topTags(JJvFC.glm, n = 5) #n="Inf"
plotMD(JJvFC.glm, legend=F)
abline(h=c(-2, 2), col="blue")
JJvFC.glmDE <- decideTestsDGE(JJvFC.glm, adjust.method = "fdr", p.value=0.1)
summary(JJvFC.glmDE)

####JJ.OnvOff####
JJ.OnvOff.glm <- glmLRT(fit.glm, contrast=glm.contrasts[,"JJ.OnvOff"])
topTags(JJ.OnvOff.glm, n = 5) #n="Inf"
plotMD(JJ.OnvOff.glm)
abline(h=c(-2, 2), col="blue")
JJ.OnvOff.glmDE <- decideTestsDGE(JJ.OnvOff.glm, adjust.method = "fdr", p.value=0.1)
summary(JJ.OnvOff.glmDE)

groupTags[["table"]] %>%   
  mutate(gene=rownames(.))%>%
  left_join(., u_protein_info, by="gene")%>%
  write_csv(., "allGenes_with_JJ_OnvOff_pvals.csv") 


#FCvBP####
FCvBP.glm <- glmLRT(fit.glm, contrast=glm.contrasts[,"FCvBP"])
topTags(FCvBP.glm, n = 5) #n="Inf"
plotMD(FCvBP.glm)
abline(h=c(-2, 2), col="blue")
FCvBP.glmDE <- decideTestsDGE(FCvBP.glm, adjust.method = "fdr", p.value=0.1)
summary(FCvBP.glmDE)
#JJvBP####
JJvBP.glm <- glmLRT(fit.glm, contrast=glm.contrasts[,"JJvBP"])
topTags(JJvBP.glm, n = 5) #n="Inf"
plotMD(JJvBP.glm)
abline(h=c(-2, 2), col="blue")
JJvBP.glmDE <- decideTestsDGE(JJvBP.glm, adjust.method = "fdr", p.value=0.1)
summary(JJvBP.glmDE)

#RCAvBP####
RCAvBP.glm <- glmLRT(fit.glm, contrast=glm.contrasts[,"RCAvBP"])
topTags(RCAvBP.glm, n = 5) #n="Inf"
plotMD(RCAvBP.glm)
abline(h=c(-2, 2), col="blue")
RCAvBP.glmDE <- decideTestsDGE(RCAvBP.glm, adjust.method = "fdr", p.value=0.1)
summary(RCAvBP.glmDE)
