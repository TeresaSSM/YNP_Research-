#Creating Biological/Gene metadata
library(tidyverse)
prots<-read.delim("swissprot.blastx.outfmt6", header=F)%>%
  select(c(1:3))%>%
  separate_wider_delim(cols=V2,delim="|", names=c("RM","ACC","ID"))%>%
  separate_wider_delim(cols=V1, delim="_i", names=c("gene","isoform"))%>%
  select(-c(RM,isoform))%>%
  unique()%>%
  group_by(gene, ACC,ID)%>%
  summarize(max_pid= max(V3))%>%
  ungroup()%>%
  group_by(gene)%>%
  slice(which.max(max_pid))%>%
  ungroup()
prots1<-prots

add_to_prots<-subset(transcript_counts$gene,
                     !(as.vector(transcript_counts$gene) %in% as.vector(prots$gene)))%>%
  data.frame(NA,NA,NA)
colnames(add_to_prots)<-c( "gene"   , "ACC" ,    "ID"  ,    "max_pid")

prots2<-rbind(prots1,add_to_prots)
write.csv(prots2, "bio_meta.csv")