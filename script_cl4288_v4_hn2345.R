rm(list=ls())

library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ReactomePA)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

setwd("/ifs/data/c2b2/ac_lab/hn2345/project_sars2classification/shared_cl4288/data")

metdat <- read.table("metadata_cumc_hn2345.txt", header = TRUE)
count <- readRDS("count_cumc_hn2345.rds")
tpm <- readRDS("tpm_cumc_hn2345.rds")


#=================================================
## Filtering for protein-coding genes in the TPM data (Reducing genes in the gene expresion matrix will save much more memory and computational time )

anno = readRDS("anno_grch38_v2.rds")
gc = intersect(rownames(tpm), anno$hgnc_symbol[anno$transcript_biotype=="protein_coding"])
length(gc)
tpm = tpm[match(gc, rownames(tpm)),]
dim(tpm)

## Filter out the genes having no expression
tpm = tpm[rowSums(tpm)>0,]
dim(tpm)

tpm <- log2(tpm + 1)
#=================================================


# Label each sample with its diagnosis

# metdat$diagnosis <- ifelse(metdat$viral.load == 'None', "Negative", ifelse(metdat$intubation == 1 | metdat$mortality == 1 | (metdat$days_hospitalized >= 7 & !is.na(metdat$days_hospitalized)), "Severe", "Mild"))

metdat$diagnosis <- ifelse(metdat$viral.load == 'None' & metdat$n_postivie==0, "Negative", ifelse(metdat$intubation == 1 | metdat$mortality == 1 | (metdat$days_hospitalized >= 7 & !is.na(metdat$days_hospitalized)), "Severe", "Mild"))




# metdat$diagnosis[metdat$diagnosis!="Negative" & metdat$n_postivie==0]=NA
metdat$diagnosis[is.na(metdat$age)]=NA
metdat$diagnosis[is.na(metdat$days_hospitalized) & metdat$diagnosis!="Negative"]=NA

metdat$diagnosis[metdat$diagnosis=="Mild" & metdat$days_hospitalized>3]=NA



#Only keep the earlist sample
## Corrected
metdat$sample.id2 = as.numeric(gsub("SS","",metdat$sample.id))
metadata_patient <- metdat %>%
  group_by(case.id) %>%
  filter(sample.id2 == min(sample.id2))


table(metadata_patient$diagnosis)
# Mild Negative   Severe 
# 103      133       83 

#=================================================

# #Only keep the earlist sample
# metadata_patient <- metdat %>%
#   group_by(case.id) %>%
#   filter(sample.id == min(sample.id))

#--------------------------------------

metadata_patient=metadata_patient[!is.na(metadata_patient$diagnosis),]
tpm = tpm[,metadata_patient$sample.id]
dim(metadata_patient)
dim(tpm)

## Corrected:

unique_diagnoses <- unique(metadata_patient$diagnosis)


# subset tpm by diagnosis
tpm_mild <- as.matrix(tpm[, which(metadata_patient$diagnosis == "Mild")])
tpm_severe <- as.matrix(tpm[, which(metadata_patient$diagnosis == "Severe")])
tpm_negative <- as.matrix(tpm[, which(metadata_patient$diagnosis == "Negative")])


#===================================
library(matrixStats) ## This R pakcage included many functions, such as rowMeans, colMeans, rowSds and colSds to compute in a matrix input

tpm_mild_mean <- rowMeans(tpm_mild)
tpm_mild_sd <- rowSds(tpm_mild)
tpm_severe_mean <- rowMeans(tpm_severe)
tpm_severe_sd <- rowSds(tpm_severe)
tpm_negative_mean <- rowMeans(tpm_negative)
tpm_negative_sd <- rowSds(tpm_negative)

dim(tpm_mild)
dim(tpm_negative)
dim(tpm_severe)
#===================================

# # combine mean expression values into a new dataframe
# tpm_by_diagnosis <- data.frame(tpm_mild_mean, tpm_severe_mean, tpm_negative_mean)
# 
# 
# # Split the TPM data into a list of dataframes, one for each diagnosis
# tpm_by_diagnosis <- lapply(unique_diagnoses, function(diag) {
#   # Get the case IDs for this diagnosis
#   case_ids <- metdat$case.id[metdat$diagnosis == diag]
#   # Get the TPM data for these case IDs
#   tpm_subset <- tpm[, colnames(tpm) %in% case_ids]
#   # Return the subset of TPM data
#   return(tpm_subset)
# })
# 


# 
# #z-score
# severe_negative = (tpm_severe_mean - t(tpm_negative_mean))/ (tpm_negative_sd + 0.00001)
# hist(severe_negative, breaks=100)
# 
# mild_negative = (tpm_mild_mean - t(tpm_negative_mean)) / (tpm_negative_sd + 0.00001)
# hist(mild_negative, breaks=100)
#================================================

## two-sample z-score
# https://statmagic.info/Content/Help-Content/two-sample-mean.html

## stdev ==0 came from 0 counts
# range(tpm_severe_mean[tpm_severe_sd==0])
# range(tpm_severe_sd)
# range(tpm_negative_sd)

tpm_negative_sd[tpm_negative_sd<0.01]=0
tpm_mild_sd[tpm_mild_sd<0.01]=0
tpm_severe_sd[tpm_severe_sd<0.01]=0

head(tpm_severe_mean)
head(tpm_negative_mean)
head(tpm_mild_mean)

## Z-score by two-sample z -test, comparing severe group and mild group
severe_negative = (tpm_severe_mean - t(tpm_negative_mean))/ sqrt(tpm_severe_sd^2/ncol(tpm_severe) + tpm_negative_sd^2/ncol(tpm_negative) )

mild_negative = (tpm_mild_mean - t(tpm_negative_mean))/ sqrt(tpm_mild_sd^2/ncol(tpm_mild) + tpm_negative_sd^2/ncol(tpm_negative) )

severe_mild = (tpm_severe_mean - t(tpm_mild_mean))/ sqrt(tpm_severe_sd^2/ncol(tpm_severe) + tpm_mild_sd^2/ncol(tpm_mild) )


## Removing NA or Inf values in the z-scores
severe_negative[is.infinite(severe_negative)]=NA
mild_negative[is.infinite(mild_negative)]=NA
severe_mild[is.infinite(severe_mild)]=NA

# Identify columns without NA values in each dataframe
severe_negative_columns <- colnames(severe_negative)[!apply(severe_negative, 2, function(x) any(is.na(x)))]
mild_negative_columns <- colnames(mild_negative)[!apply(mild_negative, 2, function(x) any(is.na(x)))]
severe_mild_columns <- colnames(severe_mild)[!apply(severe_mild, 2, function(x) any(is.na(x)))]

# Find the intersection of column names without NA values
gc3 <- Reduce(intersect, list(severe_negative_columns, mild_negative_columns, severe_mild_columns))

# Print the result
# print(gc3)

severe_negative = severe_negative[,gc3]
mild_negative = mild_negative[,gc3]
severe_mild = severe_mild[,gc3]

head(severe_negative)
head(mild_negative)
head(severe_mild)


# hist(severe_mild, breaks=100)
# hist(severe_negative, breaks=100)
# hist(mild_negative, breaks=100)
# 
# plot(severe_negative, mild_negative);
# lines(c(-10,10),c(0,0), col="blue")
# lines(c(0,0),c(-10,10), col="blue")


#================================================


ori_gene_list <- severe_negative
# print(ori_gene_list)
names(ori_gene_list) <- names(severe_negative)
gene_list_<-na.omit(ori_gene_list)
gene_list_<-sort(gene_list_, decreasing = TRUE)



ori_gene_list2 <- mild_negative
# print(ori_gene_list2)
names(ori_gene_list2) <- names(mild_negative)
gene_list_2<-na.omit(ori_gene_list2)
gene_list_2<-sort(gene_list_2, decreasing = TRUE)


ori_gene_list3 <- severe_mild
# print(ori_gene_list3)
names(ori_gene_list3) <- names(severe_mild)
gene_list_3<-na.omit(ori_gene_list3)
gene_list_3<-sort(gene_list_3, decreasing = TRUE)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"

#================================================

## Load annotation table
anno = readRDS("./anno_grch38_v2.rds")
head(anno)
## protein-coding genes:
gene.pc = intersect(names(gene_list_), anno$hgnc_symbol[anno$transcript_biotype=="protein_coding"])
gene.pc_2 = intersect(names(gene_list_2), anno$hgnc_symbol[anno$transcript_biotype=="protein_coding"])
gene.pc_3 = intersect(names(gene_list_3), anno$hgnc_symbol[anno$transcript_biotype=="protein_coding"])

## new gene list only for the protein coding genes
gene_list_ = gene_list_[match(gene.pc, names(gene_list_))]
head(gene_list_)
gene_list_2 = gene_list_2[match(gene.pc_2, names(gene_list_2))]
head(gene_list_2)
gene_list_3 = gene_list_3[match(gene.pc_3, names(gene_list_3))]
head(gene_list_3)

## convert the gene symbol to entrez ids
names(gene_list_) = anno$entrezgene_id[match(names(gene_list_), anno$hgnc_symbol)]
names(gene_list_2) = anno$entrezgene_id[match(names(gene_list_2), anno$hgnc_symbol)]
names(gene_list_3) = anno$entrezgene_id[match(names(gene_list_3), anno$hgnc_symbol)]

length(gene_list_)
length(gene_list_2)
length(gene_list_3)


# gene_list_
gene_list_<-sort(gene_list_, decreasing = TRUE)

# gene_list_2
gene_list_2<-sort(gene_list_2, decreasing = TRUE)

# gene_list_3
gene_list_3<-sort(gene_list_3, decreasing = TRUE)


#==============================================
gse <- gseGO(geneList=gene_list_, 
             ont ="BP", ############### 
             keyType = "ENTREZID", 
             nPerm = 1000, 
             minGSSize = 10, ####
             maxGSSize = 200, ####
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

gse=as.data.frame(gse)
head(gse[,c(2,5:8)],n=20)

# dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

head(gene_list_2)
gse2 <- gseGO(geneList=gene_list_2, 
              ont ="BP", ############### 
              keyType = "ENTREZID", 
              nPerm = 1000, 
              minGSSize = 10, 
              maxGSSize = 200, 
              pvalueCutoff = 1, 
              verbose = TRUE, 
              OrgDb = organism, 
              pAdjustMethod = "BH")

gse2=as.data.frame(gse2)
head(gse2[,c(2,5:8)],n=50)



#dotplot(gse2, showCategory=10, split=".sign") + facet_grid(.~.sign)

head(gene_list_3)
gse3 <- gseGO(geneList=gene_list_3, 
              ont ="BP", ############### 
              keyType = "ENTREZID", 
              nPerm = 1000, 
              minGSSize = 10, 
              maxGSSize = 200, 
              pvalueCutoff = 1, 
              verbose = TRUE, 
              OrgDb = organism, 
              pAdjustMethod = "BH")

gse3=as.data.frame(gse3)
head(gse3[,c(2,5:8)],n=20)



##############################################################################################################################
#hallmark

library(viper)

#--------------------------------------------
## Load the hall mark data
hm = read.delim("./HallMarks_GeneSets.csv",sep=",",header=F)

# Constructing a regulon object
regulon.obj = apply(hm,1,function(x){
  gset = as.character(x[-1])
  gset = gset[gset!=""]
  r = list(tfmode=rep(1,length(gset)), likelihood=rep(1,length(gset)))
  names(r$tfmode) = gset
  return(r)
})
names(regulon.obj) = hm[,1]

names(gene_list_3) = anno$hgnc_symbol[match(names(gene_list_3), anno$entrezgene_id)]

gsea.hm = aREA(eset=gene_list_3, regulon = regulon.obj)$nes
gsea.hm


##############################################################################################################################
res <- gseWP(geneList= gene_list_3,
             organism="Homo sapiens",
             # keyType = "ENTREZID",
             minGSSize = 10, 
             maxGSSize = 200, 
             # OrgDb = organism, 
             pvalueCutoff = 1)
res=as.data.frame(res)
head(res[,c(2,5:8)],n=50)

#gse_kegg

kk3 <- gseKEGG(geneList     = gene_list_3,
               organism     = 'hsa',
               minGSSize    = 10,
               pvalueCutoff = 1,
               verbose      = FALSE)
head(as.data.frame(kk3)[,c(2,7,8)],n=50)

##############################################################################################################################
#Reactome patheway enrichment
library(ReactomePA)

gse.rp = gsePathway(geneList     = gene_list_3,
                    organism     = 'human',
                    nPerm=1000,
                    minGSSize= 50,
                    maxGSSize= 500,
                    pvalueCutoff = 1,
                    verbose      = FALSE)
head(as.data.frame(gse.rp)[,c(2,7,8)],n=50)


# gene_list_3
# de <- names(gene_list_3)[abs(gene_list_3) > 1.5]
# kk4 <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
# head(kk4)

##############################################################################################################################
# Save all the results

saveRDS(gse, file = "gse_severe_negative.rds")
saveRDS(gse2, file = "gse_mild_negative.rds")
saveRDS(gse3, file = "gse_severe_mild.rds")
saveRDS(gsea.hm, file = "hallmark_severe_mild.rds")
saveRDS(kk3, file = "KEGG_severe_mild.rds")
saveRDS(kk4, file = "Reactome_severe_mild.rds")

