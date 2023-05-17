

#I've fixed the 1) and 2) you mentioned today, Friday I'll do the left part.


library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(organism, character.only = TRUE)

setwd("/Users/chenyangli/Downloads/Coviddata")
metdat <- read.table("metadata_cumc_hn2345.2.txt", header = TRUE)
count <- readRDS("count_cumc_hn2345.rds")
tpm <- readRDS("tpm_cumc_hn2345.rds")

tpm <- log2(tpm + 1)
# Label each sample with its diagnosis

metdat$diagnosis <- ifelse(metdat$viral.load == 'None', "Negative", ifelse(metdat$intubation == 1 | metdat$mortality == 1 | (metdat$days_hospitalized >= 7 & !is.na(metdat$days_hospitalized)), "Severe", "Mild"))
#Only keep the earlist sample
metadata_patient <- metdat %>%
  group_by(case.id) %>%
  filter(sample.id == min(sample.id))

metdat[1:6,]
# Get the unique diagnoses in metdat
unique_diagnoses <- unique(metdat$diagnosis)
tpm[,metdat$diagnosis == "Severe"]

# subset tpm by diagnosis
tpm_mild <- tpm[, metdat$diagnosis == "Mild"]
tpm_severe <- tpm[, metdat$diagnosis == "Severe"]
tpm_negative <- tpm[, metdat$diagnosis == "Negative"]

# calculate mean expression for each group
tpm_mild_mean <- apply(tpm_mild, 1, mean)
tpm_severe_mean <- apply(tpm_severe, 1, mean)
tpm_negative_mean <- apply(tpm_negative, 1, mean)
tpm_negative_sd <- apply(tpm_negative, 1, sd)

tpm_negative[1:5,1:10]
tpm_negative_mean[1:10]

# combine mean expression values into a new dataframe
tpm_by_diagnosis <- data.frame(tpm_mild_mean, tpm_severe_mean, tpm_negative_mean)


# Split the TPM data into a list of dataframes, one for each diagnosis
tpm_by_diagnosis <- lapply(unique_diagnoses, function(diag) {
  # Get the case IDs for this diagnosis
  case_ids <- metdat$case.id[metdat$diagnosis == diag]
  # Get the TPM data for these case IDs
  tpm_subset <- tpm[, colnames(tpm) %in% case_ids]
  # Return the subset of TPM data
  return(tpm_subset)
})

dim(tpm_mild)
dim(tpm_negative)
dim(tpm_severe)

#z-score
severe_negative = (tpm_severe_mean - t(tpm_negative_mean))/ (tpm_negative_sd + 0.00001)
mild_negative = (tpm_mild_mean - t(tpm_negative_mean)) / (tpm_negative_sd + 0.00001)

ori_gene_list <- severe_negative
print(ori_gene_list)
names(ori_gene_list) <- colnames(severe_negative)
gene_list_<-na.omit(ori_gene_list)
gene_list_<-sort(gene_list_, decreasing = TRUE)
gene_list_

ori_gene_list2 <- mild_negative
print(ori_gene_list2)
names(ori_gene_list2) <- colnames(mild_negative)
gene_list_2<-na.omit(ori_gene_list2)
gene_list_2<-sort(gene_list_2, decreasing = TRUE)
gene_list_2

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"


## Load annotation table
anno = readRDS("./anno_grch38_v2.rds")
head(anno)
## protein-coding genes:
gene.pc = intersect(names(gene_list_), anno$hgnc_symbol[anno$transcript_biotype=="protein_coding"])
gene.pc_2 = intersect(names(gene_list_2), anno$hgnc_symbol[anno$transcript_biotype=="protein_coding"])
## new gene list only for the protein coding genes
gene_list_ = gene_list_[match(gene.pc, names(gene_list_))]
gene_list_2 = gene_list_2[match(gene.pc, names(gene_list_2))]
## convert the gene symbol to entrez ids
names(gene_list_) = anno$entrezgene_id[match(names(gene_list_), anno$hgnc_symbol)]
names(gene_list_2) = anno$entrezgene_id[match(names(gene_list_2), anno$hgnc_symbol)]

length(gene_list_)
length(gene_list_2)

gene_list_2
gene_list_2<-sort(gene_list_2, decreasing = TRUE)


gse <- gseGO(geneList=gene_list_, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

head(gse)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


gse2 <- gseGO(geneList=gene_list_2, 
              ont ="ALL", 
              keyType = "ENTREZID", 
              nPerm = 10000, 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 1, 
              verbose = TRUE, 
              OrgDb = organism, 
              pAdjustMethod = "BH")

head(gse2)
dotplot(gse2, showCategory=10, split=".sign") + facet_grid(.~.sign)

saveRDS(gse, file = "gse_severe.rds")
saveRDS(gse2, file = "gse_mild.rds")
