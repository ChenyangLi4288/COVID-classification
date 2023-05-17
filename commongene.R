
library(dplyr)
library(ReactomePA)

setwd("/Users/chenyangli/Downloads/Coviddata")
gc_columbia <- readRDS("gc3.rds")
broad_count <- readRDS("/Users/chenyangli/Downloads/Coviddata/edata/broad/count_scRNA2bulk_broad_QC.rds")
broad_meta <- readRDS("/Users/chenyangli/Downloads/Coviddata/edata/broad/meta_scRNA2bulk_broad.rds")

anno = readRDS("anno_grch38_v2.rds")
gc = intersect(rownames(broad_count), anno$hgnc_symbol[anno$transcript_biotype=="protein_coding"])
length(gc)
broad_count = broad_count[match(gc, rownames(broad_count)),]
dim(broad_count)

## Filter out the genes having no expression
broad_count = broad_count[rowSums(broad_count)>0,]
dim(broad_count)

broad_count <- log2(broad_count + 1)


# Identify donor_ids for each group
COVID_negative_ids <- broad_meta$donor_id[broad_meta$SARSCoV2_PCR_Status == "neg"]
COVID_severe_ids <- broad_meta$donor_id[broad_meta$SARSCoV2_PCR_Status == "pos" & broad_meta$Cohort_Disease_WHO_Score == "COVID19_WHO_6-8"]
COVID_mild_ids <- broad_meta$donor_id[broad_meta$SARSCoV2_PCR_Status == "pos" & broad_meta$Cohort_Disease_WHO_Score == "COVID19_WHO_1-5"]

# Subset broad_count dataframe for each group
COVID_negative_counts <- broad_count[, COVID_negative_ids]
COVID_severe_counts <- broad_count[, COVID_severe_ids]
COVID_mild_counts <- broad_count[, COVID_mild_ids]


tpm_mild <- as.matrix(COVID_mild_counts)
tpm_severe <- as.matrix(COVID_severe_counts)
tpm_negative <- as.matrix(COVID_negative_counts)

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

tpm_negative_sd[tpm_negative_sd<0.01]=0
tpm_mild_sd[tpm_mild_sd<0.01]=0
tpm_severe_sd[tpm_severe_sd<0.01]=0

head(tpm_severe_mean)
head(tpm_negative_mean)
head(tpm_mild_mean)
#========================
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
print(gc3)

severe_negative = severe_negative[,gc3]
mild_negative = mild_negative[,gc3]
severe_mild = severe_mild[,gc3]

head(severe_negative)
head(mild_negative)
head(severe_mild)

