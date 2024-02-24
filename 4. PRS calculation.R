
# Run for 269-, General-, EO-, LO, merged-General-, merged-EO-, merged-LO-PRS

library(dplyr)
library(survival)
library(tibble)

# matrix processing --------------------------------------------------------------------------

# extract genotyping of corresponding SNPs from UKB general population

geno <- read.table("merged.raw", sep=" ", header=T)

phe <- readRDS("phenotype_PCa_167517.Rds")
rownames(phe) <- NULL
phe <- column_to_rownames(phe, "ID")

rownames(geno) <- geno[,1]
geno <- geno[, -c(1:6)]
library(stringr)
colnames(geno) <- str_replace(colnames(geno), "_(?=[^_]+$)", "-")
colnames(geno) <- gsub("-.*$", "", colnames(geno))
colnames(geno) <- gsub("X", "", colnames(geno))
colnames(geno) <- gsub("\\.", ":", colnames(geno))

#merge
matrix <- merge(phe, geno, by = 0)
matrix <- column_to_rownames(matrix, "Row.names")

#PRS matrix-------------------------------------------------------------------------------

# read HR
HR <- read.table("cox_result.txt", sep = "\t", header = F)
HR <- HR[, c(1,9)]
colnames(HR) <- c("SNP", "HR")

HR <- HR[HR$SNP %in% colnames(geno),]
HR$coe <- log(HR$HR)
rownames(HR) <- HR$SNP
HR <- HR[,-1]

#general population
tmatrix <- as.data.frame(t(matrix))
tmatrix <- tmatrix[-c(1:18),]
mergeAll <- merge(HR, tmatrix, by = 0)
mergeAll <- column_to_rownames(mergeAll, "Row.names")
mergeAll <- mergeAll[, -1]
mergeAll <- as.data.frame(t(mergeAll))

PRScal_noNA <- mergeAll
PRScal_noNA[is.na(PRScal_noNA)] <- 0

rm(matrix, tmatrix, geno)
gc()

b <- as.numeric(PRScal_noNA[1,])
PRScal_noNA <- PRScal_noNA[-1,]

PRS_indi <- c()
for (i in 1:ncol(PRScal_noNA)) {
  PRS_tmp <- b[i] * PRScal_noNA[,i]
  PRS_indi <- cbind(PRS_indi, PRS_tmp)
}
PRS <- apply(PRS_indi, 1, sum)

mergeAll <- mergeAll[-1,]
PRS <- cbind(mergeAll, PRS)

PRS_matrix <- merge(phe, PRS, by = 0)
rownames(PRS_matrix) <- PRS_matrix[,1]
colnames(PRS_matrix)
PRS_matrix <- PRS_matrix[, -c(1, 22:(ncol(PRS_matrix)-1))]

saveRDS(PRS_matrix, "PRS_matrix.Rds")