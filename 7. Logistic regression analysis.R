
# Run for 269-, General-, EO-, LO, merged-General-, merged-EO-, merged-LO-PRS

library(dplyr)
library(survival)
library(tibble)

#处理矩阵--------------------------------------------------------------------------

geno <- read.table("merged.raw", sep=" ", header=T)
num <- grep("_dup", colnames(geno))
geno <- geno[, -num]

phe <- readRDS("phenotype_PCa_175349.Rds")
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
PRS_matrix <- PRS_matrix[, -c(1, 20:(ncol(PRS_matrix)-1))]

#——level-----------------------------------------------------------------------------

PRScore <- PRS_matrix[PRS_matrix$Prostate_cancer_cancer_risk == 0,]$PRS
qui <- quantile(PRScore,probs = seq(0,1,0.2), name = FALSE)
first <- qui[2]
second <- qui[3]
third <- qui[4]
forth <- qui[5]

PRS_matrix$level <- 0

PRS_matrix$level[which(PRS_matrix$PRS < first)] <- 0
PRS_matrix$level[which(PRS_matrix$PRS >= first & PRS_matrix$PRS < second)] <- 1
PRS_matrix$level[which(PRS_matrix$PRS >= second & PRS_matrix$PRS < third)] <- 2
PRS_matrix$level[which(PRS_matrix$PRS >= third & PRS_matrix$PRS < forth)] <- 3
PRS_matrix$level[which(PRS_matrix$PRS >= forth)] <- 4

table(PRS_matrix$level)

#——dummy variables----------------------------------------------------------------------------
PRS_matrix$level1 <- 0
PRS_matrix$level2 <- 0
PRS_matrix$level3 <- 0
PRS_matrix$level4 <- 0

PRS_matrix$level1[which(PRS_matrix$level == 1)] <- 1
PRS_matrix$level2[which(PRS_matrix$level == 2)] <- 1
PRS_matrix$level3[which(PRS_matrix$level == 3)] <- 1
PRS_matrix$level4[which(PRS_matrix$level == 4)] <- 1

table(PRS_matrix$level, PRS_matrix$Prostate_cancer_cancer_risk)

EOPC <- PRS_matrix[which(PRS_matrix$exitage <= 55),]

LOPC <- PRS_matrix[which(PRS_matrix$exitage > 55),]

write.table(cbind(rownames(PRS_matrix), PRS_matrix), "PRS_matrix.txt", quote = F, sep = "\t", row.names = F)

# logit regression, general population-----------------------------------------------------------------------

PRS_matrix$zscore <- scale(PRS_matrix$PRS, center = T, scale = T)

logit <- glm(Prostate_cancer_cancer_risk ~ level1 + level2 + level3 + level4 +
               exitage + Assessment_centre + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             family = binomial(link = "logit"), data=PRS_matrix)
OR <- exp(summary(logit)$coefficients)[2:5,1]
LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2:5]
HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2:5]
P <- summary(logit)$coefficients[2:5,4]

res1 <- cbind(OR,LCI,HCI,P)

logit <- glm(Prostate_cancer_cancer_risk ~ zscore +
               exitage + Assessment_centre + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             family = binomial(link = "logit"), data=PRS_matrix)
OR <- exp(summary(logit)$coefficients)[2,1]
LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
P <- summary(logit)$coefficients[2,4]

res2 <- cbind(OR,LCI,HCI,P)

logit <- rbind(res2,res1)

write.table(cbind(rownames(logit),logit), "general_logit_results.txt", sep="\t", quote = F, row.names = F)

# logit regression, EOPC-------------------------------------------------------------------------------

table(EOPC$level, EOPC$Prostate_cancer_cancer_risk)
table(EOPC$level)

EOPC$zscore <- scale(EOPC$PRS, center = T, scale = T)

saveRDS(EOPC, "PRS_EOPC_matrix.Rds")
write.table(cbind(rownames(EOPC), EOPC), "PRS_EOPC_matrix.txt", quote = F, sep = "\t", row.names = F)

logit <- glm(Prostate_cancer_cancer_risk ~ level1 + level2 + level3 + level4 +
               exitage + Assessment_centre + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             family = binomial(link = "logit"), data=EOPC)
OR <- exp(summary(logit)$coefficients)[2:5,1]
LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2:5]
HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2:5]
P <- summary(logit)$coefficients[2:5,4]

res1 <- cbind(OR,LCI,HCI,P)

logit <- glm(Prostate_cancer_cancer_risk ~ zscore +
               exitage + Assessment_centre + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             family = binomial(link = "logit"), data=EOPC)
OR <- exp(summary(logit)$coefficients)[2,1]
LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
P <- summary(logit)$coefficients[2,4]

res2 <- cbind(OR,LCI,HCI,P)

logit <- rbind(res2,res1)

write.table(cbind(rownames(logit),logit), "EOPC_logit_results.txt", sep="\t", quote = F, row.names = F)

# logit regression, LOPC-------------------------------------------------------------------------------

table(LOPC$level)
table(LOPC$level, LOPC$Prostate_cancer_cancer_risk)

LOPC$zscore <- scale(LOPC$PRS, center = T, scale = T)

saveRDS(LOPC, "PRS_LOPC_matrix_10Q.Rds")
write.table(cbind(rownames(LOPC), LOPC), "PRS_LOPC_matrix_10Q.txt", quote = F, sep = "\t", row.names = F)

logit <- glm(Prostate_cancer_cancer_risk ~ level1 + level2 + level3 + level4 +
               exitage + Assessment_centre + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             family = binomial(link = "logit"), data=LOPC)
OR <- exp(summary(logit)$coefficients)[2:5,1]
LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2:5]
HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2:5]
P <- summary(logit)$coefficients[2:5,4]

res1 <- cbind(OR,LCI,HCI,P)

logit <- glm(Prostate_cancer_cancer_risk ~ zscore +
               exitage + Assessment_centre + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             family = binomial(link = "logit"), data=LOPC)
OR <- exp(summary(logit)$coefficients)[2,1]
LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
P <- summary(logit)$coefficients[2,4]

res2 <- cbind(OR,LCI,HCI,P)

logit <- rbind(res2,res1)

write.table(cbind(rownames(logit),logit), "LOPC_logit_results.txt", sep="\t", quote = F, row.names = F)

#@@****fh stratification*****@@-------------------------------------------------------------------------------

# logit regression, all -----------------------------------------------------------------------

withfh <- PRS_matrix[which(PRS_matrix$fh == 1),]
withoutfh <- PRS_matrix[which(PRS_matrix$fh == 0),]
fh <- list(withfh, withoutfh)

logitall <- c()

for (j in 1:length(fh)) {
  
  fh[[j]]$zscore <- scale(fh[[j]]$PRS, center = T, scale = T)
  
  logit <- glm(Prostate_cancer_cancer_risk ~ level1 + level2 + level3 + level4 +
                 exitage + Assessment_centre + 
                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
               family = binomial(link = "logit"), data=fh[[j]])
  OR <- exp(summary(logit)$coefficients)[2:5,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2:5]
  HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2:5]
  P <- summary(logit)$coefficients[2:5,4]
  
  res1 <- cbind(OR,LCI,HCI,P)
  
  logit <- glm(Prostate_cancer_cancer_risk ~ zscore +
                 exitage + Assessment_centre + 
                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
               family = binomial(link = "logit"), data=fh[[j]])
  OR <- exp(summary(logit)$coefficients)[2,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
  HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
  P <- summary(logit)$coefficients[2,4]
  
  res2 <- cbind(OR,LCI,HCI,P)
  
  logit <- rbind(res2,res1)
  
  logitall <- cbind(logitall, logit)
  
}

write.table(cbind(rownames(logitall),logitall),"general_logit_results_fh_stra.txt", sep="\t", quote = F, row.names = F)

#fh_EOPC--------------------------------------------------------------------------------

withfh <- EOPC[which(EOPC$fh == 1),]
withoutfh <- EOPC[which(EOPC$fh == 0),]
fh <- list(withfh, withoutfh)

logitall <- c()

for (j in 1:length(fh)) {
  
  fh[[j]]$zscore <- scale(fh[[j]]$PRS, center = T, scale = T)
  
  logit <- glm(Prostate_cancer_cancer_risk ~ level1 + level2 + level3 + level4 +
                 exitage + Assessment_centre + 
                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
               family = binomial(link = "logit"), data=fh[[j]])
  OR <- exp(summary(logit)$coefficients)[2:5,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2:5]
  HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2:5]
  P <- summary(logit)$coefficients[2:5,4]
  
  res1 <- cbind(OR,LCI,HCI,P)
  
  logit <- glm(Prostate_cancer_cancer_risk ~ zscore +
                 exitage + Assessment_centre + 
                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
               family = binomial(link = "logit"), data=fh[[j]])
  OR <- exp(summary(logit)$coefficients)[2,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
  HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
  P <- summary(logit)$coefficients[2,4]
  
  res2 <- cbind(OR,LCI,HCI,P)
  
  logit <- rbind(res2,res1)
  
  logitall <- cbind(logitall, logit)
  
}

write.table(cbind(rownames(logitall),logitall),"EOPC_logit_results_fh_stra.txt", sep="\t", quote = F, row.names = F)

#fh_LOPC--------------------------------------------------------------------------------

withfh <- LOPC[which(LOPC$fh == 1),]
withoutfh <- LOPC[which(LOPC$fh == 0),]
fh <- list(withfh, withoutfh)

logitall <- c()

for (j in 1:length(fh)) {
  
  fh[[j]]$zscore <- scale(fh[[j]]$PRS, center = T, scale = T)
  
  logit <- glm(Prostate_cancer_cancer_risk ~ level1 + level2 + level3 + level4 +
                 exitage + Assessment_centre + 
                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
               family = binomial(link = "logit"), data=fh[[j]])
  OR <- exp(summary(logit)$coefficients)[2:5,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2:5]
  HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2:5]
  P <- summary(logit)$coefficients[2:5,4]
  
  res1 <- cbind(OR,LCI,HCI,P)
  
  logit <- glm(Prostate_cancer_cancer_risk ~ zscore +
                 exitage + Assessment_centre + 
                 PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
               family = binomial(link = "logit"), data=fh[[j]])
  OR <- exp(summary(logit)$coefficients)[2,1]
  LCI <- exp(summary(logit)$coefficients[,1]-1.96*summary(logit)$coefficients[,2])[2]
  HCI <- exp(summary(logit)$coefficients[,1]+1.96*summary(logit)$coefficients[,2])[2]
  P <- summary(logit)$coefficients[2,4]
  
  res2 <- cbind(OR,LCI,HCI,P)
  
  logit <- rbind(res2,res1)
  
  logitall <- cbind(logitall, logit)
  
}

write.table(cbind(rownames(logitall),logitall),"LOPC_logit_results_fh_stra.txt", sep="\t", quote = F, row.names = F)
