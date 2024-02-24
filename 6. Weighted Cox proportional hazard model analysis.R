# Run for 269-, General-, EO-, LO, merged-General-, merged-EO-, merged-LO-PRS

PRS_matrix <- readRDS("PRS_matrix.Rds")

#——level-----------------------------------------------------------------------------

PRScore <- PRS_matrix$PRS
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

#weight-------------------------------------------------------------------------------

incidence <- read.table("White_incidence_byage/White_PCa_incidence_allage.txt", header = T, sep = "\t")
incidence$weight <- 1/incidence$Incidence*100000

PRS_matrix$weight <- 0
PRS_matrix$weight[PRS_matrix$Prostate_cancer_cancer_risk == 1] <- 1

for (i in c(1:nrow(incidence))) {
  PRS_matrix$weight[PRS_matrix$Prostate_cancer_cancer_risk == 0 & PRS_matrix$exitage == incidence[i,1]] <- incidence[i,3]
}

EOPC <- PRS_matrix[which(PRS_matrix$exitage <= 55),]

LOPC <- PRS_matrix[which(PRS_matrix$exitage > 55),]

write.table(cbind(rownames(PRS_matrix), PRS_matrix), "PRS_matrix.txt", quote = F, sep = "\t", row.names = F)

# WCoxPH, general population-----------------------------------------------------------------------

PRS_matrix$zscore <- scale(PRS_matrix$PRS, center = T, scale = T)
saveRDS(PRS_matrix, "PRS_matrix.Rds")

library(survey)

design <- svydesign(id      = ~0,
                    weights = ~weight,
                    data    = PRS_matrix)

cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ level1 + level2 + level3 + level4 +
                  age + Assessment_centre + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=PRS_matrix)
HR <- summary(cox)$coefficients[1:4,2]
CI <- summary(cox)$conf.int[1:4,c(3,4)]
P <- summary(cox)$coefficients[1:4,6]
res1 <- cbind(HR,CI,P)

cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ zscore + 
                  age + Assessment_centre + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=PRS_matrix)
HR <- summary(cox)$coefficients[1,2]
CI <- summary(cox)$conf.int[1,c(3,4)]
P <- summary(cox)$coefficients[1,6]
res2 <- c(HR,CI,P)

coxresult <- rbind(res1, res2)
rownames(coxresult)[5] <- "perSD"

write.table(cbind(rownames(coxresult), coxresult), "general_coxresults.txt", sep="\t", quote = F, row.names = F)

# WCoxPH, EOPC-------------------------------------------------------------------------------

table(EOPC$level, EOPC$Prostate_cancer_cancer_risk)
table(EOPC$level)

EOPC$zscore <- scale(EOPC$PRS, center = T, scale = T)

saveRDS(EOPC, "PRS_EOPC_matrix.Rds")

design <- svydesign(id      = ~0,
                    weights = ~weight,
                    data    = EOPC)

cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ level1 + level2 + level3 + level4 +
                  age + Assessment_centre + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=EOPC)
HR <- summary(cox)$coefficients[1:4,2]
CI <- summary(cox)$conf.int[1:4,c(3,4)]
P <- summary(cox)$coefficients[1:4,6]
res1 <- cbind(HR,CI,P)

cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ zscore + 
                  age + Assessment_centre + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=EOPC)
HR <- summary(cox)$coefficients[1,2]
CI <- summary(cox)$conf.int[1,c(3,4)]
P <- summary(cox)$coefficients[1,6]
res2 <- c(HR,CI,P)

coxresult <- rbind(res1, res2)
rownames(coxresult)[5] <- "perSD"

write.table(cbind(rownames(coxresult), coxresult), "EOPC_coxresults.txt", sep="\t", quote = F, row.names = F)

# WCoxPH, LOPC-------------------------------------------------------------------------------

table(LOPC$level)
table(LOPC$level, LOPC$Prostate_cancer_cancer_risk)

LOPC$zscore <- scale(LOPC$PRS, center = T, scale = T)

saveRDS(LOPC, "PRS_LOPC_matrix.Rds")

design <- svydesign(id      = ~0,
                    weights = ~weight,
                    data    = LOPC)

cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ level1 + level2 + level3 + level4 + 
                  age + Assessment_centre + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=LOPC)
HR <- summary(cox)$coefficients[1:4,2]
CI <- summary(cox)$conf.int[1:4,c(3,4)]
P <- summary(cox)$coefficients[1:4,6]
res1 <- cbind(HR,CI,P)

cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ zscore + 
                  age + Assessment_centre + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=LOPC)
HR <- summary(cox)$coefficients[1,2]
CI <- summary(cox)$conf.int[1,c(3,4)]
P <- summary(cox)$coefficients[1,6]
res2 <- c(HR,CI,P)

coxresult <- rbind(res1, res2)
rownames(coxresult)[5] <- "perSD"

write.table(cbind(rownames(coxresult), coxresult), "LOPC_coxresults.txt", sep="\t", quote = F, row.names = F)

#@@****fh startification*****@@-------------------------------------------------------------------------------

# WCoxPH, all -----------------------------------------------------------------------

withfh <- PRS_matrix[which(PRS_matrix$fh == 1),]
withoutfh <- PRS_matrix[which(PRS_matrix$fh == 0),]
fh <- list(withfh, withoutfh)

result <- c()

for (j in 1:length(fh)) {
  
  fh[[j]]$zscore <- scale(fh[[j]]$PRS, center = T, scale = T)
  
  design <- svydesign(id      = ~0,
                      weights = ~weight,
                      data    = fh[[j]])
  
  cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ level1 + level2 + level3 + level4 +
                    age + Assessment_centre + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=fh[[j]])
  HR <- summary(cox)$coefficients[1:4,2]
  CI <- summary(cox)$conf.int[1:4,c(3,4)]
  P <- summary(cox)$coefficients[1:4,6]
  res1 <- cbind(HR,CI,P)
  
  cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ zscore + 
                    age + Assessment_centre + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=fh[[j]])
  HR <- summary(cox)$coefficients[1,2]
  CI <- summary(cox)$conf.int[1,c(3,4)]
  P <- summary(cox)$coefficients[1,6]
  res2 <- c(HR,CI,P)
  
  coxresult <- rbind(res1, res2)
  rownames(coxresult)[5] <- "perSD"
  
  result <- cbind(result, coxresult)
  
}

write.table(cbind(rownames(result),result), "general_coxresults_fh_stra.txt",sep="\t",quote = F,row.names = F)

#fh_EOPC--------------------------------------------------------------------------------

withfh <- EOPC[which(EOPC$fh == 1),]
withoutfh <- EOPC[which(EOPC$fh == 0),]
fh <- list(withfh, withoutfh)

result <- c()

for (j in 1:length(fh)) {
  
  fh[[j]]$zscore <- scale(fh[[j]]$PRS, center = T, scale = T)
  
  design <- svydesign(id      = ~0,
                      weights = ~weight,
                      data    = fh[[j]])
  
  cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ level1 + level2 + level3 + level4 +
                    age + Assessment_centre + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=fh[[j]])
  HR <- summary(cox)$coefficients[1:4,2]
  CI <- summary(cox)$conf.int[1:4,c(3,4)]
  P <- summary(cox)$coefficients[1:4,6]
  res1 <- cbind(HR,CI,P)
  
  cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ zscore + 
                    age + Assessment_centre + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=fh[[j]])
  HR <- summary(cox)$coefficients[1,2]
  CI <- summary(cox)$conf.int[1,c(3,4)]
  P <- summary(cox)$coefficients[1,6]
  res2 <- c(HR,CI,P)
  
  coxresult <- rbind(res1, res2)
  rownames(coxresult)[5] <- "perSD"
  
  result <- cbind(result, coxresult)
  
}

write.table(cbind(rownames(result),result), "EOPC_coxresults_fh_stra.txt",sep="\t",quote = F,row.names = F)

#fh_LOPC--------------------------------------------------------------------------------

withfh <- LOPC[which(LOPC$fh == 1),]
withoutfh <- LOPC[which(LOPC$fh == 0),]
fh <- list(withfh, withoutfh)

result <- c()

for (j in 1:length(fh)) {
  
  fh[[j]]$zscore <- scale(fh[[j]]$PRS, center = T, scale = T)
  
  design <- svydesign(id      = ~0,
                      weights = ~weight,
                      data    = fh[[j]])
  
  cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ level1 + level2 + level3 + level4 +
                    age + Assessment_centre + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=fh[[j]])
  HR <- summary(cox)$coefficients[1:4,2]
  CI <- summary(cox)$conf.int[1:4,c(3,4)]
  P <- summary(cox)$coefficients[1:4,6]
  res1 <- cbind(HR,CI,P)
  
  cox <- svycoxph(Surv(Prostate_cancer_time, Prostate_cancer_cancer_risk) ~ zscore + 
                    age + Assessment_centre + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, design = design, data=fh[[j]])
  HR <- summary(cox)$coefficients[1,2]
  CI <- summary(cox)$conf.int[1,c(3,4)]
  P <- summary(cox)$coefficients[1,6]
  res2 <- c(HR,CI,P)
  
  coxresult <- rbind(res1, res2)
  rownames(coxresult)[5] <- "perSD"
  
  result <- cbind(result, coxresult)
  
}

write.table(cbind(rownames(result),result), "LOPC_coxresults_fh_stra.txt",sep="\t",quote = F,row.names = F)
