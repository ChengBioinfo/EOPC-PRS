library(gwasurvivr)

# Run for each chromosome of General-, EO- and LO-population

cov <- readRDS("phenotype_PCa_167517.Rds")
sample.ids <- as.character(cov$ID)
plinkCoxSurv(bed.file="UKBB_QC/General/chr1.bed",
             covariate.file=cov,
             sample.ids=sample.ids,
             id.column="ID",
             covariates=c("age","Assessment_centre","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
             time.to.event="Prostate_cancer_time",
             event="Prostate_cancer_cancer_risk",
             print.covs="only",
             out.file="chr1_cox",
             maf.filter=0.01,
             chunk.size = 5000,
             flip.dosage=FALSE,
             verbose=TRUE)