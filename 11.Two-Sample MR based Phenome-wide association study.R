library(data.table)
library(dplyr)
library(TwoSampleMR)

#获取工具变量SNPs-------------------------------------

load("available_outcomes.RData")

#format data ----------------------------------

data<-readRDS("cox_all_EOPC.Rds")  # generated from Genome-wide Cox regression analysis
head(data)

data$beta<-log(data$HR)
data$se<-data$beta/abs(qt(data$PVALUE/2, df = data$N))

outcome_dat <- format_data(
  dat = data,
  type = "outcome",
  # snps = NULL,
  # header = TRUE,
  # phenotype_col = "Phenotype",
  snp_col = "RSID",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "SAMP_MAF",
  effect_allele_col = "A0",
  other_allele_col = "A1",
  pval_col = "PVALUE",
  # units_col = "units",
  # ncase_col = "ncase",
  # ncontrol_col = "ncontrol",
  samplesize_col = "N",
  # gene_col = "gene",
  # id_col = "id",
  min_pval = 1e-200,
  # z_col = "z",
  # info_col = "info",
  chr_col = "CHR",
  pos_col = "POS",
  log_pval = FALSE
)

head(outcome_dat)

saveRDS(outcome_dat, "outcome_dat_EOPC_IV.Rds")

######   loop for MR as exposure and lung cancer survival as outcome

restBind<-as.data.frame(matrix(nrow = 0, ncol = 13))
colnames(restBind)<-c("id.exposure", "id.outcome","outcome", "exposure", "method", "nsnp", "b", "se",
                      "pval", "Q", "Q_df", "Q_pval", "egger_intercept")

restFix<-as.data.frame(matrix(nrow = 1, ncol = 4))
colnames(restFix)<-c("Q", "Q_df", "Q_pval", "egger_intercept")

which(grepl(id, available_outcomes$id))

############# find row count ###############
for (id in available_outcomes$id){
  test2<-(try(exposure_dat<-extract_instruments(outcomes = id,
                                                p1 = 5e-08,
                                                clump = TRUE,
                                                p2 = 5e-08,
                                                r2 = 0.001,
                                                kb = 10000,
                                                access_token = ieugwasr::check_access_token(),
                                                force_server = FALSE)))
  if(class(test2) == "NULL") {
    restFix2<-as.data.frame(matrix(nrow = 1, ncol = 13))
    colnames(restFix2)<-c("id.exposure", "id.outcome","outcome", "exposure", "method", "nsnp", "b", "se",
                          "pval", "Q", "Q_df", "Q_pval", "egger_intercept")
    restFix2[1,1]<-id
    restFix2
    restBind<-rbind(restBind, restFix2)
  }
  else{
    exposure_dat<-extract_instruments(outcomes = id,
                                      p1 = 5e-08,
                                      clump = TRUE,
                                      p2 = 5e-08,
                                      r2 = 0.001,
                                      kb = 10000,
                                      access_token = ieugwasr::check_access_token(),
                                      force_server = FALSE)
    head(exposure_dat)
    test1<-(try(dat <- harmonise_data(exposure_dat = exposure_dat,
                                      outcome_dat = outcome_dat,
                                      action = 2)))
    if (class(test1) == "try-error") {
      restFix1<-as.data.frame(matrix(nrow = 1, ncol = 13))
      colnames(restFix1)<-c("id.exposure", "id.outcome","outcome", "exposure", "method", "nsnp", "b", "se",
                            "pval", "Q", "Q_df", "Q_pval", "egger_intercept")
      restFix1[1,1]<-id
      restFix1
      restBind<-rbind(restBind, restFix1)
    }
    else{
      dat <- harmonise_data(exposure_dat = exposure_dat,
                            outcome_dat = outcome_dat,
                            action = 2)
      # MR
      test3<-(try(res <- mr(dat)))
      if (dim(test3)[1] == 0) {
        restFix3<-as.data.frame(matrix(nrow = 1, ncol = 13))
        colnames(restFix3)<-c("id.exposure", "id.outcome","outcome", "exposure", "method", "nsnp", "b", "se",
                              "pval", "Q", "Q_df", "Q_pval", "egger_intercept")
        restFix3[1,1]<-id
        restFix3
        restBind<-rbind(restBind, restFix3)
      }
      else{
        res <- mr(dat)
        het <- mr_heterogeneity(dat)
        ple <- mr_pleiotropy_test(dat)
        restBindsub<-bind_rows(res, het, ple)
        # test<-(try(restBind <- rbind(restBind, restBindsub)))
        # if (class(test) == "try-error") {
        if (dim(restBindsub)[2] != 13){
          restBindsub<-cbind(restBindsub, restFix)
          restBind<-rbind(restBind, restBindsub)
        }
        else{
          restBind<-rbind(restBind, restBindsub)
        }
      }
    }
  }
}