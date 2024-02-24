
# Run for 269-, General-, EO-, LO, merged-General-, merged-EO-, merged-LO-PRS

library(meta)

data <- read.table("metaanalysis.txt", header=T, sep = "\t")  # generated from the WCoxPH or logistic results of perSD

str(data)

data$lgOR <- log(data$OR)
data$selgOR <- (log(data$UCI) - log(data$LCI)) / (2*1.96)

meta <- metagen(TE = lgOR, seTE = selgOR, studlab = study, data = data, sm = "OR", subgroup = Group, fixed = F)

pdf("meta.pdf")
forest(meta, overall = F, text.random = "", text.w.random = "", print.subgroup.name = F, lwd = 3, hetlab = "",
       test.subgroup = F, col.study = "#926D33", print.tau2 = F, leftcols = c("studlab"))
dev.off()

het <- data.frame(I2 = meta$I2.w, P = meta$pval.Q.w)
write.table(het, "heter.txt", row.names = F, quote = F, sep = "\t")
