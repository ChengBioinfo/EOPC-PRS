
# Run for 269-, General-, EO-, LO, merged-General-, merged-EO-, merged-LO-PRS

library(dplyr)
library(survival)
library(tibble)
library(pROC)
library(ggplot2)
library(survivalROC)
library(ggthemes)

PRS_matrix <- readRDS("PRS_matrix.Rds")

#survivalROC ---------------------------------------------------------

ROC55 <- survivalROC(Stime = PRS_matrix$ExitAge,
                     status = PRS_matrix$Phe, 
                     marker = PRS_matrix$PRS,   
                     predict.time = 55,
                     method = 'KM')

PRS_matrix_de55 <- PRS_matrix[-which(PRS_matrix$ExitAge <= 55),]
ROC70 <- survivalROC(Stime = PRS_matrix_de55$ExitAge,
                     status = PRS_matrix_de55$Phe, 
                     marker = PRS_matrix_de55$PRS,   
                     predict.time = 70,
                     method = 'KM')

PRS_matrix_de70 <- PRS_matrix[-which(PRS_matrix$ExitAge <= 70),]
ROC85 <- survivalROC(Stime = PRS_matrix_de70$ExitAge,
                     status = PRS_matrix_de70$Phe, 
                     marker = PRS_matrix_de70$PRS,   
                     predict.time = 85,
                     method = 'KM')

ROC <- list(ROC55, ROC70, ROC85)

ROC_for_plot_55 <- data.frame(FP = ROC[[1]]$FP, TP = ROC[[1]]$TP, Age = "55")
ROC_for_plot_70 <- data.frame(FP = ROC[[2]]$FP, TP = ROC[[2]]$TP, Age = "70")
ROC_for_plot_85 <- data.frame(FP = ROC[[3]]$FP, TP = ROC[[3]]$TP, Age = "85")
ROC_for_plot <- rbind(ROC_for_plot_55, ROC_for_plot_70, ROC_for_plot_85)
#ROC_for_plot$Population <- factor(ROC_for_plot$Population, levels = c("All Subjects", "Positive Family History", "Negative Family History"))

pdf("survivalROC_GeneralPRS.pdf")
ggplot(ROC_for_plot, aes(x = FP, y = TP, group = Age)) +
  geom_line(aes(color = Age), size = 1.2) +
  scale_color_manual(values=c("#FCBC29", "#36907E", "#C24328"),
                     name = NULL, 
                     labels = c(paste0("55-year old\nAUC = ", round(ROC[[1]]$AUC,3)),
                                paste0("(56~70)-year old\nAUC = ", round(ROC[[2]]$AUC,3)),
                                paste0("(71~85)-year old\nAUC = ", round(ROC[[3]]$AUC,3)))) +
  ggtitle("ROC of General-PRS for PCa") +
  labs(x = "False Positive Value", y = "True Positive Value") +
  theme(axis.title.x = element_text(vjust = -1, size = 28),
        axis.title.y = element_text(vjust = 2, size = 28),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = 0, size = 24),
        axis.text.y = element_text(vjust = 0.5, size = 24),
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 28),
        panel.border = element_rect(fill = NA,color = "black", size = 1, linetype = "solid"),
        panel.grid.major = element_line(color = NA, size = .25),
        panel.grid.minor = element_line(color = NA, size = .125),
        panel.background = element_rect(fill = NA, color = NA),
        legend.position = c(0.92, 0.16),
        legend.justification = c(0.9, 0.4),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1.0,"cm")) +
  geom_abline(intercept=seq(0, 1, 1), slope=1, colour="black")
dev.off()