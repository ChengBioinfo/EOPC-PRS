
data_idw <- read.table("White_PCa_incidence_byage.txt", sep = "", header = T) #download from USCS
colnames(data_idw)[4] <- "IR" 
data_idw$Age_Group <- c(32,37,42,47,52,57,62,67,72,77,82,87)

head(data_idw)

plot(data_idw$Age_Group, data_idw$IR, col='blue', pch=19)

model_idw <- lm(data_idw$IR ~ data_idw$Age_Group, data = data_idw)

data_idw_new <- approx(data_idw$Age_Group, data_idw$IR, xout=c(30:31,33:36,38:41,43:46,48:51,53:56,58:61,63:66,68:71,73:76,78:81,83:86))
points(data_idw_new$x,data_idw_new$y,col='red',pch=19)

Incidence_rate_impute <- data.frame(Age = data_idw_new$x, Incidence = data_idw_new$y)
Incidence_rate_oringin <- data_frame(Age = data_idw$Age_Group, Incidence = data_idw$IR)
Incidence_rate <- rbind(Incidence_rate_impute, Incidence_rate_oringin)
Incidence_rate <- Incidence_rate[order(Incidence_rate$Age),]
Incidence_rate <- Incidence_rate[-c(1,2),]

write.table(Incidence_rate, file = "White_PCa_incidence_allage.txt", row.names = F, quote = F, sep = "\t")

library(ggplot2)
pdf("dotplot_White_PCa_incidence_allage.pdf", width = 8, height = 6)
theme_set(theme_bw())
ggplot(Incidence_rate, aes(x=Age, y=Incidence)) + 
  geom_line(color = "#926D33", size = 1.5) +
  labs(x = "Age", y = "Incidence (per 100,000)") +
  ylim(c(0,800)) +
  scale_x_continuous(limits=c(30,90), breaks=seq(30,90,by=10)) + 
  theme(axis.title.x = element_text(vjust = 0, size = 14),
        axis.title.y = element_text(vjust = 2, size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(margin = margin(10, 0, 10, 0), size = 18)) +
  ggtitle("Incidence of PCa")
dev.off()
