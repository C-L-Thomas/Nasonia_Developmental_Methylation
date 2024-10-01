library(readr)
library(reshape2)
library(ggplot2)

#Read in dataset
Meth_data <- read_delim("~/Whole_Development_Methylation_Percentages.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

Meth_data

Meth_data$embryo <- rowMeans(subset(Meth_data, select = c(6,7)))
Meth_data$larva <- rowMeans(subset(Meth_data, select = c(8,9,10)))
Meth_data$prepupae <- rowMeans(subset(Meth_data, select = c(11,12,13)))
Meth_data$pupae <- rowMeans(subset(Meth_data, select = c(14,15,16)))
Meth_data$adult <- rowMeans(subset(Meth_data, select = c(17,18,19)))

cols <- c(5,20,21,22,23,24)
Methylated_Genes <- Meth_data[,c(cols)]
Bi_data <- melt(Methylated_Genes, value.name = "Percentage")
colnames(Bi_data) <- c("geneID", "DStage", "Percentage")
#nrow(Bi_data)

#Read in transcriptome results
transcriptome <- read.csv("~/Normalised_Gene_Counts.csv")
transcriptome$embryo <- rowMeans(subset(transcriptome, select = c(2,3,4)))
transcriptome$larva <- rowMeans(subset(transcriptome, select = c(5,6,7)))
transcriptome$prepupae <- rowMeans(subset(transcriptome, select = c(8,9,10)))
transcriptome$pupae <- rowMeans(subset(transcriptome, select = c(11,12,13)))
transcriptome$adult <- rowMeans(subset(transcriptome, select = c(14,15,16)))

cols <- c(1,17,18,19,20,21)
RNA <- transcriptome[,c(cols)]

RNA_data <- melt(RNA, value.name = "Counts")
colnames(RNA_data) <- c("geneID", "DStage", "Counts")
#Remaining<- RNA_data[(RNA_data$geneID %in% Bi_data$geneID),]
#nrow(Remaining)

meth_exp <- merge(Bi_data, RNA_data, by=c("geneID","DStage"), all=T)
meth_exp <- meth_exp[!is.na(meth_exp$Percentage & !is.na(meth_exp$Counts)),]
meth_exp$logCounts <- log10(meth_exp$Counts)
meth_exp$DStage <- as.factor(meth_exp$DStage)
results <- na.omit(meth_exp)

ggplot(meth_exp, aes(x=Percentage, y=logCounts, colour = DStage))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Percentage Methylation")+
  ylab("log(Counts)")

####################################################
library(emmeans)

hist(results$logCounts) #is normal
model <- glm(Counts~Percentage*DStage, family=poisson,data=results, method = "glm.fit")

library(car)
#Report the interaction stats from Anova. 
Anova(model)
#summary (model)
model$deviance/model$df.residual #Overdispersion > 1.2 bad

#pseudo R2
100*((model$null.deviance-model$deviance)/model$null.deviance) #pseudo R2

#The slopes (estimates) and significant stats could be reported here
emtrends(model, pairwise~DStage, var = "Percentage")
emmip(model, DStage ~ Percentage, cov.reduce = range)

ggplot(results, aes(x = Percentage, y = log(Counts), colour = DStage)) +
  facet_wrap(~DStage, nrow=3) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(results, pred = predict(model)), aes(y = log(pred))) +
  theme(legend.position = "none")
