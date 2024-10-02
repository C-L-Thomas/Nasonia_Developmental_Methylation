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
meth_exp$meth_category[(meth_exp$Percentage >= 0 & meth_exp$Percentage <= 0.5)] <- "none"
meth_exp$meth_category[(meth_exp$Percentage > 0.5 & meth_exp$Percentage <= 30)] <- "low"
meth_exp$meth_category[(meth_exp$Percentage > 30 & meth_exp$Percentage <= 70)] <- "medium"
meth_exp$meth_category[(meth_exp$Percentage > 70 & meth_exp$Percentage <= 100)] <- "high"

v_factor_levels <- c("none", "low", "medium","high")

table(meth_exp$meth_category)

ggplot(data = meth_exp, aes(x = DStage, fill = factor(meth_category, levels = v_factor_levels))) +
  labs(x = "Developmental Stage", y = "Number of Genes") +
  geom_bar(position = "dodge") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  guides(fill=guide_legend(title="Methylation Level"))

results <- meth_exp[complete.cases(meth_exp),]
high <- subset(results, meth_category == "high")
med <- subset(results, meth_category == "medium")
low <- subset(results, meth_category == "low")

################################################################################
mean(results$Counts)
var(resultsCounts)
hist(results$Counts)

model <- glm(Counts~meth_category*DStage, family = quasipoisson, data=results)

library(car)
#Report the interaction stats from Anova. You are doing a GLM with a quasipoisson distribution
Anova(model)
summary (model)

model$deviance/model$df.residual #Overdispersion > 1.2 bad

#pseudo R2
100*((model$null.deviance-model$deviance)/model$null.deviance) #pseudo R2

library(emmeans)
emmeans(model, pairwise~meth_category*DStage )
