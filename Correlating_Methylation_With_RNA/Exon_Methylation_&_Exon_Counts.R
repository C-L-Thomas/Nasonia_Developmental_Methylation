library(readr)
library(reshape2)
library(ggplot2)

#Read in dataset
Meth_data <- read_delim("~/Total_Exons.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

Meth_data$e1 <- (Meth_data$numCs1/Meth_data$coverage1)*100
Meth_data$e2<- (Meth_data$numCs2/Meth_data$coverage2)*100
Meth_data$l1 <- (Meth_data$numCs3/Meth_data$coverage3)*100
Meth_data$l2 <- (Meth_data$numCs4/Meth_data$coverage4)*100
Meth_data$l3 <- (Meth_data$numCs5/Meth_data$coverage5)*100
Meth_data$pr1 <- (Meth_data$numCs6/Meth_data$coverage6)*100
Meth_data$pr2 <- (Meth_data$numCs7/Meth_data$coverage7)*100
Meth_data$pr3 <- (Meth_data$numCs8/Meth_data$coverage8)*100
Meth_data$pu1 <- (Meth_data$numCs9/Meth_data$coverage9)*100
Meth_data$pu2 <- (Meth_data$numCs10/Meth_data$coverage10)*100
Meth_data$pu3 <- (Meth_data$numCs11/Meth_data$coverage11)*100
Meth_data$a1 <- (Meth_data$numCs12/Meth_data$coverage12)*100
Meth_data$a2 <- (Meth_data$numCs13/Meth_data$coverage13)*100
Meth_data$a3 <- (Meth_data$numCs14/Meth_data$coverage14)*100

Meth_data$embryo <- rowMeans(subset(Meth_data, select = c(52,53)))
Meth_data$larva <- rowMeans(subset(Meth_data, select = c(54,55,56)))
Meth_data$prepupae <- rowMeans(subset(Meth_data, select = c(57,58,59)))
Meth_data$pupae <- rowMeans(subset(Meth_data, select = c(60,61,62)))
Meth_data$adult<- rowMeans(subset(Meth_data, select = c(63,64,65)))

cols <- c(50,51,66,67,68,69,70)
Methylated_Exons <- Meth_data[,c(cols)]
Methylated_Exons$unique <- paste(Methylated_Exons$geneID, sep ="_",Methylated_Exons$number)
cols <- c(8,3,4,5,6,7)
Methylated_Exons <- Methylated_Exons[,c(cols)]
Bi_data <- melt(Methylated_Exons, value.name = "Percentage")
colnames(Bi_data) <- c("ExonID", "DStage", "Percentage")
nrow(Bi_data)


#Read in transcriptome results
transcriptome <- read_delim("~/Normalised_Exons.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
transcriptome$adult<- rowMeans(subset(transcriptome, select = c(2,3,4)))
transcriptome$embryo <- rowMeans(subset(transcriptome, select = c(5,6,7)))
transcriptome$larva <- rowMeans(subset(transcriptome, select = c(8,9,10)))
transcriptome$prepupae <- rowMeans(subset(transcriptome, select = c(11,12,13)))
transcriptome$pupae <- rowMeans(subset(transcriptome, select = c(14,15,16)))

cols <- c(1,31,32,33,34,35)
RNA <- transcriptome[,c(cols)]

library(tidyverse)
here <- RNA %>% separate_wider_delim(V1,":E", names = c("gene", "exon"))
here$exon = str_remove(here$exon, "^0+")
here$unique <- paste(here$gene, sep ="_",here$exon)

cols <- c(8,4,5,6,7,3)
RNA <- here[,c(cols)]

RNA_data <- melt(RNA, value.name = "Counts")
colnames(RNA_data) <- c("ExonID", "DStage", "Counts")
#Remaining<- RNA_data[(RNA_data$geneID %in% Bi_data$geneID),]
#nrow(Remaining)

meth_exp <- merge(Bi_data, RNA_data, by=c("ExonID","DStage"), all=T)
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
var(results$Counts)
mean(results$Counts)# mean not equal to var so use quasipoisson

library(emmeans)

model <- glm(Counts~Percentage*DStage, family=quasipoisson,data=results)

library(car)
#Report the interaction stats from Anova. You are doing a GLM with a quasipoisson distribution
Anova(model)

model$deviance/model$df.residual #Overdispersion > 1.2 bad

#pseudo R2
100*((model$null.deviance-model$deviance)/model$null.deviance) #pseudo R2

#The slopes (estimates) and significant stats could be reported here
emtrends(model, pairwise~DStage, var = "Percentage")
emmip(model, DStage ~ Percentage, cov.reduce = range)

#Something like this graph might be useful
ggplot(results, aes(x = Percentage, y = log(Counts), colour = DStage)) +
  facet_wrap(~DStage, nrow=3) +
  geom_point() +
  theme_classic() +
  geom_line(data = cbind(results, pred = predict(model2)), aes(y = log(pred))) +
  theme(legend.position = "none")
