library(methylKit)

file.list <- list("DEMB5_evidence.cov", "DEMB6_evidence.cov", "DLAR4_evidence.cov",
"DLAR5_evidence.cov", "DLAR6_evidence.cov", "DPRE1_evidence.cov", "DPRE2_evidence.cov", "DPRE3_evidence.cov",
"DPUP2_evidence.cov","DPUP3_evidence.cov", "DPUP4_evidence.cov", "DADU1_evidence.cov", "DADU2_evidence.cov","DADU4_evidence.cov")

treat_conditions <- c(0,0,1,1,1,2,2,2,3,3,3,4,4,4)

raw_data <- methRead(file.list,
                     sample.id = list("EMB5","EMB6","LAR4","LAR5","LAR6","PRE1","PRE2","PRE3","PUP2","PUP3","PUP4","ADU1","ADU2","ADU4"),
                     treatment = treat_conditions,
                     assembly="Nvit_PSR_1.1",
                     dbtype = NA,
                     pipeline = "bismarkCoverage",
                     header = T, sep = "\t", mincov=1,
                     context="CpG")

filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

meth_all_data <- unite(filtered_data, destrand=FALSE ,min.per.group = 1L)
nrow(meth_all_data)

df_meth_all <- getData(meth_all_data)
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

e1 <- df_meth_all[,5:6]
e2 <- df_meth_all[,8:9]
l1 <- df_meth_all[,11:12]
l2 <- df_meth_all[,14:15]
l3 <- df_meth_all[,17:18]
pr1 <- df_meth_all[,20:21]
pr2 <- df_meth_all[,23:24]
pr3 <- df_meth_all[,26:27]
pu1 <- df_meth_all[,29:30]
pu2 <- df_meth_all[,32:33]
pu3 <- df_meth_all[,35:36]
a1 <- df_meth_all[,38:39]
a2 <- df_meth_all[,41:42]
a3 <- df_meth_all[,44:45]

e1bt <- function(a, b, p = 0.009) {binom.test(a, b, 0.009, alternative="greater")$p.value}
e2bt <- function(a, b, p = 0.009) {binom.test(a, b, 0.009, alternative="greater")$p.value}
l1bt <- function(a, b, p = 0.008) {binom.test(a, b, 0.008, alternative="greater")$p.value}
l2bt <- function(a, b, p = 0.009) {binom.test(a, b, 0.009, alternative="greater")$p.value}
l3bt <- function(a, b, p = 0.009) {binom.test(a, b, 0.009, alternative="greater")$p.value}
pr1bt <- function(a, b, p = 0.003) {binom.test(a, b, 0.003, alternative="greater")$p.value}
pr2bt <- function(a, b, p = 0.003) {binom.test(a, b, 0.003, alternative="greater")$p.value}
pr3bt <- function(a, b, p = 0.003) {binom.test(a, b, 0.003, alternative="greater")$p.value}
pu1bt <- function(a, b, p = 0.003) {binom.test(a, b, 0.003, alternative="greater")$p.value}
pu2bt <- function(a, b, p = 0.003) {binom.test(a, b, 0.003, alternative="greater")$p.value}
pu3bt <- function(a, b, p = 0.009) {binom.test(a, b, 0.009, alternative="greater")$p.value}
a1bt <- function(a, b, p = 0.003) {binom.test(a, b, 0.003, alternative="greater")$p.value}
a2bt <- function(a, b, p = 0.004) {binom.test(a, b, 0.004, alternative="greater")$p.value}
a3bt <- function(a, b, p = 0.009) {binom.test(a, b, 0.009, alternative="greater")$p.value}

colnames(e1) <- c("CT", "Ccount")
e1 <- e1[complete.cases(e1),]
e1$pVal <- mapply(e1bt, e1$Ccount, e1$CT)
e1$FDR <- p.adjust(e1$pVal, method = "BH", n = length(e1$pVal))
e1$row <- as.numeric(rownames(e1))
e1_meth_sites <- subset(e1, e1$FDR < 0.05)
allrows <- rbind(allrows, e1_meth_sites)
nrow(e1_meth_sites)

colnames(e2) <- c("CT", "Ccount")
e2 <- e2[complete.cases(e2),]
e2$pVal <- mapply(e2bt, e2$Ccount, e2$CT)
e2$FDR <- p.adjust(e2$pVal, method = "BH", n = length(e2$pVal))
e2$row <- as.numeric(rownames(e2))
e2_meth_sites <- subset(e2, e2$FDR < 0.05)
allrows <- rbind(allrows, e2_meth_sites)
nrow(e2_meth_sites)

colnames(l1) <- c("CT", "Ccount")
l1 <- l1[complete.cases(l1),]
l1$pVal <- mapply(l1bt, l1$Ccount, l1$CT)
l1$FDR <- p.adjust(l1$pVal, method = "BH", n = length(l1$pVal))
l1$row <- as.numeric(rownames(l1))
l1_meth_sites <- subset(l1, l1$FDR < 0.05)
allrows <- rbind(allrows, l1_meth_sites)
nrow(l1_meth_sites)

colnames(l2) <- c("CT", "Ccount")
l2 <- l2[complete.cases(l2),]
l2$pVal <- mapply(l2bt, l2$Ccount, l2$CT)
l2$FDR <- p.adjust(l2$pVal, method = "BH", n = length(l2$pVal))
l2$row <- as.numeric(rownames(l2))
l2_meth_sites <- subset(l2, l2$FDR < 0.05)
allrows <- rbind(allrows, l2_meth_sites)
nrow(l2_meth_sites)

colnames(l3) <- c("CT", "Ccount")
l3 <- l3[complete.cases(l3),]
l3$pVal <- mapply(l3bt, l3$Ccount, l3$CT)
l3$FDR <- p.adjust(l3$pVal, method = "BH", n = length(l3$pVal))
l3$row <- as.numeric(rownames(l3))
l3_meth_sites <- subset(l3, l3$FDR < 0.05)
allrows <- rbind(allrows, l3_meth_sites)
nrow(l3_meth_sites)

colnames(pr1) <- c("CT", "Ccount")
pr1 <- pr1[complete.cases(pr1),]
pr1$pVal <- mapply(pr1bt, pr1$Ccount, pr1$CT)
pr1$FDR <- p.adjust(pr1$pVal, method = "BH", n = length(pr1$pVal))
pr1$row <- as.numeric(rownames(pr1))
pr1_meth_sites <- subset(pr1, pr1$FDR < 0.05)
allrows <- rbind(allrows, pr1_meth_sites)
nrow(pr1_meth_sites)

colnames(pr2) <- c("CT", "Ccount")
pr2 <- pr2[complete.cases(pr2),]
pr2$pVal <- mapply(pr2bt, pr2$Ccount, pr2$CT)
pr2$FDR <- p.adjust(pr2$pVal, method = "BH", n = length(pr2$pVal))
pr2$row <- as.numeric(rownames(pr2))
pr2_meth_sites <- subset(pr2, pr2$FDR < 0.05)
allrows <- rbind(allrows, pr2_meth_sites)
nrow(pr2_meth_sites)

colnames(pr3) <- c("CT", "Ccount")
pr3 <- pr3[complete.cases(pr3),]
pr3$pVal <- mapply(pr3bt, pr3$Ccount, pr3$CT)
pr3$FDR <- p.adjust(pr3$pVal, method = "BH", n = length(pr3$pVal))
pr3$row <- as.numeric(rownames(pr3))
pr3_meth_sites <- subset(pr3, pr3$FDR < 0.05)
allrows <- rbind(allrows, pr3_meth_sites)
nrow(pr3_meth_sites)

colnames(pu1) <- c("CT", "Ccount")
pu1 <- pu1[complete.cases(pu1),]
pu1$pVal <- mapply(pu1bt, pu1$Ccount, pu1$CT)
pu1$FDR <- p.adjust(pu1$pVal, method = "BH", n = length(pu1$pVal))
pu1$row <- as.numeric(rownames(pu1))
pu1_meth_sites <- subset(pu1, pu1$FDR < 0.05)
allrows <- rbind(allrows, pu1_meth_sites)
nrow(pu1_meth_sites)

colnames(pu2) <- c("CT", "Ccount")
pu2 <- pu2[complete.cases(pu2),]
pu2$pVal <- mapply(pu2bt, pu2$Ccount, pu2$CT)
pu2$FDR <- p.adjust(pu2$pVal, method = "BH", n = length(pu2$pVal))
pu2$row <- as.numeric(rownames(pu2))
pu2_meth_sites <- subset(pu2, pu2$FDR < 0.05)
allrows <- rbind(allrows, pu2_meth_sites)
nrow(pu2_meth_sites)

colnames(pu3) <- c("CT", "Ccount")
pu3 <- pu3[complete.cases(pu3),]
pu3$pVal <- mapply(pu3bt, pu3$Ccount, pu3$CT)
pu3$FDR <- p.adjust(pu3$pVal, method = "BH", n = length(pu3$pVal))
pu3$row <- as.numeric(rownames(pu3))
pu3_meth_sites <- subset(pu3, pu3$FDR < 0.05)
allrows <- rbind(allrows, pu3_meth_sites)
nrow(pu3_meth_sites)

colnames(a1) <- c("CT", "Ccount")
a1 <- a1[complete.cases(a1),]
a1$pVal <- mapply(a1bt, a1$Ccount, a1$CT)
a1$FDR <- p.adjust(a1$pVal, method = "BH", n = length(a1$pVal))
a1$row <- as.numeric(rownames(a1))
a1_meth_sites <- subset(a1, a1$FDR < 0.05)
allrows <- rbind(allrows, a1_meth_sites)
nrow(a1_meth_sites)

colnames(a2) <- c("CT", "Ccount")
a2 <- a2[complete.cases(a2),]
a2$pVal <- mapply(a2bt, a2$Ccount, a2$CT)
a2$FDR <- p.adjust(a2$pVal, method = "BH", n = length(a2$pVal))
a2$row <- as.numeric(rownames(a2))
a2_meth_sites <- subset(a2, a2$FDR < 0.05)
allrows <- rbind(allrows, a2_meth_sites)
nrow(a2_meth_sites)

colnames(a3) <- c("CT", "Ccount")
a3 <- a3[complete.cases(a3),]
a3$pVal <- mapply(a3bt, a3$Ccount, a3$CT)
a3$FDR <- p.adjust(a3$pVal, method = "BH", n = length(a3$pVal))
a3$row <- as.numeric(rownames(a3))
a3_meth_sites <- subset(a3, a3$FDR < 0.05)
allrows <- rbind(allrows, a3_meth_sites)
nrow(a3_meth_sites)

meth_positions <- as.vector(as.numeric(unique(allrows$row)))
length(meth_positions)

subset_methBase <- methylKit::select(meth_all_data, meth_positions)
methBase_ob <- getData(subset_methBase)
write.table(methBase_ob, file="Total_Methylated_Bases.txt", quote=F, row.names = F, sep = '\t')
nrow(methBase_ob)

save.image(file = "Counted_CpGs.RData")

#------------------------------------------------------------------------------

# Generate some quick plots which can be fixed up later
pdf("CpG_Correlation.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("CpG_Cluster.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("CpG_PCA.pdf")
PCASamples(subset_methBase, adj.lim = c(1.5,0.15))
dev.off()
