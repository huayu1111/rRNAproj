library(pheatmap)
library(sva)

expData_stress <- read.table(file = "/media/wgs/InhouseData/gene_TPM_stress.txt", sep = "\t", header = T, row.names = 1)
rownames(expData_stress) <- expData_stress[,1]
expData_stress <- cbind(rowMeans(expData_stress[,c(3,4)]),rowMeans(expData_stress[,c(5,6)]),rowMeans(expData_stress[,c(7,8)]),rowMeans(expData_stress[,c(9,10)]),rowMeans(expData_stress[,c(11,12)]))
colnames(expData_stress) <- c("CX-5461","Control","Etoposide","Rapamycin","Rotonone")

expData_ES <- read.table(file = "/media/wgs/InhouseData/gene_TPM_ES.txt", sep = "\t", header = T, row.names = 1)
rownames(expData_ES) <- expData_ES[,1]
expData_ES <- cbind(rowMeans(expData_ES[,c(3,4,5)]),rowMeans(expData_ES[,c(15,16,17)]))
colnames(expData_ES) <- c("Lin28_KO","Lin28_WT")

expData_2C <- read.table(file = "/media/wgs/InhouseData/gene_TPM_2CTomato.txt", sep = "\t", header = T, row.names = 1)
rownames(expData_2C) <- expData_2C[,1]
expData_2C <- expData_2C[,seq(3,ncol(expData_2C))]
colnames(expData_2C) <- c("2CTomato+","2CTomato-")

expData_Embryo <- read.table(file="/media/wgs/EmbryoData/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_Embryo) <- expData_Embryo[,1]
expData_Embryo <- expData_Embryo[,seq(3,ncol(expData_Embryo))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_Embryo),sampInfo[,1])
colnames(expData_Embryo) <- as.character(sampInfo[matchIndexes,3])
expData_Embryo <- cbind(rowMeans(expData_Embryo[,c(1,2)]),rowMeans(expData_Embryo[,c(3,4)]),rowMeans(expData_Embryo[,c(5,6)]),rowMeans(expData_Embryo[,c(7,8)]),rowMeans(expData_Embryo[,c(9,10,11)]),rowMeans(expData_Embryo[,c(12,13,14,15)]),rowMeans(expData_Embryo[,c(16,17)]),rowMeans(expData_Embryo[,c(18,19)]))
colnames(expData_Embryo) <- c("zygote","2-cell","4-cell","8-cell","ICM","mESC","MII-Oocyte","early 2-cell")

expData_Dux <- read.table(file="/media/wgs/DuxData/RNAseqData/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_Dux) <- expData_Dux[,1]
expData_Dux <- expData_Dux[,seq(3,ncol(expData_Dux))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_Dux),sampInfo[,1])
colnames(expData_Dux) <- as.character(sampInfo[matchIndexes,3])
expData_Dux <- cbind(rowMeans(expData_Dux[,c(8,9)]),rowMeans(expData_Dux[,c(10,11)]))
colnames(expData_Dux) <- c("Dux-GFPpos","Dux-GFPneg")

expData_NELFA <- read.table(file="/media/wgs/NELFAData/RNAseqData/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_NELFA) <- expData_NELFA[,1]
expData_NELFA <- expData_NELFA[,seq(3,ncol(expData_NELFA))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_NELFA),sampInfo[,1])
colnames(expData_NELFA) <- as.character(sampInfo[matchIndexes,3])
expData_NELFA <- cbind(rowMeans(expData_NELFA[,c(9,11)]),rowMeans(expData_NELFA[,c(10,12)]))
colnames(expData_NELFA) <- c("NELFA-GFPpos","NELFA-GFPneg")

expData_CAF1 <- read.table(file="/media/wgs/CAF1Data/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_CAF1) <- expData_CAF1[,1]
expData_CAF1 <- expData_CAF1[,seq(3,ncol(expData_CAF1))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_CAF1),sampInfo[,1])
colnames(expData_CAF1) <- as.character(sampInfo[matchIndexes,3])
expData_CAF1 <- cbind(rowMeans(expData_CAF1[,c(1,2)]),rowMeans(expData_CAF1[,c(3,4)]))
colnames(expData_CAF1) <- c("CAF1_WT","CAF1_KD")

expData_Kap1 <- read.table(file="/media/wgs/Kap1Data/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_Kap1) <- expData_Kap1[,1]
expData_Kap1 <- expData_Kap1[,seq(3,ncol(expData_Kap1))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_Kap1),sampInfo[,1])
colnames(expData_Kap1) <- as.character(sampInfo[matchIndexes,3])
expData_Kap1 <- cbind(rowMeans(expData_Kap1[,c(1,2,3)]),rowMeans(expData_Kap1[,c(4,5,6)]))
colnames(expData_Kap1) <- c("Kap1_WT","Kap1_KO")

expData_Zscan4 <- read.table(file="/media/wgs/Zscan4Data/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_Zscan4) <- expData_Zscan4[,1]
expData_Zscan4 <- expData_Zscan4[,seq(3,ncol(expData_Zscan4))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_Zscan4),sampInfo[,1])
colnames(expData_Zscan4) <- as.character(sampInfo[matchIndexes,3])
expData_Zscan4 <- cbind(rowMeans(expData_Zscan4[,c(1,2)]),rowMeans(expData_Zscan4[,c(3,4)]))
colnames(expData_Zscan4) <- c("Zscan4_KO","Zscan4_WT")

expData_Dppa <- read.table(file="/media/wgs/DappaData/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_Dppa) <- expData_Dppa[,1]
expData_Dppa <- expData_Dppa[,seq(3,ncol(expData_Dppa))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_Dppa),sampInfo[,1])
colnames(expData_Dppa) <- as.character(sampInfo[matchIndexes,3])
expData_Dppa <- expData_Dppa[,c(9,10)]
colnames(expData_Dppa) <- c("Dppa4_GFPneg","Dppa4_GFPpos")

expData_LINE1 <- read.table(file="/media/wgs/LINE1Data/stringtiefile/gene_TPM.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(expData_LINE1) <- expData_LINE1[,1]
expData_LINE1 <- expData_LINE1[,seq(3,ncol(expData_LINE1))]
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(colnames(expData_LINE1),sampInfo[,1])
colnames(expData_LINE1) <- as.character(sampInfo[matchIndexes,3])
expData_LINE1 <- cbind(rowMeans(expData_LINE1[,c(1,2,3)]),rowMeans(expData_LINE1[,c(4,5,6)]))
colnames(expData_LINE1) <- c("RC ASO","LINE1 ASO")

genenames <- Reduce(intersect,list(nameStress=rownames(expData_stress),name2C=rownames(expData_2C),nameEmbryo=rownames(expData_Embryo),nameES=rownames(expData_ES),nameDux=rownames(expData_Dux),nameCAF1=rownames(expData_CAF1),nameKap1=rownames(expData_Kap1),nameZscan4=rownames(expData_Zscan4),nameDppa=rownames(expData_Dppa),nameLINE1=rownames(expData_LINE1)))

expData <- data.matrix(cbind(expData_stress[genenames,],expData_2C[genenames,],expData_Embryo[genenames,],expData_ES[genenames,],expData_Dux[genenames,],expData_NELFA[genenames,],expData_CAF1[genenames,],expData_Kap1[genenames,],expData_Zscan4[genenames,],expData_Dppa[genenames,],expData_LINE1[genenames,]))
expData <- expData[which(rowSums(expData)>0),]

expData <- log2(expData+1)
batch <- as.factor(rep(c("Stress","2C","Embryo","ES","Dux","NELFA","CAF1","Kap1","Zscan4","Dppa","LINE1"), c(ncol(expData_stress),ncol(expData_2C),ncol(expData_Embryo),ncol(expData_ES),ncol(expData_Dux),ncol(expData_NELFA),ncol(expData_CAF1),ncol(expData_Kap1),ncol(expData_Zscan4),ncol(expData_Dppa),ncol(expData_LINE1))))
modcombat = model.matrix(~1,data=batch)
combat_edata = ComBat(dat=expData, batch=batch, mod=modcombat, par.prior=TRUE)

combat_edata <- combat_edata[,c(1,2,seq(4,ncol(combat_edata)))]
pdf(file="/media/wgs/InhouseData/sample_clustering.pdf",width=20,height=15)
p <- pheatmap(combat_edata,show_rownames=FALSE,scale="row",cluster_rows=FALSE,cluster_cols=TRUE,fontsize=20)
plot(p$tree_col)
dev.off()