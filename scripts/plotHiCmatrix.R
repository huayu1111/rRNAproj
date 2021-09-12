library(Matrix)
library(pheatmap)

mapBin <- function(vec,binMat){
   vec <- unlist(strsplit(vec,"\t"))
   binIndexes <- binMat[which(binMat[,1] == vec[1] & ((binMat[,2] >= as.numeric(vec[2]) & (binMat[,2]) <= as.numeric(vec[3])) | (binMat[,3] >= as.numeric(vec[2]) & (binMat[,3]) <= as.numeric(vec[3])) | (binMat[,2] <= as.numeric(vec[2]) & (binMat[,3]) >= as.numeric(vec[3])) | (binMat[,2] >= as.numeric(vec[2]) & (binMat[,3]) <= as.numeric(vec[3])))),4]
   return(binIndexes)
}

conDataE14 <- read.table(file="/media/data/hicdata/E14/matrix/E14/kr/E14_pearson_150000_kr.matrix",sep="\t",header=FALSE,stringsAsFactors=FALSE)
freqMatrixE14 <- sparseMatrix(as.numeric(conDataE14[,1]), as.numeric(conDataE14[,2]), x = as.numeric(conDataE14[,3]))

conDataCX <- read.table(file="/media/data/hicdata/CX/matrix/CX/kr/CX_pearson_150000_kr.matrix",sep="\t",header=FALSE,stringsAsFactors=FALSE)
freqMatrixCX <- sparseMatrix(as.numeric(conDataCX[,1]), as.numeric(conDataCX[,2]), x = as.numeric(conDataCX[,3]))

mapInfoE14 <- read.table("/media/data/hicdata/E14/matrix/E14/raw/150000/E14_150000_abs.bed",sep="\t",header=FALSE,stringsAsFactors=FALSE)
geneInfo <- read.table(file="/media/data/hicdata/annodata/mouse_genes.bed",skip=7,sep="\t",header=FALSE,stringsAsFactors=FALSE)

HubData <- read.table(file="/media/data/hicdata/annodata/inactive_nucleolar_hub_mm10.bed",sep="\t",header=FALSE,stringsAsFactors=FALSE)
HubVec <- paste(HubData[,1],HubData[,2],HubData[,3],sep="\t")
HubBins <- as.numeric(unlist(sapply(HubVec,mapBin,mapInfoE14,simplify=TRUE)))
HubBins <- unique(HubBins)
HubChrs <- mapInfoE14[HubBins,1]

matchIndex <- match("ENSMUSG00000075046",geneInfo[,4])
geneVec <- paste(geneInfo[matchIndex,1],geneInfo[matchIndex,2],geneInfo[matchIndex,3],sep="\t")
DuxBins <- as.numeric(unlist(sapply(geneVec,mapBin,mapInfoE14,simplify=TRUE)))
DuxBins <- unique(DuxBins)

AllBins <- unique(c(HubBins,seq(DuxBins-4,DuxBins+5)))
annotationRow <- c(HubChrs,rep("Dux",10))
annotationRow <- data.frame(GeneClass=factor(annotationRow,levels=c("chr10","chr12","chr15","chr16","chr18","chr19","Dux")))
rownames(annotationRow) <- paste("Bin",seq(1,length(AllBins)),sep="")
annotationCol <- annotationRow

plotData <- freqMatrixE14[AllBins,AllBins] + t(freqMatrixE14[AllBins,AllBins])
rownames(plotData) <- paste("Bin",seq(1,length(AllBins)),sep="")
colnames(plotData) <- paste("Bin",seq(1,length(AllBins)),sep="")
png(file="/media/data/hicdata/E14/ABCompartmentE14DuxAll.png",units='in',width=10,height=10,res=600)
pheatmap(plotData,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=annotationRow,annotation_col=annotationRow,annotation_names_row=FALSE,annotation_names_col=FALSE, breaks = c(seq(0,0.05,0.0005)),color = colorRampPalette(c("#7676E5","white","#FF3333"))(100),fontsize=4,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

HubData <- read.table(file="/media/data/hicdata/annodata/inactive_nucleolar_hub_mm10.bed",sep="\t",header=FALSE,stringsAsFactors=FALSE)
HubVec <- paste(HubData[,1],HubData[,2],HubData[,3],sep="\t")
HubBins <- as.numeric(unlist(sapply(HubVec,mapBin,mapInfoE14,simplify=TRUE)))
HubBins <- unique(HubBins)

matchIndex <- match("ENSMUSG00000075046",geneInfo[,4])
geneVec <- paste(geneInfo[matchIndex,1],geneInfo[matchIndex,2],geneInfo[matchIndex,3],sep="\t")
DuxBins <- as.numeric(unlist(sapply(geneVec,mapBin,mapInfoE14,simplify=TRUE)))
DuxBins <- unique(DuxBins)

AllBins <- unique(c(HubBins,seq(DuxBins-4,DuxBins+5))) 
annotationRow <- c(HubChrs,rep("Dux",10))
annotationRow <- data.frame(GeneClass=factor(annotationRow,levels=c("chr10","chr12","chr15","chr16","chr18","chr19","Dux")))
rownames(annotationRow) <- paste("Bin",seq(1,length(AllBins)),sep="")
annotationCol <- annotationRow

plotData <- freqMatrixCX[AllBins,AllBins] + t(freqMatrixCX[AllBins,AllBins])
rownames(plotData) <- paste("Bin",seq(1,length(AllBins)),sep="")
colnames(plotData) <- paste("Bin",seq(1,length(AllBins)),sep="")
png(file="/media/data/hicdata/CX/ABCompartmentCXDuxAll.png",units='in',width=10,height=10,res=600)
pheatmap(plotData,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=annotationRow,annotation_col=annotationRow,annotation_names_row=FALSE,annotation_names_col=FALSE,breaks = c(seq(0,0.05,0.0005)), color = colorRampPalette(c("#7676E5","white","#FF3333"))(100),fontsize=4,show_rownames=FALSE,show_colnames=FALSE)
dev.off()