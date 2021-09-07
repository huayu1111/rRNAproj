library(reshape2)
library(ggplot2)

targetedGenes <- c("Rpl22l1","Spata5","Rps29")
resdata <- read.table("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_CX.csv",header=TRUE,row.names=1,sep=",",stringsAsFactors=FALSE)
nonnaIndexes <- which(!is.na(resdata$pvalue))
resdata <- resdata[nonnaIndexes,]
matchIndexes <- match(targetedGenes,resdata$gene_name)
resdata$classlabels <- rep("Nonsignificance",nrow(resdata))
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange > 0.5)] <- "Upregulation"
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange < -0.5)] <- "Downregulation"
resdata$size <- rep(0.5,nrow(resdata))
resdata$pvalue <- -log10(resdata$pvalue)
infIndexes <- which(is.infinite(resdata$pvalue))
noninfIndexes <- which(!is.infinite(resdata$pvalue))
maxValue <- max(resdata$pvalue[noninfIndexes])
resdata$pvalue[infIndexes] <- maxValue

png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_CX_ribogenes.png",width=12,height=10,units="cm",res=500)
p <- ggplot(resdata,aes(x=log2FoldChange,y=pvalue,colour=factor(resdata$classlabels,levels=c("Downregulation","Nonsignificance","Upregulation")))) + geom_point(size=resdata$size,shape=20,show.legend = F) + scale_colour_manual(values=c("blue","gray","red")) + geom_point(data=resdata[matchIndexes,], aes(x=log2FoldChange, y=pvalue), colour="black", size=3) # + 
p <- p + xlab(expression(log[2]("Fold change"))) + ylab(expression(-log[10]("P value")))
p <- p + theme(legend.title=element_blank())
p <- p + scale_y_continuous(limits = c(0,150))
p <- p + scale_x_continuous(limits = c(-8,8))
p <- p + annotate(geom="text", x=-7, y=145, label=paste("Down:",length(which(resdata$classlabels=="Downregulation")),sep=""), color="blue")
p <- p + annotate(geom="text", x=7, y=145, label=paste("Up:",length(which(resdata$classlabels=="Upregulation")),sep=""), color="red")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()

targetedGenes <- c("Rpl22l1","Spata5","Rps29")
resdata <- read.table("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_Etop.csv",header=TRUE,row.names=1,sep=",",stringsAsFactors=FALSE)
nonnaIndexes <- which(!is.na(resdata$pvalue))
resdata <- resdata[nonnaIndexes,]
matchIndexes <- match(targetedGenes,resdata$gene_name)
resdata$classlabels <- rep("Nonsignificance",nrow(resdata))
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange > 0.5)] <- "Upregulation"
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange < -0.5)] <- "Downregulation"
resdata$size <- rep(0.5,nrow(resdata))
resdata$pvalue <- -log10(resdata$pvalue)
infIndexes <- which(is.infinite(resdata$pvalue))
noninfIndexes <- which(!is.infinite(resdata$pvalue))
maxValue <- max(resdata$pvalue[noninfIndexes])
resdata$pvalue[infIndexes] <- maxValue

# pdf("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_Etop.pdf",width=20,height=10)
png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_Etop_ribogenes.png",width=12,height=10,units="cm",res=500)
p <- ggplot(resdata,aes(x=log2FoldChange,y=pvalue,colour=factor(resdata$classlabels,levels=c("Downregulation","Nonsignificance","Upregulation")))) + geom_point(size=resdata$size,shape=20,show.legend = F) + scale_colour_manual(values=c("blue","gray","red")) + geom_point(data=resdata[matchIndexes,], aes(x=log2FoldChange, y=pvalue), colour="black", size=3) # + geom_text(data=resdata[matchIndexes,],aes(x=log2FoldChange, y=pvalue,label=targetedGenes),colour="green",size=1.5)
p <- p + xlab(expression(log[2]("Fold change"))) + ylab(expression(-log[10]("P value")))
p <- p + theme(legend.title=element_blank())
p <- p + scale_y_continuous(limits = c(0,150))
p <- p + scale_x_continuous(limits = c(-8,8))
p <- p + annotate(geom="text", x=-7, y=145, label=paste("Down:",length(which(resdata$classlabels=="Downregulation")),sep=""), color="blue")
p <- p + annotate(geom="text", x=7, y=145, label=paste("Up:",length(which(resdata$classlabels=="Upregulation")),sep=""), color="red")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()

targetedGenes <- c("Rpl22l1","Spata5","Rps29")
resdata <- read.table("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_Rapa.csv",header=TRUE,row.names=1,sep=",",stringsAsFactors=FALSE)
nonnaIndexes <- which(!is.na(resdata$pvalue))
resdata <- resdata[nonnaIndexes,]
matchIndexes <- match(targetedGenes,resdata$gene_name)
resdata$classlabels <- rep("Nonsignificance",nrow(resdata))
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange > 0.5)] <- "Upregulation"
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange < -0.5)] <- "Downregulation"
resdata$size <- rep(0.5,nrow(resdata))
resdata$pvalue <- -log10(resdata$pvalue)
infIndexes <- which(is.infinite(resdata$pvalue))
noninfIndexes <- which(!is.infinite(resdata$pvalue))
maxValue <- max(resdata$pvalue[noninfIndexes])
resdata$pvalue[infIndexes] <- maxValue

png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_Rapa_ribogenes.png",width=12,height=10,units="cm",res=500)
p <- ggplot(resdata,aes(x=log2FoldChange,y=pvalue,colour=factor(resdata$classlabels,levels=c("Downregulation","Nonsignificance","Upregulation")))) + geom_point(size=resdata$size,shape=20,show.legend = F) + scale_colour_manual(values=c("blue","gray","red")) + geom_point(data=resdata[matchIndexes,], aes(x=log2FoldChange, y=pvalue), colour="black", size=3) # + 
p <- p + xlab(expression(log[2]("Fold change"))) + ylab(expression(-log[10]("P value")))
p <- p + theme(legend.title=element_blank())
p <- p + scale_y_continuous(limits = c(0,150))
p <- p + scale_x_continuous(limits = c(-8,8))
p <- p + annotate(geom="text", x=-7, y=145, label=paste("Down:",length(which(resdata$classlabels=="Downregulation")),sep=""), color="blue")
p <- p + annotate(geom="text", x=7, y=145, label=paste("Up:",length(which(resdata$classlabels=="Upregulation")),sep=""), color="red")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()

targetedGenes <- c("Rpl22l1","Spata5","Rps29")
resdata <- read.table("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_Rote.csv",header=TRUE,row.names=1,sep=",",stringsAsFactors=FALSE)
nonnaIndexes <- which(!is.na(resdata$pvalue))
resdata <- resdata[nonnaIndexes,]
matchIndexes <- match(targetedGenes,resdata$gene_name)
resdata$classlabels <- rep("Nonsignificance",nrow(resdata))
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange > 0.5)] <- "Upregulation"
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange < -0.5)] <- "Downregulation"
resdata$size <- rep(0.5,nrow(resdata))
resdata$pvalue <- -log10(resdata$pvalue)
infIndexes <- which(is.infinite(resdata$pvalue))
noninfIndexes <- which(!is.infinite(resdata$pvalue))
maxValue <- max(resdata$pvalue[noninfIndexes])
resdata$pvalue[infIndexes] <- maxValue

png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/CT_vs_Rote_ribogenes.png",width=12,height=10,units="cm",res=500)
p <- ggplot(resdata,aes(x=log2FoldChange,y=pvalue,colour=factor(resdata$classlabels,levels=c("Downregulation","Nonsignificance","Upregulation")))) + geom_point(size=resdata$size,shape=20,show.legend = F) + scale_colour_manual(values=c("blue","gray","red")) + geom_point(data=resdata[matchIndexes,], aes(x=log2FoldChange, y=pvalue), colour="black", size=3)
p <- p + xlab(expression(log[2]("Fold change"))) + ylab(expression(-log[10]("P value")))
p <- p + theme(legend.title=element_blank())
p <- p + scale_y_continuous(limits = c(0,150))
p <- p + scale_x_continuous(limits = c(-8,8))
p <- p + annotate(geom="text", x=-7, y=145, label=paste("Down:",length(which(resdata$classlabels=="Downregulation")),sep=""), color="blue")
p <- p + annotate(geom="text", x=7, y=145, label=paste("Up:",length(which(resdata$classlabels=="Upregulation")),sep=""), color="red")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()