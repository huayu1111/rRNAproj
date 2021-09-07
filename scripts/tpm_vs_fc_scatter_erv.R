library(ggplot2)
expdata <- read.delim(file="/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile_erv/repeat_seq_TPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
targetedGenes <- c("MERVL-int|LTR|ERVL","MT2_Mm|LTR|ERVL","GSAT_MM|Satellite|Satellite")
plotdata <- aggregate(cbind(Control.1, Control.2, CX.1, CX.2) ~ gene_name, data = expdata, FUN = "sum")
rownames(plotdata) <- plotdata$gene_name

repclasses <- unlist(lapply(strsplit(plotdata$gene_name,"\\|"),function(x){ return(x[3])}))
matchIndexes <- grep("ERV1|ERVK|ERVL",repclasses)
nomatchIndexes <- setdiff(seq(1,length(repclasses)),matchIndexes)
repclasses[nomatchIndexes] <- "Other"
repclasses[grep("ERV",repclasses)] <- "ERV"  
plotdata$Size <- rep(0.01,nrow(plotdata))
plotdata$Size <- as.numeric(plotdata$Size)
plotdata$y <- log2(rowSums(plotdata[,c(4,5)])/rowSums(plotdata[,c(2,3)]))
plotdata$x <- log10(rowSums(plotdata[,c(2,3,4,5)]))
plotdata$group <- factor(repclasses,levels=c("ERV","Other"))


plotdata <- plotdata[which(!is.na(plotdata$y) & !is.infinite(plotdata$y)),]
matchindexes <- match(targetedGenes,rownames(plotdata))

png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile_erv/CT_vs_CX_erv.png",width=12,height=10,units="cm",res=600)
p <- ggplot(plotdata,aes(x=x,y=y,color=group)) + geom_point(size=0.2) + geom_point(data=plotdata[matchindexes,], aes(x=x, y=y), colour="black", size=2) + scale_color_manual(values =c('red','grey'))
p <- p + xlab(expression(log[10]("TPM"))) + ylab(expression(log[2]("FC (CX/Control)")))
p <- p + scale_y_continuous(limits = c(-5,5))
p <- p + geom_abline(intercept=0,slope=0,colour="grey")

p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()

plotdata <- aggregate(cbind(Control.1, Control.2, Etop.1, Etop.2) ~ gene_name, data = expdata, FUN = "sum")
rownames(plotdata) <- plotdata$gene_name

repclasses <- unlist(lapply(strsplit(plotdata$gene_name,"\\|"),function(x){ return(x[3])}))
matchIndexes <- grep("ERV1|ERVK|ERVL",repclasses)
nomatchIndexes <- setdiff(seq(1,length(repclasses)),matchIndexes)
repclasses[nomatchIndexes] <- "Other"
repclasses[grep("ERV",repclasses)] <- "ERV"  
plotdata$Size <- rep(0.01,nrow(plotdata))
plotdata$Size <- as.numeric(plotdata$Size)
plotdata$y <- log2(rowSums(plotdata[,c(4,5)])/rowSums(plotdata[,c(2,3)]))
plotdata$x <- log10(rowSums(plotdata[,c(2,3,4,5)]))
plotdata$group <- factor(repclasses,levels=c("ERV","Other"))

plotdata <- plotdata[which(!is.na(plotdata$y) & !is.infinite(plotdata$y)),]
matchindexes <- match(targetedGenes,rownames(plotdata))

png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile_erv/CT_vs_Etop_erv.png",width=12,height=10,units="cm",res=600)
p <- ggplot(plotdata,aes(x=x,y=y,color=group)) + geom_point(size=0.2) + geom_point(data=plotdata[matchindexes,], aes(x=x, y=y), colour="black", size=2) + scale_color_manual(values =c('red','grey'))
p <- p + xlab(expression(log[10]("TPM"))) + ylab(expression(log[2]("FC (CX/Control)")))
p <- p + scale_y_continuous(limits = c(-5,5))
p <- p + geom_abline(intercept=0,slope=0,colour="grey")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()


plotdata <- aggregate(cbind(Control.1, Control.2, Rapa.1, Rapa.2) ~ gene_name, data = expdata, FUN = "sum")
rownames(plotdata) <- plotdata$gene_name

repclasses <- unlist(lapply(strsplit(plotdata$gene_name,"\\|"),function(x){ return(x[3])}))
matchIndexes <- grep("ERV1|ERVK|ERVL",repclasses)
nomatchIndexes <- setdiff(seq(1,length(repclasses)),matchIndexes)
repclasses[nomatchIndexes] <- "Other"
repclasses[grep("ERV",repclasses)] <- "ERV"  
plotdata$Size <- rep(0.01,nrow(plotdata))
plotdata$Size <- as.numeric(plotdata$Size)
plotdata$y <- log2(rowSums(plotdata[,c(4,5)])/rowSums(plotdata[,c(2,3)]))
plotdata$x <- log10(rowSums(plotdata[,c(2,3,4,5)]))
plotdata$group <- factor(repclasses,levels=c("ERV","Other"))

plotdata <- plotdata[which(!is.na(plotdata$y) & !is.infinite(plotdata$y)),]
matchindexes <- match(targetedGenes,rownames(plotdata))

png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile_erv/CT_vs_Rapa_erv.png",width=12,height=10,units="cm",res=600)
p <- ggplot(plotdata,aes(x=x,y=y,color=group)) + geom_point(size=0.2) + geom_point(data=plotdata[matchindexes,], aes(x=x, y=y), colour="black", size=2) + scale_color_manual(values =c('red','grey'))
p <- p + xlab(expression(log[10]("TPM"))) + ylab(expression(log[2]("FC (CX/Control)")))
p <- p + scale_y_continuous(limits = c(-5,5))
p <- p + geom_abline(intercept=0,slope=0,colour="grey")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()

plotdata <- aggregate(cbind(Control.1, Control.2, Rote.1, Rote.2) ~ gene_name, data = expdata, FUN = "sum")
rownames(plotdata) <- plotdata$gene_name

repclasses <- unlist(lapply(strsplit(plotdata$gene_name,"\\|"),function(x){ return(x[3])}))
matchIndexes <- grep("ERV1|ERVK|ERVL",repclasses)
nomatchIndexes <- setdiff(seq(1,length(repclasses)),matchIndexes)
repclasses[nomatchIndexes] <- "Other"
repclasses[grep("ERV",repclasses)] <- "ERV"  
plotdata$Size <- rep(0.01,nrow(plotdata))
plotdata$Size <- as.numeric(plotdata$Size)
plotdata$y <- log2(rowSums(plotdata[,c(4,5)])/rowSums(plotdata[,c(2,3)]))
plotdata$x <- log10(rowSums(plotdata[,c(2,3,4,5)]))
plotdata$group <- factor(repclasses,levels=c("ERV","Other"))

plotdata <- plotdata[which(!is.na(plotdata$y) & !is.infinite(plotdata$y)),]
matchindexes <- match(targetedGenes,rownames(plotdata))

png("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile_erv/CT_vs_Rote_erv.png",width=12,height=10,units="cm",res=600)
p <- ggplot(plotdata,aes(x=x,y=y,color=group)) + geom_point(size=0.2) + geom_point(data=plotdata[matchindexes,], aes(x=x, y=y), colour="black", size=2) + scale_color_manual(values =c('red','grey'))
p <- p + xlab(expression(log[10]("TPM"))) + ylab(expression(log[2]("FC (CX/Control)")))
p <- p + scale_y_continuous(limits = c(-5,5))
p <- p + geom_abline(intercept=0,slope=0,colour="grey")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()
