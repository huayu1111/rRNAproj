library(reshape2)
library(plyr)
library(ggplot2)
library(pheatmap)

expdata <- read.delim(file="/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/gene_TPM.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(expdata) <- gsub("\\.\\d+","",expdata[,1])
geneNames <- expdata[,2]
expdata[,c(1,2)] <- NULL

all_exp_des_value <- log2(rowMeans(expdata[,c(7,8)])/rowMeans(expdata[,c(3,4)]))
all_exp_des_value <- all_exp_des_value[which(!is.na(all_exp_des_value) & !is.infinite(all_exp_des_value) & all_exp_des_value > -5 & all_exp_des_value < 5)]

cluAssignFocus <- read.table(file="/public/ZhangJinLab/project_erv/stringtiefile/clusterInfo_stageFocus.txt",header=F,stringsAsFactors=F)
matchindexes <- match(cluAssignFocus[,1],rownames(expdata))

cluAssign <- cluAssignFocus[,2]
plotdata <- c()
upclusterNames <- c(1,3)
downclusterNames <- c(2,5,4)
for(classname in upclusterNames){
	matchindexes <- which(cluAssign==classname)
	matchindexes <- match(cluAssignFocus[matchindexes,1],rownames(expdata))
    matchindexes <- matchindexes[which(!is.na(matchindexes))]
	classdata <- expdata[matchindexes,]
    exp_des_value <- log2(rowMeans(classdata[,c(7,8)])/rowMeans(classdata[,c(3,4)]))
    exp_des_value <- exp_des_value[which(!is.na(exp_des_value) & !is.infinite(exp_des_value) & exp_des_value > -5 & exp_des_value < 5)]

    exp_des_label <- rep(paste("C",classname,sep=""),length(exp_des_value))
    if(length(exp_des_value) > 5){
        res = wilcox.test(exp_des_value,all_exp_des_value,paired=FALSE)
        cat(classname,length(exp_des_value),"P-value (wilcox.test):",res$p.value,mean(exp_des_value),"\n")
        if(length(plotdata)==0){
            plotdata <- cbind(exp_des_value,exp_des_label)
        }else{
            plotdata <- rbind(plotdata,cbind(exp_des_value,exp_des_label))
        }
    }
}

for(classname in downclusterNames){
	matchindexes <- which(cluAssign==classname)
	matchindexes <- match(cluAssignFocus[matchindexes,1],rownames(expdata))
    matchindexes <- matchindexes[which(!is.na(matchindexes))]
	classdata <- expdata[matchindexes,]
    exp_des_value <- log2(rowMeans(classdata[,c(7,8)])/rowMeans(classdata[,c(3,4)]))
    exp_des_value <- exp_des_value[which(!is.na(exp_des_value) & !is.infinite(exp_des_value) & exp_des_value > -5 & exp_des_value < 5)]

    exp_des_label <- rep(paste("C",classname,sep=""),length(exp_des_value))
    if(length(exp_des_value) > 5){
        res = wilcox.test(exp_des_value,all_exp_des_value,paired=FALSE)
        cat(classname,length(exp_des_value),"P-value (wilcox.test):",res$p.value,mean(exp_des_value),"\n")
        if(length(plotdata)==0){
            plotdata <- cbind(exp_des_value,exp_des_label)
        }else{
            plotdata <- rbind(plotdata,cbind(exp_des_value,exp_des_label))
        }
    }
}

all_exp_des_label <- rep("All",length(all_exp_des_value))
plotdata <- rbind(cbind(all_exp_des_value,all_exp_des_label),plotdata)

colnames(plotdata) <- c("Value","Label")
plotdata <- as.data.frame(plotdata,stringsAsFactors=F)

clusterIDs <- paste("C",c(1,3,2,5,4),sep="")
plotdata$Label=factor(plotdata$Label,levels=c("All","C1","C3","C2","C5","C4"))
plotdata$Value <- as.numeric(plotdata$Value)

plotname <- paste("/public/ZhangJinLab/project_stress/StressesRNAseqData/stringtiefile/stagegene_boxplot_Rapa","pdf",sep=".")
pdf(plotname,width=15,height=10)
p <- ggplot(data = plotdata, aes(x=Label, y=Value, fill=Label)) + geom_violin(scale="width",trim="TRUE") + geom_boxplot(width=0.1,outlier.shape = NA,color="black") + theme_classic()
p <- p + xlab("") + ylab("Log2FC") + ggtitle("")
p <- p + guides(fill=FALSE)
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),line=element_line(colour = "black",linetype=1), rect=element_rect(linetype=1), axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
print(p)
dev.off()
