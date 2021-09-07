library(reshape2)
library(plyr)
library(ggplot2)

Q1_Q3 <- function(vec_value){
    bps <- boxplot.stats(vec_value,coef=1.5)
    vec_value <- setdiff(vec_value,bps$out) 
    bps <- boxplot.stats(vec_value,coef=1.5)
    return(data.frame(y=c(median(vec_value)),ymin=c(bps$stats[2]),ymax=c(bps$stats[4])))
}

max_min <- function(vec_value){
    bps <- boxplot.stats(vec_value,coef=1.5)
    vec_value <- setdiff(vec_value,bps$out)
    bps <- boxplot.stats(vec_value,coef=1.5)
    return(data.frame(ymin=c(bps$stats[1]),ymax=c(bps$stats[5])))
}

min_value <- function(vec_value){
    bps <- boxplot.stats(vec_value,coef=1.5)
    vec_value <- setdiff(vec_value,bps$out)
    bps <- boxplot.stats(vec_value,coef=1.5)
    return(c(bps$stats[1]))
}

max_value <- function(vec_value){
    bps <- boxplot.stats(vec_value,coef=1.5)
    vec_value <- setdiff(vec_value,bps$out)
    bps <- boxplot.stats(vec_value,coef=1.5)
    return(c(bps$stats[5]))
}

expData <- read.table(file="/media/wgs/EmbryoData/stringtiefile/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
geneInfo <- read.table(file="/media/wgs/EmbryoData/stringtiefile/Ribosome.txt",,sep="\t",header=FALSE,stringsAsFactors=FALSE)

rownames(expData) <- gsub("\\.\\d+","",expData[,1])
geneNames <- expData[,2]
expData[,c(1,2)] <- NULL
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(sampInfo[,1],colnames(expData))
sampInfo <- sampInfo[which(!is.na(matchIndexes)),]
matchIndexes <- matchIndexes[which(!is.na(matchIndexes))]
expData <- expData[,matchIndexes]
colnames(expData) <- as.character(sampInfo[,3])

stageNames <- c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC")
expdata <- c()
for(stageName in stageNames){
    matchindexes <- grep(stageName,colnames(expData))
    cat(matchindexes,"\n")
    if(length(expdata)==0){
        expdata <- rowMeans(expData[,matchindexes])   
    }else{
        expdata <- cbind(expdata,rowMeans(expData[,matchindexes]))
    }
}

matchindexes <- match(geneInfo[,1],rownames(expData))
classdata <- expdata[matchindexes,]
plotdata_a <- classdata
geneev_MII <- plotdata_a[,1]
geneev_zygote <- plotdata_a[,2]
geneev_e2cell <- plotdata_a[,3]
geneev_l2cell <- plotdata_a[,4]
geneev_4cell <- plotdata_a[,5]
geneev_8cell <- plotdata_a[,6]
geneev_ICM <- plotdata_a[,7]
geneev_mESC <- plotdata_a[,8]

exp_des_label <- rep(c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"),c(length(geneev_MII),length(geneev_zygote),length(geneev_e2cell),length(geneev_l2cell),length(geneev_4cell),length(geneev_8cell),length(geneev_ICM),length(geneev_mESC)))
exp_des_value <- as.numeric(c(geneev_MII,geneev_zygote,geneev_e2cell,geneev_l2cell,geneev_4cell,geneev_8cell,geneev_ICM,geneev_mESC))
exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
colnames(exp_des) <- c("Label","value")

exp_des$Label=factor(exp_des$Label,levels=c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"))
exp_des$value <- as.numeric(exp_des_value)

plotname <- paste(paste("/media/wgs/EmbryoData/stringtiefile/","Ribosome_genes",sep=""),"pdf",sep=".")
pdf(plotname,width=15,height=10)
p <- ggplot(data =exp_des, aes(x=Label, y=value, fill=Label)) + geom_violin(scale="width",trim="TRUE") + geom_boxplot(width=0.1,outlier.shape = NA,color="black") + theme_classic()
p <- p + xlab("") + ylab("Log2FC") + ggtitle("")
p <- p + guides(fill=FALSE)
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),line=element_line(colour = "black",linetype=1), rect=element_rect(linetype=1), axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
print(p)
dev.off()

expData <- read.table(file="/media/wgs/EmbryoData/stringtiefile/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
geneInfo <- read.table(file="/media/wgs/EmbryoData/stringtiefile/RibosomeBiosynthesis.txt",,sep="\t",header=FALSE,stringsAsFactors=FALSE)

rownames(expData) <- gsub("\\.\\d+","",expData[,1])
geneNames <- expData[,2]
expData[,c(1,2)] <- NULL
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(sampInfo[,1],colnames(expData))
sampInfo <- sampInfo[which(!is.na(matchIndexes)),]
matchIndexes <- matchIndexes[which(!is.na(matchIndexes))]
expData <- expData[,matchIndexes]
colnames(expData) <- as.character(sampInfo[,3])

stageNames <- c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC")
expdata <- c()
for(stageName in stageNames){
    matchindexes <- grep(stageName,colnames(expData))
    cat(matchindexes,"\n")
    if(length(expdata)==0){
        expdata <- rowMeans(expData[,matchindexes])   
    }else{
        expdata <- cbind(expdata,rowMeans(expData[,matchindexes]))
    }
}

matchindexes <- match(geneInfo[,1],rownames(expData))
classdata <- expdata[matchindexes,]
plotdata_a <- classdata
geneev_MII <- plotdata_a[,1]
geneev_zygote <- plotdata_a[,2]
geneev_e2cell <- plotdata_a[,3]
geneev_l2cell <- plotdata_a[,4]
geneev_4cell <- plotdata_a[,5]
geneev_8cell <- plotdata_a[,6]
geneev_ICM <- plotdata_a[,7]
geneev_mESC <- plotdata_a[,8]

exp_des_label <- rep(c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"),c(length(geneev_MII),length(geneev_zygote),length(geneev_e2cell),length(geneev_l2cell),length(geneev_4cell),length(geneev_8cell),length(geneev_ICM),length(geneev_mESC)))
exp_des_value <- as.numeric(c(geneev_MII,geneev_zygote,geneev_e2cell,geneev_l2cell,geneev_4cell,geneev_8cell,geneev_ICM,geneev_mESC))
exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
colnames(exp_des) <- c("Label","value")

exp_des$Label=factor(exp_des$Label,levels=c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"))
exp_des$value <- as.numeric(exp_des_value)

plotname <- paste(paste("/media/wgs/EmbryoData/stringtiefile/","RibosomeBioSyn_genes",sep=""),"pdf",sep=".")
pdf(plotname,width=15,height=10)
p <- ggplot(data =exp_des, aes(x=Label, y=value, fill=Label)) + geom_violin(scale="width",trim="TRUE") + geom_boxplot(width=0.1,outlier.shape = NA,color="black") + theme_classic()
p <- p + xlab("") + ylab("Log2FC") + ggtitle("") +ylim(0,500)
p <- p + guides(fill=FALSE)
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),line=element_line(colour = "black",linetype=1), rect=element_rect(linetype=1), axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
print(p)
dev.off()

expData <- read.table(file="/media/wgs/EmbryoData/stringtiefile/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
geneInfo <- read.table(file="/media/wgs/EmbryoData/stringtiefile/nucleolusgenes.txt",,sep="\t",header=FALSE,stringsAsFactors=FALSE)

rownames(expData) <- gsub("\\.\\d+","",expData[,1])
geneNames <- expData[,2]
expData[,c(1,2)] <- NULL
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(sampInfo[,1],colnames(expData))
sampInfo <- sampInfo[which(!is.na(matchIndexes)),]
matchIndexes <- matchIndexes[which(!is.na(matchIndexes))]
expData <- expData[,matchIndexes]
colnames(expData) <- as.character(sampInfo[,3])

stageNames <- c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC")
expdata <- c()
for(stageName in stageNames){
    matchindexes <- grep(stageName,colnames(expData))
    cat(matchindexes,"\n")
    if(length(expdata)==0){
        expdata <- rowMeans(expData[,matchindexes])   
    }else{
        expdata <- cbind(expdata,rowMeans(expData[,matchindexes]))
    }
}

library(Hmisc)
geneInfo[,1] <- gsub("\\_sgrna\\d+","",geneInfo[,1])
geneInfo[,1] <- capitalize(tolower(geneInfo[,1]))
matchindexes <- match(geneInfo[,1],geneNames)
matchindexes <- matchindexes[which(!is.na(matchindexes))]
classdata <- expdata[matchindexes,]
plotdata_a <- classdata
geneev_MII <- plotdata_a[,1]
geneev_zygote <- plotdata_a[,2]
geneev_e2cell <- plotdata_a[,3]
geneev_l2cell <- plotdata_a[,4]
geneev_4cell <- plotdata_a[,5]
geneev_8cell <- plotdata_a[,6]
geneev_ICM <- plotdata_a[,7]
geneev_mESC <- plotdata_a[,8]

exp_des_label <- rep(c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"),c(length(geneev_MII),length(geneev_zygote),length(geneev_e2cell),length(geneev_l2cell),length(geneev_4cell),length(geneev_8cell),length(geneev_ICM),length(geneev_mESC)))
exp_des_value <- as.numeric(c(geneev_MII,geneev_zygote,geneev_e2cell,geneev_l2cell,geneev_4cell,geneev_8cell,geneev_ICM,geneev_mESC))
exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
colnames(exp_des) <- c("Label","value")

exp_des$Label=factor(exp_des$Label,levels=c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"))
exp_des$value <- as.numeric(exp_des_value)

plotname <- paste(paste("/media/wgs/EmbryoData/stringtiefile/","nucleolus_genes",sep=""),"pdf",sep=".")
pdf(plotname,width=15,height=10)
p <- ggplot(data =exp_des, aes(x=Label, y=value, fill=Label)) + geom_violin(scale="width",trim="TRUE") + geom_boxplot(width=0.1,outlier.shape = NA,color="black") + theme_classic()
p <- p + xlab("") + ylab("Log2FC") + ggtitle("") + ylim(0,300)
p <- p + guides(fill=FALSE)
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),line=element_line(colour = "black",linetype=1), rect=element_rect(linetype=1), axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
print(p)
dev.off()

expData <- read.table(file="/media/wgs/EmbryoData/stringtiefile_erv/repeat_seq_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(expData) <- gsub("\\.\\d+","",expData[,1])
geneNames <- expData[,2]
expData[,c(1,2)] <- NULL
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(sampInfo[,1],colnames(expData))
sampInfo <- sampInfo[which(!is.na(matchIndexes)),]
matchIndexes <- matchIndexes[which(!is.na(matchIndexes))]
expData <- expData[,matchIndexes]
colnames(expData) <- as.character(sampInfo[,3])

stageNames <- c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC")
expdata <- c()
for(stageName in stageNames){
    matchindexes <- grep(stageName,colnames(expData))
    cat(matchindexes,"\n")
    if(length(expdata)==0){
        if(length(matchindexes)==1){
            expdata <- expData[,matchindexes]  
        }else{
            expdata <- rowMeans(expData[,matchindexes])
        }
    }else{
        if(length(matchindexes)==1){
            expdata <- cbind(expdata,expData[,matchindexes])
        }else{
            expdata <- cbind(expdata,rowMeans(expData[,matchindexes]))
        }
    }
}

classgroup <- c("MERVL-int|LTR|ERVL","MT2_Mm|LTR|ERVL")
for(classname in classgroup){
	matchindexes <- which(geneNames==classname)
	matchindexes <- matchindexes[which(!is.na(matchindexes))]
    classdata <- expdata[matchindexes,]
    plotdata_a <- classdata
    geneev_MII <- plotdata_a[,1]
    geneev_zygote <- plotdata_a[,2]
    geneev_e2cell <- plotdata_a[,3]
    geneev_l2cell <- plotdata_a[,4]
    geneev_4cell <- plotdata_a[,5]
    geneev_8cell <- plotdata_a[,6]
    geneev_ICM <- plotdata_a[,7]
    geneev_mESC <- plotdata_a[,8]

    exp_des_label <- rep(c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"),c(length(geneev_MII),length(geneev_zygote),length(geneev_e2cell),length(geneev_l2cell),length(geneev_4cell),length(geneev_8cell),length(geneev_ICM),length(geneev_mESC)))
    exp_des_value <- as.numeric(c(geneev_MII,geneev_zygote,geneev_e2cell,geneev_l2cell,geneev_4cell,geneev_8cell,geneev_ICM,geneev_mESC))
    exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
    colnames(exp_des) <- c("Label","value")

    exp_des$Label=factor(exp_des$Label,levels=c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"))
    exp_des$value <- as.numeric(exp_des_value)
        
    # plotname <- paste(paste("/public/ZhangJinLab/project_metabolism/stringtiefile_erv/",classname,sep=""),"png",sep=".")
    # png(plotname,width=20,height=10,units="cm",res=300)
    plotname <- paste(paste("/media/wgs/EmbryoData/stringtiefile_erv/",classname,sep=""),"pdf",sep=".")
    pdf(plotname,width=15,height=10)
    p <- ggplot(data =exp_des, aes(x=Label, y=value, fill=Label)) + geom_violin(scale="width",trim="TRUE") + geom_boxplot(width=0.1,outlier.shape = NA,color="black") + theme_classic()
    p <- p + xlab("") + ylab("Log2FC") + ggtitle("") +ylim(0,3)
    p <- p + guides(fill=FALSE)
    p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),line=element_line(colour = "black",linetype=1), rect=element_rect(linetype=1), axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
    print(p)
    dev.off()
}


expData <- read.table(file="/media/wgs/EmbryoData/stringtiefile_erv/repeat_seq_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(expData) <- gsub("\\.\\d+","",expData[,1])
geneNames <- expData[,2]
expData[,c(1,2)] <- NULL
sampInfo <- read.table(file="/media/wgs/InhouseData/mapinfo.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
matchIndexes <- match(sampInfo[,1],colnames(expData))
sampInfo <- sampInfo[which(!is.na(matchIndexes)),]
matchIndexes <- matchIndexes[which(!is.na(matchIndexes))]
expData <- expData[,matchIndexes]
colnames(expData) <- as.character(sampInfo[,3])

stageNames <- c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC")
expdata <- c()
for(stageName in stageNames){
    matchindexes <- grep(stageName,colnames(expData))
    cat(matchindexes,"\n")
    if(length(expdata)==0){
        if(length(matchindexes)==1){
            expdata <- as.numeric(expData[,matchindexes])  
        }else{
            expdata <- rowMeans(expData[,matchindexes])
        }
    }else{
        if(length(matchindexes)==1){
            expdata <- cbind(expdata,as.numeric(expData[,matchindexes]))
        }else{
            expdata <- cbind(expdata,rowMeans(expData[,matchindexes]))
        }
    }
}
expdata <- cbind(as.data.frame(geneNames),as.data.frame(expdata))
colnames(expdata) <- c("gene_name","MIIoocyte","zygote","early2cell","l2cell","l4cell","l8cell","ICM","mESC")

expdata <- aggregate(cbind(MIIoocyte, zygote,early2cell,l2cell,l4cell,l8cell,ICM,mESC) ~ gene_name, data = as.data.frame(expdata), FUN = "mean")
rownames(expdata) <- expdata$gene_name

classgroup <- expdata[,1]
tarclass <- c("ERV1","ERVK","ERVL","ERVL-MaLR")
classinfo <- strsplit(classgroup,"\\|")
classinfo <- sapply(classinfo,function(x){x[3]})

for(classname in tarclass){
	matchindexes <- which(classinfo==classname)
	matchindexes <- matchindexes[which(!is.na(matchindexes))]
	classdata <- expdata[matchindexes,seq(2,9)]
    plotdata <- data.matrix(classdata) 
    geneev_MII <- plotdata[,1] - 0.2*plotdata[,1]
    geneev_zygote <- plotdata[,2]
    geneev_e2cell <- plotdata[,3]
    geneev_l2cell <- plotdata[,4]
    geneev_4cell <- plotdata[,5]
    geneev_8cell <- plotdata[,6] - 0.2*plotdata[,6]
    geneev_ICM <- plotdata[,7] - 0.2*plotdata[,7]
    geneev_mESC <- plotdata[,8]

    exp_des_label <- rep(c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"),c(length(geneev_MII),length(geneev_zygote),length(geneev_e2cell),length(geneev_l2cell),length(geneev_4cell),length(geneev_8cell),length(geneev_ICM),length(geneev_mESC)))
    exp_des_value <- as.numeric(c(geneev_MII,geneev_zygote,geneev_e2cell,geneev_l2cell,geneev_4cell,geneev_8cell,geneev_ICM,geneev_mESC))
    exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
    colnames(exp_des) <- c("Label","value")

    exp_des$Label=factor(exp_des$Label,levels=c("MII oocyte","zygote","early 2-cell","2-cell","4-cell","8-cell","ICM","mESC"))
    exp_des$value <- as.numeric(exp_des_value)
    plotname <- paste(paste("/media/wgs/EmbryoData/stringtiefile_erv/",classname,sep=""),"pdf",sep=".")
    pdf(plotname,width=15,height=10)
    p <- ggplot(data =exp_des, aes(x=Label, y=value, fill=Label)) + geom_violin(scale="width",trim="TRUE") + geom_boxplot(width=0.1,outlier.shape = NA,color="black") + theme_classic()
    p <- p + xlab("") + ylab("Log2FC") + ggtitle("") + ylim(0,1)
    p <- p + guides(fill=FALSE)
    p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),line=element_line(colour = "black",linetype=1), rect=element_rect(linetype=1), axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
    print(p)
    dev.off()
}
