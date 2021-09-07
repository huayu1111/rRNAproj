library(reshape2)
library(ggplot2)
library(RColorBrewer)
df <- read.table(file="/media/wgs/StackBarPlotData/GroupC.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
plotData <- data.frame(rows=factor(rownames(df)[row(df)],levels=c("dpc1.5","dpc2.5","dpc3.5","dpc4.5")), vars=factor(colnames(df)[col(df)],levels=c("Fragmented","X1C","X2C","X3C","X4C","X5C","X6C","X7C","X8C","Morula","BC")),values=c(as.matrix(df)))

stck <- ggplot(plotData,aes(x= rows, y=values, fill=vars)) + geom_bar(stat="identity", width=.7) + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
stck <- stck + scale_fill_manual(values=c("#A3A500",colorRampPalette(c("#2B5CD8","#EA1BE5"))(10)))
ggsave("./GroupC.pdf")


library(reshape2)
library(ggplot2)
library(RColorBrewer)
df <- read.table(file="/media/wgs/StackBarPlotData/GroupD.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE)
df <- df/15
plotData <- data.frame(rows=factor(rownames(df)[row(df)],levels=c("dpc1.5","dpc2.5","dpc3.5","dpc4.5")), vars=factor(colnames(df)[col(df)],levels=c("Fragmented","X1C","X2C","X3C","X4C","X5C","X6C","X7C","X8C","Morula","BC")),values=c(as.matrix(df)))

stck <- ggplot(plotData,aes(x= rows, y=values, fill=vars)) + geom_bar(stat="identity", width=.7) + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
stck <- stck + scale_fill_manual(values=c("#A3A500",colorRampPalette(c("#2B5CD8","#EA1BE5"))(10)))
ggsave("./GroupD.pdf")