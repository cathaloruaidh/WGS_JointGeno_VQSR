library(data.table)
library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)

inputName = args[1]
outputName = args[2]

write("Plotting Annotation Densities. \n", stderr())

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

dataALL <- fread(inputName, header = T)
write("Input file read. \n", stderr())

dataSNP <- dataALL[dataALL$TYPE == "SNP", ]
dataINDEL <- dataALL[dataALL$TYPE == "INDEL"]
dataMIXED <- dataALL[dataALL$TYPE == "MIXED"]

dataALL$FILTER_BIN <- ifelse(dataALL$FILTER == "PASS", "PASS", "FAIL")
dataSNP$FILTER_BIN <- ifelse(grepl("SNP", dataSNP$FILTER, fixed = T), "FAIL", "PASS")
dataINDEL$FILTER_BIN <- ifelse(grepl("INDEL", dataINDEL$FILTER, fixed = T), "FAIL", "PASS")


SNPFilterQD <- 2.0
SNPFilterMQ <- 40.0
SNPFilterMQRankSum <- -12.5 
SNPFilterReadPosRankSum <- -8.0
SNPFilterFS <- 60.0
SNPFilterSOR <- 3.0
SNPFilterInbreedingCoeff <- -0.8
SNPFilterDP <- mean(dataALL$DP) + 5*sd(dataALL$DP)

IndelFilterQD <- 2.0
IndelFilterReadPosRankSum <- -20.0
IndelFilterFS <- 200.0
IndelFilterSOR <- 10.0
IndelFilterInbreedingCoeff <- -0.8
IndelFilterDP <- mean(dataALL$DP) + 5*sd(dataALL$DP)

densityIntensity <- 0.3
quantileMax <- 0.999
quantileMin <- 0.001




### QD Plots 

#xMaxQDALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$QD, quantileMax, na.rm = T)
#xMaxQDSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$QD, quantileMax, na.rm = T)
#xMaxQDINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$QD, quantileMax, na.rm = T)
#xMaxQDMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$QD, quantileMax, na.rm = T)
#xMinQDALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$QD, quantileMin, na.rm = T)
#xMinQDSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$QD, quantileMin, na.rm = T)
#xMinQDINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$QD, quantileMin, na.rm = T)
#xMinQDMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$QD, quantileMin, na.rm = T)

pdf(paste0(outputName, ".recal.table.QD.pdf"))

pALL <- ggplot(dataALL, aes(x = QD, fill = FILTER_BIN)) + 
  ggtitle("All Variants") + 
  guides(fill=guide_legend(title="Filter")) + 
  geom_density(alpha=densityIntensity) 
pSNP <- ggplot(dataSNP, aes(x = QD, fill = FILTER_BIN)) + 
  ggtitle("SNPs") + 
  guides(fill=guide_legend(title="Filter")) + 
  geom_density(alpha=densityIntensity) + 
  geom_vline(xintercept = SNPFilterQD)
pINDEL <- ggplot(dataINDEL, aes(x = QD, fill = FILTER_BIN)) + 
  ggtitle("Indels") + 
  guides(fill=guide_legend(title="Filter")) + 
  geom_density(alpha=densityIntensity) + 
  geom_vline(xintercept = IndelFilterQD)
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardQD", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = QD, fill = FILTER_BIN)) + 
  ggtitle("Mixed Variants") + 
  guides(fill=guide_legend(title="Filter")) + 
  geom_density(alpha=densityIntensity) + 
  geom_vline(xintercept = IndelFilterQD)

multiplot(pALL, pINDEL, pSNP, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)

dev.off()
write("Finished QD \n", stderr())


### DP Plots

pdf(paste0(outputName, ".recal.table.DP.pdf"))

xMaxDPALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$DP, quantileMax, na.rm = T)
xMaxDPSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$DP, quantileMax, na.rm = T)
xMaxDPINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$DP, quantileMax, na.rm = T)
xMaxDPMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$DP, quantileMax, na.rm = T)
xMinDPALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$DP, quantileMin, na.rm = T)
xMinDPSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$DP, quantileMin, na.rm = T)
xMinDPINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$DP, quantileMin, na.rm = T)
xMinDPMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$DP, quantileMin, na.rm = T)

dataALL$FILTER_BIN <- ifelse(grepl("HighDP", dataALL$FILTER, fixed = T), "FAIL", "PASS")
pALL <- ggplot(dataALL, aes(x = DP, fill = FILTER_BIN)) + ggtitle("All Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinDPALL, xMaxDPALL)) + geom_density(alpha=densityIntensity) 
dataSNP$FILTER_BIN <- ifelse(grepl("HighDP", dataSNP$FILTER, fixed = T), "FAIL", "PASS")
pSNP <- ggplot(dataSNP, aes(x = DP, fill = FILTER_BIN)) + ggtitle("SNPs") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinDPSNP, xMaxDPSNP)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = SNPFilterDP)
dataINDEL$FILTER_BIN <- ifelse(grepl("HighDP", dataINDEL$FILTER, fixed = T), "FAIL", "PASS")
pINDEL <- ggplot(dataINDEL, aes(x = DP, fill = FILTER_BIN)) + ggtitle("Indels") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinDPINDEL, xMaxDPINDEL)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = IndelFilterDP)
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardDP", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = DP, fill = FILTER_BIN)) + ggtitle("Mixed Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinDPMIXED, xMaxDPMIXED)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = IndelFilterDP)

multiplot(pALL, pSNP, pINDEL, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)

dev.off()
write("Finished DP\n", stderr())


dataALL$FILTER_BIN <- ifelse(dataALL$FILTER == "PASS", "PASS", "FAIL")
dataSNP$FILTER_BIN <- ifelse(grepl("SNP", dataSNP$FILTER, fixed = T), "FAIL", "PASS")
dataINDEL$FILTER_BIN <- ifelse(grepl("INDEL", dataINDEL$FILTER, fixed = T), "FAIL", "PASS")



### FS Plots

#xMaxFSALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$FS, quantileMax, na.rm = T)
#xMaxFSSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$FS, quantileMax, na.rm = T)
#xMaxFSINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$FS, quantileMax, na.rm = T)
#xMaxFSMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$FS, quantileMax, na.rm = T)
#xMinFSALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$FS, quantileMin, na.rm = T)
#xMinFSSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$FS, quantileMin, na.rm = T)
#xMinFSINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$FS, quantileMin, na.rm = T)
#xMinFSMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$FS, quantileMin, na.rm = T)

pdf(paste0(outputName, ".recal.table.FS.pdf"))

pALL <- ggplot(dataALL, aes(x = FS, fill = FILTER_BIN)) + ggtitle("All Variants") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") 
pSNP <- ggplot(dataSNP, aes(x = FS, fill = FILTER_BIN)) + ggtitle("SNPs") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") + geom_vline(xintercept = SNPFilterFS)
pINDEL <- ggplot(dataINDEL, aes(x = FS, fill = FILTER_BIN)) + ggtitle("Indels") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") + geom_vline(xintercept = IndelFilterFS)
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardFS", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = FS, fill = FILTER_BIN)) + ggtitle("Mixed Variants") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") + geom_vline(xintercept = IndelFilterFS)

multiplot(pALL, pSNP, pINDEL, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)

dev.off()

write("Finished FS\n", stderr())

### SOR Plots

#xMaxSORALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$SOR, quantileMax, na.rm = T)
#xMaxSORSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$SOR, quantileMax, na.rm = T)
#xMaxSORINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$SOR, quantileMax, na.rm = T)
#xMaxSORMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$SOR, quantileMax, na.rm = T)
#xMinSORALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$SOR, quantileMin, na.rm = T)
#xMinSORSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$SOR, quantileMin, na.rm = T)
#xMinSORINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$SOR, quantileMin, na.rm = T)
#xMinSORMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$SOR, quantileMin, na.rm = T)

pdf(paste0(outputName, ".recal.table.SOR.pdf"))

pALL <- ggplot(dataALL, aes(x = SOR, fill = FILTER_BIN)) + ggtitle("All Variants") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") 
pSNP <- ggplot(dataSNP, aes(x = SOR, fill = FILTER_BIN)) + ggtitle("SNPs") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") + geom_vline(xintercept = SNPFilterSOR)
pINDEL <- ggplot(dataINDEL, aes(x = SOR, fill = FILTER_BIN)) + ggtitle("Indels") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") + geom_vline(xintercept = IndelFilterSOR)
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardSOR", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = SOR, fill = FILTER_BIN)) + ggtitle("Mixed Variants") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + scale_x_continuous(trans = "log10") + geom_vline(xintercept = IndelFilterSOR)

multiplot(pALL, pSNP, pINDEL, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)

dev.off()

write("Finished SOR\n", stderr())

### MQ Plots

pdf(paste0(outputName, ".recal.table.MQ.pdf"))

xMaxMQALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$MQ, quantileMax, na.rm = T)
xMaxMQSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$MQ, quantileMax, na.rm = T)
xMaxMQINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$MQ, quantileMax, na.rm = T)
xMaxMQMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$MQ, quantileMax, na.rm = T)
xMinMQALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$MQ, quantileMin, na.rm = T)
xMinMQSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$MQ, quantileMin, na.rm = T)
xMinMQINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$MQ, quantileMin, na.rm = T)
xMinMQMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$MQ, quantileMin, na.rm = T)

pALL <- ggplot(dataALL, aes(x = MQ, fill = FILTER_BIN)) + ggtitle("All Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQALL, xMaxMQALL)) + geom_density(alpha=densityIntensity) 
pSNP <- ggplot(dataSNP, aes(x = MQ, fill = FILTER_BIN)) + ggtitle("SNPs") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQSNP, xMaxMQSNP)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = SNPFilterMQ)
pINDEL <- ggplot(dataINDEL, aes(x = MQ, fill = FILTER_BIN)) + ggtitle("Indels") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQINDEL, xMaxMQINDEL)) + geom_density(alpha=densityIntensity) 
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardMQ", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = MQ, fill = FILTER_BIN)) + ggtitle("Mixed Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQMIXED, xMaxMQMIXED)) + geom_density(alpha=densityIntensity) 

multiplot(pALL, pSNP, pINDEL, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)

dev.off()

write("Finished MQ\n", stderr())

### MQRankSum Plots

pdf(paste0(outputName, ".recal.table.MQRS.pdf"))

xMaxMQRankSumALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$MQRankSum, quantileMax, na.rm = T)
xMaxMQRankSumSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$MQRankSum, quantileMax, na.rm = T)
xMaxMQRankSumINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$MQRankSum, quantileMax, na.rm = T)
xMaxMQRankSumMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$MQRankSum, quantileMax, na.rm = T)
xMinMQRankSumALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$MQRankSum, quantileMin, na.rm = T)
xMinMQRankSumSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$MQRankSum, quantileMin, na.rm = T)
xMinMQRankSumINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$MQRankSum, quantileMin, na.rm = T)
xMinMQRankSumMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$MQRankSum, quantileMin, na.rm = T)

pALL <- ggplot(dataALL, aes(x = MQRankSum, fill = FILTER_BIN)) + ggtitle("All Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQRankSumALL, xMaxMQRankSumALL)) + geom_density(alpha=densityIntensity) 
pSNP <- ggplot(dataSNP, aes(x = MQRankSum, fill = FILTER_BIN)) + ggtitle("SNPs") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQRankSumSNP, xMaxMQRankSumSNP)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = SNPFilterMQRankSum)
pINDEL <- ggplot(dataINDEL, aes(x = MQRankSum, fill = FILTER_BIN)) + ggtitle("Indels") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQRankSumINDEL, xMaxMQRankSumINDEL)) + geom_density(alpha=densityIntensity) 
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardMQRankSum", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = MQRankSum, fill = FILTER_BIN)) + ggtitle("Mixed Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinMQRankSumMIXED, xMaxMQRankSumMIXED)) + geom_density(alpha=densityIntensity) 

multiplot(pALL, pSNP, pINDEL, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)

dev.off()

write("Finished MQRankSum\n", stderr())


### ReadPosRankSum Plots

pdf(paste0(outputName, ".recal.table.RPRS.pdf"))

xMaxReadPosRankSumALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$ReadPosRankSum, quantileMax, na.rm = T)
xMaxReadPosRankSumSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$ReadPosRankSum, quantileMax, na.rm = T)
xMaxReadPosRankSumINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$ReadPosRankSum, quantileMax, na.rm = T)
xMaxReadPosRankSumMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$ReadPosRankSum, quantileMax, na.rm = T)
xMinReadPosRankSumALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$ReadPosRankSum, quantileMin, na.rm = T)
xMinReadPosRankSumSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$ReadPosRankSum, quantileMin, na.rm = T)
xMinReadPosRankSumINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$ReadPosRankSum, quantileMin, na.rm = T)
xMinReadPosRankSumMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$ReadPosRankSum, quantileMin, na.rm = T)

pALL <- ggplot(dataALL, aes(x = ReadPosRankSum, fill = FILTER_BIN)) + ggtitle("All Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinReadPosRankSumALL, xMaxReadPosRankSumALL)) + geom_density(alpha=densityIntensity) 
pSNP <- ggplot(dataSNP, aes(x = ReadPosRankSum, fill = FILTER_BIN)) + ggtitle("SNPs") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinReadPosRankSumSNP, xMaxReadPosRankSumSNP)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = SNPFilterReadPosRankSum)
pINDEL <- ggplot(dataINDEL, aes(x = ReadPosRankSum, fill = FILTER_BIN)) + ggtitle("Indels") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinReadPosRankSumINDEL, xMaxReadPosRankSumINDEL)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = IndelFilterReadPosRankSum)
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardReadPosRankSum", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = ReadPosRankSum, fill = FILTER_BIN)) + ggtitle("Mixed Variants") + guides(fill=guide_legend(title="Filter")) + scale_x_continuous(limits = c(xMinReadPosRankSumMIXED, xMaxReadPosRankSumMIXED)) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = IndelFilterReadPosRankSum)

multiplot(pALL, pSNP, pINDEL, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)

dev.off()

write("Finished ReadPosRankSum\n", stderr())

### InbreedingCoefficient Plots

#xMaxInbreedingCoeffALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$InbreedingCoeff, quantileMax, na.rm = T)
#xMaxInbreedingCoeffSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$InbreedingCoeff, quantileMax, na.rm = T)
#xMaxInbreedingCoeffINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$InbreedingCoeff, quantileMax, na.rm = T)
#xMaxInbreedingCoeffMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$InbreedingCoeff, quantileMax, na.rm = T)
#xMinInbreedingCoeffALL <- quantile(dataALL[dataALL$FILTER == "PASS",]$InbreedingCoeff, quantileMin, na.rm = T)
#xMinInbreedingCoeffSNP <- quantile(dataSNP[dataSNP$FILTER == "PASS",]$InbreedingCoeff, quantileMin, na.rm = T)
#xMinInbreedingCoeffINDEL <- quantile(dataINDEL[dataINDEL$FILTER == "PASS",]$InbreedingCoeff, quantileMin, na.rm = T)
#xMinInbreedingCoeffMIXED <- quantile(dataMIXED[dataMIXED$FILTER == "PASS",]$InbreedingCoeff, quantileMin, na.rm = T)

pdf(paste0(outputName, ".recal.table.IC.pdf"))

pALL <- ggplot(dataALL, aes(x = InbreedingCoeff, fill = FILTER_BIN)) + ggtitle("All Variants") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) 
pSNP <- ggplot(dataSNP, aes(x = InbreedingCoeff, fill = FILTER_BIN)) + ggtitle("SNPs") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = SNPFilterInbreedingCoeff)
pINDEL <- ggplot(dataINDEL, aes(x = InbreedingCoeff, fill = FILTER_BIN)) + ggtitle("Indels") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = IndelFilterInbreedingCoeff)
dataMIXED$FILTER_BIN <- ifelse(grepl("OtherHardInbreedingCoeff", dataMIXED$FILTER, fixed = T), "FAIL", "PASS")
pMIXED <- ggplot(dataMIXED, aes(x = InbreedingCoeff, fill = FILTER_BIN)) + ggtitle("Mixed Variants") + guides(fill=guide_legend(title="Filter")) + geom_density(alpha=densityIntensity) + geom_vline(xintercept = IndelFilterInbreedingCoeff)

multiplot(pALL, pSNP, pINDEL, pMIXED, cols = 2)

rm(pALL)
rm(pSNP)
rm(pINDEL)
rm(pMIXED)


dev.off()

write("Finished InbreedingCoeff\n", stderr())
