library(data.table)
library(ggplot2)
library(reshape2)


args <- commandArgs(trailingOnly = TRUE)

inputNameSNP = args[1]
inputNameINDEL = args[1]
outputName = args[3]



TranchesSNP <- fread(inputNameSNP, header = T)
TranchesINDEL <- fread(inputNameINDEL, header = T)

TranchesSNPBar <- TranchesSNP[ ,c(1,2,3)]
TranchesINDELBar <- TranchesINDEL[ ,c(1,2,3)]

TranchesSNPBar$numKnownNew <- TranchesSNPBar$numKnown
TranchesSNPBar$numNovelNew <- TranchesSNPBar$numNovel
TranchesINDELBar$numKnownNew <- TranchesINDELBar$numKnown
TranchesINDELBar$numNovelNew <- TranchesINDELBar$numNovel

TranchesSNPBar$numKnownPrev <- 0
TranchesSNPBar$numNovelPrev <- 0
TranchesINDELBar$numKnownPrev <- 0
TranchesINDELBar$numNovelPrev <- 0

for(i in 2:nrow(TranchesSNP)){
  TranchesSNPBar[i,]$numKnownNew <- TranchesSNPBar[i,]$numKnown - TranchesSNPBar[i-1,]$numKnown
  TranchesSNPBar[i,]$numKnownPrev <- TranchesSNPBar[i-1,]$numKnown
  TranchesSNPBar[i,]$numNovelNew <- TranchesSNPBar[i,]$numNovel - TranchesSNPBar[i-1,]$numNovel
  TranchesSNPBar[i,]$numNovelPrev <- TranchesSNPBar[i-1,]$numNovel
  
  TranchesINDELBar[i,]$numKnownNew <- TranchesINDELBar[i,]$numKnown - TranchesINDELBar[i-1,]$numKnown
  TranchesINDELBar[i,]$numKnownPrev <- TranchesINDELBar[i-1,]$numKnown
  TranchesINDELBar[i,]$numNovelNew <- TranchesINDELBar[i,]$numNovel - TranchesINDELBar[i-1,]$numNovel
  TranchesINDELBar[i,]$numNovelPrev <- TranchesINDELBar[i-1,]$numNovel
}


TranchesSNPBar <- TranchesSNPBar[, c(1,5,7,4,6)]
TranchesINDELBar <- TranchesINDELBar[, c(1,5,7,4,6)]

TranchesSNPBarMelt <- melt(TranchesSNPBar, id.vars = "targetTruthSensitivity")
TranchesINDELBarMelt <- melt(TranchesINDELBar, id.vars = "targetTruthSensitivity")

TranchesSNPBarMelt$targetTruthSensitivity <- factor(TranchesSNPBarMelt$targetTruthSensitivity)
TranchesINDELBarMelt$targetTruthSensitivity <- factor(TranchesINDELBarMelt$targetTruthSensitivity)


pdf(paste0(outputName, ".tranches.SNP.hist.pdf"))
ggplot(TranchesSNPBarMelt, aes(x = targetTruthSensitivity, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(limits = rev(levels(TranchesSNPBarMelt$targetTruthSensitivity))) +
  coord_flip() + 
  scale_fill_hue(l=45)

dev.off()

write("SNP histogram complete. ", stderr())


pdf(paste0(outputName, ".tranches.INDEL.hist.pdf"))

ggplot(TranchesINDELBarMelt, aes(x = targetTruthSensitivity, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(limits = rev(levels(TranchesINDELBarMelt$targetTruthSensitivity))) + 
  coord_flip()

dev.off()

write("Indel histogram complete. ", stderr())

TranchesSNPLine <- TranchesSNP[ ,c(1,4,5)]

TranchesSNPLineMelt <- melt(TranchesSNPLine, id.vars = "targetTruthSensitivity")

pdf(paste0(outputName, ".tranches.SNP.line.pdf"))

ggplot(TranchesSNPLineMelt, aes(x = targetTruthSensitivity, y = value, group = variable, col = variable)) + 
  geom_point(size = 2) + 
  geom_line(col = "dimgray") + 
  theme(axis.text.x = element_text(angle = -90)) + 
  geom_vline(xintercept = 99.9, col="#7CAE00", alpha=1) + 
  geom_vline(xintercept = TranchesSNPLineMelt[which(TranchesSNPLineMelt$value == max(TranchesSNPLineMelt[TranchesSNPLineMelt$variable == "novelTiTv",]$value)),]$targetTruthSensitivity, col="#C77CFF", alpha=1) 

dev.off()
write("SNP Ti/Tv plot complete. ", stderr())
