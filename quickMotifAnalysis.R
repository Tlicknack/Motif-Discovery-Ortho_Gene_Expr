# Quick analysis of a motifs
library(stringr)
library(ggplot2)

##### All motifs

allMotifs = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/Newest_MEME_Final_File-Background.csv", header=T)
dupKmer = allMotifs[which(duplicated(allMotifs$Kmer)),]
dupconsensus = allMotifs[which(duplicated(allMotifs$Consensus)),]

##### Motifs
aaaaaa = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_AAAAAA_withExpression.csv", header=T)
aaaaaa = aaaaaa[- (which(aaaaaa$species == "pmultimic")),]
aaaaaa_exprLevels = split(aaaaaa$FPKM, aaaaaa$species)
aaaaaa_dists = split(aaaaaa$distToStart, aaaaaa$species)

spNames = names(sapply(split(aaaaaa$FPKM, aaaaaa$species), summary))
spNames = spNames[-which(spNames == "pmultimic")]

cccccc = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_CCCCCC_withExpression.csv", header=T)
cccccc = cccccc[- (which(cccccc$species == "pmultimic")),]
cccccc_exprLevels = split(cccccc$FPKM, cccccc$species)
cccccc_dists = split(cccccc$distToStart, cccccc$species)

tttttt = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_TTTTTT_withExpression.csv", header=T)
tttttt = tttttt[- (which(tttttt$species == "pmultimic")),]
tttttt_exprLevels = split(tttttt$FPKM, tttttt$species)
tttttt_dists = split(tttttt$distToStart, tttttt$species)

gggggg = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_GGGGGG_withExpression.csv", header=T)
gggggg = gggggg[- (which(gggggg$species == "pmultimic")),]
gggggg_exprLevels = split(gggggg$FPKM, gggggg$species)
gggggg_dists = split(gggggg$distToStart, gggggg$species)


aaatcttt = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_AAATCTTT_withExpression.csv", header=T)
aaatcttt = aaatcttt[- (which(aaatcttt$species == "pmultimic")),]
aaatcttt_exprLevels = split(aaatcttt$FPKM, aaatcttt$species)
aaatcttt_dists = split(aaatcttt$distToStart, aaatcttt$species)

tcgtta = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_TCGTTA_withExpression.csv", header=T)
tcgtta = tcgtta[- (which(tcgtta$species == "pmultimic")),]
tcgtta_exprLevels = split(tcgtta$FPKM, tcgtta$species)
tcgtta_dists = split(tcgtta$distToStart, tcgtta$species)

tcccgc = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_TCCCGC_withExpression.csv", header=T)
tcccgc = tcccgc[- (which(tcccgc$species == "pmultimic")),]
tcccgc_exprLevels = split(tcccgc$FPKM, tcccgc$species)
tcccgc_dists = split(tcccgc$distToStart, tcccgc$species)

#tcagtt = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_TCAGTT_withExpression.csv", header=T)
#tcagtt_exprLevels = split(tcagtt$FPKM, tcagtt$species)
#tcagtt_dists = split(tcagtt$distToStart, tcagtt$species)

#ttaacg = read.csv("/media/tlicknac/Seagate Expansion Drive1/Motif-Discovery_Paraorthologs/motif_TTAACG_withExpression.csv", header=T)
#ttaacg_exprLevels = split(ttaacg$FPKM, ttaacg$species)
#ttaacg_dists = split(ttaacg$distToStart, ttaacg$species)

##### Plotting
#> summary(as.numeric(gcPtet))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2006  0.2572  0.2758  0.2743  0.2878  0.3781 

##### Plotting Motif Distributions
ggplot(data=aaaaaa) + geom_freqpoly(mapping = aes(x=distToStart, y=((..count..)/ sum(..count..)), group=species, color = species), binwidth=5) + 
  xlab("Distance to 5'UTR (nts)") + ylab("Proportion of Hits") + ggtitle("AAAAAA Upstream Distribution")
ggplot(data=tttttt) + geom_freqpoly(mapping = aes(x=distToStart, y=((..count..)/ sum(..count..)), group=species, color = species), binwidth=5) + 
  xlab("Distance to 5'UTR (nts)") + ylab("Proportion of Hits") + ggtitle("TTTTTT Upstream Distribution")
ggplot(data=cccccc) + geom_freqpoly(mapping = aes(x=distToStart, y=((..count..)/ sum(..count..)), group=species, color = species), binwidth=5) + 
  xlab("Distance to 5'UTR (nts)") + ylab("Proportion of Hits") + ggtitle("CCCCCC Upstream Distribution")
ggplot(data=gggggg) + geom_freqpoly(mapping = aes(x=distToStart, y=((..count..)/ sum(..count..)), group=species, color = species), binwidth=5) + 
  xlab("Distance to 5'UTR (nts)") + ylab("Proportion of Hits") + ggtitle("GGGGGG Upstream Distribution")

ggplot(data=aaatcttt) + geom_freqpoly(mapping = aes(x=distToStart, y=((..count..)/ sum(..count..)), group=species, color = species), binwidth=5) + 
  xlab("Distance to 5'UTR (nts)") + ylab("Proportion of Hits") + ggtitle("AAATCTTT Upstream Distribution")
ggplot(data=tcccgc) + geom_freqpoly(mapping = aes(x=distToStart, y=((..count..)/ sum(..count..)), group=species, color = species), binwidth=5) + 
  xlab("Distance to 5'UTR (nts)") + ylab("Proportion of Hits") + ggtitle("TCCCGC Upstream Distribution")
ggplot(data=tcgtta) + geom_freqpoly(mapping = aes(x=distToStart, y=((..count..)/ sum(..count..)), group=species, color = species), binwidth=5) + 
  xlab("Distance to 5'UTR (nts)") + ylab("Proportion of Hits") + ggtitle("TCGTTA Upstream Distribution")

##### Plotting Motif Distance vs FPKM
library(MASS)
ggplot(data=aaaaaa) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + xlab("Distance to 5' UTR") + ggtitle("AAAAAA Expression vs Distance")
ggplot(data=aaaaaa) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + stat_smooth(mapping = aes(x=distToStart, y=FPKM), method = "lm") + xlab("Distance to 5' UTR") + ggtitle("AAAAAA Expression vs Distance")
ggplot(data=aaaaaa) + geom_point(mapping = aes(x=distToStart, y=FPKM, color=species)) + facet_wrap(~ species) + xlab("Distance to 5' UTR") + ggtitle("AAAAAA Expression vs Distance")
summary(lm(aaaaaa$FPKM ~ aaaaaa$distToStart))[[8]]

ggplot(data=tttttt) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + xlab("Distance to 5' UTR") + ggtitle("TTTTTT Expression vs Distance")
ggplot(data=tttttt) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + stat_smooth(mapping = aes(x=distToStart, y=FPKM), method = "lm") + xlab("Distance to 5' UTR") + ggtitle("TTTTTT Expression vs Distance")
ggplot(data=tttttt) + geom_point(mapping = aes(x=distToStart, y=FPKM, color=species)) + facet_wrap(~ species) + xlab("Distance to 5' UTR") + ggtitle("TTTTTT Expression vs Distance")
summary(lm(tttttt$FPKM ~ tttttt$distToStart))[[8]]

ggplot(data=cccccc) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + xlab("Distance to 5' UTR") + ggtitle("CCCCCC Expression vs Distance")
ggplot(data=cccccc) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + stat_smooth(mapping = aes(x=distToStart, y=FPKM), method = "lm") + xlab("Distance to 5' UTR") + ggtitle("CCCCCC Expression vs Distance")
ggplot(data=cccccc) + geom_point(mapping = aes(x=distToStart, y=FPKM, color=species)) + facet_wrap(~ species) + xlab("Distance to 5' UTR") + ggtitle("CCCCCC Expression vs Distance")
summary(lm(cccccc$FPKM ~ cccccc$distToStart))[[8]]

ggplot(data=gggggg) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + xlab("Distance to 5' UTR") + ggtitle("GGGGGG Expression vs Distance")
ggplot(data=gggggg) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + stat_smooth(mapping = aes(x=distToStart, y=FPKM), method = "lm") + xlab("Distance to 5' UTR") + ggtitle("GGGGGG Expression vs Distance")
ggplot(data=gggggg) + geom_point(mapping = aes(x=distToStart, y=FPKM, color=species)) + facet_wrap(~ species) + xlab("Distance to 5' UTR") + ggtitle("GGGGGG Expression vs Distance")
summary(lm(gggggg$FPKM ~ gggggg$distToStart))[[8]]

ggplot(data=aaatcttt) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + xlab("Distance to 5' UTR") + ggtitle("AAATCTTT Expression vs Distance")
ggplot(data=aaatcttt) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + stat_smooth(mapping = aes(x=distToStart, y=FPKM), method = "lm") + xlab("Distance to 5' UTR") + ggtitle("AAATCTTT Expression vs Distance")
ggplot(data=aaatcttt) + geom_point(mapping = aes(x=distToStart, y=FPKM, color=species)) + facet_wrap(~ species) + xlab("Distance to 5' UTR") + ggtitle("AAATCTTT Expression vs Distance")
summary(lm(aaatcttt$FPKM ~ aaatcttt$distToStart))[[8]]

lRSquareds = list()
for(sp in spNames){
  print(sp)
  spAAATCTTT = aaatcttt[which(aaatcttt$species == sp),]
  lRSquareds[[sp]] = summary(lm(spAAATCTTT$FPKM ~ spAAATCTTT$distToStart))[[8]]
}

ggplot(data=tcccgc) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + xlab("Distance to 5' UTR") + ggtitle("TCCGC Expression vs Distance")
ggplot(data=tcccgc) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + stat_smooth(mapping = aes(x=distToStart, y=FPKM), method = "lm") + xlab("Distance to 5' UTR") + ggtitle("TCCGC Expression vs Distance")
ggplot(data=tcccgc) + geom_point(mapping = aes(x=distToStart, y=FPKM, color=species)) + facet_wrap(~ species) + xlab("Distance to 5' UTR") + ggtitle("TCCGC Expression vs Distance")
summary(lm(tcccgc$FPKM ~ tcccgc$distToStart))[[8]]

ggplot(data=tcgtta) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + xlab("Distance to 5' UTR") + ggtitle("TCGTTA Expression vs Distance")
ggplot(data=tcgtta) + geom_point(mapping = aes(x=distToStart, y=FPKM)) + stat_smooth(mapping = aes(x=distToStart, y=FPKM), method = "lm") + xlab("Distance to 5' UTR") + ggtitle("TCGTTA Expression vs Distance")
ggplot(data=tcgtta) + geom_point(mapping = aes(x=distToStart, y=FPKM, color=species)) + facet_wrap(~ species) + xlab("Distance to 5' UTR") + ggtitle("TCGTTA Expression vs Distance")
summary(lm(tcgtta$FPKM ~ tcgtta$distToStart))[[8]] #R^2=0.008

lRSquareds = list()
for(sp in spNames){
  print(sp)
  spTCGTTA = aaatcttt[which(aaatcttt$species == sp),]
  lRSquareds[[sp]] = summary(lm(spTCGTTA$FPKM ~ spTCGTTA$distToStart))[[8]]
}

##### Expression
fpkm_ptet = read.table("/media/tlicknac/Seagate Expansion Drive1/Paramecium_Transcriptomes/ptet-fpkm_prot.tab", header=T)
fpkm_ptet = fpkm_ptet[,c("gene_id", "FPKM")]
fpkm_pbi = read.table("/media/tlicknac/Seagate Expansion Drive1/Paramecium_Transcriptomes/pbi-fpkm_prot.tab", header=T)
fpkm_pbi = fpkm_pbi[,c("gene_id", "FPKM")]
fpkm_psex = read.table("/media/tlicknac/Seagate Expansion Drive1/Paramecium_Transcriptomes/psex-fpkm_prot.tab", header=T)
fpkm_psex = fpkm_psex[,c("gene_id", "FPKM")]

fpkm_ptet_noAAAAAA = fpkm_ptet[- (fpkm_ptet$gene_id %in% aaaaaa$upstreamGene), ]
fpkm_pbi_noAAAAAA =  fpkm_pbi[- (fpkm_pbi$gene_id   %in% aaaaaa$upstreamGene) ,]
fpkm_psex_noAAAAAA = fpkm_psex[- (fpkm_pbi$gene_id  %in% aaaaaa$upstreamGene) ,]
var.test(fpkm_ptet_noAAAAAA$FPKM, aaaaaa[which(aaaaaa$species == "ptet"),"FPKM"])  # p<2.2e-16  --> variances are not equal in ptet
var.test(fpkm_pbi_noAAAAAA$FPKM,  aaaaaa[which(aaaaaa$species == "pbi"), "FPKM"])  # p<2.2e-16  --> variances are not equal in ptet
var.test(fpkm_psex_noAAAAAA$FPKM,  aaaaaa[which(aaaaaa$species == "psex"), "FPKM"])  # 
t.test(aaaaaa_fpkm_ptet$FPKM, fpkm_ptet$FPKM) # p=0.00052  --> means are not equal
t.test(aaaaaa_fpkm_pbi$FPKM, fpkm_pbi$FPKM)  #p=0.06
t.test(aaaaaa_fpkm_psex$FPKM, fpkm_psex$FPKM)  #p=0.47


tttttt_fpkm_ptet = fpkm_ptet[fpkm_ptet$gene_id %in% tttttt$upstreamGene ,]
tttttt_fpkm_pbi =  fpkm_pbi[fpkm_pbi$gene_id   %in% tttttt$upstreamGene ,]
tttttt_fpkm_psex = fpkm_psex[fpkm_pbi$gene_id  %in% tttttt$upstreamGene ,]
t.test(tttttt_fpkm_ptet$FPKM, fpkm_ptet$FPKM, var.equal = F) # p=0.48
t.test(tttttt_fpkm_pbi$FPKM, fpkm_pbi$FPKM, var.equal = F)  #p=0.68
t.test(tttttt_fpkm_psex$FPKM, fpkm_psex$FPKM, var.equal = F)  #p=0.81

cccccc_fpkm_ptet = fpkm_ptet[fpkm_ptet$gene_id %in% cccccc$upstreamGene ,]
cccccc_fpkm_pbi =  fpkm_pbi[fpkm_pbi$gene_id   %in% cccccc$upstreamGene ,]
cccccc_fpkm_psex = fpkm_psex[fpkm_pbi$gene_id %in% cccccc$upstreamGene ,]
t.test(cccccc_fpkm_ptet$FPKM, fpkm_ptet$FPKM)  # p=0.36

gggggg_fpkm_ptet = fpkm_ptet[fpkm_ptet$gene_id %in% gggggg$upstreamGene ,]
gggggg_fpkm_pbi =  fpkm_pbi[fpkm_pbi$gene_id   %in% gggggg$upstreamGene ,]
gggggg_fpkm_psex = fpkm_psex[fpkm_pbi$gene_id  %in% gggggg$upstreamGene ,]
t.test(gggggg_fpkm_ptet$FPKM, fpkm_ptet$FPKM)  # p=0.31








exprTab = sapply(split(aaaaaa$FPKM, aaaaaa$species), summary)
tmpExprDf = as.data.frame(matrix(exprTab), row.names = spNames)
exprTab = str_split_fixed(tmpExprDf$V1, ",", 6)
colnames(exprTab) = colnames(distTab)
exprTab[,1] = 0.01
row.names(exprTab) = spNames
exprTab = as.data.frame(exprTab[-6,])
write.csv(exprTab, file="aaaaaa-expr.csv")
distTab = t(sapply(split(aaaaaa$distToStart, aaaaaa$species), summary))
write.csv(distTab, file="aaaaaa-dist.csv")

