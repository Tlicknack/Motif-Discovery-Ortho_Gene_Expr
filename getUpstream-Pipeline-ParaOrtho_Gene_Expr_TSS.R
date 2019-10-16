# This is the master R script that will generate a FASTA of upstream sequences for both ParaOrthologs, AllGenes, and ExpressionLevels

# This will use many Rscripts that will be sourced when appropriate:
#Part1:
# makeFPKMdir.R
# upstreamCoordPlusMinus.R
# gffToBedUpstream.R
#Part2:
# getUpstreamFromBED.R
#Part3:
#Part4:
# fixScafNamesGRanges.R
# BEDtoGRanges.R
# matchMotifFASTA.R

# This will require the following files: (No need to read in each individually... functions should do that one at a time when needed... except POFF)
# Genome Assembly: FASTA
# Genome Annotation: GFF
# ParaOrtholog Table: POFF
# RNAseq Table: FPKM

#Get List of species names to create BED files with upstream (-200,0) coordinates for each gene
source("/N/dc2/scratch/tlicknac/MotifDiscovery/Rscripts/gffToBedUpstream.R")
poffTable = read.table("~/Paramecium_POFF/all_aurelias-cut-mod.poff", header=T)
poffTable = poffTable[,-2]
speciesNames = sort(names(poffTable[-1]))
gffDir = "~/Paramecium_GFF"
gffFiles = list.files(gffDir, full.names = T)
gffFiles = sort(gffFiles[grep("*-gene.tab", gffFiles)])
gffFiles[2] = "/N/u/tlicknac/Carbonate/Paramecium_GFF/pcaud-oldgene.tab" #### HARD CODED... pcaud gene names are confusing. 
message("Make sure that the files are in the proper order!! Sort() should have done it, but double check before proceeding")
lBEDs = gffToBedUpstream(speciesNames=speciesNames, gffFiles) #lBEDs is a list with a layer of species names, and reach species name has a data frame of upstream coords for every gene
rm(list = ls()[which(ls() != "lBEDs")])
save.image(file="/N/dc2/scratch/tlicknac/MotifDiscovery/RData/upstreamBEDcoords.Rdata")
#####

#Make database for FPKM Vegetative growth expression levels
fpkmDir = "~/Paramecium_FPKM"
fpkmFiles = list.files(fpkmDir, pattern="*-xp.tab", full.names = T)
fpkmSPnames = spNames[-6]
fpkmDb = vector(mode = "list", length = length(fpkmSPnames))
names(fpkmDb) = fpkmSPnames
fpkmDb$pmultimic = NULL
for(e in 1:length(fpkmFiles)){
  spName = fpkmSPnames[e]
  if(spName != "pmultimic"){
    fpkmDb[[spName]] = makeFPKMdir(fpkmFiles[e])     
  }
}
setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/")
rm(list = ls()[which(ls() != "fpkmDb")])
save.image("fpkmDir.RData")

#####

#Make database of gene orientation and intergenic size


#####

#Assemble upstream fasta's for every gene in every genome
#load("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/upstreamBEDcoords.Rdata")
library(seqinr)
source("/N/dc2/scratch/tlicknac/MotifDiscovery/Rscripts/getUpstreamFromBED.R")
speciesNames = sort(names(poffTable[-1]))

faDir = "~/Paramecium_FASTA/StandardNames"  # make sure dir only have the 16 aurelias of this will get messed up
faFiles = sort(list.files(faDir, full.names = T))

lUpstreams = vector(mode="list", length=length(speciesNames))
names(lUpstreams) = speciesNames

setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/")
for(i in 1:length(speciesNames)){  # This takes about 9 min per species of ~40k genes !!!!!!!!!!!
  spFASTA = read.fasta(faFiles[[i]], forceDNAtolower = T, as.string = T)
  lUpstreams[[speciesNames[i]]] = getUpstreamFromBED(spBED=lBEDs[[i]], spFASTA=spFASTA) #lUpstreams is a list with a layer of species names, and each species name has every gene and its upstream region
  imageName = paste("lUpstreams_", i, "th-Species", ".RData", sep="")
  save.image(imageName)
}

for(z in 1:length(lUpstreams)){  #Fix this to remove the weird NULL
  spUpstream = lUpstreams[[z]]
  spUpstream =  spUpstream[-which(as.character(spUpstream) == "NULL")]
  lUpstreams[[z]] = spUpstream
}
rm(list = ls()[which(ls() != "lUpstreams")])
save.image("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/lUpstreams_no_null.RData")
#####


# Assemble fasta files for upstream of each species
#load("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/lUpstreams_16th-Species.RData")
setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/upstreamFasta/WholeGenome")
for(p in 1:length(lUpstreams)){  # this takes less than 30 seconds
  spUpstreams = lUpstreams[[p]]
  seqs = as.character(spUpstreams)
  seqNames = names(spUpstreams)
  fileName = paste(speciesNames[p], "_allPromoters.fasta", sep="")
  write.fasta(as.list(seqs), seqNames, file.out = fileName)
}
# Fix this: bash
# sed -e '/pattern/,+5d' file.txt
# Background File: in bash
# cd /N/dc2/scratch/tlicknac/MotifDiscovery/upstreamFasta/backgroundfile  
# cat ../WholeGenome/* > ./all_promoters_aurelias_caud_multimic.fasta
#####


#Time to assemble upstream fasta files for each WGD gene family
#load("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/lUpstreams_16th-Species.RData")

setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/upstreamFasta/ParaOrtho")
poffTable = read.table("~/Paramecium_POFF/all_aurelias-cut-mod.poff", header=T)
poffTable = poffTable[,-2]
speciesNames = sort(names(poffTable[-1]))

for(q in 1:nrow(poffTable)){
  poffRow = poffTable[q,]
  currentWGD = as.character(poffRow$WGD)
  poffRow$WGD = NULL
  vIDs = c()
  vSeqs = c()
  
  for(r in 1:length(speciesNames)){
    currentSp = speciesNames[r]
    poffRowCurrentSp = strsplit(as.character(unlist(poffRow[currentSp])), ",")
    
    for(s in 1:length(poffRowCurrentSp[[1]])){
      geneS = poffRowCurrentSp[[1]][s]
      
      if(geneS != "." & is.na(geneS) == F & geneS != "NA"){
        geneID = names(lUpstreams[[currentSp]] [which(names(lUpstreams[[currentSp]]) == geneS)])
        geneSeq = as.character( lUpstreams[[currentSp]] [which(names(lUpstreams[[currentSp]]) == geneS)] )
        
        if(identical(geneID, character(0)) == F){
          vIDs = append(vIDs, geneID)
          vSeqs = append(vSeqs, geneSeq)
        }
      }
    }
  }
  outFile = paste(currentWGD, "_upstream.fasta", sep = "")
  write.fasta(as.list(vSeqs), vIDs, file.out = outFile)
  outMessage = paste("WGD upstream fasta number", q, "has been written in ParaOrtho dir.")
  message(outMessage)
}
# Background File: in bash
# cd /N/dc2/scratch/tlicknac/MotifDiscovery/upstreamFasta/backgroundfile  
# cat ../ParaOrtho/* > ./allPOFFupstreams.fasta
#####


# Time to assemble files based on their expression levels.
library(seqinr)
setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/upstreamFasta/ExprLevel/")
#load("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/lUpstreams_16th-Species.RData")

fpkmDir = "~/Paramecium_FPKM"
fpkmFiles = list.files(fpkmDir, pattern="*-fpkm_prot.tab", full.names = T)

#for(sp in speciesNames){
sp = "ptet"
if(sp == "ptet"){
  fpkmTet = read.table(fpkmFiles[9], header=T)
  fpkmTet = as.data.frame(matrix(c(as.character(fpkmTet$gene_id), as.numeric(fpkmTet$FPKM)), ncol=2, nrow=nrow(fpkmTet)))
  colnames(fpkmTet) = c("geneID", "fpkm")
  fpkmTet$geneID = as.character(fpkmTet$geneID)
  fpkmTet$fpkm = as.numeric(as.character(fpkmTet$fpkm))
  sumExpr = summary(fpkmTet$fpkm)
  lowExpr = fpkmTet[which(fpkmTet$fpkm < sumExpr[2] & fpkmTet$fpkm > 0.05),]              #lowExpr are genes with fpkm between 0.05 and the 1st quartile
  midExpr = fpkmTet[which(fpkmTet$fpkm < sumExpr[3] & fpkmTet$fpkm > sumExpr[2]),]        #midExpr are gene with fpkm between 1st and median
  highExpr = fpkmTet[which(fpkmTet$fpkm < sumExpr[5] & fpkmTet$fpkm > sumExpr[3]),]       #midExpr are gene with fpkm between median and 3rd quartile
  veryHighExpr = fpkmTet[which(fpkmTet$fpkm < sumExpr[6] & fpkmTet$fpkm > sumExpr[5]),]   #midExpr are gene with fpkm between 3rd quartile and max
  # each category has 10,382 rows, excluding lowExpr which has 4642 due to the exlcusion of genes below fpkm > 0.05
  
  upstreamsTet = lUpstreams[["ptet"]]
  upstreamsTet = upstreamsTet[-which(as.character(upstreamsTet) == "NULL")]
  
  lowExprUpstreams = upstreamsTet[lowExpr$geneID]
  lowExprUpstreams = lowExprUpstreams[-which( as.character(lowExprUpstreams) == "NULL")]
  write.fasta(as.list(as.character(lowExprUpstreams)), names(lowExprUpstreams), file.out = "ptet_lowExpression_upstreams.fasta")
  midExprUpstreams = upstreamsTet[midExpr$geneID]
  midExprUpstreams = midExprUpstreams[-which(as.character(midExprUpstreams) == "NULL")]
  write.fasta(as.list(as.character(midExprUpstreams)), names(midExprUpstreams), file.out = "ptet_midExpression_upstreams.fasta")
  highExprUpstreams = upstreamsTet[highExpr$geneID]
  highExprUpstreams = highExprUpstreams[-which(as.character(highExprUpstreams) == "NULL")]
  write.fasta(as.list(as.character(highExprUpstreams)), names(highExprUpstreams), file.out = "ptet_highExpression_upstreams.fasta")
  veryHighExprUpstreams = upstreamsTet[veryHighExpr$geneID]
  veryHighExprUpstreams = veryHighExprUpstreams[-which(as.character(veryHighExprUpstreams) == "NULL")]
  write.fasta(as.list(as.character(veryHighExprUpstreams)), names(veryHighExprUpstreams), file.out = "ptet_veryHighExpression_upstreams.fasta")
}
#}
# Background File: already made (each spp. whole-genome file)
# Neg Control File: in bash
# cd /N/dc2/scratch/tlicknac/MotifDiscovery/ExprLevel/ 
# cat ptet_highExpression_upstreams.fasta ptet_veryHighExpression_upstreams.fasta > ptet_high_veryHigh-Expr_upstream.fasta
# cat ptet_lowExpression_upstreams.fasta ptet_midExpression_upstreams.fasta > ptet_low_mid-Expr_upstream.fasta
#####


# Now to assembly upstream fastas on the basis of TSS shape
# Taylor and I have been working on this quite a bit. So lets load some objects from that effort
setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/upstreamFasta/TssShape")
library(seqinr)
#load("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/lUpstreams_16th-Species.RData")
source("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/scripts/buildTSSDb.R")
source("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/scripts/tssToList.R")
sampleDb <- buildTSSDb(speciesNames=c("pdec","poct","pnov", "pjen"), "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Pdec.txt", "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Poct.txt", "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Pnov.txt", "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Pjen.txt")
# We can only do this in 4 spp at the moment

for(y in 1:length(sampleDb)){
  lSharp = list()
  lBroad = list()
  
  spDb = sampleDb[[y]]
  currentSp = names(sampleDb[y])
  currentUpstreams = lUpstreams[[currentSp]]
  # I have no idea how to work with this list >:( ... will just iterate through it...
  
  for(z in 1:length(spDb)){
    currentTss = spDb[z]
    
    if(is.null(currentTss[[1]]) == F){  #some TSSs are NULL
      currentID = gsub("1.G", "1.P",  names(currentTss))
      currentMSI = currentTss[[1]]$tsrMSI
      currentnTAGs = currentTss[[1]]$nTAGs
      
      outSeq = as.character(currentUpstreams[which(names(currentUpstreams) == currentID)])
      
      if(currentMSI > 0.8){
        lSharp[[currentID]] = outSeq
        outMessage = paste("Wrote to lSharp due to tsrMSI = ", currentMSI, sep="")
      }
      if(currentMSI < 0.8){
        lBroad[[currentID]] = outSeq
        outMessage = paste("Wrote to lSharp due to tsrMSI = ", currentMSI, sep="")
      }
      message(outMessage)
    } else{
      nullMessage = paste("NULL TSS: ", names(currentTss), sep="")
      message(nullMessage)
    }
  }
  outSharp = paste(currentSp, "_Sharp-promoter_upstream.fasta", sep="")
  outBroad = paste(currentSp, "_Broad-promoter_upstream.fasta", sep="")
  write.fasta(as.list(as.character(lSharp)), names(lSharp), file.out = outSharp)
  write.fasta(as.list(as.character(lBroad)), names(lBroad), file.out = outBroad)
  message("Done with species.")
}

############################################################################################################################################

# Run MEME for all types of Data

# This is all in bash: Basic Outline

# 

############################################################################################################################################

# Create an object of all possible 6mers and their relative distributions 

library(gtools)
library(seqinr)
sixmers = gtools::permutations(n=4, r=6, v=c("A", "T", "C", "G"), repeats.allowed=T)
sixMers = apply(sixmers, 1, c2s)
motifs = vector(mode="list", length=4096)
names(motifs) = sixMers


####################################################################################
# Map motifs back to the genome:

#load("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/upstreamBEDcoords.Rdata")
library(Biostrings)
library(rtracklayer)
source("/N/dc2/scratch/tlicknac/MotifDiscovery/Rscripts/fixScafNamesGRanges.R")
source("/N/dc2/scratch/tlicknac/MotifDiscovery/Rscripts/BEDtoGRanges.R")
source("/N/dc2/scratch/tlicknac/MotifDiscovery/Rscripts/matchMotifFASTA.R")
source("/N/dc2/scratch/tlicknac/MotifDiscovery/Rscripts/checkStrandGRangesReturnDistance.R")

setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/RData")

faFileDir = "/N/u/tlicknac/Carbonate/Paramecium_FASTA/StandardNames"
faFiles = list.files(faFileDir, full.names = T)

motifDir = list()

motif = "TCAGTT"

motifDir[[motif]] = vector(mode="list", length=length(spNames))  #append motifDir
names(motifDir[[motif]]) = spNames

for(m in 1:length(faFiles)){                          # THIS TAKES ABOUT 6 MIN PER 8mer!!!!
  spName = spNames[m]
  #get motif matches
  stringFASTA = readDNAStringSet(faFiles[m])
  
  if(length( strsplit( names( stringFASTA)[1], " ")[[1]]) > 1){ #fix these gosh darn scaffold names..
    names(stringFASTA) = fixScafNamesGRanges(stringFASTA)  #FUNCTION
  }    
  dfMatches = matchMotifFASTA(motif, stringFASTA)  #FUNCTION
  
  #prepare upstream BEDs to search
  spBED = lBEDs[[m]]
  spGFF = BEDtoGRanges(spBED) #FUNCTION
  
  #find motif instances in BED coordinates
  qHits = dfMatches[queryHits(suppressWarnings(findOverlaps(dfMatches, spGFF)))]
  sHits = spGFF[subjectHits(suppressWarnings(findOverlaps(dfMatches, spGFF)))]
  qHits$upstreamGene = sHits$upstreamGene
  
  #get distance to start codon
  vDistHits = c()
  for(h in 1:length(qHits)){  #arghh I have to loop GRanges objects... no lapply() functionality
    distHits = checkStrandGRangesReturnDistance(qHits[h], sHits[h])
    vDistHits = append(vDistHits, distHits)
  }
  qHits$distToStart = vDistHits
  
  #Convert to df to append to R Object
  tmpDf = data.frame(seqname=seqnames(qHits), start=start(qHits), end=end(qHits), strand=strand(qHits), upstreamGene= qHits$upstreamGene, distToStart = qHits$distToStart)
  tmpDf$species = spName
  
  motifDir[[motif]][[spName]] = tmpDf
}
# Export R data
rm(list = ls()[which(ls() != "motifDir" & ls() != "motif")])
outImage = paste("motif_", motif, ".RData", sep="")
save.image(file=outImage)

#####
# Add Expression data
library(dplyr)
motif = "TCCCGC"

load("/N/dc2/scratch/tlicknac/MotifDiscovery/RData/fpkmDb.RData")
load("../RData/motif_TCAGTT.RData")
outTable = bind_rows(motifDir[[motif]])
outTable$FPKM = NA
outTable$motif = motif

for(zzz in 1:nrow(outTable)){  # THIS TAKES 0-40MIN depending on number of motif hits.... (CCCCCC was instant, AAAAAA took 40min)
  outRow = outTable[zzz,]
  newVal = as.numeric(fpkmDb[[outRow$species]] [which( names(fpkmDb[[outRow$species]]) == outRow$upstreamGene)])
  
  if(identical(newVal, numeric(0)) == F){
    outTable[zzz,"FPKM"] = newVal 
  }
}
setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/RData")
rm(list= ls()[which(ls() != "outTable" & ls() != "motif")])
outImage2 = paste("motif_", motif, "_withExpression.RData", sep="")
save.image(outImage2)
setwd("/N/dc2/scratch/tlicknac/MotifDiscovery/motifDistributions")
outTab = paste("motif_", motif, "_withExpression.csv", sep="")
write.csv(outTable, file=outTab)
#for(spName in spNames){
#  lFpkm = fpkmDb[[spName]]
#fpkmDb[[outTable$species[1]]] [which(names(fpkmDb[[outTable$species[1] ]]) %in% outTable$upstreamGene)]
#}








