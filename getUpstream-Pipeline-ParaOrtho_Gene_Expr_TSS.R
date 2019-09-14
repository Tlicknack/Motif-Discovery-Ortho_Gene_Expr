# This is the master R script that will generate a FASTA of upstream sequences for both ParaOrthologs, AllGenes, and ExpressionLevels

# This will use many Rscripts that will be sourced when appropriate:
  # upstreamCoordPlusMinus.R
  # gffToBedUpstream.R

# This will require the following files: (No need to read in each individually... functions should do that one at a time when needed... except POFF)
  # Genome Assembly: FASTA
  # Genome Annotation: GFF
  # ParaOrtholog Table: POFF
  # RNAseq Table: FPKM

#Get List of species names to create BED files with upstream (-200,0) coordinates for each gene
source("/N/dc2/scratch/tlicknac/MotifDiscovery/gffToBedUpstream.R")
poffTable = read.table("~/Paramecium_POFF/all_aurelias-cut-mod.poff", header=T)
poffTable = poffTable[,-2]
spNames = sort(names(poffTable[-1]))
gffDir = "~/Paramecium_GFF"
gffFiles = list.files(gffDir, full.names = T)
gffFiles = sort(gffFiles[grep("*-gene.tab", gffFiles)])
message("Make sure that the files are in the proper order!! Sort() should have done it, but double check before proceeding")
lBEDs = gffToBedUpstream(speciesNames=spNames, gffFiles)
