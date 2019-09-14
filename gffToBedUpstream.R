gffToBedUpstream = function(speciesNames, speciesGFFs){
  if(speciesNames != speciesGFFs){
    stop("Names and Annotation Files differ in length.")
  }
  message("Beginning bed conversion")
  source("/media/tlicknac/Seagate Expansion Drive1/Scripts/upstreamCoordPlusMinus.R")
  #speciesNames = c("Pdec", "Pnov")
  #speciesGFFs = c("~/Paramecium_GFF/pdec-gene.tab", "~/Paramecium_GFF/pnov-gene.tab")
  
  lBEDs = list()
  
  for(i in 1:length(speciesNames)){
    currentSp = speciesNames[i]
    currentGFF = read.table(speciesGFFs[i], as.is = T)
    allGeneCoords = apply(currentGFF, 1, upstreamCoordPlusMinus)  # This returns a list(scaf, start, stop, geneID, strand)
    outName = paste(currentSp, "_Upstream.bed", sep="")
    outBED = data.frame(matrix(unlist(allGeneCoords), nrow=nrow(currentGFF), byrow=T))  # This converts the list into a BED-like Data Frame
    lBEDs[[currentSp]] = outBED
  }
  return(lBEDs)
}
