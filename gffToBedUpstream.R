gffToBedUpstream = function(speciesNames, speciesGFFs){
  if(speciesNames != speciesGFFs){
    stop("Names and Annotation Files differ in length.")
  }
  message("Beginning bed conversion")
  source("/media/tlicknac/Seagate Expansion Drive1/Scripts/upstreamCoordPlusMinus.R")
  #speciesNames = c("Pdec", "Pnov")
  #speciesGFFs = c("~/Paramecium_GFF/pdec-gene.tab", "~/Paramecium_GFF/pnov-gene.tab")
  
  for(i in 1:length(speciesNames)){
    currentSp = speciesNames[i]
    currentGFF = read.table(speciesGFFs[i], as.is = T)
    allGeneCoords = apply(currentGFF, 1, upstreamCoordPlusMinus)  # This returns a list(scaf, start, stop, geneID, strand)
  }
}