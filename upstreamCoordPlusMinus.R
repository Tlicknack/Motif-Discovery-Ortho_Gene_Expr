upstreamCoordPlusMinus = function(rowVector){
  #DONE
  scaf = as.character(unlist(rowVector[1]))
  firstCoord = as.numeric(unlist(rowVector[2]))
  secondCoord = as.numeric(unlist(rowVector[3]))
  strand = as.character(unlist(rowVector[4]))
  geneID = as.character(unlist(rowVector[5]))
  
  if(strand == "+"){
    start = as.numeric(firstCoord)-200
    stop = as.numeric(firstCoord)
  }
  if(strand == "-"){
    start = as.numeric(secondCoord)
    stop = as.numeric(secondCoord)+200
  }
  return(list(scaf, start, stop, geneID, strand))
}
