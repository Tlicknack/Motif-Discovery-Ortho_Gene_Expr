getUpstreamFromBED = function(spBED, spFASTA){
#DONE
  message("Beginning to assemble upstreams.")
  
  nRowBED = nrow(spBED)
  print(nRowBED)
  tmpListUpstreams = vector(mode="list", length=nRowBED)
  names(tmpListUpstreams) = c(as.character(unlist(spBED[4])))
  
  for(i in 1:nRowBED){
    spBEDrow = spBED[i,]
    print(spBEDrow)
    spFASTAscaf = spFASTA[which(getName(spFASTA) == as.character(unlist(spBEDrow[1])))]
    
    if(as.numeric(as.character(unlist(spBEDrow[3]))) < nchar(as.character(spFASTAscaf))){  # Only get an upstream region if the 2nd position beyond the scaffold boundaries
      upstream = substr(as.character(spFASTAscaf), start = as.numeric(as.character(unlist(spBEDrow[2]))), stop = as.numeric(as.character(unlist(spBEDrow[3]))))
      
      if(as.character(unlist(spBEDrow[5])) == "-"){
        upstream = c2s(rev(comp(s2c(upstream))))
      }
    }
    tmpListUpstreams[[as.character(unlist(spBEDrow[4]))]] = upstream
    message("tmpListUpstreams has been appended.")
  }
  return(tmpListUpstreams)
}
