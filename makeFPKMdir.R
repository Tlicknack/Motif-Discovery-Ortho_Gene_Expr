makeFPKMdir = function(fpkmName){
  fpkmTab = read.table(fpkmName, header=T)
  fpkmTab$xp = 2^(fpkmTab$xp)  # jeff stored only copy of file as log2 transformed xpression level
  fpkmOut = as.list(fpkmTab$xp)
  names(fpkmOut) = fpkmTab$gene_id
  
  return(fpkmOut)
}
