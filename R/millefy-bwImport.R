
makeCoverageMatrixBw <- function(path_bw_files, select, nbin, binsize){
  print(sprintf("Importing BigWig: %s", Sys.time()))
  
  mat <- matrix(0, ncol = nbin, nrow = length(path_bw_files))
  cov <- numeric(width(select))
  
  for(i in seq_along(path_bw_files)){
    path_bw <- path_bw_files[i]
    
    cov <- import(path_bw, selection = select, as = 'NumericList')[[1]]
    mat[i,] <- coverageIntoFixedBins(cov, binsize)
  }
  
  mat
}
