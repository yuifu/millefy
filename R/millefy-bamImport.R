#' Calculate normalization factors for visualization of BAM files
#' 
#' @param bam_files A character vector containing paths to BAM files
#' @param scaleFactor A numeric vector of normalization factors (1/(total mapped reads)*scaleFactor)
#' @return bamNormFactors numeric vector
#' @export
#'
#' @examples
#' # Gene annotation track (For faster performance, try to use \code{dt_gtf} paramter)
#' bam_files = c("example1.bam", "example2.bam")
#' normFactors = calcBamNormFactors(bam_files)
calcBamNormFactors <- function(bam_files, scaleFactor = 10^6){
  counts <- getBamCounts(bam_files)
  1/counts*scaleFactor
}

getBamCounts <- function(bam_files){
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))
  sapply(bam_files, function(x){
    countBam(x, param = param)$records
  })
}



#######################################################################################################

makeCoverageMatrix <- function(path_bam_files, select, nbin, binsize, normFactor){
  print(sprintf("Start importing BAM: %s", Sys.time()))
  
  mat <- matrix(0, ncol = nbin, nrow = length(path_bam_files))
  
  for(i in seq_along(path_bam_files)){
    path_bam <- path_bam_files[i]
    
    cov <- getCoverageFromSplicedBamImport(path_bam, select)
    cov <- coverageIntoFixedBins(cov, binsize)
    # cov[cov == 0] <- NA_integer_
    # print(dim(mat))
    # print(length(cov))
    mat[i,] <- cov * normFactor[i]
  }
  
  print(sprintf("End importing BAM: %s", Sys.time()))
  mat
}



getCoverageFromSplicedBamImport <- function(file, selection){
 gr <- splicedBamImport2(file, selection) 
 # gr <- splicedBamImport(file, selection) 
 if(length(gr) == 0){ 
    return(rep(0, end(selection) - start(selection) + 1)) 
 } 
 st <- min(start(gr)) 
 cov <- coverage(gr)[[1]] 
 rep(cov@values, cov@lengths)[start(selection):end(selection)]  # time consuming?
}

splicedBamImport2 <- function (file, selection) {
  if (!file.exists(paste(file, "bai", sep = ".")))
    stop("Unable to find index for BAM file '", file, "'. You can
         build an index using the following command:\n\t",
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  sinfo <- scanBamHeader(file)[[1]]
  res <- if (!as.character(seqnames(selection)[1]) %in%
             names(sinfo$targets)) {
    mcols(selection) <- DataFrame(score = 0)
    selection
  }else {
    param <- ScanBamParam(what = c("pos", "strand", "cigar"),
                          which = selection, flag =
                            scanBamFlag(isUnmappedQuery = FALSE))
    x <- scanBam(file, param = param)[[1]]
    # print(str(x))
    if(length(x[["pos"]]) == 0){
      return(numeric(0))
    }
    #     return(x[["cigar"]])
    
    
    lpwl <- lapply(seq_along(x[["pos"]]), function(y){returnPosAndWidthForCigar2(x[["pos"]][y], x[["cigar"]][y])})
    # lpl <- lapply(seq_along(x[["pos"]]), function(y){returnPosForCigar(x[["pos"]][y], x[["cigar"]][y])})
    # lwl <- lapply(seq_along(x[["pos"]]), function(y){returnWidthForCigar(x[["pos"]][y], x[["cigar"]][y])})
    lpl <- lapply(lpwl, function(y){y$pos})
    lwl <- lapply(lpwl, function(y){y$width})
    
    gr <- GRanges(strand=rep(x[["strand"]], sapply(lpl, length)),  # Rle object, character vector, or factor containing the strand information.
                  ranges=IRanges(start = unlist(lpl), width = unlist(lwl)),  #An IRanges object containing the ranges.
                  seqnames=seqnames(selection)[1])
    return(gr)

  }
}

splicedBamImport <- function (file, selection) {
  if (!file.exists(paste(file, "bai", sep = ".")))
    stop("Unable to find index for BAM file '", file, "'. You can
         build an index using the following command:\n\t",
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  sinfo <- scanBamHeader(file)[[1]]
  res <- if (!as.character(seqnames(selection)[1]) %in%
             names(sinfo$targets)) {
    mcols(selection) <- DataFrame(score = 0)
    selection
  }else {
    param <- ScanBamParam(what = c("pos", "strand", "cigar"),
                          which = selection, flag =
                            scanBamFlag(isUnmappedQuery = FALSE))
    x <- scanBam(file, param = param)[[1]]
    # print(str(x))
    if(length(x[["pos"]]) == 0){
      return(numeric(0))
    }
    #     return(x[["cigar"]])
    
    
    gr <- GRanges(strand=rep(x[["strand"]], returnRepForCigars(x[["pos"]], x[["cigar"]])),  # Rle object, character vector, or factor containing the strand information.
                  ranges=returnIRangesForCigars(x[["pos"]], x[["cigar"]]),  #An IRanges object containing the ranges.
                  seqnames=seqnames(selection)[1])
    return(gr)
    
    
    
    #    grs <- split(gr, strand(gr))
    #    cov <- lapply(grs[c("+", "-")], function(y)
    #      coverage(ranges(y),
    #               width=end(selection)))
    #    pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y),
    #                                                        end(y))))))
    #    if(length(pos)==0){
    #      mcols(selection) <- DataFrame(plus=0, minus=0)
    #      selection
    #    }else{
    #      GRanges(seqnames = seqnames(selection)[1],
    #              ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)),
    #              plus=as.numeric(cov[["+"]][head(pos, -1)]),
    #              minus=-as.numeric(cov[["-"]][head(pos, -1)]))
    #    }
    #  }
    #  return(res)
  }
}

returnRepForCigars <- function(poss, cigs){
  lapply(seq_along(poss), function(x){length(returnPosForCigar(poss[x], cigs[x]))}) %>% unlist
}

returnIRangesForCigars <- function(poss, cigs){
  lp <- lapply(seq_along(poss), function(x){returnPosForCigar(poss[x], cigs[x])}) %>% unlist 
  lw <- lapply(seq_along(poss), function(x){returnWidthForCigar(poss[x], cigs[x])}) %>% unlist
  IRanges(lp, width=lw)
}
# returnIRangesForCigars(poss,cigs)

returnPosAndWidthForCigar <- function(pos, cig){
  va <- strsplit(cig, "[A-Z]") %>% unlist %>% as.integer()
  al <- strsplit(cig, "[0-9]+") %>% unlist
  al <- al[al != ""]

  newpos <- pos + returnOffset(va, al)  
  width <- va
  OK <- sapply(al, returnOK)
  list(pos = newpos[OK], width = width[OK])
}

returnPosAndWidthForCigar2 <- function(pos, cig){
  va <- as.integer(unlist(strsplit(cig, "[A-Z]")))
  al <- unlist(strsplit(cig, "[0-9]+"))
  al <- al[al != ""]

  newpos <- pos + returnOffset(va, al)  
  width <- va
  OK <- sapply(al, returnOK)
  list(pos = newpos[OK], width = width[OK])
}

returnPosAndWidthForCigar3 <- function(pos, cig){
  va <- as.integer(unlist(strsplit(cig, "[A-Z]")))
  al <- unlist(strsplit(cig, "[0-9]+"))
  al <- al[al != ""]
  
  if(length(al) == 1 && al == "M"){
    list(pos = pos, witdth = va)
  }
  
  newpos <- pos + returnOffset(va, al)  
  width <- va
  OK <- sapply(al, returnOK)
  list(pos = newpos[OK], width = width[OK])
}

returnWidthForCigar <- function(pos, cig){
  
  va <- strsplit(cig, "[A-Z]") %>% unlist %>% as.integer()
  al <- strsplit(cig, "[0-9]+") %>% unlist
  al <- al[al != ""]
  
  width <- va
  OK <- sapply(al, returnOK)
  width[OK]
}

returnPosForCigar <- function(pos, cig){
  
  va <- strsplit(cig, "[A-Z]") %>% unlist %>% as.integer()
  al <- strsplit(cig, "[0-9]+") %>% unlist
  al <- al[al != ""]
  
  newpos <- pos + returnOffset(va, al)
  OK <- sapply(al, returnOK)
  newpos[OK]
}


# returnRangesForCigar <- function(pos, cig){
#   
#   va <- strsplit(cig, "[A-Z]") %>% unlist %>% as.integer()
#   al <- strsplit(cig, "[0-9]+") %>% unlist
#   al <- al[al != ""]
#   
#   newpos <- pos + returnOffset(va, al)
#   width <- va
#   OK <- sapply(al, returnOK)
#   list(newpos[OK], width[OK])
# }

returnOffset <- function(va, al){
  if(length(va) == 1){
    return(0)
  }else{
    if(al[1]=="S" || al[1] == "H"){
      va[1] = 0
    }
    c(0, sapply(1:(length(va)-1), function(x){sum(va[1:x])}))
  }
}

returnOK <- function(a){
  switch (a,
    "S" = FALSE,
    "M" = TRUE,
    "N" = FALSE,
    "I" = FALSE,
    "D" = FALSE
  )
}


# https://support.bioconductor.org/p/56047/
strandedBamImport <- function (file, selection) {
  if (!file.exists(paste(file, "bai", sep = ".")))
    stop("Unable to find index for BAM file '", file, "'. You can
         build an index using the following command:\n\t",
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  sinfo <- scanBamHeader(file)[[1]]
  res <- if (!as.character(seqnames(selection)[1]) %in%
             names(sinfo$targets)) {
    mcols(selection) <- DataFrame(score = 0)
    selection
  }else {
    param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                          which = selection, flag =
                            scanBamFlag(isUnmappedQuery = FALSE))
    x <- scanBam(file, param = param)[[1]]
    # print(str(x))
    if(length(x[["pos"]]) == 0){
      return(numeric(0))
    }
    gr <- GRanges(strand=x[["strand"]], ranges=IRanges(x[["pos"]],
                                                       width = x[["qwidth"]]), seqnames=seqnames(selection)[1])
    grs <- split(gr, strand(gr))
    cov <- lapply(grs[c("+", "-")], function(y)
      coverage(ranges(y),
               width=end(selection)))
    pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y),
                                                        end(y))))))
    if(length(pos)==0){
      mcols(selection) <- DataFrame(plus=0, minus=0)
      selection
    }else{
      GRanges(seqnames = seqnames(selection)[1],
              ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)),
              plus=as.numeric(cov[["+"]][head(pos, -1)]),
              minus=-as.numeric(cov[["-"]][head(pos, -1)]))
    }
  }
  return(res)
}
