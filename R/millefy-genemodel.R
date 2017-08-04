#' Create a data.table object that contains informaiton on exonic regions
#' 
#' @param path_gtf Path to a GTF file
#' @return dtExon data.table object
#' @export
#'
#' @examples
#' # Gene annotation track (For faster performance, try to use \code{dt_gtf} paramter)
#' path_gtf = "annotation/gencode.gtf"
#' dt_gtf_exon <- gtfToDtExon(path_gtf) 
#' geneTrack1 <- list(path_gtf = path_gtf, dt_gtf = dt_gtf_exon, label = "GENCODE")
gtfToDtExon <- function(path_gtf){
 
  dt_gtf <- fread(path_gtf, header = F) %>% subset(V3 == "exon")
  setnames(dt_gtf, c("V1", "V4", "V5", "V7"), c("chr", "start", "end", "strand"))
  # separate(dt_gtf, col = "V9", into = c("VV1", "VV2"), sep = "transcript_id \"") %>% 
  #   separate(col = "VV2", into = c("transcript_id", "VV4"), sep = "\"; ")
  dt_gtf[, transcript_id := gsub(".+transcript_id \"", "", V9)]
  dt_gtf[, transcript_id := gsub("\";.+", "", transcript_id)]
  dt_gtf[, transcript_name := gsub(".+transcript_name \"", "", V9)]
  dt_gtf[, transcript_name := gsub("\";.+", "", transcript_name)]

  dt_gtf %>% select(chr, start, end, strand, transcript_id, transcript_name)
}

plotGeneModels4 <- function(path_gtf, x_chr, x_start, x_end, dt_gtf_exon, show_transcript_id = TRUE){
  
  if(missing(dt_gtf_exon)){
    dt_gtf_exon <- gtfToDtExon(path_gtf) 
  }
  
  # transcript_ids <- dt_gtf_exon %>% subset(chr == x_chr & start <= x_end & end >= x_start) %>% {.[,unique(transcript_id)]}
  transcript_ids <- dt_gtf_exon[chr == x_chr & start <= x_end & end >= x_start, unique(transcript_id)]
  dt_gtf_exon %>% setkey(transcript_id)
  dt_gtf_exon_local <- dt_gtf_exon[transcript_ids]
  transcript_names <- dt_gtf_exon_local[, unique(transcript_name)]
  
  if(length(transcript_ids)==0){
    cat("No gene models in this region\n")
  }else{
    
    pushViewport(viewport(xscale = c(x_start, x_end), yscale = c(0,length(transcript_ids)+1), clip = "off"))
    
    
    bottom <- 1
    h <- 1
    hh <- 0.3
    len_arrow = unit(hh/2, "lines")
    for(i in 1:length(transcript_ids)){
      starts <- dt_gtf_exon_local %>% subset(transcript_id == transcript_ids[i]) %>% select(start) %>% unlist
      ends <- dt_gtf_exon_local %>% subset(transcript_id == transcript_ids[i]) %>% select(end) %>% unlist
      
      grid.lines(x=c(starts, ends), y=rep(bottom+h*(i-1),2),
                 default.units = "native", gp = gpar(col = "#554C4D"))
      
      for(j in 1:length(starts)){
        # print(starts[j]) #### debug
        grid.rect(
          x = unit((starts[j]-x_start)/(x_end - x_start + 1), "npc"),
          y = bottom+h*(i-1),
          width = unit((ends[j]-starts[j]+1)/(x_end - x_start + 1), "npc"),
          height = hh,
          gp = gpar(fill = "#ff9900"),
          default.units = "native",
          hjust = 0
        )      
      }
      
      if("strand" %in% colnames(dt_gtf_exon_local)){
        strand <- dt_gtf_exon_local %>% subset(transcript_id == transcript_ids[i]) %>% select(strand) %>% unlist %>% unique
        if(FALSE){
          angle <- ifelse(strand == "+", 140, 40)
          pos_arrow <- seq(min(starts), max(ends), ceiling((x_end-x_start+1)/100))
          # print(strand)
          grid.polyline(x=rep(pos_arrow,each=2), y=rep(bottom+h*(i-1),2*length(pos_arrow)),
                        id = rep(seq_along(pos_arrow), each = 2),
                        default.units = "native", 
                        gp = gpar(lwd = 0.5, col = "black"),
                        arrow = arrow(angle = angle, length=len_arrow))
        }
        
        if(strand == "+"){
          gene_text <- sprintf("> %s (%s)", transcript_ids[i], transcript_names[i])
        }else if (strand == "-"){
          gene_text <- sprintf("< %s (%s)", transcript_ids[i], transcript_names[i])
        }else{
          gene_text <- sprintf("%s (%s)", transcript_ids[i], transcript_names[i])
        }
      }else{
        gene_text <- sprintf("%s (%s)", transcript_ids[i], transcript_names[i])
      }
      
      xpos <- (min(x_end, max(ends)) + max(x_start, min(starts)) - 2* x_start)/2/(x_end - x_start + 1)
      
      if(show_transcript_id){
        grid.text(gene_text,
                  x = unit(xpos, "npc"), y=unit(rep(bottom+h*(i-1)+hh,2), "native"),
                  hjust = 0.5, vjust = 0,
                  gp = gpar(cex = 0.4)
        )
      }      
      
    }
    popViewport()
  }  
  
}
