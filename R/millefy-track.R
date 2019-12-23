#' Visualize read covearge in single-cell RNA-Seq data 
#' 
#' @param track_data A list of tracks.
#' @param track_type A list of track types. Track types are: "sc", "bed", "add" (bulk NGS), "avg", "gene", "title", "axis".
#' @param heights A list of track heights. Or, you can use a unit (e.g., `unit(c(1,1,12,2,1), c("null", "cm", "null", "null", "null")`).
#' @param chr A string. Chromosome name.
#' @param start An integer. Start position.
#' @param end An integer. End position.
#' @param sc_avg A logical. If TRUE (defalut), a track for averaged read coverage for every group is generated.
#' @param sc_avg_height A number. The height of the averaged read coverage track.
#' @param sc_sort_destiny 'none' (default) or 'all' or 'group'.
#' @param title A string. Title.
#' @param binsize A integer. By default bin size is automatically determined so that the number of bins is 1000.
#' @param axis A logical. If TRUE (default), axis for genomic coordinate is shown.
#' @param axis_height A number. The height of the axis track.
#' @param sc_average_mode A string. "mean" (default) or "median". How to summarise single-cell read coverage across samples in every group.
#' @param sc_avg_scale A number. Maximum value of the averaged read coverage track.
#' @param sc_avg_log A logical. If TRUE (default is FALSE), the values in the averaged read coverage track is log-transformed.
#' @param sc_type A string. "heatmap" (default) or "coverage".
#' @return Description of return values
#' @export
#'
#' @examples
#' # Path to bigWig files
#' bwfiles = Sys.glob(file.path(system.file("extdata", package="millefy"), "*.bw"))
#' 
#' # Group labels for bigWig files (same length as \\code{bwfiles})
#' groups = c("00h", "00h", "00h", "12h", "12h", "12h")
#' 
#' # Color labels for bigWig files (A named vector with the same length as the number of kinds of \\code{groups})
#' color_labels <- colorRampPalette(c("yellow", "red"))(length(unique(groups))+1)[1:length(unique(groups))]
#' names(color_labels)  <- unique(groups)
#' 
#' # Parameters
#' max_value = 7873
#' 
#' # Single cell track
#' scTrackBw <- list(path_bam_files = bwfiles, groups = groups, group_colors = color_labels, max_value = max_value, isBw=TRUE)
#' 
#' # Gene annotation track (For faster performance, try to use \\code{dt_gtf} paramter)
#' path_gtf = system.file("extdata", "example.gtf", package="millefy")
#' dt_gtf_exon <- gtfToDtExon(path_gtf)
#' geneTrack1 <- list(path_gtf = path_gtf, dt_gtf = dt_gtf_exon, label = "GENCODE")
#' 
#' # Prepare arguments for \\code{millefyPlot()}
#' tdlist <- list(scTrackBw, geneTrack1)
#' tt <- c("sc", "gene")
#' heights = c(12, 2)
#' text_main = "My plot"
#' 
#' # Location to visualize
#' chr =  "chr19" # character
#' start = 5824708 # integer
#' end = 5845478 # integer
#' 
#' ########
#' # Plot #
#' ########
#' # Plot
#' # When we don't set the sc_sort_destiny parameter (default), the order of single cells is the order of bwfiles.
#' l <- millefyPlot(track_data=tdlist, track_type=tt, heights=heights,
#'           sc_type = "heatmap",
#'           chr = chr, start = start, end = end,
#'           sc_avg = TRUE, sc_avg_height = 1,
#'           title = text_main)
#' 
#' 
#' # Replot
#' # When we set sc_sort_destiny = 'all', all single cells are reordered by diffusion maps.
#' invisible(
#'   millefyPlot(
#'         track_data=l$track_data, track_type=l$track_type, heights=l$heights,
#'         sc_type = "heatmap",
#'         chr = chr, start = start, end = end,
#'         sc_avg = TRUE, sc_avg_height = 1,
#'         title = text_main, sc_avg_scale = 10, sc_sort_destiny = 'all'
#'         )
#' )
#' 
#' # Replot
#' # When we set sc_sort_destiny = 'group', all single cells in each group are reordered by diffusion maps.
#' invisible(
#'   millefyPlot(
#'         track_data=l$track_data, track_type=l$track_type, heights=l$heights,
#'         sc_type = "heatmap",
#'         chr = chr, start = start, end = end,
#'         sc_avg = TRUE, sc_avg_height = 1,
#'         title = text_main,
#'         sc_avg_scale = 10, sc_sort_destiny = 'group'
#'         )
#' )

#' 

millefyPlot <- function(
                          track_data, track_type, heights, sc_type = c("coverage", "heatmap"), 
                          chr, start, end, binsize, title, axis = TRUE, axis_height = 1, 
                          sc_avg = TRUE, sc_avg_height = 1, sc_avg_scale=NA, sc_avg_log = FALSE, 
                          sc_average_mode = c("mean", "median"), sc_sort_destiny = c('none', 'all', 'group')
                          ){

  print(sprintf("Begin millefyPlot: %s", Sys.time()))
   
  # track_data: list of list  scale, path	(avg == TRUE, avg_height, avg_scale)
  # track_type character vector {axis, bed, title, sc, avg, add, gene}
  # heights: example: unit(c(1,1,12,2,1), c("null", "cm", "null", "null", "null")
  
  sc_type <- match.arg(sc_type, c("heatmap", "coverage"))
  sc_average_mode <- match.arg(sc_average_mode, c("mean", "median"))
  sc_sort_destiny <- match.arg(sc_sort_destiny, c('none', 'all', 'group'))
  

  if(sc_avg){
    if(!any("avg" == track_type)){
      for(i in seq_along(track_type)){
        if(track_type[i] == "sc"){
          i_after <- i
          track_data <- append(track_data, list( makeAvgTrack(track_data[[i_after]], sc_avg_scale)), after = i_after)
          track_type <- append(track_type, "avg", after = i_after)
          heights <- append(heights, calcScAvgHeight(track_data[[i_after]], sc_avg_height), after = i_after)
          break
        }
      }
    }
  }
  
  for(i in seq_along(track_type)){
    if(track_type[i] == "add"){
      if(!is.null(track_data[[i]]$trackHeight)){
        heights[i] <- track_data[[i]]$trackHeight * length(track_data[[i]]$path_bam_files)
      }
    }
  }
  
  if(!missing(title) & !any(track_type == "title")){
    track_data <- append(track_data, list(main=title), after = 0)
    track_type <- append(track_type, "title", after = 0)
    heights <- append(heights, 1, after = 0)
  }
  
  if(axis & !any(track_type == "axis")){
    track_data <- append(track_data, list("axis"))
    track_type <- append(track_type, "axis")
    heights <- append(heights, axis_height)
  }
  
  
  select <- GRanges(seqnames = chr, 
                    ranges = IRanges(start = start, end = end))
  
  if(missing(binsize)){
    nbin <- min(1000, end - start + 1)
    binsize <- ceiling((end - start + 1)/nbin)
    nbin <- calcNbin(start, end, binsize)
  }else{
    nbin <- calcNbin(start, end, binsize)
  }
  
  print(track_type)
  
  grid.newpage()
  pushViewport(plotViewport(c(1,0,1,0)))
  pushViewport(viewport(layout = grid.layout(length(track_data), 1, heights = heights)))
  
  for(i in seq_along(track_data)){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    
    if(track_type[i] == "sc"){
      if(sc_type == "heatmap"){
        track_data[[i]] <- plotScHeatmapTrack(track_data[[i]], select, nbin, binsize, sc_avg, sc_avg_log, sc_average_mode, sc_sort_destiny)
        if(sc_avg){
          # print(str(track_data[[i+1]]))
          if(is.na(sc_avg_scale) && !is.null(track_data[[i+1]]$max_value)){ # Danger!!
            sc_avg_scale <- track_data[[i+1]]$max_value
          }
          track_data[[i+1]] <- runMakeAggregationMatrix(track_data[[i]], select, nbin, binsize, sc_avg, sc_avg_log, sc_average_mode)
          if(!is.na(sc_avg_scale)){
            track_data[[i+1]]$max_value <- sc_avg_scale
          }
        }
      }else{
        track_data[[i]]$mat <- plotScCoverageTrack(track_data[[i]], select, nbin, binsize, sc_avg)
        if(sc_avg){
          track_data[[i+1]] <- runMakeAggregationMatrix(track_data[[i]], select, nbin, binsize, sc_avg)
          if(!is.na(sc_avg_scale)){
            track_data[[i+1]]$max_value <- sc_avg_scale
          }
        }
      }
    }else if(track_type[i]=="avg"){
      track_data[[i]] = plotCoverageTrack(track_data[[i]], select, nbin, binsize)
    }else if(track_type[i]=="add"){
      track_data[[i]] <- plotCoverageTrack(track_data[[i]], select, nbin, binsize)
    }else if(track_type[i]=="gene"){
      plotGeneTrack(track_data[[i]], select)
    }else if(track_type[i]=="title"){
      plotTitleTrack(track_data[[i]])
    }else if(track_type[i]=="axis"){
      plotAxisTrack(select)
    }else if(track_type[i]=="bed"){
      plotBedTrack(track_data[[i]], select)
    }
    
    popViewport()    
  }
  
  popViewport()
  popViewport()    
  
  upViewport(0)
  
  print(sprintf("Finished millefyPlot: %s", Sys.time()))
  
  millefy_plot <- list(track_data = track_data, track_type = track_type, heights = heights, 
    plot_settings = list( # Save plot setting
      sc_type=sc_type, chr=chr, start=start, end=end, binsize=binsize, title=tile, axis=axis, axis_height=axis_height,
      sc_avg=sc_avg, sc_avg_height=sc_avg_height, sc_avg_scale=sc_avg_scale, sc_avg_log=sc_avg_log,
      sc_average_mode, sc_sort_destiny
    )
  )
  set_millefy_last_plot(millefy_plot)
  return(millefy_plot)
}


# For backward compatibility
millefyPlot4 <- millefyPlot


plotBedTrack <- function(track, select){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  chr <- levels(seqnames(select))
  start <- start(select)
  end <- end(select)
  
  pushViewport(viewport(layout = grid.layout(1, length(col_ratio), widths=col_ratio)))
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  plotBed(track$path_bed, chr, start, end, track$dt_bed)
  grid.lines(x = c(0,1), y = c(0,0))
  grid.lines(x = c(0,1), y = c(1,1), gp = gpar(col="grey", lwd = 0.5))
  popViewport()
  
  if(!is.null(track$label)){
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    plotGroupLabel(track$label)
    popViewport()
  }
  
  popViewport()
}


plotTitleTrack <- function(main){
  grid.text(main)
}


plotGeneTrack <- function(track, select){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  
  pushViewport(viewport(layout = grid.layout(1, length(col_ratio), widths=col_ratio)))
  
  chr <- levels(seqnames(select))
  start <- start(select)
  end <- end(select)
  
  if(is.null(track$show_transcript_id)){
    track$show_transcript_id = TRUE
  }
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  if(is.null(track$dt_gtf)){
    plotGeneModels4(track$path_gtf, chr, start, end, show_transcript_id = track$show_transcript_id)     	
  }else{
    plotGeneModels4(track$path_gtf, chr, start, end, track$dt_gtf, show_transcript_id = track$show_transcript_id)
  }
  grid.lines(x = c(0,1), y = c(0,0), gp = gpar(col = "grey", lwd = 0.5))
  # grid.lines(x = c(0,1), y = c(1,1))
  popViewport()
  
  if(!is.null(track$label)){
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    plotGroupLabel(track$label)
    popViewport()
  }
  
  popViewport()
}



plotAxisTrack <- function(select){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  
  chr <- levels(seqnames(select))
  start <- start(select)
  end <- end(select)
  
  pushViewport(viewport(layout = grid.layout(1, length(col_ratio), widths=col_ratio)))
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  pushViewport(viewport(xscale = c(start, end), clip = "off"))
  grid.xaxis(gp = gpar(cex = 0.5))
  
  grid.text(sprintf("%s:%d-%d", chr, start, end),
            x = 0.5,
            y= 0.5,  
            hjust = 0.5,
            gp = gpar(cex = 0.7))
  
  popViewport()
  popViewport()
  
  popViewport()
}



plotCoverageTrack <- function(track, select, nbin, binsize, max_value){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  
  #     return(list(mat = mat_a, groups = groups_a, color_labels = color_labels_a))
  groups <- factor(track$groups)
 
  if(is.null(track$mat)){
    if( (!is.null(track$isBw)) && track$isBw==TRUE ){
      mat <- makeCoverageMatrixBw(track$path_bam_files, select, nbin, binsize)
    }else{
      mat <- makeCoverageMatrix(track$path_bam_files, select, nbin, binsize, track$normFactor)
    }
  }else{
    mat <- track$mat
  }
  mat[is.na(mat)] <- 0 # OK?
  
  mat2 <- mat
  
  if(! is.null(track$log) && track$log == TRUE){
    mat2 <- log10(mat2+1)
    # mat[mat == 0] <- NA
    # mat <- log10(mat) #log10(mat + min(mat,na.rm=T))
    # mat[is.na(mat)] <- 0
  }
  
  if(is.null(track$max_value)){
    max_value <- max(mat2, na.rm = T)
    if(!is.null(track$ignore_outlier) && track$ignore_outlier == TRUE){
      max_value <- maxIgnoreOutlier(mat2)
    }
  }
  
  if(missing(max_value)){
    if(is.null(track$max_value)){
      max_value <- max(mat2, na.rm = T)
    }else{
      max_value <- track$max_value
    }
  }
  if(max_value == 0){
    max_value <- 1
  }

  min_value <- 0 # min(mat2, na.rm = T)
  
  pushViewport(viewport(layout = grid.layout(1, length(col_ratio), widths=col_ratio)))
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  plotCoverage2(mat2, max_value)
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
  plotYvalue(nrow(mat2), max_value, min_value)
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  if(is.null(track$color_labels)){
    plotColorLabelNoColor(unique(groups)) 
  }else{
    color_labels <- track$color_labels
    plotColorLabel(rev(color_labels))
  }
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  plotGroupLabel(groups)
  popViewport()
  
  popViewport() 
  
  
  track$mat = mat
  track$max_value = max_value
  
  return(track)
  

}


runMakeAggregationMatrix <- function(track, select, nbin, binsize, sc_avg, sc_avg_log, sc_average_mode = "mean"){
  groups <- factor(track$groups)
  groups <- droplevels(groups)
  
  group_colors <- track$group_colors

  e <- try({
    isAppropriateColorLabels(group_colors, groups)
    }, silent=TRUE)
  if(class(e) == "try-error"){
    cat("Group colors must be a named vector.\n")
  }

  color_labels <- makeColorLabels(group_colors, groups)
  
  mat <- track$mat
  
  mat[is.na(mat)] <- 0 # OK?
  
  mat2 <- mat
  
#   if(is.null(track$max_value)){
#     if(!is.null(track$ignore_outlier) && track$ignore_outlier == TRUE){
#       max_value <- maxIgnoreOutlier(mat2)
#     }else{
#       max_value <- max(mat2, na.rm = T)
#     }
#   }else{
#     max_value <- track$max_value
#   }
#   if(max_value == 0){
#     max_value <- 1
#   }  
  
  if(sc_avg){
    mat_a <- makeAggregationMatrix(mat, groups, average_mode = sc_average_mode)
    groups_a <- makeAggregationGroups(groups, prefix = "Avg ", suffix = "")
    color_labels_a <- makeAggregationColorLabels(group_colors, groups)
    if(sc_avg_log){
      return(list(mat = mat_a, groups = groups_a, color_labels = color_labels_a, log = TRUE))
      # return(list(mat = mat_a, groups = groups_a, color_labels = color_labels_a, max_value=max_value, log = TRUE))
    }else{
      return(list(mat = mat_a, groups = groups_a, color_labels = color_labels_a))
      # return(list(mat = mat_a, groups = groups_a, color_labels = color_labels_a, max_value=max_value))
    }
  }
}


plotScHeatmapTrack <- function(track, select, nbin, binsize, sc_avg, sc_avg_log, sc_average_mode = "mean", sc_sort_destiny = 'none'){
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  
  breaks <- 100
  pal <- colorRampPalette(c('white','blue'))
  pal <- pal(breaks)
  
  groups <- factor(track$groups)
  groups <- droplevels(groups)
  
  group_colors <- track$group_colors

  e <- try({
    isAppropriateColorLabels(group_colors, groups)
    }, silent=TRUE)
  if(class(e) == "try-error"){
    cat("Group colors must be a named vector.\n")
  }

  color_labels <- makeColorLabels(group_colors, groups)
  
  if(is.null(track$mat)){
    if( (!is.null(track$isBw)) && track$isBw==TRUE ){
        mat <- makeCoverageMatrixBw(track$path_bam_files, select, nbin, binsize)  
      }else{
        mat <- makeCoverageMatrix(track$path_bam_files, select, nbin, binsize, track$normFactor)  
      }
  }else{
    mat <- track$mat
  }
  mat[is.na(mat)] <- 0 # OK?
  
  mat2 <- mat
  
  if(is.null(track$max_value)){
    if(!is.null(track$ignore_outlier) && track$ignore_outlier == TRUE){
      max_value <- maxIgnoreOutlier(mat2)
    }else{
      max_value <- max(mat2, na.rm = T)
    }
  }else{
    max_value <- track$max_value
  }
  if(max_value == 0){
    max_value <- 1
  }
  
  
  min_value <- 0 # min(mat2, na.rm = T)
  
  if(is.null(sc_sort_destiny) || sc_sort_destiny == 'none' || sum(mat2) == 0){
    color_labels2 <- color_labels
  }else{
    sel_danger <- rowSums(mat2) == 0
    sel_safe <- !(sel_danger)
    sel_safe[which(sel_danger)[1]] <- TRUE
    dc1_value <- runDiffusionMap(mat2[sel_safe,]) # dif <- DiffusionMap(as.ExpressionSet(as.data.frame(mat[sel_safe,])))
    dc1 <- numeric(nrow(mat2))
    dc1[sel_safe] <- dc1_value
    dc1[!sel_safe] <- dc1[which(sel_danger)[1]]
    
    
    if(sc_sort_destiny == 'all'){
      neworder <- order(dc1, decreasing = T)
      mat2 <- mat2[neworder,]
      color_labels2 <- color_labels[neworder]
    }else if(sc_sort_destiny == 'group'){
      neworder <- orderWithinGroups(dc1, groups)
      mat2 <- mat2[neworder,]
      color_labels2 <- color_labels[neworder]
    }else{
      message("'sort_destiny' must be set to either of 'all', 'group', 'none'")
      color_labels2 <- color_labels
    }
  }
  
  mat2[mat2 > max_value] <- max_value
  
  pushViewport(plotViewport(c(0.5,0,0.5,0)))
  
  pushViewport(viewport(layout = grid.layout(2, length(col_ratio), heights = unit(c(2,12), c("lines", "null") ), widths=col_ratio)))
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  plotHeatmap2(mat2, breaks, pal, max_value)
  popViewport()     
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  plotColorKey(pal[1:breaks], min_value, max_value)
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  plotColorLabel(rev(color_labels2))
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  if(sc_sort_destiny == 'group'){
    plotGroupLabel(groups)
  }
  popViewport()
  
  popViewport() 
  
  popViewport() 
  
  track$mat = mat
  track$max_value = max_value


  return(track)
#   if(sc_avg){
#     mat_a <- makeAggregationMatrix(mat, groups, average_mode = sc_average_mode)
#     groups_a <- makeAggregationGroups(groups, prefix = "Avg ", suffix = "")
#     color_labels_a <- makeAggregationColorLabels(group_colors, groups)
#     if(sc_avg_log){
#       return(list(mat = mat_a, groups = groups_a, color_labels = color_labels_a, max_value=max_value, log = TRUE))
#     }else{
#       return(list(mat = mat_a, groups = groups_a, color_labels = color_labels_a, max_value=max_value))
#     }
#   }
}


orderWithinGroups <- function(values, groups){
  unlist(as.vector(mapply(function(x,y){y[order(x, decreasing = T)]}, split(values, groups), split(seq_along(values), groups))))
}


makeAvgTrack <- function(track, sc_avg_scale){
  return(list(0))
}

runDiffusionMap<- function(mat){
  #   # res = hclust(dist(mat))$order
  #   # res = hclust(as.dist(cor(t(m1))), method = "median")$order
  #   res = hclust(dist(mat), method = "median")$order
  #   return(res)

  e <- try({
    dif <- DiffusionMap(as.ExpressionSet(as.data.frame(mat)));
    cat(sprintf("Eigenvalue of DC1: %f\n", eigenvalues(dif)[1]));
    res <- dif@eigenvectors[, "DC1"];
    }, silent=TRUE)

  if(class(e) == "try-error"){
    cat("There was a problem when running diffusion map. Trying PCA instead...\n")
    e <- try({
      pca =  prcomp(log10(mat+1),
       center = TRUE,
       scale. = TRUE);
      cat(sprintf("The standard deviations of PC1: %f\n", pca$sdev[1]));
      res <- pca$x[, "PC1"];
      }, silent=TRUE)

    if(class(e) == "try-error"){
      res <- rep(1, nrow(mat))
    }else{
      return(res)
    }
  }else{
    return(res)
  }
}


calcScAvgHeight <- function(track, sc_avg_height){
  groups <- factor(track$groups)
  groups <- droplevels(groups)
  sc_avg_height * length(levels(groups))
}


bedToDt <- function(path_bed){

  dt_bed <- fread(path_bed, header = F)
  setnames(dt_bed, c("V1", "V2", "V3", "V4"), c("chr", "start", "end", "name"))
  dt_bed[, start := start + 1]
  dt_bed
}
  
  
  
plotBed <- function(path_bed, x_chr, x_start, x_end, dt_bed){
  
  if(missing(dt_bed)){
    dt_bed <- bedToDt(path_bed) 
  }
  
  dt_bed_local <- dt_bed %>% subset(chr == x_chr & start <= x_end & end >= x_start)
  
  h <- 0.5
  
  if(nrow(dt_bed_local)>0){
    
    pushViewport(viewport(xscale = c(x_start, x_end), yscale = c(0,1), clip = "on"))
    
    for(i in 1:nrow(dt_bed_local)){
      grid.rect(
        x = unit((dt_bed_local[i,start]-x_start)/(x_end - x_start + 1), "npc"),
        y = 0,
        width = unit((dt_bed_local[i,end]-dt_bed_local[i,start]+1)/(x_end - x_start + 1), "npc"),
        height = h,
        gp = gpar(fill = "#B32E3D"),
        default.units = "native",
        hjust = 0
      )
    }
    if("name" %in% colnames(dt_bed_local)){
      for(i in 1:nrow(dt_bed_local)){
        xpos <- (min(x_end, dt_bed_local[i,end]) + max(x_start, dt_bed_local[i,start]) - 2* x_start)/2/(x_end - x_start + 1)
        # print(xpos)
        grid.text(dt_bed_local[i,name],
                  x = unit(xpos, "npc"), y=unit(h, "native"),
                  hjust = 0, vjust = 0,
                  rot = 15,
                  gp = gpar(cex = 0.3)
            )

     }
    
    popViewport()
    }
  }
}
