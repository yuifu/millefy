.millefy_plot_store <- function() {
  .millefy_last_plot <- NULL
  
  list(
    get = function() .millefy_last_plot,
    set = function(value) .millefy_last_plot <<- value
  )
}
.millefy_store <- .millefy_plot_store()

#' Set the last Millefy plot to be fetched by millefy_last_plot()
#'
#' @seealso [millefy_last_plot()]
#' @export
#' @keywords internal
set_millefy_last_plot <- function(value) .millefy_store$set(value)


#' Retrieve the last Millefy plot to be modified or created.
#' @export
#' @keywords internal
millefy_last_plot <- function() .millefy_store$get()

#' Adjust the last Millefy plot
#' 
#' @param max_value A number. The maximum value of heatmap color scaling.
#' @param heights A list of track heights. Or, you can use a unit (e.g., `unit(c(1,1,12,2,1), c("null", "cm", "null", "null", "null")`).
#' @param sc_type A string. "heatmap" (default) or "coverage".
#' @param title A string. Title.
#' @param axis A logical. If TRUE (default), axis for genomic coordinate is shown.
#' @param axis_height A number. The height of the axis track.
#' @param sc_avg A logical. If TRUE (defalut), a track for averaged read coverage for every group is generated.
#' @param sc_avg_height A number. The height of the averaged read coverage track.
#' @param sc_avg_scale A number. Maximum value of the averaged read coverage track.
#' @param sc_avg_log A logical. If TRUE (default is FALSE), the values in the averaged read coverage track is log-transformed.
#' @param sc_average_mode A string. "mean" (default) or "median". How to summarise single-cell read coverage across samples in every group.
#' @param sc_sort_destiny 'none' (default) or 'all' or 'group'.
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
#' path_gtf <- system.file("extdata", "example.gtf", package="millefy")
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
#' 
#' # Plot
#' l <- millefyPlot(track_data=tdlist, track_type=tt, heights=heights,
#'           sc_type = "heatmap",
#'           chr = chr, start = start, end = end,
#'           sc_avg = TRUE, sc_avg_height = 1,
#'           title = text_main)
#' 
#' # Replot ajusting max value of the heatmap
#' millefy_adjust(max_value = 100)
millefy_adjust <- function(
        max_value=NA_real_,
        heights=NA,
        sc_type = NA_character_,
        title = NA_character_,
        axis = NA,
        axis_height = NA_real_,
        sc_avg = NA,
        sc_avg_height = NA_real_,
        sc_avg_scale = NA_real_,
        sc_avg_log = NA_real_,
        sc_average_mode = NA,
        sc_sort_destiny = NA
){
    l = millefy_last_plot()
    if(is.null(l)) stop('No Millefy plot available.')

    # Adjust max_value
    if(!is.na(max_value)) l$track_data[[2]]$max_value <- max_value
    # Adjust heights
    if(!is.na(heights) && length(heights)==length(l$heights)) l$heights <- heights
    # Adjust parameters (except for chr, start, end, track_data, track_type)
    if(!is.na(sc_type)) l$plot_settings$sc_type <- sc_type
    if(!is.na(title)) l$plot_settings$title <- title
    if(!is.na(axis)) l$plot_settings$axis <- axis
    if(!is.na(axis_height)) l$plot_settings$axis_height <- axis_height
    if(!is.na(sc_avg)) l$plot_settings$sc_avg <- sc_avg
    if(!is.na(sc_avg_height)) l$plot_settings$sc_avg_height <- sc_avg_height
    if(!is.na(sc_avg_scale)) l$plot_settings$sc_avg_scale <- sc_avg_scale
    if(!is.na(sc_avg_log)) l$plot_settings$sc_avg_log <- sc_avg_log
    if(!is.na(sc_average_mode)) l$plot_settings$sc_average_mode <- sc_average_mode
    if(!is.na(sc_sort_destiny)) l$plot_settings$sc_sort_destiny <- sc_sort_destiny

    # Replot
    invisible(
        millefyPlot(track_data = l$track_data, 
                    track_type = l$track_type, 
                    heights = l$heights, 
                    sc_type = l$plot_settings$sc_type,
                    chr = l$plot_settings$chr,
                    start = l$plot_settings$start,
                    end = l$plot_settings$end,
                    binsize = l$plot_settings$binsize,
                    title = l$plot_settings$tile,
                    axis = l$plot_settings$axis,
                    axis_height = l$plot_settings$axis_height,
                    sc_avg = l$plot_settings$sc_avg,
                    sc_avg_height = l$plot_settings$sc_avg_height,
                    sc_avg_scale = l$plot_settings$sc_avg_scale,
                    sc_avg_log = l$plot_settings$sc_avg_log,
                    sc_average_mode = l$plot_settings$sc_average_mode,
                    sc_sort_destiny = l$plot_settings$sc_sort_destiny
                    )
    )
    
}

