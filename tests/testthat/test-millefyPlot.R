context("millefyPlot")

test_that("millefyPlot works.",{

	bwfiles = Sys.glob(file.path(system.file("extdata", package="millefy"), "*.bw"))
	groups = c("00h", "00h", "00h", "12h", "12h", "12h")
	color_labels <- colorRampPalette(c("yellow", "red"))(length(unique(groups))+1)[1:length(unique(groups))]
	names(color_labels)  <- unique(groups)

	max_value = 7873

	scTrackBw <- list(path_bam_files = bwfiles, groups = groups, group_colors = color_labels, max_value = max_value, isBw=TRUE)

	path_gtf = system.file("extdata", "example.gtf", package="millefy")
	dt_gtf_exon <- gtfToDtExon(path_gtf)
	geneTrack1 <- list(path_gtf = path_gtf, dt_gtf = dt_gtf_exon, label = "example genes")

	tdlist <- list(scTrackBw, geneTrack1)
	tt <- c("sc", "gene")
	heights = c(12, 2)
	text_main = "My plot"

	# Defining location to visualize
	chr =  "chr19" # character
	start = 5824708 # integer
	end = 5845478 # integer

	

	l <- millefyPlot(track_data=tdlist, track_type=tt, heights=heights,
          sc_type = "heatmap",
          chr = chr, start = start, end = end,
          sc_avg = TRUE, sc_avg_height = 1,
          title = text_main)
	

	expect_equal(length(bwfiles), nrow(l$track_data[[2]]$mat))
})
