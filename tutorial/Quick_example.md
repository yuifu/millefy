# Quick example

## Prepare tracks and parameters

```
# Path to bigWig files
bwfiles = Sys.glob(file.path(system.file("extdata", package="millefy"), "*.bw"))

# Group labels for bigWig files (same length as \\code{bwfiles})
groups = c("00h", "00h", "00h", "12h", "12h", "12h")

# Color labels for bigWig files (A named vector with the same length as the number of kinds of \\code{groups})
color_labels <- colorRampPalette(c("yellow", "red"))(length(unique(groups))+1)[1:length(unique(groups))]
names(color_labels)  <- unique(groups)

# Parameters
max_value = 7873

# Single cell track
scTrackBw <- list(path_bam_files = bwfiles, groups = groups, group_colors = color_labels, max_value = max_value, isBw=TRUE)

# Gene annotation track (For faster performance, try to use \\code{dt_gtf} paramter)
path_gtf = system.file("extdata", "example.gtf", package="millefy")
dt_gtf_exon <- gtfToDtExon(path_gtf)
geneTrack1 <- list(path_gtf = path_gtf, dt_gtf = dt_gtf_exon, label = "GENCODE")

# Prepare arguments for \\code{millefyPlot()}
tdlist <- list(scTrackBw, geneTrack1)
tt <- c("sc", "gene")
heights = c(12, 2)
text_main = "My plot"

# Location to visualize
chr =  "chr19" # character
start = 5824708 # integer
end = 5845478 # integer

```

## Plot

When we don't set the `sc_sort_destiny` parameter (default), the order of single cells is the order of `bwfiles`.

```
# Plot
l <- millefyPlot(track_data=tdlist, track_type=tt, heights=heights,
          sc_type = "heatmap",
          chr = chr, start = start, end = end,
          sc_avg = TRUE, sc_avg_height = 1,
          title = text_main)
```
<img src="../img/millefy_plot_example_default.png" width="40%" />

<!-- ## Replot using `millefy_adjust()`
You can adjust the Millefy plot using `millefy_adjust()`.

For example, when we set `sc_sort_destiny = 'all'`, all single cells are reordered by diffusion maps.

```
# Replot
millefy_adjust(sc_avg_scale = 10, sc_sort_destiny = 'all')
```


When we set `sc_sort_destiny = 'group'`, all single cells in each group are reordered by diffusion maps.

 -->
When we set `sc_sort_destiny = 'all'`, all single cells are reordered by diffusion maps.

```
l <- millefyPlot(track_data=tdlist, track_type=tt, heights=heights,
          sc_type = "heatmap",
          chr = chr, start = start, end = end,
          sc_avg = TRUE, sc_avg_height = 1,
          title = text_main)
```

<img src="../img/millefy_plot_example_all.png" width="40%" />

You can do the same thing using `millefy_adjust()`. Use of `millefy_adjust()` is faster since it reuses read coverage data in the last Millefy plot.

```
# Replot
millefy_adjust(sc_sort_destiny = 'all')
```

When we set `sc_sort_destiny = 'group'`, all single cells in each group are reordered by diffusion maps.

```
# Replot
millefy_adjust(sc_sort_destiny = 'group')
```

<img src="../img/millefy_plot_example_group.png" width="40%" />
