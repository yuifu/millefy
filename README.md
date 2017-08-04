# Millefy
Millefy: Genome browser-like visualization of single-cell RNA-seq dataset

## Example of millefy plot

TODO

## Requirements

- R (version 3.2.2 or higher)

## Dependency

- Rsamtools
- GenomicRanges
- data.table
- dplyr
- magrittr
- grid
- destiny
- rtracklayer
- IRanges
- tidyr

## Installation

```
devtools::install_github("yuifu/millefy")
```

## Usage

### Example

```
# Path to bigWig files
bwfiles = c("bw1.bw", "bw2.bw", "bw3.bw", "bw4.bw", "bw5.bw")
# Group labels for bigWig files (same length as \\code{bwfiles})
groups = c("A", "A", "A", "B", "B")
# Color labels for bigWig files (same length as \\code{bwfiles})
color_labels <- colorRampPalette(c("yellow", "red"))(length(unique(groups))+1)[1:length(unique(groups))]
# Parameters
max_value = 7873
# Single cell track
scTrackBw <- list(path_bam_files = bwfiles, groups = groups, group_colors = color_labels, max_value = max_value, isBw=TRUE)
# Gene annotation track (For faster performance, try to use \\code{dt_gtf} paramter)
path_gtf = "annotation/gencode.gtf"
dt_gtf_exon <- gtfToDtExon(path_gtf)
geneTrack1 <- list(path_gtf = path_gtf, dt_gtf = dt_gtf_exon, label = "GENCODE")
# Prepare arguments for \\code{millefyPlot4()}
tdlist <- list(scTrackBw, geneTrack1)
tt <- c("sc", "gene")
heights = c(12, 2)
text_main = "My plot"

# Plot
l <- millefyPlot4(track_data=tdlist, track_type=tt, heights=heights,
          sc_type = "heatmap",
          chr = chr, start = start, end = end,
          sc_avg = TRUE, sc_avg_height = 1,
          title = text_main)

# Replot
invisible(
  millefyPlot4(
        track_data=l$track_data, track_type=l$track_type, heights=l$heights,
        sc_type = "heatmap",
        chr = chr, start = start, end = end,
        sc_avg = TRUE, sc_avg_height = 1,
        title = text_main, sc_avg_scale = 10, sc_sort_destiny = 'all'
        )
)

# Replot
invisible(
  millefyPlot4(
        track_data=l$track_data, track_type=l$track_type, heights=l$heights,
        sc_type = "heatmap",
        chr = chr, start = start, end = end,
        sc_avg = TRUE, sc_avg_height = 1,
        title = text_main,
        sc_avg_scale = 10, sc_sort_destiny = 'group'
        )
)
```

### Test

```
path_gtf = system.file("extdata", "example.gtf", package="millefy")
bwfiles = Sys.glob(file.path(system.file("extdata", "example_bw", package="millefy"), "*.bw"))
```



## License

Copyright (c) 2017 Haruka Ozaki Released under the MIT license
