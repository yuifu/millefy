plotHeatmap2 <- function(mat, breaks, pal, fixed_max){
  nr <- nrow(mat)
  grid.raster(matrix(pal[as.numeric(cut2(mat, breaks, fixed_max = fixed_max))], nrow = nr),
              interpolate = F, width =1, height = 1)
}

plotHeatmap <- function(mat, breaks, pal){
  nr <- nrow(mat)
  grid.raster(matrix(pal[as.numeric(cut(mat, breaks))], nrow = nr),
              interpolate = F, width =1, height = 1)
}

plotColorKey <- function(colors, min_color, max_color){
  pushViewport(viewport(layout = grid.layout(
    1,2,
    widths = unit(c(4,1), c("null","null"))
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  

  grid.raster(matrix(colors,nrow=1), y = 0.25, interpolate = F, width = 1, height = 0.5, vjust = 0)
  grid.rect(x=0, y = 0.25, width=1, height = 0.5, gp = gpar(fill = NA, col="black"), vjust = 0, hjust = 0)
  
  grid.text(min_color, x=0, y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.5))
  grid.text(sprintf("%.3f", max_color), x=1, y=0.24, hjust = 1, vjust = 1, gp = gpar(cex = 0.5))
  
  popViewport()
  popViewport()
}



cut2 <- function (x, breaks, labels = NULL, include.lowest = FALSE, right = TRUE, 
                  dig.lab = 3L, ordered_result = FALSE, fixed_max, ...) 
{
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (length(breaks) == 1L) {
    if (is.na(breaks) || breaks < 2L) 
      stop("invalid number of intervals")
    nb <- as.integer(breaks + 1)
    dx <- diff(rx <- range(x, na.rm = TRUE))
    rx[2L] <- fixed_max
    if (dx == 0) {
      dx <- abs(rx[1L])
      breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, 
                        length.out = nb)
    }
    else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + 
                               dx/1000)
    }
  }
  else nb <- length(breaks <- sort.int(as.double(breaks)))
  if (anyDuplicated(breaks)) 
    stop("'breaks' are not unique")
  codes.only <- FALSE
  if (is.null(labels)) {
    for (dig in dig.lab:max(12L, dig.lab)) {
      ch.br <- formatC(0 + breaks, digits = dig, width = 1L)
      if (ok <- all(ch.br[-1L] != ch.br[-nb])) 
        break
    }
    labels <- if (ok) 
      paste0(if (right) 
        "("
        else "[", ch.br[-nb], ",", ch.br[-1L], if (right) 
          "]"
        else ")")
    else paste("Range", seq_len(nb - 1L), sep = "_")
    if (ok && include.lowest) {
      if (right) 
        substr(labels[1L], 1L, 1L) <- "["
      else substring(labels[nb - 1L], nchar(labels[nb - 
                                                     1L], "c")) <- "]"
    }
  }
  else if (is.logical(labels) && !labels) 
    codes.only <- TRUE
  else if (length(labels) != nb - 1L) 
    stop("lengths of 'breaks' and 'labels' differ")
  code <- .bincode(x, breaks, right, include.lowest)
  if (codes.only) 
    code
  else factor(code, seq_along(labels), labels, ordered = ordered_result)
}

# cut2(1:6, 3, fixed_max = 10)