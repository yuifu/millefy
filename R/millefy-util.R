coverageIntoFixedBins <- function(cov, binsize){
  # cat(sprintf("length(cov): %d", length(cov))) # debug
  # cat(sprintf("binsize: %d\n", binsize)) # debug
  # cat(sprintf("nbin: %d\n", length(seq(1,length(cov), by = binsize)))) # debug
  # sapply(seq(1,length(cov), by = binsize), calcPartialMeanOfCoverage)
  sapply(seq(1,length(cov), by = binsize), function(x){
    mean(cov[x:min((x+binsize-1),length(cov))], na.rm=T)
  })
}


calcNbin <- function(start, end, binsize){
  ceiling((end - start+1)/binsize)
}

plotColorLabelNoColor <- function(groups){
  
  pushViewport(viewport(yscale = c(1, length(groups)+1)))
  
  for(i in seq_along(groups)){
    grid.rect(x=unit(0, "npc"), y=i, width = unit(1, "npc"), height = 1, 
              default.units = "native", vjust = 0, hjust = 0,
              gp = gpar(fill = NA, col = "grey"))
  }
  
  popViewport()
}

plotColorLabel <- function(color_labels){
  
  pushViewport(viewport(yscale = c(1, length(color_labels)+1)))
  
  for(i in seq_along(color_labels)){
    grid.rect(x=unit(0, "npc"), y=i, width = unit(1, "npc"), height = 1, 
              default.units = "native", vjust = 0, hjust = 0,
              gp = gpar(fill = color_labels[i], col = NA))
  }
  
  popViewport()
}


plotGroupLabel <- function(groups){
  
  lv <- levels(factor(groups))
  
  # pushViewport(viewport(yscale = c(1, length(lv)+1)))
  pushViewport(viewport(layout = grid.layout(
    length(lv),1,
    heights = unit(sapply(lv, function(x){table(factor(groups))[x]}), "null"),
    widths = unit(1, "null")
  )))
  
  for(i in seq_along(lv)){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))    
    
    grid.text(lv[i],
              x = 0.9,
              y= 0.5,  
              hjust = 1,
              gp = gpar(cex = 0.7))
    
    popViewport()
  }
  
  popViewport()
}

# Maybe unnessesary
addAggregationRows <- function(mat, groups){
  levels_genuine <- levels(groups)[levels(groups) %in% groups]
  
  rbind(mat,
        tapply(1:nrow(mat), groups, function(x){
          
          if(length(x) == 1){
            return(mat[x,])
          }
          colMeans(mat[x,], na.rm = T)
          
        }) %>% unlist %>% matrix(nrow = length(levels_genuine), byrow = T)
  )
}

makeAggregationMatrix <- function(mat, groups, average_mode){
  levels_genuine <- levels(groups)[levels(groups) %in% groups]
  if(average_mode == "mean"){
    tapply(1:nrow(mat), groups, function(x){
      
      if(length(x) == 1){
        return(mat[x,])
      }
      colMeans(mat[x,], na.rm = T)
      
    }) %>% unlist %>% matrix(nrow = length(levels_genuine), byrow = T)
  }else{
    tapply(1:nrow(mat), groups, function(x){
      
      if(length(x) == 1){
        return(mat[x,])
      }
      apply(mat[x,], 2, function(y){median(y,na.rm=TRUE)})
      
    }) %>% unlist %>% matrix(nrow = length(levels_genuine), byrow = T)
  }
}

makeAggregationGroups <- function(groups, prefix = "Avg ", suffix = ""){
  levels_genuine <- levels(groups)[levels(groups) %in% groups]
  factor(paste0(prefix, levels_genuine, suffix), levels = paste0(prefix, levels_genuine, suffix))
}

# Maybe unnessesary
addAggregationGroups <- function(groups, prefix = "Avg ", suffix = ""){
  levels_genuine <- levels(groups)[levels(groups) %in% groups]
  factor(c(as.character(groups), paste0(prefix, levels_genuine, suffix)), levels = c(levels(groups), paste0(prefix, levels_genuine, suffix)))
}

makeColorLabels <- function(group_colors, groups){
  group_colors[levels(groups)[groups]]
}

makeAggregationColorLabels <- function(group_colors, groups){
  group_colors[levels(groups)]
}

plotYvalue <- function(nr, max_color, min_color){
  
  pushViewport(viewport(layout = grid.layout(nr,1)))
  for(i in 1:nr){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.text(sprintf("%.3f", max_color),
              x = 0.05, y= 0.95,
              hjust = 0,
              vjust = 1,
              gp = gpar(cex = 0.3)
    )
    grid.text(sprintf("%.3f", min_color),
              x = 0.05, y= 0.05,
              hjust = 0,
              vjust = 0,
              gp = gpar(cex = 0.3)
    )
    popViewport()
  }
  popViewport()
}


maxIgnoreOutlier <- function(mat){
  # i_max <- which.max(colSums(mat, na.rm = T))
  # i_max_of_max <- which.max(mat[,i_max])
  # first_max <- max(mat[,i_max], na.rm = T) 
  # second_max <- max(mat[,i_max][-i_max_of_max], na.rm = T)
  first_max <- max(mat, na.rm = T) 
  second_max <- max(mat[mat<=quantile(mat, 0.999, na.rm = T)], na.rm = T)
  if(first_max > 10*second_max){
    second_max
  }else{
    first_max
  }
}

