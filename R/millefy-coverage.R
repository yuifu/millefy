plotCoverage2 <- function(mat, max_color){
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  pushViewport(viewport(layout = grid.layout(nr,1)))
  
  for(i in 1:nr){
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1, xscale = c(1,nc), yscale = c(0, max_color), clip = "on"))
    
    grid.polygon(
      x = c(1,1:nc,nc,1),
      y = c(0,mat[i, ],0,0),
      gp=gpar(fill="blue", col = NA, alpha = 0.5),
      default.units = "native"
    )
    
    popViewport()
  }
  
  popViewport()
}


plotCoverage <- function(mat, max_color){
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  pushViewport(viewport(xscale = c(1,nc), yscale = c(0,(nr * max_color))))
  # print(nr * max_color)
  #       grid.yaxis()
  #       grid.xaxis()
  
  
  grid.polygon(
    x = rep(c(1,1:nc,nc,1), nr),
    y = sapply(nr:1, function(x){
      c(0,mat[nr - x + 1, ],0,0) + (x-1) * max_color
    }),
    id = rep(1:nr, each = nc+3),
    gp=gpar(fill="blue", col = NA, alpha = 0.5),
    default.units="native"
  )
  
  popViewport()
}