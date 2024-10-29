## METADATA ===============================================================
## Description: Generate volcanic island with 2d-gausians with decay over time
## 
## R version: 4.0.2 for Windows
## Date: 2022-05-30 23:48:22
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##


# LIBRARIES ----------
# lib <- c("raster", "lattice", "rayimage", "viridis")
# sapply(lib, require, character.only = TRUE, quietly = TRUE, warn.conflicts = TRUE)
library(raster)

# FUNCTIONS ----------

# `plotting

plot_persp <- function(x, y, z, lcol=topo.colors, ...) {
  ## getting the value of the midpoint
  zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
  ## calculating the breaks
  breaks <- hist(zz, plot=FALSE)$breaks
  ## cutting up zz
  cols <- lcol(length(breaks)-1)
  zzz <- cut(zz, breaks=breaks, labels=cols)
  ## plotting
  persp(x, y, z, col=as.character(zzz),
        phi=30, theta=-25, ltheta = -70,
        expand = 0.5, border = NA, box = FALSE, shade = 0.75, ...)
  ## return breaks 
  list(breaks=breaks)#, colors=cols) #...and colors for the legend
  
}

plot_multiple_persp <- function(l, lcol=topo.colors, scol="darkblue", ...){
  #l is the list of z values...
  # merge lists for ranges as a large vector..
  lv <- do.call(c, l)
  breaks <- hist(lv, plot=FALSE)$breaks
  ## cutting up zz
  cols <- c(scol, lcol(length(breaks)-2))
  # make sure par exits cleanly
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  times <- length(l)
  par(mfrow=c(1,times))
  for (i in 1:length(l)){
    zi <- l[[i]]
    # zifacet <- zi[-1, -1] + zi[-1, -ncol(zi)] + zi[-nrow(zi), -1] + zi[-nrow(zi), -ncol(zi)]
    zz <- (zi[-1,-1] + zi[-1,-ncol(zi)] + zi[-nrow(zi),-1] + zi[-nrow(zi),-ncol(zi)])/4
    zzz <- cut(zz, breaks=breaks, labels=cols)
    ## plotting
    persp(x, y, zi,
          zlim=c(min(lv, na.rm=T), max(lv, na.rm=T)),
          col=as.character(zzz),
          phi=30, theta=-25, ltheta = -70,
          expand = 0.5, border = NA, box = FALSE, shade = 0.25, ...)
    if(!any(is.na(names(l)))){
      title(names(l)[i])
    }
    
    
    # if last
    if (i == length(l)){
      # [] add legend
      # https://www.rdocumentation.org/packages/fields/versions/13.3/topics/colorbar.plot
    }
  }
}

# 'shape

drop <- function(x,y, a=3, ...) { 
  r <- sqrt(x^2+y^2) 
  z <- a * sin(r)/r 
}

gaussian_2d <- function(x,y, x0, y0, sdx, sdy, a) { 
  r <- (((x-x0)^2)/(2*sdx^2))+(((y-y0)^2)/(2*sdy^2))
  z <- a*exp(-r)
  return(z)
}

# 'generator

random_lanscape <- function(x,y, n=2){
  zf <- 0
  for (i in 1:n){
    zf <- zf+outer(x, y, gaussian_2d, x0=sample(x,1), y0=sample(y,1), sdx=rnorm(1,2,2), sdy=rnorm(1,2,2), a=rnorm(1,2,2))
  }
  return(zf)
}

set_landscape_t <- function(x=seq(-10, 10, length=60),y=NA, t=-4:0, 
                            x0t=function(t){rep(0,length(t))}, 
                            y0t=function(t){rep(0,length(t))}, 
                            sdxt=function(t){2/(abs(t)+1)}, sdyt=NA, 
                            at=function(t){abs(t)+1}){
  if (any(is.na(y))){y=x}
  if (any(is.na(sdyt))){sdyt=sdxt}
  l <- list()
  for (i in 1:length(t)){
    ti <- t[i]
    l[[i]] <- outer(x, y, gaussian_2d, x0=x0t(ti), y0=y0t(ti), sdx=sdxt(ti), sdy=sdyt(ti), a=at(ti))
    names(l)[i] <- paste0("x0=",x0t(ti), "|y0=",y0t(ti), "|sdx=", sdxt(ti),"|sdy=",sdyt(ti),"|a=", at(ti))
  }
  return(l)
}



# PLOTTING ----------

##### Example plots

### single simple -----
x <- seq(-10, 10, length=60)
y <- x
z <- outer(x, y, gaussian_2d, x0=0, y0=0, sdx=2, sdy=2, a=10)
plot_persp(x, y, z)


### single simple drop -----
x <- seq(-10, 10, length=60)
y <- x
z <- outer(x, y, drop, x0=0, y0=0, a=10)
plot_persp(x, y, z)




### multiple -----
x <- seq(-10, 10, length=60)
y <- x
l <- list()
for (ti in 1:4){
  l[[ti]] <- outer(x, y, gaussian_2d, x0=0, y0=0, sdx=ti, sdy=ti, a=ti^-0.5)
  #l[[ti]][l[[ti]]<0.1] <- NA
}
plot_multiple_persp(l, lcol=terrain.colors)


### Two islands  -----

x <- seq(-10, 10, length=60)
y <- x
z1 <- outer(x, y, gaussian_2d, x0=-5, y0=-5, sdx=2, sdy=2, a=10)
plot_persp(x, y, z1)
z2 <- outer(x, y, gaussian_2d, x0=4, y0=5, sdx=3, sdy=4, a=5)
plot_persp(x, y, z2)
zf <- z1+z2 # add z2
plot_persp(x, y, zf)

# create a list through time
l <- list(z1=zf,
          z2=zf+outer(x, y, gaussian_2d, x0=-8, y0=-7, sdx=2, sdy=1, a=4),
          z3=zf+outer(x, y, gaussian_2d, x0=-8, y0=-7, sdx=2, sdy=1, a=13))
plot_multiple_persp(l, lcol=terrain.colors)

# converting to rasters
rb <- brick(lapply(l, raster))
plot(rb)

# considering sea levels 
selev <- c(1,2,1)
ls <- l
for (zi in 1:3) {
  ms <- ls[[zi]]<=selev[zi]
  ls[[zi]][ms] <- NA
}
rb <- brick(lapply(ls, raster))
plot(rb)
lapply(list, function)

## example generate landscape
# single simple
x <- seq(-10, 10, length=60)
y <- x
# few
plot_persp(x, y, random_lanscape(x,y,10))
# many
plot_persp(x, y, random_lanscape(x,y,100), lcol=terrain.colors)


###### Example for in time
l <- set_landscape_t()
plot_multiple_persp(l, lcol=terrain.colors)
