#' @importFrom MASS kde2d
#' @export get_density

get_density <- \(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- base::findInterval(x, dens$x)
  iy <- base::findInterval(y, dens$y)
  ii <- base::cbind(ix, iy)
  return(dens$z[ii])
}
