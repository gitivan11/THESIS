### This file contains
# PointsSubset

PointsSubset <- function(x, y, subset, ...) {
  if(any(subset))
    points(x[subset], y[subset], ... )
  return(invisible(NULL))
}
