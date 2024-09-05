### Function to calculate distance to 1:1 line ----
dist_point_line <- function(a, slope, intercept) {
  b = c(1, intercept+slope)
  c = c(-intercept/slope,0)       
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  return(det(m)/sqrt(sum(v1*v1)))
}