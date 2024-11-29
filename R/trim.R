trim <- function(x, minimum = 0, maximum = 1){
  max(min(maximum, x), minimum)
}
