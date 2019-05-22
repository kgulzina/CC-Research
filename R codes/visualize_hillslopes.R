library("tidyverse")

set.seed(1)
d = bind_rows(
  data.frame(hillslope = 1,
             soiltype = "A",
             subhillslope = 1:12,
             rise = runif(12),
             run = runif(12),
             stringsAsFactors = FALSE),
  data.frame(hillslope = 2,
             soiltype = "B",
             subhillslope = 1:12,
             rise = runif(12),
             run = runif(12),
             stringsAsFactors = FALSE),
  data.frame(hillslope = 3,
             soiltype = "C",
             subhillslope = 1:12,
             rise = runif(12),
             run = runif(12),
             stringsAsFactors = FALSE)
)


calc_height = function(d) {
  
}