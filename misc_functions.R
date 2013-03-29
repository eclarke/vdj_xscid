orderByCopy <- function(gtsp.df) {
  return(gtsp.df[order(gtsp.df$copy), ])
}