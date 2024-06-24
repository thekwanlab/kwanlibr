#' Save a figure to both PNG and PDF
#' 
#' ggsave_vector_raster(filename, ...) is a wrapper to ggplot2::ggsave. Provide 
#' a base file name and this will
#' save the most recently displayed figure to both a raster (default PNG) and
#' vector (default PDF) format in the same location
#' 
#' @param filename The base filename of the images you want to save. DO NOT 
#' include an extension
#' @param width Figure width passed to ggsave. Default NA
#' @param height Figure height passed to ggsave. Default NA
#' @param vector_device Device to use for the vector image. See ggsave 
#' documentation for options. Default 'pdf'
#' @param raster_device Device to use for the raster image. See ggsave 
#' documentation for options. Default 'png'
#' @param ... Additional options passed to ggsave
#' @return None
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#' saved_plot = ggplot2::ggplot(mtcars, aes(cyl, mpg)) + ggplot2::geom_point()
#' saved_plot
#' 
#' ggsave_vector_raster('figure_1A')
#' ggsave_vector_raster('figure_2A', width=8, height=6, vector_device='svg', raster_device='jpeg')
#' 
#' # Pass extra arguments to ggsave() through ...
#' ggsave_vector_raster('figure_1A', plot=saved_plot) 
ggsave_vector_raster <- function(
  filename, 
  width = NA, 
  height = NA, 
  vector_device = 'pdf', 
  raster_device = 'png',
  ...
) {
  
  ggplot2::ggsave(
    filename = paste0(filename, '.', vector_device),
    device = vector_device,
    width = width,
    height = height,
    ...
  )
  ggplot2::ggsave(
    filename = paste0(filename, '.', raster_device),
     device = raster_device,
     width = width,
     height = height,
     ...
  )
}

#' clamp
#' 
#' Bound a value or array of values
#' 
#' @param x Value or array of values to clamp
#' @param lower lower bound of the clamp
#' @param upper upper bound of the clamp
#' 
#' @export
#' 
#' @examples
#' clamp(1:5) # returns c(1,2,3,4,5)
#' clamp(1:5, lower=3) # returns c(3,3,3,4,5)
clamp <- function(
  x,
  lower = -Inf,
  upper = Inf
) {
  if (lower > upper) stop("Lower bound cannot be greater than upper bound")
  ifelse(x > upper, 
         upper,
         ifelse(x < lower, lower, x))
}