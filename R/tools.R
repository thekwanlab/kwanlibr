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
#' ggsave_vector_raster('figure_1A')
#' ggsave_vector_raster('figure_2A', width=8, height=6, vector_device='svg', raster_device='jpeg')
#' 
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
