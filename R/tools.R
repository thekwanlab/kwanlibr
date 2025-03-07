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
#' ggsave_vector_raster('figure_2A', width=6, height=6, dpi=600, vector_device='svg', raster_device='jpeg')
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
  if (lower > upper) stop(paste("Lower bound",
                                lower,
                                "cannot be greater than upper bound",
                                upper))
  ifelse(x > upper,
         upper,
         ifelse(x < lower, lower, x))
}

#' draw_PCA()
#'
#' run PCA and draw PCA plot
#'
#' @param df A data-matrix or data-frame
#' @param label_list A vector that contains assigned groups of the sample, the 
#' length of the vector should be the same as the number of samples in the df 
#' @param color_list A vector that contains colors of the choice, the length of 
#' the vector should be equal to the number of unique elements in label_list. By 
#' default, use ggplot default colors
#' @param legend Boolean values, if True, draw the legend on the PCA, if false, 
#' no legend. Default is set to TRUE
#' 
#' @import ggplot2
#' 
#' @export
#' 
#' @examples
#' draw_PCA(iris[-5], label_list = iris$Species, color_list = c("green", "blue", "yellow"))
#' draw_PCA(iris[-5], label_list = iris$Species)
#' draw_PCA(iris[-5], legend=FALSE)
#' draw_PCA(iris[-5]) if no label_list, no legend is produced

draw_PCA <- function(df, 
                     label_list=NULL, 
                     color_list=NULL, 
                     legend=TRUE){
  #Check if df has non-numerical columns
  if (!all(sapply(df, is.numeric))) {
    return("Check the structure of the df, make sure df only contains numerical inputs.")
  }
  
  # Validate length of label_list if provided
  if (!is.null(label_list) && length(label_list) != nrow(df)) {
    return("Check the label_list, make sure the length matches the number of samples in df.")
  }
  
  # Validate color_list if label_list is provided
  if (!is.null(label_list) && !is.null(color_list) && length(color_list) != length(unique(label_list))) {
    return("Check the color_list, make sure the number of distinct colors matches 
           the number of unique labels in label_list.")
  }
  pca_result <- prcomp(df,center = TRUE)
  
  #Build a data frame
  pcaData <- as.data.frame(pca_result$x)
  
  if (!is.null(label_list)) {
    pcaData$Group <- label_list # add group to df
  }
  
  # Calculate % Variance Explained
  PoV <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  VoPC1 <- sprintf("%.2f%%", PoV[1] * 100)
  VoPC2 <- sprintf("%.2f%%", PoV[2] * 100)
  
  # Create the plot
  if (is.null(label_list)) {
    p <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
      geom_point(size = 2)
  } else {
    p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 2)
  }
  p <- p + xlab(paste0("PC1: ", VoPC1)) +
    ylab(paste0("PC2: ", VoPC2)) +
    theme_minimal()
  
  if (!is.null(color_list)) {
    p <- p + scale_color_manual(values=color_list)
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}
