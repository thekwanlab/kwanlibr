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
#' @param data A data-matrix or data-frame, where the rows of the data should be
#' observations, and columns should be numerical variables
#' @param label A vector that contains assigned groups of the sample, the
#' length of the vector should be the same as the number of samples in the df
#' @param color A vector that contains colors of the choice, the length of
#' the vector should be equal to the number of unique elements in label. By
#' default, use ggplot default colors
#' @param legend Boolean values, if True, draw the legend on the PCA, if false,
#' no legend. Default is set to TRUE
#' @return A ggplot object of PCA
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' draw_PCA(iris[-5], label = iris$Species, color = c("green", "blue", "yellow"))
#' draw_PCA(iris[-5], label_list = iris$Species)
#' draw_PCA(iris[-5], legend=FALSE)
#' draw_PCA(iris[-5]) # if no label_list, no legend is produced

draw_PCA <- function(data,
                     label=NULL,
                     color=NULL,
                     legend=TRUE){
  #Check if data has non-numerical columns
  if (!all(sapply(data, is.numeric))) {
    stop("Check the structure of the data, make sure data only contains numerical inputs.")
  }

  # Validate length of label if provided
  if (!is.null(label) && length(label) != nrow(data)) {
    stop("Check the label, make sure the length matches the number of samples in data.")
  }

  # Validate color if label is provided
  if (!is.null(label) && !is.null(color) && length(color) != length(unique(label))) {
    stop("Check the color, make sure the number of distinct colors matches the
         number of unique labels in label.")
  }
  pca_result <- prcomp(data,center = TRUE)

  #Build a data frame
  pcaData <- as.data.frame(pca_result$x)

  if (!is.null(label)) {
    pcaData$Group <- label # add group to df
  }

  # Calculate % Variance Explained
  PoV <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  VoPC1 <- sprintf("%.2f%%", PoV[1] * 100)
  VoPC2 <- sprintf("%.2f%%", PoV[2] * 100)

  # Create the plot
  if (is.null(label)) {
    p <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
      geom_point(size = 2)
  } else {
    p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 2)
  }
  p <- p + xlab(paste0("PC1: ", VoPC1)) +
    ylab(paste0("PC2: ", VoPC2)) +
    theme_minimal()

  if (!is.null(color)) {
    p <- p + scale_color_manual(values=color)
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

#' Create and return a volcano plot
#'
#' draw_volcano_general(data, colors, size, alpha, criteria, xdiff, ymax) creates
#' and returns a volcano plot. This function allows flexible coloring based on a
#' user-defined grouping column in the dataset or applies a uniform color if no
#' grouping is specified.
#'
#' @param data A \code{data.frame} containing the input data. Must include at least
#' the columns \code{Fold} (x-axis) and \code{negLogFDR} (y-axis).
#' @param colors Character. Either:
#'   \itemize{
#'     \item A single unnamed color string if \code{criteria} is \code{NULL} that is
#'     used for all points.
#'     \item A named character vector of colors, where names correspond to grouping
#'     labels in the specified
#'     \code{criteria} column. For example: \code{setNames(c("grey", "tomato"),
#'     c("NS", "Sig"))}.
#'     }
#' @param size Numeric. Point size in the plot. Default is \code{1}.
#' @param alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque). Default is \code{1}.
#' @param criteria Optional. Character string. The name of a column in \code{data}
#' that defines groupings for coloring. If \code{NULL}, all points use the same color.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#'
#' @return A \code{ggplot} object representing the volcano plot.
#' @export
#' @examples
#' draw_volcano_general(data=data, colors = "grey")

draw_volcano_general <- function(
    data,
    colors,
    size=1,
    alpha=1,
    criteria=NULL,
    xdiff = 5,
    ymax = 40
) {
  if (is.null(data) || !is.data.frame(data)) {
    stop("`data` must be a valid data frame.")
  }

  p <- ggplot()

  if (is.null(criteria)) {
    if (length(colors) != 1 || !is.character(colors)) {
      stop("When `criteria` is NULL, `colors` must be a single unnamed color string.")
    }
    p <- p + geom_point(
      data=data,
      aes(x=Fold, y=negLogFDR),
      color=colors,
      size=size,
      alpha=alpha
    )
  } else {
    if (!(criteria %in% colnames(data))) {
      stop(paste0("Column '", criteria, "' not found in the data."))
    }

    group_levels <- names(colors)
    if (is.null(group_levels) || any(group_levels == "")) {
      stop("`colors` must be a named vector with group names when `criteria` is used.")
    }
    for (grp in group_levels) {
      p <- p + geom_point(
        data = dplyr::filter(data, .data[[criteria]] == grp),
        aes(x = Fold, y = negLogFDR, color = .data[[criteria]]),
        size = size,
        alpha = alpha
      )
    }
    p <- p + scale_color_manual(values = colors)
  }

  # plot styling
  p <- p +
    xlim(-xdiff, xdiff) +
    ylim(0, ymax) +
    labs(
      x = TeX('$\\log_2$ FC'),
      y = TeX('$-\\log_{10}$ FDR')) +
    theme(
      aspect.ratio = 1,
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )

  return(p)
}
