% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_diffbind_analyzer.R
\name{make_diffbind_PCA_plot}
\alias{make_diffbind_PCA_plot}
\title{Create and save a PCA plot}
\usage{
make_diffbind_PCA_plot(
  dba_object,
  figure_title_nocontrast,
  figure_title_contrast,
  save_name,
  save_directory = NULL,
  point_size = 8,
  width = 20,
  height = 9,
  color = c("darkmagenta", "aquamarine4")
)
}
\arguments{
\item{dba_object}{DBA object (RDS format) generated by \code{perform_diffbind()}.}

\item{figure_title_nocontrast}{Character. Title for the PCA plot using all consensus
regions.}

\item{figure_title_contrast}{Character. Title for the PCA plot using only differentially
bound regions.}

\item{save_name}{Character. Base name (without extension) for the saved plot files.}

\item{save_directory}{Character. Directory path where output plots will be saved.}

\item{point_size}{Numeric. Size of the points, default is set to 8.}

\item{width}{Numeric. Width dimensions of saved plot in inches. Default is 20.}

\item{height}{Numeric. Height dimensions of saved plot in inches. Default is 9.}

\item{color}{Character vector of length 2. Colors to use for the two biological conditions
in PCA plots. Default is \code{c("darkmagenta", "aquamarine4")}.}
}
\value{
A \code{gridExtra} object containing a side-by-side layout of two PCA plots.
}
\description{
make_diffbind_PCA_plot(dba_object, figure_title_nocontrast, figure_title_contrast,
save_directory, save_name, point_size, width, height, color) Creates side-by-side PCA
plots to compare the variance structure of all consensus binding regions vs differentially
bound (DB) regions from a DBA object. Both plots are saved as PNG and PDF files under
designated file path.
}
\examples{
make_diffbind_PCA_plot(
  dba_object = dba.obj,
  figure_title_nocontrast = "PCA of All Consensus Regions",
  figure_title_contrast = "PCA of Differential Binding Regions",
  save_directory = "diffbind_figures",
  save_name = "diffbind_pca_plots")
}
\keyword{PCA}
\keyword{plot}
