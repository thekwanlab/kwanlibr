% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_diffbind_analyzer.R
\name{merge_diffbind_with_DE}
\alias{merge_diffbind_with_DE}
\title{Merge Differential Binding and Expression Data}
\usage{
merge_diffbind_with_DE(
  dba_object,
  DE_file_path,
  tss_file_path,
  save_directory,
  contrast_number = 1
)
}
\arguments{
\item{dba_object}{DBA object (RDS format) generated by \code{perform_diffbind()}.}

\item{DE_file_path}{Character. Path to the CSV file containing bulk RNA-seq differential
expression results. This file must include the columns \code{gene_name}, \code{logFC},
and \code{FDR}.}

\item{tss_file_path}{Character. Transcription Start Site file path to be used in annotation
with Differential binding sites.}

\item{save_directory}{Character. Directory path for storing intermediate files, including the TSS
BED file, sorted BED file, and nearest gene lookup TSV.}

\item{contrast_number}{Integer. Index indicating which contrast to extract from the
DBA object. Defeault is set to \code{1}.}
}
\value{
A data frame containing merged binding and expression data. Each row corresponds to
a binding site and includes:
\itemize{
\item Site genomic coordinates
\item Site ID
\item Nearest gene name
\item \code{DE.logFC} and \code{DE.FDR} values
\item \code{DB.logFC} and \code{DB.FDR} values
}
}
\description{
merge_diffbind_with_DE(dba_object, DE_file_path, tss_file_path, save_directory, contrast_number)
integrates differential binding sites with differential expression data (from bulk RNA-seq)
by joining them based on their nearest gene TSS and returns a data frame where each row
represents a differential binding site along with its nearest gene's differential expression.
}
\examples{
merge_diffbind_with_DE(
  dba_object = dba.obj,
  DE_file_path = path_to_DE_csv_file,
  tss_file_path = path_to_tss_file,
  save_directory = "diffbind_files")
}
