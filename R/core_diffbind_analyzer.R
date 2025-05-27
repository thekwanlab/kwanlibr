#=============================================================================
# kwanlibr is a package that contains wrapper functions
# for DiffBind functions.

#---------
# TO-DO
#---------
# - create testing suite?
# - changelog?
# - create runable examples in documentation?
#=============================================================================

#install dependencies

#DiffBind dependencies
#' @import DiffBind
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#' @import latex2exp
#' @import gtools
#' @import tibble

#---------------------
# Functions
#---------------------

#' @title Perform Diffbind Analysis
#'
#' @name perform_diffbind
#'
#' @description
#' perform_diffbind(sample_sheet_path, control_level, min_members, analysis_method,
#' data_type, fdr_threshold, normalization) performs a complete DiffBind analysis
#' workflow using a provided sample sheet, which includes blacklist filtering, counting,
#' normalization, and contrast setup.
#'
#' @param sample_sheet_path Character. The path to the sample CSV file.
#' Normally the file is named 'diffbind_samples.csv' and placed in the current
#' working directory. It must contain the following columns:
#' \itemize{
#'  \item \strong{"Sample_ID"}
#'  \item \strong{"Tissue"}
#'  \item \strong{"Condition"}
#'  \item \strong{"bamReads"}
#'  \item \strong{"Peaks"}
#'  \item \strong{"PeakCaller"}
#'  }
#' @param control_level Character. Reference group name (e.g., "HET", "Ctrl").
#' @param min_members Integer. Minimum number of replicates per condition (default is \code{NULL}).
#' Set to 2 if any condition has only 2 replicates, otherwise set to \code{NULL}.
#' @param analysis_method The analysis method, either \code{DBA_DESEQ2} (default) or \code{DBA_EDGER}.
#' @param data_type Output format of peak data (default is \code{DBA_DATA_FRAME}).
#' Other options are \code{DBA_DATA_GRANGES}, \code{DBA_DATA_RANGEDDATA}.
#' @param fdr_threshold Numeric. FDR cutoff for significance, used in both filtering
#' and visualization (Default is 0.05).
#' @param normalization Normalization method (default is \code{DBA_NORM_LIB}). Other options
#' include \code{DBA_NORM_RLE}, \code{DBA_NORM_TMM}, etc.
#' @return analyzed DBA object
#' @keywords diffbind
#' @export
#' @examples
#' perform_diffbind(
#'   sample_sheet_path = "diffbind_samples.csv",
#'   control_level = "HET",
#'   min_members = 2)

perform_diffbind <- function(
    sample_sheet_path,
    control_level,
    min_members = NULL,
    analysis_method = DBA_DESEQ2,
    data_type = DBA_DATA_FRAME,
    fdr_threshold = 0.05,
    normalization = DBA_NORM_LIB
){
  if (is.null(sample_sheet_path) || !file.exists(sample_sheet_path)) {
    stop("The 'sampleSheetPath' is invalid or the file does not exist.")
  }

  samplesheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)
  required_columns <- c("Sample_ID", "Tissue", "Condition", "bamReads", "Peaks", "PeakCaller")

  missing_columns <- setdiff(required_columns, colnames(samplesheet))

  if (length(missing_columns) > 0) {
    stop("The sample sheet is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  if (!analysis_method %in% c(DBA_DESEQ2, DBA_EDGER)) {
    stop("Invalid 'analysisMethod'. Choose either 'DBA_DESEQ2' or 'DBA_EDGER'.")
  }
  if (!control_level %in% samplesheet$Condition) {
    stop("the specified reference'", control_level, "'is not found in the Condition")
  }

  #Generate dba object
  dba.obj.samples <- DiffBind::dba(sampleSheet = sample_sheet_path)
  dba.obj.samples$config$th <- fdr_threshold
  dba.obj.samples$config$DataType <- data_type
  dba.obj.samples$config$AnalysisMethod <- analysis_method

  dba.obj.blacklist <- DiffBind::dba.blacklist(dba.obj.samples)
  dba.obj.count <- DiffBind::dba.count(dba.obj.blacklist, summits = TRUE)
  dba.obj.normalize <- DiffBind::dba.normalize(dba.obj.count, normalize = normalization)

  if (min_members == 2){
    dba.obj.contrast <- DiffBind::dba.contrast(dba.obj.normalize,
                                              minMembers = min_members,
                                              reorderMeta = list(Condition = control_level))
  } else if (is.null(min_members)) {
    dba.obj.contrast <- DiffBind::dba.contrast(dba.obj.normalize,
                                              reorderMeta = list(Condition = control_level))
  }

  dba.obj <- DiffBind::dba.analyze(dba.obj.contrast)
  return(dba.obj)
}

#' Save DBA object
#' It returns nothing
#' save_diffbind_object(dba_object, file_suffix, dbaobj_save_path) saves the
#' DBA object under the specified file path.
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param file_suffix Character. Suffix for the output DBA object file.
#' @param dbaobj_save_path Character. Directory path to save the resulting
#' DBA object (default is \code{NULL}). If \code{NULL}, current working
#' directory is used.
#' @export
#' @examples
#' save_diffbind_object(
#'   dba_object = dba.obj,
#'   file_suffix = "h3k27me3",
#'   dbaobj_save_path = "test_results")

save_diffbind_object <- function(
    dba_object,
    file_suffix = NULL,
    dbaobj_save_path = NULL
){
  if (is.null(file_suffix)){
    stop("'file_suffix' is required to uniquely identify the saved RDS file.")
  }
  if (is.null(dbaobj_save_path)) {
    dbaobj_save_path <- getwd()
  }

  #Save the DBA object
  saveRDS(dba_object, file=file.path(dbaobj_save_path, paste0("dba_obj_", file_suffix, ".RDS")))

  message("DiffBind analysis completed successfully. Results saved to:",
          file.path(dbaobj_save_path, paste0("dba_obj_", file_suffix, ".RDS")))
}

#' Retrieve Differentially Bind Sites from a DBA object
#'
#' get_DBsites(dba_object, fdr_threshold, contrast_number) Extracts differentially
#' bound (DB) regions from a DBA object and filters regions based on a user-defined
#' FDR threshold, finally returns a dataframe that contains sites information.
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param fdr_threshold Numeric. Significance cutoff for filtering DB sites.
#' Default is \code{0.05}.
#' @param contrast_number Integer. Index specifying which contrast to extract
#' results from. Defaults to \code{1} corresponding to the first contrast.
#' @return A dataframe that contains filtered sites information.
#' @export
#' @examples
#' get_DBsites(dba_object = dba.obj)

get_DBsites <- function(
    dba_object,
    fdr_threshold = 0.05,
    contrast_number = 1
){
  if (contrast_number <= 0 || contrast_number > length(dba.obj$contrasts)) {
    stop("Invalid contrast_number. Please provide a valid contrast index.")
  }

  sites <- as.data.frame(dba.report(dba_object, contrast = contrast_number, th = 1))

  if (is.null(sites)) {
    stop("No sites found in the report. Please check your contrast setup.")
  }

  #filter sites
  filtered.sites <- sites %>% filter(FDR <= fdr_threshold)

  print(paste('Positive logFC:', sum(filtered.sites$Fold > 0)))
  print(paste('Negative logFC:', sum(filtered.sites$Fold < 0)))

  return(filtered.sites)
}

#' save different sites files
#' It doesn't return anything
#' save_sites(dba_object, save_path, fdr_threshold, contrast_number) saves site
#' dataframes in bed format based on different fdr_threshold, which includes allDB
#' sites, upDB sites, and downDB sites.
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param save_path Character. Output directory where results are written.
#' @param fdr_threshold Numeric. If it is 1, then it saves all DB sites, if
#' it is 0.05, then it save upDB and downDB sites separately.
#' @param contrast_number Integer. Index specifying which contrast to extract
#' results from. Defaults to \code{1} corresponding to the first contrast.
#' @export
#' @examples
#' save_sites(
#'   dba_object = dba.obj,
#'   save_path = "diffbind_files")

save_sites <- function(
    dba_object,
    save_path,
    fdr_threshold = 0.05,
    contrast_number = 1
){
  if (is.null(save_path)) {
    save_path <- getwd()
  }
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  sig.sites <- get_DBsites(dba_object,
                           fdr_threshold = fdr_threshold,
                           contrast_number = contrast_number)

  #UpDB.sites sorted by FDR
  upDB.sites <- sig.sites %>%
    filter(Fold > 0) %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    arrange(FDR) %>%
    select(Chr, Start, End, site.id, Score) %>%
    mutate(Position = paste0(Chr, ':', Start, "-", End))

  write.table(upDB.sites, file = file.path(save_path, 'upDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)

  #DownDB.sites sorted by FDR
  DownDB.sites <- sig.sites %>%
    filter(Fold < 0) %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    arrange(FDR) %>%
    select(Chr, Start, End, site.id, Score) %>%
    mutate(Position = paste0(Chr, ':', Start, "-", End))

  write.table(DownDB.sites, file = file.path(save_path, 'downDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)

  all.sites <- get_DBsites(dba_object,
                          fdr_threshold = 1,
                          contrast_number = contrast_number)

  #all.sites order by chromosomes
  Consensus.sites <- all.sites %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    select(Chr, Start, End, site.id, FDR, Score) %>%
    mutate(Position = paste0(Chr, ':', Start, "-", End)) %>%
    slice(gtools::mixedorder(Chr), Start) %>%
    arrange(Chr, Start)

  write.table(Consensus.sites, file = file.path(save_path,'consensus_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)
}

#' Create and save a density plot
#'
#' make_density_plot(dba_object, figure_title, file_name, fdr_threshold, figure_save_path,
#' contrast_number, width, height, color) creates a density plot from a DBA object
#' with specified figure title and saves it in pdf and png format under the designated
#' file path.
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param figure_title Character. Title for density plot.
#' @param file_name Character. Base name used for saving the plot files.
#' @param fdr_threshold Numeric. FDR cutoff for significance (Default is 1).
#' @param figure_save_path Character. Directory path where the plot files will be saved.
#' Default is set to \code{NULL}.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object. Defeault is set to \code{1}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param color Character. Color to fill the density plot. Default is \code{"aquamarine4"}.
#' @return A ggplot object of density plot.
#' @keywords Density Plot
#' @export
#' @examples
#' make_density_plot(
#'   dba_object = dba.obj,
#'   figure_title = "Distribution of Differential Binding Regions",
#'   file_name = "df_fold_density_h3k27ac",
#'   figure_save_path = "diffbind_figures")

make_density_plot <- function(
    dba_object,
    figure_title,
    file_name,
    fdr_threshold = 1,
    figure_save_path=NULL,
    contrast_number=1,
    width=8,
    height=6,
    color='aquamarine4'
){
  if (is.null(figure_save_path)) {
    figure_save_path <- getwd()
  }
  if (is.null(file_name)) {
    stop("A valid file name is required to be specified.")
  }

  sites <- get_DBsites(dba_object = dba_object,
                      fdr_threshold = fdr_threshold,
                      contrast_number = contrast_number)

  #draw density plot
  p <- sites %>%
    ggplot(aes(x=Fold)) +
    geom_density(aes(y=after_stat(count)), fill=color) +
    labs(x = TeX('$\\log_2$ FC'),
         title = TeX(paste0(figure_title,
                            ' [ FDR < ', fdr_threshold, ' ]')))

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(figure_save_path, file_name),
    width = width, height = height, dpi=600,
    plot = p
  )

  message("Density Plot is generated successfully. Results saved to: \n",
          file.path(figure_save_path, paste0(file_name ,'.png')), " and \n",
          file.path(figure_save_path, paste0(file_name ,'.pdf')))
  return(p)
}

#' Create and save a PCA plot
#'
#' make_PCA_plot_diffbind(dba_object, figure_title_nocontrast, figure_title_contrast,
#' figure_save_path, file_name, size, width, height, color) Creates side-by-side PCA
#' plots to compare the variance structure of all consensus binding regions vs differentially
#' bound (DB) regions from a DBA object. Both plots are saved as PNG and PDF files under
#' designated file path.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param figure_title_nocontrast Character. Title for the PCA plot using all consensus
#' regions.
#' @param figure_title_contrast Character. Title for the PCA plot using only differentially
#' bound regions.
#' @param figure_save_path Character. Directory path where output plots will be saved.
#' @param file_name Character. Base name (without extension) for the saved plot files.
#' @param size Numeric. Size of the points, default is set to 8.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 20.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 9.
#' @param color Character vector of length 2. Colors to use for the two biological conditions
#' in PCA plots. Default is \code{c("darkmagenta", "aquamarine4")}.
#' @return A \code{gridExtra} object containing a side-by-side layout of two PCA plots.
#' @keywords PCA plot
#' @import gridExtra
#' @export
#' @examples
#' make_PCA_plot_diffbind(
#'   dba_object = dba.obj,
#'   figure_title_nocontrast = "PCA of All Consensus Regions",
#'   figure_title_contrast = "PCA of Differential Binding Regions",
#'   figure_save_path = "diffbind_figures",
#'   file_name = "diffbind_pca_plots")

make_PCA_plot_diffbind <- function(
    dba_object,
    figure_title_nocontrast,
    figure_title_contrast,
    figure_save_path,
    file_name,
    size=8,
    width=20,
    height=9,
    color=c('darkmagenta','aquamarine4')
){
  if (is.null(figure_save_path)) {
    figure_save_path = getwd()
  }

  # For all consensus regions:
  df_nocontrast_pre <- data.frame(dba.peakset(dba_object, bRetrieve = TRUE))
  df_nocontrast <- t(df_nocontrast_pre[,6:ncol(df_nocontrast_pre)])
  df_label <- dba.obj$samples$Condition
  PoV_nocontrast_rank <- prcomp(df_nocontrast, center = TRUE)$sdev^2
  sum_pc_nocontrast <- sum(PoV_nocontrast_rank[1:2])
  legend_label <- paste0("Sum of variance of 2 PCs: ", sum_pc_nocontrast)
  p1 <- kwanlibr::draw_PCA(df_nocontrast, label=df_label, color=color) +
    geom_point(size = size) +
    labs(caption = legend_label) +
    ggtitle(figure_title_nocontrast) +
    theme(aspect.ratio = 1) +
    theme_bw()

  # For all DB regions:
  df_contrast_pre <- dba.report(dba_object, bCounts = TRUE)
  df_contrast <- t(df_contrast_pre[,10:ncol(df_contrast_pre)])
  rownames(df_contrast) <- NULL
  PoV_contrast_rank <- prcomp(df_contrast, center = TRUE)$sdev^2
  sum_pc_contrast <- sum(PoV_contrast_rank[1:2])
  legend_label <- paste0("Sum of variance of 2 PCs: ", sum_pc_contrast)
  p2 <- kwanlibr::draw_PCA(df_contrast, label = df_label, color=color) +
    geom_point(size = size) +
    labs(caption = legend_label) +
    ggtitle(figure_title_contrast) +
    theme(aspect.ratio = 1) +
    theme_bw()

  #combine two plots into one
  p <- gridExtra::grid.arrange(p1, p2, nrow=1)

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(figure_save_path, file_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  message("PCA Plot is generated successfully. Results saved to: \n",
          file.path(figure_save_path, paste0(file_name ,'.png')), " and \n",
          file.path(figure_save_path, paste0(file_name ,'.pdf')))
  return(p)
}

#' Generate volcano-plot-ready data from DiffBind results
#'
#' make_volcano_sites(dba_object, contrast_number, fdr_threshold, xdiff, ymax) makes
#' volcano-plottable table and applies filtering based on FDR threshold and clamps
#' values.
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object. Defeault is set to \code{1}.
#' @param fdr_threshold Numeric. FDR cutoff for significance. Default is set to \code{0.05}.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @return A dataframe that contains volcano-plottable data.
#' @export
#' @examples
#' make_volcano_sites(dba_object = dba.obj)

make_volcano_sites <- function(
    dba_object,
    contrast_number=1,
    fdr_threshold=0.05,
    xdiff=5,
    ymax=40
){
  sites <- get_DBsites(dba_object,
                      fdr_threshold = 1,
                      contrast_number = contrast_number)

  sig.sites <- get_DBsites(dba_object,
                          fdr_threshold = fdr_threshold,
                          contrast_number = contrast_number)

  downsampled_nonsig_sites <- sites %>%
    filter(FDR > fdr_threshold) %>%
    sample_frac(1, replace=FALSE)

  volcano.sites <- sig.sites %>%
    bind_rows(downsampled_nonsig_sites) %>%
    mutate(negLogFDR = -log10(FDR)) %>%
    mutate(Fold = kwanlibr::clamp(Fold, -xdiff, xdiff),
           negLogFDR = kwanlibr::clamp(negLogFDR, 0, ymax)) %>%
    mutate(FDR = if_else(FDR <= fdr_threshold,
                         paste('<=',fdr_threshold),
                         paste('>',fdr_threshold)))
  return(volcano.sites)
}

#' Create and save a volcano plot
#'
#' make_volcano_plot_diffbind(dba_object, figure_title, figure_save_path, file_name,
#' size, alpha, width, height, xdiff, ymax, contrast_number, fdr_threshold, color)
#' generates volcano plot and returns the ggplot object and it also saves the file
#' under the designated figure file path.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param figure_title Character. Figure title for volcano plot, ex. "h3k27me3 cKO vs cHET"
#' @param figure_save_path Character. Directory path where the output plots will
#' be saved. Default is \code{NULL}.
#' @param file_name Character. Base name for the output volcano plot figures.
#' @param size Numeric. Size of the points in volcano plots, default is set to 1.
#' @param alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque). Default is \code{1}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @param contrast_number Integer. Index indicating which contrast to extract from
#' the DBA object. Defeault is set to \code{1}.
#' @param fdr_threshold Numeric. The FDR threshold used for grouping, default is \code{0.05}.
#' @param color Character vector of length 2. By default the color is \code{c('aquamarine4', 'grey')}
#' @return A ggplot object of volcano plot
#' @keywords volcano plot
#' @export
#' @examples
#' make_volcano_plot_diffbind(
#'   dba_object = dba.obj,
#'   figure_title = "h3k27ac cKO vs cHET",
#'   figure_save_path = "diffbind_figures",
#'   file_name = "db_volcano")

make_volcano_plot_diffbind <- function(
    dba_object,
    figure_title,
    figure_save_path=NULL,
    file_name,
    size=1,
    alpha=1,
    width=8,
    height=6,
    xdiff=5,
    ymax=40,
    contrast_number=1,
    fdr_threshold=0.05,
    color=c('aquamarine4', 'grey')
){
  if (is.null(figure_save_path)) {
    figure_save_path = getwd()
  }

  volcano.sites <- make_volcano_sites(dba.obj,
                                      xdiff = xdiff,
                                      ymax = ymax)

  sig.sites <- get_DBsites(dba_object = dba_object,
                           fdr_threshold = fdr_threshold,
                           contrast_number = contrast_number)

  p <-draw_volcano_general(volcano.sites,
                           criteria="FDR",
                           colors = setNames(c(color[1], color[2]),
                                             c(paste('<=', fdr_threshold),
                                               paste('>', fdr_threshold))),
                           size=size,
                           alpha=alpha,
                           xdiff = xdiff,
                           ymax = ymax)

  figure_title <- paste0(figure_title, " [", dim(sig.sites)[1], ' DB Regions FDR <= ',
                  fdr_threshold,
                  ']')

  p <- p + labs(title = figure_title)

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(figure_save_path, file_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  message("Volcano Plot is generated successfully. Results saved to: \n",
          file.path(figure_save_path, paste0(file_name ,'.png')), " and \n",
          file.path(figure_save_path, paste0(file_name ,'.pdf')))
  return(p)
}

#' Merge Differential Binding and Expression Data
#'
#' merge_sites_with_exp(dba_object, RNAfile_path, join_data_sites, save_path, contrast_number,
#' fdr_threshold, bulk_fdr_cutoff) integrates differential binding sites with differential
#' expression data (from bulk RNA-seq) by joining them based on their nearest gene TSS and
#' returns a data frame where each row represents a differential binding site along with its
#' nearest gene's differential expression.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param RNAfile_path Character. Path to the CSV file containing bulk RNA-seq differential
#' expression results. This file must include the columns \code{gene_name}, \code{logFC},
#' and \code{FDR}.
#' @param join_data_sites dataframe (e.g., volcano_sites) to be merged with DE data.
#' @param save_path Character. Directory path for storing intermediate files, including the TSS
#' BED file, sorted BED file, and nearest gene lookup TSV.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object. Defeault is set to \code{1}.
#' @param fdr_threshold Numeric. Significance threshold for FDR in differential binding dataframe.
#' Default is 0.05.
#' @param bulk_fdr_cutoff Numeric. Significance threshold for RNA-seq differential expression.
#' Defualt is 0.05.
#' @return A data frame containing merged binding and expression data. Each row corresponds to
#' a binding site and includes:
#' \itemize{
#'   \item Site genomic coordinates
#'   \item Site ID
#'   \item Nearest gene name
#'   \item \code{bulk.logFC} and \code{bulk.FDR} values
#' }
#' @export
#' @examples
#' merge_sites_with_exp(
#'   dba_object = dba.obj,
#'   RNAfile_path = "path_to_RNAseq_csv_file",
#'   join_data_sites = volcano.sites,
#'   save_path = "diffbind_files")

merge_sites_with_exp <- function(
    dba_object,
    RNAfile_path,
    join_data_sites,
    save_path,
    contrast_number=1,
    fdr_threshold=0.05,
    bulk_fdr_cutoff=0.05
){
  if (is.null(RNAfile_path) || !file.exists(RNAfile_path)) {
    stop("`RNAfile_path` must be a valid path to the RNA-seq results CSV file.")
  }
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  sig.sites <- get_DBsites(dba_object = dba_object,
                           fdr_threshold = fdr_threshold,
                           contrast_number = contrast_number)

  allDB.sites = sig.sites %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    select(Chr, Start, End, site.id, Score)

  write.table(allDB.sites, file = file.path(save_path, 'allDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)

  # run bash script to retrieve the nearest gene ID
  allDBsite_path <- file.path(save_path, "allDB_sites.bed")
  tss_bed <- file.path(save_path, "TSS.bed")
  sorted_bed <- file.path(save_path, "allDB_sites.sorted.bed")
  lookup_tsv <- file.path(save_path, "DB_site_nearest_gene_lookup.tsv")

  script <- system.file("scripts", "db_annotate_regions.sh", package = "kwanlibr")
  if (script == "") stop("Annotation script not found in kwanlibr package.")
  system2("bash", args = c(script, allDBsite_path, tss_bed, sorted_bed, lookup_tsv))

  # read the annotated gene expression RNA sequencing file
  db.site.gene.lookup = read.table(lookup_tsv, col.names = c('site.id','gene_name'))
  bulk.DE = read.csv(RNAfile_path) %>% distinct(gene_name, .keep_all = TRUE)

  # inner join sites and the diff expressed data
  db.site.bulk.data = db.site.gene.lookup %>%
    inner_join(bulk.DE, by = 'gene_name') %>%
    select(site.id, gene_name, logFC, FDR) %>%
    dplyr::rename(bulk.logFC = logFC) %>%
    dplyr::rename(bulk.FDR = FDR)

  db.site.bulk.data.meets.FDR = db.site.bulk.data %>% filter(bulk.FDR < bulk_fdr_cutoff)

  data.sites.bulk.color = join_data_sites %>%
    tibble::rownames_to_column('site.id') %>%
    mutate(site.id = as.numeric(site.id)) %>%
    left_join(db.site.bulk.data.meets.FDR, by = 'site.id')

  return(data.sites.bulk.color)
}

#' Make Volcano Plot from DB and DE Data (Color by bulk logFC direction)
#'
#' make_volcano_plot_from_merged(merged_df, figure_title, figure_save_path, file_name,
#' color, xdiff, ymax, size, alpha, width, height, contrast_number, bulk_fdr_cutoff)
#' merges differential binding (DB) sites with RNA-seq differential expression (DE)
#' data using \code{merge_sites_with_exp()}, then creates a volcano plot. All sites
#' are shown in gray; DB sites mapped to significant DE genes are colored by direction
#' of bulk logFC (Up/Down).
#'
#' @param merged_df Merged data generated from running \code{merge_sites_with_exp()}
#' @param figure_title Character. The title of the volcano plot.
#' @param figure_save_path Character. The path to save the volcano plot.
#' @param file_name Character. Output figure file_name (e.g., "volcano.png").
#' @param color Character vector of length 3, in order of background color, down FC color,
#' and up FC color. E.x. \code{c('grey','steelblue', 'tomato')}.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @param size Numeric. Point size in the plot. Default is \code{1}.
#' @param alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque). Default is \code{1}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param contrast_number Integer. Index indicating which contrast to extract from
#' the DBA object. Defeault is set to \code{1}.
#' @param bulk_fdr_cutoff Numeric. Significance threshold for RNA-seq differential expression.
#' Defualt is 0.05.
#' @return A ggplot object of DB and DE merged volcano plot
#' @export
#' @examples
#' make_volcano_plot_from_merged(
#'   merged_df = merged_volcano,
#'   figure_title = 'cKO vs cHet Differential Binding sites',
#'   figure_save_path = "diffbind_figures",
#'   file_name = "db_volcano_bulk_color_binary",
#'   color = c('grey','steelblue', 'tomato'))

make_volcano_plot_from_merged <- function(
    merged_df,
    figure_title,
    figure_save_path,
    file_name,
    color,
    xdiff=5,
    ymax=40,
    size=1,
    alpha=1,
    width = 8,
    height = 6,
    contrast_number=1,
    bulk_fdr_cutoff=0.05
){
  # Add criteria column
  merged_df <- merged_df %>%
    mutate(
      bulk.logFC.direction = case_when(
        !is.na(bulk.FDR) & bulk.FDR < bulk_fdr_cutoff & bulk.logFC > 0 ~ "Up",
        !is.na(bulk.FDR) & bulk.FDR < bulk_fdr_cutoff & bulk.logFC < 0 ~ "Down",
        TRUE ~ "NM"),  # Not Matched
      # Set explicit factor levels to control draw order (NS first)
      bulk.logFC.direction = factor(bulk.logFC.direction, levels = c("NM", "Down", "Up"))
    )

  # Plot
  p <- draw_volcano_general(
    data = merged_df,
    colors = setNames(c(color[1], color[2], color[3]), c("NM", "Down", "Up")),
    criteria = "bulk.logFC.direction",
    size = size,
    alpha = alpha,
    xdiff = xdiff,
    ymax = ymax
  )

  # Add titles etc
  p <- p +
    # exclude "NM" in legend, keep only "down" and "up"
    scale_color_manual(name = "bulkRNAseq FC",
                       values = setNames(c(color[1], color[2], color[3]), c("NM", "Down", "Up")),
                       breaks = c("Down", "Up")) +
    xlim(-xdiff, xdiff) +
    ylim(0, ymax) +
    labs(title = paste0(figure_title,
                        '\n[ color by nearest TSS bulkRNAseq logFC, filtered by bulkRNAseq FDR < ',
                        bulk_fdr_cutoff, ' ]'),
         color = "bulk RNA-seq FC",
         x = TeX('$\\log_2$( differential binding FC )'),
         y = TeX('$-\\log_{10}$( differential binding FDR )'))

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(figure_save_path, file_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  return(p)
}

#' Create and save the scattor plot from DB and DE Data
#'
#' make_scatter_plot_from_merged(merged_df, figure_title, figure_save_path, file_name,
#' width, height, size, regression, contrast_number, point_color, line_color) create
#' scatter plot from merged differential binding sites with differential expression
#' data and save the plot under designated directory.
#'
#' @param merged_df Merged data generated from running \code{merge_sites_with_exp()}
#' @param figure_title Character. The scatter plot title.
#' @param figure_save_path Character. Directory to save the figure output.
#' @param file_name Character. The base name of the scatter plot file.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param size Numeric. The size of the point in scatter plot. Default is 1.
#' @param regression logical value. Default is set to 'TRUE', a regression line is
#' drawn to the scatter plot. otherwise no regression line.
#' @param contrast_number Integer. Index indicating which contrast to extract from
#' the DBA object. Defeault is set to \code{1}.
#' @param point_color Character Value. The color of the point. Default is \code{'black'}.
#' @param line_color Character Value. The color of the regression line. Default is
#' \code{'aquamarine4'}.
#' @export
#' @examples
#' make_scatter_plot_from_merged(
#'   merged_df = merged_sigsites,
#'   figure_title = "DB vs DE log2FC",
#'   figure_save_path = "diffbind_figures",
#'   file_name = "db_bulk_scatter")

make_scatter_plot_from_merged <- function(
    merged_df,
    figure_title,
    figure_save_path,
    file_name,
    width=8,
    height=6,
    size=1,
    regression=TRUE,
    contrast_number=1,
    point_color='black',
    line_color='aquamarine4'
){
  # exclude rows with NA Values in any column
  merged_df <- merged_df %>% na.omit()

  p <- merged_df %>%
    ggplot(aes(x=Fold, y=bulk.logFC)) +
    geom_point(color=point_color,
               size=size)

  if (regression) {
    db.bulk.scatter.regression = glm(bulk.logFC ~ Fold, data = merged_df)
    beta = summary(db.bulk.scatter.regression)$coefficients['Fold', 1]
    p.val = summary(db.bulk.scatter.regression)$coefficients['Fold', 4]
    p <- p + geom_smooth(method = glm, color = line_color)+
      xlim(c(-1,1) * max(abs(merged_df$Fold))) +
      ylim(c(-1,1) * max(abs(merged_df$bulk.logFC))) +
      labs(x = TeX('Differential Binding Region $\\log_2$FC'),
           y = TeX('Nearby Gene RNAseq $\\log_2$FC'),
           title = TeX(paste0(figure_title, " ",
                              '$\\beta_{DB.region} =', round(beta, 2), '$,',
                              '  p value = ', signif(p.val, digits=2))))
  } else {
    p <- p +
      xlim(c(-1,1) * max(abs(merged_df$Fold))) +
      ylim(c(-1,1) * max(abs(merged_df$bulk.logFC))) +
      labs(x = TeX('Differential Binding Region $\\log_2$FC'),
           y = TeX('Nearby Gene RNAseq $\\log_2$FC'),
           title = TeX(paste0(figure_title)))
  }

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(figure_save_path, file_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  return(p)
}
