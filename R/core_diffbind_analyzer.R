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

  if (is.null(min_members)) {
    dba.obj.contrast <- DiffBind::dba.contrast(dba.obj.normalize,
                                               reorderMeta = list(Condition = control_level))
  } else {
    dba.obj.contrast <- DiffBind::dba.contrast(dba.obj.normalize,
                                               minMembers = min_members,
                                               reorderMeta = list(Condition = control_level))
  }

  dba.obj <- DiffBind::dba.analyze(dba.obj.contrast)

  message("DiffBind analysis completed successfully.")

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
    file_suffix,
    dbaobj_save_path = NULL
){
  if (is.null(dbaobj_save_path)) {
    dbaobj_save_path <- getwd()
  }

  #Save the DBA object
  saveRDS(dba_object, file=file.path(dbaobj_save_path, paste0("dba_obj_", file_suffix, ".RDS")))
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

  sites <- as.data.frame(dba.report(dba_object, contrast = contrast_number, th = fdr_threshold))

  if (is.null(sites)) {
    stop("No sites found in the report. Please check your contrast setup.")
  }

  print(paste('Positive logFC:', sum(sites$Fold > 0)))
  print(paste('Negative logFC:', sum(sites$Fold < 0)))

  return(sites)
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
#' contrast_number, xdiff, width, height, color) creates a density plot from a DBA object
#' with specified figure title and saves it in pdf and png format under the designated
#' file path.
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param figure_title Character. Title for density plot.
#' @param file_name Character. Base name used for saving the plot files.
#' @param fdr_threshold Numeric. FDR cutoff for DB significance (Default is 1).
#' @param figure_save_path Character. Directory path where the plot files will be saved.
#' Default is set to \code{NULL}.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object. Defeault is set to \code{1}.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is NULL.
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
    xdiff=NULL,
    width=8,
    height=6,
    color='aquamarine4'
){
  if (is.null(file_name)) {
    stop("A valid file name is required to be specified.")
  }

  #strip off extensions in the filename
  if (tools::file_ext(file_name) != "") {
    file_name <- tools::file_path_sans_ext(file_name)
    warning("file name should not include any extensions, file name changed to: ", file_name)
  }

  if (is.null(figure_save_path)) {
    figure_save_path <- getwd()
  }

  sites <- get_DBsites(dba_object = dba_object,
                      fdr_threshold = fdr_threshold,
                      contrast_number = contrast_number)

  #draw density plot
  p <- sites %>%
    ggplot(aes(x=Fold)) +
    geom_density(aes(y=after_stat(count)), fill=color) +
    labs(x = TeX('$\\log_2$ FC'), title = figure_title)

  # Add xlim conditionally
  if (!is.null(xdiff)) {
    p <- p + xlim(-xdiff, xdiff)
  }

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
#' figure_save_path, file_name, point_size, width, height, color) Creates side-by-side PCA
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
#' @param point_size Numeric. Size of the points, default is set to 8.
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
    point_size=8,
    width=20,
    height=9,
    color=c('darkmagenta','aquamarine4')
){
  #strip off extensions in the filename
  if (tools::file_ext(file_name) != "") {
    file_name <- tools::file_path_sans_ext(file_name)
    warning("file name should not include any extensions, file name changed to: ", file_name)
  }

  if (is.null(figure_save_path)) {
    figure_save_path = getwd()
  }

  if (length(color) != 2) {
    stop("There should be 2 color choices.")
  }

  # For all consensus regions:
  df_nocontrast_pre <- data.frame(dba.peakset(dba_object, bRetrieve = TRUE))
  df_nocontrast <- t(df_nocontrast_pre[,6:ncol(df_nocontrast_pre)])
  df_label <- dba.obj$samples$Condition
  PoV_nocontrast_rank <- prcomp(df_nocontrast, center = TRUE)$sdev^2
  sum_pc_nocontrast <- sum(PoV_nocontrast_rank[1:2])
  legend_label <- paste0("Sum of variance of 2 PCs: ", sum_pc_nocontrast)
  p1 <- kwanlibr::draw_PCA(df_nocontrast, label=df_label, color=color) +
    geom_point(size = point_size) +
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
    geom_point(size = point_size) +
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
#'
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
  print("For all sites:")
  sites <- get_DBsites(dba_object,
                      fdr_threshold = 1,
                      contrast_number = contrast_number)
  print("For significant sites:")
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
#' make_volcano_plot_diffbind(dba_object, figure_title, file_name, figure_save_path,
#' point_size, point_alpha, width, height, xdiff, ymax, contrast_number, fdr_threshold, color)
#' generates volcano plot and returns the ggplot object and it also saves the file
#' under the designated figure file path.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param figure_title Character. Figure title for volcano plot, ex. "h3k27me3 cKO vs cHET"
#' @param file_name Character. Base name for the output volcano plot figures.
#' @param figure_save_path Character. Directory path where the output plots will
#' be saved. Default is \code{NULL}.
#' @param point_size Numeric. Size of the points in volcano plots, default is set to 1.
#' @param point_alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque). Default is \code{1}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @param contrast_number Integer. Index indicating which contrast to extract from
#' the DBA object. Defeault is set to \code{1}.
#' @param fdr_threshold Numeric. The FDR threshold used for grouping, default is \code{0.05}.
#' @param color Character vector of length 2. One for insignificant DB sites, one for significant
#' DB sites. Usually 'grey' is for insignificant sites. By default, it is \code{c('aquamarine4', 'grey')}
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
    file_name,
    figure_save_path=NULL,
    point_size=1,
    point_alpha=1,
    width=8,
    height=6,
    xdiff=5,
    ymax=40,
    contrast_number=1,
    fdr_threshold=0.05,
    color=c('aquamarine4', 'grey')
){
  #strip off extensions in the filename
  if (tools::file_ext(file_name) != "") {
    file_name <- tools::file_path_sans_ext(file_name)
    warning("file name should not include any extensions, file name changed to: ", file_name)
  }

  if (is.null(figure_save_path)) {
    figure_save_path = getwd()
  }

  if (length(color) != 2) {
    stop("There should be 2 color choices.")
  }

  volcano.sites <- make_volcano_sites(dba.obj,
                                      contrast_number = contrast_number,
                                      fdr_threshold = fdr_threshold,
                                      xdiff = xdiff,
                                      ymax = ymax)

  p <- ggplot(data = volcano.sites, aes(x = Fold, y = negLogFDR, color = FDR)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = color) +
    xlim(-xdiff, xdiff) +
    ylim(0, ymax) +
    labs(x = TeX('$\\log_2$ FC'),
         y = TeX('$-\\log_{10}$ FDR'),
         title = figure_title)

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
#' merge_sites_with_DE(dba_object, DE_file_path, tss_file_path, save_path, contrast_number)
#' integrates differential binding sites with differential expression data (from bulk RNA-seq)
#' by joining them based on their nearest gene TSS and returns a data frame where each row
#' represents a differential binding site along with its nearest gene's differential expression.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param DE_file_path Character. Path to the CSV file containing bulk RNA-seq differential
#' expression results. This file must include the columns \code{gene_name}, \code{logFC},
#' and \code{FDR}.
#' @param tss_file_path Character. Transcription Start Site file path to be used in annotation
#' with Differential binding sites.
#' @param save_path Character. Directory path for storing intermediate files, including the TSS
#' BED file, sorted BED file, and nearest gene lookup TSV.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object. Defeault is set to \code{1}.
#' @return A data frame containing merged binding and expression data. Each row corresponds to
#' a binding site and includes:
#' \itemize{
#'   \item Site genomic coordinates
#'   \item Site ID
#'   \item Nearest gene name
#'   \item \code{DE.logFC} and \code{DE.FDR} values
#'   \item \code{DB.logFC} and \code{DB.FDR} values
#' }
#' @export
#' @examples
#' merge_sites_with_DE(
#'   dba_object = dba.obj,
#'   DE_file_path = path_to_DE_csv_file,
#'   tss_file_path = path_to_tss_file,
#'   save_path = "diffbind_files")

merge_sites_with_DE <- function(
    dba_object,
    DE_file_path,
    tss_file_path,
    save_path,
    contrast_number=1
){
  if (is.null(DE_file_path) || !file.exists(DE_file_path)) {
    stop("`DE_file_path` must be a valid path to the RNA-seq results CSV file.")
  }

  if (!file.exists(tss_file_path)) {
    stop("`tss_file_path` must be a valid path to the TSS BED file.")
  }

  #check if DE file has required column
  DE = read.csv(DE_file_path)
  required_cols <- c("gene_name", "logFC", "FDR")
  missing_cols <- setdiff(required_cols, colnames(DE))
  if (length(missing_cols) > 0) {
    stop("The DE file is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  sites <- get_DBsites(dba_object = dba_object,
                       fdr_threshold = 1,
                       contrast_number = contrast_number)

  allDB.sites = sites %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    select(Chr, Start, End, site.id, Score)

  write.table(allDB.sites, file = file.path(save_path, 'allDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)

  # run bash commands to retrieve the nearest gene ID
  allDB_bed <- file.path(save_path, "allDB_sites.bed")
  lookup_tsv <- file.path(save_path, "DB_site_nearest_gene_lookup.tsv")
  sorted_bed <- tempfile(pattern = "sorted_allDB_sites_", fileext = ".bed")

  cmd <- sprintf(
    'module load Bioinformatics bedops/2.4.41 && sort-bed %s > %s && closest-features --closest --delim "\t" %s %s | awk \'{ print $4, $9 }\' > %s',
    shQuote(allDB_bed),
    shQuote(sorted_bed),
    shQuote(sorted_bed),
    shQuote(tss_file_path),
    shQuote(lookup_tsv)
  )
  system2("bash", args = c("-l", "-c", cmd))

  # read the annotated gene expression RNA sequencing file
  db_site_gene_lookup = read.table(lookup_tsv,
                                   fill = TRUE,
                                   col.names = c('site.id','gene_name'))

  DE = DE %>% distinct(gene_name, .keep_all = TRUE)

  # inner join sites and the diff expressed data
  db_site_DE_data = db_site_gene_lookup %>%
    inner_join(DE, by = 'gene_name') %>%
    select(site.id, gene_name, logFC, FDR) %>%
    dplyr::rename(DE.logFC = logFC) %>%
    dplyr::rename(DE.FDR = FDR)

  full_sites_and_DE_data = sites %>%
    tibble::rownames_to_column('site.id') %>%
    mutate(site.id = as.numeric(site.id)) %>%
    dplyr::rename(DB.logFC = Fold) %>%
    dplyr::rename(DB.FDR = FDR) %>%
    left_join(db_site_DE_data, by = 'site.id')

  return(full_sites_and_DE_data)
}

#' Make Volcano Plot from DB and DE Data (Color by DE logFC direction)
#'
#' make_volcano_plot_from_merged(merged_df, figure_title, figure_save_path, file_name,
#' xdiff, ymax, point_size, point_alpha, width, height, DB_fdr_cutoff, DE_fdr_cutoff,
#' color) creates and saves a volcano plot based on merged DB and DE data. All sites
#' are shown in gray; DB sites mapped to significant DE genes are colored by direction
#' of DE logFC (Up/Down).
#'
#' @param merged_df Merged data generated from running \code{merge_sites_with_DE()}
#' @param figure_title Character. The title of the volcano plot.
#' @param figure_save_path Character. The path to save the volcano plot.
#' @param file_name Character. Output figure file_name (e.g., "volcano.png").
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @param point_size Numeric. Point size in the plot. Default is \code{1}.
#' @param point_alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque). Default is \code{1}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param DB_fdr_cutoff Numeric. Significance threshold for FDR in differential binding
#' dataframe. Default is 0.05.
#' @param DE_fdr_cutoff Numeric. Significance threshold for FDR in RNA-seq differential
#' expression. Defualt is 0.05.
#' @param color Character vector of length 2, in order of down FC color, then up FC
#' color. default is set to \code{c('steelblue', 'tomato')}.
#' @return A ggplot object of DB and DE merged volcano plot
#' @export
#' @examples
#' make_volcano_plot_from_merged(
#'   merged_df = merged_volcano,
#'   figure_title = 'cKO vs cHet Differential Binding sites',
#'   figure_save_path = "diffbind_figures",
#'   file_name = "db_volcano_DE_color_binary")

make_volcano_plot_from_merged <- function(
    merged_df,
    figure_title,
    figure_save_path,
    file_name,
    xdiff=5,
    ymax=40,
    point_size=1,
    point_alpha=1,
    width = 8,
    height = 6,
    DB_fdr_cutoff=0.05,
    DE_fdr_cutoff=0.05,
    color=c("steelblue", "tomato")
){
  # Check for required columns
  required_cols <- c("DB.logFC", "DB.FDR", "DE.logFC", "DE.FDR")
  missing_cols <- setdiff(required_cols, colnames(merged_df))
  if (length(missing_cols) > 0) {
    stop("The merged data is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  #strip off extensions in the filename
  if (tools::file_ext(file_name) != "") {
    file_name <- tools::file_path_sans_ext(file_name)
    warning("file name should not include any extensions, file name changed to: ", file_name)
  }

  if (length(color) != 2) {
    stop("There should be 2 color choices.")
  }

  #filtering by FDR
  sig_merged_df <- merged_df %>% filter(DB.FDR < DB_fdr_cutoff)

  downsampled_nonsig_merged_df <- merged_df %>%
    filter(DB.FDR > DB_fdr_cutoff) %>%
    sample_frac(1, replace = FALSE)

  volcano_merged_df <- sig_merged_df %>% bind_rows(downsampled_nonsig_merged_df)
  volcano_merged_df <- merged_df %>% mutate(negLogDBFDR = -log10(DB.FDR))

  # Add criteria column
  volcano_DE_groups <- volcano_merged_df %>%
    filter(DB.FDR < DB_fdr_cutoff) %>%
    filter(DE.FDR < DE_fdr_cutoff) %>%
    mutate(DE.logFC.direction = if_else(DE.logFC > 0, 'Up', 'Down'))

  p <- ggplot() +
    geom_point(data = volcano_merged_df,
               aes(x=DB.logFC, y=negLogDBFDR),
               size = 1,
               alpha = 1,
               color = 'grey') +
    geom_point(data = volcano_DE_groups,
               aes(x=DB.logFC, y=negLogDBFDR, col=DE.logFC.direction),
               size = point_size,
               alpha = point_alpha) +
    scale_colour_manual(values = color, na.value='grey') +
    xlim(-xdiff, xdiff) +
    ylim(0, ymax) +
    labs(x = TeX('$\\log_2$( differential binding FC )'),
         y = TeX('$-\\log_{10}$( differential binding FDR )'),
         colour = 'bulkRNAseq FC',
         title = figure_title)

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
#' width, height, point_size, regression, DB_fdr_cutoff, DE_fdr_cutoff, point_color,
#' line_color) create and saves scatter plot from merged DB and DE data and save the
#' plot under designated directory.
#'
#' @param merged_df Merged data generated from running \code{merge_sites_with_DE()}
#' @param figure_title Character. The scatter plot title.
#' @param figure_save_path Character. Directory to save the figure output.
#' @param file_name Character. The base name of the scatter plot file.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param point_size Numeric. The size of the point in scatter plot. Default is 1.
#' @param regression logical value. Default is set to 'TRUE', a regression line is
#' drawn to the scatter plot. otherwise no regression line.
#' @param DB_fdr_cutoff Numeric. Significance threshold for FDR in differential binding
#' dataframe. Default is 0.05.
#' @param DE_fdr_cutoff Numeric. Significance threshold for FDR in RNA-seq differential
#' expression. Defualt is 0.05.
#' @param point_color Character Value. The color of the point. Default is \code{'black'}.
#' @param line_color Character Value. The color of the regression line. Default is
#' \code{'aquamarine4'}.
#' @import tools
#' @return ggplot object of scatter plot
#' @export
#' @examples
#' make_scatter_plot_from_merged(
#'   merged_df = merged_df,
#'   figure_title = "DB vs DE log2FC",
#'   figure_save_path = "diffbind_figures",
#'   file_name = "db_DE_scatter")

make_scatter_plot_from_merged <- function(
    merged_df,
    figure_title,
    figure_save_path,
    file_name,
    width=8,
    height=6,
    point_size=1,
    regression=TRUE,
    DB_fdr_cutoff=0.05,
    DE_fdr_cutoff=1,
    point_color='black',
    line_color='aquamarine4'
){
  # Check for required columns
  required_cols <- c("DB.logFC", "DB.FDR", "DE.logFC", "DE.FDR")
  missing_cols <- setdiff(required_cols, colnames(merged_df))
  if (length(missing_cols) > 0) {
    stop("The merged data is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  #strip off extensions in the filename
  if (tools::file_ext(file_name) != "") {
    file_name <- tools::file_path_sans_ext(file_name)
    warning("file name should not include any extensions, file name changed to: ", file_name)
  }

  # exclude rows with NA Values in any column
  sig_merged_df <- merged_df %>%
    na.omit() %>%
    filter(DB.FDR < DB_fdr_cutoff) %>%
    filter(DE.FDR < DE_fdr_cutoff)

  p <- sig_merged_df %>%
    ggplot(aes(x = DB.logFC, y = DE.logFC)) +
    geom_point(color = point_color,
               size = point_size)

  if (regression) {
    db.DE.scatter.regression = glm(DE.logFC ~ DB.logFC, data = sig_merged_df)
    beta = summary(db.DE.scatter.regression)$coefficients['DB.logFC', 1]
    p.val = summary(db.DE.scatter.regression)$coefficients['DB.logFC', 4]
    p <- p + geom_smooth(method = glm, color = line_color)+
      xlim(c(-1,1) * max(abs(merged_df$DB.logFC))) +
      ylim(c(-1,1) * max(abs(merged_df$DE.logFC))) +
      labs(x = TeX('Differential Binding Region $\\log_2$FC'),
           y = TeX('Nearby Gene RNAseq $\\log_2$FC'),
           title = figure_title)

    cat(paste0("beta = ", beta, "\np value = ", p.val, "\n"))

  } else {
    p <- p +
      xlim(c(-1,1) * max(abs(sig_merged_df$DB.logFC))) +
      ylim(c(-1,1) * max(abs(sig_merged_df$DE.logFC))) +
      labs(x = TeX('Differential Binding Region $\\log_2$FC'),
           y = TeX('Nearby Gene RNAseq $\\log_2$FC'),
           title = figure_title)
  }

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(figure_save_path, file_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  return(p)
}
