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
#' @importFrom tibble rownames_to_column
#' @importFrom latex2exp TeX
#' @importFrom gtools mixedorder

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
#'  \item \strong{"SampleID"}
#'  \item \strong{"Condition"}
#'  \item \strong{"bamReads"}
#'  \item \strong{"Peaks"}
#'  \item \strong{"PeakCaller"}
#'  }
#' @param control_level Name of reference group (e.g., "HET", "Ctrl").
#' @param min_members Minimum number of replicates per condition.
#' Set to 2 if any condition has only 2 replicates, otherwise set to \code{NULL}.
#' @param analysis_method The analysis method, either \code{DBA_DESEQ2} or \code{DBA_EDGER}.
#' @param fdr_threshold FDR cutoff for significance, used in both filtering and visualization.
#' @param normalization Normalization method. One of \code{DBA_NORM_LIB}, \code{DBA_NORM_RLE},
#' \code{DBA_NORM_TMM}, etc. See DiffBind documentation.
#' @return analyzed diffbind object
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
    fdr_threshold = 0.05,
    normalization = DBA_NORM_LIB
){
  if (is.null(sample_sheet_path) || !file.exists(sample_sheet_path)) {
    stop("'sample_sheet_path' is invalid or the file does not exist.")
  }

  samplesheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)
  required_columns <- c("SampleID", "Condition", "bamReads", "Peaks", "PeakCaller")

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
  dba.obj.samples$config$DataType <- DBA_DATA_FRAME
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
#'
#' save_diffbind_object(dba_object, file_suffix, save_directory) saves the
#' DBA object under the specified file path. It returns nothing
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param file_suffix Character. Suffix for the output DBA object file.
#' @param save_directory Character. Directory path to save the resulting
#' DBA object. If \code{NULL}, current working directory is used.
#' @export
#' @examples
#' save_diffbind_object(
#'   dba_object = dba.obj,
#'   file_suffix = "h3k27me3",
#'   save_directory = "test_results")

save_diffbind_object <- function(
    dba_object,
    file_suffix,
    save_directory = NULL
){
  if (is.null(save_directory)) {
    save_directory <- getwd()
  }

  file_path = file.path(save_directory, paste0("dba_obj_", file_suffix, ".RDS"))
  message(paste('Saving diffbind object to', file_path))

  #Save the DBA object
  saveRDS(dba_object, file=file_path)
}

#' Retrieve Differentially Bind Sites from a DBA object
#'
#' get_diffbind_sites(dba_object, fdr_threshold, contrast_number) Extracts differentially
#' bound (DB) regions from a DBA object and filters regions based on a user-defined
#' FDR threshold, finally returns a dataframe that contains sites information.
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param fdr_threshold Numeric. Significance cutoff for filtering DB sites.
#' @param contrast_number Integer. Index specifying which contrast to extract
#' results from.
#' @param verbose Boolean. If true, it prints information on number of DB sites.
#' Silent if false.
#' @return A dataframe that contains filtered sites information.
#' @export
#' @examples
#' get_diffbind_sites(dba_object = dba.obj)

get_diffbind_sites <- function(
    dba_object,
    fdr_threshold = 0.05,
    contrast_number = 1,
    verbose = TRUE
){
  if (contrast_number <= 0 || contrast_number > length(dba_object$contrasts)) {
    stop("Invalid contrast_number. Please provide a valid contrast index.")
  }

  # Initialize a flag variable to FALSE, used to track if the specific DiffBind warning
  # "No sites above threshold" was triggered.
  warning_triggered <- FALSE

  # Call dba.report() and convert its result to a data frame, while intercepting warnings.
  # The withCallingHandlers function allows us to "catch" warnings when they are issued:
  sites <- withCallingHandlers({
    as.data.frame(dba.report(dba_object, contrast = contrast_number, th = fdr_threshold))
  }, warning = function(w) {
    if (grepl("No sites above threshold", w$message)) {
      warning_triggered <<- TRUE           # set flag for custom warning later
      invokeRestart("muffleWarning")       # suppress default warning
    }
  })

  # If no sites found or fdr threshold is too stringent, produce the error
  if (is.null(sites) || warning_triggered) {
    stop(paste(
      "No significant sites were found for your comparison.",
      "This can happen if your FDR cutoff is too strict or if the contrast settings need adjustment.",
      "Try increasing the FDR threshold or reviewing your experimental contrast."
    ))
  }

  if (verbose) {
    print(paste('Positive logFC:', sum(sites$Fold > 0)))
    print(paste('Negative logFC:', sum(sites$Fold < 0)))
  }

  return(sites)
}

#' save differential binding sites files
#'
#' save_diffbind_sites(dba_object, save_directory, fdr_threshold, contrast_number) saves site
#' dataframes in bed format based on different fdr_threshold, which includes allDB
#' sites, upDB sites, and downDB sites. It doesn't return anything
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param save_directory Character. Output directory where results are written.
#' @param fdr_threshold Numeric. FDR cutoff for significantly differntially
#' bound up and down sites.
#' @param contrast_number Integer. Index specifying which contrast to extract
#' results from.
#' @export
#' @examples
#' save_diffbind_sites(
#'   dba_object = dba.obj,
#'   save_directory = "diffbind_files")

save_diffbind_sites <- function(
    dba_object,
    save_directory = NULL,
    fdr_threshold = 0.05,
    contrast_number = 1
){
  old.scipen = options('scipen')
  options(scipen = 99)

  if (is.null(save_directory)) {
    save_directory <- getwd()
  }
  if (!dir.exists(save_directory)) {
    dir.create(save_directory, recursive = TRUE)
  }

  tryCatch({
    sig.sites <- get_diffbind_sites(dba_object,
                            fdr_threshold = fdr_threshold,
                            contrast_number = contrast_number,
                            verbose = FALSE)

    #UpDB.sites sorted by FDR
    upDB.sites <- sig.sites %>%
      filter(Fold > 0) %>%
      mutate(Score = -log10(FDR)) %>%
      tibble::rownames_to_column('site.id') %>%
      arrange(FDR) %>%
      select(Chr, Start, End, site.id, Score) %>%
      mutate(Position = paste0(Chr, ':', Start, "-", End))

    upDB_sites_file_name = file.path(save_directory, 'upDB_sites.bed')
    write.table(upDB.sites, file = upDB_sites_file_name,
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

    downDB_sites_file_name = file.path(save_directory, 'downDB_sites.bed')
    write.table(DownDB.sites, file = downDB_sites_file_name,
                quote = FALSE, sep = '\t',
                row.names = FALSE, col.names = FALSE)

    message(paste(
      'Files saved:',
      upDB_sites_file_name,
      downDB_sites_file_name,
      sep='\n'))
  }, error = function(message) {
    message(message)
    message("No differential binding sites detected. Omitting files.")
  })

  all.sites <- get_diffbind_sites(dba_object,
                          fdr_threshold = 1,
                          contrast_number = contrast_number,
                          verbose = FALSE)

  # all.sites order by chromosomes and base pair
  # sorting both Chr and Start is challenging due to
  # a) default sorting of Chr yields chr1 -> chr10 .... -> chr2
  # b) dplyr::arrange discards the current group order
  Consensus.sites <- all.sites %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    select(Chr, Start, End, site.id, FDR, Score) %>%
    mutate(Position = paste0(Chr, ':', Start, "-", End)) %>%
    slice(gtools::mixedorder(Chr)) %>%
    group_by(Chr) %>%
    mutate(across(everything(), ~.[order(Start)]))

  consensus_sites_file_name = file.path(save_directory,'consensus_sites.bed')
  write.table(Consensus.sites, file = consensus_sites_file_name,
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)

  message(paste(
    'Files saved:',
    consensus_sites_file_name,
    sep='\n'))
  options(scipen = old.scipen)
}

#' Create and save a density plot
#'
#' make_diffbind_density_plot(dba_object, figure_title, save_name, fdr_threshold, save_directory,
#' contrast_number, xdiff, width, height, color) creates a density plot from a DBA object
#' with specified figure title and saves it in pdf and png format under the designated
#' file path.
#'
#' @param dba_object DBA object obtained after \code{perform_diffbind()}.
#' @param figure_title Character. Title for density plot.
#' @param save_name Character. Base name used for saving the plot files.
#' @param fdr_threshold Numeric. FDR cutoff for DB significance.
#' @param save_directory Character. Directory path where the plot files will be saved.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Setting to
#' \code{NULL} frees the limits.
#' @param width Numeric. Width dimensions of saved plot in inches.
#' @param height Numeric. Height dimensions of saved plot in inches.
#' @param color Character. Color to fill the density plot.
#' @return A ggplot object of density plot.
#' @keywords Density Plot
#' @export
#' @examples
#' make_diffbind_density_plot(
#'   dba_object = dba.obj,
#'   figure_title = "Distribution of Differential Binding Regions",
#'   save_name = "df_fold_density_h3k27ac",
#'   save_directory = "diffbind_figures")

make_diffbind_density_plot <- function(
    dba_object,
    figure_title,
    save_name,
    fdr_threshold = 1,
    save_directory = NULL,
    contrast_number=1,
    xdiff=3,
    width=8,
    height=6,
    color='aquamarine4'
){
  if (is.null(save_name)) {
    stop("A valid file name is required to be specified.")
  }

  if (is.null(save_directory)) {
    save_directory <- getwd()
  }

  sites <- get_diffbind_sites(dba_object = dba_object,
                      fdr_threshold = fdr_threshold,
                      contrast_number = contrast_number,
                      verbose = FALSE)

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
    filename = file.path(save_directory, save_name),
    width = width, height = height, dpi=600,
    plot = p
  )

  return(p)
}

#' Create and save a PCA plot
#'
#' make_diffbind_PCA_plot(dba_object, figure_title_nocontrast, figure_title_contrast,
#' save_directory, save_name, point_size, width, height, color) Creates side-by-side PCA
#' plots to compare the variance structure of all consensus binding regions vs differentially
#' bound (DB) regions from a DBA object. Both plots are saved as PNG and PDF files under
#' designated file path.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param figure_title_nocontrast Character. Title for the PCA plot using all consensus
#' regions.
#' @param save_directory Character. Directory path where output plots will be saved.
#' @param save_name Character. Base name (without extension) for the saved plot files.
#' @param figure_title_contrast Character. Title for the PCA plot using only differentially
#' bound regions.
#' @param point_size Numeric. Size of the points.
#' @param width Numeric. Width dimensions of saved plot in inches.
#' @param height Numeric. Height dimensions of saved plot in inches.
#' @param color Character vector of length 2. Colors to use for the two biological conditions
#' in PCA plots.
#' @return A \code{gridExtra} object containing a side-by-side layout of two PCA plots.
#' @keywords PCA plot
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' make_diffbind_PCA_plot(
#'   dba_object = dba.obj,
#'   figure_title_nocontrast = "PCA of All Consensus Regions",
#'   figure_title_contrast = "PCA of Differential Binding Regions",
#'   save_directory = "diffbind_figures",
#'   save_name = "diffbind_pca_plots")

make_diffbind_PCA_plot <- function(
    dba_object,
    figure_title_nocontrast,
    save_name,
    figure_title_contrast=NULL,
    save_directory = NULL,
    point_size=8,
    width=20,
    height=9,
    color=c('darkmagenta','aquamarine4')
){
  if (is.null(figure_title_contrast)) {
    figure_title_contrast <- paste(figure_title_nocontrast, "Differential Binding Sites")
  }

  if (is.null(save_directory)) {
    save_directory = getwd()
  }

  if (length(color) != 2) {
    stop("There should be 2 color choices.")
  }

  # For all consensus regions:
  df_nocontrast_pre <- data.frame(dba.peakset(dba_object, bRetrieve = TRUE))
  df_nocontrast <- t(df_nocontrast_pre[,6:ncol(df_nocontrast_pre)])
  df_label <- dba_object$samples$Condition
  p1 <- kwanlibr::draw_PCA(df_nocontrast, label=df_label, color=color) +
    theme_bw() +
    geom_point(size = point_size) +
    ggtitle(figure_title_nocontrast) +
    theme(plot.title = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size=20))

  # For all DB regions:
  df_contrast_pre <- dba.report(dba_object, bCounts = TRUE)
  df_contrast <- t(df_contrast_pre[,10:ncol(df_contrast_pre)])
  rownames(df_contrast) <- NULL
  p2 <- kwanlibr::draw_PCA(df_contrast, label = df_label, color=color) +
    theme_bw() +
    geom_point(size = point_size) +
    ggtitle(figure_title_contrast) +
    theme(plot.title = element_text(size=20),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size=20))

  #combine two plots into one
  p <- gridExtra::grid.arrange(p1, p2, nrow=1)

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(save_directory, save_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  return(p)
}

#' Generate volcano-plot-ready data from DiffBind results
#'
#' get_diffbind_volcano_data(dba_object, contrast_number, fdr_threshold, xdiff, ymax) makes
#' volcano-plottable table and applies filtering based on FDR threshold and clamps
#' values.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object.
#' @param fdr_threshold Numeric. FDR cutoff for significance.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}.
#' @return A dataframe that contains volcano-plottable data.
#' @export
#' @examples
#' get_diffbind_volcano_data(dba_object = dba.obj)

get_diffbind_volcano_data <- function(
    dba_object,
    contrast_number=1,
    fdr_threshold=0.05,
    xdiff=3,
    ymax=20
){
  sites <- get_diffbind_sites(dba_object,
                      fdr_threshold = 1,
                      contrast_number = contrast_number,
                      verbose = FALSE)
  sig.sites <- tryCatch(
    get_diffbind_sites(
      dba_object,
      fdr_threshold = fdr_threshold,
      contrast_number = contrast_number,
      verbose = FALSE)
    , error = function(message) {
      message(message)
      return(head(sites, 0)) # If no signifcant sites, stack an empty dataframe
    }
  )

  downsampled_nonsig_sites <- sites %>%
    filter(FDR > fdr_threshold) %>%
    sample_frac(1, replace=FALSE)

  volcano.sites <- sig.sites %>%
    bind_rows(downsampled_nonsig_sites) %>%
    mutate(negLogFDR = -log10(FDR)) %>%
    mutate(Fold = kwanlibr::clamp(Fold, -xdiff, xdiff),
           negLogFDR = kwanlibr::clamp(negLogFDR, 0, ymax)) %>%
    mutate(FDR_is_significant = if_else(FDR <= fdr_threshold,
                                        paste('<=',fdr_threshold),
                                        paste('>',fdr_threshold))) %>%
    mutate(FDR_is_significant = factor(FDR_is_significant, 
                                       levels = paste(c('>','<='), fdr_threshold)))
  return(volcano.sites)
}

#' Create and save a volcano plot
#'
#' make_diffbind_volcano_plot(dba_object, figure_title, save_name, save_directory,
#' point_size, point_alpha, width, height, xdiff, ymax, contrast_number, fdr_threshold, color)
#' generates volcano plot and returns the ggplot object and it also saves the file
#' under the designated figure file path.
#'
#' @param dba_object DBA object (RDS format) generated by \code{perform_diffbind()}.
#' @param figure_title Character. Figure title for volcano plot, ex. "h3k27me3 cKO vs cHET"
#' @param save_name Character. Base name for the output volcano plot figures.
#' @param save_directory Character. Directory path where the output plots will be saved.
#' @param point_size Numeric. Size of the points in volcano plots.
#' @param point_alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque).
#' @param width Numeric. Width dimensions of saved plot in inches.
#' @param height Numeric. Height dimensions of saved plot in inches.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}.
#' @param contrast_number Integer. Index indicating which contrast to extract from the DBA object.
#' @param fdr_threshold Numeric. The FDR threshold used for grouping.
#' @param color Character vector of length 2. First element is the color for insignificant DB sites.
#' Second element is the color for significant DB sites.
#' @return A ggplot object of volcano plot
#' @keywords volcano plot
#' @export
#' @examples
#' make_diffbind_volcano_plot(
#'   dba_object = dba.obj,
#'   figure_title = "h3k27ac cKO vs cHET",
#'   save_directory = "diffbind_figures",
#'   save_name = "db_volcano")

make_diffbind_volcano_plot <- function(
    dba_object,
    figure_title,
    save_name,
    save_directory=NULL,
    point_size=1,
    point_alpha=1,
    width=8,
    height=6,
    xdiff=3,
    ymax=20,
    contrast_number=1,
    fdr_threshold=0.05,
    color=c('grey', 'aquamarine4')
){

  if (is.null(save_directory)) {
    save_directory = getwd()
  }

  if (length(color) != 2) {
    stop("There should be 2 color choices.")
  }

  volcano.sites <- get_diffbind_volcano_data(dba_object,
                                      contrast_number = contrast_number,
                                      fdr_threshold = fdr_threshold,
                                      xdiff = xdiff,
                                      ymax = ymax)

  p <- ggplot(data = volcano.sites, aes(x = Fold, y = negLogFDR, color = FDR_is_significant)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = color) +
    xlim(-xdiff, xdiff) +
    ylim(0, ymax) +
    theme(aspect.ratio = 1.0) +
    labs(x = TeX('$\\log_2$ FC'),
         y = TeX('$-\\log_{10}$ FDR'),
         title = figure_title)

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(save_directory, save_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  return(p)
}

#' Merge Differential Binding and Expression Data
#'
#' merge_diffbind_with_DE(dba_object, DE_file_path, tss_file_path, save_directory, contrast_number)
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
#' @param save_directory Character. Directory path for storing intermediate files, including the TSS
#' BED file, sorted BED file, and nearest gene lookup TSV.
#' @param contrast_number Integer. Index indicating which contrast to extract from the
#' DBA object.
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
#' merge_diffbind_with_DE(
#'   dba_object = dba.obj,
#'   DE_file_path = path_to_DE_csv_file,
#'   tss_file_path = path_to_tss_file,
#'   save_directory = "diffbind_files")

merge_diffbind_with_DE <- function(
    dba_object,
    DE_file_path,
    tss_file_path,
    save_directory,
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

  if (!dir.exists(save_directory)) {
    dir.create(save_directory, recursive = TRUE)
  }

  sites <- get_diffbind_sites(dba_object = dba_object,
                       fdr_threshold = 1,
                       contrast_number = contrast_number,
                       verbose = FALSE)

  allDB.sites = sites %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    select(Chr, Start, End, site.id, Score)

  old.scipen = options('scipen')
  options(scipen = 99) # (practically) disable scientific notation in output
  write.table(allDB.sites, file = file.path(save_directory, 'allDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)
  options(scipen = old.scipen)

  # run bash commands to retrieve the nearest gene ID
  allDB_bed <- file.path(save_directory, "allDB_sites.bed")
  lookup_tsv <- file.path(save_directory, "DB_site_nearest_gene_lookup.tsv")
  sorted_bed <- tempfile(pattern = "sorted_allDB_sites_", tmpdir = save_directory, fileext = ".bed")

  system(paste("module load Bioinformatics bedops/2.4.41 &&",
               "sort-bed", allDB_bed, ">", sorted_bed,
               "&& closest-features --closest --delim '\t'", sorted_bed,
               tss_file_path,
               "| awk '{ print $4 \"\\t\" $9 }' >", lookup_tsv), intern = TRUE)

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
#' make_diffbind_volcano_plot_from_merged(merged_df, figure_title, save_directory, save_name,
#' xdiff, ymax, point_size, point_alpha, width, height, DB_fdr_cutoff, DE_fdr_cutoff,
#' color) creates and saves a volcano plot based on merged DB and DE data. All sites
#' are shown in gray; DB sites mapped to significant DE genes are colored by direction
#' of DE logFC (Up/Down).
#'
#' @param merged_df Merged data generated from running \code{merge_diffbind_with_DE()}
#' @param figure_title Character. The title of the volcano plot.
#' @param save_directory Character. The path to save the volcano plot.
#' @param save_name Character. Output figure save_name (e.g., "volcano.png").
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}.
#' @param point_size Numeric. Point size in the plot.
#' @param point_alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque).
#' @param width Numeric. Width dimensions of saved plot in inches.
#' @param height Numeric. Height dimensions of saved plot in inches.
#' @param DB_fdr_cutoff Numeric. Significance threshold for FDR in differential binding
#' dataframe.
#' @param DE_fdr_cutoff Numeric. Significance threshold for FDR in RNA-seq differential
#' expression.
#' @param color Character vector of length 2, in order of down FC color, then up FC color.
#' @return A ggplot object of DB and DE merged volcano plot.
#' @export
#' @examples
#' make_diffbind_volcano_plot_from_merged(
#'   merged_df = merged_volcano,
#'   figure_title = 'cKO vs cHet Differential Binding sites',
#'   save_directory = "diffbind_figures",
#'   save_name = "db_volcano_DE_color_binary")

make_diffbind_volcano_plot_from_merged <- function(
    merged_df,
    figure_title,
    save_directory,
    save_name,
    xdiff=3,
    ymax=20,
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
    theme(aspect.ratio = 1.0) +
    labs(x = TeX('$\\log_2$( differential binding FC )'),
         y = TeX('$-\\log_{10}$( differential binding FDR )'),
         colour = 'bulkRNAseq FC',
         title = figure_title)

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(save_directory, save_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  return(p)
}

#' Create and return a GLM object from DB and DE Data
#'
#' DB_DE_regression_model(merged_df) fits a generalized linear model (GLM) to test
#' the association between DE and DB data. It regresses \code{DE.logFC} against
#' \code{DB.logFC} in a merged data frame.
#'
#' @param merged_df Merged data generated from running \code{merge_diffbind_with_DE()}.
#' @return An object of class inheriting from 'glm'.
#' @details The summary of the regression is printed to the console.
#' @export
#' @examples
#' DB_DE_regression_model(merged_df = merged_df)

DB_DE_regression_model <- function(
    merged_df
){
  regression <- glm(DE.logFC ~ DB.logFC, data = merged_df)

  print(summary(regression))
  return(regression)
}

#' Create and save the scattor plot from DB and DE Data
#'
#' make_scatter_plot_from_merged(merged_df, figure_title, save_directory, save_name,
#' width, height, point_size, regression, DB_fdr_cutoff, DE_fdr_cutoff, point_color,
#' line_color) create and saves scatter plot from merged DB and DE data and save the
#' plot under designated directory.
#'
#' @param merged_df Merged data generated from running \code{merge_diffbind_with_DE()}
#' @param figure_title Character. The scatter plot title.
#' @param save_directory Character. Directory to save the figure output.
#' @param save_name Character. The base name of the scatter plot file.
#' @param width Numeric. Width dimensions of saved plot in inches.
#' @param height Numeric. Height dimensions of saved plot in inches.
#' @param point_size Numeric. The size of the point in scatter plot.
#' @param regression logical value of whether a regression line is requested on the scatter plot.
#' @param DB_fdr_cutoff Numeric. Significance threshold for FDR in differential binding
#' dataframe.
#' @param DE_fdr_cutoff Numeric. Significance threshold for FDR in RNA-seq differential
#' expression. when setting it to 1, it will include all the genes so we can look for overall
#' trends.
#' @param point_color Character Value. The color of the point.
#' @param line_color Character Value. The color of the regression line.
#' @import tools
#' @return ggplot object of scatter plot
#' @export
#' @examples
#' make_scatter_plot_from_merged(
#'   merged_df = merged_df,
#'   figure_title = "DB vs DE log2FC",
#'   save_directory = "diffbind_figures",
#'   save_name = "db_DE_scatter")

make_scatter_plot_from_merged <- function(
    merged_df,
    figure_title,
    save_directory,
    save_name,
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
    # Since there is no straightforward way to access the model object within geom_smooth(),
    # we fit the model externally to retrieve summary details from the glm object.
    DB.DE.scatter.regression = DB_DE_regression_model(sig_merged_df)
    p <- p + geom_smooth(method = glm, color = line_color)
  }

  p <- p +
    xlim(c(-1,1) * max(abs(sig_merged_df$DB.logFC))) +
    ylim(c(-1,1) * max(abs(sig_merged_df$DE.logFC))) +
    theme(aspect.ratio = 1.0) +
    labs(x = TeX('Differential Binding Region $\\log_2$FC'),
         y = TeX('Nearby Gene RNAseq $\\log_2$FC'),
         title = figure_title)

  #save the plot
  kwanlibr::ggsave_vector_raster(
    filename = file.path(save_directory, save_name),
    width = width, height = height, dpi = 600,
    plot = p
  )

  return(p)
}
