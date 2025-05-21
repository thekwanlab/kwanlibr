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
#' perform_diffbind(samplesheetpath, reference, minMembers, analysismethod, datatype,
#' FDRthreshold, normalization, filesuffix, filePath) performs a complete DiffBind
#' analysis workflow using a provided sample sheet. It includes blacklist filtering,
#' counting, normalization, and contrast setup. The resulting DBA object is saved as
#' an RDS file under specific file path.
#'
#' @param samplesheetpath Character. The path to the sample CSV file.
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
#' @param reference Character. Reference group name (e.g., "HET", "Ctrl").
#' @param minMembers Integer. Minimum number of replicates per condition (default is \code{NULL}).
#' Set to 2 if any condition has only 2 replicates, otherwise set to \code{NULL}.
#' @param analysismethod The analysis method, either \code{DBA_DESEQ2} (default) or \code{DBA_EDGER}.
#' @param datatype Output format of peak data (default is \code{DBA_DATA_FRAME}).
#' Other options are \code{DBA_DATA_GRANGES}, \code{DBA_DATA_RANGEDDATA}.
#' @param FDRthreshold Numeric. FDR cutoff for significance, used in both filtering
#' and visualization (Default is 0.05).
#' @param normalization Normalization method (default is \code{DBA_NORM_LIB}). Other options
#' include \code{DBA_NORM_RLE}, \code{DBA_NORM_TMM}, etc.
#' @param filesuffix Character. Suffix for the output DBA object file.
#' @param filePath Character. Directory path to save the resulting DBA object (default is \code{NULL}).
#' If \code{NULL}, current working directory is used.
#' @export
#' @examples
#' perform_diffbind(samplesheetpath = "diffbind_samples.csv",
#'                  reference = "HET",
#'                  minMembers = 2,
#'                  filesuffix = "h3k27me3",
#'                  filepath = "diffbind_files")

perform_diffbind <- function(
    samplesheetpath,
    reference,
    minMembers=NULL,
    analysismethod=DBA_DESEQ2,
    datatype = DBA_DATA_FRAME,
    FDRthreshold = 0.05,
    normalization = DBA_NORM_LIB,
    filesuffix=NULL,
    filePath=NULL
) {
  if (is.null(samplesheetpath) || !file.exists(samplesheetpath)) {
    stop("The 'samplesheetpath' is invalid or the file does not exist.")
  }
  samplesheet <- read.csv(samplesheetpath, stringsAsFactors = FALSE)

  required_columns <- c("Sample_ID", "Tissue", "Condition", "bamReads",
                        "Peaks", "PeakCaller")
  missing_columns <- setdiff(required_columns, colnames(samplesheet))

  if (length(missing_columns) > 0) {
    stop("The sample sheet is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  if (!analysismethod %in% c(DBA_DESEQ2, DBA_EDGER)) {
    stop("Invalid 'analysismethod'. Choose either 'DBA_DESEQ2' or 'DBA_EDGER'.")
  }
  if (!reference %in% samplesheet$Condition) {
    stop("the specified reference'", reference, "'is not found in the Condition")
  }
  if (is.null(filePath)) {
    filePath <- getwd()
  }
  if (is.null(filesuffix)){
    stop("'filesuffix' is required to uniquely identify the saved RDS file.")
  }

  #Generate dba object
  dba.obj.samples = DiffBind::dba(sampleSheet = samplesheetpath)
  dba.obj.samples$config$th = FDRthreshold
  dba.obj.samples$config$DataType = datatype
  dba.obj.samples$config$AnalysisMethod = analysismethod

  dba.obj.blacklist = DiffBind::dba.blacklist(dba.obj.samples)
  dba.obj.count = DiffBind::dba.count(dba.obj.blacklist, summits = TRUE)
  dba.obj.normalize = DiffBind::dba.normalize(dba.obj.count, normalize=normalization)

  if (minMembers == 2){
    dba.obj.contrast = DiffBind::dba.contrast(dba.obj.normalize,
                                              minMembers = minMembers,
                                              reorderMeta = list(Condition=reference))
  } else if (is.null(minMembers)) {
    dba.obj.contrast = DiffBind::dba.contrast(dba.obj.normalize,
                                              reorderMeta = list(Condition=reference))
  }
  dba.obj = DiffBind::dba.analyze(dba.obj.contrast)

  #Save the DBA object
  saveRDS(dba.obj, file=paste0(filePath, "dba_obj_", filesuffix, ".RDS"))

  message("DiffBind analysis completed successfully. Results saved to:",
          paste0(filePath, paste0("dba_obj_", filesuffix, ".RDS")))
}

#' Retrieve Differentially Bind Sites from a DBA object
#'
#' get_DBsites(dbaobjpath, FDR_THRESHOLD, contrastnumber, filepath) Extracts differentially
#' bound (DB) regions from a DBA object and writes the results into BED-formatted files. The
#' function filters regions based on a user-defined FDR threshold and saves separate files for
#' all significant sites, upregulated, and downregulated sites.
#'
#' @param dbaobjpath Character. File path to the DBA object (RDS file) created from \code{perform_diffbind()}
#' @param FDR_THRESHOLD Numeric. Significance cutoff for filtering DB sites. Default is \code{0.05}.
#' @param contrastnumber Integer. Index specifying which contrast to extract results from. Defaults to \code{1},
#' corresponding to the first contrast.
#' @param filepath Character. Output directory where results are written.
#' @export
#' @examples
#' get_DBsites(dbaobjpath=dbaobjpath)

get_DBsites <- function(
    dbaobjpath,
    FDR_THRESHOLD=0.05,
    contrastnumber=1,
    filepath=NULL
) {
  if (is.null(dbaobjpath)|| !file.exists(dbaobjpath)){
    stop("check again on dbaobj to see if it is a valid path and file.")
  }

  dba.obj <- tryCatch({
    readRDS(dbaobjpath)
  }, error = function(e) {
    stop("Failed to read the DBA object from path: ", dbaobjpath, "\nError: ", e$message)
  })

  if (contrastnumber <= 0 || contrastnumber > length(dba.obj$contrasts)) {
    stop("Invalid contrastnumber. Please provide a valid contrast index.")
  }
  if (is.null(filepath)) {
    filepath <- getwd()
  }

  #read DBA object
  dba.obj = readRDS(dbaobjpath)
  sites = as.data.frame(dba.report(dba.obj, contrast = contrastnumber, th = 1))
  if (is.null(sites)) {
    stop("No sites found in the report. Please check your contrast setup.")
  }

  #sig.sites
  sig.sites = sites %>% filter(FDR <= FDR_THRESHOLD)
  write.table(sig.sites, file = file.path(filepath,'sig_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = TRUE, col.names = TRUE)
  print(head(sig.sites))
  print(paste('Positive logFC:', sum(sig.sites$Fold > 0)))
  print(paste('Negative logFC:', sum(sig.sites$Fold < 0)))
  if (!dir.exists(filepath)) {
    dir.create(filepath, recursive = TRUE)
  }
  #allDB.sites order by chromosomes
  allDB.sites = sig.sites %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    select(Chr, Start, End, site.id, FDR, Score) %>%
    mutate(Position = paste0(Chr, ':', Start, "-", End)) %>%
    slice(gtools::mixedorder(Chr), Start)
  write.table(allDB.sites, file = file.path(filepath,'allDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)
  #UpDB.sites sorted by FDR
  upDB.sites = sig.sites %>%
    filter(Fold > 0) %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    arrange(FDR) %>%
    select(Chr, Start, End, site.id, Score) %>%
    mutate(Position = paste0(Chr, ':', Start, "-", End))
  write.table(upDB.sites, file = file.path(filepath, 'upDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)
  #DownDB.sites sorted by FDR
  DownDB.sites = sig.sites %>%
    filter(Fold < 0) %>%
    mutate(Score = -log10(FDR)) %>%
    tibble::rownames_to_column('site.id') %>%
    arrange(FDR) %>%
    select(Chr, Start, End, site.id, Score) %>%
    mutate(Position = paste0(Chr, ':', Start, "-", End))
  write.table(DownDB.sites, file = file.path(filepath, 'downDB_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = FALSE, col.names = FALSE)
}

#' Create and save a density plot
#'
#' make_density_plot(dbaobjpath, figure_title, filename, figure_dir, contrast_number,
#' FDR_THRESHOLD, width, height, DENSITY_color) creates a density plot from a DBA
#' object with specified figure title and saves it in pdf and png format under the
#' designated file path.
#'
#' @param dbaobjpath Character. Path to the saved DBA object (RDS format) from \code{perform_diffbind()}.
#' @param figure_title Character. Title displayed at the top of the plot.
#' @param filename Character. Base name (without extension) used for saving the plot files.
#' @param figure_dir Character. Directory path where the plot files will be saved. Default
#' is set to \code{NULL}.
#' @param contrastnumber Integer. Index indicating which contrast to extract from the DBA
#' object. Defeault is set to \code{1}.
#' @param FDR_THRESHOLD Numeric. The FDR threshold used for annotation purposes in the plot
#' title. Default is \code{0.05}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param DENSITY_color Character. Color to fill the density plot. Default is \code{"aquamarine4"}.
#' @return A \code{ggplot2} object containing the generated density plot.
#' @keywords Density Plot
#' @export
#' @examples
#' make_density_plot(
#'   dbaobjpath=dbaobjpath,
#'   figure_title="Distribution of Differential Binding Regions",
#'   figure_dir="diffbind_figures"),
#'   filename="df_fold_density_h3k27ac")

make_density_plot <- function(
    dbaobjpath,
    figure_title,
    filename,
    figure_dir=NULL,
    contrastnumber=1,
    FDR_THRESHOLD=0.05,
    width=8,
    height=6,
    DENSITY_color='aquamarine4'
) {
  if (is.null(dbaobjpath)|| !file.exists(dbaobjpath)){
    stop("check again on dbaobj to see if it is a valid path and file.")
  }
  if (is.null(figure_dir)) {
    figure_dir <- getwd()
  }
  #read DBA object
  dba.obj = readRDS(dbaobjpath)
  sites = as.data.frame(DiffBind::dba.report(dba.obj, contrast=contrastnumber, th=1))

  #draw density plot
  p <- sites %>%
    ggplot(aes(x=Fold)) +
    geom_density(aes(y=after_stat(count)), fill=DENSITY_color) +
    labs(x = TeX('$\\log_2$ FC'),
         title = TeX(paste0(figure_title,
                            ' [ FDR < ', FDR_THRESHOLD, ' ]')))

  #save the plot
  ggplot2::ggsave(filename = file.path(figure_dir, paste0(filename ,'.png')),
                  width=width, height=height)
  ggplot2::ggsave(filename = file.path(figure_dir, paste0(filename ,'.pdf')),
                  width=width, height=height)

  message("Density Plot is generated successfully. Results saved to: \n",
          file.path(figure_dir, paste0(filename ,'.png')), " and \n",
          file.path(figure_dir, paste0(filename ,'.pdf')))
  return(p)
}

#' Create and save a PCA plot
#'
#' make_PCA_plot_diffbind(dbaobjpath, figure_title_nocontrast, figure_title_contrast,
#' figure_dir, filename, size, width, height, color) Creates side-by-side PCA plots to
#' compare the variance structure of all consensus binding regions vs differentially
#' bound (DB) regions from a DBA object. Both plots are saved as PNG and PDF files.
#'
#' @param dbaobjpath Character. Path to the saved DBA object (RDS format) generated by
#' \code{perform_diffbind()}.
#' @param figure_title_nocontrast Character. Title for the PCA plot using all consensus
#' regions.
#' @param figure_title_contrast Character. Title for the PCA plot using only differentially
#' bound regions.
#' @param figure_dir Character. Directory path where output plots will be saved.
#' @param filename Character. Base name (without extension) for the saved plot files.
#' @param size Numeric. Size of the points, default is set to 8.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 20.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 9.
#' @param color Character vector of length 2. Colors to use for the two biological conditions
#' in PCA plots. Default is \code{c("darkmagenta", "aquamarine4")}.
#' @return A \code{gridExtra} object containing a side-by-side layout of two PCA plots.
#' @import gridExtra
#' @export
#' @examples
#' make_PCA_plot_diffbind(
#'   dbaobjpath = dbaobjpath,
#'   figure_title_nocontrast = "PCA of All Consensus Regions",
#'   figure_title_contrast = "PCA of Differential Binding Regions",
#'   figure_dir = "diffbind_figures",
#'   filename = "diffbind_pca_plots",
#'   color = c("blue", "red"))

make_PCA_plot_diffbind <- function(
    dbaobjpath,
    figure_title_nocontrast,
    figure_title_contrast,
    figure_dir,
    filename,
    size=8,
    width=20,
    height=9,
    color=c('darkmagenta','aquamarine4')
){
  if (is.null(dbaobjpath)|| !file.exists(dbaobjpath)){
    stop("check again on dbaobj to see if it is a valid path and file.")
  }
  if (is.null(figure_dir)) {
    figure_dir <- getwd()
  }

  dba.obj = readRDS(dbaobjpath)
  # For all consensus regions:
  df_nocontrast_pre <- data.frame(dba.peakset(dba.obj, bRetrieve = TRUE))
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
  df_contrast_pre <- dba.report(dba.obj, bCounts = TRUE)
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
  ggplot2::ggsave(plot = p, filename = file.path(figure_dir, paste0(filename, '.png')),
                  width=width, height=height)
  ggplot2::ggsave(plot = p, filename = file.path(figure_dir, paste0(filename, '.pdf')),
                  width=width, height=height)

  message("PCA Plot is generated successfully. Results saved to: \n",
          file.path(figure_dir, paste0(filename ,'.png')), " and \n",
          file.path(figure_dir, paste0(filename ,'.pdf')))
  return(p)
}

#' Generate volcano-plot-ready data from DiffBind results
#'
#' make_volcano_sites(dbaobjpath, sigsitespath, filepath, FDR_THRESHOLD, xdiff, ymax)
#' makes volcano-plottable table and applies filtering based on FDR threshold and clamps
#' values.
#' It does not return the data.
#' @param dbaobjpath Character. Path to the saved DBA object (RDS format) generated by
#' \code{perform_diffbind()}.
#' @param sigsitespath Character. Path to the sig_sites.bed generated by \code{get_DBsites()}.
#' @param filepath Character. Output directory where volcano data is written.
#' @param FDR_THRESHOLD Numeric. FDR cutoff for significance. Default is set to 0.05.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @export
#' @examples
#' make_volcano_sites(dbaobjpath=dbaobjpath,
#'                    sigsitespath=sigsitespath,
#'                    filepath="diffbind_files")

make_volcano_sites <- function(dbaobjpath,
                               sigsitespath,
                               filepath,
                               FDR_THRESHOLD=0.05,
                               xdiff=5,
                               ymax=40) {
  if (is.null(dbaobjpath)|| !file.exists(dbaobjpath)){
    stop("check again on dbaobj to see if it is a valid path and file.")
  }
  if (is.null(filepath)) {
    filepath <- getwd()
  }

  dba.obj = readRDS(dbaobjpath)
  sites = dba.report(dba.obj, contrast = 1, th = 1)

  sig.sites <- read.table(sigsitespath,
                           header=TRUE, sep="\t",
                           stringsAsFactors=FALSE, quote="")

  downsampled_nonsig_sites = sites %>%
    filter(FDR > FDR_THRESHOLD) %>%
    sample_frac(1, replace=FALSE)

  volcano.sites = sig.sites %>%
    bind_rows(downsampled_nonsig_sites) %>%
    mutate(negLogFDR = -log10(FDR)) %>%
    mutate(Fold = kwanlibr::clamp(Fold, -xdiff, xdiff),
           negLogFDR = kwanlibr::clamp(negLogFDR, 0, ymax)) %>%
    mutate(FDR = if_else(FDR <= FDR_THRESHOLD,
                         paste('<=',FDR_THRESHOLD),
                         paste('>',FDR_THRESHOLD)))

  write.table(volcano.sites, file = file.path(filepath, 'volcano_sites.bed'),
              quote = FALSE, sep = '\t',
              row.names = TRUE, col.names = TRUE)
}

#' Create and save a volcano plot
#'
#' make_volcano_plot_diffbind(volcanopath, sigsitespath, figure_title, figure_dir,
#' filename, size, alpha, width, height, xdiff, ymax, FDR_THRESHOLD, color) generates
#' volcano plot from data in volcano_sites bed files and returns the ggplot object
#' also saves the file under the designated figure file path.
#'
#' @param volcanopath Character. Path to the volcano_sites.bed generated by \code{make_volcano_sites()}.
#' @param sigsitespath Character. Path to the sig_sites.bed generated by \code{get_DBsites()}.
#' @param figure_title Character. Figure title for volcano plot, ex. h3k27me3 cKO vs cHET
#' @param figure_dir Character. Directory path where the output plots will be saved.
#' @param filename Character. Base name for the output figures.
#' @param size Numeric. Size of the points, default is set to 1.
#' @param alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque). Default is \code{1}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @param FDR_THRESHOLD Numeric. The FDR threshold used for grouping, default is \code{0.05}.
#' @param color Character vector of length 2. By default the color is \code{c('aquamarine4', 'grey')}
#' @return ggplot object
#' @export
#' @examples
#' make_volcano_plot_diffbind(volcanopath=volcanopath,
#'                            sigsitespath=sigsitespath,
#'                            figure_title="h3k27ac cKO vs cHET",
#'                            figure_dir="diffbind_figures"),
#'                            filename="db_volcano")

make_volcano_plot_diffbind <- function(
    volcanopath,
    sigsitespath,
    figure_title,
    figure_dir,
    filename,
    size=1,
    alpha=1,
    width=8,
    height=6,
    xdiff=5,
    ymax=40,
    FDR_THRESHOLD=0.05,
    color=c('aquamarine4', 'grey')
    ){
  if (is.null(volcanopath)|| !file.exists(volcanopath)){
    stop("check again on dbaobj to see if it is a valid path and file.")
  }

  volcano.sites <- read.table(volcanopath,
                              header=TRUE, sep="\t",
                              stringsAsFactors=FALSE, quote="")
  sig.sites <- read.table(sigsitespath,
                          header=TRUE, sep="\t",
                          stringsAsFactors=FALSE, quote="")

  p <-draw_volcano_general(volcano.sites,
                           criteria="FDR",
                           colors = setNames(c(color[1], color[2]),
                                             c(paste('<=', FDR_THRESHOLD),
                                               paste('>', FDR_THRESHOLD))),
                           size=size,
                           alpha=alpha,
                           xdiff = xdiff,
                           ymax = ymax)
  figure_title <- paste0(figure_title, " [", dim(sig.sites)[1], ' DB Regions FDR <= ',
                  FDR_THRESHOLD,
                  ']')
  p <- p + labs(title = figure_title)
  ggplot2::ggsave(filename = file.path(figure_dir, paste0(filename, '.png')),
                  width=width, height=height)
  ggplot2::ggsave(filename = file.path(figure_dir, paste0(filename, '.pdf')),
                  width=width, height=height)
  message("Volcano Plot is generated successfully. Results saved to: \n",
          file.path(figure_dir, paste0(filename ,'.png')), " and \n",
          file.path(figure_dir, paste0(filename ,'.pdf')))
  return(p)
}

#' Merge Differential Binding and Expression Data
#'
#' merge_sites_with_exp(allDBsitepath, RNAfilepath, joindatapath, filepath, BULK_FDR_CUTOFF)
#' integrates differential binding sites with differential expression data (from bulk RNA-seq)
#' by joining them based on their nearest gene TSS and returns a data frame where each row
#' represents a differential binding site along with its nearest gene's differential expression.
#'
#' @param allDBsitepath Character. Path to the BED file containing all differentially
#' bound sites. This file is typically generated by \code{get_DBsites()}.
#' @param RNAfilepath Character. Path to the CSV file containing bulk RNA-seq differential
#' expression results. This file must include the columns \code{gene_name}, \code{logFC},
#' and \code{FDR}.
#' @param joindatapath Character. Path to the BED-like table (e.g., volcano_sites.bed) to
#' be annotated with expression data. This is the backbone data for visualization.
#' @param filepath Character. Directory path for storing intermediate files, including the TSS
#' BED file, sorted BED file, and nearest gene lookup TSV.
#' @param BULK_FDR_CUTOFF Numeric. Significance threshold for RNA-seq differential expression
#' (e.g., 0.05). Defualt is set as 0.05.
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
#' merge_sites_with_exp(allDBsitepath = "diffbind_files/allDB_sites.bed",
#'                      volcanositepath = "diffbind_files/volcano_sites.bed",
#'                      RNAfilepath = "path_to_RNAseq_csv_file",
#'                      filepath = "diffbind_files")

merge_sites_with_exp <- function(
    allDBsitepath,
    RNAfilepath,
    joindatapath,
    filepath,
    BULK_FDR_CUTOFF=0.05){
  if (is.null(allDBsitepath) || !file.exists(allDBsitepath)) {
    stop("`allDBsitepath` must be a valid path to the DB sites BED file.")
  }
  if (is.null(RNAfilepath) || !file.exists(RNAfilepath)) {
    stop("`RNAfilepath` must be a valid path to the RNA-seq results CSV file.")
  }
  if (is.null(joindatapath) || !file.exists(joindatapath)) {
    stop("`joindatapath` must be a valid path to the specified data file.")
  }
  # run bash script to retrieve the nearest gene ID
  tss_bed <- file.path(filepath, "TSS.bed")
  sorted_bed <- file.path(filepath, "allDB_sites.sorted.bed")
  lookup_tsv <- file.path(filepath, "DB_site_nearest_gene_lookup.tsv")

  script <- system.file("scripts", "db_annotate_regions.sh", package = "kwanlibr")
  if (script == "") stop("Annotation script not found in kwanlibr package.")
  system2("bash", args = c(script, allDBsitepath, tss_bed, sorted_bed, lookup_tsv))

  # read the annotated gene expression RNA sequencing file
  db.site.gene.lookup = read.table(lookup_tsv, col.names = c('site.id','gene_name'))
  bulk.DE = read.csv(RNAfilepath) %>% distinct(gene_name, .keep_all = TRUE)

  # inner join sites and the diff expressed data
  db.site.bulk.data = db.site.gene.lookup %>%
    inner_join(bulk.DE, by = 'gene_name') %>%
    select(site.id, gene_name, logFC, FDR) %>%
    dplyr::rename(bulk.logFC = logFC) %>%
    dplyr::rename(bulk.FDR = FDR)

  db.site.bulk.data.meets.FDR = db.site.bulk.data %>% filter(bulk.FDR < BULK_FDR_CUTOFF)

  data.sites <- read.table(joindatapath,
                           header=TRUE, sep="\t",
                           stringsAsFactors=FALSE, quote="")

  data.sites.bulk.color = data.sites %>%
    tibble::rownames_to_column('site.id') %>%
    mutate(site.id = as.numeric(site.id)) %>%
    left_join(db.site.bulk.data.meets.FDR, by = 'site.id')

  return(data.sites.bulk.color)
}

#' Make Volcano Plot from DB and DE Data (Color by bulk logFC direction)
#'
#' make_volcano_plot_from_merged(allDBsitepath, RNAfilepath, volcanopath, filepath,
#' figure_title, filename, color, xdiff, ymax, size, alpha, width, height, BULK_FDR_CUTOFF)
#' merges differential binding (DB) sites with RNA-seq differential expression (DE)
#' data using \code{merge_sites_with_exp()}, then creates a volcano plot. All sites
#' are shown in gray; DB sites mapped to significant DE genes are colored by direction
#' of bulk logFC (Up/Down). Plot is saved and returned.
#'
#' @param allDBsitepath Character. Path to BED file with differential binding sites.
#' @param RNAfilepath Character. Path to RNA-seq results CSV (must include gene_name, logFC, FDR).
#' @param volcanopath Character. Path to BED-like file of all sites (for background).
#' @param filepath Character. Directory to save intermediate and output files.
#' @param figure_title Character. The title of the volcano plot.
#' @param filename Character. Output figure filename (e.g., "volcano.png").
#' @param color Character vector of length 3, in order of background color, down FC color,
#' and up FC color. E.x. \code{c('grey','steelblue', 'tomato')}.
#' @param xdiff Numeric. Limits the x-axis range to \code{[-xdiff, xdiff]}. Default is \code{5}.
#' @param ymax Numeric. Limits the y-axis range to \code{[0, ymax]}. Default is \code{40}.
#' @param size Numeric. Point size in the plot. Default is \code{1}.
#' @param alpha Numeric. Transparency level of the points. Ranges from 0 (fully
#' transparent) to 1 (fully opaque). Default is \code{1}.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param BULK_FDR_CUTOFF Numeric FDR threshold to filter DE genes. Default is set
#' to 0.05.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' make_volcano_plot_from_merged(allDBsitepath = "diffbind_files/allDB_sites.bed",
#'                               RNAfilepath="",
#'                               volcanopath="diffbind_files/volcano_sites.bed",
#'                               filepath="diffbind_files",
#'                               figure_title='cKO vs cHet Differential Binding sites',
#'                               figure_path="diffbind_figures",
#'                               filename="db_volcano_bulk_color_binary",
#'                               color=c('grey','steelblue', 'tomato'))

make_volcano_plot_from_merged <- function(
    allDBsitepath,
    RNAfilepath,
    volcanopath,
    filepath,
    figure_title,
    figure_path,
    filename,
    color,
    xdiff=5,
    ymax=40,
    size=1,
    alpha=1,
    width = 8,
    height = 6,
    BULK_FDR_CUTOFF=0.05
) {
  merged_df <- merge_sites_with_exp(allDBsitepath = allDBsitepath,
                                    RNAfilepath = RNAfilepath,
                                    joindatapath = volcanopath,
                                    filepath = filepath,
                                    BULK_FDR_CUTOFF = BULK_FDR_CUTOFF)
  # Add criteria column
  merged_df <- merged_df %>%
    mutate(
      negLogFDR = -log10(FDR),
      bulk.logFC.direction = case_when(
        !is.na(bulk.FDR) & bulk.FDR < 0.05 & bulk.logFC > 0 ~ "Up",
        !is.na(bulk.FDR) & bulk.FDR < 0.05 & bulk.logFC < 0 ~ "Down",
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
                        BULK_FDR_CUTOFF, ' ]'),
         color = "bulk RNA-seq FC",
         x = TeX('$\\log_2$( differential binding FC )'),
         y = TeX('$-\\log_{10}$( differential binding FDR )'))

  ggplot2::ggsave(filename = file.path(figure_path, paste0(filename, ".png")),
                  plot = p,
                  width = width,
                  height = height)
  ggplot2::ggsave(filename = file.path(figure_path, paste0(filename, ".pdf")),
                  plot = p,
                  width = width,
                  height = height)
  return(p)
}

#' Create and save the scattor plot from DB and DE Data
#'
#' make_scatter_plot_from_merged(allDBsitepath, RNAfilepath, sigsitespath, filepath,
#' figure_title, figure_path, filename, width, height, size, regression, BULK_FDR_CUTOFF,
#' point_color, line_color) create scatter plot from merged differential binding sites
#' with differential expression data and save the plot under designated directory.
#'
#' @param allDBsitepath Character. Path to BED file of all differential binding sites.
#' @param RNAfilepath Character. Path to RNA-seq results CSV (must include gene_name, logFC, FDR).
#' @param sigsitespath Character. Path to BED file of significant sites.
#' @param filepath Character. Directory path for storing intermediate files, including
#' the TSS BED file, sorted BED file, and nearest gene lookup TSV.
#' @param figure_title Character. The scatter plot title.
#' @param figure_path Character. Directory to save the figure output.
#' @param filename Character. The name of the file.
#' @param width Numeric. Width dimensions of saved plot in inches. Default is 8.
#' @param height Numeric. Height dimensions of saved plot in inches. Default is 6.
#' @param size Numeric. The size of the point in scatter plot. Default is 1.
#' @param regression logical value. Default is set to 'TRUE', a regression line is
#' drawn to the scatter plot. otherwise no regression line.
#' @param BULK_FDR_CUTOFF Numeric FDR threshold to filter DE genes. Default is 0.05.
#' @param point_color Character Value. The color of the point. Default is \code{'black'}.
#' @param line_color Character Value. The color of the regression line. Default is
#' \code{'aquamarine4'}.
#' @export
#' @examples
#' make_scatter_plot_from_merged()

make_scatter_plot_from_merged <- function(
    allDBsitepath,
    RNAfilepath,
    sigsitespath,
    filepath,
    figure_title,
    figure_path,
    filename,
    width=8,
    height=6,
    size=1,
    regression=TRUE,
    BULK_FDR_CUTOFF=0.05,
    point_color='black',
    line_color='aquamarine4'){

  merged_df <- merge_sites_with_exp(allDBsitepath = allDBsitepath,
                                    RNAfilepath = RNAfilepath,
                                    joindatapath = sigsitespath,
                                    filepath = filepath,
                                    BULK_FDR_CUTOFF = BULK_FDR_CUTOFF)

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
  ggplot2::ggsave(filename = file.path(figure_path, paste0(filename, '.png')),
                  width=width, height=height)
  ggplot2::ggsave(filename = file.path(figure_path, paste0(filename, '.pdf')),
                  width=width, height=height)
  return(p)
}
