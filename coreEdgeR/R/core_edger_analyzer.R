#=============================================================================
# coreEdgeR is a package that contains wrapper functions
# for edgeR functions.
# They were originally for [...]

#---------
# TO-DO
#---------
# - create testing suite?
# - changelog?
# - create runable examples in documentation?
#=============================================================================

# install dependencies
if (length(find.package("librarian", quiet=TRUE)) == 0) {
  install.packages("librarian")
}

librarian::shelf(
  edgeR,
  rtracklayer,
  org.Mm.eg.db,
  data.table,
  ggrepel,
  quiet=FALSE,
  ask=FALSE
)

#constants
IMPORTANT_TAGS <- c("gene_id", "gene_type", "gene_name")


#---------------------
# Functions
#---------------------

#' Get reference
#'
#' kwanlibr_get_gtf(gtf, verbose) returns the reference specified by
#' the path gtf. 
#' @param gtf absolute path to the reference
#' @param verbose set to TRUE to print number of rows. Default is FALSE
#' @return A dataframe of the reference genome 
#' @keywords reference
#' @export
#' @examples
#' kwanlibr_get_gtf(gtf='/nfs/turbo/umms-kykwan/projects/reference/gtf/gencode.vM14.primary_assembly.ERCC.annotation.gtf')

kwanlibr_get_gtf <- function(
  gtf='/nfs/turbo/umms-kykwan/projects/reference/gtf/gencode.vM14.primary_assembly.annotation.gtf',
  verbose=FALSE
) {
  
  ## Read in reference genome for location info
  gtf <- rtracklayer::readGFF(
    gtf,
    columns=c("seqid", "start", "end"),
    tags=IMPORTANT_TAGS,
    filter=list(type=c("gene"))
  )
  
  if (verbose) {
    print(paste("Number of Rows:", nrow(gtf)))
  }
  
  return(gtf)
}


#' Perform edgeR
#'
#' kwanlibr_perform_edger(
#' sampleTable, fileCol, idCol, condCol, batchCol, filePrefix, gtf, saveName) 
#' performs edger analysis on sampleTable
#' @param sampleTable a dataframe with sample data
#' @param fileCol header that refers to file names
#' @param idCol header that refers to the unique identifier of each sample
#' @param condCol header that refers to the condition of each sample
#' @param batchCol header that refers to the batch of the sample
#' @param filePrefix prefix that the sample files start with. Default is NULL
#' @param gtf a dataframe of the reference genome
#' @param saveName filename to save edger tables under
#' @return A dataframe of the edger analysis
#' @keywords edger
#' @export
#' @examples
#' kwanlibr_perform_edger(smc3, idCol="Pool.Name", gtf=gtf, saveName="tables/Smc3EdgeR")

kwanlibr_perform_edger <- function(
  sampleTable,
  fileCol="Filename", 
  idCol="Sample", 
  condCol="condition",
  batchCol=NULL,
  filePrefix=NULL, gtf=NULL,
  saveName=NULL
) {
  if (class(filePrefix) == "NULL") {
    filePrefix <- getwd()
  }
  filePrefix <- gsub("/+$", "", filePrefix)
  
  if (class(gtf) == "NULL") {
    gtf <- kwanlibr_get_gtf()
  }
    
  if(class(batchCol) == "NULL"){
    sampleTable <- data.frame(
      sampleName = sampleTable[[idCol]],
      files = file.path(filePrefix, sampleTable[[fileCol]]),
      condition = sampleTable[[condCol]]
    )
  }
  else{
    sampleTable <- data.frame(
      sampleName = sampleTable[[idCol]],
      files = file.path(filePrefix, sampleTable[[fileCol]]),
      condition = sampleTable[[condCol]],
      batch = sampleTable[[batchCol]]
    )
  }
  
  ## Perform analysis
  dge <- edgeR::readDGE(sampleTable)
  ## remove last rows to avoid problems with edger as they are meta tags, not genes
  dge$counts <- head(dge$count, -5)
  ## Add in gene names to use later on
  merged_dge_genes <- merge(
    dge$counts,
    gtf,
    by.x = "row.names",
    by.y = "gene_id"
  )
  ## Column name gets lost in merge to become "Row.names"
  colnames(merged_dge_genes)[1] <- "gene_id"
  dge$genes <- merged_dge_genes[,IMPORTANT_TAGS]
  ## Remove variable to save on memory
  remove(merged_dge_genes)
  dge <- edgeR::calcNormFactors(dge)
  if(class(batchCol) == "NULL"){
    design <- model.matrix(
      ~ droplevels(sampleTable[["condition"]])
    )
  } else {
    design <- model.matrix(
      ~ sampleTable[["condition"]] + sampleTable[["batch"]]
    )
  }
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmFit(dge, design)
  lrt <- edgeR::glmLRT(fit, coef=2)
  
  print(summary(limma::decideTests(lrt, p.value = 0.001)))
  
  if (class(saveName) != "NULL") {
    final_table <- as.data.frame(edgeR::topTags(lrt, n = nrow(dge$counts)))
    dir.create(
      dirname(saveName),
      recursive = TRUE,
      showWarnings = FALSE
    )
    write.csv(
      final_table,
      file=paste(saveName, ".csv", sep="")
    )
    write.csv(
      edgeR::topTags(lrt, n = 5000),
      file=paste(saveName, "_top5k.csv", sep="")
    ) 
  }
  
  return(lrt)
}


#' Get a subset of table rows
#'
#' kwanlibr_subset(sampleTable, col, ...) selects a subset of rows from 
#' sampleTable where the value of col is in ... 
#' @param sampleTable table we want to subset
#' @param col subset sampleTable based on value of col
#' @return A dataframe with the subset of selected rows
#' @keywords subset
#' @export
#' @examples
#' kwanlibr_subset(smc3, "Pool.Name", "Smc")

kwanlibr_subset <- function(sampleTable, col, ...) {
  return(
    sampleTable[grepl(paste(c(...), collapse="|"),
    sampleTable[[col]], perl=TRUE),]
  )
}


#' Label samples as ko (knockout) or CON (control)
#'
#' kwanlibr_label_con(sampleTable, col, pattern) labels a sample CON, if the
#' value in col matches pattern
#' @param sampleTable table we want to subset
#' @param col subset sampleTable based on value of col
#' @param pattern regex of which samples are CON
#' @return A copy of sampleTable with new column 'Condition' where the value is either ko or CON
#' @keywords label
#' @export
#' @examples
#' kwanlibr_label_con(smc3, "Condition", "ctrl")

kwanlibr_label_con <- function(sampleTable, col, pattern) {
  # Note: CON is purposefully all uppercase in order to ensure that it is before
  #  ko alphabetically.Otherwise, fold changes will be flipped with
  #  controls being treated as ko.
  # The levels argument specifies the order of the factors, but 
  #   manipulations outside of this function may revert the table back to
  #   default ordering, hence the uppercasing of CON.
  sampleTable$condition <- factor(
    ifelse(
      grepl(pattern, sampleTable[[col]]),
      "CON",
      "ko"
    ),
    levels = c("CON", "ko")
  )
  return(sampleTable)
}


#' Create a volcano plot
#'
#' kwanlibr_make_volcano(lrt, figure_title, filename, figure_dir, fdr, xdiff,
#' ymax, intersect, intersect_only, label_genes) creates a volcano plot from the 
#' data in lrt with title volcano_figure_title and saves it in a pdf and png 
#' under the folder figure_dir 
#' @param lrt a dataframe of the egdger data to plot
#' @param figure_title name of the plot
#' @param filename name the plot is saved under
#' @param figure_dir directory to save the plots under
#' @param fdr false discovery rate threshold. Default is 0.01
#' @param xdiff difference in the x-axis. Default is 5
#' @param ymax upper bound of y-axis. Default is 40
#' @param intersect Default is NUMM
#' @param intersect_only Default is "magenta"
#' @param label_genes Default is NULL
#' @keywords volcano plot
#' @export
#' @examples
#' kwanlibr_make_volcano(smc3, figure_title="Smc3", filename="edger_smc3",
#' figure_dir=figure_dir)

kwanlibr_make_volcano <- function(
  lrt,
  figure_title,
  filename,
  figure_dir,
  fdr=0.01,
  xdiff=5,
  ymax=40,
  intersect=NULL,
  intersect_only="magenta",
  label_genes=NULL
) {
  dir.create(figure_dir)

  # Desperately needs refactoring
  volcano_df <- edgeR::topTags(lrt, n = nrow(lrt$table))$table
  volcano_df$negLogPval <- -log10(volcano_df$PValue)
  cko_DEGs <- volcano_df$FDR < fdr
  up_DEGs <- volcano_df$logFC > 0
  
  if (intersect_only == FALSE) {
    volcano_df$color <- ifelse(
      cko_DEGs,
      ifelse(up_DEGs, "red", "blue"),
      "gray50"
    )
  }

  if (!is.null(intersect)) {
    intersect_cko_DEGs <- (volcano_df$gene_name %in% intersect) & (cko_DEGs)
    
    if (intersect_only == FALSE) {
      volcano_df$alpha <- ifelse(
        intersect_cko_DEGs,
        1,
        ifelse(cko_DEGs, deg_only_alpha, 0.1)
      )
    } else {
      volcano_df$color <- ifelse(intersect_cko_DEGs, intersect_only, "gray50")
      volcano_df$alpha <- ifelse(intersect_cko_DEGs, 1, 0.1)
    }
    
  } else {
    volcano_df$color <- ifelse(
      cko_DEGs,
      ifelse(up_DEGs, "red", "blue"),
      "gray50"
    )
    volcano_df$alpha <- ifelse(cko_DEGs, 1, 0.1)
  }
  

  # Create ceilings
  volcano_df$logFC <- ifelse(
    (volcano_df$logFC > xdiff) | (volcano_df$logFC < -xdiff),
    sign(volcano_df$logFC) * xdiff,
    volcano_df$logFC
  )
  volcano_df$negLogPval <- ifelse(
    volcano_df$negLogPval > ymax,
    ymax,
    volcano_df$negLogPval
  )
  
  # Plot significant points on top of insignificant points
  gray <- subset(volcano_df, color == "gray50")
  notgray <- subset(volcano_df, color != "gray50")
  
  p <- ggplot2::ggplot(volcano_df, aes(logFC, negLogPval)) +
    ggplot2::geom_point(data=gray, col=gray$color, alpha=gray$alpha) +
    ggplot2::geom_point(data=notgray, col=notgray$color, alpha=notgray$alpha) +
    ggplot2::xlim(-xdiff, xdiff) +
    ggplot2::ylim(0, ymax) +
    ggplot2::ggtitle(bquote(italic(.(figure_title))~"cKO")) +
    ggplot2::xlab(expression(log[2]*"FoldChange")) +
    ggplot2::ylab(expression(-log[10]*"("*italic("P")*"-value)")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 45),
          axis.title = ggplot2::element_text(size = 30),
          axis.text = ggplot2::element_text(size=23))

  if (!is.null(label_genes)) {
    label_df <- volcano_df[cko_DEGs | (volcano_df$gene_name %in% label_genes),]
    label_genes <- label_df$gene_name %in% label_genes
    label_df$gene_name <- ifelse(
      label_genes,
      paste0(label_df$gene_name),
      ""
    )
    p <- p +
       ggplot2::geom_label_repel(
        data = label_df,
        aes(logFC, negLogPval, label = gene_name),
        max.overlaps = Inf,
        box.padding = 0.5,
        # force = 4,
        size = 5,
        min.segment.length = 0,
        fill = "white",
        fontface = "italic"
        # parse = TRUE,
      ) + 
       ggplot2::geom_point(data=label_df[label_genes,], col="cyan")
  }
  
  ggplot2::ggsave(
    paste0("volcano_", filename, ".pdf"),
    path = figure_dir,
    plot=p,
    width = 6, height = 6
  )
   ggplot2::ggsave(
    filename = paste0("volcano_", filename, ".png"),
    path = figure_dir,
    plot=p,
    width = 6, height = 6
  )
  
  return(p)
}


#' Old volcano plot
#'
#' kwanlibr_make_volcano_adg(lrt, figure_title, filename, figure_dir, fdr, xdiff,
#' ymax, intersect, intersect_only) creates a volcano plot from the 
#' data in lrt with title volcano_figure_title and saves it in a pdf and png 
#' under the folder figure_dir. THIS PLOT HAS OLD POINT SIZE BEHAVIOUR.
#' @param lrt a dataframe of the egdger data to plot
#' @param figure_title name of the plot
#' @param filename name the plot is saved under
#' @param figure_dir directory to save the plots under
#' @param fdr false discovery rate threshold. Default is 0.01
#' @param xdiff difference in the x-axis. Default is 5
#' @param ymax upper bound of y-axis. Default is 40
#' @param intersect Default is NUMM
#' @param intersect_only Default is "magenta"
#' @keywords volcano plot
#' @export
#' @examples
#' kwanlibr_make_volcano_adg(smc3, figure_title="Smc3", filename="edger_smc3",
#' figure_dir=figure_dir)

# Delete this to restore old point size behavior
kwanlibr_make_volcano_adg <- function(
  lrt,
  figure_title,
  filename,
  figure_dir,
  fdr=0.01,
  xdiff=5,
  ymax=40,
  intersect=NULL,
  intersect_only="magenta"
) {
  volcano_df <- edgeR::topTags(lrt, n = nrow(lrt$table))$table
  volcano_df$negLogPval <- -log10(volcano_df$PValue)
  cko_DEGs <- volcano_df$FDR < fdr
  up_DEGs <- volcano_df$logFC > 0
  
  if (intersect_only == FALSE) {
    volcano_df$color <- ifelse(
      cko_DEGs,
      ifelse(up_DEGs, "red", "blue"),
      "gray50"
    )
  }
  
  if (!is.null(intersect)) {
    intersect_cko_DEGs <- (volcano_df$gene_name %in% intersect) & (cko_DEGs)
    
    if (intersect_only == FALSE) {
      volcano_df$alpha <- ifelse(
        intersect_cko_DEGs,
        1,
        ifelse(cko_DEGs, deg_only_alpha, 0.1)
      )
    } else {
      volcano_df$color <- ifelse(intersect_cko_DEGs, intersect_only, "gray50")
      volcano_df$alpha <- ifelse(intersect_cko_DEGs, 1, 0.1)
    }
    
  } else {
    volcano_df$alpha <- ifelse(cko_DEGs, 1, 0.1)
  }
  
  # Create ceiling
  volcano_df$logFC <- ifelse(
    (volcano_df$logFC > xdiff) | (volcano_df$logFC < -xdiff),
    sign(volcano_df$logFC) * xdiff,
    volcano_df$logFC
  )
  volcano_df$negLogPval <- ifelse(
    volcano_df$negLogPval > ymax,
    ymax,
    volcano_df$negLogPval
  )
  
  gray <- subset(volcano_df, color == "gray50")
  notgray <- subset(volcano_df, color != "gray50")
  
  p <- ggplot(volcano_df, aes(logFC, negLogPval)) +
    geom_point(data=gray, col=gray$color, alpha=gray$alpha) +
    geom_point(data=notgray, col=notgray$color, alpha=notgray$alpha, size = 2.5) +
    xlim(-xdiff, xdiff) +
    ylim(0, ymax) +
    ggtitle(bquote(italic(.(figure_title))~"cKO")) +
    xlab(expression(log[2]*"FoldChange")) +
    ylab(expression(-log[10]*"("*italic("P")*"-value)")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 45),
          axis.title = element_text(size = 30),
          axis.text = element_text(size=23))
  # geom_text_repel(
  #   data = volcano_df[cko_DEGs,],
  #   aes(label = symbol),
  #   #size = 5,
  #   #min.segment.length = unit(0, 'lines')
  # )
  
  ggsave(
    paste0(figure_dir, "/volcano_", filename, ".pdf"),
    plot=p,
    width = 6, height = 6
  )
  ggsave(
    paste0(figure_dir, "/volcano_", filename, ".png"),
    plot=p,
    width = 6, height = 6
  )
  
  return(p)
}


#test function for documentation generation

#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

cat_function <- function(love=TRUE){
  if(love==TRUE){
    print("I love cats!")
  }
  else {
    print("I am not a cool person.")
  }
}

