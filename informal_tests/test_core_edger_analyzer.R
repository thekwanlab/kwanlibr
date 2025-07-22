config = config::get(file = 'informal_tests/test_suite_config.yaml')
devtools::document()
devtools::load_all(config$paths$package)

# ==============
# Test get_gtf
# ==============

gtf = kwanlibr::get_gtf(config$paths$gtf)

samples = read.table(config$paths$gene_expression_samples_example, header=TRUE, sep=',')
Filename = c("1395-ML-1_AGTCAA-TAGATC_S156.tsv",
             "1395-ML-1_AGTTCC-TAGATC_S157.tsv",
             "1395-ML-1_ATGTCA-TAGATC_S158.tsv",
             "1395-ML-1_CCGTCC-TAGATC_S159.tsv",
             "1395-ML-1_GTAGAG-TAGATC_S160.tsv",
             "1395-ML-1_GTCCGC-TAGATC_S161.tsv",
             "1395-ML-1_GTGAAA-TAGATC_S162.tsv",
             "1395-ML-1_GTGGCC-TAGATC_S163.tsv")


# ===================
# Test subset_samples
# ===================
exp_samples = kwanlibr::subset_samples(samples, "Pool.Name", "Smc")

exp_samples$Filename = Filename

# ===================
# Test label_control_samples
# ===================
exp_samples = kwanlibr::label_control_samples(exp_samples, 'Condition', 'ctrl')

# ===================
# Test perform_edger
# ===================

# GTF object
lrt = kwanlibr::perform_edger(
  exp_samples,
  filePrefix = config$paths$gene_expression_gene_counts_dir,
  idCol = 'Pool.Name',
  gtf=gtf,
  saveName = paste0(config$paths$test_results_dump, '/', 'edger_lrt')
)

# GTF file path
lrt = kwanlibr::perform_edger(
  exp_samples,
  filePrefix = config$paths$gene_counts_dir,
  idCol = 'Pool.Name',
  gtf=config$paths$gtf,
  saveName = paste0(config$paths$test_results_dump, '/', 'edger_lrt')
)

# no GTF
lrt = kwanlibr::perform_edger(
  exp_samples,
  filePrefix = config$paths$gene_counts_dir,
  idCol = 'Pool.Name',
  saveName = paste0(config$paths$test_results_dump, '/', 'edger_lrt')
)

# ===================
# Test make_volcano_dataframe
# ===================

volcano_df = kwanlibr::make_volcano_dataframe(lrt)
head(volcano_df)

# ===================
# Test draw_volcano
# ===================

p = kwanlibr::draw_volcano(lrt, figure_title = 'Volcano Time')
p

# ===================
# Test label_volcano_genes
# ===================

p_labeled = label_volcano_plot_genes(p, lrt, c('Cars','Cars2','Cars3','Cars4'))
p_labeled

# ===================
# Test make_volcano
# ===================

kwanlibr::make_volcano(
  lrt = lrt,
  figure_title="Gene X cKO vs Control",
  filename = 'test_make_volcano_labels',
  figure_dir = config$paths$test_results_dump,
  label_genes = c('Cars', 'Cars2', 'Cux1', 'Cux2')
)

# ===================
# Test make_MA_plot
# ===================

kwanlibr::make_MA_plot(
  lrt =lrt,
  figure_title="logCPM vs logFC",
  filename = "test_make_ma_plot",
  figure_dir = config$paths$test_results_dump
)

kwanlibr::make_MA_plot(
  lrt =lrt,
  figure_title="logCPM vs logFC",
  filename = "test_make_ma_plot",
  figure_dir = config$paths$test_results_dump,
  xmax = 10,
  smooth_line = TRUE
)

# ==============
# Test make_PCA_plot
# ==============

make_PCA_plot(lrt, color = c("blue", "red"), figure_title = "Principal Component Analysis",
              figure_dir = config$paths$test_results_dump, filename = "plot_assigned_color_and_legend")
make_PCA_plot(lrt, legend = FALSE, figure_title = "Principal Component Analysis",
              figure_dir = config$paths$test_results_dump, filename = "plot_default_color_no_legend")
make_PCA_plot(lrt, figure_title = "Principal Component Analysis",
              figure_dir = config$paths$test_results_dump, filename = "plot_default_color_and_legend")
make_PCA_plot(lrt, label = TRUE, figure_title = "Principal Component Analysis",
              figure_dir = config$paths$test_results_dump, filename = "plot_labels_default_color_and_legend")
make_PCA_plot(lrt, figure_title = "Principal Component Analysis",
              figure_dir = config$paths$test_results_dump, filename = "plot_no_labels_default_color_and_legend")
