config = config::get(file = 'informal_tests/test_suite_config.yaml')
devtools::load_all(config$paths$package)

# ==============
# Test get_gtf
# ==============

gtf = kwanlibr::get_gtf(config$paths$gtf)

samples = read.table(config$paths$samples_example, header=TRUE, sep=',')
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
lrt = kwanlibr::perform_edger(
  exp_samples,
  filePrefix = config$paths$gene_counts_dir,
  idCol = 'Pool.Name',
  gtf=gtf,
  saveName = paste0(config$paths$test_results_dump, 'edger_lrt')
)

# ===================
# Test make_volcano
# ===================
kwanlibr::make_volcano(
  lrt = lrt,
  figure_title="Gene X cKO vs Control",
  filename = 'test_make_volcano',
  figure_dir = config$paths$test_results_dump,
  
)
