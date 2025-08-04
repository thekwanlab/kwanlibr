config = config::get(file = 'informal_tests/test_suite_config.yaml')
devtools::document()
devtools::load_all(config$paths$package)

# ==============
# Test perform_diffbind
# ==============
samplecsv <- config$paths$binding_sites_samples_example
dba.obj <- perform_diffbind(sample_sheet_path = samplecsv,
                            control_level ="HET",
                            min_members = 2,
                            normalization = DBA_NORM_RLE)

# ==============
# Test save_diffbind_object
# ==============
save_diffbind_object(dba_object=dba.obj,
                     file_suffix = "h3k27ac",
                     save_directory = config$paths$test_results_dump)

# ==============
# Test get_diffbind_sites
# ==============
sig.sites <- get_diffbind_sites(dba_object=dba.obj, verbose = FALSE)
sig.sites <- get_diffbind_sites(dba_object=dba.obj, fdr_threshold = 0.0)
all.sites <- get_diffbind_sites(dba_object=dba.obj, fdr_threshold = 1)

# ==============
# Test save_diffbind_sites
# ==============
save_diffbind_sites(dba.obj,
                    save_directory = file.path(config$paths$test_results_dump, "diffbind_files"))
save_diffbind_sites(dba.obj,
                    save_directory = file.path(config$paths$test_results_dump, "diffbind_files"),
                    fdr_threshold = 0.0)

#===============
# Test make_diffbind_density_plot
#===============
make_diffbind_density_plot(dba.obj,
                           figure_title="Distribution of Differential Binding Regions",
                           save_name="df_fold_density_h3k27ac",
                           save_directory = file.path(config$paths$test_results_dump, "diffbind_figures"))

make_diffbind_density_plot(dba.obj,
                           figure_title="Distribution of Differential Binding Regions",
                           save_name="df_fold_density_h3k27ac",
                           save_directory = file.path(config$paths$test_results_dump, "diffbind_figures"),
                           xdiff = 1)

make_diffbind_density_plot(dba.obj,
                           figure_title="Distribution of Differential Binding Regions",
                           save_name="df_fold_density_h3k27ac",
                           save_directory = file.path(config$paths$test_results_dump, "diffbind_figures"),
                           fdr_threshold = 0.0)

#===============
# Test make_diffbind_PCA_plot
#===============
make_diffbind_PCA_plot(dba.obj,
                       figure_title_nocontrast="All consensus Peaks: h3k27ac HET vs KO",
                       save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
                       save_name="PCA_consensus_vs_diffbind_peaks_kwanlib_h3k27ac")

make_diffbind_PCA_plot(dba.obj,
                       figure_title_nocontrast="All consensus Peaks: h3k27ac",
                       figure_title_contrast="'Differential Binding Peaks: h3k27ac HET vs KO'",
                       save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
                       save_name="PCA_consensus_vs_diffbind_peaks_kwanlib_h3k27ac")

#===============
# Test get_diffbind_volcano_data
#===============
volcano.sites <- get_diffbind_volcano_data(dba.obj, xdiff=1, ymax=10)
volcano.sites <- get_diffbind_volcano_data(dba.obj, xdiff=1, ymax=10, fdr_threshold = 0.0)

#===============
# Test make_diffbind_volcano_plot
#===============
make_diffbind_volcano_plot(dba.obj,
                           figure_title="h3k27ac cKO vs cHET",
                           save_name="db_volcano_h3k27ac",
                           save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
                           point_size = 2,
                           point_alpha = 0.8)

make_diffbind_volcano_plot(dba.obj,
                           figure_title="h3k27ac cKO vs cHET",
                           save_name="db_volcano_h3k27ac",
                           save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
                           fdr_threshold = 0.0)
#================
# Test merge_diffbind_with_DE
#================
options(scipen = 15)
merged_data <- merge_diffbind_with_DE(dba.obj,
                                      DE_file_path = config$paths$binding_sites_gene_expression_example,
                                      tss_file_path = config$paths$binding_sites_transcription_start_site_example,
                                      save_directory = file.path(config$paths$test_results_dump, "diffbind_figures"))
#=================
# Test make_diffbind_volcano_plot_from_merged
#=================
make_diffbind_volcano_plot_from_merged(
  merged_df = merged_data,
  figure_title=paste0('cKO vs cHet Differential Binding sites', '\n[ color by nearest TSS bulkRNAseq logFC, filtered by bulkRNAseq FDR < 0.05 ]'),
  save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
  save_name="db_volcano_bulk_color_binary_h3k4me1",
  xdiff = 2,
  ymax = 17.5,
  point_size = 2,
  point_alpha = 1,
  width = 8,
  height = 6)

make_diffbind_volcano_plot_from_merged(
  merged_df = merged_data,
  figure_title=paste0('cKO vs cHet Differential Binding sites', '\n[ color by nearest TSS bulkRNAseq logFC, filtered by bulkRNAseq FDR < 0.05 ]'),
  save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
  save_name="db_volcano_bulk_color_binary_h3k4me1",
  xdiff = 2,
  ymax = 17.5,
  point_size = 2,
  point_alpha = 1,
  width = 8,
  height = 6,
  DE_fdr_cutoff = 0.0,
  DB_fdr_cutoff = 0.0)


# ==============
# Test add_DB_DE_regression
# ==============
sig_merged_data <- merged_data %>%
  na.omit() %>%
  filter(DB.FDR < 0.05) %>%
  filter(DE.FDR < 1)

DB_DE_regression <- DB_DE_regression_model(merged_df = sig_merged_data)

#=================
# Test make_scatter_plot_from_merged
#=================
make_scatter_plot_from_merged(merged_df = merged_data,
                              figure_title ="DB sites log2FC vs Nearby Gene RNAseq log2FC",
                              save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
                              save_name="db_bulk_scatter_h3k4me1",
                              width = 8,
                              height = 6,
                              regression = TRUE)

make_scatter_plot_from_merged(merged_df = merged_data,
                              figure_title ="DB sites log2FC vs Nearby Gene RNAseq log2FC",
                              save_directory=file.path(config$paths$test_results_dump, "diffbind_figures"),
                              save_name="db_bulk_scatter_no_regression_line_h3k4me1",
                              width = 8,
                              height = 6,
                              regression = FALSE)
