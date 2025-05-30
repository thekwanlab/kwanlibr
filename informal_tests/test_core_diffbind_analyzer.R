config = config::get(file = 'informal_tests/test_suite_config.yaml')
devtools::document()
devtools::load_all(config$paths$package)

# ==============
# Test perform_diffbind
# ==============
samplecsv <- config$paths$binding_sites_samples_example
dba.obj <- perform_diffbind(sample_sheet_path = samplecsv,
                            min_members = 2,
                            control_level ="HET",
                            fdr_threshold = 0.05,
                            normalization = DBA_NORM_RLE)

# ==============
# Test save_diffbind_object
# ==============
save_diffbind_object(dba_object=dba.obj,
                     file_suffix = "h3k27ac",
                     dbaobj_save_path = "informal_tests/test_results")

# ==============
# Test get_DBsites
# ==============
sig.sites <- get_DBsites(dba.obj, fdr_threshold = 0.05, contrast_number = 1)

all.sites <- get_DBsites(dba.obj, fdr_threshold = 1, contrast_number = 1)

# ==============
# Test save_sites
# ==============
save_sites(dba.obj,
           save_path = paste0(config$paths$test_results_dump, "diffbind_files"),
           fdr_threshold = 0.05,
           contrast_number = 1)

#===============
# Test make_density_plot
#===============
make_density_plot(dba.obj,
                  fdr_threshold = 1,
                  color='aquamarine4',
                  figure_title="Distribution of Differential Binding Regions",
                  figure_save_path = paste0(config$paths$test_results_dump, "diffbind_figures"),
                  x_lim = c(-1, 1),
                  file_name="df_fold_density_h3k27ac")

#===============
# Test make_PCA_plot_diffbind
#===============
make_PCA_plot_diffbind(dba.obj,
                       figure_title_nocontrast="All consensus Peaks: h3k27ac HET vs KO",
                       figure_title_contrast="'Differential Binding Peaks: h3k27ac HET vs KO'",
                       figure_save_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                       file_name="PCA_consensus_vs_diffbind_peaks_kwanlib_h3k27ac")

#===============
# Test make_volcano_sites
#===============
volcano.sites <- make_volcano_sites(dba.obj, xdiff=1, ymax=10)

#===============
# Test make_volcano_plot_diffbind
#===============
make_volcano_plot_diffbind(dba.obj,
                           figure_title="h3k27ac cKO vs cHET",
                           figure_save_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                           file_name="db_volcano_h3k27ac",
                           point_size = 2,
                           point_alpha = 0.8,
                           contrast_number = 1,
                           xdiff = 2,
                           ymax = 10)
#================
# Test merge_sites_with_exp
#================
dba.obj <- readRDS("/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/dba_obj_mll_me1.RDS")
options(scipen = 15)
merged_data <- merge_sites_with_exp(dba.obj,
                                    DE_file_path = "/nfs/turbo/umms-kykwan/projects/mll/bulk_rna_seq/tables/mll_edgeR.csv",
                                    tss_file_path = "informal_tests/test_results/diffbind_files_me1/TSS.bed",
                                    save_path = paste0(config$paths$test_results_dump, "diffbind_files_me1"))
#=================
# Test make_volcano_plot_from_merged
#=================
make_volcano_plot_from_merged(merged_df = merged_data,
                              figure_title='cKO vs cHet Differential Binding sites',
                              figure_save_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                              file_name=paste0('cKO vs cHet Differential Binding sites',
                                               '\n[ color by nearest TSS bulkRNAseq logFC, filtered by bulkRNAseq FDR < 0.05 ]'),
                              xdiff = 2,
                              ymax = 17.5,
                              point_size=2,
                              point_alpha=1,
                              width = 8,
                              height = 6)

#=================
# Test make_scatter_plot_from_merged
#=================
make_scatter_plot_from_merged(merged_df = merged_data,
                              figure_title ="DB sites log2FC vs Nearby Gene RNAseq log2FC",
                              file_name="db_bulk_scatter_h3k4me1",
                              figure_save_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                              width = 8,
                              height = 6,
                              regression = TRUE)

make_scatter_plot_from_merged(merged_df = merged_data,
                              figure_title ="DB sites log2FC vs Nearby Gene RNAseq log2FC",
                              file_name="db_bulk_scatter_no_regression_line_h3k4me1",
                              figure_save_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                              width = 8,
                              height = 6,
                              regression = FALSE)
