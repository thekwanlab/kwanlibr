config = config::get(file = 'informal_tests/test_suite_config.yaml')
devtools::document()
devtools::load_all(config$paths$package)

# ==============
# Test perform_diffbind
# ==============
samplecsv = config$paths$binding_sites_samples_example
perform_diffbind(samplesheetpath = samplecsv,
                 minMembers = 2,
                 reference="HET",
                 FDRthreshold = 0.05,
                 filesuffix="h3k27ac",
                 normalization = DBA_NORM_RLE,
                 filePath=config$paths$test_results_dump)

# ==============
# Test get_DBsites
# ==============
dbapath = paste0(config$paths$test_results_dump, "dba_obj_h3k27ac.RDS")
get_DBsites(dbapath,filepath = paste0(config$paths$test_results_dump, "diffbind_files"))

#===============
# Test make_density_plot
#===============
dbapath = paste0(config$paths$test_results_dump, "dba_obj_h3k27ac.RDS")
make_density_plot(dbapath,
                  DENSITY_color='aquamarine4',
                  figure_title="Distribution of Differential Binding Regions",
                  figure_dir = paste0(config$paths$test_results_dump, "diffbind_figures"),
                  filename="df_fold_density_h3k27ac")

#===============
# Test make_PCA_plot_diffbind
#===============
dbapath = paste0(config$paths$test_results_dump, "dba_obj_h3k27ac.RDS")
make_PCA_plot_diffbind(dbapath,
                       figure_title_nocontrast="All consensus Peaks: h3k27ac HET vs KO",
                       figure_title_contrast="'Differential Binding Peaks: h3k27ac HET vs KO'",
                       figure_dir=paste0(config$paths$test_results_dump, "diffbind_figures"),
                       filename="PCA_consensus_vs_diffbind_peaks_kwanlib_h3k27ac")

#===============
# Test make_volcano_sites
#===============
dbapath = paste0(config$paths$test_results_dump, "dba_obj_h3k27ac.RDS")
make_volcano_sites(dbapath,
                   sigsitespath=paste0(config$paths$test_results_dump, "diffbind_files/", "sig_sites.bed"),
                   filepath=paste0(config$paths$test_results_dump, "diffbind_files"),
                   xdiff=1,
                   ymax=10)

#===============
# Test make_volcano_plot_diffbind
#===============
make_volcano_plot_diffbind(volcanopath = paste0(config$paths$test_results_dump, "diffbind_files/", "volcano_sites.bed"),
                           sigsitespath = paste0(config$paths$test_results_dump, "diffbind_files/", "sig_sites.bed"),
                           figure_title="h3k27ac cKO vs cHET",
                           figure_dir=paste0(config$paths$test_results_dump, "diffbind_figures"),
                           filename="db_volcano",
                           xdiff = 1,
                           ymax = 10)
#================
# Test merge_sites_with_exp
#================
merged_data <- merge_sites_with_exp(allDBsitepath = "/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/allDB_sites.bed",
                     joindatapath = "/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/volcano_sites.bed",
                     RNAfilepath = "/nfs/turbo/umms-kykwan/projects/mll/bulk_rna_seq/tables/mll_edgeR.csv",
                     filepath = paste0(config$paths$test_results_dump, "diffbind_files"),
                     BULK_FDR_CUTOFF = 0.05)

#=================
# Test make_volcano_plot_from_merged
#=================
make_volcano_plot_from_merged(allDBsitepath = "/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/allDB_sites.bed",
                              RNAfilepath="/nfs/turbo/umms-kykwan/projects/mll/bulk_rna_seq/tables/mll_edgeR.csv",
                              volcanopath="/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/volcano_sites.bed",
                              filepath=paste0(config$paths$test_results_dump, "diffbind_files"),
                              BULK_FDR_CUTOFF=0.05,
                              figure_title='cKO vs cHet Differential Binding sites',
                              filename="db_volcano_bulk_color_binary",
                              figure_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                              color=c('grey','steelblue', 'tomato'),
                              xdiff = 2,
                              ymax = 17.5,
                              size=1,
                              alpha=1,
                              width = 8,
                              height = 6)

#=================
# Test make_scatter_plot_from_merged
#=================
make_scatter_plot_from_merged(allDBsitepath = "/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/allDB_sites.bed",
                              RNAfilepath="/nfs/turbo/umms-kykwan/projects/mll/bulk_rna_seq/tables/mll_edgeR.csv",
                              sigsitespath="/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/sig_sites.bed",
                              filepath=paste0(config$paths$test_results_dump, "diffbind_files"),
                              BULK_FDR_CUTOFF=1,
                              figure_title ="DB sites log2FC vs Nearby Gene RNAseq log2FC",
                              filename="db_bulk_scatter",
                              figure_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                              width = 8,
                              height = 6,
                              regression = TRUE)

make_scatter_plot_from_merged(allDBsitepath = "/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/allDB_sites.bed",
                              RNAfilepath="/nfs/turbo/umms-kykwan/projects/mll/bulk_rna_seq/tables/mll_edgeR.csv",
                              sigsitespath="/nfs/turbo/umms-kykwan/projects/mll/h3k4_me1_cutntag/diffbind_files/sig_sites.bed",
                              filepath=paste0(config$paths$test_results_dump, "diffbind_files"),
                              BULK_FDR_CUTOFF=1,
                              figure_title ="DB sites log2FC vs Nearby Gene RNAseq log2FC",
                              filename="db_bulk_scatter",
                              figure_path=paste0(config$paths$test_results_dump, "diffbind_figures"),
                              width = 8,
                              height = 6,
                              regression = FALSE)
