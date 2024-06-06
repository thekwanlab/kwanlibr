config = config::get(file = 'informal_tests/test_suite_config.yaml')
devtools::load_all(config$paths$package)

# =========================
# test ggsave_vector_raster
# =========================
df = data.frame(a=rnorm(100), b=runif(100))
p_hist = ggplot(df, aes(x=a)) +
  geom_histogram()

# pass arguments via ..., save plot from object
ggsave_vector_raster(file.path(config$paths$test_results_dump,
                               'test_from_object_ggsave_vector_raster'),
                     plot = p_hist)

# save from most recently displayed plot, change dims
p_hist
ggsave_vector_raster(file.path(config$paths$test_results_dump,
                               'test_from_display_ggsave_vector_raster'),
                     width=8, height=2)
