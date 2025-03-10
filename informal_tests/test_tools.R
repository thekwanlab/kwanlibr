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

# =========================
# test clamp
# =========================

x = 1:20
clamp(x)
clamp(x, lower=4, upper=15)
clamp(x, lower=10, upper=10)
clamp(x, lower=10, upper=2) # throw error

# =========================
# test draw_PCA
# =========================

draw_PCA(iris[-5], label = iris$Species, color = c("green", "blue", "yellow"))
draw_PCA(iris[-5], label = iris$Species)
draw_PCA(iris[-5], label = iris$Species, legend = FALSE)
draw_PCA(iris[-5], legend=FALSE)
draw_PCA(iris[-5])
draw_PCA(mtcars, label = mtcars$vs, color = c("coral", "cyan"))
draw_PCA(mtcars, label = mtcars$vs, color = c("coral", "cyan"), legend = FALSE)
draw_PCA(mtcars, label = mtcars$carb, legend = FALSE)
draw_PCA(mtcars, label = mtcars$gear)

library(palmerpenguins)
df <- na.omit(penguins)
draw_PCA(df[-c(1,2,7)], label = df$species, color = c("cyan", "coral", "orange"))
draw_PCA(df[-c(1,2,7)], label = df$sex, color = c("cyan", "violet"))
draw_PCA(df[-c(1,2,7)], label = df$island, color = c("cyan", "violet", "lightblue4"))
draw_PCA(df[-c(1,2,7)], label = df$year)
draw_PCA(df[-c(1,2,7)], label = df$year, legend = FALSE)
