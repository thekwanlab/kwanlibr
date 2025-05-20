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
draw_PCA(mtcars, label = as.factor(mtcars$am), color = c("coral", "cyan"))
draw_PCA(mtcars, label = as.factor(mtcars$vs), color = c("coral", "cyan"), legend = FALSE)
draw_PCA(mtcars, label = mtcars$carb, legend = FALSE)
draw_PCA(mtcars, label = mtcars$gear)

# Error-Handling Tests
draw_PCA(iris, label=iris$Species, color=c("green", "blue", "yellow"))
draw_PCA(iris[-5], label=iris$Species[-1], color=c("green", "blue", "yellow"))
draw_PCA(iris[-5], label=append(iris$Species, "2"), color=c("green", "blue", "yellow"))
draw_PCA(iris[-5], label=iris$Species, color=c("green", "blue"))
draw_PCA(iris[-5], label=iris$Species, color=c("green", "blue", "coral", "lightblue3"))

mtcars_data <- mtcars
mtcars_data$car_type <- rownames(mtcars)
draw_PCA(mtcars_data, label = as.factor(mtcars$vs), color = c("coral", "cyan"))
draw_PCA(mtcars, label = as.factor(mtcars$vs[-2]), color = c("coral", "cyan"))
draw_PCA(mtcars, label = as.factor(append(mtcars$vs, 1)), color = c("coral", "cyan"))
draw_PCA(mtcars, label = mtcars$carb, color = c("coral"))
draw_PCA(mtcars, label = mtcars$carb, color = c("coral", "cyan", "darkseagreen2"))

# =========================
# test draw_volcano_general()
# =========================

# Generate mock data
set.seed(123)
num_genes <- 1000
FDR_THRESHOLD <- 0.05
genes <- paste0("Gene", 1:num_genes)
expression_values <- rnorm(num_genes, mean = 0, sd = 2)
p_values <- runif(num_genes, min = 1e-10, max = 6)

data <- data.frame(Gene = genes, Fold = expression_values, FDR = p_values) %>%
  mutate(negLogFDR = -log10(FDR))

data_with_criteria <- data %>%
  mutate(FDR = ifelse(FDR <= FDR_THRESHOLD, paste('<=',FDR_THRESHOLD), paste('>',FDR_THRESHOLD)))

draw_volcano_general(data, colors = "grey", ymax = 2)
draw_volcano_general(data_with_criteria,
                     criteria = "FDR",
                     colors = setNames(c("aquamarine2", "grey"), c('<= 0.05', '> 0.05')),
                     ymax = 2)
