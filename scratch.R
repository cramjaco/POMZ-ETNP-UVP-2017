moo <- calc_small_psd(unbinned)
debug(calc_small_psd)

iris %>% ggplot(aes(x = Sepal.Width, y = Sepal.Length, fill = Species, shape = Species)) + geom_point(size = 2, alpha = 0.75) +
  scale_shape_manual(values = 21:25) +
  scale_fill_manual(values = c("black", "red", "blue"), breaks = c("setosa", "versicolor", "virginica"))
