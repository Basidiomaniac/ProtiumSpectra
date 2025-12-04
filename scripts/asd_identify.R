# Script to conduct linear discriminant analysis (LDA) on spectral dataframes
# Code adapted from appetizer_course_dm20251103.R (author: Niko Darci-Maher)
# Adapted by Laurel Miller, Nov 2025
# Input: Spectral dataframe from asd_cleanup.R (nir)
# Output: LDA model capable of discriminating species at 99% accuracy

library(viridis)
library(klaR)
library(lme4)
library(randomForest)
library(caret)
library(class)

# plot mean wavelength per sample
nir_sample <- nir %>%
  group_by(sample) %>%
  summarise(
    across(starts_with("nm"), mean),
    across(-starts_with("nm"),
    ~ first(.x), .names = "{col}"))
nir_sample_long = nir_sample %>% pivot_longer(cols = starts_with('nm'),
               names_to = 'wavelength',
               names_prefix = 'nm',
               values_to = 'reflectance_mean') %>%
              mutate(wavelength = as.numeric(wavelength),
              sample = as.character(sample))

line_plot_sample = ggplot(nir_sample_long, aes(x = wavelength, y = reflectance_mean, color = species)) +
  geom_line(aes(group = sample)) 
# line_plot_sample

# means for each species, plot
nir_species = nir %>% group_by(species) %>%
  summarize_at(vars(starts_with('nm')), mean)
nir_species = nir_species %>% pivot_longer(cols = starts_with('nm'),
                                           names_to = 'wavelength',
                                           names_prefix = 'nm',
                                           values_to = 'reflectance_mean') %>%
  mutate(wavelength = as.numeric(wavelength))
line_plot_species = ggplot(nir_species, aes(x = wavelength, y = reflectance_mean, color = species)) +
  geom_line(aes(group = species)) +
  scale_color_viridis_d(option = "plasma")
# line_plot_species

line_plot_species = ggplot() +
  geom_line(data=nir_species, aes(x=wavelength, y=reflectance_mean), alpha = 0.5) +
  geom_line(data = nir_species, aes(x=wavelength, y=reflectance_mean, color=population_lineage, group = population_lineage)) +
  scale_color_viridis_d(option = "plasma")
# line_plot_species

# plot with density by species
nir_long = lengthen(nir)
line_plot_species = ggplot(nir_long, aes(x = wavelength, y = intensity, color = species)) +
  geom_line(aes(group = sample), alpha = 0.2) +
  scale_color_viridis_d(option = "plasma")
# line_plot_species

# check the normality of all the wavelengths
# NOTE: When running all test data, multimodal distributions show up. Why??
# ANSWER: Each date had unique data skews. Day-to-day calibration MATTERS.
# In future sampling, should find a way to standardize results. Standard sample?
nir_norm = nir %>% mutate_at(vars(starts_with('nm')), scale)
nir_long = nir_norm %>% dplyr::select(sample, replicate, species, date, starts_with('nm')) %>%
  pivot_longer(cols = starts_with('nm'),
               names_to = 'wavelength',
               values_to = 'reflectance')
normaldist_plot = ggplot(nir_long, aes(x = reflectance, color = wavelength)) + 
  geom_density(alpha = 0.5) + 
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ date)
# normaldist_plot

# hmm, measurements across dates seem to vary, but this is also because we 
# sampled different species on different days. Let's compare within species
# across dates.
# currently crashes
# date_lmer = lmer(intensity ~ date + (1 | species), data = nir_long)
# summary(date_lmer)
# no surprise: species and date are sig. predictors of intensity,


# separate out wavelength data, produce LDA model!
# pull the first for each species as test data, reserve the rest as training data
nir_training = nir %>% group_by(species) %>%
                  filter(sample != min(sample)) # all but first sample per species clade
nir_testing = nir %>% group_by(species) %>%
                  filter(sample == min(sample)) # ONLY first sample per species clade
lda_model = lda(x = nir_training[, wavelengths], grouping = nir_training$species)

# make the species prediction! note that we grab data by species clades, but only
# predict to species level. this avoids issues if clades are spectrally distinct
# but we aren't quite ready to predict clades outright, so we try species
species_predict = MASS:::predict.lda(lda_model, newdata = nir_testing[, wavelengths])

# how did we do on actual species accuracy?
accuracy = mean(species_predict$class == nir_testing$species)
accuracy # for all test data, the accuracy comes out at 89%, not bad

# let's average all the points for each sample and see if that helps
# note we use nir_sample, which is already averaged to sample
sample_avg_train = nir_sample %>% group_by(species) %>%
  filter(sample != (min(sample))) # all but first sample per species
sample_avg_test = nir_sample %>% group_by(species) %>%
  filter(sample == (min(sample))) # ONLY first sample per species
sample_lda_model = lda(x = sample_avg_train[, wavelengths], grouping = sample_avg_train$species)
sample_avg_predictions = MASS:::predict.lda(sample_lda_model, newdata = sample_avg_test[, wavelengths])
mean(sample_avg_predictions$class == sample_avg_test$species)
mean(sample_avg_predictions$posterior)
# cool! accuracy bumps up from 80% to 90% when averaged to sample

# Get LD1 and LD2 scores for species
lda_data <- data.frame(
  Species = nir_testing$species,
  Sample = nir_testing$sample,
  LD1 = species_predict$x[, 1],
  LD2 = species_predict$x[, 2] # Only the first two LDs are used
)

# get LD1 and LD2 for averaged data
sample_lda_data = data.frame(
  Species = sample_avg_test$species,
  Sample = sample_avg_test$sample,
  LD1 = sample_avg_predictions$x[, 1],
  LD2 = sample_avg_predictions$x[, 2]
)

# plot LD in 2D space
lda_space_plot <- ggplot(sample_lda_data, aes(x = LD1, y = LD2, color = Species)) +
  geom_point(alpha = 0.9) +
  stat_ellipse(geom = "polygon", type = "norm", aes(fill = Species), alpha = 0.1, show.legend = TRUE) +
  labs(
    title = "LDA Coverage of Protium spp.",
    subtitle = paste0(
      "LD1 explains ", round(sample_lda_model$svd[1]^2 / sum(sample_lda_model$svd^2) * 100, 1), 
      "% of variance\n", " LD2 explains ", round(sample_lda_model$svd[2]^2 / sum(sample_lda_model$svd^2) * 100, 1),
      "% of variance"
    ),
    x = "Linear Discriminant 1 (LD1)",
    y = "Linear Discriminant 2 (LD2)"
  )
# lda_space_plot # looks like the LDA covers a decent amount of variance, species are well-separated

# let's see how good ALL our LD's are
ld_variance = round(sample_lda_model$svd^2 / sum(sample_lda_model$svd^2) * 100, 1) %>% data.frame() 
colnames(ld_variance) = "explains"
rownames(ld_variance) = c("LD1", "LD2", "LD3", "LD4", "LD5", "LD6", "LD7", "LD8", "LD9")
ld_variance$ld = rownames(ld_variance)
lda_bar_plot = ggplot(ld_variance, aes(x = ld, y = explains)) +
  geom_col(fill = "#13536a", color = "#08191f") + 
  labs(
    title = "Linear Discriminant Effectiveness",
    x = "Linear Discriminant",
    y = "% of Variance Explained")
# lda_bar_plot

# Predict species (and get LD's) for ALL data points, not averaged by sample
# This is for visualization
full_data_predict = MASS:::predict.lda(sample_lda_model, newdata = nir_testing[, wavelengths])
full_lda_data = data.frame(
  Species = nir_testing$species,
  Sample = nir_testing$sample,
  LD1 = full_data_predict$x[, 1],
  LD2 = full_data_predict$x[, 2] # Only the first two LDs are used
)

# urchin plot - a function to visualize samples and their averages in LD space.
# provide your lda data as above, plus x and y limits for the plot.
urchin_plot = function(lda_data, lda_points, xmin, xmax, ymin, ymax) {
  # map boundaries. these will need to be adjusted to fit your values.
  x_min = xmin
  x_max = xmax
  y_min = ymin
  y_max = ymax
  grid = expand.grid(LD1 = seq(x_min, x_max, by = 0.1),
                     LD2 = seq(y_min, y_max, by = 0.1))  
  grid$Species = knn(
    train = lda_data[, c("LD1", "LD2")],
    test  = grid[, c("LD1", "LD2")],
    cl    = lda_data$Species,
    k = 1
  )
  centers = lda_points %>% group_by(Sample) %>% 
    summarize(LD1_cent = mean(LD1, na.rm = TRUE),
              LD2_cent = mean(LD2, na.rm = TRUE))
  full_lda_centers = lda_points %>%  left_join(centers, by = "Sample")
  tile_plot = ggplot() + 
    geom_tile(data = grid, aes(x = LD1, y = LD2, fill = Species), alpha = 0.8) +
    geom_segment(data = full_lda_centers, 
                 aes(x = LD1, y = LD2,
                     xend = LD1_cent, yend = LD2_cent),
                 alpha = 0.1, linewidth = 0.5) +
    geom_point( # plot sample averages (centers)
      data = full_lda_centers,
      aes(x = LD1_cent, y = LD2_cent, color = Species, size = 2)
    ) +
    labs(
      subtitle = paste0(
        "LD1 explains ", round(lda_model$svd[1]^2 / sum(lda_model$svd^2) * 100, 1), "% of variance\n",
        "LD2 explains ", round(lda_model$svd[2]^2 / sum(lda_model$svd^2) * 100, 1), "% of variance"
      ),
      x = "Linear Discriminant 1",
      y = "Linear Discriminant 2"
    ) +
    scale_color_viridis_d(begin = 0, end = 0.9, option = "A") + # colors points
    scale_fill_viridis_d(option = "C",) + # colors background grid
    geom_point(data = full_lda_centers, # this adds the full 12 points per sample
      aes(x = LD1, y = LD2, color = Species, size = 1),
      alpha = 0.1)
  tile_plot
}

# visualizing sample and their subpoints on LD space (10 samples, 12 subpoints each)
urchin_plot(sample_lda_data, full_lda_data, -50, 25, -50, 25)

# urchin plot of ALL points in the dataset
nir_predict = MASS:::predict.lda(sample_lda_model, newdata = nir[, wavelengths])
mean(nir_predict$class == nir$species) # 93% accuracy on full data
nir_lda_data = data.frame(
  Species = nir$species,
  Sample = nir$sample,
  LD1 = nir_predict$x[, 1],
  LD2 = nir_predict$x[, 2] # Only the first two LDs are used
)
# urchin_plot(sample_lda_data, nir_lda_data, -50, 25, -50, 25)
# looks messy, even though classification is 93% accurate

# let's look at loading of each wavelengths on the model LD's
lda_loading = sample_lda_model$scaling %>%
  as.data.frame() %>%
  mutate(wavelength = wavelengths) %>%
  pivot_longer(cols = starts_with("LD"), 
               names_to = "LD", 
               values_to = "loading")
loading_plot = ggplot(lda_loading, aes(x = wavelength, y = loading, group = LD, color = LD)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  labs(x = "Wavelength (nm)",
       y = "LDA loading",
       title = "LDA Loadings by Wavelength") +
  scale_color_viridis_d(option = "plasma") +
  theme(legend.position = "none") +
  facet_wrap(~ LD, ncol = 1) +
  theme_void()
loading_plot

# Summary: LDA works really well for species ID! Just collapsing our samples
# by avg. for each wavelength brings us to about 90% accuracy.
#
# NEXT: open asd_analysis.R to produce statistics on wavelength variation and
# potential connections to the metabolome, environment, etc.
