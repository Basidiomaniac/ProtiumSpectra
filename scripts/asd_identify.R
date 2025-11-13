# Script to conduct linear discriminant analysis (LDA) on spectral dataframes
# Code adapted from appetizer_course_dm20251103.R (author: Niko Darci-Maher)
# Adapted by Laurel Miller, Nov 2025
# Input: Spectral dataframe from asd_cleanup.R (nir)
# Output: LDA model capable of discriminating species at 99% accuracy

library(viridis)
library(klaR)
library(lme4)

# color by species and plot all sample spectra (mult. per leaf and species)
nir_sample = nir %>% group_by(sample) %>%
  summarize_at(vars(starts_with('nm')), mean)
nir_sample = nir_sample %>% pivot_longer(cols = starts_with('nm'),
               names_to = 'wavelength',
               names_prefix = 'nm',
               values_to = 'reflectance_mean') %>%
              left_join(nir 
              %>% dplyr::select(sample, species) %>% unique, by = 'sample') %>%
              mutate(wavelength = as.numeric(wavelength),
              sample = as.character(sample))

line_plot_sample = ggplot(nir_sample, aes(x = wavelength, y = reflectance_mean, color = species)) +
  geom_line(aes(group = sample))
line_plot_sample


# means for each species, plot
nir_species = nir %>% group_by(species) %>%
  summarize_at(vars(starts_with('nm')), mean)
nir_species = nir_species %>% pivot_longer(cols = starts_with('nm'),
                                           names_to = 'wavelength',
                                           names_prefix = 'nm',
                                           values_to = 'reflectance_mean') %>%
  mutate(wavelength = as.numeric(wavelength))
line_plot_species = ggplot(nir_species, aes(x = wavelength, y = reflectance_mean, color = species)) +
  geom_line(aes(group = species))
line_plot_species

# plot with density by species
nir_long = lengthen(nir)
line_plot_species = ggplot(nir_long, aes(x = wavelength, y = intensity, color = species)) +
  geom_line(aes(group = sample), alpha = 0.2) +
  facet_wrap(~ date)
line_plot_species

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
normaldist_plot

# hmm, measurements across dates seem to vary, but this is also because we 
# sampled different species on different days. Let's compare within species
# across dates.
date_lmer = lmer(intensity ~ date + (1 | species), data = nir_long)
summary(date_lmer)

# no surprise: species and date are sig. predictors of intensity,


# separate out wavelength data, produce LDA model!
# pull the first for each species as test data, reserve the rest as training data
nir_training = nir %>% group_by(species, date) %>%
                  filter(sample != min(sample)) # all but first sample per species
nir_testing = nir %>% group_by(species, date) %>%
                  filter(sample == min(sample)) # ONLY first sample per species
wavelengths = c("date_num", grep("^nm", names(nir_training), value = TRUE))
lda_model = lda(x = nir_training[, wavelengths], grouping = nir_training$species)

# make the species prediction!
species_predict = MASS:::predict.lda(lda_model, newdata = nir_testing[, wavelengths])

# how did we do on actual species accuracy?
accuracy = mean(species_predict$class == nir_testing$species)
accuracy # for test dataset, accuracy is 99.08%!

# Get LD1 and LD2 scores for species
lda_data <- data.frame(
  Species = nir_testing$species,
  Date = nir_testing$date_num,
  LD1 = species_predict$x[, 1],
  LD2 = species_predict$x[, 2] # Only the first two LDs are used
)

# plot LD in 2D space
lda_space_plot <- ggplot(lda_data, aes(x = LD1, y = LD2, color = Species)) +
  geom_point(alpha = 0.9) +
  stat_ellipse(geom = "polygon", type = "norm", aes(fill = Species), alpha = 0.1, show.legend = TRUE) +
  labs(
    title = "LDA Coverage of Protium spp.",
    subtitle = paste0(
      "LD1 explains ", round(lda_model$svd[1]^2 / sum(lda_model$svd^2) * 100, 1), 
      "% of variance\n", " LD2 explains ", round(lda_model$svd[2]^2 / sum(lda_model$svd^2) * 100, 1),
      "% of variance"
    ),
    x = "Linear Discriminant 1 (LD1)",
    y = "Linear Discriminant 2 (LD2)"
  ) +
  scale_x_continuous(limits = c(-50, 40)) +
  scale_y_continuous(limits = c(-50, 40)) +
  facet_wrap(~ Date)
lda_space_plot # looks like the LDA covers a decent amount of variance, species are well-separated

# let's see how good ALL our LD's are
ld_variance = round(lda_model$svd^2 / sum(lda_model$svd^2) * 100, 1) %>% data.frame() 
colnames(ld_variance) = "explains"
rownames(ld_variance) = c("LD1", "LD2", "LD3", "LD4", "LD5", "LD6", "LD7", "LD8")
ld_variance$ld = rownames(ld_variance)
lda_bar_plot = ggplot(ld_variance, aes(x = ld, y = explains)) +
  geom_col(fill = "#13536a", color = "#08191f") + 
  labs(
    title = "Linear Discriminant Effectiveness",
    x = "Linear Discriminant",
    y = "% of Variance Explained")
lda_bar_plot

# map boundaries. these will need to be adjusted to fit your values.
x_min = -50
x_max = 40
y_min = -50
y_max = 40
grid = expand.grid(LD1 = seq(x_min, x_max, by = 0.1),
                    LD2 = seq(y_min, y_max, by = 0.1))

# make a NEW LDA just trained on LD1 and LD2 from original model.
# this "flattens" our wavelengths so we can predict boundaries
lda_boundary_model <- lda(Species ~ LD1 + LD2, data = lda_data) 
grid$Species = MASS:::predict.lda(lda_boundary_model, newdata = grid)$class
tile_plot = ggplot() + 
  geom_tile(data = grid, aes(x = LD1, y = LD2, fill = Species), alpha = 0.8) +
  geom_point(
    data = lda_data,
    aes(x = LD1, y = LD2, color = Species)
  ) +
  labs(
    subtitle = paste0(
      "LD1 explains ", round(lda_model$svd[1]^2 / sum(lda_model$svd^2) * 100, 1), "% of variance\n",
      "LD2 explains ", round(lda_model$svd[2]^2 / sum(lda_model$svd^2) * 100, 1), "% of variance"
    ),
    x = "Linear Discriminant 1",
    y = "Linear Discriminant 2"
  ) +
  scale_color_viridis_d(begin = 0, end = 0.9, option = "A") +
  scale_fill_viridis_d(option = "C",) +
  facet_wrap(~ Date) +
  theme_classic()
tile_plot
# for the test data, we see contours that encapsulate our data's LD's pretty well.
# there are a couple samples from paniculatum, ferrugineum, and subserratum that
# land outside their species' bounds.
# NEXT: open asd_analysis.R to produce statistics on wavelength variation and
# potential connections to the metabolome, environment, etc.
