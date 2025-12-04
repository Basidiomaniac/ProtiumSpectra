# Script to analyze patterns in NIR spectra from samples. Analysis is 
# intra- and inter- specific, plus connecting our most informative wavelengths
# to metabolite data(?)
# Last updated November 2025, Laurel Miller

library(sf)
library(ggmap)
library(leaflet)

# test PCA of all species across all wavelengths
nir_num = nir %>% dplyr::select(starts_with('nm'))
pca_all = prcomp(nir_num, center = T, scale. = T)
nir = cbind(nir, pca_all$x[,1:5])
pcaplot = ggplot(nir, aes(x = PC1, y = PC2, color = species)) + 
  geom_point() +
  scale_color_viridis_d(option = "A")
pcaplot

# # try identification on clade level
# clade_sample = nir_sample[nir_sample$specieslineage != "",]
# clade_avg_train = nir_sample %>% group_by(specieslineage) %>%
#   filter(sample != (min(sample))) # all but first sample per species
# clade_avg_test = nir_sample %>% group_by(specieslineage) %>%
#   filter(sample == (min(sample))) # ONLY first sample per species
# clade_lda_model = lda(x = clade_avg_train[, wavelengths], grouping = clade_avg_train$specieslineage)
# clade_avg_predictions = MASS:::predict.lda(clade_lda_model, newdata = clade_avg_test[, wavelengths])
# mean(clade_avg_predictions$class == clade_avg_test$specieslineage)

# looking back at the subserratum species complex, let's try our hand at
# clustering in LDA space
# first, spectral clustering, to see a spatial representation of similarity
subserratum_lda = nir_lda_data %>% filter(Species == "subserratum")

# similarity matrix across LD1 and LD2
similarity_matrix <- exp(-dist(subserratum[, 3:4])^2 / (2 * 1^2))

#Compute Eigenvalues and Eigenvectors
eigen_result <- eigen(similarity_matrix)
eigenvalues <- eigen_result$values
eigenvectors <- eigen_result$vectors

#Choose the First k Eigenvectors
k = 9 # 9 clades were assigned genomically in Misiewicz et al. 2023
selected_eigenvectors <- eigenvectors[, 1:k]

# k-mean cluster, add to main dataframe
cluster_assignments <- kmeans(selected_eigenvectors, centers = k)$cluster
subserratum$Cluster <- factor(cluster_assignments)

# visualize clusters
ggplot(subserratum, aes(LD1, LD2, color = Cluster, label = Species)) +
  geom_point() +
  scale_color_viridis_d(option = "C")

# looking at lineages
subserratum = nir[nir$species == "subserratum",]
subserratum = subserratum %>%
  filter(genome == TRUE)
sub_lineage = subserratum %>% group_by(population_lineage) %>%
  summarize_at(vars(starts_with('nm')), mean)
sub_lineage = sub_lineage %>% pivot_longer(cols = starts_with('nm'),
                                           names_to = 'wavelength',
                                           names_prefix = 'nm',
                                           values_to = 'reflectance_mean') %>%
  mutate(wavelength = as.numeric(wavelength))
line_lineage = ggplot() +
  geom_line(data = sub_lineage, aes(x=wavelength, y=reflectance_mean, 
                color=population_lineage, group = population_lineage)) +
  scale_color_viridis_d(option = "plasma")
line_lineage # plot out spectra by lineage (genome confirmed)

# remake model for subserratum lineage, plot out
subserratum_training = subserratum %>% group_by(population_lineage) %>%
  filter(sample != min(sample))
subserratum_testing = subserratum %>% group_by(population_lineage) %>%
  filter(sample == min(sample)) 
lineage_model = lda(x = subserratum_training[, wavelengths], 
                    grouping = subserratum_training$population_lineage)

# model, predict lineage
species_predict = MASS:::predict.lda(lineage_model, newdata = subserratum_testing[, wavelengths])
accuracy = mean(species_predict$class == subserratum_testing$population_lineage)
accuracy

# Get LD1 and LD2 scores for species
lineage_lda_data <- data.frame(
  Species = subserratum_testing$population_lineage,
  Sample = subserratum_testing$sample,
  LD1 = species_predict$x[, 1],
  LD2 = species_predict$x[, 2] # Only the first two LDs are used
)

urchin_plot(lineage_lda_data, lineage_lda_data, -10, 10, -10, 10)
# meh, 36% accurate. averaging actually brings success DOWN here.
# spectra might just be too messy for lineage level data?

subserratum$whitesand = subserratum$population_lineage
whitesand_train = subserratum %>% group_by(whitesand) %>%
  filter(sample != min(sample))
whitesand_test = subserratum %>% group_by(whitesand) %>%
  filter(sample == min(sample)) 
lineage_model = lda(x = subserratum_training[, wavelengths], 
                    grouping = subserratum_training$population_lineage)

# which lineages are most identifiable when considered as binary (A vs. not A, etc.)
# model, predict lineage
lineage_results = list()
lineages = unique(confirmed$population_lineage)
lineages = lineages[-5]
confirmed$population_lineage[confirmed$population_lineage == ""] = "ferrugineum"
accuracy_vector = c()
for (lineage in lineages) {
  confirmed$binary_label = ifelse(confirmed$population_lineage == lineage, lineage, "OTHER")
  lineage_train = confirmed %>% group_by(population_lineage) %>%
    filter(sample != min(sample))
  lineage_test = confirmed %>% group_by(population_lineage) %>%
    filter(sample == min(sample)) 
  lineage_model = lda(x = lineage_train[, wavelengths], 
                      grouping = lineage_train$binary_label)
  species_predict = MASS:::predict.lda(lineage_model, newdata = lineage_test[, wavelengths])
  accuracy = mean(species_predict$class == lineage_test$binary_label)
  accuracy_vector = c(accuracy_vector, accuracy)
}

# summarize success by lineage in binary classification
names(accuracy_vector) = lineages
accuracy_vector

# map coords of data for the Misiewicz et al. dataset
lin_color <- colorFactor(palette = 'RdYlGn', nir$population_lineage)
confirmed = nir %>% filter(genome == TRUE)
leaflet() %>% 
  addTiles() %>% 
  addCircleMarkers(data = confirmed, lat = ~lat, lng = ~long,
                   color = ~lin_color(population_lineage),
                   popup = confirmed$population_lineage) %>% 
  addLegend('bottomright', pal = lin_color, values = confirmed$population_lineage,
            title = 'Population Lineages',
            opacity = 1)
