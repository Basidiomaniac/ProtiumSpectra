# Script to analyze patterns in NIR spectra from samples. Analysis is 
# intra- and inter- specific, plus connecting our most informative wavelengths
# to metabolite data(?)
# Last updated November 2025, Laurel Miller

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
subserratum = nir_lda_data %>% filter(Species == "subserratum")

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


# next, let's try hierarchical clustering, which we can compare to the 
# group's phylogenetic tree