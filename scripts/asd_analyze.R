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