# Script to import and prepare .asd files from MODEL field spectrometer
# Code adapted from appetizer_course_dm20251103.R

library(tidyverse)
library(asdreader)
library(viridis)

# functions

# provide a df containing spectral data, rows as samples, columns as wavelengths
lengthen = function(df) {
  spectra_long <- df %>%
    pivot_longer(
      cols = starts_with("nm"),   # select wavelength columns
      names_to = "wavelength",
      values_to = "intensity"
    ) %>%
    mutate(
      wavelength = as.numeric(sub("nm", "", wavelength)) # convert nm350 → 350
    )
  return(spectra_long)
}
spectra_plot = function(df_long, category) {
  plot = ggplot(df_long, aes(x = wavelength, y = intensity, color = category, group = category)) +
    geom_line(size = 0.8, alpha = 0.8) +
    theme_minimal() +
    labs(x = "Wavelength (nm)", y = "Reflectance",) +
    theme(legend.position = "none") +
    if (is.numeric(df_long$category)) {
      plot <- plot + scale_color_viridis_c(option = "plasma")
    } else {
      plot <- plot + scale_color_viridis_d(option = "plasma")
    }
  plot
}

# FILL IN DIRECTORIES
asd_dir = '/home/girlmunculus/spectraProj/asd'
github_dir = '/home/girlmunculus/GitHub/ProtiumSpectra/'

# some basic test commands using provided test file
test_file = asd_file()
test_spectra = get_spectra(test_file)
head(test_spectra)
plot(as.numeric(test_spectra), type = 'l')

# read .asd files
setwd(asd_dir)
filepaths = list.files(pattern = "*.asd", 
                       recursive = TRUE, full.names = TRUE)
spectra = get_spectra(filepaths, "reflectance")
plot(as.numeric(spectra), type = 'l')

# clean nir data up
colnames(spectra) <- paste("nm", colnames(spectra), sep = "")
rownames(spectra) <- str_extract(rownames(spectra), "Protium.*")
spectra_clean = as.data.frame(spectra) %>% 
                mutate(sample = str_split_i(rownames(spectra), 'Protium', 2) %>% 
                str_split_i(., 'x', 1) %>%
                as.numeric,
              replicate = str_split_i(rownames(spectra), 'x', 2) %>% 
                str_split_i(., '.asd', 1) %>%
                as.numeric) %>%
                select(sample, replicate, starts_with('nm')) %>%
                filter(!if_any(everything(), ~ str_detect(.x, "Inf"))) %>%
                filter(!is.na(sample)) # remove white reading NA's

# turn spectra into long df for plotting
spectra_long = lengthen(spectra_clean)

# plot spectral chart (this takes a while)
# you can change 'color =' to groupings of interest (species, etc.)
spectra_plot(spectra_long, spectra_long$sample)

# read in the metadata on each sample
setwd(github_dir)
meta = read.csv('subserratum_NIR_project.csv')
species = meta %>% mutate(sample = str_split_i(ViewSpecFile.Folder, 'Protium', 2) %>%
                            as.numeric) %>%
  rename(population_lineage = population.lineage,
         geographic_location = geographic.location) %>%
  select(sample, species, population_lineage, geographic_location)

# combine species metadata with nir data
nir = spectra_clean %>% left_join(species, by = 'sample') %>%
  select(sample, replicate, species, population_lineage, geographic_location, starts_with('nm'))

# PCA of all species across all wavelengths
nir_num = nir %>% select(starts_with('nm'))
pca_all = prcomp(nir_num, center = T, scale. = T)
nir = cbind(nir, pca_all$x[,1:5])
pcaplot = ggplot(nir, aes(x = PC1, y = PC2, color = geographic_location)) + geom_point()
pcaplot

# plot spectra by species
nir_long = lengthen(nir)
spectra_plot(nir_long, nir_long$species)
