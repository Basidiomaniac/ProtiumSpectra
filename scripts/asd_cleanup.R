# Script to import and prepare .asd files from field spectrometer
# Code adapted from appetizer_course_dm20251103.R
# Adapted by Laurel Miller, Nov 2025
# Inputs: folder(s) full of .asd files, metadata in .csv format
# Outputs: Dataframes of spectral data with metadata, functions to lengthen
# for plotting and a plotting function using ggplot2.

library(tidyverse)
library(asdreader)
library(viridis)
library(prospectr)
library(lubridate)

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

spectra_cleanup = function(df) {
  spectra = df
  colnames(spectra) <- paste("nm", colnames(spectra), sep = "")
  rownames(spectra) <- sub(".*(?=Protium)", "", rownames(spectra), perl = TRUE)
  spectra_clean = as.data.frame(spectra) %>% 
    mutate(sample = str_split_i(rownames(spectra), 'Protium', 2) %>% 
             str_split_i(., 'x', 1) %>%
             as.numeric,
           replicate = str_split_i(rownames(spectra), 'x', 2) %>% 
             str_split_i(., '.asd', 1) %>%
             as.numeric) %>%
    filter(!if_any(everything(), ~ str_detect(.x, "Inf"))) %>%
    filter(!is.na(sample)) # remove white reading NA's
  return(spectra)
}

# FILL IN DIRECTORIES FOR YOUR OWN FILES
asd_dir = '/home/girlmunculus/spectraProj/asd/' # location of .asd files
github_dir = '/home/girlmunculus/GitHub/ProtiumSpectra/' # location of GitHub project

# some basic test commands using provided test file
test_file = asd_file()
test_spectra = get_spectra(test_file)
head(test_spectra)
plot(as.numeric(test_spectra), type = 'l')

# read .asd files
setwd(asd_dir)
filepaths = list.files(pattern = "\\.asd$", 
                       recursive = TRUE, full.names = TRUE)

# go through files, adding collection date and merging to main
spectra = get_spectra(filepaths, "reflectance")
# plot(as.numeric(spectra), type = 'l')

# clean nir data up
colnames(spectra) <- paste("nm", colnames(spectra), sep = "")
rownames(spectra) <- sub(".*(?=Protium)", "", rownames(spectra), perl = TRUE)
spectra_clean = as.data.frame(spectra) %>% 
                filter(!if_any(everything(), ~ str_detect(.x, "Inf"))) %>%
                filter(if_any(everything(), ~ . < 1)) %>%
                filter(!is.na(sample)) # remove white reading NA's

# splicing to fix steps in measurements caused by spectrometer's three scanners
# NOTE: This currently only fixes the 1000nm gap, not the 1801nm one.
spectra_spliced = spliceCorrection(X = spectra_clean, wav = 350:2500)
spectra_cleaner = spectra_spliced %>%
        as.data.frame() %>% 
        mutate(sample = str_split_i(rownames(spectra_clean), 'Protium', 2) %>% 
         str_split_i(., 'x', 1) %>%
         as.numeric,
       replicate = str_split_i(rownames(spectra_clean), 'x', 2) %>% 
         str_split_i(., '.asd', 1) %>%
         as.numeric) %>%
    dplyr::select(sample, replicate, starts_with('nm')) %>% 
  mutate(sample = str_pad(sample, width = 3, side = "left", pad = "0"))

# turn spectra into long df for plotting (takes a bit)
# spectra_long = lengthen(spectra_clean)

# plot spectral chart (this takes a while)
# you can change 'color =' to groupings of interest (species, etc.)
# spectra_plot(spectra_long, spectra_long$sample)

# stepclass here


# read in the metadata on each sample
setwd(github_dir)
meta = read.csv('subserratum_NIR_project.csv')
species = meta %>% mutate(sample = str_split_i(ViewSpecFile.Folder, 'Protium', 2) %>%
                            as.numeric) %>%
  rename(population_lineage = population.lineage,
         geographic_location = geographic.location,
         date = DataFolder) %>%
  dplyr::select(sample, species, population_lineage, geographic_location, date) %>% 
  mutate(sample = str_pad(sample, width = 3, side = "left", pad = "0"))

# combine species metadata with nir data
nir = spectra_cleaner %>% left_join(species, by = 'sample') %>%
  dplyr::select(sample, replicate, species, population_lineage, 
                geographic_location, date, starts_with('nm')) %>% 
  mutate(date_parsed = parse_date_time(date, orders = "B y")) %>%
  mutate(date_num = as.numeric(date_parsed))  # make date numeric for LDA later
# test plot spectra by species
nir_long = lengthen(nir)
spectra_plot(nir_long, nir_long$species)

# now, open asd_identify.R to start linear discriminant analysis (LDA)
# for species identification.
