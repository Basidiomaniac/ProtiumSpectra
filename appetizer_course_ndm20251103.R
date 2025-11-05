library(tidyverse)

theme_set(theme_bw())

setwd('~/src/protium_spectra/')

# read in the march20 data
nir = read.csv('data/march20_allsamples_reflectance.csv')
# nir1 = read.csv('data/march20_allsamples_reflectance_deriv1.csv')

# clean up the dataframe
summary(nir[,1:5])
nir = nir %>% mutate(Wavelength = paste0('nm', as.character(Wavelength)))
rownames(nir) = nir$Wavelength
nir = nir %>% select(-Wavelength)
nir = nir %>% t %>% data.frame
nir = nir %>% mutate(sample = str_split_i(rownames(nir), 'Protium', 2) %>% 
                         str_split_i(., 'x', 1) %>%
                         as.numeric,
                     replicate = str_split_i(rownames(nir), 'x', 2) %>% 
                         str_split_i(., '.asd', 1) %>%
                         as.numeric) %>%
    select(sample, replicate, starts_with('nm')) %>%
    filter(!is.na(sample)) # remove the white reference reading

# read in the metadata on each sample
meta = read.csv('data/subserratum_NIR_project.csv')
species = meta %>% mutate(sample = str_split_i(ViewSpecFile.Folder, 'Protium', 2) %>%
                              as.numeric) %>%
                rename(population_lineage = population.lineage,
                       geographic_location = geographic.location) %>%
                select(sample, species, population_lineage, geographic_location)

# merge species info onto nir info
nir = nir %>% left_join(species, by = 'sample') %>%
    select(sample, replicate, species, population_lineage, geographic_location, starts_with('nm'))

# take a peek at the relative variance of the different spectra within and between species

# do pca on all the wavelengths
nir_num = nir %>% select(starts_with('nm'))
pcres = prcomp(nir_num, center = T, scale. = T)
nir = cbind(nir, pcres$x[,1:5])

pcaplot = ggplot(nir, aes(x = PC1, y = PC2, color = species)) + geom_point()
pcaplot
ggsave('fig/protium_appetizer_pca_allwavelengths.png', pcaplot,
       width = 8, height = 6, units = 'in', dpi = 800)

# get the means for each sample and species
nir_sample = nir %>% group_by(sample) %>%
                summarize_at(vars(starts_with('nm')), mean)
nir_sample = nir_sample %>% pivot_longer(cols = starts_with('nm'),
                                           names_to = 'wavelength',
                                           names_prefix = 'nm',
                                           values_to = 'reflectance_mean') %>%
    left_join(nir %>% select(sample, species) %>% unique, by = 'sample') %>%
    mutate(wavelength = as.numeric(wavelength),
           sample = as.character(sample))

line_plot_sample = ggplot(nir_sample, aes(x = wavelength, y = reflectance_mean, color = species)) +
    geom_line(aes(group = sample))
line_plot_sample
ggsave('fig/protium_appetizer_spectra_bysample.png', line_plot_sample,
       width = 8, height = 6, units = 'in', dpi = 800)

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
ggsave('fig/protium_appetizer_spectra_byspecies.png', line_plot_species,
       width = 8, height = 6, units = 'in', dpi = 800)

# build an LDA model to predict species based on spectra

# check the normality of all the wavelengths
nir_norm = nir %>% mutate_at(vars(starts_with('nm')), scale)
nir_long = nir_norm %>% select(sample, replicate, species, starts_with('nm')) %>%
                    pivot_longer(cols = starts_with('nm'),
                                names_to = 'wavelength',
                                values_to = 'reflectance')
normaldist_plot = ggplot(nir_long, aes(x = reflectance, color = wavelength)) + 
    geom_density(alpha = 0.5) + 
    theme(legend.position = 'none')
normaldist_plot
ggsave('fig/protium_appetizer_normalized_reflectance_density.png', normaldist_plot,
       width = 8, height = 6, units = 'in', dpi = 800)

library(caret)

trainvect = createDataPartition(nir_norm$species, p = 0.9, list = F)
train = nir_norm[trainvect,] %>% select(species, starts_with('nm'))
test = nir_norm[-trainvect,] %>% select(species, starts_with('nm'))

fit_control = trainControl(method = 'repeatedcv',
                           number = 3,
                           repeats = 3)

lda_fit = train(species ~ ., data = train,
                method = 'lda',
                trControl = fit_control)
lda_fit

test_pred = predict(lda_fit, newdata = test)

sum(test_pred == test$species) / length(test$species)

# use the model to predict the species of the mystery plant (Protium159)
mystery = read.csv('data/april23_100_to_169_reflectance.csv')
# clean up the dataframe
summary(mystery[,1:5])
mystery = mystery %>% mutate(Wavelength = paste0('nm', as.character(Wavelength)))
rownames(mystery) = mystery$Wavelength
mystery = mystery %>% select(-Wavelength)
mystery = mystery %>% t %>% data.frame
mystery = mystery %>% mutate(sample = str_split_i(rownames(mystery), 'Protium', 2) %>% 
                         str_split_i(., 'x', 1) %>%
                         as.numeric,
                     replicate = str_split_i(rownames(mystery), 'x', 2) %>% 
                         str_split_i(., '.asd', 1) %>%
                         as.numeric) %>%
    select(sample, replicate, starts_with('nm')) %>%
    filter(!is.na(sample)) %>% # remove the white reference reading
    filter(sample == 159)

# all zeros??!!? hmmmm




