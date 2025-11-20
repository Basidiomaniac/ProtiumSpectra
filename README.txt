This project contains code for importing, cleaning, and analyzing NIR spectral data in the form of .asd files.
Scripts provided for producing LDA models to discriminate between taxonomic categories based on NIR data.
Code for this project was written by Niko Darci-Maher, Laurel Miller
Last updated 2025

No installation required, just install libraries and run R scripts.

Required input:
-.asd files containing NIR spectral data. You should edit the filepaths in asd_cleanup.R to match your file locations.
asd_cleanup.R searches recursively, so it's okay if they're in subfolders. The main requirement here is .asd files,
which contain one row per reading, and one column per wavelength measured. We configured our spectrometer to name
files as "path/Date/Protium###x######.asd" representing sample number and reading number per sample. So
Protium003x000011.asd is the 11th reading on sample 3. To group multiple readings per sample, you'll need to name
your rows of data appropriately.
-metadata. You can see an example of ours in NIR_test_data_metadata.csv. Line 113 (species = ...) allocates 
metadata from the meta dataframe into your main spectral data. For any metadata you have, just fill in the
rename portion of the command.

Scripts:
asd_cleanup.R - Code for importing, cleaning, and combining spectral data files with metadata. After this code,
you're ready to begin building models and analyzing your data.
asd_identify.R - Code to build LDA models, test success of averaged versus non-averaged models, and plot
results on an abstracted "LD space" of your wavelength data.
asd_analyze.R - Code to analyze wave data. Still in progress.

Example data files:
March 20 - Contains Protium samples 1-52. Can be used with the "appetizer course" code for testing species ID
April 23 - Protium 100 - 196 
April 30 - Protium 197 - 291 (210 has no leaves, so no sample!)
May 7 - Protium 292 - 450
These are all of the herbarium specimens from the sheets. Not all readings have 12 measurements.

This project has two courses: An Appetizer Course and a Main Course.
The Appetizer Course goes with March20th data folder (plus one reading in the April 23rd folder).
The Main Course goes with google sheet “subserratum NIR project” and data folders April23, April30, and May7.
The Appetizer Course is to see if we can use discriminant analysis LDA to separate out Protium species. 
We set up a training data set of 10 different species.  These 10 species are listed in the google sheet. 
Each species includes 5-6 individuals.  We measured each individual 12 times (3 leaves, 4 measures per leaf, 2 on bottom side, 2 on top side). 
We can see how much variation there is within measurements of the same leaf, the same individual, and among individuals of the same species. 
We can then use the “mystery plant” (and also there are other repeats of these 10 species in the other data folders) to see if the computer can 
accurately place them into the right species. The Main Course is to see how FT-NIR data will vary among lineages of a species complex, 
and if there is phylogenetic signal in spectral data. See the Misiewicz et al. 2023, Molecular Ecology paper to see a phylogeny of these lineages.
You can see in the datasheet in column F and G the species and lineages of all the samples.

Guide to the Data Folders and Google Sheets
Each sample has 12 readings.  These are four scans from 3 leaflets each (one on the bottom side base, one on the bottom side tip, 
one on the top side base, one on the top side tip).  The sample number corresponds with column 

Note: Protium159 contains a mystery plant that corresponds to one of the species from Protium 1-52-- if you guys can figure out what species 
Protium159 from the spectral data – then we will know our method works!  (Appetizer Course) 


