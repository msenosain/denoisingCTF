# denoisingCTF: a CyTOF denoising pipeline

Here we present a noise removal pipeline for CyTOF data. You can either use the 
main function `rm_noise` to remove the noise applying the current trained models,
or train your own models and then call them into the main function.

## Installation
You will need to install the `devtools` CRAN package. All other dependencies 
will be automatically installed with the package installation.

```R
# Install devtools
install.packages("devtools")

# Install denoisingCTF
devtools::install_github("msenosain/denoisingCTF")
```

## Usage

### Denoising CyTOF data

This function will remove zeros, beads and debris from the FCS files in the 
current working directory using previously trained models. 

```R
library(denoisingCTF)
# Set working directory in folder with files to be denoised
rm_noise(file_type = '.fcs|.FCS', rm_beads = TRUE, rm_debris = TRUE,
    use.current.model.beads = TRUE, use.current.model.debris = TRUE,
    alg_db = 'RF', alg_bd = 'RF')

# Mandatory: 52
# Optional: 20 29 46 34 44 47 15 37 40 31 49
# Beads: 14 26 28 41 51
# GP_dc: 2 52 73 74 75 76

```

### Preprocessing datasets and model training

If you want to train your own classification models for beads and debris, you
can do it following these steps:

1. Beads Model:
```R
# Unsupervised beads detection and training/testing data sets generation.
Beads_TrainTest(sample_size = 40, method = 'k_means', 
    bsample = 5000, class_col = 'BeadsSmp_ID')

# Train model
TrainModel(train_set, test_set, alg = 'all', class_col = 'BeadsSmp_ID', 
    seed = 40, name_0 = 'cells', name_1 = 'beads', label = 'beads',
    allowParallel = T, free_cores = 2)

```

2. Debris model:
``` R
# Removal of zeros, beads and addition of row ID column
pre_gate(sample_size=30, model_beads=model_beads, alg_bd = 'RF') 

# Manually gate debris and dead cells using Gaussian Parameters and a live/dead cell marker.

# Comparison of pre-gated and post-gated files, noise labeling and training/testing data sets generation.
post_gate(bsample = 5000, path_pregated = '../') 

# Train model
TrainModel(train_set, test_set, alg = 'all', class_col = 'GP_Noise', 
    seed = 40, name_0 = 'cells', name_1 = 'debris', label = 'debris', 
    allowParallel = T, free_cores = 2)

```

## Tutorial

For a detailed example, see our tutorial.


## Debris and Beads model training and evaluation

For further details on these models, click here.





