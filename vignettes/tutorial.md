# denoisingCTF Tutorial

This tutorial will guide you

## Getting started

## Preparing your data

## Denoising CyTOF data

This function will remove zeros, beads and debris from the FCS files in the 
current working directory using previously trained models. 

First, a prompt will ask the user for the "Mandatory markers", which can be any
marker that _should_ be expressed in the cells (e.g. post-fixation DNA intercalator).
Any cell/event/row that has zero expression of this marker is removed.

Then, a second prompt will ask the user for the "Optional markers", which can 
be markers used for cell type identification. In this case, any cell/event/row
that has zero expression of _all_ this markers, meaning that it cannot be 
identified as any existing cell type, is removed.



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

## Preprocessing datasets and model training
If you want to train your own classification models for beads and debris, you
can do it following these steps:

1. Beads Model:
```R
# Unsupervised beads detection
# 'Good' data is saved and training and test datasets are generated
Beads_TrainTest(sample_size = 40, method = 'k_means', 
    bsample = 5000, class_col = 'BeadsSmp_ID')

# Train model
TrainModel(train_set, test_set, alg = 'all', class_col = 'BeadsSmp_ID', 
    seed = 40, name_0 = 'cells', name_1 = 'beads', label = 'beads',
    allowParallel = T, free_cores = 2)

```

2. Debris model:
``` R
# Remove zeros and add row ID column
pre_gate(sample_size=30, model_beads=model_beads, alg_bd = 'RF') 

# Gate debris using Gaussiam Parameters
# Set working directory in the gated files folder
# Compares files and saves training and test datasets
post_gate(bsample = 5000, path_pregated = '../') 

# Train model
TrainModel(train_set, test_set, alg = 'all', class_col = 'GP_Noise', 
    seed = 40, name_0 = 'cells', name_1 = 'debris', label = 'debris', 
    allowParallel = T, free_cores = 2)

```

