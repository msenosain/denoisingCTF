#' Preprocess data before debris manual gating
#'
#' This function selects a random sample and reads the FCS files, remove zeros 
#' and beads and then adds a row ID column and writes a new FCS file. It uses 
#' files on current working directory. Files should have same parameters 
#' desc/column names. 
#'
#' @param sample_size A numeric vector. Number of samples (files) to be used.
#' @param use.current.model Logical. If \code{FALSE}, an object of class \code{train} can 
#'        be passed to the argument \code{model_beads}. Default is \code{TRUE}.
#' @param model_beads Object of class \code{train}. Needed if \code{use.current.model} 
#'        is \code{FALSE}.
#' @param alg_bd A character vector with the name of the algorithm used to train
#'        \code{model_beads}. It can be 'RF' for Random Forest or 'XGB' for XGBoost.
#' @export pre_gate()
pre_gate <- function(sample_size = 20, use.current.model = TRUE, model_beads=model_beads, 
    alg_bd = 'RF'){
    # Read model
    if(use.current.model){
        model_beads.file <- system.file('data', 'model_rf_beads.rds', package='denoisingCTF',
                                mustWork=TRUE)
        model_beads <- readRDS("model_rf_beads.rds")
    }


    # Read files list
    files_list <- list.files(pattern='.FCS|.fcs')
    print('File list successfully created!')

    # Create output directory
    output_path <- file.path(getwd(), 'toy_debris')
    if (!dir.exists(output_path)){
      dir.create(output_path, recursive = TRUE)
      message("\nOutput directory created\n")
    }

    # Read first file and get column order to correct later
    smp <- flowCore::read.FCS(files_list[1], transformation = FALSE)
    col_nms <- colnames(data.frame(smp@exprs))
    desc_new <- smp@parameters@data$desc
    cn <- as.matrix(paste0(col_nms, '_', desc_new))

    # Obtain descriptions and add new
    length_desc <- length(desc_new)
    names_desc <- names(desc_new)
    desc_new <- c('ID', desc_new)
    names(desc_new) <- c(names_desc, paste('$P', (length_desc+1), 'S', sep = '' ))


    # Ask for channels to remove zeros    
    print(cn)
    prompt1 <- "Enter the column INDICES of the 'mandatory' markers (separated by single space only, no comas allowed) \n"
    prompt2 <- "Enter the column INDICES of the 'optional' markers (separated by single space only, no comas allowed) \n"
    mand_idx <- as.numeric(strsplit(readline(prompt1), " ")[[1]])
    opt_idx <- as.numeric(strsplit(readline(prompt2), " ")[[1]]) 

    # Ask for channels to remove beads
    prompt3 <- "Enter the column INDICES of the beads channels Ce140, Eu151, Eu153, Ho165, Lu175 (separated by single space only, no comas allowed) \n"
    ft_beads <- col_nms[as.numeric(strsplit(readline(prompt3), " ")[[1]])]
    
    # Extract sample of files list
    files_list <- sample(files_list, size = sample_size)

    # Process files iteratively
    for (i in 1:length(files_list)){
        # Read exprs from FCS file
        dt <- data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)

        # Order columns
        dt <- dt[col_nms]

        # Remove zeros
        dt <- rm_zeros(dt, # 52 ||  20 29 46 34 44 47 15 37 40 31 49
            mand_idx = mand_idx, opt_idx = opt_idx)

        # Remove beads
        k <- predict_cl(dt, model = model_beads, features = ft_beads, 
                alg = alg_bd, label = 'beads') # 14 26 28 41 51
        dt <- dt[-k,]

        # Add row ID
        dt <- tibble::rowid_to_column(dt, var = 'ID')

        # Write new FCS
        filenm <- tools::file_path_sans_ext(basename(files_list[i]))
        cytofCore::cytofCore.write.FCS(as.matrix(dt), 
            file = file.path(output_path, paste0(filenm, '_tdb.FCS')),
            what = "numeric", channelDescriptions = desc_new)
    }

}



#' Post-gate debris labeling
#'
#' This function compares the pre-gated files with the gated files and adds the
#' label of noise and cells to the events in each file. Then it uses 
#' \code{BalancedSample} and \code{TrainTest} to generate a class-balanced 
#' training and test sets. The working directory should be the folder that 
#' contains the gated files and the pre-gated files should be one directory up.
#'
#' @param bsample A numeric vector. Size of the sample from each class for 
#'        to be passed to \code{BalancedSample} function.
#' @param path_pregated A character vector with full path names. Default is the 
#'        one directory up.
#' @export post_gate()
post_gate <- function (bsample = 5000, path_pregated = '../'){

    # Create output folder
    output_path <- file.path(getwd(), 'labeled')
    if (!dir.exists(output_path)){
      dir.create(output_path, recursive = TRUE)
      message("\nOutput directory created\n")
    }

    # Read pre-gated files list
    preg_flist <- list.files(path = path_pregated, pattern = '.FCS|.fcs')
    print('Pre-gate file list successfully created!')

    # Read post-gated files list
    postg_flist <- list.files(pattern='.FCS|.fcs')
    print('Post-gate file list successfully created!')

    # gsub('_.*', '',x) to select only ptID
    # names(preg_flist) <- gsub('_.*', '',preg_flist)
    # names(postg_flist) <- gsub('_.*', '',postg_flist)

    # Order list elements
    preg_flist <- preg_flist[order(preg_flist)]
    postg_flist <- postg_flist[order(postg_flist)]

    # Read first file and get column order to correct later
    smp <- flowCore::read.FCS(paste0('../', preg_flist[1]), transformation = FALSE)
    col_nms <- colnames(data.frame(smp@exprs))
    desc_new <- smp@parameters@data$desc

    # Obtain descriptions and add new
    length_desc <- length(desc_new)
    names_desc <- names(desc_new)
    desc_new <- c(desc_new[2:length_desc], 'GP_Noise')
    names(desc_new) <- names_desc


    dt_ls <- list()
    for (i in 1:length(preg_flist)){
        
        # Read files
        dt <- data.frame(flowCore::read.FCS(paste0('../', preg_flist[i]), transformation = FALSE)@exprs)
        dt_gated <- data.frame(flowCore::read.FCS(postg_flist[i], transformation = FALSE)@exprs)
        
        # Order columns
        dt <- dt[col_nms]
        dt_gated <- dt_gated[col_nms]

        # indexes for rows gated out (noise)
        idx <- which(!(dt$ID %in% dt_gated$ID)) 
        dt['GP_Noise'] <- 0 #'cells'
        dt[idx, 'GP_Noise'] <- 1 #'noise'

        # Remove ID column
        dt['ID'] <- NULL

        # Write new FCS
        filenm <- tools::file_path_sans_ext(basename(preg_flist[i]))
        cytofCore::cytofCore.write.FCS(as.matrix(dt), 
            file = file.path(output_path, paste0(filenm, '_lb.FCS')),
            what = "numeric", channelDescriptions = desc_new)

        # Get balanced sample
        dt <- BalancedSample(dt, sample_size = bsample, class_col = 'GP_Noise')
        dt_ls[[i]] <- dt
   
    }
    # Add names to list
    names(dt_ls) <- tools::file_path_sans_ext(basename(preg_flist)) 

    # Compute train/test sets
    TrainTest(dt_ls, output_path = output_path, label = 'debris', class_col = 'GP_Noise')
}
