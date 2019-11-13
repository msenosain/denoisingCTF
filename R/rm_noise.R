#' Noise removal function
#'
#' This function 
#'
#' @param file_type A regular expression for the file extension. Only file names 
#'        which match the regular expression will be returned. FCS only.
#' @param rm_beads Logical. If \code{TRUE} beads are removed. Default is \code{TRUE}.
#' @param rm_debris Logical. If \code{TRUE} debris is removed. Default is \code{TRUE}.
#' @param use.current.model.beads Logical. If \code{FALSE}, an object of class \code{train} can 
#'        be passed to the argument \code{model_beads}. Default is \code{TRUE}.
#' @param use.current.model.debris Logical. If \code{FALSE}, an object of class \code{train} can 
#'        be passed to the argument \code{model_beads}. Default is \code{TRUE}.
#' @param model_beads Object of class \code{train}. Needed if \code{use.current.model.beads} 
#'        is \code{FALSE}.
#' @param model_debris Object of class \code{train}. Needed if \code{use.current.model.debris} 
#'        is \code{FALSE}.
#' @param alg_db A character vector with the name of the algorithm used to train
#'        \code{model_beads}. It can be 'RF' for Random Forest or 'XGB' for XGBoost.
#' @param alg_bd A character vector with the name of the algorithm used to train
#'        \code{model_beads}. It can be 'RF' for Random Forest or 'XGB' for XGBoost.
rm_noise <- function(file_type = '.fcs|.FCS', rm_beads = TRUE, rm_debris = TRUE,
                     use.current.model.beads = TRUE, use.current.model.debris = TRUE,
                     model_beads = model_beads, model_debris = model_debris, 
                     alg_db = 'RF', alg_bd = 'RF'){

    # Read modelS
    if(use.current.model.beads){
        model_beads.file <- system.file('data', 'model_rf_beads.rds', package='denoisingCTF',
                                mustWork=TRUE)
        model_beads <- readRDS(model_beads.file)
    }

    if(use.current.model.debris){
        model_debris.file <- system.file('data', 'model_rf_debris.rds', package='denoisingCTF',
                                mustWork=TRUE)
        model_debris <- readRDS(model_debris.file)
    }
    
    # Read files list
    fcs <- c('.FCS', '.fcs', '.fcs|.FCS')
    if (!(file_type %in% fcs)){
      stop("Data type not supported!")
    }
    files_list <- list.files(pattern=file_type)
    print('File list successfully created!')

    # Create output folder
    output_path <- file.path(getwd(), 'output')
    if (!dir.exists(output_path)){
      dir.create(output_path, recursive = TRUE)
      message("\nOutput directory created\n")
    }

    # Read first file to correct column order later
    smp <- flowCore::read.FCS(files_list[1], transformation = FALSE)
    descrp <- smp@parameters@data$desc
    smp <- data.frame(smp@exprs)
    col_nms <- colnames(smp)

    # Ask for channels to remove zeros
    cn <- as.matrix(paste0(col_nms, '_', descrp))
    if(rm_debris){
        print(cn)
        prompt1 <- "Enter the column INDICES of the 'mandatory' markers 
        (separated by single space only, no comas allowed) \n"
        prompt2 <- "Enter the column INDICES of the 'optional' markers 
        (separated by single space only, no comas allowed) \n"
        mand_idx <- as.numeric(strsplit(readline(prompt1), " ")[[1]])
        opt_idx <- as.numeric(strsplit(readline(prompt2), " ")[[1]])        
    }


    # channels to remove beads 
    if(rm_beads){
        #print(cn)
        prompt <- "Enter the column INDICES of the beads channels Ce140, Eu151, 
        Eu153, Ho165, Lu175 (separated by single space only, no comas allowed) \n"
        ft_beads <- col_nms[as.numeric(strsplit(readline(prompt), " ")[[1]])]
    }

    # channels to remove debris
    if(rm_debris){
        #print(cn)
        prompt <- "Enter the column INDICES of the gaussian parameters channels 
        'Event_length', 'Center', 'Offset', 'Residual', 'Width' and intact-cells 
        marker channel (separated by single space only, no comas allowed) \n"
        ft_debris <- col_nms[as.numeric(strsplit(readline(prompt), " ")[[1]])]

    }


    # Iterate between files to remove noise
    for(i in 1:length(files_list)){
        # Read exprs from FCS file
        dt <- data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)

        # Order columns
        dt <- dt[col_nms]

        k <- c()
        numcols <- ncol(dt)
        
        # Remove beads
        if(rm_beads){
            k2 <- predict_cl(dt, features = ft_beads, model = model_beads, 
                alg = alg_bd, label = 'beads')
            k <- c(k, k2)
            dt <- cbind(dt, beads = 0)
            dt[k2, 'beads'] <- 1
        }

        # Remove debris
        if(rm_debris){
            dt <- rm_zeros(dt, mand_idx = mand_idx, opt_idx = opt_idx)
            k1 <- predict_cl(dt, model = model_debris, features = ft_debris, 
                alg = alg_db, label = 'debris')
            k <- c(k, k1)
            dt <- cbind(dt, debris = 0)
            dt[k1, 'debris'] <- 1
        }

        # Create output folder for CSV files
        output_path2 <- file.path(output_path, 'noiseCL')
        if (!dir.exists(output_path2)){
          dir.create(output_path2, recursive = TRUE)
        }
        # Write csv file
        filename <- tools::file_path_sans_ext(basename(files_list[i]))
        write.csv(dt, file = file.path(output_path2, paste0(filename, 
            '_noiseCL.csv')), row.names=FALSE)
        
        # remove k rows from dt
        k <- unique(k)
        dt <- dt[-k,1:numcols]

        # Write new FCS
        cytofCore::cytofCore.write.FCS(as.matrix(dt), 
            file = file.path(output_path, paste0(filename, '_denoised.FCS')),
            what = "numeric", channelDescriptions = descrp)

    }

}