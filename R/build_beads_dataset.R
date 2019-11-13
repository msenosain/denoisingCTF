#' Merge beads-only and cells-only
#'
#' This function labels and merges beads-only and cells-only files collected 
#' on the same date separately. Column names in these files should be edited as 
#' needed to ensure are the same across files prior using this function. 
#' Working directory should be the folder that contain the files of interest. 
#' Files corresponding to the same date will be merged and new FCS files with 
#' a class column will be outputted. After the merge, files can be normalized
#' using a bead normalization software.
#'
BeadsIDmerge <- function(){
    files_list <- list.files(pattern='.FCS|.fcs')
    smp <- flowCore::read.FCS(files_list[1], transformation = FALSE)

    # Obtain descriptions and add new
    desc_new <- smp@parameters@data$desc
    length_desc <- length(desc_new)
    names_desc <- names(desc_new)
    desc_new <- c(desc_new, 'beads_ID')
    names(desc_new) <- c(names_desc, paste('$P', (length_desc+1), 'S', sep = '' ))

    # Extract data as list 
    dt <- extract_data() #function in basic_functions.R 

    # Column order
    colorder <- c(colnames(data.frame(smp@exprs)))
    dt <- dt %>% purrr::map(function(x) x<-x[colorder])

    # Add identifier
    for (i in names(dt)){
        if(any(grepl("control+", i))){
            dt[[i]]['beads_ID'] <- 0 #'cells'
        } else {
            dt[[i]]['beads_ID'] <- 1 #'beads'
        }
    }
    
    # Merge files
    bd_l <- list()
    for (i in unique(gsub("_.*$",  "",names(dt)))){ #files with same date are merged
        k <- which(gsub("_.*$",  "",names(dt)) %in% i)
        bd_l[[i]] <- dplyr::bind_rows(dt[k])
    }
    
    # Create output directory
    output_path <- file.path(getwd(), 'merged_beadscells')
    if (!dir.exists(output_path)){
      dir.create(output_path, recursive = TRUE)
      message("\nOutput directory created\n")
    }
    
    # Write FCS files
    bd_l %>%
      names(.) %>%
      purrr::map(~ cytofCore::cytofCore.write.FCS(as.matrix(bd_l[[.]]),
        file = paste0(output_path, '/', ., '_bc.FCS'), 
        what = "numeric", channelDescriptions = desc_new))

}




#' Beads unsupervised classification
#'
#' This function uses clustering to classify the events into beads or cells. 
#' Beads are expected to be a group of very little dispersion, so the 
#' classification is assessed by looking into the Coefficient of Variation of
#' the rows classified as beads. If this is > \code{CV_max}, the classification 
#' failed and the \code{df} is converted to 0. Otherwise, \code{df} with a new
#' class column is outputted. The summary stats are written in CSV files either
#' case.
#'
#' @param df An object of class \code{data.frame}.
#' @param method A character vector with the name of the clustering algorithm 
#'        to be used. Options are Gaussian Mixture Models 'GMM' or K-means
#'        'k_means'.
#' @param beads_ch A numeric vector with the column indexes for the beads channels.
#' @param n_clusters A numeric vector indicating the number of clusters to be targeted.
#' @param CV_max A numeric vector indicating the maximum value for the coefficient 
#'        of variation (CV) to evaluate beads cluster. 
#' @param filename A character vector with the name of the file that originated
#'        the data.frame.
#' @param output_path A character vector with full path names.
#' @param FCSdesc FCS column descriptions (extracted from \code{flowFrame} object)
#'
#' @return If CV< \code{CV_max} returns a \code{data.frame} with new class column
#'         for beads or cells. Otherwise, returns a vector = 0. 
BeadsUnsup <- function(df, method = c('GMM', 'k_means'), beads_ch = beads_ch,
    n_clusters = 2, CV_max = 0.05, filename = filename, 
    output_path = output_path, FCSdesc = FCSdesc){
  
    method <- match.arg(method)
  
    df_beads <- t_asinh(df[,beads_ch])
    if (method == 'GMM'){
        cl_data <- mclust::Mclust(df_beads, G = n_clusters)$classification
        df_beads <- as.data.frame(cbind(df_beads, class = cl_data))
    }
    if (method == 'k_means'){
        cl_data <- stats::kmeans(df_beads, centers = n_clusters)$cluster
        df_beads <- as.data.frame(cbind(df_beads, class = cl_data))
    } 

    # Generate summary statistics
    summary_c <- psych::describeBy(df_beads[, 1:length(beads_ch)],
                               group  = df_beads$class,
                               digits = 6)

    # Determine to which group (1 or 2) the beads belong. If not '2', force it
    # Beads = 1 , cells = 0
    if (sum(summary_c[[2]]$median) < sum(summary_c[[1]]$median)){
        k_0 <- c(which(df_beads$class == 1))
        k_1 <- c(which(df_beads$class == 2))
        df_beads$class[k_0] <- 1
        df_beads$class[k_1] <- 0
        summary_c <- psych::describeBy(df_beads[, 1:length(beads_ch)],
                                 group  = df_beads$class,
                                 digits = 6)
    } else {
        k_0 <- c(which(df_beads$class == 1))
        k_1 <- c(which(df_beads$class == 2))
        df_beads$class[k_0] <- 0
        df_beads$class[k_1] <- 1    
    }

    # Write FCS
    df['BeadsSmp_ID'] <- df_beads$class
    cytofCore::cytofCore.write.FCS(as.matrix(df), 
        file = file.path(output_path, paste0(filename, '_bd.FCS')),
        what = "numeric", channelDescriptions = FCSdesc)    

    # Create output folder
    output_path <- file.path(output_path, 'stats')
    if (!dir.exists(output_path)){
      dir.create(output_path, recursive = TRUE)
      message("\nOutput directory created\n")
    }


    # computing coefficient of variation (CV) https://www.westgard.com/lesson34.htm
    if(sum(summary_c[[2]]$sd/summary_c[[2]]$mean) < CV_max*length(beads_ch)){
        print(filename)
        print('PASSED CV checkup')
        
        # Write csv file with stats
        write.csv(summary_c[[2]], file = file.path(output_path, 
            paste0(filename, '_beadsstats.csv')))
      
    } else {
        print(filename)
        print('Clustering FAILED to detect beads.')
        df <- 0

        # Write csv file with stats
        write.csv(summary_c[[2]], file = file.path(output_path, 
            paste0(filename, '_FAILEDbeadsstats.csv')))
    }

    return(df)
}




#' Wrapper function to generate datasets for model training
#'
#' This function generates training and test data sets for bead classification
#' model training. First it uses \code{BeadsUnsup} to select for 'good' samples.
#' Then it uses \code{BalancedSample} and \code{TrainTest} to generate a
#' class-balanced training and test sets.
#'
#' @param sample_size A numeric vector. Number of samples (files) to be used.
#' @param method A character vector. Clustering method for \code{BeadsUnsup}.
#' @param bsample A numeric vector. Size of the sample from each class for 
#'        to be passed to \code{BalancedSample} function.
#' @param class_col A character vector with to label the column that identify 
#'        the classes.
Beads_TrainTest <- function(sample_size = 30, method = method, 
    bsample = 5000, class_col = 'BeadsSmp_ID'){

    # Read files list
    files_list <- list.files(pattern='.FCS|.fcs')
    print('File list successfully created!')

    # Create output directory
    output_path <- file.path(getwd(), 'beads_samplesDS')
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
    desc_new <- c(desc_new, class_col)
    names(desc_new) <- c(names_desc, paste0('$P', (length_desc+1), 'S'))


    # Ask for beads channels  
    print(cn)
    prompt <- "Enter the column INDICES of the beads channels (separated by 
               single space only, no comas allowed) \n"
    beads_ch <- as.numeric(strsplit(readline(prompt), " ")[[1]]) 
    
    # Extract sample of files list
    files_list <- sample(files_list, size = sample_size)

    # Process files iteratively
    dt_ls <- list()
    for (i in 1:length(files_list)){
        # Read exprs from FCS file
        dt <- data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)

        # Order columns
        dt <- dt[col_nms]

        # Classify beads
        filename <- tools::file_path_sans_ext(basename(files_list[i]))
        dt <- BeadsUnsup(dt, method = method, beads_ch = beads_ch, 
            n_clusters = 2, CV_max = 0.05, filename = filename, 
            output_path = output_path, FCSdesc = desc_new)

        # Get balanced sample
        if (is.data.frame(dt)){
            dt <- BalancedSample(dt, sample_size = bsample, class_col = class_col) #function from basic_utils
        }
        dt_ls[[i]] <- dt
    }

    # Add names to list
    names(dt_ls) <- tools::file_path_sans_ext(basename(files_list)) 

    # Remove dt != df from dt_ls
    dt_ls <- Filter(is.data.frame, dt_ls)

    # Compute train/test sets
    TrainTest(dt_ls, output_path = output_path, label = 'beads')
    
}