
#' Arcsinh transformation 
#'
#' This function takes a data.frame and applies basic \code{asinh} function 
#' using a cofactor.
#'
#' @param df An object of class \code{data.frame} containing columns to be 
#'        scaled. The function ignores non-numeric columns.
#' @param cofactor The cofactor for the scaling. This could be any number >0:Inf.
#'
#' @return Returns a \code{data.frame} with the scaled data.
t_asinh <- function(df, cofactor = 5) {
    nums <- vapply(df, is.numeric, FUN.VALUE = logical(1)) # only operates in numeric
    df[,nums] <- asinh(df[,nums]/cofactor)
    (df)
}



#' Extract data (fcs, csv or txt)
#'
#' This function reads files of current folder into a \code{list} of 
#' \code{data.frames}
#'
#' @param file_type A regular expression for the file extension. Only file names 
#'        which match the regular expression will be returned. Data types suported
#'        are .FCS, .CSV and .TXT
#' @param sampling Logical. Should a sample of size \code{sample_size} extracted 
#'        from the list of files? If \code{FALSE} the function will extract all 
#'        files of the list. Default is \code{FALSE}
#' @param sample_size Numeric. If \code{sampling=TRUE} this would be the sample 
#'        size.
#' @param SpecifyList Logical. If \code{TRUE}, a vector of file names shold be 
#'        given in \code{f_list}. Default is \code{FALSE}
#' @param f_list If \code{SpecifyList=TRUE}, this should be a vector of file names.
#'
#' @return Returns a \code{list} of \code{data.frames} of the extracted data.
extract_data <- function (file_type = '.fcs|.FCS',
                          sampling = FALSE, sample_size = 10, 
                          SpecifyList = FALSE, f_list) {
  fcs <- c('.FCS', '.fcs', '.fcs|.FCS')
  csv <- c('.CSV', '.csv', '.csv|.CSV')
  txt <- c('.TXT', '.txt', '.txt|.TXT')
  
  if (!(file_type %in% c(fcs,csv,txt))){
    stop("Data type not supported!")
  }

  if (SpecifyList){
    files_list <- f_list
  } else {
    files_list <- list.files(pattern=file_type)
  }
  
  if (sampling == T){
        files_list <- sample(files_list, size = sample_size)
  }
  exprs_data <- list()
  if (file_type %in% fcs){
    for (i in 1:length(files_list)){
      exprs_data[[i]] <-as.data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)
    }}
  if (file_type %in% csv){
    for (i in 1:length(files_list)){
      exprs_data[[i]] <-as.data.frame(read.csv(files_list[i]))
    }}
  if (file_type %in% txt){
    for (i in 1:length(files_list)){
      exprs_data[[i]] <-as.data.frame(read.csv(files_list[i], sep = '\t'))
    }}
  
  names(exprs_data) <- tools::file_path_sans_ext(basename(files_list))
  return(exprs_data)
}