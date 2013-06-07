#' Generate the annotation table that will be used for the quality 
#' string
#'
#' @title Generate the annotation table that will be used for the 
#' quality string.
#' @param input.file character string containing the full path to
#' the input file (a Bowtie output file preprocessed by 
#' `preprocess_bowtie_output.R`)
#' @param chip character, the type of ALMAC array under
#' consideration. It can be either \code{adxcrc} or \code{xcel}.
#' The default is \code{adxcrc}.
#' @param output.file character string containing the full path to 
#' the output file where the probe sets and their quality string 
#' will be saved. 
#' @param by.refseq logical, indicates whether the number of
#' mismatches should be tabulated against the refseq_ids rather
#' than against the gene_ids. Default is FALSE
#' @param rm.affx logical, should the probe sets starting with
#' \code{AFFX} (usually controls) be removed? Default is TRUE.
#' @param rm.asense logical. Should the probes mapping on the 
#' opposite strand be removed? Default is FALSE, in which case 
#' the gene/refseq ID is preceded by \code{ANTISENSE}.
#' @examples
#' \dontrun{
#' input.file <- "/export/big/Data/almac_annotations/new_annotation_2013/output/formatted_xcel_subset.txt"
#' output.file <- "~/Desktop/pset_and_qscores.RData"
#' add_quality_string(input.file = input.file, chip = "xcel", output.file = output.file)
#' }
#' @return an (invisible) data frame containing the probe set ID and
#' the associated quality string.
#' @author Giovanni d'Ario
#' @export
add_quality_string <- function(
    input.file=NULL,
    chip=c("adxcrc", "xcel"),
    rm.affx=TRUE,
    rm.asense=FALSE,
    by.refseq=FALSE,
    output.file=NULL) {
    
    if(is.null(input.file))
        stop("Please enter the bowtie Data file name")
    if(is.null(output.file))
        stop("Please enter the output file name")
    
    chip <- match.arg(chip)
    
    ## All the columns are characters with the exception of the
    ## number of mismatches (the last column).
    colClasses <- c(rep("character", 
                        ifelse(chip == "adxcrc", 5, 6)),
                        "integer")
    message("Reading the bowtie Data file...")
    Data <- read.delim(file = input.file,
                       header = TRUE,
                       colClasses = colClasses)
    
    ## if rm.asense = TRUE, remove everything that is mapped
    ## to the reverse strand, otherwise, prepend ANTISENSE to the
    ## corresponding gene_id and refseq
    idx_antisense <- Data$strand %in% "-"
    
    if(rm.asense) {
        Data <- Data[!idx_antisense, ]
    } else {
        Data$gene_id[idx_antisense] <- paste("ANTISENSE", 
                                             Data$gene_id[idx_antisense],
                                             sep = "_")
        Data$refseq_id[idx_antisense] <- paste("ANTISENSE", 
                                               Data$refseq_id[idx_antisense],
                                               sep = "_")
    }
    
    ## If rm.affx, remove the control probes that begin by "AFFX"
    if(rm.affx) {
        idx_affx <- grep("AFFX", Data[, 1])
        ## Make sure that idx_affx is longer than zero
        if(length(idx_affx) > 0)
            Data <- Data[-idx_affx, ]
    }
    
    if(by.refseq) {
        ## Create a column that contains both the Gene ID and the Refseq ID
        Data$gene_refseq <- paste(Data$gene_id, 
                                  Data$refseq_id, sep = ":")
        # Data$refseq_id <- Data$gene_id <- NULL
    } else {
        Data$refseq_id <- NULL
        Data <- unique(Data)
    }
    
    ## Sort by probe set and by probe
    message("Sorting by probe set ID")
    idx_sort <- order(Data$probeset_id)
    Data <- Data[idx_sort, ]
    
    ## Split the matrix by probe set
    message("Splitting by probe set...")
    by_pset <- split(Data, Data$probeset_id)
    
    ## Compute the cross tabulation of the number of mismatches
    ## by refseq_id or (default) by gene_id
    message("Cross tabulating mismatches and Gene IDs...")
    if(by.refseq) {
        message("Cross tabulating mismatches and Refseq IDs...")
        mm_table <- lapply(by_pset, function(x)
            xtabs(~ n_mismatches + gene_refseq, data = x))  
    } else {
        mm_table <- lapply(by_pset, function(x)
            xtabs(~ n_mismatches + gene_id, data = x))
    }    
    
    ## Generate the quality strings
    message("Generating the quality strings...")
    
    quality_string <- lapply(mm_table, table_to_quality_string)
    probe_set_annotation <- data.frame(probeset_id = names(quality_string),
                                       quality_string = unlist(quality_string),
                                       stringsAsFactors = FALSE)
    rownames(probe_set_annotation) <- NULL
    
    message(paste("Saving the output in ", 
                  output.file, "...", 
                  sep = "" ))
    save(probe_set_annotation, file = output.file,
         compress = TRUE)
    return(invisible(probe_set_annotation))
}
