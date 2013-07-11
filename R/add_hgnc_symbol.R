#' @title Add the HGNC gene symbols to an annotation table
#' @description This function adds the HGNC gene symbol to an annotation
#' table obtained by applying \code{assign_score} to the output of a 
#' Bowtie alignment.
#' @details The \code{add_hgnc_symbol} function takes an annotation table 
#' produced by \code{assign_score} and adds two columns, one with 
#' the HGNC gene symbol, and one specifying whether the probe is 
#' mapping the forward or the reverse strand of the gene.
#' @param annot data frame containing the probeset ID, the quality 
#' string and the best match. It should be the output of 
#' \code{assign_score} (thus it should also include the final score).
#' @param best_match_col character, the name of the column containing
#' the best matching gene. Default is \code{best_match}.
#' @examples
#' data(psets_and_qstrings)
#' annotation_with_score <- assign_score(x=psets_and_qstrings)
#' annotation_with_hgnc <- add_hgnc_symbol(annot=annotation_with_score)
#' @return a data frame containing the annotation table with the 
#' final score, the best matches, the HGNC symbol and the strand 
#' information.
#' @seealso assign_score
#' @author Giovanni d'Ario
#' @export
add_hgnc_symbol <- function(annot,
                           best_match_col="best_match") {
    if(!require(AnnotationDbi))
        stop("AnnotationDbi is required by this function")
    if(!require(org.Hs.eg.db))
        stop("The 'org.Hs.eg.db' package does not seem to be installed")
    best_match <- as.character(annot[[best_match_col]])
    
    # Identify the probesets mapping the reverse strand
    idx_asense <- grep("ANTISENSE", best_match)
    
    # Remove the `ANTISENSE` tag from the ENTREZ IDs
    best_match <- sub("ANTISENSE_", "", best_match)
    
    # Extract the HGNC from the ENTREZ IDs
    hgnc_symbol <- AnnotationDbi::mget(best_match, 
                                       org.Hs.egSYMBOL,
                                       ifnotfound = NA)
    
    # In case of multiple HGNC symbols, return only the first one
    hgnc_symbol <- sapply(hgnc_symbol, function(x) return(x[1]))
    # Check that the length of the HGNC symbols have the same length
    # as the elements in the annotation table
    if(length(hgnc_symbol) != nrow(annot))
        stop("Matching error: debug add_hgnc_score")
    
    # Add the information on the strand
    strand <- rep("forward", nrow(annot))
    strand[idx_asense] <- "reverse"
    
    annot$hgnc_symbol <- hgnc_symbol
    annot$strand  <- strand

    return(annot)
}