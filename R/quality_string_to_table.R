##' Transform a qualityString into a mismatch-by-GeneID matrix 
##'
##' This function takes a qualityString and returns a matrix where
##' the columns names are the mapped Gene IDs and the row names are
##' the number
##' of mismatches. The entries in the table are the number of probes
##' falling in each combination.
##' @title Transform a qualityString into a mismatch-by-GeneID matrix.
##' @param s a qualityString
##' @return a matrix where the column names are the Gene IDs, the
##' rownames are the number of mismatches and the entries are the
##' number of probes falling into each combination.
##' @author Giovanni d'Ario
##' @export
quality_string_to_table <- function(s) {
    tmp <- unlist(strsplit(s, split = "//"))
    tmp <- strsplit(tmp, split = "/")
    gids <- sapply(tmp, function(x) return(x[1]))
    
    rn <- get_row_names(tmp)
    counts <- sapply(tmp, get_counts)
    ## These are characters, and we want numbers
    counts <- matrix(as.numeric(counts), ncol = length(gids))
    colnames(counts) <- gids
    rownames(counts) <- rn
    return(counts)
}

#' Converts a probe set to a quality table
#' 
#' \code{probeset_to_table} takes a probe set ID and an annotation
#' object (containing the quality string information) and returns 
#' the relative quality table.
#' @param p character. A probe set ID.
#' @param annot data frame. The annotation table containing (at 
#' least) the probe set ID and the quality string.
#' @param pcol character. The name of the column containing the probe
#' set IDs. Default is \code{probeset_id}.
#' @param qcol character. The name of the column containing the quality 
#' string. Default is \code{quality_string}.
#' @examples
#' data(psets_and_qstrings)
#' p <- "200003_s_at"
#' probeset_to_table(p = p, annot = psets_and_qstrings)
#' @return a quality table
#' @author Giovanni d'Ario
#' @export
probeset_to_table <- function(p=NULL, 
                              annot=NULL, 
                              pcol="probeset_id", 
                              qcol="quality_string") {
    if(is.null(p))
        stop("I need a probeset ID to proceed")
    if(is.null(annot))
        stop("I need an annotation data frame to proceed")
    
    idx_pset <- match(p, annot[[pcol]])
    s <- annot[[qcol]][idx_pset]
    
    quality_string_to_table(s)
}