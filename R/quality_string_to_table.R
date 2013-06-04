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
##' 
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

