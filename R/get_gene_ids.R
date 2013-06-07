#' Extract the Gene Ids from a quality string
#' 
#' Accessor function to extract the Entrez/Refseq gene IDs from a 
#' quality string. This function is not supposed to be used by the
#' end user.
#' @param s character. A quality string.
#' @param simplify logical. Should the output be simplified?
#' @author Giovanni d'Ario
get_gene_ids <- function(s, simplify=FALSE) {
    tmp <- unlist(strsplit(s, split = "//"))
    tmp <- strsplit(tmp, split = "/")
    gidrefs <- sapply(tmp, function(x) return(x[1]))
    if(!simplify) {
        gidrefs <- sapply(gidrefs, strsplit, split = ":")
        gidrefs <- t(as.data.frame(gidrefs))
        colnames(gidrefs) <- c("gene_id", "refseq_id")
        rownames(gidrefs) <- NULL
    }
    return(gidrefs)
}


get_row_names <- function(tmp) {
    x <- tmp[[1]]
    mms <- unlist(strsplit(x[2], split = ";"))
    mms <- strsplit(mms, split = ":")
    ## Rownames, i.e. the number of mismatches
    rn <- sapply(mms, function(x) return(x[1]))
    rn
}

get_counts <- function(x) {
    tmp1 <- strsplit(x, split = ";")
    tmp2 <- strsplit(tmp1[[2]], split = ":")
    cnt <- sapply(tmp2, function(x) return(x[2]))
    cnt
}
