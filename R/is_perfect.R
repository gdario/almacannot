##' Check if the probe set is "perfect"
##'
##' A "perfect" probe set is a probe set where all the 11 probes are
##' unequivocally mapping the same transcript, on the positive 
##' strand and with no mismatches. This function returns TRUE if 
##' the quality string of a probe set is such that the probe set 
##' turns out to be "perfect"
##' @title Check if the quality string tells that the probe set is
##' "perfect"
##' @param a quality string
##' @return A logical value. If TRUE the probe set is "perfect".
##' @author Giovanni d'Ario
##' @export
is_perfect <- function(s) {
    x <- quality_string_to_table(s)
    idx1 <- (nrow(x) == 1) & (ncol(x) == 1) # maps one single gene?
    idx2 <- x[1,1] == 11
    return(idx1 & idx2)
}

##' Check if all the probe sets associated with a particular
##' combination of Gene ID and Refseq ID are perfect
##'
##' @title Check if all the probe sets associated with a particular
##' combination of Gene ID and Refseq ID are perfect
##' @param generef a string of the format geneid:refseqid
##' @param annot an annotation data frame, obtained by applying
##' \code{almac_reannotation} to the output of bowtie
##' @return a logical value: TRUE if all the probe sets associate with
##' the generef are 'perfect'.
##' @author Giovanni d'Ario
all_perfect_probes <- function(generef=NULL,
                               annot=NULL) {
    idx <- grep(generef, annot$quality_string)
    tmp <- annot[idx, ]
    out <- sapply(tmp$quality_string, is_perfect)
    return(all(out))
}
