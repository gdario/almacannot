##' Check if the probe set is "perfect"
##'
##' A "perfect" probe set is a probe set where all the 11 probes are
##' unequivocally mapping the same transcript, on the positive 
##' strand and with no mismatches. This function returns TRUE if 
##' the quality string of a probe set is such that the probe set 
##' turns out to be "perfect"
##' @title Check if the quality string tells that the probe set is
##' "perfect"
##' @param s a quality string
##' @param n_probes integer. The number of probes composing a probeset.
##' \code{n_probes} is 11 in the old ADX platform and 12 in the new 
##' Excel platform.
##' @return A logical value. If TRUE the probe set is "perfect".
##' @author Giovanni d'Ario
##' @export
is_perfect <- function(s, n_probes=NULL) {
    if(is.null(n_probes))
        stop("Please specify the number of probes composing a probeset")
    x <- quality_string_to_table(s)
    idx1 <- (nrow(x) == 1) & (ncol(x) == 1) # maps one single gene?
    idx2 <- x[1,1] == n_probes
    return(idx1 & idx2)
}

##' Check if all the probe sets associated with a particular
##' combination of Gene ID and Refseq ID are perfect
##'
##' @title Check if all the probe sets associated with a particular
##' combination of Gene ID and Refseq ID are perfect
##' @param gene_id a string of the format geneid:refseqid
##' @param annot an annotation data frame, obtained by applying
##' @param n_probes integer, the number of probes composing a probeset.
##' \code{n_probes} is 11 in the old ADX platform and 12 in the new 
##' Excel platform.
##' \code{almac_reannotation} to the output of bowtie
##' @return a logical value: TRUE if all the probe sets associate with
##' the gene_id are 'perfect'.
##' @author Giovanni d'Ario
are_all_psets_perfect <- function(gene_id=NULL,
                               annot=NULL,
                               n_probes=NULL) {
    if(is.null(n_probes))
        stop("Please specify the number of probes composing a probeset")
    
    idx <- grep(gene_id, annot$quality_string)
    tmp <- annot[idx, ]
    out <- sapply(tmp$quality_string, is_perfect, n_probes = n_probes)
    return(all(out))
}
