#' Add the best matching gene and the quality score to an annotation
#' table.
#' 
#' The \code{assign_score} function takes a data frame containing
#' (at least) the probe set ID and the quality scores and adds to
#' it two columns containing the ID of the best matching gene/refseq
#' and the overall quality score.
#' @param x data frame containing (at least) the probe set IDs and
#' the quality string computed by \code{add_quality_string}
#' @param qs_col character, name of the column containing the quality
#' string. Defaults to \code{quality_string}.
#' @examples
#' data(psets_and_qstrings)
#' new_annot <- assign_score(psets_and_qstrings)
#' @return a data frame obtaining from the binding of the original
#' input \code{x} and of the two additional columns with the best 
#' match and the final score.
#' @author Giovanni d'Ario
#' @export
assign_score <- function(x, qs_col="quality_string") {
    quality_strings <- x[[qs_col]]
    best_match_and_score <- sapply(quality_strings,
                                   compute_score,
                                   simplify = FALSE)
    names(best_match_and_score) <- NULL
    best_match_and_score <- do.call(rbind, best_match_and_score)
    if(nrow(best_match_and_score) == nrow(x))
        out <- cbind(x, best_match_and_score)
    return(out)
}

#' Assign a score to each gene of a quality table
#' 
#' \code{compute_score} takes a table produced by 
#' \code{quality_string_to_table} assigns a numeric score to each
#' column (gene or refseq, depending on how the table has been 
#' created). The details of the score are described below.
#' The idea behind this function is to provide the
#' user with a quick measure of the qualityt of the probe set.
#' However, when a low score probe set is encountered, we recommend
#' to call \code{quality_string_to_table} in order to have a more
#' complete picture.
#' The scoring method used here assign one point to each
#' zero-mismatch probe, -1/2 to each 1-mismatch, -1/3 to 2-mismatch 
#' and -1/4 to 3-mismatch probes. The gene with the highest score 
#' is then selected. The final score is obtained subtacting the 
#' absolute value of the second best score from the first one. The 
#' absolute value is due to the fact that some second hits could be
#' only mismatches, and an overall negative score. 
#' The rationale behind this choice is probably clearer 
#' looking at some examples.
#' \enumerate{
#' \item if all the probes are perfectly matching and no mismatch or
#' cross-hybridization problem occurs, the score equals the number
#' of matching probes.
#' \item If two genes are perfectly matched by all the probes, the 
#' final score is zero (it is not possible to say which gene is 
#' being mapped). This would be not possible if, for example, we 
#' penalized the probes by a weighted sum of mismatches, using
#' for example, the number of mismatches as weights. In such a case,
#' two perfectly matching genes would both have zero score.
#' \item No distinction is made between sense and antisense strand.
#' If the best match is on an antisense strand, the reported best
#' match will be the relative Entrez ID preceded by `ANTISENSE`.
#' \item If a probe set has many spurious matches, its score could
#' be negative.
#' }
#' @param x either a string or a table produced by 
#' \code{quality_string_to_table}. If \code{x} is a string, it is 
#' silently converted into a table by calling 
#' \code{quality_string_to_table}.
#' @param w numeric weights to be given to the number of mismatches.
#' @return a data frame with the ID of the best matching gene/refseq
#' and the relative final score.
#' @examples 
#' data("psets_and_qstrings")
#' s <- psets_and_qstrings$quality_string[1]
#' x <- quality_string_to_table(s)
#' @author Giovanni d'Ario <giovanni.dario@gmail.com>
compute_score <- function(x, w=c(1, -1/2, -1/3, -1/4)) {
    
    if(is.character(x))
        x <- quality_string_to_table(x)
    ## The rownames contain (as characters) the number of mismatches
    mismatches <- as.numeric(rownames(x))
    ## Assign the weigths (using only the existing mismatches)
    weights <- w[mismatches + 1]
    ## Compute the scores
    scores  <- weights %*% x
    
    ### order according to the score
    idx <- order(scores, decreasing=TRUE)
    scores <- scores[, idx, drop = FALSE]
    
    ## Identify the best score and the best gene
    best_match <- colnames(scores)[1]
    
    ## Compute the final score as the best score minus the second
    ## best score
    ## N.B.: we take the absolute value of the second score since 
    ## it could be negative. Consider, for example, the case where
    ## the second best match has no perfect match but several probes
    ## with one mismatch.
    if(ncol(scores) > 1) {
        final_score <- scores[1] - abs(scores[2])
    } else {
        final_score <- scores[1]
    }
    return(data.frame(best_match = best_match, 
                      final_score = final_score,
                      stringsAsFactors = FALSE))
}
