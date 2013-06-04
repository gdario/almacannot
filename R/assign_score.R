#' Assign a score to each gene of a quality table
#' 
#' \code{assign_score} takes a table produced by 
#' \code{quality_string_to_table} assigns a numeric score to each
#' column (gene or refseq, depending on how the table has been 
#' created). The details of the score are described in the 
#' details section. The idea behind this function is to provide the
#' user with a quick measure of the qualityt of the probe set.
#' However, when a low score probe set is encountered, we recommend
#' to call \code{quality_string_to_table} in order to have a more
#' complete picture.
#' @details The scoring method used here assign one point to each
#' zero-mismatch probe, -1/2 to each 1-mismatch, -1/3 to 2-mismatch 
#' and -1/4 to 3-mismatch probes. The gene with the highest score 
#' is then selected while the score of all the other matches is 
#' sign reversed. The final score is obtained summing all the partial
#' scores. The rationale behind this choice is probably clearer 
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
#' @return a numeric vector containing the scores for each gene in
#' the quality table
#' @examples 
#' data("data/psets_and_qstrings.RData")
#' s <- probe_set_annotation$quality_string[1]
#' x <- quality_string_to_table(x)
#' @author Giovanni d'Ario <giovanni.dario@gmail.com>
assign_score <- function(x, w=c(1, -1/2, -1/3, -1/4)) {
    
    if(is.character(x))
        x <- quality_string_to_table(x)
    ## The rownames contain (as characters) the number of mismatches
    mismatches <- as.numeric(rownames(x))
    ## Assign the weigths (using only the existing mismatches)
    weights <- w[mismatches + 1]
    ## Compute the scores
    scores  <- weights %*% x
    
    ## Identify the best score and the best gene
    idx_best <- which.max(scores)
    best_match <- colnames(scores)[idx_best]
    
    ## Compute the final score as the best score minus the absolute
    ## value of the other scores
    final_score <- scores[idx_best] - sum(abs(scores[-idx_best]))
    return(list(best_match = best_match, final_score = final_score))
}

