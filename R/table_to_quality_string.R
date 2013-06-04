#' Convert a quality string into a table
#' 
#' \code{table_to_quality_string} takes a quality string and returns a table
#' with the number of mismatches in the rows and the Entrez IDs or
#' the Refseq IDs, depending on how the quality string was obtained,
#' in the columns.
#' @param x: character. The quality string
#' @return a table with the number of mismatches in the rows and the
#' Entrez/Refseq IDs in the columns
#' @examples
#' data("psets_and_qstrings")
#' s <- psets_and_qstrings$quality_string[1]
#' quality_string_to_table(s)
#' @author Giovanni d'Ario
#' @export
table_to_quality_string <- function(x) {
    mm.nmm <- apply(x, 2, function(z) 
        paste(rownames(x), 
              z, sep = ":", 
              collapse = ";"))
    refs.mm.nmm <- paste(paste(colnames(x), mm.nmm, 
                               sep = "/"), 
                         collapse = "//")
    return(refs.mm.nmm)
}
