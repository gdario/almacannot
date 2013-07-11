context("Checks on assign_score")

has_max <- function(maxval) {
    function(x) {
        expectation(
        max(x) == maxval,
        paste("The maximum value is not", maxval)
        )
    }
}

has_dim <- function(nr, nc) {
    function(x) {
        expectation(
            (nrow(x) == nr) & (ncol(x) == nc),
            "Dimensions are not the specified ones"
        )        
    }
}

test_that("compute_score works on one row matrices", {
    ## Perfect probe set
    s <- "10594/0:11"
    score <- almacannot:::compute_score(s, w = c(1, -1/2, -1/3, -1/4))
    expect_that(score, is_a("data.frame"))
    expect_that(score, has_dim(nr = 1, nc = 2))
    expect_that(score[[2]], equals(11))
})

test_that("assign_score returns sensible scores", {
    data(psets_and_qstrings)
    new_annot <- assign_score(psets_and_qstrings)
    
    expect_that(new_annot, is_a("data.frame"))
    expect_that(new_annot$final_score, has_max(11))
})
