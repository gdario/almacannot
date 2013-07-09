context("Test the is_perfect function")

test_that("is_perfect identifies perfect probe sets", {
    s <- "10594/0:11"
    expect_that(is_perfect(s, n_probes=11), is_true())
    s <- "825/0:0;1:0;3:1//826/0:10;1:1;3:0//84290/0:0;1:0;3:1"
    expect_that(is_perfect(s, 11), is_false())
})