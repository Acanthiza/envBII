
  tests <- c("test_01.R", "test_02.R")

  test_fun <- function(x) {

    purrr::map(tests
               , source
               )

  }

  ls_size <- c(1, 3, 40)

  purrr::map(ls_size
             , ~test_fun(x = .)
             )
