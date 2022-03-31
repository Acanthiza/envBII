
  library(magrittr)

  if(!exists("run_from")) run_from <- 0
  if(!exists("run_to")) run_to <- 100

  commit_notes <- "Square-root sr as outcome"

  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= run_to & as.numeric(names(.)) >= if(run_from == 0) 1 else run_from]) %>%
    purrr::walk(source
                , verbose = TRUE
                )
