
  library(magrittr)

  if(!exists("run_from")) run_from <- 0
  if(!exists("run_to")) run_to <- 40

  # Overall settings
  if(!exists("max_cores")) max_cores <- 14

  if(!exists("ls_size")) ls_size <- 200

  if(!exists("current")) current <- "../envEco/out/SA_50_current"

  if(!exists("cv_method")) cv_method <- "cv_normal"

  dir() %>%
    grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
    setNames(stringr::str_extract(.,"\\d{4}")) %>%
    `[` (names(.)[as.numeric(names(.)) <= run_to & as.numeric(names(.)) >= if(run_from == 0) 1 else run_from]) %>%
    purrr::walk(source
                , verbose = TRUE
                )


