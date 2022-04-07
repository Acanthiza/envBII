
  library(magrittr)

  if(!exists("run_from")) run_from <- 0
  if(!exists("run_to")) run_to <- 100

  do_run <- function(x) {

    ls_size <<- x

    commit_notes <<- paste0("Aggregation landscape at "
                           , ls_size
                           , " m"
                           )

    dir() %>%
      grep("^\\d{4}_.*\\.R$",.,value=TRUE) %>%
      setNames(stringr::str_extract(.,"\\d{4}")) %>%
      `[` (names(.)[as.numeric(names(.)) <= run_to & as.numeric(names(.)) >= if(run_from == 0) 1 else run_from]) %>%
      purrr::walk(source
                  , verbose = TRUE
                  )

  }

  sizes <- c(200, 400, 800, 1600, 3200)

  purrr::map(sizes
             , ~do_run(.)
             )


  #------Compare results--------

  sr_bii_tifs <- fs::dir_ls(recurse = TRUE
                            , regexp = "sr_bii"
                            ) %>%
    tibble::enframe(name = NULL
                    , value = "path"
                    ) %>%
    dplyr::mutate(agg_size = as.numeric(basename(dirname(path)))) %>%
    dplyr::filter(agg_size %in% sizes) %>%
    dplyr::arrange(agg_size) %>%
    dplyr::mutate(r = purrr::map(path, terra::rast)
                  , cells = purrr::map_dbl(r, terra::ncell)
                  , sr_bii_extract = purrr::map(r
                                                , terra::extract
                                                , y = terra::vect(ibra_sub)
                                                )
                  ) %>%
    tidyr::unnest(cols = sr_bii_extract) %>%
    dplyr::group_by(ID, agg_size) %>%
    dplyr::summarise(cells = n()
                    , sr_bii = mean(lyr1, na.rm = TRUE)
                    , sr_bii_se = sd(lyr1) / sqrt(cells)
                    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cells > 100) %>%
    dplyr::left_join(ibra_sub %>%
                       sf::st_set_geometry(NULL) %>%
                       dplyr::mutate(ID = row_number()) %>%
                       dplyr::select(ID, contains("_N"))
                     )


  ggplot(sr_bii_tifs, aes(sr_bii, as.factor(agg_size))) +
    geom_col() +
    geom_errorbarh(aes(xmin = sr_bii - sr_bii_se, xmax = sr_bii + sr_bii_se)) +
    facet_wrap(~IBRA_SUB_N)
