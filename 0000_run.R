
  run_parallel <- TRUE

  library(magrittr)

  source("0010_setup.R")

  if(!exists("run_from")) run_from <- 20
  if(!exists("run_to")) run_to <- 100

  do_run <- function(x, y) {

    ls_size <<- x

    current <<- y

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

  sizes <- data.frame(size = c(200, 400, 800, 1600, 3200))

  aois <- data.frame(aoi = c("../envEco/out/KI_50_current"
                             , "../envEco/out/HF_50_current"
                             )
                     ) %>%
    dplyr::mutate(aoi = forcats::fct_inorder(aoi))

  runs <- sizes %>%
    dplyr::left_join(aois, by = character()) %>%
    dplyr::arrange(aoi)


  #-------parallel------

  if(run_parallel) {

    library(furrr)

    # Cores to use for any parallel processing
    use_cores <- if(nrow(runs) > 8) 8 else nrow(runs)

    # Plan for any furrr functions
    future::plan(sequential)

    future::plan(multisession
                 , workers = use_cores
                 )

  }



  #-----RUN-------

  furrr::future_walk2(runs$size
                      , runs$aoi
                      , do_run
                      )


  #------Compare results--------

  # model fit

  sr_fits <- fs::dir_ls(recurse = TRUE
                        , regexp = "fit"
                        ) %>%
    grep("old", ., value = TRUE, invert = TRUE) %>%
    tibble::enframe(name = NULL
                    , value = "path"
                    ) %>%
    dplyr::mutate(agg_size = as.numeric(basename(dirname(path)))) %>%
    dplyr::filter(agg_size %in% sizes$size) %>%
    dplyr::arrange(agg_size) %>%
    dplyr::mutate(obj = purrr::map(path, rio::import)
                  , mod_plot = purrr::map(obj, "mod_plot")
                  , resid_plot = purrr::map(obj, "resid_plot")
                  , imp_plot = purrr::map(obj, "imp_plot")
                  )

  mod_plots <- patchwork::wrap_plots(purrr::map2(sr_fits$mod_plot
                                                 , sr_fits$agg_size
                                                 , ~ .x + labs(title = .y)
                                                 )
                                     )

  resid_plots <- patchwork::wrap_plots(purrr::map2(sr_fits$resid_plot
                                                 , sr_fits$agg_size
                                                 , ~ .x + labs(title = .y)
                                                 )
                                     )

  imp_plots <- patchwork::wrap_plots(purrr::map2(sr_fits$imp_plot
                                                 , sr_fits$agg_size
                                                 , ~ .x + labs(title = .y)
                                                 )
                                     )




  # model results

  sr_bii_tifs <- fs::dir_ls(recurse = TRUE
                            , regexp = "sr_bii"
                            ) %>%
    grep("old", ., value = TRUE, invert = TRUE) %>%
    tibble::enframe(name = NULL
                    , value = "path"
                    ) %>%
    dplyr::mutate(agg_size = as.numeric(basename(dirname(path)))) %>%
    dplyr::filter(agg_size %in% sizes$size) %>%
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
    dplyr::summarise(cells = dplyr::n()
                    , sr_bii = mean(lyr1, na.rm = TRUE)
                    , sr_bii_se = sd(lyr1, na.rm = TRUE) / sqrt(cells)
                    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cells > 100) %>%
    dplyr::left_join(ibra_sub %>%
                       sf::st_set_geometry(NULL) %>%
                       dplyr::mutate(ID = dplyr::row_number()) %>%
                       dplyr::select(ID, contains("_N"))
                     )


  ggplot(sr_bii_tifs, aes(sr_bii, as.factor(agg_size))) +
    geom_col() +
    geom_errorbarh(aes(xmin = sr_bii - sr_bii_se, xmax = sr_bii + sr_bii_se)) +
    facet_wrap(~IBRA_SUB_N)
