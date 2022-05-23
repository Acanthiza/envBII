

#--------model results--------

  extract_best_met <- function(race_results
                               , metric = "rmse"
                               , col = "mean"
                               ) {

    bm <- race_results %>%
      workflowsets::rank_results() %>%
      dplyr::filter(.metric == metric) %>%
      dplyr::slice(1) %>%
      dplyr::pull(!!rlang::ensym(col))

    bm

  }

  fit_info <- tibble(fit = list(fit)) %>%
    dplyr::mutate(imp_plot = purrr::map(fit, "imp_plot")
                  , wf_results = purrr::map(fit, "wf_results")
                  , func = purrr::map_chr(wf_results
                                          , ~extract_best_met(.
                                                              , "rmse"
                                                              , "model"
                                                              )
                                          )
                  , rmse = purrr::map_dbl(wf_results
                                          , ~extract_best_met(.
                                                              , "rmse"
                                                              )
                                          )
                  , rsq = purrr::map_dbl(wf_results
                                         , ~extract_best_met(.
                                                             , "rsq"
                                                             )
                                         )
                  , test = purrr::map(fit
                                            , ~yardstick::metrics(.$test_res
                                                                  , truth = sr
                                                                  , estimate = .pred
                                                                  ) %>%
                                              tidyr::pivot_wider(names_from = .metric
                                                                 , values_from = .estimate
                                                                 , names_prefix = "test_"
                                                                 )
                                            )
                  ) %>%
    tidyr::unnest(cols = c(test)) %>%
    dplyr::select(contains("test_")
                  , rmse, rsq
                  , func
                  )

  sr_bii_res <- tibble(r = list(sr_bii)) %>%
    dplyr::mutate(tif_cells = purrr::map_dbl(r, terra::ncell)
                  , sr_bii_extract = purrr::map(r
                                                , terra::extract
                                                , y = terra::vect(lsa)
                                                )
                  ) %>%
    tidyr::unnest(cols = sr_bii_extract) %>%
    dplyr::left_join(lsa %>%
                       sf::st_set_geometry(NULL) %>%
                       dplyr::mutate(ID = dplyr::row_number()) %>%
                       dplyr::select(ID, where(is.character))
                     )

  sr_bii_res_summary <- sr_bii_res %>%
    dplyr::group_by(ID, tif_cells, across(any_of(names(fit_info))), LSA) %>%
    dplyr::summarise(cells = dplyr::n()
                     , sr_bii = mean(lyr1, na.rm = TRUE)
                     , sr_bii_se = sd(lyr1, na.rm = TRUE) / sqrt(cells)
                     ) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(sr_bii_res %>%
                       dplyr::summarise(cells = dplyr::n()
                                        , sr_bii = mean(lyr1, na.rm = TRUE)
                                        , sr_bii_se = sd(lyr1, na.rm = TRUE) / sqrt(cells)
                                        ) %>%
                       dplyr::mutate(LSA = "State")
                     )


  sw_sr_res <- ggplot(sr_bii_res_summary
         , aes(sr_bii
               , LSA
               , fill = LSA
               )
         ) +
    geom_col() +
    geom_errorbarh(aes(xmin = sr_bii - sr_bii_se, xmax = sr_bii + sr_bii_se)
                   , height = 0.2
                   ) +
    geom_vline(aes(xintercept = 0)
               , linetype = 2
               , colour = "red"
               ) +
    lsa_pal_fill



  plot_ras <- function(star_ras, bks = 5, mp = 0) {

    # bks - classes in continuous legend = bks + 2
    # mp - midpoint for continuous legend.

    n_cells <- stars::st_dimensions(star_ras)$x$to * stars::st_dimensions(star_ras)$y$to

    sr_cuts <- unique(c(-1, hist(star_ras, bks, plot = FALSE)$breaks, 1))

    tmap::tm_shape(star_ras
                   , raster.downsample = n_cells > 42468400/4
    ) +
      tmap::tm_raster(breaks = sr_cuts
                      , palette = "viridis"
                      , midpoint = mp
      )

  }

  tif <- terra::sources(sr_bii) %>%
    stars::read_stars()

  map_tif <- plot_ras(tif, 8, 0)
