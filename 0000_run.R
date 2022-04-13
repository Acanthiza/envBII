
  run_parallel <- 6

  library(magrittr)

  if(!exists("run_from")) run_from <- 0
  if(!exists("run_to")) run_to <- 100

  do_run <- function(x, y, z) {

    ls_size <<- x

    current <<- as.character(y)

    cv_method <<- as.character(z)

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

  sizes <- data.frame(size = c(100, 200, 400, 800, 1600
                               , 3200
                               #, 6400
                               )
                      )

  aois <- data.frame(aoi = c("../envEco/out/KI_50_current"
                             , "../envEco/out/HF_50_current"
                             )
                     ) %>%
    dplyr::mutate(aoi = forcats::fct_inorder(aoi))

  cv_methods <- data.frame(cv_method = c("cv_normal"
                                     , "cv_spatial"
                                     )
                          )

  runs <- sizes %>%
    dplyr::left_join(aois, by = character()) %>%
    dplyr::left_join(cv_methods, by = character()) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(aoi) %>%
    dplyr::mutate(check = fs::path(gsub("envEco", "envBII", aoi)
                                   , cv_method
                                   , size
                                   , "sr_bii.tif"
                                   )
                  , done = file.exists(check)
                  )


  #-------parallel------

  if(!isFALSE(run_parallel)) {

    library(furrr)

    # Cores to use for any parallel processing
    use_cores <- if(nrow(runs[!runs$done,]) == 0) {

        0

      } else if(nrow(runs[!runs$done,]) >= run_parallel) {

        run_parallel

      } else {

        nrow(runs[!runs$done,])

      }

    # Plan for any furrr functions

    future::plan(sequential)

    if(use_cores > 0) {

      future::plan(multisession
                   , workers = use_cores
                   )

    }

  }


  #-----RUN-------

  furrr::future_pwalk(list(runs$size[!runs$done]
                           , runs$aoi[!runs$done]
                           , runs$cv_method[!runs$done]
                           )
                      , do_run
                      )


  #------Compare results--------

  # model fit

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

  sr_fits <- runs %>%
    dplyr::mutate(aoi = gsub(".*\\/|_\\d{1,}.*", "", aoi)
                  , path = gsub("sr_bii\\.tif", "fit.rds", check)
                  , done = file.exists(path)
                  ) %>%
    dplyr::select(size, aoi, cv_method, done, path) %>%
    dplyr::filter(size %in% sizes$size) %>%
    dplyr::arrange(size) %>%
    dplyr::filter(done) %>%
    dplyr::mutate(fit = purrr::map(path, rio::import)
                  #, mod_plot = purrr::map(fit, "mod_plot")
                  #, resid_plot = purrr::map(fit, "resid_plot")
                  , imp_plot = purrr::map(fit, "imp_plot")
                  , wf_results = purrr::map(fit
                                            , function(x) if(is.null(x$wf_results)) x$race_results else x$wf_results
                                            )
                  , func = purrr::map_chr(wf_results
                                          , ~extract_best_met(.
                                                              , "rmse"
                                                              , "model"
                                                              )
                                          )
                  # , pack = purrr::map_chr(fit
                  #                         , ~ .$fit$fit$fit$fit$call %>% `[[`(1) %>% `[[`(2) %>% as.character()
                  #                         )
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
                  , full_model = purrr::map(fit
                                            , ~yardstick::metrics(.$final_fit
                                                                  , truth = sr
                                                                  , estimate = .pred
                                                                  ) %>%
                                              tidyr::pivot_wider(names_from = .metric
                                                                 , values_from = .estimate
                                                                 , names_prefix = "full_model"
                                                                 )
                                            )
                  ) %>%
    tidyr::unnest(cols = c(full_model))

  sr_fits_best <- sr_fits %>%
    dplyr::group_by(aoi, cv_method) %>%
    dplyr::mutate(best = full_modelrmse == min(full_modelrmse)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = paste0(aoi
                                 , ": "
                                 , size
                                 , " ("
                                 , cv_method
                                 , ")"
                                 )
                  )

  sr_fits_short <- sr_fits_best %>%
    dplyr::select(size, aoi, cv_method
                  , contains("full")
                  , rmse, rsq
                  , func
                  #, pack
                  , best, label
                  )

  library(ggplot2)
  library(patchwork)

  mod_plots <- sr_fits_best %>%
    dplyr::mutate(final_fit = purrr::map(fit, "final_fit")) %>%
    dplyr::select(any_of(names(sr_fits_short)), final_fit) %>%
    tidyr::unnest(cols = c(final_fit)) %>%
    ggplot(aes(sr
               , .pred^2
               , colour = best
               )
           ) +
      geom_jitter(shape = ".") +
      geom_abline(color = "gray50", lty = 2) +
      geom_smooth() +
      facet_grid(size ~ paste0(aoi, ": ", cv_method))


  resid_plots <- sr_fits_best %>%
    dplyr::mutate(final_fit = purrr::map(fit, "final_fit")) %>%
    dplyr::select(any_of(names(sr_fits_short)), final_fit) %>%
    tidyr::unnest(cols = c(final_fit)) %>%
    ggplot(aes(sr, resid, colour = best)) +
      geom_jitter(shape = ".") +
      geom_hline(aes(yintercept = 0), color = "gray50", lty = 2) +
      geom_smooth() +
      facet_grid(size ~ paste0(aoi, ": ", cv_method))

  imp_plots <- patchwork::wrap_plots(purrr::map2(sr_fits_best$imp_plot[sr_fits_best$best]
                                                 , sr_fits_best$label[sr_fits_best$best]
                                                 , ~ .x + labs(title = .y)
                                                 )
                                     )


  #-------maps--------

  # IBRA Sub
  ibra_sub <- sf::st_read(fs::path("..","envEco", "out","shp","ibra_sub.shp")) %>%
    sf::st_transform(crs = 4283) %>%
    sf::st_make_valid()

  # LSA
  lsa <- sf::st_read(fs::path("..","envEco", "out","shp","lsa.shp")) %>%
    sf::st_transform(crs = 4283) %>%
    sf::st_make_valid()


  #--------model results--------

  sr_bii_tifs <- runs %>%
    dplyr::mutate(aoi = gsub(".*\\/|_\\d{1,}.*", "", aoi)
                  , path = check
                  , done = file.exists(path)
                  ) %>%
    dplyr::select(size, aoi, cv_method, done, path) %>%
    dplyr::filter(size %in% sizes$size) %>%
    dplyr::arrange(size) %>%
    dplyr::filter(done) %>%
    dplyr::mutate(r = purrr::map(path, terra::rast)
                  , tif_cells = purrr::map_dbl(r, terra::ncell)
                  , sr_bii_extract = purrr::map(r
                                                , terra::extract
                                                , y = terra::vect(lsa)
                                                )
                  ) %>%
    dplyr::left_join(sr_fits_short) %>%
    tidyr::unnest(cols = sr_bii_extract) %>%
    dplyr::left_join(lsa %>%
                       sf::st_set_geometry(NULL) %>%
                       dplyr::mutate(ID = dplyr::row_number()) %>%
                       dplyr::select(ID, where(is.character))
                     )

  sr_bii_tifs_summary <- sr_bii_tifs %>%
    dplyr::group_by(ID, tif_cells, across(any_of(names(sr_fits_short))), LSA, path) %>%
    dplyr::summarise(cells = dplyr::n()
                    , sr_bii = mean(lyr1, na.rm = TRUE)
                    , sr_bii_se = sd(lyr1, na.rm = TRUE) / sqrt(cells)
                    ) %>%
    dplyr::ungroup()


  ggplot(sr_bii_tifs_summary %>%
           dplyr::filter(aoi == LSA)
         , aes(sr_bii
               , as.factor(size)
               , fill = best
               #, colour = best
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
    facet_grid(cv_method ~ LSA)


  #------outputs-------

  final_tifs <- sr_bii_tifs_summary %>%
    dplyr::group_by(aoi) %>%
    dplyr::filter(aoi == LSA
                  , best
                  , rmse == min(rmse)
                  ) %>%
    dplyr::ungroup()

  hf_tif <- final_tifs %>%
    dplyr::filter(aoi == "HF") %>%
    dplyr::pull(path) %>%
    stars::read_stars()

  ki_tif <- final_tifs %>%
    dplyr::filter(aoi == "KI") %>%
    dplyr::pull(path) %>%
    stars::read_stars()

  #-----visualise-------

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



  hf_map <- plot_ras(hf_tif, 8, 0)

  hf_map
  hist(hf_tif)


  ki_map <- plot_ras(ki_tif, 8, 0)

  ki_map
  hist(ki_tif)


