
  run_parallel <- 14

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

  sizes <- data.frame(size = c(200, 400, 800, 1600, 3200, 6400, 12800))

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
      dplyr::pull(!!ensym(col))

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
                  , mod_plot = purrr::map(fit, "mod_plot")
                  , resid_plot = purrr::map(fit, "resid_plot")
                  , imp_plot = purrr::map(fit, "imp_plot")
                  , func = purrr::map_chr(fit
                                          , ~extract_best_met(.$race_results
                                                              , "rmse"
                                                              , "model"
                                                              )
                                          )
                  , pack = purrr::map_chr(fit
                                          , ~ .$fit$fit$fit$fit$call %>% `[[`(1) %>% `[[`(2) %>% as.character()
                                          )
                  , rmse = purrr::map_dbl(fit
                                           , ~extract_best_met(.$race_results
                                                               , "rmse"
                                                               )
                                           )
                  , rsq = purrr::map_dbl(fit
                                          , ~extract_best_met(.$race_results
                                                              , "rsq"
                                                              )
                                          )
                  ) %>%
    dplyr::group_by(aoi) %>%
    dplyr::mutate(best = rmse == min(rmse)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = paste0(aoi
                                 , ": "
                                 , size
                                 , " ("
                                 , cv_method
                                 , ")"
                                 )
                  )

  sr_fits_short <- sr_fits %>%
    dplyr::select(size, aoi, cv_method, rmse, rsq, func, pack, best, label)

  library(ggplot2)
  library(patchwork)

  mod_plots <- sr_fits %>%
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


  resid_plots <- sr_fits %>%
    dplyr::mutate(final_fit = purrr::map(fit, "final_fit")) %>%
    dplyr::select(any_of(names(sr_fits_short)), final_fit) %>%
    tidyr::unnest(cols = c(final_fit)) %>%
    ggplot(aes(sr, resid, colour = best)) +
      geom_jitter(shape = ".") +
      geom_hline(aes(yintercept = 0), color = "gray50", lty = 2) +
      geom_smooth() +
      facet_grid(size ~ paste0(aoi, ": ", cv_method))

  imp_plots <- patchwork::wrap_plots(purrr::map2(sr_fits$imp_plot
                                                 , sr_fits$label
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
                  , cells = purrr::map_dbl(r, terra::ncell)
                  , sr_bii_extract = purrr::map(r
                                                , terra::extract
                                                , y = terra::vect(lsa)
                                                )
                  ) %>%
    dplyr::left_join(sr_fits_short) %>%
    tidyr::unnest(cols = sr_bii_extract) %>%
    dplyr::group_by(ID, across(any_of(names(sr_fits_short)))) %>%
    dplyr::summarise(cells = dplyr::n()
                    , sr_bii = mean(lyr1, na.rm = TRUE)
                    , sr_bii_se = sd(lyr1, na.rm = TRUE) / sqrt(cells)
                    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cells > 100) %>%
    dplyr::left_join(lsa %>%
                       sf::st_set_geometry(NULL) %>%
                       dplyr::mutate(ID = dplyr::row_number()) %>%
                       dplyr::select(ID, where(is.character))
                     )


  ggplot(sr_bii_tifs, aes(sr_bii
                          , as.factor(size)
                          , fill = best
                          )
         ) +
    geom_col() +
    geom_errorbarh(aes(xmin = sr_bii - sr_bii_se, xmax = sr_bii + sr_bii_se)) +
    geom_vline(aes(xintercept = 0)
               , linetype = 2
               , colour = "red"
               ) +
    facet_grid(LSA ~ paste0(aoi, ": ", cv_method))
