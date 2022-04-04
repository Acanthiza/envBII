
  #-------env data------

  rect_around_point <- function(x,xsize,ysize){
    bbox <- sf::st_bbox(x)
    bbox <- bbox +c(xsize/2,ysize/2,-xsize/2,-ysize/2)
    return(sf::st_as_sfc(bbox))
  }


  sr_rectangles <- sr_data %>%
    sf::st_as_sf(coords = c("long", "lat")
                 , crs = 4283
                 , remove = FALSE
                 ) %>%
    sf::st_transform(crs = 3577) %>%
    dplyr::mutate(ID = row_number()
                  , sf = purrr::map(geometry
                           , rect_around_point
                           , xsize = agg_cells*30/2
                           , ysize = agg_cells*30/2
                           )
                  ) %>%
    tidyr::unnest(cols = sf) %>%
    sf::st_set_geometry(NULL) %>%
    sf::st_as_sf(crs = 3577) %>%
    sf::st_transform(crs = 4283) %>%
    terra::vect()

  flor_env_dynamic <- purrr::map(epochs$r
                                , terra::extract
                                , y = sr_rectangles
                                , na.rm = TRUE
                                ) %>%
    stats::setNames(epochs$name) %>%
    dplyr::bind_rows(.id = "name") %>%
    tibble::as_tibble() %>%
    dplyr::filter(!is.na(LC_NAME)) %>%
    dplyr::left_join(sr_rectangles %>%
                       sf::st_as_sf() %>%
                       sf::st_set_geometry(NULL)
                     ) %>%
    dplyr::add_count(across(any_of(names(sr_data))), name, name = "cells") %>%
    dplyr::count(across(any_of(names(sr_data))), name, LC_NAME, cells, name = "count") %>%
    dplyr::mutate(p = count / cells) %>%
    tidyr::pivot_wider(id_cols = any_of(c("name", names(sr_data)))
                       , names_from = "LC_NAME"
                       , values_from = p
                       , values_fill = 0
                       )

  flor_env <- sr_data %>%
    dplyr::mutate(ID = row_number()) %>%
    dplyr::left_join(luyear) %>%
    dplyr::select(ID, lat, long, sr, name) %>%
    dplyr::inner_join(flor_env_dynamic) %>%
    dplyr::select(-ID, -name, -any_of(names(sr_data)), sr) %>%
    na.omit()

  out_file <- fs::path(out_dir, "fit.rds")

  if(!file.exists(out_file)) {

    make_and_fit_sr_model <- function(env_df
                                      , folds = 2
                                      , reps = 2
                                      , tune_size = 5
                                      , context = names(sr_data)
                                      ) {

      res <- list()

      res$original_data <- env_df

      env_vars <- names(env_df)[!names(env_df) %in% context]

      # Splits
      res$data_split <- initial_split(env_df
                                  , strata = sr
                                  )

      res$data_train <- training(res$data_split)
      res$data_test  <- testing(res$data_split)

      # Recipe
      res$rec <- recipe(env_df) %>%
        update_role(any_of(context)
                    , new_role = "id"
                    ) %>%
        update_role(sr
                    , new_role = "outcome"
                    ) %>%
        update_role(any_of(env_vars)
                 , new_role = "predictor"
                 ) %>%
        step_corr(all_predictors()
                  , threshold = 0.95
                  ) %>%
        step_sqrt(all_outcomes()
                 , skip = TRUE
                 )

      # Resampling
      res$folds <- vfold_cv(res$data_train, v = folds, repeats = reps)

      # Models
      res$rf_mod <- rand_forest(mode = "regression") %>%
        parsnip::set_args(mtry = tune()
                           , trees = tune()
                          ) %>%
        set_engine("randomForest")

      res$glm_mod <- linear_reg(mode = "regression") %>%
        parsnip::set_args(penalty = tune()
                          , mixture = tune()
                          ) %>%
        parsnip::set_engine("glmnet")

      res$xgb_mod <- boost_tree(mode = "regression") %>%
        parsnip::set_args(mtry = tune()
                          , trees = tune()
                          , min_n = tune()
                          , tree_depth = tune()
                          , learn_rate = tune()
                          , loss_reduction = tune()
                          , sample_size = tune()
                          , stop_iter = tune()
                          ) %>%
        parsnip::set_engine("xgboost")

      # Workflow
      res$mod_wf <- workflow_set(preproc = list(rec = res$rec)
                                 , models = list(glm = res$glm_mod
                                                 , rf = res$rf_mod
                                                 , xgb = res$xgb_mod
                                                 )
                                 )


      res$race_ctrl <- control_race(save_pred = TRUE
                                    , parallel_over = "everything"
                                    , save_workflow = TRUE
                                    )



      res$race_results <- res$mod_wf %>%
        workflow_map("tune_race_anova"
                     , seed = 8591
                     , resamples = res$folds
                     , grid = tune_size
                     , control = res$race_ctrl
                     )

      res$wf_best <- res$race_results %>%
        rank_results() %>%
        filter(.metric == "rmse") %>%
        dplyr::slice(1) %>%
        dplyr::pull(wflow_id)

      res$mod_best <- res$race_results %>%
        extract_workflow_set_result(res$wf_best) %>%
        select_best(metric = "rmse")

      # Final model
      res$final_wf <- res$mod_wf %>%
        extract_workflow(res$wf_best) %>%
        finalize_workflow(res$mod_best)

      # Fit final model to whole data set
      res$fit <- res$final_wf %>%
        fit(data = res$original_data)

      # Final results
      res$final_fit <- res$fit %>%
        augment(res$original_data)


      # Model plot
      res$mod_plot <- res$final_fit %>%
        ggplot(aes(sr, .pred^2)) +
        geom_point() +
        geom_abline(color = "gray50", lty = 2)

      # Residual plot
      res$resid_plot <- res$final_fit %>%
        dplyr::mutate(resid = sr - .pred^2) %>%
        ggplot(aes(sr, resid)) +
        geom_point() +
        geom_hline(aes(yintercept = 0)
                   , colour = "gray50"
                   , lty = 2
                   )

      # Importance plot
      res$imp_plot <- res$fit %>%
        extract_fit_parsnip() %>%
        vip::vip()

      return(res)

    }

    # fit <- make_and_fit_sr_model(flor_env, context = names(sr_data))

    fit <- make_and_fit_sr_model(flor_env
                                 , folds = 5
                                 , reps = 10
                                 , tune_size = 30
                                 , context = names(sr_data)
                                 )

    rio::export(fit, out_file)

  } else fit <- rio::import(out_file)


  #------predict--------

  purrr::walk2(env_stack
               , names(env_stack)
              , function(x, y) {

                 out_file <- fs::path(out_dir, paste0(y, ".tif"))

                 if(!file.exists(out_file)) {

                   terra::predict(x
                                  , fit$fit
                                  , na.rm = TRUE
                                  , type = "raw"
                                  , filename = out_file
                                  , overwrite = TRUE
                                  , cores = 1
                                  )

                 }

               }

               )


  sr_tifs <- purrr::map(names(env_stack)
                        , ~terra::rast(fs::path(out_dir, paste0(., ".tif")))
                        ) %>%
    stats::setNames(paste0(names(env_stack)))


  sr_tifs <- purrr::map(sr_tifs
                        , function(x) x^2
                        )

  #-------sr bii---------

  out_file <- fs::path(out_dir, "sr_bii.tif")

  if(!file.exists(out_file)) {

    fvi <- function(ref, current){

      v <- current / ref

      # z <- ref / current
      #
      # v <- ifelse(v <= 1 | is.na(v)
      #             , -1 * v
      #             , z
      #             )

      return(v)

      }

    sr_bii <- terra::lapp(c(sr_tifs[[1]], sr_tifs[[length(sr_tifs)]])
                          , fun = fvi
                          , filename = out_file
                          , overwrite = TRUE
                          )

  }

  sr_bii <- stars::read_stars(out_file)


  #-----visualise-------

  sr_ras_current_star <- stars::read_stars(out_file)

  plot_ras <- function(star_ras) {

    n_cells <- stars::st_dimensions(star_ras)$x$to * stars::st_dimensions(star_ras)$y$to

    sr_cuts <- c(hist(star_ras, 5, plot = FALSE)$breaks, Inf)

    tmap::tm_shape(star_ras
                   , raster.downsample = n_cells > 42468400/10
                   ) +
      tmap::tm_raster(breaks = sr_cuts
                      , palette = "viridis"
                      , midpoint = 1
                      )

  }

  plot_ras(sr_ras_current_star)

  hist(sr_ras_current_star)

