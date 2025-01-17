
  #------env stack-------

  epoch_stack <- epochs %>%
    dplyr::select(name, agg_path) %>%
    dplyr::mutate(s = purrr::map(agg_path, terra::rast)
                  , s = purrr::map(s, function(x) x %>% stats::setNames(lulandcover$LC_NAME))
                  )

  static_stack <- terra::rast(statics$agg_path) %>%
    stats::setNames(statics$name)


  #-------env data------

  sr_points <- sr_data %>%
    sf::st_as_sf(coords = c("long", "lat")
                 , crs = 4283
                 , remove = FALSE
                 ) %>%
    dplyr::mutate(ID = row_number()) %>%
    sf::st_transform(crs = 7844) %>%
    terra::vect()

  out_file <- fs::path(out_dir, "sr_data_epochs.rds")

  if(!file.exists(out_file)) {

    sr_data_epochs <- purrr::map(epoch_stack$s
                                  , terra::extract
                                  , y = sr_points
                                  ) %>%
      stats::setNames(epochs$name) %>%
      dplyr::bind_rows(.id = "name") %>%
      tibble::as_tibble() %>%
      dplyr::inner_join(sr_points %>%
                         sf::st_as_sf() %>%
                         sf::st_set_geometry(NULL)
                       ) %>%
      dplyr::select(any_of(names(sr_data)), any_of(names(epoch_stack$s[[1]])))

    rio::export(sr_data_epochs
                , out_file
                )


  } else sr_data_epochs <- rio::import(out_file)


  out_file <- fs::path(out_dir, "sr_data_statics.rds")

  if(!file.exists(out_file)) {

    sr_data_statics <- terra::extract(static_stack
                                      , y = sr_points
                                      ) %>%
      tibble::as_tibble() %>%
      dplyr::inner_join(sr_points %>%
                         sf::st_as_sf() %>%
                         sf::st_set_geometry(NULL)
                       ) %>%
      dplyr::select(any_of(names(sr_data)), any_of(names(static_stack)))

    rio::export(sr_data_statics
                , out_file
                )

  } else sr_data_statics <- rio::import(out_file)


  sr_env <- sr_data %>%
    dplyr::inner_join(sr_data_epochs) %>%
    dplyr::inner_join(sr_data_statics) %>%
    na.omit()


  #-------Model--------

  out_file <- fs::path(out_dir, "fit.rds")

  if(!file.exists(out_file)) {

    make_and_fit_sr_model <- function(env_df
                                      , folds = 2
                                      , reps = 2
                                      , tune_size = 3
                                      , cv_method = "cv_normal"
                                      , context = names(sr_data)
                                      ) {

      res <- list()

      res$original_data <- env_df

      env_vars <- names(env_df)[!names(env_df) %in% c(context)]

      # Splits
      res$data_split <- initial_split(env_df
                                  , strata = sr
                                  )

      res$data_train <- training(res$data_split)
      res$data_test  <- testing(res$data_split)

      # Recipe
      res$rec <- recipe(sr ~ .
                        , env_df
                        ) %>%
        update_role(any_of(context)
                    , new_role = "context"
                    ) %>%
        update_role(-any_of(context)
                    , new_role = "predictor"
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
      res$folds <- if(cv_method == "cv_sptaial") {

        spatial_clustering_cv(res$data_train
                              , coords = c(lat, long)
                              , v = folds
                              #, repeats = reps
                              )


      } else {

        vfold_cv(res$data_train, v = folds, repeats = reps)

      }


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

      # res$xgb_mod <- boost_tree(mode = "regression") %>%
      #   parsnip::set_args(mtry = tune()
      #                     , trees = tune()
      #                     , min_n = tune()
      #                     , tree_depth = tune()
      #                     , learn_rate = tune()
      #                     , loss_reduction = tune()
      #                     , sample_size = tune()
      #                     , stop_iter = tune()
      #                     ) %>%
      #   parsnip::set_engine("xgboost")

      # Workflow
      res$mod_wf <- workflow_set(models = list(glm = res$glm_mod
                                                 , rf = res$rf_mod
                                                 #, xgb = res$xgb_mod
                                                 )
                                 , preproc = list(rec = res$rec)
                                 )


      res$race_ctrl <- control_race(save_pred = TRUE
                                    , parallel_over = "everything"
                                    , save_workflow = TRUE
                                    )



      res$wf_results <- res$mod_wf %>%
        workflowsets::workflow_map(if(cv_method == "cv_sptaial") "tune_grid" else "tune_race_anova"
                                   , seed = 8591
                                   , resamples = res$folds
                                   , grid = tune_size
                                   , control = res$race_ctrl
                                   )

      res$wf_plot <- autoplot(res$wf_results)

      res$wf_best <- res$wf_results %>%
        workflowsets::rank_results(select_best = TRUE) %>%
        dplyr::filter(.metric == "rmse") %>%
        dplyr::top_n(-1, mean) %>%
        dplyr::pull(wflow_id)

      res$mod_best <- res$wf_results %>%
        workflowsets::extract_workflow_set_result(res$wf_best) %>%
        tune::select_best(metric = "rmse")

      # Final model
      res$final_wf <- res$mod_wf %>%
        workflowsets::extract_workflow(res$wf_best) %>%
        tune::finalize_workflow(res$mod_best)

      # Final fit
      res$fit <- res$final_wf %>%
        fit(data = res$original_data)

      # Fit final model to test data set to get final metrics
      res$test <- res$final_wf %>%
        tune::last_fit(split = res$data_split)

      res$test_res <- res$test$.predictions[[1]] %>%
        dplyr::mutate(pred = .pred ^2
                      , resid = sr - pred
                      )

      # Importance plot
      res$imp_plot <- res$test %>%
        tune::extract_fit_parsnip() %>%
        vip::vip()


      return(res)

    }

    # fit <- make_and_fit_sr_model(flor_env, context = names(sr_data))

    fit <- make_and_fit_sr_model(sr_env
                                 , folds = 5
                                 , reps = 5
                                 , tune_size = 5
                                 , context = names(sr_data)
                                 )

    rio::export(fit, out_file)

  } else fit <- rio::import(out_file)

  mod_plot <- fit$test_res %>%
    ggplot(aes(sr
               , pred
               )
           ) +
    geom_jitter(shape = ".") +
    geom_abline(color = "gray50", lty = 2) +
    geom_smooth()


  resid_plot <- fit$test_res %>%
    ggplot(aes(sr, resid)) +
      geom_jitter(shape = ".") +
      geom_hline(aes(yintercept = 0), color = "gray50", lty = 2) +
      geom_smooth()

  imp_plot <- fit$imp_plot


  #------predict--------

  purrr::walk2(epoch_stack$s
               , epoch_stack$name
                , function(x, y) {

                 out_file <- fs::path(out_dir, paste0(y, ".tif"))

                 if(!file.exists(out_file)) {

                   # These vars are in the recipe fed to the tidymodels fit but
                    # were not predictors. Need to provide these to get
                    # terra::predict to work as it expects all the same
                    # columns as were provided to the fit
                   context_vars <- fit$rec$var_info %>%
                     dplyr::filter(role == "context") %>%
                     dplyr::select(variable, type)

                   context_vars_names <- as.list(context_vars$variable)

                   names(context_vars_names) <- context_vars$variable

                   context_vars_names[context_vars$type == "numeric"] <- 1
                   context_vars_names[context_vars$type != "numeric"] <- "1"

                   context <- as.data.frame(context_vars_names)

                   this_stack <- c(x, static_stack)

                   # Predict fit across the raster(s) in env_stack
                   terra::predict(this_stack
                                  , fit$fit
                                  , na.rm = TRUE
                                  , type = "raw"
                                  , filename = out_file
                                  , const = context
                                  , overwrite = TRUE
                                  , cores = 1
                                  )

                 }

               }

               )


  sr_tifs <- purrr::map(epoch_stack$name
                        , ~terra::rast(fs::path(out_dir, paste0(., ".tif")))
                        ) %>%
    stats::setNames(paste0(epoch_stack$name))


  sr_tifs <- purrr::map(sr_tifs
                        , function(x) x^2
                        )


  #-------sr bii---------

  out_file <- fs::path(out_dir, "sr_bii.tif")

  if(!file.exists(out_file)) {

    fvi <- function(ref, current){

      v <- current / ref

      z <- ref / current

      v <- ifelse(v <= 1 | is.na(v)
                  , -1 * (1 - v)
                  , 1 - z
                  )

      return(v)

      }

    sr_bii <- terra::lapp(c(sr_tifs[[1]], sr_tifs[[length(sr_tifs)]])
                          , fun = fvi
                          , filename = out_file
                          , overwrite = TRUE
                          )

  } else sr_bii <- terra::rast(out_file)



