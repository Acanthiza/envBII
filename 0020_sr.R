
  #-------env pca------

  flor_env_static_long <- purrr::map(env_stack
                                     , envRaster::get_env_data
                                     , df = contexts
                                     )

  flor_env_static <- purrr::map(flor_env_static_long
                                , function(x) x %>%
                                  dplyr::mutate(name = gsub("_aggregated", "", name)) %>%
                                  dplyr::filter(!is.na(value)) %>%
                                  dplyr::select(any_of(visit_cols), name, value) %>%
                                  tidyr::pivot_wider(names_from = name
                                                     , values_from = value
                                                     )
                                )

  flor_env <- purrr::map(flor_env_static
                         , function(x) contexts %>%
                           dplyr::left_join(x) %>%
                           dplyr::select(-lat, -long) %>%
                           na.omit()
                         )


  #------sr data--------

  # sr_data <- purrr::map(env_pca
  #                       , function(x) contexts %>%
  #                         dplyr::left_join(x$pca_res_cell) %>%
  #                         na.omit()
  #                       )
  #
  # sr_data_sf <- purrr::map(sr_data
  #                          , function(x) x %>%
  #                            sf::st_as_sf(coords = c("long", "lat")
  #                                         , crs = 4283
  #                                         )
  #                          )


  # ------sr model-------

  # pcs <- grep("^pc", names(sr_data[[1]]), value = TRUE)


  out_file <- "fit.rds"

  if(!file.exists(out_file)) {

    fit <- purrr::map(flor_env
                      , function(x) {

                        res <- list()

                        res$original_data <- x

                        env_vars <- names(x)[!names(x) %in% names(contexts)]

                        # Splits
                        res$data_split <- initial_split(x
                                                    , strata = sr
                                                    )

                        res$data_train <- training(res$data_split)
                        res$data_test  <- testing(res$data_split)

                        # Recipe
                        res$rec <- recipe(x) %>%
                          update_role(any_of(names(contexts))
                                      , new_role = "id"
                                      ) %>%
                          update_role(year
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
                                    )

                        # Resampling
                        res$folds <- vfold_cv(res$data_train, v = 5, repeats = 20)

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
                                       , grid = 30
                                       , control = res$race_ctrl
                                       )

                        res$wf_best <- res$race_results %>%
                          rank_results() %>%
                          filter(.metric == "rsq") %>%
                          dplyr::slice(1) %>%
                          dplyr::pull(wflow_id)

                        res$mod_best <- res$race_results %>%
                          extract_workflow_set_result(res$wf_best) %>%
                          select_best(metric = "rmse")

                        # Final model
                        res$final_mod <- res$mod_wf %>%
                          extract_workflow(res$wf_best) %>%
                          finalize_workflow(res$mod_best)

                        # Fit final model to whole data set
                        res$fit <- res$final_mod %>%
                          fit(data = res$original_data)

                        # Final results
                        res$final_fit <- res$fit %>%
                          augment(res$original_data)


                        # Model plot
                        res$mod_plot <- res$final_fit %>%
                          ggplot(aes(sr, .pred)) +
                          geom_point() +
                          geom_abline(color = "gray50", lty = 2)

                        # Residual plot
                        res$resid_plot <- res$final_fit %>%
                          dplyr::mutate(resid = sr - .pred) %>%
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

                      )


    rio::export(fit, out_file)

  } else fit <- rio::import(out_file)

  mods <- map(fit, "fit")

  #------predict--------

  purrr::pwalk(list(env_stack
                    , mods
                    , names(mods)
                    )
               , function(x, y, z) {

                 out_file <- fs::path(paste0("sr_", z, ".tif"))

                 if(!file.exists(out_file)) {

                   terra::predict(x
                                  , y
                                  , na.rm = TRUE
                                  , const = data.frame(year = 2020L)
                                  , type = "raw"
                                  , filename = out_file
                                  , overwrite = TRUE
                                  , cores = 1
                                  )

                 }

               }

               )



  sr_potential <- terra::rast("sr_potential.tif")
  sr_current <- terra::rast("sr_current.tif")

  # Align rasters
  sr_current <- terra::project(sr_current
                             , sr_potential
                             , filename = "temp.tif"
                             )

  fs::file_delete("sr_current.tif")
  fs::file_copy("temp.tif", "sr_current.tif")

  sr_current <- terra::rast("sr_current.tif")





  #-------sr bii---------

  out_file <- fs::path("sr_bii.tif")

  if(!file.exists(out_file)) {

    fvi <- function(ref, current){

      v <- current / ref

      #v <- ifelse(v <= 1 | is.na(v), v, 1)

      return(v)

      }

    sr_bii <- terra::lapp(c(sr_potential, sr_current)
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
                      )

  }

  plot_ras(sr_ras_current_star)

  hist(sr_ras_current_star)

