

library(tidymodels)

set.seed(4595)

splits <- initial_split(use_data, strata = sr)

data_other <- training(splits)
data_test  <- testing(splits)

val_set <- validation_split(data_other
                            , strata = sr
                            , prop = 0.80
                            )

folds <- vfold_cv(data_train, v = 10)

rf_mod <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
  set_engine("ranger", num.threads = 2) %>%
  set_mode("regression")

rf_wf <- workflow() %>%
  add_model(rf_mod) %>%
  add_formula(sr ~ .)

rf_res <- rf_wf %>%
  tune_grid(val_set
            , grid = 25
            , control = control_grid(save_pred = TRUE)
            , metrics = metric_set(rmse)
            )

autoplot(rf_res)

rf_best <- rf_res %>%
  select_best(metric = "rmse")

resid_plot <- rf_res %>%
  collect_predictions(parameters = rf_best) %>%
  dplyr::mutate(resid = .pred - sr) %>%
  ggplot(aes(sr, resid)) +
    geom_point() +
    labs(y = "Residual", x = "Species richness") +
    # Scale and size the x- and y-axis uniformly:
    coord_obs_pred()
