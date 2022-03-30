
  potential <- "D:/env/projects/envEco/out/KI_5_potential"
  current <- "D:/env/projects/envEco/out/KI_50_current"

  #-------overall settings-------

  visit_cols <- c("lat", "long", "year")

  # How many pca axes to use?
  pca_axes <- 4

  # What style to use to break up pca axes?
  find_interval_style <- "quantile"

  # How many cuts along pc 1?
  pca_cuts <- 20
  # cuts on pc 2 will be 20/2
  # cuts on pc3 will be 20/3
  # etc.


  tmap::tmap_mode("view")

  tmap::tmap_options(basemaps = c("OpenStreetMap.Mapnik"
                                  , "Esri.WorldImagery"
                                  )
                     , limits = c(facets.plot = 100)
                     )

  #-------packages-------

  library(magrittr)
  library(terra)
  library(envRaster)
  library(caret)
  library(dplyr)
  library(purrr)

  library(tidymodels)
  library(finetune)


  #------import-------

  stuff <- fs::dir_info(c(current, potential)
                     , recurse = TRUE
                     ) %>%
    dplyr::mutate(out_dir = stringr::str_extract(path, "\\w{2,5}_\\d{1,3}_\\w{6,20}")
                  , out_time = stringr::str_extract(path, "\\d{4}-\\d{2}-\\d{2}-\\d{4}")
                  , out_time = lubridate::ymd_hm(out_time)
                  , obj = gsub("\\..*", "", basename(path))
                  ) %>%
    tidyr::separate(out_dir
                    , into = c("aoi", "buffer", "cur_pot")
                    , sep = "_"
                    ) %>%
    dplyr::select(path, type, out_time, aoi, buffer, cur_pot, size, obj) %>%
    dplyr::group_by(aoi, buffer, cur_pot) %>%
    dplyr::mutate(tot_size = sum(size)) %>%
    dplyr::filter(type != "directory"
                  , tot_size > 1e+9
                  , out_time == max(out_time, na.rm = TRUE)
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(new_obj = paste0(obj, "_", cur_pot))

  # tifs <- stuff %>%
  #   dplyr::filter(grepl("\\.tif$", path)) %>%
  #   dplyr::mutate(ras = purrr::map(path, terra::rast))
  #
  # eco_ras <- tifs %>%
  #   dplyr::filter(grepl("eco_pred", path)) %>%
  #   dplyr::pull(ras) %>%
  #   `[[`(1)
  #
  # sr_ras <- tifs %>%
  #   dplyr::filter(grepl("sr_ras", path)) %>%
  #   dplyr::pull(ras) %>%
  #   `[[`(1)
  #
  # sr_ras_star <- tifs %>%
  #   dplyr::filter(grepl("sr_ras", path)) %>%
  #   dplyr::pull(path) %>%
  #   stars::read_stars()

  env_pca <- stuff %>%
    dplyr::filter(grepl("env_pca", obj)) %>%
    dplyr::select(cur_pot, new_obj, path) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    dplyr::pull(obj, new_obj)

  taxa_taxonomy <- stuff %>%
    dplyr::filter(grepl("taxa_taxonomy\\.rds", path)
                  , cur_pot == "potential"
                  ) %>%
    dplyr::select(cur_pot, new_obj, path) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    dplyr::pull(obj, new_obj)

  flor_tidy <- stuff %>%
    dplyr::filter(grepl("filt_summary\\.rds", path)) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    tidyr::unnest(cols = c(obj)) %>%
    dplyr::filter(name == "flor_tidy"
                  , cur_pot == "potential"
                  ) %>%
    dplyr::mutate(new_obj = paste0(name, "_", cur_pot)) %>%
    dplyr::pull(obj, new_obj) %>%
    purrr::map(function(x) x %>%
                 dplyr::inner_join(taxa_taxonomy$taxa_taxonomy_potential %>%
                                      dplyr::select(taxa, ind) %>%
                                      dplyr::filter(ind != "N")
                                    )
               )

  contexts <- flor_tidy$flor_tidy_potential %>%
    dplyr::count(across(any_of(visit_cols)), name = "sr")

  flor_env <- stuff %>%
    dplyr::filter(grepl("flor_env", obj)) %>%
    dplyr::select(cur_pot, new_obj, path) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    dplyr::pull(obj, new_obj) %>%
    purrr::map(function(x) x %>%
                 dplyr::inner_join(contexts %>%
                                     dplyr::select(-sr)
                                   )
               )

  lc_tidy <-  stuff %>%
    dplyr::filter(grepl("filt_summary_lc\\.rds", path)) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    tidyr::unnest(cols = c(obj)) %>%
    dplyr::filter(name == "lc_tidy"
                  , cur_pot == "current"
                  ) %>%
    dplyr::mutate(new_obj = paste0(name, "_", cur_pot)) %>%
    dplyr::pull(obj, new_obj)


  #--------current env-------

  # Static are used for predict (and train, if they are only static)
  rasters_aligned <- "../../data/raster/aligned/ki_50"

  use_ras_type <- unique(envEcosystems::env$process[envEcosystems::env$provider != "KIDTM1m"])

  no_ras_type <- c("^_"
                   , "isocluster"
                   , "twi"
                   , "test"
                   , "90-"
                   , "00-"
                   , "median"
                   , "winter"
                   , "spring"
                   #, "autumn"
                   , "sd"
                   #, "min"
                   #, "max"
                   #, "03\\d{4}05" # autumn
                   , "06\\d{4}08" # winter
                   , "09\\d{4}11" # spring
                   )

  statics <- envRaster::parse_env_tif(rasters_aligned) %>%
    dplyr::filter(process %in% use_ras_type) %>%
    dplyr::filter(!grepl(paste0(no_ras_type,collapse = "|")
                         , path
                         )
                  ) %>%
    envFunc::filter_test_func() %>%
    dplyr::mutate(name = gsub("_aligned\\.tif","",basename(path))
                  , agg = paste0(name, "_aggregated.tif")
                  , agg_path = fs::path(gsub("aligned", "aggregated", dirname(path))
                                        , agg
                                        )
                  , r = map(path, terra::rast)
                  , agg_exists = purrr::map_lgl(agg_path, file.exists)
                  )

  agg_dir <- dirname(statics$agg_path[[1]])

  if(!file.exists(agg_dir)) fs::dir_create(agg_dir)

  purrr::walk2(statics$r[!statics$agg_exists]
               , statics$agg_path[!statics$agg_exists]
               , ~terra::aggregate(.x
                                   , fact = 30
                                   , fun = "mean"
                                   , na.rm = TRUE
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )

  env_stack <- list()

  env_stack$current <- statics %>%
    dplyr::pull(agg_path) %>%
    terra::rast() %>%
    stats::setNames(statics$name)


  #--------potential env-------

  # Static are used for predict (and train, if they are only static)
  rasters_aligned <- "../../data/raster/aligned/ki_5"

  # Potential
  use_ras_type <- unique(envEcosystems::env$process[envEcosystems::env$group != "satellite"])

  no_ras_type <- c("^_"
                   , "isocluster"
                   , "twi"
                   , "test"
                   , "90-"
                   , "00-"
                   , "median"
                   , "winter"
                   , "spring"
                   #, "autumn"
                   , "sd"
                   #, "min"
                   #, "max"
                   #, "03\\d{4}05" # autumn
                   , "06\\d{4}08" # winter
                   , "09\\d{4}11" # spring
                   )

  statics <- envRaster::parse_env_tif(rasters_aligned) %>%
    dplyr::filter(process %in% use_ras_type) %>%
    dplyr::filter(!grepl(paste0(no_ras_type,collapse = "|")
                         , path
                         )
                  ) %>%
    envFunc::filter_test_func() %>%
    dplyr::mutate(name = gsub("_aligned\\.tif","",basename(path))
                  , agg = paste0(name, "_aggregated.tif")
                  , agg_path = fs::path(gsub("aligned", "aggregated", dirname(path))
                                        , agg
                                        )
                  , r = map(path, terra::rast)
                  , agg_exists = purrr::map_lgl(agg_path, file.exists)
                  )

  agg_dir <- dirname(statics$agg_path[[1]])

  if(!file.exists(agg_dir)) fs::dir_create(agg_dir)

  purrr::walk2(statics$r[!statics$agg_exists]
               , statics$agg_path[!statics$agg_exists]
               , ~terra::aggregate(.x
                                   , fact = 30
                                   , fun = "mean"
                                   , na.rm = TRUE
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )

  env_stack$potential <- statics %>%
    dplyr::pull(agg_path) %>%
    terra::rast() %>%
    stats::setNames(statics$name)


