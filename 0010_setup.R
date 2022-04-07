
  #------Project-------

  if(!exists("current")) current <- "../envEco/out/KI_50_current"

  if(!exists("ls_size")) ls_size <- 1600

  agg_cells <- floor(sqrt((ls_size * ls_size) / (30 * 30)))

  out_dir <- fs::path("out"
                      , basename(current)
                      , ls_size
                      )

  data_dir <- fs::path("data"
                       , basename(current)
                       )

  if(!file.exists(out_dir)) fs::dir_create(out_dir)

  #-------packages-------

  library(magrittr)
  library(terra)
  library(envRaster)
  library(caret)
  library(dplyr)
  library(purrr)
  library(tmap)

  library(tidymodels)
  library(finetune)



  #-------overall settings-------

  visit_cols <- c("lat", "long", "year")

  tmap::tmap_mode("view")

  tmap::tmap_options(basemaps = c("OpenStreetMap.Mapnik"
                                  , "Esri.WorldImagery"
                                  )
                     , limits = c(facets.plot = 100)
                     )


  options(scipen = 999)


  #-------maps--------

  # IBRA Sub
  ibra_sub <- sf::st_read(fs::path("..","envEco", "out","shp","ibra_sub.shp")) %>%
    sf::st_transform(crs = 4283) %>%
    sf::st_make_valid()


  #-------lookups-----------
  # epochs
  luep <- rio::import("luEp.csv") %>%
    dplyr::select(!contains("short")) %>%
    dplyr::rename(name = landcover) %>%
    dplyr::mutate(name = gsub("-", "_", name))

  # epoch : year
  luyear <- envFunc::make_epochs(epoch_overlap = FALSE) %>%
    tidyr::unnest(cols = c(year)) %>%
    dplyr::left_join(luep %>%
                       dplyr::rename(end = epochEnd) %>%
                       dplyr::select(-epoch)
                     )


  # landcover
  lulandcover <- rio::import("luLandcover.csv") %>%
    tibble::as_tibble()

  #------import-------

  stuff <- fs::dir_info(current
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

  taxa_taxonomy <- stuff %>%
    dplyr::filter(grepl("taxa_taxonomy\\.rds", path)) %>%
    dplyr::select(cur_pot, new_obj, path) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    dplyr::select(obj) %>%
    tidyr::unnest(cols = c(obj))

  flor_tidy <- stuff %>%
    dplyr::filter(grepl("filt_summary\\.rds", path)) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    tidyr::unnest(cols = c(obj)) %>%
    dplyr::filter(name == "flor_tidy") %>%
    dplyr::select(obj) %>%
    tidyr::unnest(cols = c(obj)) %>%
    dplyr::inner_join(taxa_taxonomy %>%
                        dplyr::select(taxa, ind) %>%
                        dplyr::filter(ind == "Y")
                      )

  sr_data <- flor_tidy %>%
    dplyr::count(across(any_of(visit_cols)), name = "sr") %>%
    dplyr::left_join(luyear %>%
                       dplyr::select(year, name)
                     )

  aoi <- stuff %>%
    dplyr::filter(grepl("aoi\\.rds", path)) %>%
    dplyr::select(cur_pot, new_obj, path) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    dplyr::select(obj) %>%
    tidyr::unnest(cols = c(obj)) %>%
    sf::st_as_sf(crs = 3577) %>%
    sf::st_transform(4283)


  #-------Prepare landcover--------

  # based on code here: https://stackoverflow.com/questions/69546515/how-to-aggregate-categorical-spatraster

  landcover_folder <- "../../data/raster/landcover"

  epochs <- fs::dir_info(landcover_folder
                          , recurse = TRUE
                          ) %>%
    dplyr::select(path) %>%
    dplyr::filter(grepl("tif$", path)) %>%
    envFunc::filter_test_func() %>%
    dplyr::mutate(name = gsub("\\.tif","",basename(path))
                  , name = gsub("LANDCOVER", "Landcover", name)
                  ) %>%
    dplyr::left_join(luep) %>%
    dplyr::mutate(r = purrr::map(path, terra::rast)
                  , raw_path = fs::path(data_dir
                                      , "landcover"
                                      , "raw"
                                      , paste0(name, ".tif")
                                      )
                  , raw_exists = file.exists(raw_path)
                  , seg_path = fs::path(data_dir
                                        , "landcover"
                                        , "segregated"
                                        , paste0(name, ".tif")
                                        )
                  , seg_exists = file.exists(seg_path)
                  , agg_path = fs::path(data_dir
                                        , "landcover"
                                        , "aggregated"
                                        , ls_size
                                        , paste0(name, ".tif")
                                        )
                  , agg_exists = file.exists(agg_path)
                  )

  #------aoi raw--------

  raw_dir <- dirname(epochs$raw_path[[1]])

  if(!file.exists(raw_dir)) fs::dir_create(raw_dir)

  purrr::walk2(epochs$r[!epochs$raw_exists]
               , epochs$raw_path[!epochs$raw_exists]
               , ~terra::crop(.x
                              , y = vect(aoi)
                              , filename = .y
                              , overwrite = TRUE
                              )
               )


  #------aoi seg------

  seg_dir <- dirname(epochs$seg_path[[1]])

  if(!file.exists(seg_dir)) fs::dir_create(seg_dir)

  epochs <- epochs %>%
    dplyr::mutate(r_raw = purrr::map(raw_path, terra::rast))

  purrr::walk2(epochs$r_raw[!epochs$seg_exists]
               , epochs$seg_path[!epochs$seg_exists]
               , ~terra::segregate(.x
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )


  #------aggregation-------

  agg_dir <- dirname(epochs$agg_path[[1]])

  if(!file.exists(agg_dir)) fs::dir_create(agg_dir)

  epochs <- epochs %>%
    dplyr::mutate(r_seg = purrr::map(seg_path, terra::rast))

  agg_func <- function(x) {

    sum_x <- sum(x, na.rm = TRUE)
    length_x <- sum(!is.na(x))

    p <- sum_x / length_x

    return(p)

  }

  purrr::walk2(epochs$r_seg[!epochs$agg_exists]
               , epochs$agg_path[!epochs$agg_exists]
               , ~terra::aggregate(.x
                                   , fact = agg_cells
                                   , fun = agg_func
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )


  #-------Prepare satellite---------

  rasters_aligned <- fs::path("../../data/raster/aligned"
                              , gsub("_current|_potential", "", basename(current))
                              )

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
    dplyr::mutate(r = map(path
                          , terra::rast
                          )
                  , name = gsub("_aligned|\\.tif", "", basename(path))
                  , rep_path = fs::path(data_dir
                                        , "statics"
                                        , "reproject"
                                        , gsub("_aligned", "", basename(path))
                                        )
                  , rep_exists = file.exists(rep_path)
                  , agg_path = fs::path(data_dir
                                      , "statics"
                                      , "aggregate"
                                      , ls_size
                                      , gsub("_aligned", "", basename(path))
                                      )
                  , agg_exists = purrr::map_lgl(agg_path, file.exists)
                  )


  # reproject

  rep_dir <- dirname(statics$rep_path[[1]])

  if(!file.exists(rep_dir)) fs::dir_create(rep_dir)

  target_ras <- terra::rast(epochs$raw_path[[1]])

  purrr::walk2(statics$r[!statics$rep_exists]
               , statics$rep_path[!statics$rep_exists]
               , ~terra::project(.x
                                 , y = target_ras
                                 , filename = .y
                                 , overwrite = TRUE
                                 )
               )


  # aggregate

  agg_dir <- dirname(statics$agg_path[[1]])

  if(!file.exists(agg_dir)) fs::dir_create(agg_dir)

  statics <- statics %>%
    #envFunc::filter_test_func(test_col = "rep_path") %>%
    dplyr::mutate(r_rep = map(rep_path, terra::rast))

  target_ras <- terra::rast(epochs$agg_path[[1]])

  purrr::walk2(statics$r_rep[!statics$agg_exists]
               , statics$agg_path[!statics$agg_exists]
               , ~terra::resample(.x
                                  , y = target_ras
                                   , method = "average"
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )


  #------env stack-------

  epoch_stack <- epochs %>%
    dplyr::select(name, agg_path) %>%
    dplyr::mutate(s = purrr::map(agg_path, terra::rast)
                  , s = purrr::map(s, function(x) x %>% stats::setNames(lulandcover$LC_NAME))
                  )

  static_stack <- terra::rast(statics$agg_path) %>%
    stats::setNames(statics$name)




