
  current <- "D:/env/projects/envEco/out/KI_50_current"

  aggregation_landscape <- 2000 # m

  agg_cells <- floor(sqrt((aggregation_landscape * aggregation_landscape) / (30 * 30)))

  out_dir <- fs::path("out"
                      , basename(current)
                      , aggregation_landscape
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
                        dplyr::filter(ind != "N")
                      )

  sr_data <- flor_tidy %>%
    dplyr::count(across(any_of(visit_cols)), name = "sr")


  #-------Prepare landcover--------
  # based on code here: https://stackoverflow.com/questions/69546515/how-to-aggregate-categorical-spatraster

  raster_folder <- "../../data/raster/landcover"

  epochs <- fs::dir_info(raster_folder
                          , recurse = TRUE
                          ) %>%
    dplyr::select(path) %>%
    dplyr::filter(grepl("tif$", path)
                  , grepl("raw", path)
                  ) %>%
    envFunc::filter_test_func() %>%
    dplyr::mutate(name = gsub("1\\.tif","",basename(path))) %>%
    dplyr::left_join(luep) %>%
    dplyr::mutate(seg_path = fs::path(gsub("raw", "segregated", dirname(path))
                                      , basename(path)
                                      )
                  , r = map(path, terra::rast)
                  , seg_exists = purrr::map_lgl(seg_path, file.exists)
                  , agg_path = fs::path(gsub("raw", "aggregated", dirname(path))
                                        , aggregation_landscape
                                        , basename(path)
                                        )
                  , r = map(path, terra::rast)
                  , agg_exists = purrr::map_lgl(agg_path, file.exists)
                  )

  #------segregation------

  seg_dir <- dirname(epochs$seg_path[[1]])

  if(!file.exists(seg_dir)) fs::dir_create(seg_dir)

  purrr::walk2(epochs$r[!epochs$seg_exists]
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
    dplyr::mutate(seg = purrr::map(seg_path, terra::rast))

  agg_func <- function(x) {

    sum_x <- sum(x, na.rm = TRUE)
    length_x <- sum(!is.na(x))

    p <- sum_x / length_x

    return(p)

  }

  purrr::walk2(epochs$seg[!epochs$agg_exists]
               , epochs$agg_path[!epochs$agg_exists]
               , ~terra::aggregate(.x
                                   , fact = agg_cells
                                   , fun = agg_func
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )

  #------env stack-------

  env_stack <- epochs %>%
    dplyr::mutate(s = purrr::map(agg_path, terra::rast)
                  , s = purrr::map(s, function(x) x %>% stats::setNames(lulandcover$LC_NAME))
                  ) %>%
    dplyr::pull(s) %>%
    stats::setNames(epochs$name)




