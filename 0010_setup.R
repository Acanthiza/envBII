
  #------Project-------

  if(!exists("current")) current <- "../envEco/out/KI_50_current"

  if(!exists("ls_size")) ls_size <- 1600

  if(!exists("cv_method")) cv_method <- "cv_normal"

  agg_cells <- floor(sqrt((ls_size * ls_size) / (30 * 30)))


  # dir for data inputs
  data_dir <- fs::path("data"
                       , basename(current)
                       )

  if(!file.exists(data_dir)) fs::dir_create(data_dir)


  # dir for analysis outputs
  out_dir <- fs::path("out"
                      , basename(current)
                      , cv_method
                      , ls_size
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
  library(spatialsample)
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
