
  #------Project-------

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

  library(tidyverse)

  library(terra)
  library(envRaster)
  library(caret)
  library(tmap)
  library(future)

  library(tidymodels)
  library(spatialsample)
  library(finetune)



  #-------overall settings-------

  toi <- "class"

  visit_cols <- c("lat", "long", "year", "toi")

  tmap::tmap_mode("view")

  tmap::tmap_options(basemaps = c("OpenStreetMap.Mapnik"
                                  , "Esri.WorldImagery"
                                  )
                     , limits = c(facets.plot = 100)
                     )


  options(scipen = 999)


  future::plan(sequential)

  future::plan(multisession
               , workers = if(parallel::detectCores() > max_cores) max_cores else parallel::detectCores() - 1
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


  # generic data
  xfun::in_dir(fs::path("..", "RC")
               , source(fs::path("common", "genericData.R"))
               )
