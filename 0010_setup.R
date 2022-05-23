
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

  library(magrittr)
  library(terra)
  library(envRaster)
  library(caret)
  library(dplyr)
  library(purrr)
  library(tmap)
  library(future)

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



  #-------maps--------

  # IBRA Sub
  ibra_sub <- sf::st_read(fs::path("..","envEco", "out","shp","ibra_sub.shp")) %>%
    sf::st_transform(crs = 4283) %>%
    sf::st_make_valid()

  # LSA
  lsa <- sf::st_read(fs::path("..","envEco", "out","shp","lsa.shp")) %>%
    sf::st_transform(crs = 4283) %>%
    sf::st_make_valid()


  #------palettes------

  lulsa <- rio::import(fs::path("data", "luLSA.csv"))

  # Set colours for LSAs - LSARegion
  lsa_palette <- mapply(FUN = function(red,green,blue,alpha) rgb(red
                                                                 , green
                                                                 , blue
                                                                 , alpha
                                                                 , maxColorValue = 255
                                                                 )
                        , red = lulsa$Red
                        , green = lulsa$Green
                        , blue = lulsa$Blue
                        , alpha = lulsa$Alpha
                        )


  names(lsa_palette) <- lulsa$LSA
  lsa_pal_fill <- scale_fill_manual(name = "LSA", values = lsa_palette, drop = TRUE)
  lsa_pal_col <- scale_colour_manual(name = "LSA", values = lsa_palette, drop = TRUE)
