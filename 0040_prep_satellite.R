
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


  #-------static reproject-------

  rep_dir <- dirname(statics$rep_path[[1]])

  if(!file.exists(rep_dir)) fs::dir_create(rep_dir)

  target_ras <- epochs$r[[1]]

  purrr::walk2(statics$r[!statics$rep_exists]
               , statics$rep_path[!statics$rep_exists]
               , ~terra::project(.x
                                 , y = target_ras
                                 , filename = .y
                                 , overwrite = TRUE
                                 )
               )


  #------static agg---------

  agg_dir <- dirname(statics$agg_path[[1]])

  if(!file.exists(agg_dir)) fs::dir_create(agg_dir)

  statics <- statics %>%
    #envFunc::filter_test_func(test_col = "rep_path") %>%
    dplyr::mutate(r_rep = map(rep_path, terra::rast))

  purrr::walk2(statics$r_rep[!statics$agg_exists]
               , statics$agg_path[!statics$agg_exists]
               , ~terra::aggregate(x = .x
                                   , fact = agg_cells
                                   , fun = "mean"
                                   , na.rm = TRUE
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )
