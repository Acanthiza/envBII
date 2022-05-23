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


  #------aoi seg------

  seg_dir <- dirname(epochs$seg_path[[1]])

  if(!file.exists(seg_dir)) fs::dir_create(seg_dir)

  purrr::walk2(epochs$r[!epochs$seg_exists]
               , epochs$seg_path[!epochs$seg_exists]
               , ~terra::segregate(.x
                                   , filename = .y
                                   , overwrite = FALSE
                                   )
               )


  #------aoi agg-------

  agg_dir <- dirname(epochs$agg_path[[1]])

  if(!file.exists(agg_dir)) fs::dir_create(agg_dir)

  epochs <- epochs %>%
    dplyr::mutate(r_seg = purrr::map(seg_path, terra::rast))

  purrr::walk2(epochs$r_seg[!epochs$agg_exists]
               , epochs$agg_path[!epochs$agg_exists]
               , ~terra::aggregate(.x
                                   , fact = agg_cells
                                   , fun = "mean"
                                   , na.rm = TRUE
                                   , filename = .y
                                   , overwrite = TRUE
                                   )
               )
