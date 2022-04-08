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


  #------aoi agg-------

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
