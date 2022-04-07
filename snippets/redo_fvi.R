
  # make sure fvi from 0020_sr.R is available



  # Assemble tifs and remake sr_bii.tif using fvi

  tifs <- fs::dir_ls(recurse = TRUE
                     , regexp = "\\d{4}\\.tif$"
                     ) %>%
    tibble::enframe(name = NULL
                    , value = "path"
                    ) %>%
    dplyr::mutate(agg_size = basename(dirname(path))
                  , name = gsub("\\.tif", "", basename(path))
                  ) %>%
    dplyr::filter(agg_size %in% sizes) %>%
    dplyr::arrange(agg_size) %>%
    tidyr::nest(data = -agg_size) %>%
    dplyr::mutate(r_list = purrr::map(data
                               , function(x) purrr::map2(x$path
                                                         , x$name
                                                         , ~terra::rast(.x) %>%
                                                           stats::setNames(.y)
                                                         )
                               )
                  , r_list = purrr::map(r_list
                                        , stats::setNames
                                        , nm = names(env_stack)
                                        )
                  , out_path = purrr::map_chr(agg_size
                                              , ~fs::path(dirname(out_dir)
                                                          , .
                                                          , "sr_bii.tif"
                                                          )
                                              )
                  , sr_bii = purrr::map2(r_list
                                         , out_path
                                        , ~terra::lapp(c(.x[[1]], .x[[length(.x)]])
                                                       , fun = fvi
                                                       , filename = .y
                                                       , overwrite = TRUE
                                                       )
                                        )
                  )

  # Copy tifs to network

  base_network_path <- fs::path("//env.sa.gov.au/dfsroot/IST/DEHProjects/Landscapes/envBII")

  tifs_to_copy <- fs::dir_ls(regexp = "sr_bii"
                            , recurse = TRUE
                            ) %>%
    tibble::enframe(name = NULL, value = "from_path") %>%
    dplyr::mutate(to_path = purrr::map_chr(from_path
                                       , ~fs::path(base_network_path
                                                  , .
                                                  )
                                       )
                  )


  purrr::walk2(tifs_to_copy$from_path
               , tifs_to_copy$to_path
               , file.copy
               , overwrite = TRUE
               )

