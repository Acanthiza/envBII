

#------import-------

  stuff <- fs::dir_info(current
                        , recurse = TRUE
                        ) %>%
    dplyr::filter(type == "file"
                  , !grepl("landcover|random", path)
                  ) %>%
    dplyr::select(path, size) %>%
    dplyr::mutate(out_dir = stringr::str_extract(path, "\\w{2,5}_\\d{1,3}_\\w{6,20}")
                  , toi = stringr::str_extract(path, "(?<=SA_50_current\\/)[[:alpha:]]*(?=\\/)")
                  , date = lubridate::ymd_hm(stringr::str_extract(path, "\\d{4}-\\d{2}-\\d{2}-\\d{4}"))
                  ) %>%
    tidyr::separate(out_dir
                    , into = c("aoi", "buffer", "cur_pot")
                    , sep = "_|\\/"
                    ) %>%
    dplyr::group_by(toi, aoi, buffer, cur_pot) %>%
    dplyr::filter(date == max(date, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-size, - date)

  data <- stuff %>%
    dplyr::filter(grepl("bio_tidy|\\/taxa\\.rds", path)) %>%
    dplyr::mutate(file = gsub("\\..*|bio_", "", basename(path))) %>%
    tidyr::pivot_wider(names_from = file, values_from = path) %>%
    dplyr::mutate(tidy = purrr::map(tidy, rio::import)
                  , taxa = purrr::map(taxa
                                      , rio::import
                                      )
                  , tidy_unthemed = purrr::map(tidy
                                               , . %>%
                                                 dplyr::left_join(theme_prep) %>%
                                                 dplyr::filter(is.na(theme)) %>%
                                                 dplyr::count(taxa)
                                               )
                  , sr = purrr::map(tidy
                                    , . %>%
                                      dplyr::left_join(theme_prep) %>%
                                      dplyr::group_by(across(any_of(c(visit_cols, "theme")))) %>%
                                      dplyr::summarise(sr = n()) %>%
                                      dplyr::ungroup()
                                    )
                  )


  aoi <- stuff %>%
    dplyr::filter(grepl("aoi\\.rds", path)) %>%
    dplyr::mutate(obj = map(path, rio::import)) %>%
    dplyr::select(obj) %>%
    dplyr::distinct() %>%
    tidyr::unnest(cols = c(obj)) %>%
    sf::st_as_sf(crs = 3577) %>%
    sf::st_transform(7844)
