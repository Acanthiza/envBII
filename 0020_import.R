

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
