

  #-----themes------

  themes <- c("IW", "CM", "TE")

  theme <- read_csv("data/ThSpp_DataPrepared.csv") %>%
    dplyr::select(original_name = Spp
                  , kingdom
                  , !!ensym(toi)
                  #, tsn
                  , any_of(themes)
                  ) %>%
    dplyr::distinct() %>%
    tidyr::pivot_longer(any_of(themes)
                        , names_to = "theme"
                        , values_to = "value"
                        ) %>%
    dplyr::filter(value > 0.2) %>%
    dplyr::filter(!is.na(kingdom))

  themeGBIF <- theme %>%
    tidyr::nest(data = -kingdom) %>%
    dplyr::mutate(gbif = purrr::map2(data
                                     , kingdom
                                     , ~envClean::get_gbif_tax(.x
                                                               , taxa_col = "original_name"
                                                               , out_file = fs::path("out"
                                                                                     , paste0(.y
                                                                                              , "_themeGBIF.csv"
                                                                                              )
                                                                                     )
                                                               , king_type = .y
                                                               )
                                     )
                  ) %>%
    dplyr::select(gbif) %>%
    tidyr::unnest(cols = c(gbif))

  theme_prep <- theme %>%
    dplyr::left_join(themeGBIF %>%
                       dplyr::distinct(original_name, taxa)
                     ) %>%
    dplyr::select(taxa, everything(), -original_name) %>%
    dplyr::distinct()

  rio::export(theme_prep
              , fs::path("out"
                         , "theme_prep.rds"
                         )
              )

  theme_missing <- bio_tidy %>%
    dplyr::left_join(theme_prep) %>%
    dplyr::filter(is.na(value)) %>%
    dplyr::count(taxa)

