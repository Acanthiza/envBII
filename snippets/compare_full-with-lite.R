
  recent <- 2019
  reference <- recent - 20

  spp <- c("Stipiturus malachurus"
           , "Melanodryas cucullata"
           , "Neophema chrysostoma"
           , "Stagonopleura guttata"
           )

  full_dir <- fs::path(".."
                       , "envTre"
                       , "out"
                       , "SA_0_class"
                       , "runs"
                       , "2022-02-07-0728"
                       )

  simple_dir <- fs::path(".."
                       , "envISI"
                       , "out"
                       , "mods"
                       )

  #------ plot function-------

  make_spp_plot <- function(mod_type
                            , taxa
                            , common
                            , mod_res
                            ) {

    plot_title <- bquote(~italic(.(taxa))*":" ~ .(common)* ". List-length corrected reporting rate")

    plot_subtitle <- paste0(stringr::str_to_title(mod_type)
                            , " model"
                            )

    plot_caption <- paste0("At median list length: "
                          , unique(mod_res$pred$list_length)
                          , " species"
                          )

    res <- list()

    temp <- mod_res$year_diff_df %>%
      dplyr::mutate(taxa = taxa
                    , common = common
                    )

    diff_l <- temp %>%
      dplyr::group_by(across(contains("IBRA"))) %>%
      dplyr::summarise(lower = sum(diff < 0) / n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(likelihood = map(lower
                                     , ~cut(.
                                            , breaks = c(0,lulikelihood$maxVal)
                                            , labels = lulikelihood$likelihood
                                            , include.lowest = TRUE
                                            )
                                     )
                    ) %>%
      tidyr::unnest(cols = c(likelihood)) %>%
      dplyr::mutate(likelihood = forcats::fct_expand(likelihood
                                                     , levels(lulikelihood$likelihood)
                                                     )
                    )

    diff <- temp %>%
      dplyr::left_join(diff_l)

    res$year_diff_plot <- ggplot(diff
                                 , aes(diff, fill = likelihood)
                                 ) +
      geom_density() +
      geom_vline(aes(xintercept = 0)
                 , linetype = 2
                 , colour = "red"
                 ) +
      scale_fill_viridis_d(drop = FALSE) +
      facet_wrap(~IBRA_SUB_N
                 , scales = "free"
                 ) +
      theme(axis.text.y = element_blank()
            , axis.ticks.y = element_blank()
            ) +
      labs(title = plot_title
          , fill = "Likelihood of decline"
          , subtitle = plot_subtitle
          , caption = plot_caption
          )

    res$plot_line <- mod_res$pred %>%
      ggplot(aes(x = year, y = pred)) +
      geom_line(aes(group = .draw)
                , alpha = 0.1
                ) +
      geom_vline(xintercept = c(reference, recent)
                 , linetype = 2
                 , colour = "red"
                 ) +
      facet_wrap(~IBRA_SUB_N
                 , scales = "free_y"
                 ) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title = plot_title
           , subtitle = plot_subtitle
           , caption = plot_caption
           )

    return(res)

  }


  # #-------Full--------
  # mod_type <- "Full"
  #
  # full_mod_path <- dir_ls(full_dir
  #                        , regexp = spp
  #                        ) %>%
  #   grep("list-length_mod", ., value = TRUE)
  #
  # full_mod <- rio::import(full_mod_path)
  #
  # full_mod_res <- make_mod_res(full_mod_path)
  #
  # full_mod_plots <- make_spp_plot(full_mod_res)
  #
  # #--------simple--------
  # mod_type <- "simple"
  #
  # simple_mod_path <- dir_ls(simple_dir
  #                        , regexp = spp
  #                        ) %>%
  #   grep("_summary", ., value = TRUE, invert = TRUE)V
  #
  # simple_mod <- rio::import(simple_mod_path)
  #
  # simple_mod_res <- make_mod_res(simple_mod_path)
  #
  # simple_mod_plots <- make_spp_plot(simple_mod_res)

  #----- comparisons-------

  comparisons <- taxa_taxonomy %>%
    dplyr::filter(taxa %in% spp) %>%
    dplyr::select(taxa, common) %>%
    dplyr::mutate(full_path = fs::path(full_dir
                                       , paste0("list-length_mod_", taxa, ".rds")
                                       )
                  , simple_path = fs::path(simple_dir
                                         , paste0(taxa, ".rds")
                                         )
                  ) %>%
    tidyr::pivot_longer(contains("_path")
                        , names_to = "mod_type"
                        , values_to = "mod_path"
                        ) %>%
    dplyr::mutate(mod_res = purrr::map(mod_path
                                       , make_mod_res
                                       , reference = reference
                                       , recent = recent
                                       , draws = 500
                                       )
                  , mod_type = gsub("_path", "", mod_type)
                  , mod_plot = pmap(list(mod_type
                                         , taxa
                                         , common
                                         , mod_res
                                         )
                                   , make_spp_plot
                                   )
                  )


  comp_plots <- comparisons %>%
    dplyr::select(-matches("_path|_res")) %>%
    tidyr::pivot_wider(names_from = "mod_type"
                       , values_from = "mod_plot"
                       ) %>%
    dplyr::mutate(patch_plot = purrr::map2(full
                                           , simple
                                           , ~patchwork::wrap_plots(.x[[1]]
                                                                    , .y[[1]]
                                                                    , .x[[2]]
                                                                    , .y[[2]]
                                                                    , guides = "collect"
                                                                    )
                                           )
                  )


  comp_plots$patch_plot[[4]]




