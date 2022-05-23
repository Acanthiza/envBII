

#-----visualise-------

  plot_ras <- function(star_ras, bks = 5, mp = 0) {

    # bks - classes in continuous legend = bks + 2
    # mp - midpoint for continuous legend.

    n_cells <- stars::st_dimensions(star_ras)$x$to * stars::st_dimensions(star_ras)$y$to

    sr_cuts <- unique(c(-1, hist(star_ras, bks, plot = FALSE)$breaks, 1))

    tmap::tm_shape(star_ras
                   , raster.downsample = n_cells > 42468400/4
                   ) +
      tmap::tm_raster(breaks = sr_cuts
                      , palette = "viridis"
                      , midpoint = mp
                      )

  }

  tif <- terra::sources(sr_bii) %>%
    stars::read_stars()

  map_tif <- plot_ras(tif, 8, 0)

  hist(tif)



