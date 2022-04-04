
 #---------To github--------

 now <- format(Sys.time(), "%Y-%m-%d %H:%M")

  if(nrow(gert::git_status() > 1)) {

    envFunc::git_commit_env(paste0(commit_notes
                                   , ". Successful run: "
                                   , now
                                   )
                            )

  }


  #------To network-------

  base_network_path <- fs::path("//env.sa.gov.au/dfsroot/IST/DEHProjects/Landscapes/envBII")

  filesRegex <- c("gpkg$" # geopackage
                  , "tif$" # raster
                  , "xml$"
                  , "report.docx"
                  , "ecosystemsDesc.csv"
                  , "slides.html"
                  , "short.docx"
                  )

  files_to_copy <- tibble::tibble(path = fs::dir_ls(regexp = paste0(filesRegex,collapse = "|"))) %>%
    dplyr::filter(!grepl("~",path))


  fs::dir_create(base_network_path)

  purrr::walk(files_to_copy$path
              , ~file.copy(from = .
                           , to = fs::path(base_network_path, .)
                           , overwrite = TRUE
                           )
              )




