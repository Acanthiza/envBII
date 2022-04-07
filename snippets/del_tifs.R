

  del_tifs <- fs::dir_ls(fs::path("data", "KI_50_current", "landcover", "raw")) %>%
    tibble::enframe(name = NULL, value = "path") %>%
    dplyr::filter(grepl("\\d{4}1\\.", path))

  if(FALSE) {

    fs::file_delete(del_tifs$path)

  }
