fetch_vix <- function() {
  csv_path <- "data/vix.csv"

  if (file.exists(csv_path)) {
    vix_data <- readr::read_csv(
      csv_path,
      col_types = readr::cols(
        date = readr::col_date(),
        vix_close = readr::col_double()
      )
    )
  } else {
    vix_xts <- quantmod::getSymbols(
      "^VIX",
      src = "yahoo",
      from = "1990-01-02",
      auto.assign = FALSE
    )

    vix_data <- tibble::tibble(
      date = zoo::index(vix_xts),
      vix_close = as.numeric(quantmod::Cl(vix_xts))
    )

    dir.create("data", showWarnings = FALSE)
    readr::write_csv(vix_data, csv_path)
  }

  vix_data |>
    dplyr::mutate(vix_log_return = log(vix_close / dplyr::lag(vix_close))) |>
    dplyr::filter(!is.na(vix_log_return))
}
