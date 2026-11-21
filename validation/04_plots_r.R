# Module 4 validation (graphs and maps), R side. Run on the Mac (R console),
# from the root of densnv-py. Writes the tables behind the figures (aggregated
# counts) to validation/out/ (local, not in git).
path_txt <- "/Users/fdzul/Library/CloudStorage/Dropbox/dataset/dge/denv/2026/DENGUE2_.txt"
dir.create("validation/out", showWarnings = FALSE)
dataset <- data.table::fread(path_txt, quote = "", fill = TRUE, encoding = "Latin-1")
w <- function(x, f) write.csv(as.data.frame(x), file.path("validation/out", f),
                              row.names = FALSE, fileEncoding = "UTF-8")
year <- 2026

# plt_heatmap: data of the plotly trace
h1 <- plotly::plotly_data(densnv::plt_heatmap(dataset, state = "YUCATAN", year = year,
                                              status = c(1, 2), hospitalizados = FALSE))
w(h1[, c("DES_MPO_REP", "week", "value")], "heatmap_state_r.csv")
h2 <- plotly::plotly_data(densnv::plt_heatmap(dataset, state = NULL, year = year,
                                              status = 2, hospitalizados = TRUE))
w(h2[, c("DES_MPO_REP", "week", "value")], "heatmap_country_r.csv")

# plt_pyramid: data of the ggplot
w(densnv::plt_pyramid(dataset)$data[, c("age_class", "sex", "count")], "pyramid_r.csv")

# plt_treemap: data of the ggplot
t1 <- densnv::plt_treemap(dataset, country = TRUE, year = year)$data
w(t1[, c("DES_EDO_REP", "DES_DIAG_FINAL", "value")], "treemap_country_r.csv")
t2 <- densnv::plt_treemap(dataset, country = FALSE, cve_edo = "31", year = year)$data
w(t2[, c("DES_MPO_REP", "DES_DIAG_FINAL", "value")], "treemap_state_r.csv")

# mp_bumpmap: the polygons with their value and position in the ranking
b1 <- densnv::mp_bumpmap(dataset, year = year, state = TRUE)
w(sf::st_drop_geometry(b1$layers[[2]]$data)[, c("NOMGEO", "value", "fine_cap")], "bump_country_r.csv")
b2 <- densnv::mp_bumpmap(dataset, year = year, state = FALSE, cve_edo = "31")
w(sf::st_drop_geometry(b2$layers[[2]]$data)[, c("NOMGEO", "value", "fine_cap")], "bump_state_r.csv")

# plt_serotype and mp_serotype return plots/maps without their tables: the same
# aggregations as their code
conf <- dataset |> dplyr::filter(ANO == year, ESTATUS_CASO == 2)
w(conf |> dplyr::filter(DES_DIAG_FINAL %in% c("DENGUE CON SIGNOS DE ALARMA", "DENGUE NO GRAVE",
                                               "DENGUE GRAVE"), DES_EDO_REP != "") |>
    dplyr::count(DES_EDO_REP, DES_DIAG_FINAL), "serotype_diag_r.csv")
w(conf |> dplyr::filter(DENGUE_SER_TRIPLEX %in% 1:4, DES_EDO_REP != "") |>
    dplyr::count(DES_EDO_REP, DENGUE_SER_TRIPLEX), "serotype_ser_r.csv")
ms <- dataset |> dplyr::filter(ESTATUS_CASO == 2, DENGUE_SER_TRIPLEX %in% 1:4,
                               sprintf("%02d", as.integer(CVE_EDO_REP)) == "31") |>
  dplyr::count(CVE_MPO_REP = sprintf("%03d", as.integer(CVE_MPO_REP)), DENGUE_SER_TRIPLEX)
w(ms, "mp_serotype_r.csv")
