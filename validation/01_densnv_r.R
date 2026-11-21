# Module 1 validation, R side. Run on the Mac (R console), from the root of
# densnv-py. Writes aggregated counts only (no case rows) to validation/out/.
path_txt <- "/Users/fdzul/Library/CloudStorage/Dropbox/dataset/dge/denv/2026/DENGUE2_.txt"
dir.create("validation/out", showWarnings = FALSE)

loc <- densnv::cases_by_loc(path = path_txt, snv_dataset = "DENV") |>
  sf::st_drop_geometry()
write.csv(loc[, c("cvegeo", "probable", "confirmado", "DNG", "DCSA", "DG",
                  "DENV1", "DENV2", "DENV3", "DENV4")],
          "validation/out/loc_r.csv", row.names = FALSE)

agem <- densnv::cases_by_agem(path_vbd = dirname(path_txt), vbd = "DENGUE",
                              status_case = c(1, 2), country = TRUE) |>
  sf::st_drop_geometry()
write.csv(agem[, c("CVE_ENT", "CVE_MUN", "n")], "validation/out/agem_r.csv",
          row.names = FALSE)

st <- densnv::read_dengue_dataset(path = path_txt, spatial_resolution = "country",
                                  status_caso = c(1, 2))
write.csv(as.data.frame(table(CVE_EDO_RES = sprintf("%02d", as.integer(st$CVE_EDO_RES)))),
          "validation/out/country_r.csv", row.names = FALSE)
