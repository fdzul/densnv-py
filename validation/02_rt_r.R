# Module 2 validation (Rt), R side. Run on the Mac (R console), from the root
# of densnv-py. Writes only Rt estimates (no case rows) to validation/out/.
path_txt <- "/Users/fdzul/Library/CloudStorage/Dropbox/dataset/dge/denv/2026/DENGUE2_.txt"
dir.create("validation/out", showWarnings = FALSE)
dataset <- data.table::fread(path_txt, quote = "", fill = TRUE)

areas <- list(
  mx      = list(level = "country"),
  yuc     = list(level = "state", cve_edo = "31"),
  merida  = list(level = "municipality", cve_edo = "31", cve_mpo = "050"),
  merida_loc = list(level = "locality", cve_edo = "31", cve_mpo = "050", cve_loc = "0001"),
  acapulco = list(level = "locality", cve_edo = "12", cve_mpo = "001", cve_loc = "0001"),
  yuc_res = list(level = "state", cve_edo = "31", place = "RES")
)
win <- list(); curve <- list()
for (k in names(areas)) {
  a <- areas[[k]]
  r <- do.call(densnv::rt_dengue,
               c(list(dataset = dataset, plot = FALSE, window_days = 28), a))
  win[[k]] <- cbind(area = k, r$rt_ventana[, c("casos", "rt", "rt_inf", "rt_sup", "fase")])
  if (!is.null(r$rt)) {
    R <- r$rt$R
    curve[[k]] <- data.frame(area = k, t_end = R$t_end, mean = R[["Mean(R)"]],
                             q025 = R[["Quantile.0.025(R)"]],
                             q975 = R[["Quantile.0.975(R)"]])
  }
}
write.csv(do.call(rbind, win), "validation/out/rt_window_r.csv", row.names = FALSE)
write.csv(do.call(rbind, curve), "validation/out/rt_curve_r.csv", row.names = FALSE)
