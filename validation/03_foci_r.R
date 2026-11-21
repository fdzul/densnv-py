# Module 3 validation (mp_active_foci), R side. Run on the Mac (R console),
# from the root of densnv-py. Writes to validation/out/ (local, not in git):
# the geocoded cases needed by Python and the results of R.
load("/Users/fdzul/Documents/geocoding_mx_2026_dengue/4.RData/dengue_mx_2026.RData")  # z
dir.create("validation/out", showWarnings = FALSE)

cols <- intersect(c("long", "lat", "ESTATUS_CASO", "FEC_INI_SIGNOS_SINT",
                    "CVE_EDO_REP", "CVE_MPO_REP", "CVE_LOC_REP", "DES_LOC_REP",
                    "CVE_EDO_RES", "CVE_MPO_RES", "CVE_LOC_RES", "DES_LOC_RES"), names(z))
zz <- as.data.frame(z)[, cols]
zz$FEC_INI_SIGNOS_SINT <- format(as.Date(zz$FEC_INI_SIGNOS_SINT), "%Y-%m-%d")
write.csv(zz, "validation/out/z_geocoded.csv", row.names = FALSE, fileEncoding = "UTF-8")

runs <- expand.grid(loc = c("31|Merida", "12|Acapulco de Juarez", "26|Hermosillo",
                            "27|Villahermosa", "03|La Paz"),
                    method = c("window", "st-dbscan"), stringsAsFactors = FALSE)
foci <- list(); diag <- list()
for (k in seq_len(nrow(runs))) {
  edo <- sub("\\|.*", "", runs$loc[k]); loc <- sub(".*\\|", "", runs$loc[k])
  r <- tryCatch(densnv::mp_active_foci(geocoded_dataset = z, cve_edo = edo, locality = loc,
                                       temporal_scale = 28, buffer_m = 0,
                                       method = runs$method[k], n_sim = 999, seed = 123,
                                       print_map = FALSE),
                error = function(e) { message(loc, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  id <- paste(edo, loc, runs$method[k], sep = "|")
  if (!is.null(r$foci)) {
    f <- sf::st_drop_geometry(r$foci)
    foci[[id]] <- data.frame(run = id, casos = f$casos, casos_recientes = f$casos_recientes,
                             confirmados = f$confirmados, knox_pares = f$knox_pares,
                             manzanas = f$manzanas, area_transmision_km2 = f$area_transmision_km2,
                             p_holm = f$p_holm)
  }
  t <- r$diagnostico$tabla
  diag[[id]] <- data.frame(run = id, knox = r$knox$knox, pares_espacio = r$knox$pares_espacio,
                           n_casos = r$knox$n_casos, knox_p = r$knox$p_value,
                           t[, c("casos_periodo", "casos_ventana", "rt", "fase", "incremento",
                                 "cobertura_casos", "pct_casos_en_focos", "pct_area_focos",
                                 "focos", "manzanas", "escenario")])
}
write.csv(do.call(rbind, foci), "validation/out/foci_r.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(do.call(rbind, diag), "validation/out/diag_r.csv", row.names = FALSE, fileEncoding = "UTF-8")
