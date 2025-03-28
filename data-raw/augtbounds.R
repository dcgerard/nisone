## code to prepare `augtbounds` dataset goes here
df <- expand.grid(n = seq(2, 100), alpha = c(0.01, 0.05, 0.1, 0.2))
df$eta <- NA_real_
df$nu <- NA_real_
for (i in seq_len(nrow(df))) {
  eout <- eta_alpha(alpha = df$alpha[[i]], n = df$n[[i]])
  df$eta[[i]] <- eout[["eta"]]
  df$nu[[i]] <- eout[["nu"]] ## nu is (A - mu) / (sigma / sqrt(n)), which is different than in paper
}
augtbounds <- df
usethis::use_data(augtbounds, overwrite = TRUE, internal = TRUE)
