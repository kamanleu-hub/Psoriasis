
#Arthor: Chi-Erh Chiu
# Study-level table
studies <- data_csv %>%
  mutate(
    drug = if_else(!is.na(Roflumilast), "Roflumilast", "Tapinarof"),
    study_id = trimws(sub("\\{.*", "", ID))
  ) %>%
  select(
    refid = Refid,
    study_id,
    trial_id = trial_ID,
    drug,
    n_total = size,
    followup_months = followup,
    pct_female = percent_female,
    rob_overall = RoB
  ) %>%
  arrange(drug, study_id)
print(studies, n = Inf)


# Make study IDs unique
make_unique <- function(d) {
  d %>% mutate(study_id = paste0(study_id, " [", refid, "]"))
}

# Dose labels 
dose_map <- tibble(
  refid = c(1, 5, 59, 948, 3174, 3434, 3491, 3583, 3641, 3642, 3647, 3648),
  dose  = c("0.3% foam", "0.3% cream", "1% cream", "0.3% foam",
            "1% cream",  "1% cream",  "1% cream", "1% cream",
            "1% cream",  "1% cream",  "0.3% cream", "0.3% cream")
)

add_dose <- function(d) {
  d %>%
    left_join(dose_map, by = "refid") %>%
    mutate(study_id = paste0(study_id, " (", dose, ")")) %>%
    select(-dose)
}

#RoB
table(studies$rob_overall)
studies %>% 
  select(study_id, drug, rob_overall) %>% 
  arrange(rob_overall, drug) %>%
  print(n = Inf)

# Build outcome datasets
# IGA/PGA success 
iga <- data_csv %>%
  mutate(
    drug = if_else(!is.na(Roflumilast), "Roflumilast", "Tapinarof"),
    study_id = trimws(sub("\\{.*", "", ID)),
    event_int = if_else(Refid == 3434,
                        clinical_severity_counts_comp_categ,
                        clinical_severity_counts_int_categ),
    n_int     = if_else(Refid == 3434,
                        clinical_severity_n_comp_categ,
                        clinical_severity_n_int_categ)
  ) %>%
  select(
    refid = Refid, study_id, drug,
    event.e = event_int, n.e = n_int,
    event.c = clinical_severity_counts_ctrl_categ,
    n.c     = clinical_severity_n_ctrl_categ
  ) %>%
  arrange(drug, study_id) %>%
  make_unique() %>%
  add_dose()

print(iga, n = Inf)


# Any AE 
ae <- data_csv %>%
  mutate(
    drug = if_else(!is.na(Roflumilast), "Roflumilast", "Tapinarof"),
    study_id = trimws(sub("\\{.*", "", ID)),
    event_int = if_else(Refid == 3434,
                        nAEparticipants_counts_comp_categ,
                        nAEparticipants_counts_int_categ),
    n_int     = if_else(Refid == 3434,
                        nAEparticipants_n_comp_categ,
                        nAEparticipants_n_int_categ)
  ) %>%
  select(
    refid = Refid, study_id, drug,
    event.e = event_int, n.e = n_int,
    event.c = nAEparticipants_counts_ctrl_categ,
    n.c     = nAEparticipants_n_ctrl_categ
  ) %>%
  arrange(drug, study_id) %>%
  make_unique() %>%
  add_dose()


# GI AEs
gi <- data_csv %>%
  mutate(
    drug = if_else(!is.na(Roflumilast), "Roflumilast", "Tapinarof"),
    study_id = trimws(sub("\\{.*", "", ID))
  ) %>%
  select(
    refid = Refid, study_id, drug,
    event.e = GI_counts_int_categ, n.e = GI_n_int_categ,
    event.c = GI_counts_ctrl_categ, n.c = GI_n_ctrl_categ
  ) %>%
  arrange(drug, study_id) %>%
  make_unique() %>%
  add_dose()


# Irritation
irritation <- data_csv %>%
  mutate(
    drug = if_else(!is.na(Roflumilast), "Roflumilast", "Tapinarof"),
    study_id = trimws(sub("\\{.*", "", ID))
  ) %>%
  select(
    refid = Refid, study_id, drug,
    event.e = irritation_counts_int_categ, n.e = irritation_n_int_categ,
    event.c = irritation_counts_ctrl_categ, n.c = irritation_n_ctrl_categ
  ) %>%
  arrange(drug, study_id) %>%
  make_unique() %>%
  add_dose()


# PSD change from baseline (continuous)
psd <- data_csv %>%
  mutate(
    drug = if_else(!is.na(Roflumilast), "Roflumilast", "Tapinarof"),
    study_id = trimws(sub("\\{.*", "", ID))
  ) %>%
  select(
    refid = Refid, study_id, drug,
    n.e    = patient_symptoms_n_int_cont,
    mean.e = patient_symptoms_mean_int_cont,
    sd.e   = patient_symptoms_SD_int_cont,
    n.c    = patient_symptoms_n_ctrl_cont,
    mean.c = patient_symptoms_mean_ctrl_cont,
    sd.c   = patient_symptoms_SD_ctrl_cont
  ) %>%
  arrange(drug, study_id)

# Organon SE to SD correction 
psd_fixed <- psd %>%
  mutate(
    se_to_sd_corrected = refid %in% c(3641, 3642),
    sd.e = if_else(se_to_sd_corrected, sd.e * sqrt(n.e), sd.e),
    sd.c = if_else(se_to_sd_corrected, sd.c * sqrt(n.c), sd.c)
  ) %>%
  make_unique() %>%
  add_dose()

print(psd_fixed %>% filter(!is.na(mean.e)), n = Inf)


# Network meta-analysis 
build_contrasts_bin <- function(d_wide) {
  d_wide %>%
    filter(!is.na(event.e), !is.na(n.e), !is.na(event.c), !is.na(n.c)) %>%
    mutate(treat1 = drug, treat2 = "Vehicle")
}

build_contrasts_cont <- function(d_wide) {
  d_wide %>%
    filter(!is.na(mean.e), !is.na(sd.e), !is.na(n.e),
           !is.na(mean.c), !is.na(sd.c), !is.na(n.c)) %>%
    mutate(treat1 = drug, treat2 = "Vehicle")
}

# NMA: IGA 
pw_iga <- pairwise(
  treat = list(treat1, treat2),
  event = list(event.e, event.c),
  n     = list(n.e, n.c),
  studlab = study_id,
  data    = build_contrasts_bin(iga),
  sm = "RR"
)
nma_iga <- netmeta(TE = TE, seTE = seTE,
                   treat1 = treat1, treat2 = treat2, studlab = studlab,
                   data = pw_iga, sm = "RR",
                   common = FALSE, random = TRUE,
                   reference.group = "Vehicle")
summary(nma_iga)

# NMA: Any AE 
pw_ae <- pairwise(
  treat = list(treat1, treat2),
  event = list(event.e, event.c),
  n     = list(n.e, n.c),
  studlab = study_id,
  data    = build_contrasts_bin(ae),
  sm = "RR"
)
nma_ae <- netmeta(TE = TE, seTE = seTE,
                  treat1 = treat1, treat2 = treat2, studlab = studlab,
                  data = pw_ae, sm = "RR",
                  common = FALSE, random = TRUE,
                  reference.group = "Vehicle")
summary(nma_ae)

# NMA: GI
pw_gi <- pairwise(
  treat = list(treat1, treat2),
  event = list(event.e, event.c),
  n     = list(n.e, n.c),
  studlab = study_id,
  data    = build_contrasts_bin(gi),
  sm = "RR", allstudies = TRUE
)
nma_gi <- netmeta(TE = TE, seTE = seTE,
                  treat1 = treat1, treat2 = treat2, studlab = studlab,
                  data = pw_gi, sm = "RR",
                  common = FALSE, random = TRUE,
                  reference.group = "Vehicle")
summary(nma_gi)

# NMA: Irritation
pw_irr <- pairwise(
  treat = list(treat1, treat2),
  event = list(event.e, event.c),
  n     = list(n.e, n.c),
  studlab = study_id,
  data    = build_contrasts_bin(irritation),
  sm = "RR", allstudies = TRUE
)
nma_irr <- netmeta(TE = TE, seTE = seTE,
                   treat1 = treat1, treat2 = treat2, studlab = studlab,
                   data = pw_irr, sm = "RR",
                   common = FALSE, random = TRUE,
                   reference.group = "Vehicle")
summary(nma_irr)

# NMA: PSD 
pw_psd <- pairwise(
  treat = list(treat1, treat2),
  mean  = list(mean.e, mean.c),
  sd    = list(sd.e, sd.c),
  n     = list(n.e, n.c),
  studlab = study_id,
  data    = build_contrasts_cont(psd_fixed),
  sm = "MD"
)
nma_psd <- netmeta(TE = TE, seTE = seTE,
                   treat1 = treat1, treat2 = treat2, studlab = studlab,
                   data = pw_psd, sm = "MD",
                   common = FALSE, random = TRUE,
                   reference.group = "Vehicle")
summary(nma_psd)


# Network plot and NMA forest plots
pdf("network_plot.pdf", width = 6, height = 6)
netgraph(nma_iga, plastic = FALSE,
         thickness = "number.of.studies",
         points = TRUE, cex.points = 5, cex = 1.2)
dev.off()

pdf("forest_nma_iga.pdf", width = 8, height = 4)
forest(nma_iga, reference.group = "Vehicle", sortvar = TE,
       smlab = "IGA/PGA success (RR)")
dev.off()

pdf("forest_nma_any_ae.pdf", width = 8, height = 4)
forest(nma_ae, reference.group = "Vehicle", sortvar = TE,
       smlab = "Any AE (RR)")
dev.off()



# Rankings
cat("\n=== Rankings: IGA/PGA success (higher = more efficacious) ===\n")
print(netrank(nma_iga, small.values = "undesirable"))

cat("\n=== Rankings: Any AE (higher = safer) ===\n")
print(netrank(nma_ae, small.values = "desirable"))

# Headline table — head-to-head NMA results
get_h2h <- function(nma, t1 = "Roflumilast", t2 = "Tapinarof", is_log = TRUE) {
  TE_val <- nma$TE.random[t1, t2]
  se_val <- nma$seTE.random[t1, t2]
  lo <- TE_val - 1.96 * se_val
  hi <- TE_val + 1.96 * se_val
  p  <- 2 * (1 - pnorm(abs(TE_val / se_val)))
  if (is_log) c(est = exp(TE_val), lo = exp(lo), hi = exp(hi), p = p)
  else        c(est = TE_val, lo = lo, hi = hi, p = p)
}

iga_h2h <- get_h2h(nma_iga)
ae_h2h  <- get_h2h(nma_ae)
gi_h2h  <- get_h2h(nma_gi)
irr_h2h <- get_h2h(nma_irr)
psd_h2h <- get_h2h(nma_psd, is_log = FALSE)

headline <- tibble(
  outcome  = c("IGA/PGA success", "Any AE", "GI AE", "Irritation AE", "PSD change"),
  measure  = c("RR", "RR", "RR", "RR", "MD"),
  estimate = c(iga_h2h["est"], ae_h2h["est"], gi_h2h["est"], irr_h2h["est"], psd_h2h["est"]),
  ci_lower = c(iga_h2h["lo"],  ae_h2h["lo"],  gi_h2h["lo"],  irr_h2h["lo"],  psd_h2h["lo"]),
  ci_upper = c(iga_h2h["hi"],  ae_h2h["hi"],  gi_h2h["hi"],  irr_h2h["hi"],  psd_h2h["hi"]),
  p_value  = c(iga_h2h["p"],   ae_h2h["p"],   gi_h2h["p"],   irr_h2h["p"],   psd_h2h["p"])
) %>%
  mutate(across(c(estimate, ci_lower, ci_upper), ~round(., 2)),
         p_value = round(p_value, 3))

cat("\n=== HEADLINE: Roflumilast vs Tapinarof (UPDATED data) ===\n")
print(headline)
write_csv(headline, "headline_results_v2.csv")

cat("\n========== Pipeline complete ==========\n")
cat("Files saved to:", getwd(), "\n")


#NMA forest plots with strength, k, and N labels
#Compute Ns per outcome to verify
cat("\n=== Ns per outcome (verify before using in labels) ===\n")
cat("\nIGA:\n")
iga %>%
  filter(!is.na(event.e), !is.na(n.e), !is.na(event.c), !is.na(n.c)) %>%
  group_by(drug) %>% summarise(k = n(), N = sum(n.e + n.c, na.rm = TRUE)) %>% print()
cat("\nAny AE:\n")
ae %>%
  filter(!is.na(event.e), !is.na(n.e), !is.na(event.c), !is.na(n.c)) %>%
  group_by(drug) %>% summarise(k = n(), N = sum(n.e + n.c, na.rm = TRUE)) %>% print()
cat("\nGI:\n")
gi %>%
  filter(!is.na(event.e), !is.na(n.e), !is.na(event.c), !is.na(n.c)) %>%
  group_by(drug) %>% summarise(k = n(), N = sum(n.e + n.c, na.rm = TRUE)) %>% print()
cat("\nIrritation:\n")
irritation %>%
  filter(!is.na(event.e), !is.na(n.e), !is.na(event.c), !is.na(n.c)) %>%
  group_by(drug) %>% summarise(k = n(), N = sum(n.e + n.c, na.rm = TRUE)) %>% print()
cat("\nPSD:\n")
psd_fixed %>%
  filter(!is.na(mean.e), !is.na(sd.e), !is.na(n.e),
         !is.na(mean.c), !is.na(sd.c), !is.na(n.c)) %>%
  group_by(drug) %>% summarise(k = n(), N = sum(n.e + n.c, na.rm = TRUE)) %>% print()

# Evidence text per treatment, keyed by name
kn_data <- function(nma, vals) {
  out <- unname(vals[nma$trts])
  out[is.na(out)] <- ""
  data.frame(kN = out, row.names = nma$trts, check.names = FALSE)
}

kn_iga <- c(Roflumilast = "(5 RCTs, N = 1701)", Tapinarof = "(6 RCTs, N = 1836)", Vehicle = "")
kn_ae  <- c(Roflumilast = "(5 RCTs, N = 1831)", Tapinarof = "(7 RCTs, N = 1856)", Vehicle = "")
kn_gi  <- c(Roflumilast = "(4 RCTs, N = 1615)", Tapinarof = "(5 RCTs, N = 1679)", Vehicle = "")
kn_irr <- c(Roflumilast = "(1 RCT, N = 302)",   Tapinarof = "(4 RCTs, N = 1603)", Vehicle = "")
kn_psd <- c(Roflumilast = "(4 RCTs, N = 1345)", Tapinarof = "(2 RCTs, N = 1025)", Vehicle = "")

trt_labels <- c("Roflumilast 0.3%", "Tapinarof 1%", "Vehicle")
draw_v2 <- function(nma, vals, file, smlab) {
  pdf(file, width = 10, height = 4)
  forest(nma, reference.group = "Vehicle", sortvar = TE,
         smlab     = smlab,
         labels    = trt_labels,
         add.data  = kn_data(nma, vals),
         leftcols  = c("studlab", "kN"),
         leftlabs  = c("Treatment", ""),      # name the column if you prefer
         rightcols = c("effect", "ci"),
         fs.study.labels   = 11,              # drug names
         fs.study          = 9,               # evidence column + RR + CI
         just.studlab      = "left",
         just.addcols.left = "left",
         spacing = 1.2)
  dev.off()
}

draw_v2(nma_iga, kn_iga, "forest_nma_iga_v2.pdf",        "IGA/PGA success (RR)")
draw_v2(nma_ae,  kn_ae,  "forest_nma_any_ae_v2.pdf",     "Any AE (RR)")
draw_v2(nma_gi,  kn_gi,  "forest_nma_gi.pdf",         "GI AE (RR)")
draw_v2(nma_irr, kn_irr, "forest_nma_irritation.pdf", "Irritation AE (RR)")
draw_v2(nma_psd, kn_psd, "forest_nma_psd.pdf",        "PSD change (MD)")

cat("\n========== Done ==========\n")
