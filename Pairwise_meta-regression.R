raw_dat <- raw_dat %>%
  rename(
    baseline_bsa_text = `What is baseline Body Surface Area (BSA)? Enter any information - sometimes this is a maximum, sometimes a mean, sometimes categories (i.e. mild, medium, or severe).`,
    patient_itch_mean_int_cont = `patient itch_mean_int_cont`,
    patient_itch_n_int_categ = `patient itch_n_int_categ`
  ) %>%
  mutate(
    asian_group = if_else(
      `% Asian_comment` == 100,
      "Asian Country",
      "Western Country",
      missing = "Western Country"
    ),
    asian_group = factor(asian_group, levels = c("Asian Country", "Western Country"))
  ) %>%
  mutate(
    vehicle = coalesce(
      `What "vehicle" or inactive carrier is used with the active ingredient? -> Cream`,
      `What "vehicle" or inactive carrier is used with the active ingredient? -> Foam`
    )
  )

prep_bsa_dat <- function(raw_dat) {
  
  extract_num <- function(x, pattern) {
    as.numeric(str_match(x, regex(pattern, ignore_case = TRUE))[, 2])
  }
  
  raw_dat %>%
    select(
      ID,
      Refid,
      intervention_description,
      baseline_bsa_text,
      Roflumilast,
      Tapinarof,
      control,
      everything()
    ) %>%
    mutate(
      concentration = str_extract(intervention_description, "\\b\\d+\\.?\\d*%"),
      
      bsa_int = extract_num(
        baseline_bsa_text,
        "intervention[^0-9]*([0-9]+\\.?[0-9]*)"
      ),
      bsa_int_sd = extract_num(
        baseline_bsa_text,
        "intervention[^0-9]*[0-9]+\\.?[0-9]*\\s*\\(([0-9]+\\.?[0-9]*)\\)"
      ),
      
      bsa_ctrl = extract_num(
        baseline_bsa_text,
        "(?:vehicle|control|comparator)[^0-9]*([0-9]+\\.?[0-9]*)"
      ),
      bsa_ctrl_sd = extract_num(
        baseline_bsa_text,
        "(?:vehicle|control|comparator)[^0-9]*[0-9]+\\.?[0-9]*\\s*\\(([0-9]+\\.?[0-9]*)\\)"
      ),
      
      overall_bsa = extract_num(
        baseline_bsa_text,
        "mean bsa[^0-9]*([0-9]+\\.?[0-9]*)"
      ),
      overall_bsa_sd = extract_num(
        baseline_bsa_text,
        "mean bsa[^0-9]*[0-9]+\\.?[0-9]*\\s*\\(([0-9]+\\.?[0-9]*)\\)"
      ),
      
      bsa_int = coalesce(bsa_int, overall_bsa),
      bsa_ctrl = coalesce(bsa_ctrl, overall_bsa),
      bsa_int_sd = coalesce(bsa_int_sd, overall_bsa_sd),
      bsa_ctrl_sd = coalesce(bsa_ctrl_sd, overall_bsa_sd)
    ) %>%
    select(-overall_bsa, -overall_bsa_sd)
}
sev_dat <- prep_bsa_dat(raw_dat)

swap_comp_to_int <- function(dat, refid = 3434) {
  
  comp_cols <- names(dat)[grepl("_comp_", names(dat))]
  
  for (comp_col in comp_cols) {
    
    int_col <- sub("_comp_", "_int_", comp_col)
    
    if (int_col %in% names(dat)) {
      
      dat[[int_col]] <- ifelse(
        dat$Refid == refid,
        dat[[comp_col]],
        dat[[int_col]]
      )
      
    }
  }
  
  dat
}
sev_dat <- swap_comp_to_int(sev_dat, refid = 3434)
sev_dat <- sev_dat %>%
  mutate(across(contains("_cont"), as.numeric))

run_bsa_metareg <- function(dat,
                            drug_col,
                            event_e_col,
                            n_e_col,
                            event_c_col,
                            n_c_col,
                            asian_subgroup = FALSE,
                            vehicle_subgroup = FALSE) {
  dat2 <- dat %>%
    filter(!is.na(.data[[drug_col]])) %>%
    filter(
      !is.na(.data[[event_e_col]]),
      !is.na(.data[[n_e_col]]),
      !is.na(.data[[event_c_col]]),
      !is.na(.data[[n_c_col]])
    ) %>%
    mutate(
      across(
        all_of(c(event_e_col, n_e_col, event_c_col, n_c_col)),
        as.numeric
      ),
      
      baseline_bsa =
        (bsa_int * .data[[n_e_col]] +
           bsa_ctrl * .data[[n_c_col]]) /
        (.data[[n_e_col]] + .data[[n_c_col]]),
      
      subgroup_eth = factor(
        if_else(
          `% Asian_comment` == 100,
          "East Asian Country",
          "Western Country",
          missing = "Western Country"
        ),
        levels = c("East Asian Country", "Western Country")
      ),
      
      subgroup_vehicle = factor(
        vehicle,
        levels = c("Foam", "Cream")
      ),
      
      total_n = .data[[n_e_col]] + .data[[n_c_col]],
      
      sample_size = sprintf(
        "%d (%.1f)",
        total_n,
        100 * total_n / sum(total_n)
      )
    )
  
  dat2$baseline_bsa[is.nan(dat2$baseline_bsa)] <- NA_real_
  
  subgroup_var <- NULL
  subgroup_name <- ""
  
  if (asian_subgroup) {
    subgroup_var <- dat2$subgroup_eth
  } else if (vehicle_subgroup) {
    subgroup_var <- dat2$subgroup_vehicle
  }
  
  meta_sub <- meta::metabin(
    event.e = dat2[[event_e_col]],
    n.e     = dat2[[n_e_col]],
    event.c = dat2[[event_c_col]],
    n.c     = dat2[[n_c_col]],
    studlab = dat2$ID,
    data    = dat2,
    sm      = "RR",
    random  = TRUE,
    common  = FALSE,
    method.random.ci = "HK",
    subgroup = subgroup_var,
    subgroup.name = subgroup_name
  )
  
  meta_nosub <- meta::metabin(
    event.e = dat2[[event_e_col]],
    n.e     = dat2[[n_e_col]],
    event.c = dat2[[event_c_col]],
    n.c     = dat2[[n_c_col]],
    studlab = dat2$ID,
    data    = dat2,
    sm      = "RR",
    random  = TRUE,
    common  = FALSE,
    method.random.ci = "HK"
  )
  
  mr <- meta::metareg(meta_sub, ~ baseline_bsa)
  
  list(
    data = dat2,
    meta = meta_sub,
    meta_nosub = meta_nosub,
    metareg = mr
  )
}
plot_metareg <- function(res, vehicle_subgroup = FALSE) {
  
  if (vehicle_subgroup) {
    leftcols <- c("studlab", "sample_size")
    leftlabs <- c("Study", "Sample size (%)")
  } else {
    leftcols <- c("studlab", "vehicle", "sample_size")
    leftlabs <- c("Study", "Vehicle", "Sample size (%)")
  }
  
  forest(
    res$meta,
    studlab = TRUE,
    leftcols = leftcols,
    leftlabs = leftlabs,
    rightcols = c("effect", "ci"),
    print.subgroup.name = FALSE
  )
  
  summary(res$meta)
  print(coef(summary(res$metareg)))
  
  bubble(
    res$metareg,
    studlab = res$data$Refid,
    xlab = "Mean baseline BSA (%)",
    ylab = "Log RR"
  )
}

#C.S. Roflumilast----
res_R_cs <- run_bsa_metareg(
  sev_dat,
  drug_col    = "Roflumilast",
  event_e_col = "clinical_severity_counts_int_categ",
  n_e_col     = "clinical_severity_n_int_categ",
  event_c_col = "clinical_severity_counts_ctrl_categ",
  n_c_col     = "clinical_severity_n_ctrl_categ",
  vehicle_subgroup = TRUE
)

plot_metareg(res_R_cs,vehicle_subgroup = TRUE)

#C.S. Tapinarof----
res_T_cs <- run_bsa_metareg(
  sev_dat,
  drug_col    = "Tapinarof",
  event_e_col = "clinical_severity_counts_int_categ",
  n_e_col     = "clinical_severity_n_int_categ",
  event_c_col = "clinical_severity_counts_ctrl_categ",
  n_c_col     = "clinical_severity_n_ctrl_categ",
  asian_subgroup = TRUE
)
plot_metareg(res_T_cs)
summary(res_T_cs$meta)
#A.E. Roflumilast----
res_R_ae <- run_bsa_metareg(
  sev_dat,
  drug_col    = "Roflumilast",
  event_e_col = "nAEparticipants_counts_int_categ",
  n_e_col     = "nAEparticipants_n_int_categ",
  event_c_col = "nAEparticipants_counts_ctrl_categ",
  n_c_col     = "nAEparticipants_n_ctrl_categ",
  vehicle_subgroup = TRUE
)
plot_metareg(res_R_ae,vehicle_subgroup = TRUE)
summary(res_R_ae$meta)
#A.E. Tapinarof----
res_T_ae <- run_bsa_metareg(
  sev_dat,
  drug_col    = "Tapinarof",
  event_e_col = "nAEparticipants_counts_int_categ",
  n_e_col     = "nAEparticipants_n_int_categ",
  event_c_col = "nAEparticipants_counts_ctrl_categ",
  n_c_col     = "nAEparticipants_n_ctrl_categ",
  asian_subgroup = TRUE
)
plot_metareg(res_T_ae)


#Itch Roflumilast  binary----
res_R_itch <- run_bsa_metareg(
  sev_dat,
  drug_col   = "Roflumilast",
  event_e_col = "patient_itch_counts_int_categ",
  n_e_col     = "patient_itch_n_int_categ",
  event_c_col = "patient_itch_counts_ctrl_categ",
  n_c_col     = "patient_itch_n_ctrl_categ",
  vehicle_subgroup = TRUE
)

plot_metareg(res_R_itch,  vehicle_subgroup = TRUE)
#Itch Tapinarof itch continuous----
tap_itch <- sev_dat %>%
  filter(!is.na(Tapinarof)) %>%
  filter(
    !is.na(patient_itch_mean_int_cont),
    !is.na(patient_itch_SD_int_cont),
    !is.na(patient_itch_n_int_cont),
    !is.na(patient_itch_mean_ctrl_cont),
    !is.na(patient_itch_SD_ctrl_cont),
    !is.na(patient_itch_n_ctrl_cont)
  ) %>%
  mutate(
    across(
      c(
        patient_itch_n_int_cont,
        patient_itch_mean_int_cont,
        patient_itch_SD_int_cont,
        patient_itch_n_ctrl_cont,
        patient_itch_mean_ctrl_cont,
        patient_itch_SD_ctrl_cont
      ),
      as.numeric
    ),
    baseline_bsa = (
      bsa_int * patient_itch_n_int_cont +
        bsa_ctrl * patient_itch_n_ctrl_cont
    ) / (
      patient_itch_n_int_cont + patient_itch_n_ctrl_cont
    ),
    total_n = patient_itch_n_int_cont + patient_itch_n_ctrl_cont,
    sample_size = sprintf(
      "%d (%.1f)",
      total_n,
      100 * total_n / sum(total_n)
    )
  )

tap_itch$baseline_bsa[is.nan(tap_itch$baseline_bsa)] <- NA_real_

m_tap_itch <- metacont(
  n.e = patient_itch_n_int_cont,
  mean.e = patient_itch_mean_int_cont,
  sd.e = patient_itch_SD_int_cont,
  n.c = patient_itch_n_ctrl_cont,
  mean.c = patient_itch_mean_ctrl_cont,
  sd.c = patient_itch_SD_ctrl_cont,
  studlab = ID,
  data = tap_itch,
  sm = "SMD",
  random = TRUE,
  common = FALSE,
  method.random.ci = "HK"
)

summary(m_tap_itch)

forest(
  m_tap_itch,
  leftcols = c("studlab", "vehicle", "sample_size"),
  leftlabs = c("Study", "Vehicle", "Sample size (%)"),
  rightcols = c("effect", "ci")
)

mr_tap_itch <- metareg(
  m_tap_itch,
  ~ baseline_bsa
)

summary(mr_tap_itch)

bubble(
  mr_tap_itch,
  studlab = mr_tap_itch$data$Refid,
  xlab = "Mean baseline BSA (%)",
  ylab = "SMD"
)
#folliculitis Tapinarof----
res_T_f <- run_bsa_metareg(
  sev_dat,
  drug_col    = "Tapinarof",
  event_e_col = "folliculitis_counts_int_categ",
  n_e_col     = "folliculitis_n_int_categ",
  event_c_col = "folliculitis_counts_ctrl_categ",
  n_c_col     = "folliculitis_n_ctrl_categ",
  asian_subgroup = TRUE
)

plot_metareg(res_T_f)



#PSD Roflumilast----
rof_symptoms <- sev_dat %>%
  filter(!is.na(Roflumilast)) %>%
  filter(
    !is.na(patient_symptoms_n_int_cont),
    !is.na(patient_symptoms_mean_int_cont),
    !is.na(patient_symptoms_SD_int_cont),
    !is.na(patient_symptoms_n_ctrl_cont),
    !is.na(patient_symptoms_mean_ctrl_cont),
    !is.na(patient_symptoms_SD_ctrl_cont)
  ) %>%
  mutate(
    across(
      c(
        patient_symptoms_n_int_cont,
        patient_symptoms_mean_int_cont,
        patient_symptoms_SD_int_cont,
        patient_symptoms_n_ctrl_cont,
        patient_symptoms_mean_ctrl_cont,
        patient_symptoms_SD_ctrl_cont
      ),
      as.numeric
    ),
    baseline_bsa = (
      bsa_int * patient_symptoms_n_int_cont +
        bsa_ctrl * patient_symptoms_n_ctrl_cont
    ) / (
      patient_symptoms_n_int_cont + patient_symptoms_n_ctrl_cont
    ),
    subgroup_vehicle = factor(
      vehicle,
      levels = c("Foam", "Cream")
    ),
    total_n = patient_symptoms_n_int_cont + patient_symptoms_n_ctrl_cont,
    sample_size = sprintf(
      "%d (%.1f)",
      total_n,
      100 * total_n / sum(total_n)
    )
  )

rof_symptoms$baseline_bsa[is.nan(rof_symptoms$baseline_bsa)] <- NA_real_

m_rof_symptoms <- metacont(
  n.e = patient_symptoms_n_int_cont,
  mean.e = patient_symptoms_mean_int_cont,
  sd.e = patient_symptoms_SD_int_cont,
  n.c = patient_symptoms_n_ctrl_cont,
  mean.c = patient_symptoms_mean_ctrl_cont,
  sd.c = patient_symptoms_SD_ctrl_cont,
  studlab = ID,
  data = rof_symptoms,
  sm = "SMD",
  random = TRUE,
  common = FALSE,
  method.random.ci = "HK",
  subgroup = subgroup_vehicle
)

m_rof_symptoms_nosub <- metacont(
  n.e = patient_symptoms_n_int_cont,
  mean.e = patient_symptoms_mean_int_cont,
  sd.e = patient_symptoms_SD_int_cont,
  n.c = patient_symptoms_n_ctrl_cont,
  mean.c = patient_symptoms_mean_ctrl_cont,
  sd.c = patient_symptoms_SD_ctrl_cont,
  studlab = ID,
  data = rof_symptoms,
  sm = "SMD",
  random = TRUE,
  common = FALSE,
  method.random.ci = "HK"
)

summary(m_rof_symptoms)

forest(
  m_rof_symptoms,
  leftcols = c("studlab", "sample_size"),
  leftlabs = c("Study", "Sample size (%)"),
  rightcols = c("effect", "ci"),
  print.subgroup.name = FALSE
)

mr_rof_symptoms <- metareg(
  m_rof_symptoms,
  ~ baseline_bsa
)

summary(mr_rof_symptoms)

bubble(
  mr_rof_symptoms,
  studlab = mr_rof_symptoms$data$Refid,
  xlab = "Mean baseline BSA (%)",
  ylab = "SMD"
)
#Bias test----
run_bias_tests <- function(x, outcome_name, k_min = 3) {
  stopifnot(inherits(x, "meta"))
  
  cat("\n==============================\n")
  cat(outcome_name, "\n")
  cat("==============================\n\n")
  
  cat("Begg's test:\n")
  print(meta::metabias(x, method.bias = "Begg", k.min = k_min))
  
  cat("\nEgger's test:\n")
  print(meta::metabias(x, method.bias = "Egger", k.min = k_min))
}

run_bias_tests(res_R_cs$meta_nosub, "Clinical Severity - Roflumilast")
run_bias_tests(res_T_cs$meta_nosub, "Clinical Severity - Tapinarof")

run_bias_tests(res_R_ae$meta_nosub, "Adverse Events - Roflumilast")
run_bias_tests(res_T_ae$meta_nosub, "Adverse Events - Tapinarof")

run_bias_tests(res_R_itch$meta_nosub, "Itch - Roflumilast")
run_bias_tests(m_tap_itch, "Itch - Tapinarof")

run_bias_tests(res_T_f$meta_nosub, "Folliculitis - Tapinarof")

run_bias_tests(m_rof_symptoms_nosub, "Patient Diary Symptoms - Roflumilast")
