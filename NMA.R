library(dplyr)
library(readr)
library(stringr)
library(metafor)
library(netmeta)
library(tidyr)
library(ggplot2)
raw_dat <- read_csv("C:/Users/10531/Downloads/distillersr-PCORI_Psoriasis KQ1 +2_2026-05-10 (1).csv")

col_df <- data.frame(
  column_number = seq_along(names(raw_dat)),
  column_name = names(raw_dat)
)

raw_dat <- raw_dat %>%
  rename(
    baseline_bsa_text = `What is baseline Body Surface Area (BSA)? Enter any information - sometimes this is a maximum, sometimes a mean, sometimes categories (i.e. mild, medium, or severe).`
  )

extract_after_keyword <- function(x, keyword) {
  m <- str_match(
    x,
    regex(
      paste0(keyword, "[^0-9]*([0-9]+\\.?[0-9]*)"),
      ignore_case = TRUE
    )
  )
  as.numeric(m[, 2])
}

# Helper: extract a single overall mean BSA
extract_mean_bsa <- function(x) {
  m <- str_match(
    x,
    regex("mean bsa[^0-9]*([0-9]+\\.?[0-9]*)", ignore_case = TRUE)
  )
  as.numeric(m[, 2])
}

sev_dat <- raw_dat %>%
  select(
    ID,
    clinical_severity,
    clinical_severity_measure_categ,
    clinical_severity_measure_cont,
    clinical_severity_n_int_categ,
    clinical_severity_counts_int_categ,
    clinical_severity_n_ctrl_categ,
    clinical_severity_counts_ctrl_categ,
    clinical_severity_n_comp_categ,
    clinical_severity_counts_comp_categ,
    followup,
    intervention_description,
    baseline_bsa_text,
    Roflumilast,
    Tapinarof,
    control
  ) %>%
  mutate(
    concentration = str_extract(
      intervention_description,
      "\\b\\d+\\.?\\d*%"
    ),
    
    bsa_int = as.numeric(
      str_match(
        baseline_bsa_text,
        regex("intervention[^0-9]*([0-9]+\\.?[0-9]*)",
              ignore_case = TRUE)
      )[,2]
    ),
    
    bsa_ctrl = as.numeric(
      str_match(
        baseline_bsa_text,
        regex("(?:vehicle|control|comparator)[^0-9]*([0-9]+\\.?[0-9]*)",
              ignore_case = TRUE)
      )[,2]
    ),
    
    overall_bsa = as.numeric(
      str_match(
        baseline_bsa_text,
        regex("mean bsa[^0-9]*([0-9]+\\.?[0-9]*)",
              ignore_case = TRUE)
      )[,2]
    ),
    
    # if only one mean reported, copy to both arms
    bsa_int = coalesce(bsa_int, overall_bsa),
    bsa_ctrl = coalesce(bsa_ctrl, overall_bsa)
  ) %>%
  select(-overall_bsa)


col_df <- data.frame(
  column_number = seq_along(names(sev_dat)),
  column_name = names(sev_dat)
)

sev_es <- escalc(
  measure = "RR",
  ai = clinical_severity_counts_int_categ,
  bi = clinical_severity_n_int_categ - clinical_severity_counts_int_categ,
  ci = clinical_severity_counts_ctrl_categ,
  di = clinical_severity_n_ctrl_categ - clinical_severity_counts_ctrl_categ,
  data = sev_dat
)

sev_dat2 <- sev_dat %>%
  mutate(
    treat = coalesce(Roflumilast, Tapinarof)
  )

arm_dat <- bind_rows(
  sev_dat2 %>%
    transmute(
      study = ID,
      treat = treat,
      event = clinical_severity_counts_int_categ,
      n = clinical_severity_n_int_categ,
      baseline_bsa = bsa_int
    ),
  sev_dat2 %>%
    transmute(
      study = ID,
      treat = control,
      event = clinical_severity_counts_ctrl_categ,
      n = clinical_severity_n_ctrl_categ,
      baseline_bsa = bsa_ctrl
    )
) %>%
  filter(!is.na(study), !is.na(treat), !is.na(event), !is.na(n)) %>%
  mutate(
    treat = str_trim(as.character(treat)),
    study = as.character(study)
  )

#NMA----
pw <- pairwise(
  treat = treat,
  event = event,
  n = n,
  studlab = study,
  data = arm_dat,
  sm = "RR"
)

net1 <- netmeta(
  TE = pw$TE,
  seTE = pw$seTE,
  treat1 = pw$treat1,
  treat2 = pw$treat2,
  studlab = pw$studlab,
  sm = "RR",
  random = TRUE,
  reference.group = "No intervention"
)

#NMA output
summary(net1)
forest(
  net1,
  ref = "No intervention",
  text.random = "RR (Random Effects)",
  smlab = "Treatment vs No intervention"
)
netgraph(net1)
netrank(net1)
netleague(net1)