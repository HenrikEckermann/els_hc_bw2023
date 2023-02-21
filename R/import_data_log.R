library(mice)
library(tidyverse)
library(here)
library(glue)


# set seed for imputation later
set.seed(1)
df <- readxl::read_excel(here("data/Dataforanalysis Argentina FINAL.xlsx"))

# rename/prepare variables for analysis
df <- df %>%
  select(
    id = DEM.PNUM,
    birthweight = MD.WEIGHTF,
    head_cf = MD.HEADF,
    hc_0_20 = "Cortisol week 0-20 (pg/mg)",
    hc_21_40 = "Cortisol week 21-40 (pg/mg)",
    age_mother = "Age at participation",
    length_mother = DEM.LENGTHM,
    weight_mother = DEM.WEIGHTM,
    smoke_preg = DEM.SMOKE,
    alcohol_preg = DEM.ALCOHOL,
    edu = DEM.EDUM,
    n_inhab = DEM.LIVINGP,
    n_rooms = DEM.LIVINGR,
    life_ev_sum = LEI.sev,
    life_ev_freq = LEI.freq,
    jtv_sum = JTV.sum,
    ace_sum = ACE.sum,
    neo_female = MD.SEXF,
    gest_age = DM.GEST,
    high_bp_preg = MD.HBP,
    threat_premat = MD.TPB,
    csection = MD.TYPEBIRTH,
    haircolor = PELO.1,
    hairtreat_6m = PELO.2,
    hairwashes_per_week = PELO.4,
    hairproducts = PELO.6,
    corticost = PELO.8,
    life_event_3m = PELO.11
  )


# note missing values as NA
df <- df %>% mutate_all(function(value) ifelse(value == 99999, NA, value)) %>%
  mutate(
    csection = as.factor(csection),
    threat_premat = as.factor(threat_premat),
    high_bp_preg = as.factor(high_bp_preg),
    neo_female = as.factor(neo_female),
    alcohol_preg = as.factor(alcohol_preg),
    smoke_preg = as.factor(smoke_preg),
    hairproducts = as.factor(hairproducts),
    corticost = as.factor(corticost),
    life_event_3m = as.factor(life_event_3m)
    )
# code factors
fac_list <- list(
  "csection", "threat_premat", "high_bp_preg", "neo_female", 
  "alcohol_preg", "smoke_preg", "hairproducts", "corticost", "life_event_3m")
map(fac_list, ~levels(df[[.x]]))

# select variables for analyses
predictors <- c(
  "life_ev_sum" ,
  "life_ev_freq",
  "jtv_sum",
  "ace_sum",
  "age_mother",
  "length_mother",
  "weight_mother",
  "edu",
  "neo_female",
  "n_inhab",
  "n_rooms"
)
haircortisol <- c(
  "hc_0_20",
  "hc_21_40"
)

hc_covariates <- c(
  "haircolor",
  "hairtreat_6m",
  "hairwashes_per_week",
  "hairproducts",
  "corticost",
  "life_event_3m"
)
outcomes <- c("head_cf", "birthweight", "gest_age")

df_hc <- select(df, id, all_of(haircortisol), all_of(hc_covariates))
df <- select(df, id, all_of(predictors), all_of(haircortisol), all_of(outcomes), all_of(hc_covariates))

mlr::summarizeColumns(df)
missing_ids <- filter(df, is.na(hc_21_40)) %>% .$id

# inspect edu
df <- mutate(df, edu = as.integer(as.factor(edu)))
count(df, edu)
ggplot(df, aes(edu)) +
  geom_density()

ggplot(df, aes(edu, life_ev_sum)) +
  geom_jitter() +
  geom_smooth()

# 69 total but 45 complete cases
dim(df)
dim(na.omit(df))
unique(df$id) %>% length()

df_mis <- as.data.frame(df)

# change hw to categorical
count(df_mis, hairwashes_per_week)
df_mis <- df_mis %>%
            mutate(
              hw = as.factor(ifelse(hairwashes_per_week < 4, "<4",
                                ifelse(hairwashes_per_week >=4, ">4", NA))),
                   haircolor = as.factor(ifelse(haircolor == 5, 3, haircolor)),
                   hairtreat_6m = as.factor(
                          ifelse(hairtreat_6m == 1, 0,
                            ifelse(hairtreat_6m > 1 & hairtreat_6m < 7, 1, NA))
                          ),
                    across(matches("hc_\\d+_\\d+"), function(x) log(x)),
              id = as.factor(id)
                              ) %>%
          select(-hairwashes_per_week)
ggplot(df_mis, aes(as.factor(life_event_3m), hc_0_20)) +
  geom_point()

count(df_mis, life_event_3m)

# I need this df in order to scale the two hc vars together
hctemp <- select(df_mis, contains("hc"), id) %>%
  pivot_longer(
    -id,
    names_to = "varname",
    values_to = "hc") %>%
    mutate(hc = scale(hc)[, 1]) %>%
    pivot_wider(names_from = "varname", values_from = "hc")

d <- df_mis %>%
     select(-contains("hc_")) %>% # we want to scale them together
     mutate(els = scale(jtv_sum)[, 1] + scale(ace_sum)[, 1]) %>%
     mutate(across(where(is.numeric), function(x) scale(x)[, 1])) %>%
     left_join(hctemp, by = "id") %>%
     select(
      lei = life_ev_sum,
      cc = corticost,
      age = age_mother,
      hc1 = hc_0_20,
      hc2 = hc_21_40,
      length = length_mother,
      weight = weight_mother,
      ga = gest_age,
      cf = head_cf,
      bw = birthweight,
      everything(),
      -life_ev_freq
  ) %>%
  select(-id)

# edu treated differently here so that imputation as ordinal works below
d2 <- df_mis %>%
     select(-contains("hc_")) %>% # we want to scale them together
     mutate(
       els = scale(jtv_sum)[, 1] + scale(ace_sum)[, 1],
       edu = as.factor(edu)
     ) %>%
     mutate(across(where(is.numeric), function(x) scale(x)[, 1])) %>%
     left_join(hctemp, by = "id") %>%
     select(
      lei = life_ev_sum,
      cc = corticost,
      age = age_mother,
      hc1 = hc_0_20,
      hc2 = hc_21_40,
      length = length_mother,
      weight = weight_mother,
      ga = gest_age,
      cf = head_cf,
      bw = birthweight,
      everything(),
      -life_ev_freq
  ) %>% mutate(edu = as.integer(edu)) %>%
  select(-id)


if (!file.exists(here::here("data/data_cl_log_sc.Rds"))) {
  # I inspected imputations
  df_imp <- mice(d, m = 50, method = "pmm")
  # for the ordinal model I need edu as non-standardized
  df_imp_edu <- mice(d2, m = 50, method = "pmm")

  save(
    df,
    df_imp,
    df_hc,
    missing_ids,
    df_imp_edu,
    file = here::here("data/data_cl_log_sc.Rds")
  )
 } else {
  load(here::here("data/data_cl_log_sc.Rds"))
}


# descriptive statistics table
dftemp <- df %>%
  mutate(
    haircolor = as.factor(ifelse(haircolor == 1, "black", ifelse(haircolor == 2, "brown", ifelse(haircolor == 3, "blond", ifelse(haircolor == 4, "red", ifelse(haircolor == 5,"grey", NA)))))),
    hc_0_20 = log(hc_0_20),
    hc_21_40 = log(hc_21_40),
    ci = head_cf/birthweight
  )
tbl <- mlr::summarizeColumns(dftemp)
new_name <- c(
  "empty",
  "Life Events",
  "empty",
  "Childhood Trauma Questionnaire",
  "Adverse Childhood Events",
  "Hair Cortisol 0 - 20 wks",
  "Hair Cortisol 21 - 40 wks",
  "Maternal Age",
  "Maternal Height",
  "Maternal Weight",
  "Educational Level",
  "empty",
  "Gestational Age",
  "empty",
  "empty",
  "Head Circumference",
  "Birth weight",
  "empty",
  "empty",
  "Hair washes per Week",
  "empty",
  "empty",
  "empty",
  "Cephalization index "
  )
tbl$name <- new_name
tbl <- filter(tbl, name != "empty") %>%
 select(Variable = name, Missing = na, Mean = mean, SD = disp, Median = median, MAD = mad, Min = min, Max = max) %>%
 mutate(across(where(is.numeric), function(x) round(x, 2))) %>%
 arrange(Variable)

# Carolina requested a correlation table:
cor_tbl <- mutate(dftemp, across(everything(), function(x) as.numeric(x))) %>%
  select(-c(
    id, life_ev_freq, neo_female, n_inhab, n_rooms, haircolor, hairtreat_6m, 
    hairproducts, corticost, life_event_3m)
    ) 

new_name <- c(
  "Life Events",
  "Childhood Trauma Questionnaire",
  "Adverse Childhood Events",
  "Hair Cortisol 0 - 20 wks",
  "Hair Cortisol 21 - 40 wks",
  "Maternal Age",
  "Maternal Height",
  "Maternal Weight",
  "Educational Level",
  "Gestational Age",
  "Head Circumference",
  "Birth weight",
  "Hair washes per Week",
  "Cephalization index "
  )

colnames(cor_tbl) <- new_name
cor_tbl <- corx::corx(
  as.data.frame(cor_tbl),
  triangle = "lower",
  stars = c(0.05, 0.01, 0.001),
  describe = NULL
  )

save(tbl, cor_tbl, file = here::here("data/descr_tbl.Rds"))



tbl
df_mis$head_cf
df_mis %>%
  mutate(
    els = scale(ace_sum)[, 1] + scale(jtv_sum)[, 1],
    lei = scale(life_ev_sum)[, 1],
    ci = head_cf/birthweight
  ) %>%
  ggplot(aes(hc_21_40, ci)) +
  geom_point() +
  geom_smooth()
df_mis %>%
  mutate(
    els = scale(ace_sum)[, 1] + scale(jtv_sum)[, 1],
    lei = scale(life_ev_sum)[, 1],
    ci = head_cf/birthweight
  ) %>%
  ggplot(aes(lei, hc_0_20)) +
  geom_point() +
  geom_smooth()
df_mis %>%
  mutate(
    els = scale(ace_sum)[, 1] + scale(jtv_sum)[, 1],
    lei = scale(life_ev_sum)[, 1],
    ci = head_cf/birthweight
  ) %>%
  ggplot(aes(log(hc_0_20), log(hc_21_40))) +
  geom_point() 
mutate(df_mis, across(where(is.numeric), function(x) scale(x)[, 1])) %>% filter(corticost == 1)
ggplot(df, aes(log(hc_0_20))) +
  geom_density()



# Fig 3
#remotes::install_github('jorvlan/raincloudplots')
library(raincloudplots)

dfcloud <- df %>%
  mutate(id = 1:dim(df_mis)[1]) %>%
  mutate(across(contains("hc_"), function(x) log(x)))

df_1x1 <- data_1x1(
  array_1 = dfcloud$hc_0_20,
  array_2 = dfcloud$hc_21_40
)

f3 <- raincloud_1x1_repmes(
  data = df_1x1,
  colors = (c('#C00000', '#C00000')),
  fills = (c('#C00000', '#C00000')),
  line_color = 'gray',
  line_alpha = .5,
  size = 5,
  alpha = .6,
  align_clouds = FALSE) +
  scale_x_continuous(breaks=c(1,2), labels=c("1", "2"), limits=c(0, 3)) +
  xlab("Sample") +
  ylab("Hair Cortisol (log)") +
  theme_bw(base_size = 20)
f3
save(f3, file = here::here("data/figures.Rds"))
