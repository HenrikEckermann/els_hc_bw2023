options(repr.plot.width = 15, repr.plot.height = 15, repr.plot.res = 100)

library(tidyverse)
library(dagitty)
library(brms)


load(here::here("data/data_cl_log_sc.Rds"))


# For data driven model structuring I will first proceed with one example dataset
# and in case of close decisions during
# CI I will test several datasets. For model fitting I will use them all.
d <- mice::complete(df_imp, 1)


g <- dagitty("dag{
    {hw cc} -> hc1 <- els -> lei-> {ga bw cf} <- hc1
    {hw cc} -> hc2 <- els -> lei-> {ga bw cf} <- hc2
    els -> {ga bw cf}
    lei -> {hc1 -> hc2}
    els -> {length -> weight} -> {ga bw cf hc1 hc2}
    age -> {hc1 hc2 ga bw <-> cf lei hw edu}
    ga -> {bw cf}
    weight -> lei
    els -> edu -> {lei weight hc1 hc2 bw ga cf}
    {els lei edu age} -> cc -> {weight bw ga cf}
}")
coordinates(g) <- list(
  x = c(
    els = 0, hc1 = 0.75, hc2 = 1.25, lei = 1,
    hw = 1, cc = 0.75, ga = 2, bw = 2, cf = 2, length = 0.25, weight = 0.25, age = 0.35,
    edu = 0.25, nishta = 2.5),
  y = c(
    els = 0, lei = 1, hc1 = -1.25, hc2 = -1.25, hw = -1.5, cc = -1.5,
    bw = -0.75, ga = 0, cf = 0.75, length = -0.75, weight = -0.5, age = -1,
    edu = 0.5)
)

rethinking::drawdag(g)

# for the publication I use upper cases
g_uc <- dagitty("dag{
  {HW CC} -> HC1 <- ELS -> LEI-> {GA BW CI} <- HC1
  {HW CC} -> HC2 <- ELS -> LEI-> {GA BW CI} <- HC2
  ELS -> {GA BW CI}
  LEI -> {HC1 -> HC2}
  ELS -> {LENGTH -> WEIGHT} -> {GA BW CI HC1 HC2}
  AGE -> {HC1 HC2 GA BW <-> CI LEI HW EDU}
  GA -> {BW CI}
  WEIGHT -> LEI
  ELS -> EDU -> {LEI WEIGHT HC1 HC2 BW GA CI}
  {ELS LEI EDU AGE} -> CC -> {WEIGHT BW GA CI}
  nishta}")
  
coordinates(g_uc) <- list(
  x = c(
    ELS = -0.25, HC1 = 0.75, HC2 = 1.25, LEI = 1,
    HW = 1, CC = 0.75, GA = 2, BW = 2, CI = 2, LENGTH = 0.23, WEIGHT = 0.23, AGE = 0.35, EDU = 0.23, nishta = 2.5),
  y = c(
    ELS = 0, LEI = 1, HC1 = -1.25, HC2 = -1.25, HW = -1.5, CC = -1.5,
    BW = -0.75, GA = 0, CI = 0.75, LENGTH = -0.75, WEIGHT = -0.5, AGE = -1,
    EDU = 0.5, nishta = 0))

rethinking::drawdag(g_uc)
impliedConditionalIndependencies(g)


# lets first make sure that the data agrees with above dag
# age _||_ edu

# we need to add edge between age and edu as there are also younger women

# cc _||_ hw | age
summary(glm(cc ~ hw + age, d, family = binomial()))
# edu _||_ hw
summary(lm(edu ~ hw + age, d))
# edu _||_ lngth | els
summary(lm(edu ~ length + els, d))
# age _||_ els
cor.test(d$age, d$els)
# age _||_ hw
t.test(age ~ hw, d)
ggplot(d, aes(hw, age)) +
  geom_jitter(width = 0.1)
# added an arrow from age to hw, younger women wash more often
# age _||_ lngt
cor.test(d$age, d$length)
# age _||_ wght
cor.test(d$age, d$weight)
# bw _||_ hw | age, cc, els, hc_1, hc_2, lei, lngt, wght, edu
summary(lm(bw ~ hw + age + cc + els + hc1 + hc2 + lei + length + weight + edu, d))
# cf _||_ hw | age, cc, els, hc_1, hc_2, lei, lngt, wght, edu
summary(lm(cf/bw ~ hw + age + els + hc1 + hc2 + lei + length + weight + edu + cc, d))
# els _||_ hw
summary(lm(els ~ hw, d))
# ga _||_ hw | age, cc, els, hc_1, hc_2, lei, lngt, wght, edu
summary(lm(ga ~ hw + age + cc + els + hc1 + hc2 + lei + length + weight + edu, d))
# hw _||_ lei | age
summary(glm(hw ~ lei + age, d, family = binomial()))
# hw _||_ lngt
summary(glm(hw ~ length, d, family = binomial()))
# hw _||_ wght + age
summary(glm(hw ~ weight + age, d, family = binomial()))
# lei _||_ lngt | els
summary(lm(lei ~ length + els, d))
# lei _||_ wght | els
summary(lm(lei ~ weight + els, d))
# added arrow from weight to lei (weight was measured before)
# lngt _||_ wght | els
summary(lm(length ~ weight + els, d))
# obviously added length -> weight
# lei _||_ lngt | cc, els, wght, edu
summary(lm(lei ~ length + cc + els + weight + edu, d))
# the final DAG is consistent with the data and my state of knowledge.
# I go ahead and check if we can estimate
# the paths coefficients of interest

# we are interested only in specific paths. I will check whether we can
# determine those path coefficients unbiased whereas we use hc as proxy for c
# els -> hc1
adjustmentSets(g, exposure = "els", outcome = "hc1", effect = "total")
adjustmentSets(g, exposure = "els", outcome = "hc1", effect = "direct")
# els -> hc2
adjustmentSets(g, exposure = "els", outcome = "hc2", effect = "direct")
# els -> lei
adjustmentSets(g, exposure = "els", outcome = "lei", effect = "direct")
# els -> ga
adjustmentSets(g, exposure = "els", outcome = "ga", effect = "direct")
# els -> bw
adjustmentSets(g, exposure = "els", outcome = "bw", effect = "direct")
# els -> cf
adjustmentSets(g, exposure = "els", outcome = "cf", effect = "direct")

# lei -> hc1
adjustmentSets(g, exposure = "lei", outcome = "hc1", effect = "direct")
# lei -> hc2
adjustmentSets(g, exposure = "lei", outcome = "hc2", effect = "direct")
# lei -> ga
adjustmentSets(g, exposure = "lei", outcome = "ga", effect = "direct")
# lei -> bw
adjustmentSets(g, exposure = "lei", outcome = "bw", effect = "direct")
# lei -> cf
adjustmentSets(g, exposure = "lei", outcome = "cf", effect = "direct")


# hc1 -> ga
adjustmentSets(g, exposure = "hc1", outcome = "ga", effect = "direct")
# hc1 -> cf
adjustmentSets(g, exposure = "hc1", outcome = "cf", effect = "direct")
# hc1 -> bw
adjustmentSets(g, exposure = "hc1", outcome = "bw", effect = "direct")

# hc2 -> ga
adjustmentSets(g, exposure = "hc2", outcome = "ga", effect = "direct")
# hc2 -> cf
adjustmentSets(g, exposure = "hc2", outcome = "cf", effect = "direct")
# hc2 -> bw
adjustmentSets(g, exposure = "hc2", outcome = "bw", effect = "direct")
# looks good, we can identify all paths without bias under the assumed DAG.
# Note however, that we use a proxy for cortisol. Thus, there is bias from wash out.
# the less measurement error from the hair cortisol measurement procedure, the 
# less biased estimates will

# for simplicity I will keep standard priors for now unless models do not
# converge.

form1 <- bf(lei ~ els + edu + weight)
get_prior(form1, complete(df_imp, 3), family = student())
if (!file.exists(here::here("data/coeftbl_cc_log_sc.Rds"))) {
  # els -> hc1
  adjustmentSets(g, exposure = "els", outcome = "hc1", effect = "total")
  els.hc1.t <- brm_multiple(
    family = student(),
    formula = hc1 ~ els,
    data = df_imp,
    file = "data/m_els.hc1.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    ),
  )
  # extract posterior
  b.els.hc1.t <- posterior_summary(els.hc1.t, variable = "b_els")
  adjustmentSets(g, exposure = "els", outcome = "hc1", effect = "direct")
  els.hc1 <- brm_multiple(
    family = student(),
    formula = hc1 ~ els + age + cc + edu + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_els.hc1",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    ),
  )
  # extract posterior
  b.els.hc1 <- posterior_summary(els.hc1, variable = "b_els")
  # els -> hc2
  adjustmentSets(g, exposure = "els", outcome = "hc2", effect = "direct")
  els.hc2 <- brm_multiple(
    family = student(),
    formula = hc2 ~ els + age + cc  + edu + hc1 + hw + lei + length + weight,
    data = df_imp,
    file = "data/m_els.hc2",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.hc2 <- posterior_summary(els.hc2, variable = "b_els")
  # els -> lei
  adjustmentSets(g, exposure = "els", outcome = "lei", effect = "direct")
  # just to doublecheck that nothing is open
  sum(paths(g, "els", "lei")$open)
  els.lei <- brm_multiple(
    family = student(),
    formula = lei ~ els + edu + weight,
    data = df_imp,
    file = "data/m_els.lei",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.lei <- posterior_summary(els.lei, variable = "b_els")
  # els -> ga
  adjustmentSets(g, exposure = "els", outcome = "ga", effect = "total")
  els.ga.t <- brm_multiple(
    family = student(),
    formula = ga ~ els,
    data = df_imp,
    file = "data/m_els.ga.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.ga.t <- posterior_summary(els.ga.t, variable = "b_els")
  adjustmentSets(g, exposure = "els", outcome = "ga", effect = "direct")
  els.ga <- brm_multiple(
    family = student(),
    formula = ga ~ els + age + cc + edu + hc1 + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_els.ga",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.ga <- posterior_summary(els.ga, variable = "b_els")
  # els -> bw
  adjustmentSets(g, exposure = "els", outcome = "bw", effect = "total")
  els.bw.t <- brm_multiple(
    family = student(),
    formula = bw ~ els,
    data = df_imp,
    file = "data/m_els.bw.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.bw.t <- posterior_summary(els.bw.t, variable = "b_els")
  adjustmentSets(g, exposure = "els", outcome = "bw", effect = "direct")
  els.bw <- brm_multiple(
    family = student(),
    formula = bw ~ els + age + cc + edu + ga + hc1 + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_els.bw",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.bw <- posterior_summary(els.bw, variable = "b_els")
  b.els.bw
  # els -> cf
  adjustmentSets(g, exposure = "els", outcome = "cf", effect = "total")
  els.cf.t <- brm_multiple(
    family = student(),
    formula = cf/bw ~ els,
    data = df_imp,
    file = "data/m_els.cf.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.cf.t <- posterior_summary(els.cf.t, variable = "b_els")
  adjustmentSets(g, exposure = "els", outcome = "cf", effect = "direct")
  els.cf <- brm_multiple(
    family = student(),
    formula = cf/bw ~ els + age + cc + edu + ga + hc1 + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_els.cf",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.els.cf <- posterior_summary(els.cf, variable = "b_els")

  # lei -> hc1
  adjustmentSets(g, exposure = "lei", outcome = "hc1", effect = "direct")
  lei.hc1 <- brm_multiple(
    family = student(),
    formula = hc1 ~ lei + age + cc + edu + els + weight + hw,
    data = df_imp,
    file = "data/m_lei.hc1",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.hc1 <- posterior_summary(lei.hc1, variable = "b_lei")

  # lei -> hc2 total
  adjustmentSets(g, exposure = "lei", outcome = "hc2", effect = "total")
  paths(g, "lei", "hc2")$paths
  sum(paths(g, "lei", "hc2")$open)
  lei.hc2.t <- brm_multiple(
    family = student(),
    formula = hc2 ~ lei,
    data = df_imp,
    file = "data/m_lei.hc2.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.hc2.t <- posterior_summary(lei.hc2.t, variable = "b_lei")

  # lei -> hc2
  adjustmentSets(g, exposure = "lei", outcome = "hc2", effect = "direct")
  lei.hc2 <- brm_multiple(
    family = student(),
    formula = hc2 ~ lei + age + cc + edu + els + hc1 + hw + length + weight,
    data = df_imp,
    file = "data/m_els.hc2",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.hc2 <- posterior_summary(lei.hc2, variable = "b_lei")
  # lei -> ga
  adjustmentSets(g, exposure = "lei", outcome = "ga", effect = "total")
  sum(paths(g, "lei", "ga")$open)
  lei.ga.t <- brm_multiple(
    family = student(),
    formula = ga ~ lei,
    data = df_imp,
    file = "data/m_lei.ga.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.ga.t <- posterior_summary(lei.ga.t, variable = "b_lei")
  adjustmentSets(g, exposure = "lei", outcome = "ga", effect = "direct")
  lei.ga <- brm_multiple(
    family = student(),
    formula = ga ~ lei + age + cc + edu + els + hc1 + hc2 + hw + length + weight,
    data = df_imp,
    file = "data/m_lei.ga",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.ga <- posterior_summary(lei.ga, variable = "b_lei")
  # lei -> bw
  adjustmentSets(g, exposure = "lei", outcome = "bw", effect = "total")
  # this path cannot be estimated in an unbiased manner
  sum(paths(g, "lei", "bw")$open)
  lei.bw.t <- brm_multiple(
    family = student(),
    formula = bw ~ lei,
    data = df_imp,
    file = "data/m_lei.bw.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.bw.t <-   posterior_summary(lei.bw.t, variable = "b_lei")
  adjustmentSets(g, exposure = "lei", outcome = "bw", effect = "direct")
  lei.bw <- brm_multiple(
    family = student(),
    formula = bw ~ lei + age + cc + edu + els + ga + hc1 + hc2 + hw + length + weight,
    data = df_imp,
    file = "data/m_lei.bw",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.bw <-   posterior_summary(lei.bw, variable = "b_lei")
  # lei -> cf
  adjustmentSets(g, exposure = "lei", outcome = "cf", effect = "total")
  sum(paths(g, "lei", "cf")$open)
  lei.cf.t <- brm_multiple(
    family = student(),
    formula = cf/bw ~ lei,
    data = df_imp,
    file = "data/m_lei.cf.t",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.cf.t <- posterior_summary(lei.cf.t, variable = "b_lei")
  adjustmentSets(g, exposure = "lei", outcome = "cf", effect = "direct")
  lei.cf <- brm_multiple(
    family = student(),
    formula = cf/bw ~ lei + age + cc + edu + els + ga + hc1 + hc2 + hw + length + weight,
    data = df_imp,
    file = "data/m_lei.cf",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.lei.cf <- posterior_summary(lei.cf, variable = "b_lei")



  # hc1 -> ga
  adjustmentSets(g, exposure = "hc1", outcome = "ga", effect = "direct")
  hc1.ga <- brm_multiple(
    family = student(),
    formula = ga ~ hc1 + age + cc + edu + els + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_hc1.ga",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.hc1.ga <- posterior_summary(hc1.ga, variable = "b_hc1")
  # hc1 -> cf
  adjustmentSets(g, exposure = "hc1", outcome = "cf", effect = "direct")
  hc1.cf <- brm_multiple(
    family = student(),
    formula = cf/bw ~ hc1 + age + cc + edu + els + ga + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_hc1.cf",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.hc1.cf <- posterior_summary(hc1.cf, variable = "b_hc1")
  # hc1 -> bw
  adjustmentSets(g, exposure = "hc1", outcome = "bw", effect = "direct")
  hc1.bw <- brm_multiple(
    family = student(),
    formula = bw ~ hc1 + age + cc + edu + els + ga + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_hc1.bw",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.hc1.bw <- posterior_summary(hc1.bw, variable = "b_hc1")

  # hc2 -> ga
  adjustmentSets(g, exposure = "hc2", outcome = "ga", effect = "direct")
  hc2.ga <- brm_multiple(
    family = student(),
    formula = ga ~ hc2 + age + cc + edu + els + hc1 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_hc2.ga",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.hc2.ga <- posterior_summary(hc2.ga, variable = "b_hc2")
  # hc2 -> cf
  adjustmentSets(g, exposure = "hc2", outcome = "cf", effect = "direct")
  hc2.cf <- brm_multiple(
    family = student(),
    formula = cf/bw ~ hc2 + age + cc + edu + els + ga + hc1 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_hc2.cf",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.hc2.cf <- posterior_summary(hc2.cf, variable = "b_hc2")
  # hc2 -> bw
  adjustmentSets(g, exposure = "hc2", outcome = "bw", effect = "direct")
  hc2.bw <- brm_multiple(
    family = student(),
    formula = bw ~ hc2 + age + cc + edu + els + ga + hc1 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_hc2.bw",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.hc2.bw <- posterior_summary(hc2.bw, variable = "b_hc2")

  # hc1 -> hc2
  adjustmentSets(g, exposure = "hc1", outcome = "hc2", effect = "direct")
  hc1.hc2 <- brm_multiple(
    family = student(),
    formula = hc2 ~ hc1 + age + cc + edu + els + hw + lei + length + weight,
    data = df_imp,
    file = "data/m_hc1.hc2",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.hc1.hc2 <- posterior_summary(hc1.hc2, variable = "b_hc1")

  # els -> edu
  # sum(paths(g, from = "els", to = "edu")$open)
  adjustmentSets(g, exposure = "els", outcome = "edu", effect = "direct")
  sum(paths(g, "els", "edu")$open)
  # normal model with edu as outcome does not converge. edu is an ordinal
  # variable and modeling it as continuous as predictor works usually well
  # but in this particular case I need to fit an ordinal model
  str(complete(df_imp))
  els.edu <- brm_multiple(data = df_imp_edu,
  file = "data/m_els.edu",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b")
    ),
      family = cumulative,
      edu ~ 1 + els
    )


  b.els.edu <- posterior_summary(els.edu, variable = "b_els")
  # edu -> lei
  adjustmentSets(g, exposure = "edu", outcome = "lei", effect = "direct")
  sum(paths(g, "edu", "lei")$open)
  edu.lei <- brm_multiple(
    family = student(),
    formula = lei ~ edu + age + els + weight,
    data = df_imp,
    file = "data/m_edu.lei",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.edu.lei <- posterior_summary(edu.lei, variable = "b_edu")


  # ga -> bw
  adjustmentSets(g, exposure = "ga", outcome = "bw", effect = "direct")
  ga.bw <- brm_multiple(
    family = student(),
    formula = bw ~ ga + age + cc + edu + els + hc1 + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_ga.bw",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.ga.bw <- posterior_summary(ga.bw, variable = "b_ga")

  # ga -> cf
  adjustmentSets(g, exposure = "ga", outcome = "cf", effect = "direct")
  ga.cf <- brm_multiple(
    family = student(),
    formula = cf/bw ~ ga + age + cc + edu + els + hc1 + hc2 + lei + length + weight + hw,
    data = df_imp,
    file = "data/m_ga.cf",
    prior = c(
      prior(normal(0, 10), "Intercept"),
      prior(normal(0, 0.5), "b"),
      prior(exponential(1), "sigma"),
      prior(gamma(2, 0.1), "nu")
    )
  )
  b.ga.cf <- posterior_summary(ga.cf, variable = "b_ga")

  coeftbl <- tibble(
    exposure = c(rep("els", 10), rep("lei", 8), rep("hc1", 3), rep("hc2", 3), "hc1", "els", "edu", "ga", "ga"),
    outcome = c(
      "hc1", "hc1", "hc2", "lei", "ga", "ga", "bw", "bw", "cf", "cf", "hc1",
      "hc2", "ga", "ga", "bw", "bw", "cf", "cf", "ga", "bw", "cf", "ga", "bw", "cf", "hc2", "edu", "lei", "bw", "cf"
    ),
    effect = c(
      "total", rep("direct", 4), "total", "direct", "total", "direct",
      "total", "direct", "direct", "total", "direct", "total", "direct", "total", rep("direct", 12)
    ),
    b = c(
      b.els.hc1.t[1, 1],
      b.els.hc1[1, 1],
      b.els.hc2[1, 1],
      b.els.lei[1, 1],
      b.els.ga[1, 1],
      b.els.ga.t[1, 1],
      b.els.bw[1, 1],
      b.els.bw.t[1, 1],
      b.els.cf[1, 1],
      b.els.cf.t[1, 1],
      b.lei.hc1[1, 1],
      b.lei.hc2[1, 1],
      b.lei.ga.t[1, 1],
      b.lei.ga[1, 1],
      b.lei.bw.t[1, 1],
      b.lei.bw[1, 1],
      b.lei.cf.t[1, 1],
      b.lei.cf[1, 1],
      b.hc1.ga[1, 1],
      b.hc1.bw[1, 1],
      b.hc1.cf[1, 1],
      b.hc2.ga[1, 1],
      b.hc2.bw[1, 1],
      b.hc2.cf[1, 1],
      b.hc1.hc2[1, 1],
      b.els.edu[1, 1],
      b.edu.lei[1, 1],
      b.ga.bw[1, 1],
      b.ga.cf[1, 1]
    ),
    lower = c(
      b.els.hc1.t[1, 3],
      b.els.hc1[1, 3],
      b.els.hc2[1, 3],
      b.els.lei[1, 3],
      b.els.ga[1, 3],
      b.els.ga.t[1, 3],
      b.els.bw[1, 3],
      b.els.bw.t[1, 3],
      b.els.cf[1, 3],
      b.els.cf.t[1, 3],
      b.lei.hc1[1, 3],
      b.lei.hc2[1, 3],
      b.lei.ga.t[1, 3],
      b.lei.ga[1, 3],
      b.lei.bw.t[1, 3],
      b.lei.bw[1, 3],
      b.lei.cf.t[1, 3],
      b.lei.cf[1, 3],
      b.hc1.ga[1, 3],
      b.hc1.bw[1, 3],
      b.hc1.cf[1, 3],
      b.hc2.ga[1, 3],
      b.hc2.bw[1, 3],
      b.hc2.cf[1, 3],
      b.hc1.hc2[1, 3],
      b.els.edu[1, 3],
      b.edu.lei[1, 3],
      b.ga.bw[1, 3],
      b.ga.cf[1, 3]
    ),
    upper = c(
      b.els.hc1.t[1, 4],
      b.els.hc1[1, 4],
      b.els.hc2[1, 4],
      b.els.lei[1, 4],
      b.els.ga[1, 4],
      b.els.ga.t[1, 4],
      b.els.bw[1, 4],
      b.els.bw.t[1, 4],
      b.els.cf[1, 4],
      b.els.cf.t[1, 4],
      b.lei.hc1[1, 4],
      b.lei.hc2[1, 4],
      b.lei.ga.t[1, 4],
      b.lei.ga[1, 4],
      b.lei.bw.t[1, 4],
      b.lei.bw[1, 4],
      b.lei.cf.t[1, 4],
      b.lei.cf[1, 4],
      b.hc1.ga[1, 4],
      b.hc1.bw[1, 4],
      b.hc1.cf[1, 4],
      b.hc2.ga[1, 4],
      b.hc2.bw[1, 4],
      b.hc2.cf[1, 4],
      b.hc1.hc2[1, 4],
      b.els.edu[1, 4],
      b.edu.lei[1, 4],
      b.ga.bw[1, 4],
      b.ga.cf[1, 4]
    ))

    # add probablities of directions
    effects <- c(
      "b_els",
      "b_els",
      "b_els",
      "b_els",
      "b_els",
      "b_els",
      "b_els",
      "b_els",
      "b_els",
      "b_els",
      "b_lei",
      "b_lei",
      "b_lei",
      "b_lei",
      "b_lei",
      "b_lei",
      "b_lei",
      "b_lei",
      "b_hc1",
      "b_hc1",
      "b_hc1",
      "b_hc2",
      "b_hc2",
      "b_hc2",
      "b_hc1",
      "b_els",
      "b_edu",
      "b_ga",
      "b_ga"
    )
    models <- list(
      els.hc1.t,
      els.hc1,
      els.hc2,
      els.lei,
      els.ga,
      els.ga.t,
      els.bw,
      els.bw.t,
      els.cf,
      els.cf.t,
      lei.hc1,
      lei.hc2,
      lei.ga.t,
      lei.ga,
      lei.bw.t,
      lei.bw,
      lei.cf.t,
      lei.cf,
      hc1.ga,
      hc1.bw,
      hc1.cf,
      hc2.ga,
      hc2.bw,
      hc2.cf,
      hc1.hc2,
      els.edu,
      edu.lei,
      ga.bw,
      ga.cf
    )

    ps <- map2(models, effects, function(model, effect) {
      post <- posterior_samples(model, variable = effect)
      p <- round(mean(post[[effect]] > 0), 3)
      # ifelse(p > 0.5, p, 1 - p)
    }) %>% as.numeric()
    coeftbl$p <- ps
  save(coeftbl, file = here::here("data/coeftbl_cc_log_sc.Rds"))
  save(models, file = here::here("data/models_cc_log_sc.Rds"))
  } else {
  load(here::here("data/coeftbl_cc_log_sc.Rds"))
  load(here::here("data/models_cc_log_sc.Rds"))
}


coeftbl <- coeftbl %>%  mutate(across(where(is.numeric), function(x) round(x, 2))) %>%
  arrange(exposure, outcome)
coeftbl
count(d, cc)
colnames(coeftbl) <- str_to_title(colnames(coeftbl))
coeftbl %>% filter(!(Exposure == "lei" & Outcome == "bw" & Effect == "total"))

# interaction ELS * LEI for hc1
adjustmentSets(g, exposure = "lei", outcome = "hc1", effect = "direct")
lei.hc1.x <- brm_multiple(
  family = student(),
  formula = hc1 ~ lei * els + age + cc + edu + weight + hw,
  data = df_imp,
  file = here::here("data/lei.hc1.x"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
models[[length(models) + 1]] <- lei.hc1.x
# in that case we would expect that the effect of lei on hc increases among
# the mothers that had no ELS.
h2 <- c(h2 = "lei:els > 0")
hypothesis(lei.hc1.x, h2)
mean(posterior_samples(lei.hc1.x, pars = "b_lei:els")[, 1] > 0)
quantile(complete(df_imp)$els, 0.75)
summary(lei.hc1.x)

# to explore it extensively lets make a category upper quantile els
dfimp_alt <- map(1:50, function(m) {
  complete(df_imp, m = m)[[1]] %>%
  mutate(els_cat = as.factor(ifelse(els < 0.653, 0, 1)))
})


lei.hc.1.x.cat <- brm_multiple(
  family = student(),
  formula = hc1 ~ lei * els_cat + age + cc + edu + weight + hw,
  data = dfimp_alt,
  file = here::here("data/lei.hc.1.x.cat"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(lei.hc.1.x.cat)
models[[length(models) + 1]] <- lei.hc.1.x.cat

# interaction ELS * LEI for hc2
adjustmentSets(g, exposure = "lei", outcome = "hc2", effect = "direct")
lei.hc2.x <- brm_multiple(
  family = student(),
  formula = hc2 ~ lei * els + hc1 + age + cc + edu + weight + hw,
  data = df_imp,
  file = here::here("data/lei.hc2.x"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
# in that case we would expect that the effect of lei on hc increases among
# the mothers that had no ELS.
mean(posterior_samples(lei.hc2.x, pars = "b_lei:els")[, 1] > 0)
summary(lei.hc2.x)
models[[length(models) + 1]] <- lei.hc2.x

# we could derive the following outcomes from the existing models but this 
# is easier to interpret for most:
# birthoutcomes ~ change over time in HC.
dfimp_alt <- map(1:50, function(m) {
  complete(df_imp, m = m)[[1]] %>%
  mutate(hcdiff = hc2 - hc1, hcdiff = scale(hcdiff)[, 1])
})

# hcdiff -> bw
hcdiff.bw <- brm_multiple(
  family = student(),
  formula = bw ~ hcdiff + age + cc + edu + els + ga + lei + length + weight + hw,
  data = dfimp_alt,
  file = here::here("data/hcdiff.bw"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)

b.hcdiff.bw <- posterior_summary(hcdiff.bw, variable = "b_hcdiff")
b.hcdiff.bw.post <- as_draws_df(hcdiff.bw, variable = "b_hcdiff") %>%
  summarise(mean(b_hcdiff < 0))
b.hcdiff.bw.post
b.hcdiff.bw
summary(hcdiff.bw)

models[[length(models) + 1]] <- hcdiff.bw


colnames(complete(df_imp))
adjustmentSets(g, exposure = "hc2", outcome = "cf", effect = "direct")
# hcdiff -> ci
hcdiff.cf <- brm_multiple(
  family = student(),
  formula = cf ~ hcdiff + age + cc + edu + els + ga + lei + length + weight + hw,
  data = dfimp_alt,
  file = here::here("data/hcdiff.cf"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)

summary(hcdiff.cf)
models[[length(models) + 1]] <- hcdiff.cf


# hcdiff -> ga
hcdiff.ga <- brm_multiple(
  family = student(),
  formula = ga ~ hcdiff + age + cc + edu + els + lei + length + weight + hw,
  data = dfimp_alt,
  file = here::here("data/hcdiff.ga"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)

summary(hcdiff.ga)
models[[length(models) + 1]] <- hcdiff.ga


# conclusion: stronger increase in HC during pregnancy seems to be associated 
# with lower birthweight AFTER accoutning for gestational age and other
# factors. We could have also derived that from this model:
hc1.bw <- brm_multiple(
  family = student(),
  formula = bw ~ hc1 + age + cc + edu + els + ga + hc2 + lei + length + weight + hw,
  data = df_imp,
  file = here::here("data/hc1.bw"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(hc1.bw)

# Multi level model to test whether increase over time is 
# statistically significant 
dfimp_time <- map(1:50, function(m) {
  complete(df_imp, m = m)[[1]] %>%
  mutate(id = 1:69) %>%
  pivot_longer(contains("hc"), names_to = "time", values_to = "hc") %>%
  mutate(time = as.factor(as.numeric(str_extract(time, "\\d"))))
})

dftime <- complete(df_imp) %>% 
  mutate(id = 1:69) %>%
  pivot_longer(contains("hc"), names_to = "time", values_to = "hc") %>%
  mutate(time = as.factor(as.numeric(str_extract(time, "\\d"))))


formtime = bf(hc ~ time + age + cc + edu + els + lei + length + weight + hw + (1 + time | id))

get_prior(formtime, data = dftime)
time.hc <- brm_multiple(
  family = student(),
  formula = hc ~ time + age + cc + edu + els + lei + length + weight + hw + (1 + time | id),
  data = dfimp_time,
  file = here::here("data/time.hc"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu"),
    prior(exponential(1), "sd"),
    prior(lkj(1), "cor")
  ),
  # warmup = 1000,
  # iter = 2000,
  #control = list(adapt_delta = 0.9, max_treedepth = 10)
)
summary(time.hc)

pred <- predict(time.hc, re_formula = NULL) %>% as.data.frame()
resid_df <- bind_cols(pred, dfimp_time[[1]]) %>%
  mutate(resid = hc - Estimate) %>%
  select(resid, time)

# residuals a homogeneous when accouting for individual variability
ggplot(resid_df, aes(time, resid)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.1)

models[[length(models) + 1]] <- time.hc

# time.hc2 <- brm_multiple(
#   family = student(),
#   formula = bf(hc ~ time + age + cc + edu + els + lei + length + weight + hw + (1 + time | id), sigma ~ time), # some factors maybe affect C depending on when they were present most OR affect distal hair part more
#   data = dfimp_time,
#   file = here::here("data/time.hc2"),
#   prior = c(
#     prior(normal(0, 10), "Intercept"),
#     prior(normal(0, 0.5), "b"),
#     #prior(exponential(1), "sigma"),
#     prior(gamma(2, 0.1), "nu"),
#     prior(exponential(1), "sd"),
#     prior(lkj(1), "cor")
#   ),
#   # warmup = 1000,
#   # iter = 2000,
#   #control = list(adapt_delta = 0.9, max_treedepth = 10)
# )



# time.hc.2 <- brm(
#   family = student(),
#   formula = hc ~ time + cc + hw + (1 + time | id),
#   data = dftime,
#   prior = c(
#     prior(normal(0, 10), "Intercept"),
#     prior(normal(0, 0.5), "b"),
#     prior(exponential(1), "sigma"),
#     prior(gamma(2, 0.1), "nu"),
#     prior(exponential(1), "sd")
#   )
# )
# summary(time.hc.2)

# are there variables that explain the change in HC over time?
x.hcdiff <- brm_multiple(
  family = student(),
  formula = hcdiff ~  age + cc + edu + els * lei + length + weight + hw,
  data = dfimp_alt,
  file = here::here("data/x.hcdiff"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(x.hcdiff)

models[[length(models) + 1]] <- x.hcdiff

x.hcdiff2 <- brm_multiple(
  family = student(),
  formula = hcdiff ~  age + cc + edu + els + lei + length + weight + hw,
  data = dfimp_alt,
  file = here::here("data/x.hcdiff2"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(x.hcdiff2)

models[[length(models) + 1]] <- x.hcdiff.2

temp <- as_draws_df(time.hc)
slopes <- select(temp, matches("r_id\\[\\d+,time2\\]")) %>% 
  pivot_longer(everything(), names_to = "id", values_to = "slopes") %>%
  group_by(id) %>%
  summarise(slope = mean(slopes), sdslope = sd(slopes))  %>%
  mutate(id = as.integer(str_match(id, "\\[(\\d+),time2\\]")[, 2])) %>%
  arrange(id)
dfimp_slopes <- map(dfimp_time, function(d) {
  d %>%
    pivot_wider(
      names_from = time, 
      values_from = hc, 
      names_glue = "{.value}{time}") %>%
    mutate(hcdiff = hc2 - hc1, scale(hcdiff)[, 1]) %>%
    left_join(slopes, by = "id")}) 

x.hcdiff3 <- brm_multiple(
  family = student(),
  formula = slope ~  age + edu + els + lei + length + weight,
  data = dfimp_slopes,
  file = here::here("data/x.hcdiff3"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(x.hcdiff3)


x.hcdiff4 <- brm_multiple(
  family = student(),
  formula = bf(slope | se(sdslope, sigma = TRUE) ~  age + edu + els + lei + length + weight),
  data = dfimp_slopes,
  file = here::here("data/x.hcdiff4"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(x.hcdiff4)


x.hcdiff5 <- brm_multiple(
  family = student(),
  formula = slope ~  age + edu + els + lei + length + weight + hw + cc,
  data = dfimp_slopes,
  file = here::here("data/x.hcdiff5"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(x.hcdiff5)

# hcdiff -> bw
hcdiff.bw2 <- brm_multiple(
  family = student(),
  formula = bw ~ me(slope, sdslope) + age + cc + edu + els + ga + lei + length + weight + hw,
  data = dfimp_slopes,
  file = here::here("data/hcdiff.bw2"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(hcdiff.bw2)


# hcdiff -> bw
hcdiff.bw3 <- brm_multiple(
  family = student(),
  formula = bw ~ slope + age + cc + edu + els + ga + lei + length + weight + hw,
  data = dfimp_slopes,
  file = here::here("hcdiff.bw3"),
  prior = c(
    prior(normal(0, 10), "Intercept"),
    prior(normal(0, 0.5), "b"),
    prior(exponential(1), "sigma"),
    prior(gamma(2, 0.1), "nu")
  )
)
summary(hcdiff.bw3)

colnames(temp)
# plot random effects 
re <- select(
  temp, 
  matches("r_id\\[\\d+,time2\\]"), 
  matches("r_id\\[\\d+,Intercept\\]")) %>% 
  pivot_longer(everything(), names_to = "effect", values_to = "values") %>%
  group_by(effect) %>%
  summarise(estimate = median(values)) %>%
  ungroup %>%
  mutate(
    id = as.integer(str_match(effect, "\\[(\\d+),.+\\]")[, 2]),
    effect = ifelse(str_detect(effect, "Intercept"), "Intercepts", "Slopes")
  ) %>%
  pivot_wider(names_from = effect, values_from = estimate)
p_re <- ggplot(re, aes(Intercepts, Slopes)) +
  geom_point(size = 5) +
  theme_bw(base_size = 20)
save(p_re, file = here::here("data/p_re.Rds"))  


mutate(re, Slopes = Slopes + 1.16) %>%
  filter(Slopes <=0)
# how to present the other effect?!
summary(x.hcdiff3)


# I will make the supplementary tables in an extra script
save(models, file = "data/models_for_supp.Rds")
