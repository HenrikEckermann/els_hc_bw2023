library(tidyverse)
library(brms)
library(glue)

# Bayesian helper functions for after model fitting 
get_lower <- function(x, level = 0.025) {quantile(x, level)}
get_upper <- function(x, level = 0.975) {quantile(x, level)}

print_estimate <- function(model, effect){
  tab <- posterior_samples(model, effect) %>%
    select(all_of(effect)) %>%
    summarise_all(
      list(beta = median, lower = get_lower, upper = get_upper)
    ) %>% mutate(across(everything(), function(x) format(round(x, 3), nsmall = 3)))
  return(glue("$\beta$ = {tab$beta}, 95% CI [{tab$lower}, {tab$upper}]"))
}

# to create tables
brms_tbl <- function(model, exclude_id = TRUE) {
  tb <- posterior_samples(model) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    group_by(Parameter) %>%
    summarise(
      Estimate = median(Value),
      SD = sd(Value),
      lower = quantile(Value, 0.025),
      upper = quantile(Value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), function(x) format(round(x, 2), nsmall = 2))) %>%
    mutate("95% CI" = glue("{lower} : {upper}")) %>%
    select(Parameter, Estimate, SD, "95% CI")

  if (exclude_id) {
    tb <- tb  %>%
      filter(!grepl("^r_id", Parameter), Parameter != "lp__")
  }
  return(tb)
}

load(file = "data/models_for_supp.Rds")


# to avoid redundant tables I retain ids only of models that we do not yet have
# in f_ids. First needs to standardize the forms in order to compare them:
forms <- map_chr(models, ~as.character(.x$formula)[1])
forms <- map(forms, function(form) {
  f <- str_replace(form, "~", "+")
  set <- sort(str_split_1(f, " \\+ "))
  con <- str_flatten(set, "_")
  return(con)
})

# now only retain unique forms
f_ids <- c(1)
for (i in seq_along(forms)) {
  if (!all(forms[[i]] %in% forms[f_ids])) {
    f_ids <- c(f_ids, i)
  }
}
f_ids

table_nr <- seq_along(models[f_ids])

# print header followed by table for each model 
file.create(here::here("Rmd/table_supplement.Rmd"))
write_lines(
  '---\noutput:\n  pdf_document:\n    toc: false\n    toc_depth : 2\n---', file = here::here("Rmd/table_supplement.Rmd"), append = FALSE)

pmap(list(models[f_ids], table_nr), function(model, nr) {
  # write_lines(
  #   x = glue("\n\n### Supplementary Table {nr}\n \n"),
  #   file = here::here("article/table_supplement.Rmd"), append = TRUE)
  
  write_lines(
    x = knitr::kable(
      brms_tbl(model, exclude_id = TRUE), # ID must be excluded to guarantee anonymity
      caption = glue('Coefficients for the model using formula: {as.character(model$formula)[1]}\n \n')
      ),
    file = here::here("Rmd/table_supplement.Rmd"), append = TRUE)
})


# before rendering, please manually escape "|" in the multilevel model form
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
rmarkdown::render(here::here("Rmd/table_supplement.Rmd"))
