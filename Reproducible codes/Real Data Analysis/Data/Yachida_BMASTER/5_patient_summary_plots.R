library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
setwd("U:/BMASTER/Data raw/Yachida_BMASTER")
# ---- Load & coerce ----
df <- read_csv("meta_used.csv", show_col_types = FALSE) |>
  mutate(
    Age            = suppressWarnings(as.numeric(Age)),
    BMI            = suppressWarnings(as.numeric(BMI)),
    `Brinkman.Index` = suppressWarnings(as.numeric(`Brinkman.Index`)),
    Alcohol        = suppressWarnings(as.numeric(Alcohol))
  )

# ---- Categories (same thresholds as before) ----
smoking_cat <- function(x){
  ifelse(is.na(x), NA_character_,
         ifelse(x <= 0, "Never",
                ifelse(x < 400, "Light", "Heavy")))
}
drinking_cat <- function(x){
  ifelse(is.na(x), NA_character_,
         ifelse(x <= 0, "Abstainer",
                ifelse(x <= 100, "Low",
                       ifelse(x <= 250, "Moderate", "High"))))
}

df <- df |>
  mutate(
    Smoking.Category  = factor(smoking_cat(`Brinkman.Index`),
                               levels = c("Never","Light","Heavy")),
    Drinking.Category = factor(drinking_cat(Alcohol),
                               levels = c("Abstainer","Low","Moderate","High")),
    Gender            = str_to_title(as.character(Gender))
  )

# ---- Style helper (takes a ggplot object 'p') ----
finalize_plot <- function(p, title, xlab = NULL, ylab = NULL) {
  p +
    ggplot2::ggtitle(title) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )
}

save_plot <- function(p, filename, w = 6, h = 4, dpi = 300) {
  ggsave(filename, p, width = w, height = h, dpi = dpi)
  print(p)
}

# ---- 1) Age ----
age <- df |> filter(!is.na(Age)) |> pull(Age)
bins_age <- min(20, max(8, floor(sqrt(length(age)) + 1)))
p_age <- ggplot(df, aes(x = Age)) +
  geom_histogram(bins = bins_age, color = "black", fill = "orange")
p_age <- finalize_plot(p_age, "Age (years)", "Age", "Count")
save_plot(p_age, "Supp_fig_age_hist.png")

# ---- 2) BMI ----
bmi <- df |> filter(!is.na(BMI)) |> pull(BMI)
bins_bmi <- min(20, max(8, floor(sqrt(length(bmi)) + 1)))
p_bmi <- ggplot(df, aes(x = BMI)) +
  geom_histogram(bins = bins_bmi, color = "black", fill = "orange")
p_bmi <- finalize_plot(p_bmi, "Body Mass Index (kg/m²)", "BMI", "Count")
save_plot(p_bmi, "Supp_fig_bmi_hist.png")

# ---- 3) Smoking ----
smk_counts <- df |>
  count(Smoking.Category, .drop = FALSE) |>
  mutate(Smoking.Category = factor(Smoking.Category, levels = c("Never","Light","Heavy")))
p_smk <- ggplot(smk_counts, aes(x = Smoking.Category, y = n)) +
  geom_col(color = "black", fill = "orange") +
  geom_text(aes(label = n), vjust = -0.2)
p_smk <- finalize_plot(p_smk, "Smoking (Brinkman Index categories)", "Category", "N")
save_plot(p_smk, "Supp_fig_smoking_bar.png")

# ---- 4) Drinking ----
drk_counts <- df |>
  count(Drinking.Category, .drop = FALSE) |>
  mutate(Drinking.Category = factor(Drinking.Category,
                                    levels = c("Abstainer","Low","Moderate","High")))
title_drink <- "Drinking (Alcohol categories)\n(0=Abstainer; 0–100 Low; 100–250 Moderate; >250 High)"
p_drk <- ggplot(drk_counts, aes(x = Drinking.Category, y = n)) +
  geom_col(color = "black", fill = "orange") +
  geom_text(aes(label = n), vjust = -0.2)
p_drk <- finalize_plot(p_drk, title_drink, "Category", "N")
save_plot(p_drk, "Supp_fig_drinking_bar.png")

# ---- 5) Sex ----
sex_counts <- df |> mutate(Gender = factor(Gender)) |> count(Gender)
p_sex <- ggplot(sex_counts, aes(x = Gender, y = n)) +
  geom_col(color = "black", fill = "orange") +
  geom_text(aes(label = n), vjust = -0.2)
p_sex <- finalize_plot(p_sex, "Sex", "Sex", "N")
save_plot(p_sex, "Supp_fig_sex_bar.png")
