vaers_COV2020_2025 <- list(Merged_final20, Merged_final21, Merged_final22, Merged_final23, Merged_final24, Merged_final25) |>  bind_rows()
View(vaers_COV2020_2025)
write.csv( vaers_COV2020_2025, file = "VAERS_COVIDCASES_2020_2025.csv", row.names = FALSE)
#Let the WoRk BeGin
library(tidyverse)
library(lubridate)
library(fasttime) #---------fasttttt
# date problems─
parse_date_fast <- function(x) {
  # Try ISO format first, then MDY
  out <- as.Date(x, format = "%Y-%m-%d")
  bad <- is.na(out)
  if (any(bad)) {
    out[bad] <- as.Date(x[bad], format = "%m/%d/%Y")
  }
  out
}
report_year_vec <- as.integer(vaers_COV2020_2025$REPORT_YEAR)
vax_year_vec    <- year(parse_date_fast(vaers_COV2020_2025$VAX_DATE))

# ── Pre-compute grepl vectors (once each, vectorised) ─────────────────────
sex_upper  <- toupper(vaers_COV2020_2025$SEX)
is_pfizer  <- grepl("PFIZER|BIONTECH|BNT162",vaers_COV2020_2025$VAX_MANU, ignore.case = TRUE) |grepl("PFIZER|BNT162", vaers_COV2020_2025$VAX_TYPE, ignore.case = TRUE)
is_moderna <- grepl("MODERNA|mRNA-1273",vaers_COV2020_2025$VAX_MANU,ignore.case = TRUE) |grepl("MODERNA|mRNA-1273", vaers_COV2020_2025$VAX_TYPE, ignore.case = TRUE)
is_janssen <- grepl("JANSSEN|J&J|Ad26|COV2",vaers_COV2020_2025$VAX_MANU, ignore.case = TRUE) |grepl("JANSSEN|Ad26", vaers_COV2020_2025$VAX_TYPE, ignore.case = TRUE)
#Preparing the tableee
df_prepared <- vaers_COV2020_2025 |>
  mutate(
    report_year = report_year_vec,
    vax_year    = vax_year_vec,
    age_group = factor(
      case_when(
        is.na(AGE_YRS)                 ~ "Unknown",
        AGE_YRS >  0  & AGE_YRS <= 10  ~ "0-10",
        AGE_YRS >= 11 & AGE_YRS <= 20  ~ "11-20",
        AGE_YRS >= 21 & AGE_YRS <= 30  ~ "21-30",
        AGE_YRS >= 31 & AGE_YRS <= 40  ~ "31-40",
        AGE_YRS >= 41 & AGE_YRS <= 50  ~ "41-50",
        AGE_YRS >= 51 & AGE_YRS <= 60  ~ "51-60",
        AGE_YRS >= 61 & AGE_YRS <= 70  ~ "61-70",
        AGE_YRS >= 71 & AGE_YRS <= 80  ~ "71-80",
        AGE_YRS >= 81 & AGE_YRS <= 90  ~ "81-90",
        AGE_YRS >= 91 & AGE_YRS <= 100 ~ "91-100",
        TRUE                            ~ "Other"
      ),
      levels = c("0-10","11-20","21-30","31-40","41-50", "51-60","61-70","71-80","81-90","91-100","Other")),
    gender = factor(
      case_when(
        sex_upper %in% c("F","FEMALE") ~ "Female",
        sex_upper %in% c("M","MALE")   ~ "Male",
        TRUE                           ~ "Unknown"
      ),
      levels = c("Female","Male","Unknown")),
    vax_type = factor(
      case_when(
        is_pfizer  ~ "BNT162b2 (Pfizer)",
        is_moderna ~ "mRNA-1273 (Moderna)",
        is_janssen ~ "Ad26.COV2.S (Janssen)",
        TRUE       ~ "Other"      ),
      levels = c("BNT162b2 (Pfizer)","mRNA-1273 (Moderna)",
                 "Ad26.COV2.S (Janssen)","Other")),
    admin_dose = factor(
      case_when(
        is.na(VAX_DOSE_SERIES)   ~ "Unknown",
        VAX_DOSE_SERIES == 1     ~ "Dose 1",
        VAX_DOSE_SERIES == 2     ~ "Dose 2",
        VAX_DOSE_SERIES == 3     ~ "Dose 3",
        VAX_DOSE_SERIES == 4     ~ "Dose 4",
        VAX_DOSE_SERIES == 5     ~ "Dose 5",
        VAX_DOSE_SERIES >= 6     ~ "Dose 6+",
        TRUE                     ~ "Unknown"
      ),
      levels = c("Dose 1","Dose 2","Dose 3","Dose 4",
                 "Dose 5","Dose 6+")),
    onset_interval = factor(
      case_when(
        suppressWarnings(as.numeric(DAYS_TO_ONSET)) == 0   ~ "Day 0",
        suppressWarnings(as.numeric(DAYS_TO_ONSET)) <= 7   ~ "Days 1-7",
        suppressWarnings(as.numeric(DAYS_TO_ONSET)) <= 14  ~ "Days 8-14",
        suppressWarnings(as.numeric(DAYS_TO_ONSET)) <= 21  ~ "Days 15-21",
        suppressWarnings(as.numeric(DAYS_TO_ONSET)) <= 28  ~ "Days 22-28",
        TRUE                                               ~ "Days > 28"
      ),
      levels = c("Day 0","Days 1-7","Days 8-14", "Days 15-21","Days 22-28","Days > 28","Unknown"),
      ordered = TRUE),
    died           = DIED     == "Y",
    life_threat    = L_THREAT == "Y",
    hospitalized   = HOSPITAL == "Y",
    prolonged_hosp = X_STAY   == "Y",
    er_or_urgent   = ER_VISIT == "Y" | X_STAY == "Y",
    disabled       = DISABLE  == "Y",
    outcome_severity = factor(
      case_when(
        DIED     == "Y"                     ~ "Died",
        L_THREAT == "Y"                     ~ "Life-threatening",
        HOSPITAL == "Y"                     ~ "Hospitalized",
        ER_VISIT == "Y" | X_STAY == "Y"     ~ "ER / Urgent care",
        TRUE                                ~ "Other / Not serious"
      ),
      levels = c("Died","Life-threatening","Hospitalized", "ER / Urgent care","Other / Not serious")),
    has_is  = replace_na(has_is,  FALSE),
    has_hm  = replace_na(has_hm,  FALSE),
    stroke = as.integer(has_is | has_hm),
   stroke_type = factor(
      case_when(
        has_hm  ~ "Haemorrhagic",
        has_is  ~ "Ischaemic",
        TRUE    ~ "No Stroke"
      ),
      levels = c("Ischaemic","Haemorrhagic","No Stroke")
    )) |>
  select( VAERS_ID, report_year, vax_year, AGE_YRS, age_group,gender, vax_type,admin_dose, onset_interval, died, life_threat, hospitalized, er_or_urgent, prolonged_hosp, disabled, outcome_severity, has_is, has_hm, stroke, stroke_type)
glimpse(df_prepared)
cat("Columns now:", ncol(df_prepared), "\n")   # should be 21
cat("Rows       :", nrow(df_prepared), "\n")
view(df_prepared)
cat("Ischaemic only    :", sum( df_prepared$has_is & !df_prepared$has_hm), "\n")
cat("Ischaemic only    :", sum( df_prepared$has_is), "\n")
cat("Haemorrhagic only :", sum(!df_prepared$has_is &  df_prepared$has_hm), "\n")
cat("Overlap IS + HM   :", sum( df_prepared$has_is &  df_prepared$has_hm),  "\n")
cat("Any stroke        :", sum( df_prepared$stroke), "\n")
cat("Stroke rate (%):", round(mean(df_prepared$stroke) * 100, 4), "\n")

# Distribution check
df_prepared |> count(stroke_type)
df_prepared |> count(report_year)
df_prepared |> count(vax_type,   sort = TRUE)
df_prepared |> count(age_group,  sort = TRUE)
df_prepared |> count(onset_interval)
library(gtsummary)
df_prepared <- df_prepared |>
  mutate(
    died  = !is.na(died) & died == "Y",
    life_threat  = !is.na(life_threat) & life_threat == "Y",
    hospitalized = !is.na(hospitalized) & hospitalized == "Y",
    er_or_urgent = (!is.na(er_or_urgent) & er_or_urgent == "Y") |
      (!is.na(prolonged_hosp)   & prolonged_hosp   == "Y"),
    disabled     = !is.na(disabled)  & disabled  == "Y"
  )
library(gtsummary)
library(gt)

table1 <- df_prepared |>
  tbl_summary(
    by = stroke_type,
    include = c(AGE_YRS,age_group,gender,vax_type,admin_dose,onset_interval,outcome_severity),
    label = list(AGE_YRS ~ "Mean Age (yrs)",
                 age_group ~"Age",
                 gender ~ "Sex", 
                 vax_type ~ "Type of Vaccine",
                 admin_dose ~ "Dosage", 
                 onset_interval ~ "Onset interval Post-Vaccination", 
                 outcome_severity ~ "Outcome"),
    statistic = list(
      AGE_YRS ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) |>
  add_p(
    test = list(
      all_categorical() ~ "chisq.test"
    ),
    test.args = all_categorical() ~ list(simulate.p.value = TRUE)
  ) |>  bold_labels() |> as_flex_table() |>  save_as_docx(path = "Table1_Stroke_Comparison.docx")
#|>as_gt() |> gtsave(filename = "Table1_Stroke_Comparison.png", vwidth = 1800, vheight = 1200)
#GG magic
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
#Total AE Reports for COVID-19 vaccines by Year
totalreports_plot <- df_prepared |>  count(report_year)
ggplot(totalreports_plot, aes(x = factor(report_year), y = n)) +
  geom_col(fill = "steelblue") +  scale_y_continuous(labels = comma) +
  labs(title = "Total AE Reports for COVID-19 vaccines by Year",
       x = "Year",
       y = "Number of Reports") + theme_minimal()
ggsave(filename = "Total AE Reports for COVID-19 vaccines by Year.png")
totalstrokereports_plot <- df_prepared |> filter(stroke == 1) |> count(report_year)
ggplot(totalstrokereports_plot, aes(x = factor(report_year), y = n)) +
  geom_col(fill = "firebrick") +
  labs(title = "Stroke Reports for COVID-19 vaccines by Year",
       x = "Year",
       y = "Number of Stroke Reports") + theme_minimal()
ggsave(filename = "Stroke Reports for COVID-19 vaccines by Year.png")
#OR WE COULD MAKE THINGS MORE PRESENTABLE
library(dplyr)
library(ggplot2)
library(scales)
theme_set( theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()))
setwd("F:/VAERS Project/Main Files")
#1. Total AE Reports by Year
totalreports_plot <- df_prepared |>  count(report_year)
PICTURE1 <- ggplot(totalreports_plot, aes(x = factor(report_year), y = n)) +
  geom_col(fill = "#2C7FB8", width = 0.7) +  scale_y_continuous(labels = comma) +
  labs(title = "Total AE Reports by Year",
       x = "Year",
       y = "Number of Reports") +   theme(plot.title = element_text(hjust = 0.5))
ggsave("Total AE Reports by Year.png", plot = PICTURE1, width = 7, height = 5, dpi = 300)

# 2. STROKE REPORTS
totalstrokereports_plot <- df_prepared |>  filter(stroke == 1) |> count(report_year)

PICTURE2 <- ggplot(totalstrokereports_plot, aes(x = factor(report_year), y = n)) +
  geom_col(fill = "#CB3E38", width = 0.7) +
  scale_y_continuous(labels = comma) +
  labs(title = "Stroke Reports by Year",
       x = "Year",
       y = "Number of Stroke Reports")+  theme(plot.title = element_text(hjust = 0.5))
ggsave("Stroke Reports by Year.png", plot = PICTURE2, width = 7, height = 5, dpi = 300)

#3. Stroke subtypes by Year
stroke_subtype_data <- df_prepared |>
  filter(stroke == 1,  stroke_type != "No Stroke") |> count(report_year, stroke_type)
PICTURE3 <- ggplot(stroke_subtype_data, aes(x = factor(report_year), y = n, fill = stroke_type)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_y_continuous(labels = comma) +
  labs(title = "Stroke Subtypes by Year",
       x = "Year",
       y = "Number of Cases",
       fill = "Subtype") + theme(plot.title = element_text(hjust = 0.5))
ggsave("Stroke Subtypes by Year.png",  plot = PICTURE3, width = 7, height = 5, dpi = 300)
#4. Age-wise distribution
age_data <- df_prepared |>  count(age_group, stroke)
PICTURE4 <- ggplot(age_data, aes(x = age_group, y = n, fill = factor(stroke))) +
  geom_col(position = "dodge", width = 0.7) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c("0" = "#F8766D", "1" = "#00BFC4"),
    labels = c("0" = "No Stroke", "1" = "Stroke"),
    name = "Group") +
  labs(title = "Age Group Distribution",
       x = "Age Group",
       y = "Number of Reports") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Age Group Distribution.png", plot = PICTURE4, width = 7, height = 5, dpi = 300)
#5. Gender-based distribution
gender_data <- df_prepared |> filter(stroke == 1, gender != "Unknown") |> count(stroke_type, gender)
PICTURE5 <- ggplot(gender_data, aes(x = stroke_type, y = n, fill = gender)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_y_continuous(labels = comma) +
  labs(title = "Sex Distribution by Stroke Type",
       x = "Stroke Type",
       y = "Number of Cases",
       fill = "Sex") + theme(plot.title = element_text(hjust = 0.5))
ggsave("Sex Distribution by Stroke Type.png", plot = PICTURE5, width = 7, height = 5, dpi = 300)
#6. Stroke reports compared across gender and age
gender_age_data <- df_prepared |> filter(stroke == 1) |>  filter(gender %in% c("Female", "Male")) |>  count(age_group, gender)
PLOT_SEX_AGE <- ggplot(gender_age_data, aes(x = age_group, y = n, fill = gender)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) + scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c("Female" = "#E15759",
               "Male" = "#4E79A7"),
    name = "Sex") +
  labs(title = "Stroke Cases by Sex Across Age Groups",
       x = "Age Group",
       y = "Number of Stroke Reports") +
  theme_minimal(base_size = 13) +  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Stroke Cases by Sex Across Age Groups.png",PLOT_SEX_AGE, width = 9, height = 5, dpi = 300)
#7. Median days to onset
onset_data <- df_prepared |> filter(stroke == 1) |> mutate(onset = as.numeric(onset_interval)) |>
  filter(!is.na(onset), stroke_type != "No Stroke") |> group_by(stroke_type, vax_type) |>
  summarise(median_days = median(onset), .groups = "drop")
PICTURE6 <- ggplot(onset_data, aes(x = vax_type, y = median_days, fill = stroke_type)) +
  geom_col(position = position_dodge(width = 0.75),  width = 0.7) +  coord_flip() +   
  labs(title = "Median Days to Stroke Onset by Vaccine Type",
       x = "Vaccine Type",
       y = "Median Days",
       fill = "Stroke Type") + theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top")
ggsave("Median Days to Onset.png", PICTURE6, width = 8, height = 6, dpi = 300)
#8 Comparing Stroke Reporting rate across vaccines
vaccine_data <- df_prepared |>  group_by(vax_type) |> summarise(rate = mean(stroke) * 100, .groups = "drop")
PICTURE7 <- ggplot(vaccine_data, aes(x = vax_type, y = rate)) +
  geom_col(fill = "#2C7FB8", width = 0.7) +
  labs(title = "Stroke Reporting Rate by Vaccine",
       x = "Vaccine Type",
       y = "Reporting Rate (%)") + theme(plot.title = element_text(hjust = 0.5))
ggsave("Stroke Rate by Vaccine.png", plot = PICTURE7, width = 7, height = 5, dpi = 300)
###THE END