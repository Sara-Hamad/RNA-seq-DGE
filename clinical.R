library(dplyr)
library(ggplot2)
library(patchwork)

setwd("E:/Courses/RNA-seq/Bakry/clinical/info_clinical_normal")
clinical_data_normal <- read.delim("clinical.tsv")
 
#pre processing normal clinical data
clinical_data_normal <- clinical_data_normal %>%
  mutate(across(everything(), as.character)) %>%  # Convert all columns to character
  mutate(across(everything(), ~na_if(.x, "--")))  # Replace "--" with NA


#identify non numiric data  
non_numeric <- clinical_data_normal$days_to_death[!grepl("^[0-9.-]+$", clinical_data_normal$days_to_death)]
print(non_numeric)
  

#removing the "--" from the hybird coulmn 
clinical_data_normal <- clinical_data_normal %>%
  mutate(days_to_death = gsub("--", NA, days_to_death))

#setting coulmns with numbers into numeric 
clinical_data_normal <- clinical_data_normal %>%
  mutate(
    age_at_index = as.numeric(age_at_index),
    days_to_birth = as.numeric(days_to_birth),
    days_to_death = as.numeric(days_to_death),
    year_of_birth = as.numeric(year_of_birth),
    age_at_diagnosis = as.numeric(age_at_diagnosis),
    days_to_diagnosis = as.numeric(days_to_diagnosis),
    days_to_last_follow_up = as.numeric(days_to_last_follow_up)
  )

# Convert 'age_at_diagnosis' to years if it's currently in days
clinical_data_normal$age_at_diagnosis_years <- clinical_data_normal$age_at_diagnosis / 365.25



# Define your plots
p1 <- ggplot(clinical_data_normal, aes(x = age_at_diagnosis_years)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black") +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  labs(title = "Distribution of Age at Diagnosis _sample taken", x = "Age at Diagnosis (years)", y = "Number of Patients") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 16, face = "bold"))

p2 <- ggplot(clinical_data_normal, aes(x = factor(3), y = age_at_diagnosis_years)) +
  geom_point(stat = "identity", position = position_jitter(width = 0.1)) +
  labs(y = "Age at Diagnosis (years)", x = "", title = "Dot Plot of Age at Diagnosis") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

p3 <- ggplot(clinical_data_normal, aes(x = race, y = age_at_diagnosis_years)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  labs(title = "Age at Diagnosis by Race", x = "Race", y = "Age at Diagnosis (years)") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), plot.title = element_text(hjust = 0.5))

spacer <- plot_spacer()
# Combine plots with space between them
p_combined <- p1 + spacer + p2 + spacer + p3

# Open a PNG device
png("normal_clinical_plots_horizontal.png", width = 24, height = 8, units = 'in', res = 300)
print(p_combined)
dev.off()

#empty the workspace
rm(list =ls())

#repeat steps for tumor 
setwd("E:/Courses/RNA-seq/Bakry/clinical/info_clinical_tumor")
clinical_data_tumor <- read.delim("clinical.tsv")


clinical_data_tumor <- clinical_data_tumor %>%
  mutate(across(everything(), as.character)) %>%  # Convert all columns to character
  mutate(across(everything(), ~na_if(.x, "--")))  # Replace "--" with NA


#identify non numiric data  
non_numeric <- clinical_data_tumor$days_to_death[!grepl("^[0-9.-]+$", clinical_data_tumor$days_to_death)]
print(non_numeric)



clinical_data_tumor <- clinical_data_tumor %>%
  mutate(days_to_death = gsub("--", NA, days_to_death))%>%
  mutate(days_to_last_follow_up  = gsub("--", NA, days_to_last_follow_up ))


clinical_data_tumor <- clinical_data_tumor %>%
  mutate(
    age_at_index = as.numeric(age_at_index),
    days_to_birth = as.numeric(days_to_birth),
    days_to_death = as.numeric(days_to_death),
    year_of_birth = as.numeric(year_of_birth),
    age_at_diagnosis = as.numeric(age_at_diagnosis),
    days_to_diagnosis = as.numeric(days_to_diagnosis),
    days_to_last_follow_up = as.numeric(days_to_last_follow_up)
  )

# Convert 'age_at_diagnosis' to years if it's currently in days
clinical_data_tumor$age_at_diagnosis_years <- clinical_data_tumor$age_at_diagnosis / 365.25


#histogram of tumor samples 
p1<-ggplot(clinical_data_tumor, aes(x = age_at_diagnosis_years)) +
  geom_histogram(
    binwidth = 5, 
    fill = "blue",
    color = "black"
  ) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  labs(
    title = "Distribution of Age at Diagnosis _samples were taken_",
    x = "Age at Diagnosis (years)",
    y = "Number of Patients"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

#dot_plot 
p2 <-ggplot(clinical_data_tumor, aes(x = factor(1), y = age_at_diagnosis_years)) + 
  geom_point(stat = "identity", position = position_jitter(width = 0.1)) +
  labs(y = "Age at Diagnosis (years)", x = "", title = "Dot Plot of Age at Diagnosis") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p3 <- ggplot(data = clinical_data_tumor, aes(x =vital_status , y =age_at_diagnosis_years )) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',position = position_jitter(width = 0.1)) +
  labs(
    title = "Age at Diagnosis _samples were taken_",
    x = "vital_status",
    y = "Age at Diagnosis (years)"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

spacer <- plot_spacer()
# Combine plots with space between them
p_combined <- p1 + spacer + p2 + spacer + p3

# Open a PNG device
png("_tumor_clinical_plots_.png", width = 24, height = 8, units = 'in', res = 300)
print(p_combined)
dev.off()

