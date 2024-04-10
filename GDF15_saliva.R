# GDF15 in Saliva

# Preparation ------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(jpeg)
library(png)
library(cowplot)
library(corrplot)
library(psych)
library(ggsignif)

layout_600_300 <- theme(
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 12),
  legend.text = element_text(size = 12)
)

layout_facet <- theme(
  axis.text = element_text(size = 8),
  axis.title = element_text(size = 10),
  legend.text = element_text(size = 10)
)

# Data load  ------------------------------------------------------------------
Saliva_Data <- read.csv(here::here("Data", "MiSBIE_Saliva_Data_ALL.csv"))  ## Loading all GDF15 saliva ELISA data


Saliva_Data <- Saliva_Data %>%
  group_by(ID) %>%
  mutate(Replicate = row_number()) %>%
  separate(ID, into = c("ID", "BloodDraw")) %>%
  group_by(ID, BloodDraw) %>%
  mutate(mean = mean(GDF15, na.rm = T)) %>%
  mutate(sd = sd(GDF15, na.rm = T)) %>%
  mutate(CV = sd/mean*100) %>%
  na.omit() %>%
  mutate(BloodDraw = case_when(
    BloodDraw == "D1F" ~ "Fasting1", 
    BloodDraw == "D2F" ~ "Fasting2", 
    BloodDraw == "D1AW" ~ "Awakening1",
    BloodDraw == "D1T30" ~ "Day1T30",
    BloodDraw == "D1T45" ~ "Day1T45",
    BloodDraw == "D1BT" ~ "Day1Bedtime",
    BloodDraw == "D2AW" ~ "Awakening2",
    BloodDraw == "D2T30" ~ "Day2T30",
    BloodDraw == "D2T45" ~ "Day2T45",
    BloodDraw == "D2BT" ~ "Day2Bedtime",
    BloodDraw == "D3AW" ~ "Awakening3",
    BloodDraw == "D3T30" ~ "Day3T30",
    BloodDraw == "D3T45" ~ "Day3T45",
    BloodDraw == "D3BT" ~ "Day3Bedtime",
    BloodDraw == "CP" ~ "ColdPressor",
    BloodDraw == "MRI1" ~ "MRI1",
    BloodDraw == "MRI2" ~ "MRI2",
    BloodDraw == "MRI3" ~ "MRI3",
    BloodDraw == "SR1" ~ "Minus5",
    BloodDraw == "SR2" ~ "5",
    BloodDraw == "SR3" ~ "10",
    BloodDraw == "SR4" ~ "20",
    BloodDraw == "SR5" ~ "30",
    BloodDraw == "SR6" ~ "60",
    BloodDraw == "SR7" ~ "90",
    BloodDraw == "SR8" ~ "120"
  )) %>%
  mutate(Time = case_when(
    BloodDraw == "Awakening1" ~ "210",
    BloodDraw == "Day1T30" ~ "240",
    BloodDraw == "Day1T45" ~ "255",
    BloodDraw == "Day1Bedtime" ~ "280",
    BloodDraw == "Awakening2" ~ "300",
    BloodDraw == "Day2T30" ~ "330",
    BloodDraw == "Day2T45" ~ "345",
    BloodDraw == "Day2Bedtime" ~ "370",
    BloodDraw == "Awakening3" ~ "390",
    BloodDraw == "Day3T30" ~ "420",
    BloodDraw == "Day3T45" ~ "435",
    BloodDraw == "Day3Bedtime" ~ "460",
    BloodDraw == "ColdPressor" ~ "150", 
    BloodDraw == "MRI1" ~ "175", 
    BloodDraw == "MRI2" ~ "185",
    BloodDraw == "MRI3" ~ "195",
    BloodDraw == "Fasting1" ~ "-20", 
    BloodDraw == "Fasting2" ~ "160",
    BloodDraw == "Minus5" ~ "-5",
    BloodDraw == "5" ~ "5",
    BloodDraw == "10" ~ "10",
    BloodDraw == "20" ~ "20",
    BloodDraw == "30" ~ "30",
    BloodDraw == "60" ~ "60",
    BloodDraw == "90" ~ "90",
    BloodDraw == "120" ~ "120"
  )) %>%
  mutate(Time = as.numeric(Time))

meta <- read.csv(here::here("Data","MiSBIE_Redcap_data_2024-02-22.csv")) %>%
  mutate(ID = as.character(ID))
Saliva_Data <-
  full_join(Saliva_Data, meta) %>%
  mutate(Sex = case_when(Sex =="1"~"Male",
                         Sex =="3"~"Trans",
                         Sex == "2" ~ "Female"))%>%
  mutate(Condition = case_when(Condition =="0"~"Control",
                               Condition == "1" ~ "3243A>G",
                               Condition == "2" ~ "Deletion",
                               Condition == "3" ~ "MELAS"))

Pre_age_Saliva <- Saliva_Data %>%
  select(ID, BloodDraw, mean, Age, Sex, Condition, Disease_Severity) %>%
  unique()


Data_wide <- Pre_age_Saliva %>%
  select(ID, Sex, Condition, Age, BloodDraw, mean)%>%
  unique()%>%
  pivot_wider(names_from = BloodDraw, names_prefix = "t_", values_from = mean) # format in time

write.csv(Pre_age_Saliva, "Processed_Data/GDF15_Saliva_Pre_Age_Correction_ALL.csv")
write.csv(Data_wide, "Processed_Data/GDF15_Saliva_Pre_Age_Correction_wide_ALL.csv")
rm(Data_wide)

# Exclude data below minimum detection level: 2.0pg/mL

Saliva_Data_min_detect <- Saliva_Data %>%
  select(ID, BloodDraw, mean, Age, Sex, Condition, Disease_Severity) %>%
  unique()%>%
  filter(mean >= 2.0)

Data_wide <- Saliva_Data_min_detect %>%
  select(ID, Sex, Condition, Age, BloodDraw, mean)%>%
  unique()%>%
  pivot_wider(names_from = BloodDraw, names_prefix = "t_", values_from = mean) # format in time


write.csv(Saliva_Data_min_detect, "Processed_Data/GDF15_Saliva_Pre_Age_Correction.csv")
write.csv(Data_wide, "Processed_Data/GDF15_Saliva_Pre_Age_Correction_wide.csv")