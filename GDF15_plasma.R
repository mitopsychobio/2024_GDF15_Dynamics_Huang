# GDF15 in Plasma

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
Plasma_Data <- read.csv(here::here("Data", "MiSBIE_Plasma_Data_ALL.csv"))  ## Loading all GDF15 saliva ELISA data

Plasma_Data <- Plasma_Data %>%
  group_by(ID) %>%
  mutate(Replicate = row_number()) %>%
  separate(ID, into = c("ID", "BloodDraw")) %>%
  group_by(ID, BloodDraw) %>%
  mutate(mean = mean(GDF15, na.rm = T)) %>%
  mutate(sd = sd(GDF15, na.rm = T)) %>%
  na.omit() %>%
  mutate(BloodDraw = case_when(
    BloodDraw == "F" ~ "Fasting1", 
    BloodDraw == "D1F" ~ "Fasting1",
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
    BloodDraw == "Fasting1" ~ "-20", 
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
Plasma_Data <-
  full_join(Plasma_Data, meta) %>%
  mutate(Sex = case_when(Sex =="1"~"Male",
                         Sex =="3"~"Trans",
                         Sex == "2" ~ "Female"))%>%
  mutate(Condition = case_when(Condition =="0"~"Control",
                               Condition == "1" ~ "3243A>G",
                               Condition == "2" ~ "Deletion",
                               Condition == "3" ~ "MELAS"))

Pre_age_Plasma <- Plasma_Data %>%
  select(ID, BloodDraw, mean, Age, Sex, Condition, Disease_Severity) %>%
  unique()


Data_wide <- Pre_age_Plasma %>%
  select(ID, Sex, Condition, Age, BloodDraw, mean)%>%
  unique()%>%
  pivot_wider(names_from = BloodDraw, names_prefix = "t_", values_from = mean) # format in time

write.csv(Pre_age_Plasma, "Processed_Data/GDF15_Plasma_Pre_Age_Correction.csv")
write.csv(Data_wide, "Processed_Data/GDF15_Plasma_Pre_Age_Correction_wide.csv")
rm(Data_wide)

#Age Correction
## Step 1: Get residuals from linear model: Residual = Measured_GDF15 - Predicted_GDF15

Data_Controls_Fasting <- Plasma_Data %>%
  ungroup %>%
  filter(Condition %in% "Control") %>%
  filter(BloodDraw %in% "Fasting1") %>%
  unite(ID, ID, BloodDraw, Replicate) %>%
  column_to_rownames("ID") %>%
  select(c(GDF15, Age))

Data_Controls_Fasting %>%
  ggplot(aes(x = Age, y =GDF15)) +
  geom_point() +
  xlim(18,60) +
  ylim(0,500) +
  stat_cor(aes(label = paste(..r.label..,..rr.label.., ..p.label.., sep = "~`,`~")), method = "spearman", label.x = 20, label.y = 1000, size =5) +
  geom_smooth(method = "lm") +
  ylab("Circulating GDF15 level (pg/mL)")

Data_Controls <- Plasma_Data %>%
  ungroup %>%
  filter(Condition %in% "Control") %>%
  unite(ID, ID, BloodDraw, Replicate) %>%
  column_to_rownames("ID") %>%
  select(c(Age))

Data_MD <- Plasma_Data %>%
  ungroup %>%
  filter(!Condition %in% "Control") %>%
  unite(ID, ID, BloodDraw, Replicate) %>%
  column_to_rownames("ID") %>%
  select(c(Age))

# Setup the model
mod <- lm(GDF15~Age, Data_Controls_Fasting)
summary(mod)

# Predict CTRL, all other GDF15 timepoints 
predicted_GDF15 <- predict(mod, newdata = Data_Controls)
predicted_GDF15_CTRL <- as.data.frame(predicted_GDF15) %>%
  rownames_to_column("ID")
rm(predicted_GDF15)

# Predict MD
predicted_GDF15 <- predict(mod, newdata = Data_MD)
predicted_GDF15_MD <- as.data.frame(predicted_GDF15) %>%
  rownames_to_column("ID")
rm(predicted_GDF15)

# Combine both
predicted_GDF15 = bind_rows(predicted_GDF15_MD, predicted_GDF15_CTRL)
rm(predicted_GDF15_MD, predicted_GDF15_CTRL)

tmp <- Plasma_Data %>%
  unite(ID, ID, BloodDraw, Replicate) 
final <- full_join(tmp , predicted_GDF15) %>%
  mutate(residual = GDF15 - predicted_GDF15) %>%
  separate(ID, c("ID", "BloodDraw"))
rm(Data_Controls_Fasting,  mod,  predicted_GDF15, predicted_GDF15_CTRL, predicted_GDF15_MD, Data_Controls, Data_MD,  tmp, Fasting_Data)

## Step 2: GDF15_corrected = Median_GDF15_Predicted + residual

final2 <- final %>%
  mutate(GDF15_Median = median(predicted_GDF15, na.rm = T)) %>%
  mutate(age_corrected_GDF15 = residual + GDF15_Median) #%>%
rm(final)
#select(-mean_GDF15)
# select(-mean) %>%
# rename(mean = residual)
final2 %>%
  filter(BloodDraw %in% "Fasting1") %>%
  #filter(Condition == "Control") %>%
  rename(measured_GDF15 = GDF15) %>%
  ggplot(aes(x = Age, color = ID)) +
  geom_hline(yintercept = 362, linetype = "dashed") +
  geom_point(aes(y = measured_GDF15, color = "Measured GDF15")) +
  geom_point(aes(y = age_corrected_GDF15, color = "Age-corrected GDF15")) +
  geom_point(aes(y=predicted_GDF15, color= "Predicted GDF15")) +
  scale_linetype_manual(values = 2) +
  scale_color_manual(values = c('red', 'green2', 'blue'), guide = 'legend', name = "Legend") +
  ylab("GDF15 [pg/ml]") -> p



plotly::ggplotly(p)

colnames(final2)
Plasma_Data <- final2 %>%
  select(-c(GDF15_Median, predicted_GDF15)) %>%
  rename(measured_GDF15 = GDF15) %>%
  group_by(ID, BloodDraw) %>%
  mutate(mean_measured_GDF15 = mean(measured_GDF15, na.rm = T)) %>%
  mutate(sd_measured_GDF15 = sd(measured_GDF15, na.rm = T)) %>%
  mutate(CV = sd_measured_GDF15/mean_measured_GDF15*100) %>%
  mutate(mean_residual = mean(residual, na.rm = T)) %>%
  mutate(sd_residual = sd(residual, na.rm = T)) %>%
  mutate(mean_age_corrected_GDF15 = mean(age_corrected_GDF15, na.rm = T)) %>%
  mutate(sd_age_corrected_GDF15 = sd(age_corrected_GDF15, na.rm = T)) %>%
  select(-c(measured_GDF15, age_corrected_GDF15, residual)) %>%
  unique() %>%
  #na.omit() %>%
  select(-c(mean, sd, CV))
# rm(final2)
write.csv(Plasma_Data, "Processed_Data/GDF15_age_corrected_Plasma_ALL.csv")

Data_wide <- read.csv(here::here("Processed_Data", "GDF15_age_corrected_Plasma_ALL.csv"))%>%
  select(c(ID, Sex, Condition, Age, Time, mean_age_corrected_GDF15))%>%
  unique()%>%
  pivot_wider(names_from = Time, names_prefix = "t_", values_from = mean_age_corrected_GDF15) # format in time
write.csv(Data_wide, "Processed_Data/GDF15_age_corrected_Plasma_ALL_wide.csv")


