library(tibble)
dt <- read.table("re_timepoint.txt", header = T, sep = "\t")

data <- read.table("meta.analysis_0.0002_0.25_unique.tsv",header = T,sep = "\t")
tax <- read.table("abundance_Group_tax.txt", header = T,  sep = "\t")
feature <- merge(tax,data,by = "feature")
feature = as.data.frame(feature)
feature <- feature%>%
  select(feature,species,enriched)

feature1 <- feature[feature$enriched == "Control",]
feature1$enriched[feature1$enriched == "Control"] <- "Donor"
feature2 <- feature[feature$enriched == "Case",]
feature2$enriched[feature2$enriched == "Case"] <- "pre-FMT"

sample_map = read.table("Sample_timepoint.txt", header = T,  sep = "\t",stringsAsFactors = F)
sample_map$Timepoint[sample_map$Timepoint == "Recipient before FMT"] <- "pre-FMT"

sample_map <- sample_map %>%
  mutate(
    Day = case_when(
      Timepoint %in% c("Donor", "pre-FMT") ~ Timepoint,
      TRUE ~ as.character(Day)
    )
  )

merge_data <- merge(feature,dt,by.x = "feature",by.y = "name")
#write.table(merge_data,file="re_keystone_timepoint.txt",sep="\t", quote = F,row.names = F)

###########post-FMT和donor的差异
library(Maaslin2)
library(dplyr)

dt <- read.table("re_keystone_timepoint.txt", header = T, sep = "\t")
dt <- dt[,-c(2:3)]
dt <- column_to_rownames(dt,var="feature")

sampf <- unique(subset(sample_map, Sample %in% colnames(dt), c("PRJ", "Sample", "Day")))
rownames(sampf) <- sampf$Sample
dtf <- as.data.frame(t(dt[, sampf$Sample]))
rownames(dtf) <- rownames(sampf) 

output_dir <- "maaslin_PRJ"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fit_data <- Maaslin2(
  input_data = dtf,
  input_metadata = sampf,
  output = output_dir,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = c("PRJ", "Day"), 
  random_effects = c(), 
  reference = c("Day,Donor", "PRJ,PRJNA674880"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

save_path <- file.path(output_dir, "maaslin_results_PRJ.RData")
save(fit_data, file = save_path)

###########
load("maaslin_PRJ/maaslin_results_PRJ.RData")
results <- fit_data$results 

all_day_results <- results[results$metadata=="Day",]

significant_results <- all_day_results %>%
  mutate(
    Significance = ifelse(qval < 0.05, "Significant", "Not significant")
  )

#############Donor组富集
Donor_dt <- significant_results %>%
  filter(feature %in% feature1$feature)

Donor_dt <- Donor_dt %>%
  mutate(
    enriched = case_when(
      coef < 0 ~ "Control",
      coef > 0 ~ "Case",
      TRUE ~ "Natural"
    )
  )

Donor_control <- Donor_dt[Donor_dt$enriched == "Control",]
donor_signif_summary <- Donor_control %>%
  group_by(value, Significance) %>%
  summarise(n = n(), .groups = 'drop')

donor_wide <- donor_signif_summary %>%
  pivot_wider(
    names_from = value,
    values_from = n,
    values_fill = 0
  )%>%
  select(Significance, `pre-FMT`, `7`, `30`, `60`, `90`, `180`, `360`)

donor_wide <- donor_wide %>%
  mutate(
    Significance = recode(
      Significance,
      "Not significant" = "Control_enriched_Not_significant",
      "Significant" = "Control_enriched_Significant"
    )
  )
donor_wide_control <- donor_wide

###########
Donor_case <- Donor_dt[Donor_dt$enriched == "Case",]
donor_signif_summary <- Donor_case %>%
  group_by(value, Significance) %>%
  summarise(n = n(), .groups = 'drop')

donor_wide <- donor_signif_summary %>%
  pivot_wider(
    names_from = value,
    values_from = n,
    values_fill = 0
  )
donor_wide <- donor_wide %>%
  mutate(`pre-FMT` = 0)

donor_wide <- donor_wide %>%
  mutate(
    Significance = recode(
      Significance,
      "Not significant" = "Case_enriched_Not_significant",
      "Significant" = "Case_enriched_Significant"
    )
  )
donor_wide_case <- donor_wide

combined_dt1 <- rbind(donor_wide_control,donor_wide_case)

###################pre组富集
pre_dt <- significant_results %>%
  filter(feature %in% feature2$feature)

pre_dt <- pre_dt %>%
  mutate(
    enriched = case_when(
      coef < 0 ~ "Control",
      coef > 0 ~ "Case",
      TRUE ~ "Natural"
    )
  )

pre_control <- pre_dt[pre_dt$enriched == "Control",]
pre_signif_summary <- pre_control %>%
  group_by(value, Significance) %>%
  summarise(n = n(), .groups = 'drop')

pre_wide <- pre_signif_summary %>%
  pivot_wider(
    names_from = value,
    values_from = n,
    values_fill = 0
  )%>%
  select(Significance, `pre-FMT`, `7`, `30`, `60`, `90`, `180`, `360`)

pre_wide <- pre_wide %>%
  mutate(
    Significance = recode(
      Significance,
      "Not significant" = "Control_enriched_Not_significant",
      "Significant" = "Control_enriched_Significant"
    )
  )
pre_wide_control <- pre_wide

###########
pre_case <- pre_dt[pre_dt$enriched == "Case",]
pre_signif_summary <- pre_case %>%
  group_by(value, Significance) %>%
  summarise(n = n(), .groups = 'drop')

pre_wide <- pre_signif_summary %>%
  pivot_wider(
    names_from = value,
    values_from = n,
    values_fill = 0
  )%>%
  select(Significance, `pre-FMT`, `7`, `30`, `60`, `90`, `180`, `360`)

pre_wide <- pre_wide %>%
  mutate(
    Significance = recode(
      Significance,
      "Not significant" = "Case_enriched_Not_significant",
      "Significant" = "Case_enriched_Significant"
    )
  )
pre_wide_case <- pre_wide

combined_dt2 <- rbind(pre_wide_control,pre_wide_case)

desired_order <- c(
  "Control_enriched_Significant",
  "Control_enriched_Not_significant",
  "Case_enriched_Not_significant",
  "Case_enriched_Significant"
)

combined_dt_ordered1 <- combined_dt1 %>%
  mutate(Significance = factor(Significance, levels = desired_order)) %>%
  arrange(Significance)
combined_dt_ordered2 <- combined_dt2 %>%
  mutate(Significance = factor(Significance, levels = desired_order)) %>%
  arrange(Significance)

combined_long1 <- combined_dt_ordered1 %>%
  pivot_longer(cols = `pre-FMT`:`360`, names_to = "Day", values_to = "Count")

combined_long2 <- combined_dt_ordered2 %>%
  pivot_longer(cols = `pre-FMT`:`360`, names_to = "Day", values_to = "Count")

ggplot(combined_long1, aes(x = factor(Day, levels = c("pre-FMT", "7", "30", "60", "90", "180", "360")), y = Count, fill = Significance)) +
  geom_bar(stat = "identity") +
  labs(
    title = "",
    x = "Time point",
    y = "Count",
    fill = "Significance"
  ) +
  scale_fill_manual(
    name = "Group",
    values = c(
      "Control_enriched_Significant" = "#0077B5",
      "Control_enriched_Not_significant" = "lightblue",
      "Case_enriched_Not_significant" = "#ffcccc",
      "Case_enriched_Significant" = "#DC143C"
    ),
    breaks = desired_order
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 17)
    axis.text.y = element_text(size = 17),
    axis.ticks = element_line(color = "black"), 
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20), 
    aspect.ratio = 1,
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )


##############
significant_results1 <- significant_results[,c(1,3,4)]#coef
significant_results2 <- significant_results[,c(1,3,8)]#qvalue

tax <- read.table("../../../02.taxon/taxonomy.txt", header = T,  sep = "\t")

##########
merged_data <- merge(significant_results1,tax[,c(1,8)], by.x="feature",by.y="name")
merged_data <- merged_data%>%
  select(feature,value,coef)
merged_data_wide <- merged_data %>%
  pivot_wider(
    names_from = value,
    values_from = coef
  )
merged_data_wide <- merged_data_wide %>%
  select(feature, `pre-FMT`, `7`, `30`, `60`, `90`, `180`, `360`)
#write.table(merged_data_wide,file="keystone_coef.txt", sep="\t", quote = F,row.names = F)

#########
merged_data <- merge(significant_results2,tax[,c(1,8)], by.x="feature",by.y="name")
merged_data <- merged_data%>%
  select(feature,value,qval)
value_order <- c("pre-FMT", "7", "30", "60", "90", "180", "360")
merged_data$value <- factor(merged_data$value, levels = value_order, ordered = TRUE)
sorted_data <- merged_data %>%
  arrange(feature, value)

sorted_data$signif <- ifelse(
  sorted_data$qval < 0.05,
  '3,10,#000000,1,-1',
  NA_character_ 
)
#write.table(sorted_data,file="keystone_qval_sorted.txt", sep="\t", quote = F,row.names = F)