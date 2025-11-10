library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(ggpubr)
library(tibble)

sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

numFunc = function(x){
  if(x < 0.001){formatC(x, digits = 1, width = 1, format = "e", flag = "0")}
  else if(x<0.05){formatC(x, digits = 3, width = 1, format = "f", flag = "0")}
  else{NA}
}

zy_alpha = function(dt=NA, sample_map=NA, group="Group", ID="Sample", # 必须参数
                    index="shannon", # 计算参数
                    sample.color=NA, # 美化参数
                    box_width=0.5, # 箱式图宽度
                    title="alpha diversity", # 文字参数,
                    violin = F
){
  # pvalue给的是非精确计算exact=F
  ## colors 
  if (any(is.na(sample.color))){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  message(paste(length(sample.color), "of groups to plot"))
  
  ## align dt and group
  dt = dt[,sample_map[,ID]]
  dt = dt[rowSums(dt)!=0,]
  
  #alpha
  if(tolower(index) == "obs"){
    alpha = data.frame(alpha=colSums((dt>0)+0))
  }else{
    alpha = data.frame(alpha = vegan::diversity(t(dt),index=index))
  }
  
  dm = merge(alpha,sample_map, by.x='row.names', by.y=ID)
  comp = combn(as.character(unique(dm[,group])),2,list)
  
  p = ggplot(dm, aes(x=.data[[group]], y=alpha,fill=.data[[group]]))
  if(isTRUE(violin)){
    p <- p+
      geom_violin()+
      geom_boxplot(width=box_width, fill="white",
                   position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.colour = NA)
  }else{
    p <- p+ 
      geom_boxplot(position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.color="#c1c1c1")
  }
  
  ylabs = structure(c("Number of OTUs","Shannon index", "1 - Simpson index", "Invsimpson index"),
                    names=c("obs", "shannon", "simpson","invsimpson"))
  ylab = ylabs[tolower(index)]
  
  
  p <- p+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values=sample.color)+
    geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=numFunc)+
    # geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1)+
    labs(title=title, y = ylab, x=NULL)
  
  p
}

align_dt_sample <- function(dt, sample_map, ID=NA){
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dt = dt[,sample_map[,ID]] %>% filter(rowSums(.) !=0)
  list(dt=dt, sample_map=sample_map)
}

calculate_alpha_by_sample <- function(dt, sample_map) {
  # 确保只保留丰度表和 sample_map 中都存在的样本
  common_samples <- intersect(colnames(dt), sample_map$Sample)
  dt <- dt[, common_samples]
  sample_map <- sample_map[sample_map$Sample %in% common_samples, , drop = FALSE]
  
  alpha <- list()
  alpha$obs <- colSums(dt > 0)                        
  alpha$shannon <- diversity(t(dt), index = "shannon")  
  alpha$simpson <- diversity(t(dt), index = "simpson")  
  alpha$invsimpson <- diversity(t(dt), index = "invsimpson") 
  
  alpha_df <- data.frame(
    Sample = names(alpha$obs),
    obs = alpha$obs,
    shannon = alpha$shannon,
    simpson = alpha$simpson,
    invsimpson = alpha$invsimpson
  )
  
  alpha_df <- merge(alpha_df, sample_map, by.x = "Sample", by.y = "Sample")
  
  alpha_df <- alpha_df[, c("Sample", "PRJ", "Group", "obs", "shannon", "simpson", "invsimpson")]
  
  return(alpha_df)
}

###########数据处理
dt <- read.table("profile.txt", header = T, sep = "\t")

sample_map = read.table("sample.info.txt", header = T,  sep = "\t",stringsAsFactors = F)
sample_map$Timepoint[sample_map$Timepoint == "Recipient before FMT"] <- "pre-FMT"
sample_map <- sample_map %>%
  mutate(
    Day = case_when(
      Timepoint %in% c("Donor", "pre-FMT") ~ Timepoint,
      TRUE ~ as.character(Day)
    )
  )
sample_map_cleaned <- sample_map %>%
  filter(!is.na(Day))
sample_map_cleaned <- distinct(sample_map_cleaned, Sample, .keep_all = TRUE)

result <- align_dt_sample(dt2, sample_map_cleaned, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample
#write.table(dt3,file="re_timepoint.txt",sep="\t", quote = F,row.names = F)
#write.table(aligned_sample,file="Sample_timepoint.txt",sep="\t", quote = F,row.names = F)

############combat校正
dt <- read.table("re_timepoint.txt", header = T,  sep = "\t")
dt <- column_to_rownames(dt,var="name") 
ped = min(setdiff(as.numeric(as.matrix(dt)),0))/100
dt = dt+ped

sample_map = read.table("Sample_timepoint.txt", header = T,  sep = "\t",stringsAsFactors = F)
dt_log <-  asin(sqrt(dt))

mod <- model.matrix(~ Group, data = sample_map)

combat_data <- ComBat(
  dat = dt_log,
  batch = sample_map$PRJ,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)
#write.table(combat_data,file="re_timepoint_combat.txt",sep="\t", quote = F,row.names = T)

################α多样性
dt <- read.table("re_timepoint.txt", header = T, sep = "\t")
dt <- column_to_rownames(dt,var="name")
sample_map = read.table("Sample_timepoint.txt", header = T,  sep = "\t",stringsAsFactors = F)

bac_dt <- dt %>%
  rownames_to_column('row_name') %>% 
  filter(grepl("^MGY", row_name)) %>% 
  column_to_rownames('row_name') 

fun_dt <- dt %>%
  rownames_to_column('row_name') %>% 
  filter(grepl("^F", row_name)) %>% 
  column_to_rownames('row_name') 

vir_dt <- dt %>%
  rownames_to_column('row_name') %>% 
  filter(grepl("^V", row_name)) %>%
  column_to_rownames('row_name') 

alpha_result <- calculate_alpha_by_sample(dt, sample_map)
alpha_result <- alpha_result[,1:5]
write.csv(alpha_result, "raw_all_alpha_summary.csv",sep="\t", quote = F,row.names = FALSE)

#################
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(metafor)
library(readxl)
library(ggplot2)
library(pheatmap) 
library(permute)
library(lattice)
library(vegan)
all_dt <- read.table("raw_all_alpha_summary.csv", header = T,  sep = ",")
bac_dt <- read.table("raw_bac_alpha_summary.csv", header = T,  sep = ",")
fun_dt <- read.table("raw_fun_alpha_summary.csv", header = T,  sep = ",")
vir_dt <- read.table("raw_vir_alpha_summary.csv", header = T,  sep = ",")

all_dt  <- merge(all_dt,sample_map[,c(2,6)],by="Sample")
colnames(all_dt)[colnames(all_dt) == "obs"] <- "all_obs"
colnames(all_dt)[colnames(all_dt) == "shannon"] <- "all_shannon"
bac_dt  <- merge(bac_dt,sample_map[,c(2,6)],by="Sample")
colnames(bac_dt)[colnames(bac_dt) == "obs"] <- "bac_obs"
colnames(bac_dt)[colnames(bac_dt) == "shannon"] <- "bac_shannon"
fun_dt  <- merge(fun_dt,sample_map[,c(2,6)],by="Sample")
colnames(fun_dt)[colnames(fun_dt) == "obs"] <- "fun_obs"
colnames(fun_dt)[colnames(fun_dt) == "shannon"] <- "fun_shannon"
vir_dt  <- merge(vir_dt,sample_map[,c(2,6)],by="Sample")
colnames(vir_dt)[colnames(vir_dt) == "obs"] <- "vir_obs"
colnames(vir_dt)[colnames(vir_dt) == "shannon"] <- "vir_shannon"

merge_dt <- cbind(all_dt,bac_dt,fun_dt,vir_dt)
merge_dt <- merge_dt[,c(1,4,5,10,11,16,17,22,23)]

##########不较正
dt <-  merge_dt
dt <- column_to_rownames(dt,var="Sample")
dt <- t(dt)

library(Maaslin2)
library(dplyr)
sampf <- unique(subset(sample_map, Sample %in% colnames(dt), c("PRJ", "Sample", "Day")))
row.names(sampf) <- sampf$Sample
dtf <- as.data.frame(t(dt))

output_dir <- "maaslin_raw"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fit_data <- Maaslin2(
  input_data = dtf,
  input_metadata = sampf,
  output = output_dir,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = c( "Day"),            
  random_effects = c(),                     
  reference = c("Day,Donor"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

save_path <- file.path(output_dir, "all_maaslin_results.RData")
save(fit_data, file = save_path)

###########较正
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
  reference = c("Day,Donor", "PRJ,Langdon_2021"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

save_path <- file.path(output_dir, "all_maaslin_results.RData")
save(fit_data, file = save_path)


###################绘图
raw_dt <- read.table("maaslin_raw/all_results.tsv", header = T,  sep = "\t")
prj_dt <- read.table("maaslin_PRJ/all_results.tsv", header = T,  sep = "\t")
prj_dt <- prj_dt[prj_dt$metadata != "PRJ",]

df <- prj_dt %>%
  mutate(star = ifelse(qval < 0.05, "*** ", ""))

df <- df %>%
  mutate(feature = paste(feature, value, sep = "_"))

desired_order <- c("all_obs_pre-FMT","all_obs_7","all_obs_30","all_obs_60","all_obs_90" ,"all_obs_180","all_obs_360",
                   "bac_obs_pre-FMT","bac_obs_7","bac_obs_30","bac_obs_60","bac_obs_90" ,"bac_obs_180","bac_obs_360",
                   "fun_obs_pre-FMT","fun_obs_7","fun_obs_30","fun_obs_60","fun_obs_90" ,"fun_obs_180","fun_obs_360",
                   "vir_obs_pre-FMT","vir_obs_7","vir_obs_30","vir_obs_60","vir_obs_90" ,"vir_obs_180","vir_obs_360",
                   "all_shannon_pre-FMT","all_shannon_7","all_shannon_30","all_shannon_60","all_shannon_90" ,"all_shannon_180","all_shannon_360",
                   "bac_shannon_pre-FMT","bac_shannon_7","bac_shannon_30","bac_shannon_60","bac_shannon_90" ,"bac_shannon_180","bac_shannon_360",
                   "fun_shannon_pre-FMT","fun_shannon_7","fun_shannon_30","fun_shannon_60","fun_shannon_90" ,"fun_shannon_180","fun_shannon_360",
                   "vir_shannon_pre-FMT","vir_shannon_7","vir_shannon_30","vir_shannon_60","vir_shannon_90" ,"vir_shannon_180","vir_shannon_360")
df$feature <- factor(df$feature, levels = rev(desired_order))

library(ggplot2)
# 绘图
ggplot(df, aes(x = feature, y = coef)) +
  geom_col(aes(fill = ifelse(coef >= 0, "Positive", "Negative")), width = 0.8) +
  scale_fill_manual(values = c("Negative" = "#0077B5", "Positive" = "#DC143C"), guide = "none") + 
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr),
                width = 0.2, color = "gray30") +
  geom_text(aes(y = coef + sign(coef) * stderr * 2, label = star),
            size = 6, fontface = "bold", color = "red") +
  coord_flip() +
  labs(
    title = "",
    x = "Feature",
    y = "Coefficient value"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 19),
    axis.text.y = element_text(size = 17),
    axis.text.x = element_text(size = 17),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
    axis.ticks = element_line(color = "black", size = 0.5), # 添加刻度线
    axis.ticks.length = unit(0.2, "cm"), # 调整刻度线长度
    aspect.ratio = 8
  )

#################PCoA
zy_pcoa <- function(dt=NA, sample_map=NA, group=NA, ID=NA, sample.color=NULL,
                    ado_method="bray", pca_method="bray",
                    levels=0.95, star_plot=F, ellipse_plot=T,
                    cut_rate = NA, cut_num=NA,
                    title="PCoA", x=1, y=2, ados=T, mydist=NULL){
  
  if(typeof(mydist) != "NULL"){
    fmt_profile = align_dist_sample(mydist, sample_map, ID=ID)
    mydist = fmt_profile$dist
  }else{
    # 对齐profile和分组的样本名称
    fmt_profile = align_dt_sample(dt, sample_map, ID=ID)
    dt = fmt_profile$dt
  }
  sample_map = fmt_profile$sample_map
  
  if(is.finite(cut_rate) || is.finite(cut_num)){
    sample_map = filter_group(sample_map, group=group, cut_rate = cut_rate, cut_num=cut_num)
  }
  
  ## colors 
  if ( typeof(sample.color) == "NULL" ){
    sample.color = rainbow(length(unique(sample_map[,group]))) 
  }
  
  # 统计每个分组各有多少,作为新的图例
  group_summ <- sample_map %>%
    dplyr::select(all_of(group)) %>%
    dplyr::group_by(across({{group}})) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(new_label=paste(!!sym(group), " (", count, ")", sep=""))
  
  # 确保 new_label 和 group 名称长度一致
  stopifnot(length(group_summ$new_label) == length(group_summ[[group]]))
  new_label <- setNames(as.character(group_summ$new_label), as.character(group_summ[[group]]))
  
  message(paste(length(unique(sample_map[,group])), "of groups to plot"))
  
  if(typeof(mydist) == "NULL"){
    mydist = vegdist(t(dt), method = pca_method)
  }else{
    mydist = as.dist(mydist)
  }
  
  ado_r2 = ado_p = NA
  if (isTRUE(ados)){
    if(length(unique(sample_map[,group])) > 1){
      ## adonis
      ado = adonis2(mydist ~ sample_map[,group])
      ado_r2 = round(ado$R2[1], digits = 4)
      ado_p = ado$`Pr(>F)`[1]
    }
  }
  
  ## PCoA
  pcoa = cmdscale(mydist, k=min(10, dim(mydist)[1]-1), eig=T)
  eigs = signif(pcoa$eig/sum(pcoa$eig), 4)*100
  point = pcoa$points
  
  colnames(point) = paste("pcoa.", 1:ncol(point),sep="")
  
  xlab = paste("PCoA", x, " (",eigs[x],"%)", sep="")
  ylab = paste("PCoA", y, " (",eigs[y],"%)", sep="")
  substitle <- paste0("'R'^2~'='~'", ado_r2, "'~~italic('p')~'='~'", ado_p, "'") %>% 
    as.formula() %>% 
    eval()
  
  dm = merge(point, sample_map, by.x='row.names', by.y=ID)
  
  # 定义形状：圆形（16），三角形（17），方形（15）
  unique_groups <- unique(dm$Group)
  shapes <- setNames(c(16, 17, 15)[seq_along(unique_groups)], as.character(unique_groups))
  
  p1 <- ggscatter(
    data = dm,
    x = paste0("pcoa.", x),
    y = paste0("pcoa.", y),
    color = "Day",
    shape = "Group",
    star.plot = star_plot,
    ellipse.level = levels,
    ellipse = ellipse_plot
  ) +
    theme_bw() +
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
    theme(
      panel.grid = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      panel.border = element_rect(colour = "black", linewidth = 0.25)
    ) +
    scale_color_manual(values = sample.color, labels = new_label) +
    scale_shape_manual(values = shapes) + # 设置形状
    labs(x = xlab, y = ylab, title = title, subtitle = substitle)
  
  list(plot = p1, new_label = new_label, dm = dm, mydist = mydist)
}

result <- zy_pcoa(
  dt = dt,
  sample_map = sample_map,
  group = "Day",
  ID = "Sample",
  pca_method = "bray",
  ellipse_plot = FALSE,
  sample.color = NULL  
)

dm <- result$dm
p <- result$plot

day_colors <- c(
  "Donor" = "#0077B5",         
  "pre-FMT" = "#fde0dd",       
  "7" = "#fcbba1",
  "30" = "#fc9272",
  "60" = "#fb6a4a",
  "90" = "#ef3b2c",
  "180" = "#cb181d",
  "360" = "#a50f15"
)

new_label <- result$new_label
dm$Day <- factor(dm$Day, levels = c("Donor", "pre-FMT", "7", "30", "60", "90", "180", "360"))

p_new <- p +
  scale_color_manual(
    name = "Day",
    values = day_colors,
    labels = new_label[levels(factor(dm$Day))]
  ) +
  scale_fill_manual(
    name = "Day",
    values = day_colors,
    labels = new_label[levels(factor(dm$Day))]
  ) +
  coord_fixed() +
  theme(
    aspect.ratio = 1,
    text = element_text(size = 17),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

print(p_new)

#############与供体的距离
distance_matrix <- as.matrix(result$mydist)
donor_samples <- sample_map[sample_map$Day == "Donor", "Sample"]
non_donor_samples <- sample_map[sample_map$Day != "Donor", "Sample"]
sub_dist <- distance_matrix[rownames(distance_matrix) %in% non_donor_samples, 
                            colnames(distance_matrix) %in% donor_samples]
mean_distances <- apply(sub_dist, 1, mean)
mean_distances <- as.data.frame(mean_distances)
mean_distances$Sample <- row.names(mean_distances)
row.names(sample_map) <- sample_map$Sample
merge_dt <- left_join(sample_map,mean_distances,by="Sample")

plot_data <- merge_dt %>%
  filter(!is.na(mean_distances))

total_color1 <- c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                  "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
                  "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
                  "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
                  "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4",
                  "#b2df8a", "#33a02c", "#fb9a99")

desired_order <- c("Donor", "pre-FMT", "7", "30", "60", "90", "180", "360")
custom_colors <- total_color1[1:length(desired_order)]

p <- ggplot(plot_data, aes(x = factor(Day, levels = desired_order), y = mean_distances, fill = factor(Day, levels = desired_order))) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") + 
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Bray-Curtis Distance to Donor by Group",
    x = "Group",
    y = "Average Bray-Curtis
Distance to Donor"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19), 
    plot.title = element_text(size = 16, face = "bold"), 
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
    aspect.ratio = 1
  )

print(p)

#############不同people
dt <- read.table("re_timepoint.txt", header = T, sep = "\t")

sample_map  <- read.table("../../timepoint-v2.txt", header = T,  sep = "\t",stringsAsFactors = F)

sample_map$Timepoint[sample_map$Timepoint == "Recipient before FMT"] <- "pre-FMT"
sample_map <- sample_map %>%
  mutate(
    Day = case_when(
      Timepoint %in% c("Donor", "pre-FMT") ~ Timepoint,
      TRUE ~ as.character(Day)
    )
  )

merge_dt <- merge(sample_map,dt,by = "Sample")

merge_dt1 <- merge_dt[merge_dt$PRJ == "PRJNA637878",]
merge_dt2 <- merge_dt[merge_dt$PRJ == "PRJNA674880",]
merge_dt3 <- merge_dt[merge_dt$PRJ == "PRJNA701961",]
merge_dt4 <- merge_dt[merge_dt$PRJ == "PRJNA705895",]

##############多个Donor求平均值
merge_dt2.1 <- merge_dt2[merge_dt2$People != "A1", ]
subset_rows <- merge_dt2[merge_dt2$People == "A1", ]

numeric_cols <- subset_rows %>%
  select(9:ncol(subset_rows)) %>%
  mutate(across(everything(), ~ as.numeric(.)))

mean_values <- colMeans(numeric_cols, na.rm = TRUE)

metadata <- subset_rows[1, 1:8, drop = FALSE]
mean_df <- cbind(metadata, as.data.frame(t(mean_values)))

names(mean_df)[9:ncol(mean_df)] <- names(merge_dt2)[9:ncol(merge_dt2)]
mean_df_new <- mean_df[rep(1, 3), ]
row.names(mean_df_new) <- NULL
mean_df_new$People <- 2:4
mean_df_new$Sample <- 2:4

merge_dt2.2 <- rbind(merge_dt2.1,mean_df_new)

############单个Donor分配
merge_dt3.1 <- merge_dt3 %>%
  group_by(Sample) %>%
  mutate(
    Sample = if (n() > 1) paste0(Sample, ".", 1:n()) else Sample
  ) %>%
  ungroup()

#############合并
data <- rbind(merge_dt1,merge_dt2.2,merge_dt3.1,merge_dt4)

data[, 9:ncol(data)] <- lapply(data[, 9:ncol(data)], function(x) as.numeric(as.character(x)))

dt <- data[,-c(2:8)]
row.names(dt) <- NULL
dt <- column_to_rownames(dt, var="Sample")
dt <- as.data.frame(t(dt))

sample_map <- as.data.frame(data[,1:8])

result <- zy_pcoa(
  dt = dt,
  sample_map = sample_map,
  group = "Day",
  ID = "Sample",
  pca_method = "bray",
  ellipse_plot = FALSE,
  sample.color = NULL
)

distance_matrix <- as.matrix(result$mydist)
donor_samples <- sample_map[sample_map$Day == "Donor", "Sample"]
non_donor_samples <- sample_map[sample_map$Day != "Donor", "Sample"]

sub_dist <- distance_matrix[rownames(distance_matrix) %in% non_donor_samples, 
                            colnames(distance_matrix) %in% donor_samples]

mean_distances <- apply(sub_dist, 1, mean)
mean_distances <- as.data.frame(mean_distances)
mean_distances$Sample <- row.names(mean_distances)

row.names(sample_map) <- sample_map$Sample

merge_dt <- left_join(sample_map,mean_distances,by="Sample")

plot_data <- merge_dt %>%
  filter(!is.na(mean_distances))
plot_data$Day[plot_data$Day == "Recipient before FMT"] <- "pre-FMT"

total_color1 <- c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                  "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
                  "#bebada", "#fb8072", "#80b1d3",
                  "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
                  "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4",
                  "#b2df8a", "#33a02c", "#fb9a99")

desired_order <- c("pre-FMT", "7", "30", "60", "90", "180","360")

p <- ggplot(plot_data, aes(x = factor(Day, levels = desired_order), y = mean_distances, group = People, color = as.factor(People))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = setNames(total_color1[1:length(unique(plot_data$People))], unique(plot_data$People))) + # 使用 total_color1 颜色
  labs(
    x = "Day",
    y = "Mean Distance",
    color = "People ID",
    title = "Mean Distance by Day for Each Person"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(fill = NA, color = "black", size = 1),
    
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    
    plot.title = element_text(size = 16, face = "bold"),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    
    legend.position = "right",
    aspect.ratio = 1
  )

#############绘图
unique_prjs <- unique(plot_data$PRJ)
#color_palette <- c("#cab2d6", "#fccde5", "#a6cee3","#fb9a99")
color_palette <- c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")

p <- ggplot(plot_data, aes(
  x = factor(Day, levels = desired_order),
  y = mean_distances,
  group = People,
  color = PRJ
)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = color_palette) +
  labs(
    x = "Day",
    y = "Average Bray-Curtis
Distance to Donor",
    color = "PRJ",
    title = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 19),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.position = "right",
    aspect.ratio = 1
  )
