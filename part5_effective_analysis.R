###########有效组和无效组在FMT前和后的差异菌株
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

#############数据处理
dt <- read.table("profile.txt", header = T, sep = "\t")

sample_map = read.table("multiple_FMT_Sample.txt", header = T,  sep = "\t",stringsAsFactors = F)
sample_map = sample_map[,c(2,7,9)]
sample_map <- sample_map%>%
  filter(sample_map$Type != "Multiple_2")
sample_map1 <- sample_map%>%
  filter(sample_map$Group != "post-FMT")
sample_map2 <- sample_map%>%
  filter(sample_map$Group != "pre-FMT")

result <- align_dt_sample(dt, sample_map1, ID = "Sample")
aligned_sample <- result$sample_map
aligned_dt <- result$dt
dt2 <- aligned_dt %>%
  rownames_to_column(var = "name")
#write.table(dt2,file="re_multiple_pre.txt",sep="\t", quote = F,row.names = F)
#write.table(aligned_sample,file="Sample_multiple_pre.txt",sep="\t", quote = F,row.names = F)

result <- align_dt_sample(dt, sample_map2, ID = "Sample")
aligned_sample <- result$sample_map
aligned_dt <- result$dt
dt2 <- aligned_dt %>%
  rownames_to_column(var = "name")
#write.table(dt2,file="re_multiple_post.txt",sep="\t", quote = F,row.names = F)
#write.table(aligned_sample,file="Sample_multiple_post.txt",sep="\t", quote = F,row.names = F)

################差异分析，以pre-FMT数据为例，post-FMT数据是一样的
dt <- read.table("re_multiple_pre.txt", header = T, sep = "\t")
dt <- column_to_rownames(dt,var="name")
sample_map = read.table("Sample_multiple_pre.txt", header = T,  sep = "\t",stringsAsFactors = F)

sample_map2 <- sample_map%>%
  filter(sample_map$Group != "Donor")

result <- align_dt_sample(dt, sample_map2, ID = "Sample")
aligned_sample <- result$sample_map
aligned_dt <- result$dt
dt_filtered <- aligned_dt
sample_map <- aligned_sample

library(Maaslin2)
library(dplyr)

sampf <- unique(subset(aligned_sample, Sample %in% colnames(aligned_dt), c("Sample", "Group", "Type")))
row.names(sampf) <- sampf$Sample
dtf <- as.data.frame(t(aligned_dt[, sampf$Sample]))
rownames(dtf) <- rownames(sampf) 

output_dir <- "./"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fit_data <- Maaslin2(
  input_data = dtf,
  input_metadata = sampf,
  output = output_dir,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = c("Type"),               
  random_effects = c(),                         
  reference = c("Type,Single_1"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

save_path <- file.path(output_dir, "maaslin_results_pre.RData")
save(fit_data, file = save_path)

results <- fit_data$results
significant_results <- all_day_results %>%
  mutate(
    Significance = ifelse(qval < 0.05, "Significant", "Not significant")
  )
results_sorted <- results[order(-abs(results$coef), results$pval), ]
results_sorted <- results_sorted %>%
  mutate(keystone = ifelse(feature %in% feature1$feature, "keystone", ""))

##################随机森林分析
#############不同类别的数据处理
dt <- read.table("re_multiple_pre.txt", header = T, sep = "\t")
sample_map = read.table("Sample_multiple_pre.txt", header = T,  sep = "\t",stringsAsFactors = F)

bac_dt <- aligned_dt %>%
  rownames_to_column('row_name') %>% 
  filter(grepl("^MGY", row_name)) %>% 
  column_to_rownames('row_name')  

result <- align_dt_sample(bac_dt, sample_map, ID = "Sample")
aligned_sample <- result$sample_map
aligned_dt <- result$dt

write.table(aligned_dt,file="Bacteria_abundance.txt",sep="\t", quote = F,row.names = T)
write.table(aligned_sample,file="Bacteria_group.txt",sep="\t", quote = F,row.names = F)

fun_dt <- dt_filtered %>%
  rownames_to_column('row_name') %>%
  filter(grepl("^F", row_name)) %>% 
  column_to_rownames('row_name') 

result <- align_dt_sample(fun_dt, sample_map, ID = "Sample")
aligned_sample <- result$sample_map
aligned_dt <- result$dt

write.table(aligned_dt,file="Fungi_abundance.txt",sep="\t", quote = F,row.names = T)
write.table(aligned_sample,file="Fungi_group.txt",sep="\t", quote = F,row.names = F)

vir_dt <- dt_filtered %>%
  rownames_to_column('row_name') %>%
  filter(grepl("^V", row_name)) %>%
  column_to_rownames('row_name')

result <- align_dt_sample(vir_dt, sample_map, ID = "Sample")
aligned_sample <- result$sample_map
aligned_dt <- result$dt

write.table(aligned_dt,file="Virus_abundance.txt",sep="\t", quote = F,row.names = T)
write.table(aligned_sample,file="Virus_group.txt",sep="\t", quote = F,row.names = F)

############用python脚本进行随机森林分析
############用shell脚本计算AUC值
############绘图
library(ggplot2)
library(pROC)
library(dplyr)
library(data.table)

my_join <- function(dt, group, join_str="zysplitstr"){
  dt[,group] = lapply(dt[,group],as.character)
  dt$zy_tmp_group = apply(dt[,group], 1, function(x)paste(x, collapse = join_str))
  return(dt)
}

my_split <- function(dt, column, into, sep="zysplitstr"){
  splits <- strsplit(dt[[column]], sep)
  split_dt <- do.call(rbind, splits)
  colnames(split_dt) <- into
  dt[, into] <- split_dt
  return(dt)
}

map_name <- function(roc.list){
  # 返回映射后的新名字
  oc = c() # old name
  nc = c() # new name
  for(rb in names(roc.list)){
    # b = signif(ci(roc.list[[rb]], of="auc")*100, digits=3)
    b = sprintf("%0.1f",ci(roc.list[[rb]], of="auc")*100)
    c = paste(rb, " (",b[2], "%)\t95% CI: " , b[1],"%-",b[3],"%", sep="")
    oc = c(oc,rb)
    nc = c(nc,c)
  }
  names(nc) = oc
  nc
}

calc_auc <- function(dt, pred=NA, true=NA, group=NULL,acc=F,levels=NA,
                     boot_n=2000){
  # group 支持向量
  roc.list = list()
  
  if(typeof(group) == "NULL"){
    if(sum(is.na(levels)) == 1){
      levels = unique(dt[,true])
    }
    roc.list['AUC'] = list( roc(dt[,true], dt[,pred], levels=levels))
    grps = "AUC"
    
  }else{
    ### 如果分组是向量，那就连接起来
    if(length(group) > 1){
      old_group = group
      dt = my_join(dt, group)
      group = "zy_tmp_group"
    }
    grps = unique(dt[,group])
    if(length(grps) == 1){
      if(sum(is.na(levels)) == 1){
        levels = unique(dt[,true])
      }
      roc.list['AUC'] = list(roc(dt[,true], dt[,pred], levels=levels))
    }else{
      for(g in grps){
        temp_dt = dt[dt[,group]==g,]
        if(sum(is.na(levels)) == 1){
          levels = unique(temp_dt[,true])
        }
        roc.list[as.character(g)] = list(  roc(temp_dt[,true], temp_dt[,pred], levels=levels))
      }
    }
  }
  names_ = names(roc.list)
  result_auc = matrix(NA, ncol=6,nrow=length(grps),
                      dimnames=list(names_, c("low","auc","high","low_acc","acc","high_acc")))
  for(i in names_){
    b = sprintf("%0.4f",ci(roc.list[[i]], of="auc")*100)
    result_auc[i,c(1,2,3)] = as.numeric(b)
    if(isTRUE(acc)){
      ac = ci.coords(roc.list[[i]], x="best", ret="accuracy", transpose=F)
      ac = sprintf("%0.4f",unlist(ac)*100)
      ac[2] = sprintf("%0.4f",coords(roc.list[[i]], x="best", ret="accuracy", transpose=F)*100)
      result_auc[i,c(4,5,6)] = ac
    }
  }
  
  result_auc = as.data.frame(result_auc)
  if(exists("old_group")){
    result_auc$zy_tmp_group = rownames(result_auc)
    result_auc = my_split(result_auc, group, old_group)
    result_auc$zy_tmp_group <- NULL
    rownames(result_auc) = NULL
  } else if( typeof(group) != "NULL" ){
    result_auc[,group] = rownames(result_auc)
    rownames(result_auc) = NULL
  }
  list(table=result_auc)
}

plot_roc <- function(dt, pred=NA, true=NA, group=NULL, levels=NA,
                     fill=FALSE,
                     cols = NA, conf_level=0.95, boot_n=2000){
  old_scipen = getOption("scipen")
  old_digits = getOption("digits")
  options(scipen=0)
  options(digits=7)
  
  roc.list = list()
  if(typeof(group) == "NULL"){
    if(sum(is.na(cols)) == 1){cols = "darkblue"}
    if(sum(is.na(levels)) == 1){levels = unique(dt[,true])}
    roc.list['AUC'] = list(roc(dt[,true], dt[,pred], levels=levels))
  }else{
    if(length(group) > 1){
      old_group = group
      dt = my_join(dt, group, " - ")
      group = "zy_tmp_group"
    }
    grps = unique(dt[,group])
    if(sum(is.na(cols)) == 1){
      cols=c(1:length(grps))
    }
    for(g in grps){
      if(sum(is.na(levels)) == 1){levels = unique(dt[,true])}
      temp_dt = dt[dt[,group]==g,] %>% droplevels
      roc.list[as.character(g)] = list(  roc(temp_dt[,true], temp_dt[,pred], levels=levels))
    }
  }
  new_name_map <- map_name(roc.list)
  p <- ggroc(roc.list)+
    theme_bw()+
    geom_segment(data = data.frame(x = 0, y = 1),
                 aes(x = x, y = y, xend = 1, yend = 0),
                 color = "#d9d9d9", lwd = .4, inherit.aes = F)+
    scale_color_manual(values=cols, labels=new_name_map)+
    theme(
      panel.grid.minor = element_blank()
      ,panel.grid = element_line(linetype="dashed", color="black", linewidth = 0.2)
      ,panel.border = element_rect(color="black", linewidth = 0.5)
      ,axis.ticks = element_line(color="black", linewidth = 0.5)
      ,axis.ticks.length = unit(2,"mm")
      ,axis.text  = element_text(color="black")
    )
  
  if(isTRUE(fill)){
    ci.list <- lapply(roc.list, function(rocobj)
      setDT(
        data.frame(
          ci.se(rocobj, specificities=seq(0, 1, 0.1)), check.names=F)
        ,keep.rownames = T)
    )
    data_ci <- bind_rows(ci.list, .id="plot_group")
    data_ci$rn = as.numeric(data_ci$rn)
    p <- p+
      geom_ribbon(data=data_ci,aes(x=rn, ymin=`2.5%`, ymax=`97.5%`, fill=plot_group),
                  alpha=.3,
                  inherit.aes = F)+ # 必须有参数inherit.aes
      scale_fill_manual(values=cols, labels=new_name_map)
  }
  
  options(scipen = old_scipen)
  options(digits = old_digits)
  list(plot=p, ROC=roc.list, labels=new_name_map)
}

##################数据处理
auc1 <- read.table('auc/All_auc.txt', sep = '\t', header = TRUE, fill = TRUE)
auc2 <- read.table('auc/Bacteria_auc.txt', sep = '\t', header = TRUE, fill = TRUE)
auc3 <- read.table('auc/Fungi_auc.txt', sep = '\t', header = TRUE, fill = TRUE)
auc4 <- read.table('auc/Virus_auc.txt', sep = '\t', header = TRUE, fill = TRUE)

auc_list <- list(
  auc1 = auc1,
  auc2 = auc2,
  auc3 = auc3,
  auc4 = auc4
)

############species和AUC曲线
summarize_by_X <- function(data) {
  data %>%
    group_by(X) %>%
    summarise(
      avg_auc = mean(auc, na.rm = TRUE),
      avg_high = mean(high, na.rm = TRUE),
      avg_low = mean(low, na.rm = TRUE)
    )
}

summary_list <- lapply(auc_list, summarize_by_X)
summary_list_named <- lapply(names(summary_list), function(name) {
  summary_list[[name]]$source <- name
  summary_list[[name]]
})

all_data <- bind_rows(summary_list_named)

sources_map <- c("auc1" = "All", 
                 "auc2" = "Bacteria", 
                 "auc3" = "Fungi",
                 "auc4" = "Virus")

all_data <- all_data %>%
  mutate(source = case_when(
    source %in% names(sources_map) ~ sources_map[source],
    TRUE ~ source 
  ))

keep_x <- c(1, seq(5, 100, by = 5))
filtered_data <- all_data %>%
  filter(X %in% keep_x)

filtered_data <- filtered_data %>% mutate(X = factor(X))

my_colors <- c("#49B192", "#E64B35", "#4DBBD5", "#3C5488")

filtered_data$source <- factor(filtered_data$source, levels = unique(filtered_data$source))

p <- ggplot(filtered_data, aes(x = X, y = avg_auc, color = source, group = source)) +
  
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  
  scale_color_manual(values = my_colors) +
  
  ylim(50, 100) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    aspect.ratio = 0.5
  ) +
  
  labs(
    title = "OTUs vs AUC",
    x = "Number of OTUs",
    y = "%AUC mean",
    color = "Project"
  )

print(p)

#################species重要性排序
data <- read.table("All",header = T,sep = "\t")

data2 <- data[1:30,-2]

tax <- read.table("taxonomy.txt", header = T,  sep = "\t")
data3 <- left_join(data2,tax ,by = c("Feature" = "name"))
data3 <- data3 %>%
  dplyr::select(species,Mean_Importance,Std_Importance)
names(data3)[1] <- "Feature"

data4 <- data3 %>%
  arrange(desc(Mean_Importance)) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))

ggplot(data4, aes(x = Mean_Importance, y = Feature)) +
  geom_col(aes(fill = Mean_Importance), show.legend = FALSE) +
  labs(
    title = "Top Features: Random Forest Importance",
    x = "Mean Importance",
    y = "Feature"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
    axis.ticks = element_line(colour = "black"),
    aspect.ratio = 2
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))