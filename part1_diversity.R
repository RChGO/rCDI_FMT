library(tibble)

############样本处理
align_dt_sample <- function(dt, sample_map, ID=NA){
  dt = dt[, colSums(dt, na.rm=T) != 0]
  
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dt = dt[,sample_map[,ID]] %>% filter(rowSums(., na.rm=T) !=0)
  list(dt=dt, sample_map=sample_map)
}

dt <- read.table("all.profile.txt", header = T,  sep = "\t")
sample_map <- read.table("sample.info.txt", header = T,  sep = "\t")
sample_map <- sample_map[sample_map$Group != "post-FMT",]
result <- align_dt_sample(dt2, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map
#write.table(aligned_dt,file="Donor_pre-FMT_profile.txt",sep="\t", quote = F,row.names = T)
#write.table(aligned_sample,file="Donor_pre-FMT_sample.info.txt",sep="\t", quote = F,row.names =F)

####donor pre-FMT 样本
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

# 主函数
zy_alpha <- function(dt = NA, sample_map = NA, group = "Group", ID = "Sample",
                     index = "shannon", 
                     sample.color = NA,
                     box_width = 0.5, 
                     title = "Alpha Diversity",
                     violin = FALSE, 
                     show_significance = TRUE, 
                     show_points = TRUE) {  
  
  common_samples <- intersect(rownames(dt), sample_map[[ID]])
  dt <- dt[common_samples, ]
  sample_map <- sample_map[sample_map[[ID]] %in% common_samples, ]
  
  if(tolower(index) == "shannon"){
    alpha <- data.frame(alpha = colSums((dt > 0) + 0))
  } else {
    alpha <- data.frame(alpha = vegan::diversity(t(dt), index = tolower(index)))
  }
  
  dm <- merge(alpha, sample_map, by.x = row.names, by.y = ID)
  names(dm)[1] <- "Sample"
  
  comp <- combn(as.character(unique(dm[[group]])), 2, simplify = FALSE)
  
  if (is.null(sample.color) || any(is.na(sample.color))) {
    sample.color <- brewer.pal(n = length(unique(dm[[group]])), name = "Set3")
  }
  
  message(paste(length(sample.color), "组别用于绘图"))
  
  p <- ggplot(dm, aes_string(x = group, y = "alpha", fill = group)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = sample.color)
  
  if (isTRUE(violin)) {
    p <- p +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = box_width, fill = "white",
                   position = position_dodge2(preserve = 'single'),
                   outlier.shape = NA, outlier.colour = NA)
  } else {
    p <- p + 
      geom_boxplot(position = position_dodge2(preserve = 'single'),
                   outlier.shape = 21, outlier.fill = NA, outlier.color = "#c1c1c1")
  }
  
  if (isTRUE(show_points)) {
    p <- p + geom_jitter(aes(color = .data[[group]]), width = 0.2, size = 2.5, alpha = 0.7)
  }
  
  ylabs <- structure(c("Number of OTUs","Shannon index", "1 - Simpson index", "Invsimpson index"),
                     names=c("obs", "shannon", "simpson","invsimpson"))
  ylab <- ylabs[tolower(index)]
  
  if (isTRUE(show_significance)) {
    p <- p +
      stat_compare_means(comparisons = comp,
                         method = "wilcox.test",
                         label = label_func,
                         vline.xjust = 0.75,
                         step_increase = 0.08) +
      stat_compare_means(label.y = max(dm$alpha) * 1.1,
                         method = "wilcox.test",
                         label = label_func)
  }
  
  p <- p +
    labs(title = title, y = ylab, x = NULL)
  
  print(p)
}

calculate_alpha_by_sample <- function(dt, sample_map) {
  common_samples <- intersect(colnames(dt), sample_map$Sample)
  dt <- dt[, common_samples]
  sample_map <- sample_map[sample_map$Sample %in% common_samples, , drop = FALSE]
  
  alpha <- list()
  alpha$obs <- colSums(dt > 0)                          # 观测 OTU 数
  alpha$shannon <- diversity(t(dt), index = "shannon")   # 香农指数
  alpha$simpson <- diversity(t(dt), index = "simpson")   # 辛普森指数
  alpha$invsimpson <- diversity(t(dt), index = "invsimpson")  # 逆辛普森指数
  
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

#############计算多样性指数
dt <- read.table("Donor_pre-FMT_profile.txt", header = T,  sep = "\t")
sample_map = read.table("Donor_pre-FMT_sample.info.txt", header = T,  sep = "\t")

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
#write.csv(alpha_result, "all_alpha_summary.csv",sep="\t", quote = F,row.names = FALSE)

##########绘图
all_dt <- read.table("all_alpha_summary.csv", header = T,  sep = ",")
bac_dt <- read.table("bac_alpha_summary.csv", header = T,  sep = ",")
fun_dt <- read.table("fun_alpha_summary.csv", header = T,  sep = ",")
vir_dt <- read.table("vir_alpha_summary.csv", header = T,  sep = ",")

all_dt  <- merge(all_dt,sample_map[,c(5,8)],by="Sample")
colnames(all_dt)[colnames(all_dt) == "obs"] <- "all_obs"
colnames(all_dt)[colnames(all_dt) == "shannon"] <- "all_shannon"
bac_dt  <- merge(bac_dt,sample_map[,c(5,8)],by="Sample")
colnames(bac_dt)[colnames(bac_dt) == "obs"] <- "bac_obs"
colnames(bac_dt)[colnames(bac_dt) == "shannon"] <- "bac_shannon"
fun_dt  <- merge(fun_dt,sample_map[,c(5,8)],by="Sample")
colnames(fun_dt)[colnames(fun_dt) == "obs"] <- "fun_obs"
colnames(fun_dt)[colnames(fun_dt) == "shannon"] <- "fun_shannon"
vir_dt  <- merge(vir_dt,sample_map[,c(5,8)],by="Sample")
colnames(vir_dt)[colnames(vir_dt) == "obs"] <- "vir_obs"
colnames(vir_dt)[colnames(vir_dt) == "shannon"] <- "vir_shannon"

merge_dt <- cbind(all_dt,bac_dt,fun_dt,vir_dt)
merge_dt2 <- merge_dt[,c(1,2,3,4,5,10,11,16,17,22,23)]

##############α多样性指数meta分析
myread <- function (inf){
  dt = fread(inf, sep="\t", header=T, check.names=F, data.table = F)
  row.names(dt) = dt[,1]
  dt = dt[,-1]
  dt
}


align_dt = function(x,y){
  a = intersect(rownames(x), y)
  a1 = x[a,]
  
  b = setdiff(y, a)
  
  z =x[b,]
  rownames(z) = b
  z[is.na(z)] = 0
  
  rbind(a1,z)
}

library(dplyr)
library(metafor)
library(foreach)
library(doParallel)

meta_metafor_parallel <- function(dt, group = "group", group_pair = c("Disease", "Control"),
                                  proj = "proj", measure = "SMD", method = "REML",
                                  ncpus = 1) {
  dt <- dt %>%
    rename(proj = all_of(proj), group = all_of(group))
  
  feature <- setdiff(colnames(dt), c("sample", "group", "proj"))
  
  results_list <- parallel::mclapply(feature, function(i) {
    cat("\rProcessing: ", i)
    
    tib <- dt %>%
      subset(group %in% group_pair[1]) %>%
      dplyr::select(all_of(i), proj) %>%
      rename(index = all_of(i)) %>%
      group_by(proj) %>%
      summarise(d_Mean = mean(index), d_Sd = sd(index), d_N = n(), .groups = 'drop')
    
    tib2 <- dt %>%
      subset(group %in% group_pair[2]) %>%
      dplyr::select(all_of(i), proj) %>%
      rename(index = all_of(i)) %>%
      group_by(proj) %>%
      summarise(c_Mean = mean(index), c_Sd = sd(index), c_N = n(), .groups = 'drop')
    
    meta_in <- merge(tib, tib2, by = "proj")
    
    smd_meta <- escalc(measure = measure, data = meta_in, append = TRUE,
                       m1i = d_Mean, m2i = c_Mean,
                       sd1i = d_Sd, sd2i = c_Sd,
                       n1i = d_N, n2i = c_N)
    
    non_na <- smd_meta %>% dplyr::filter(!is.na(yi)) %>% nrow()
    if (non_na != 0) {
      smd_rma <- tryCatch({
        rma(yi, vi, method = method, data = smd_meta)
      }, error = function(e){
        message("默认参数不能收敛，重新设置 for ", i)
        rma(yi, vi, method = method, data = smd_meta, control = list(stepadj = 0.5, maxiter = 1000))
      })
      
      out_df <- smd_rma$data %>%
        add_column(
          measure = measure,
          model = "RM",
          method_tau2 = method,
          val_tau2 = as.numeric(smd_rma$tau2),
          I2 = paste0(round(smd_rma$I2, digits = 2), "%"),
          Q = smd_rma$QE,
          Q_pval = smd_rma$QEp,
          feature = i,
          estimate = as.numeric(smd_rma$beta),
          ci_lb = smd_rma$ci.lb,
          ci_ub = smd_rma$ci.ub,
          pval = smd_rma$pval,
          .before = 1
        )
      return(out_df)
    } else {
      message("\nWarning: ", i, " is not suitable.")
      return(NULL)
    }
  }, mc.cores = ncpus)
  
  meta_outp <- do.call(rbind, Filter(Negate(is.null), results_list))
  
  cat(paste0("== estimate > 0, ==> ", group_pair[1]," ==\n== estimate < 0, ==> ",group_pair[2], " =="))
  return(meta_outp)
}

res = meta_metafor_parallel(merge_dt2, group="Group", group_pair=c("Donor", "pre-FMT"), proj="PRJ", measure = "SMD",method = "REML",ncpus = 40)
res.qval <- res %>%
  dplyr::select(feature, pval) %>%
  unique() %>%
  mutate(qval = p.adjust(pval, method="BH"))

resf <- res %>%
  merge(res.qval[,c("feature","qval")], by='feature') %>%
  mutate(I2 = as.numeric(gsub("%","", I2)))
#write.table(resf, "meta.analysis_α.tsv",sep = "\t", row.names = F, quote = FALSE)

forest_data <- resf %>%
  select(feature, proj, estimate, ci_lb, ci_ub, yi, vi, qval)

forest_plot <- forest_data %>%
  ggplot(aes(x = yi, 
             y = factor(feature, levels = rev(c("all_obs", "bac_obs", "fun_obs", "vir_obs",
                                                "all_shannon", "bac_shannon", "fun_shannon", "vir_shannon"))))) +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub),
                 height = 0.15, color = "black", size = 0.5) +
  geom_point(aes(shape = proj), size = 3, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  theme_minimal() +
  labs(title = "Forest Plot",
       x = "Hedges' g (Standardized Mean Difference)",
       y = "Feature") +
  scale_shape_manual(values = 1:9, name = "Study") +
  theme(
    axis.text.y = element_text(size = 19, hjust = 1),
    axis.text.x = element_text(size = 19),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 14), 
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1, linetype = "dashed")
  )+
  coord_fixed(ratio = 1.5)

forest_plot 


#####################################β多样性
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(readxl)
options(warn=-1)

filter_group <- function(sample_map, group=NA, cut_rate = NA, cut_num=NA){
  cut_off = cut_num
  if(is.finite(cut_rate)){
    if(cut_rate>1 || cut_rate < 0){
      stopifnot("cut_rate should range(0,1)" = 1)
    }else{
      cut_off = sample_map %>%
        dplyr::select(!!sym(group)) %>%
        dplyr::group_by(!!sym(group)) %>%
        dplyr::summarise(count=n()) %>%
        pull(count) %>%
        max() * cut_rate
      cut_off = as.integer(cut_off)
    }
  }
  
  select_grps <- sample_map %>%
    dplyr::select(!!sym(group)) %>%
    dplyr::group_by(!!sym(group)) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::filter(count > cut_off) %>%
    dplyr::pull(!!sym(group))
  message("???????????")
  sample_map %>%
    dplyr::mutate( !!sym(group) := ifelse( !!sym(group) %in% select_grps, !!sym(group), paste("less_",cut_off, sep="")))
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

align_dist_sample <- function(dist, sample_map, ID=NA){
  dist = as.matrix(dist)
  intersect_id = intersect(sample_map[,ID], colnames(dist))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dist = dist[ sample_map[,ID], sample_map[,ID] ]
  list(dist=as.dist(dist), sample_map=sample_map)
}

zy_pcoa <- function(dt=NA, sample_map=NA, group=NA, ID=NA, sample.color=NULL,
                    ado_method="bray", pca_method="bray",
                    levels=0.95, star_plot=F, ellipse_plot=T,
                    cut_rate = NA, cut_num=NA,
                    title="PCoA", x=1, y=2, ados=T, mydist=NULL){
  
  if(typeof(mydist) != "NULL"){
    fmt_profile = align_dist_sample(mydist, sample_map, ID=ID)
    mydist = fmt_profile$dist
  }else{
    fmt_profile = align_dt_sample(dt, sample_map, ID=ID)
    dt = fmt_profile$dt
  }
  sample_map = fmt_profile$sample_map
  
  if(is.finite(cut_rate) || is.finite(cut_num)){
    sample_map = filter_group(sample_map, group=group, cut_rate = cut_rate, cut_num=cut_num)
  }
  
  if ( typeof(sample.color) == "NULL" ){
    sample.color = rainbow(length(unique(sample_map[,group]))) 
  }
  
  group_summ <- sample_map %>%
    dplyr::select(all_of(group)) %>%
    dplyr::group_by(across({{group}})) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(new_label=paste(!!sym(group), " (", count, ")", sep=""))
  
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
  
  if (typeof(sample.color) == "NULL"){
    sample.color = rainbow(length(unique(sample_map[,group]))) 
  } 
  unique_prjs <- unique(dm$PRJ)
  n_prj <- length(unique_prjs)
  
  base_shapes <- c(15,0, 17,2,16, 1, 5, 6, 3, 4) 
  
  if (n_prj > length(base_shapes)) {
    stop("Too many PRJ groups! Only support up to 10 different shapes.")
  }
  selected_shapes <- base_shapes[seq_len(n_prj)]
  shapes <- setNames(selected_shapes, as.character(unique_prjs))
  
  p1 <- ggscatter(
    data = dm,
    x = paste0("pcoa.", x),
    y = paste0("pcoa.", y),
    color = "Group",  
    shape = "PRJ",  
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
    scale_shape_manual(values = shapes) +  
    labs(x = xlab, y = ylab, title = title, subtitle = substitle)
  
  list(plot = p1, new_label = new_label, dm = dm, dist_matrix = mydist)
}

result <- zy_pcoa(
  dt = dt,
  sample_map = sample_map,
  group = "PRJ",
  ID = "Sample",
  pca_method = "bray",
  ellipse_plot = FALSE,
  sample.color = total_color1
)
save(result,file="all_PCoA.Rdata")

# 假设 result$plot 是一个 ggplot 对象
p_main <- result$plot +
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        aspect.ratio = 1 )+
  labs(
    color = "Group", 
    shape = "Project" 
  )+
  scale_shape_manual(values = 1:7) 

###############adonis
unique_prjs <- unique(sample_map$PRJ)

results_list <- list()

for (prj in unique_prjs) {
  prj_indicator <- dm$PRJ == prj
  n_prj <- sum(prj_indicator)
  if (n_prj < 2) {
    warning(paste("PRJ =", prj, "has < 2 samples, skipped."))
    next
  }
  
  dist_sub <- dist_matrix[prj_indicator, prj_indicator]
  group_sub <- dm$Group[prj_indicator]  # 注意：Group 是分组变量（如 Treatment）
  
  if (length(unique(group_sub)) < 2) {
    warning(paste("PRJ =", prj, "has only one Group level, skipped."))
    next
  }
  
  ado_prj <- adonis2(dist_sub ~ group_sub, permutations = 999)
  
  # 提取结果
  r2_prj <- ado_prj$R2[1]  
  p_prj <- ado_prj$`Pr(>F)`[1]
  df_model <- ado_prj$Df[1] 
  
  adj_r2_prj <- RsquareAdj(r2_prj, m = df_model, n = n_prj)
  
  results_list[[prj]] <- list(
    R_squared = r2_prj,
    Adjusted_R_squared = adj_r2_prj,
    P_value = p_prj
  )
}

results_df <- tibble(
  PRJ = names(results_list),
  R_Squared = sapply(results_list, function(x) x$R_squared),
  Adj_R2 = sapply(results_list, function(x) x$Adjusted_R_squared),
  P_value = sapply(results_list, function(x) x$P_value)
)
print(results_df)

my_colors <- c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
               "#6a3d9a", "#ffff99", "#b15928","#8dd3c7")

if(length(my_colors) < length(unique(results_df$PRJ))) {
  stop("Color vector is not long enough to cover all PRJ groups.")
}
named_colors <- setNames(my_colors, unique(as.character(results_df$PRJ)))

ggplot(results_df, aes(x = PRJ, y = Adj_R2, fill = PRJ)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0("P = ", sprintf("%.3f", P_value))), 
            vjust = -0.5, color = "red") +
  labs(
    title = "Adjusted R² by PRJ",
    x = "Project (PRJ)",
    y = "Adjusted R-squared"
  ) +
  ylim(0, 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 14, hjust = 0.5),  
    aspect.ratio = 1,
    legend.position = "none" ,     
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()  
  ) +
  scale_fill_manual(values = named_colors)

###############adonis计算R2和adjR2
load("all_PCoA.Rdata")
dist_matrix <- result$dist_matrix
dist_matrix <- as.matrix(dist_matrix)
dm <- result$dm 
rownames(dm) <- NULL
dm <- column_to_rownames(dm,var="Row.names")
common_samples <- intersect(rownames(dist_matrix), rownames(dm))
common_samples <- sort(common_samples)
dist_matrix <- dist_matrix[common_samples, common_samples]
dm <- dm[common_samples, ]
identical(rownames(dist_matrix), rownames(dm))
dm <- rownames_to_column(dm, var="Row.names") 

ado_group <- adonis2(dist_matrix[sample_map$Sample,sample_map$Sample] ~ Group, data = sample_map, permutations = 999, distance = matrix)
r2 <- ado_group$R2[1] 
df_model <- ado_group$Df[1]
n <- nrow(dm)
adj_r2 <- RsquareAdj(r2, m = df_model, n = n)

ado_prj <- adonis2(dist_matrix[sample_map$Sample,sample_map$Sample] ~ PRJ, data = sample_map, permutations = 999, distance = matrix)
r2 <- ado_prj$R2[1] 
df_model <- ado_prj$Df[1]
n <- sum(ado_prj$Df) 
n <- nrow(dm) 
adj_r2 <- RsquareAdj(r2, m = df_model, n = n)
