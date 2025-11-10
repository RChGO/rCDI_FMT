###############数据处理
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

dt <- read.table("Donor_pre-FMT_profile.txt", header = T,  sep = "\t")
sample_map <- read.table("Donor_pre-FMT_sample.info.txt", header = T,  sep = "\t")

dt2 <- aligned_dt %>%
  rownames_to_column(var = "name")

dt_long <- dt %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Abundance")

dt_grouped <- dt_long %>%
  left_join(sample_map, by = c("Sample" = "Sample"))

mean_abundance <- dt_grouped %>%
  group_by(name, Group) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

mean_abundance_wide <- mean_abundance %>%
  pivot_wider(names_from = Group, values_from = MeanAbundance) %>%
  replace(is.na(.), 0) 

filtered_names1 <- mean_abundance_wide %>%
  filter(Donor > 0.0002 | `pre-FMT` > 0.0002) %>%
  pull(name)

prevalence <- dt_grouped %>%
  group_by(name, Group) %>%
  summarise(
    total_samples = n(),
    present_samples = sum(Abundance > 0, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(total_samples, present_samples),
    names_sep = "_"
  ) %>%
  mutate(
    Donor_prevalence = present_samples_Donor / total_samples_Donor,
    preFMT_prevalence = `present_samples_pre-FMT` / `total_samples_pre-FMT`
  )

tax <- read.table("taxonomy.txt", header = T,  sep = "\t")
prevalence_dt <- merge(prevalence,tax,by="name")
prevalence_dt2 <- prevalence_dt %>%
  filter(Donor_prevalence > 0.25 | `preFMT_prevalence` > 0.25)
filtered_names2 <- prevalence_dt2$name
filtered_names <- intersect(filtered_names1, filtered_names2)
dt_filtered <- dt2 %>%
  filter(name %in% filtered_names)
#write.table(dt_filtered,file = "re_filtered.0.0002_0.25.txt", sep = "\t", row.names = TRUE, quote = FALSE)

###############meta.analysis
library(data.table)
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

genome = read.table("taxonomy.txt", sep="\t", header=T, check.names=F)
dt <- read.table("re_filtered.0.0002_0.25.txt",header = T,row.names = 1,sep = "\t")
rownames(dt) <- NULL
dt <- column_to_rownames(dt,var="name")
dt_names <- rownames(dt)
genome <- genome[genome$name %in% dt_names, ]
genome <- as.data.frame(genome)

dt.tmp = align_dt(dt, genome$name)
sample = read.table("Donor_pre-FMT_sample.info.txt", header = T,  sep = "\t")

sample = subset(sample, Sample %in% colnames(dt.tmp))

numeric_cols <- sapply(dt.tmp, is.numeric)

dtt <-  t(asin(sqrt(dt.tmp)))
dm = merge(dtt, sample[,c("Sample", "Group", "PRJ")], by.x='row.names', by.y='Sample')
rownames(dm) <- make.unique(as.character(dm[,1]))
dm = dm[,-1]

res = meta_metafor_parallel(dm, group="Group", group_pair=c("Donor", "pre-FMT"), proj="PRJ", measure = "SMD",method = "REML",ncpus = 40)

library(dplyr)

res.qval <- res %>%
  dplyr::select(feature, pval) %>%
  unique() %>%
  mutate(qval = p.adjust(pval, method="BH"))

resf <- res %>%
  merge(res.qval[,c("feature","qval")], by='feature') %>%
  mutate(I2 = as.numeric(gsub("%","", I2)))

resf$feature = factor(resf$feature, levels=unique(resf$feature))
#write.table(resf, "meta.analysis_0.0002_0.25.tsv",sep = "\t", row.names = F, quote = FALSE)

resf_filter <- resf %>%
  filter(I2 <= 50, qval <= 0.05) %>%
  group_by(feature) %>%
  filter((ci_lb > 0 & ci_ub > 0) | (ci_lb < 0 & ci_ub < 0))
#write.table(resf, "meta.analysis_0.0002_0.25_filter",sep = "\t", row.names = F, quote = FALSE)

############wilcoxon
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(parallel)

zy_pvalue_one_vs_others = function(dt=NA, sample_map=NA, group="Group", ID="Sample",p.method="wilcox.test"){
  grps = unique(sample_map[,group])
  result = rbind()
  for(g in grps){
    sample_map = sample_map %>% mutate(zy_temp_group=ifelse(get(`group`)==g,get(`group`),"other"))
    temp_result = zy_pvalue(dt, sample_map, group="zy_temp_group", ID=ID, p.method=p.method)
    temp_result$qvalue = p.adjust(temp_result$pvalue, method = "BH")
    result = rbind(result,temp_result)
  }
  
  return(result)
}

my_siglabs <- function(dt, name='name', lab1='g1', lab2='g2'){
  features = unique(dt[,name])
  res = rbind()
  for( feature in features){
    dtf = dt[dt[,name] == feature, ]
    labs = unique(unlist(dtf[,c(lab1, lab2)]))
    for(lab in labs){
      sig_labs = unique(unlist(dtf[dtf[,lab1]== lab | dtf[,lab2]==lab, c(lab1, lab2)]))
      tolabs = paste(setdiff(sig_labs, lab), collapse=",")
      tmp <- data.frame(name=feature, from=lab, to=tolabs)
      res = rbind(res, tmp)
    }
  }
  res
}

zy_pvalue = function(dt=NA, sample_map=NA, group="Group", ID="Sample", p.method="wilcox.test"){    
  test.arg = c(wilcox.test, t.test)
  names(test.arg) = c("wilcox.test", "t.test")
  my_test = test.arg[[p.method]]
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }else{
    message("\033[31mInfo\t数据和分组完全匹配\033[0m")
  }
  dt = dt[, sample_map[,ID]]
  
  raw_ncol = ncol(dt)
  raw_nrow = nrow(dt)
  dt = dt[rowSums(dt)!=0,]
  f_ncol = ncol(dt)
  f_nrow = nrow(dt)
  
  if(f_ncol != raw_ncol || raw_nrow != f_nrow){
    message(paste("delete all items is 0 -> columns:", raw_ncol - f_ncol, " rows:" ,raw_nrow - f_nrow, sep=""))
  }
  
  grps = unique(sample_map[,group])
  com = t(combn(grps,2))
  nspecies = nrow(dt)
  names_ = rownames(dt)
  result = data.frame(matrix(NA,nrow = nrow(com)*nspecies, ncol = 19,
                             dimnames = list(NULL,c("name","g1","g2","Avg.g1","Avg.g2","fold_change","enriched",
                                                    "all.avg","all.var","pvalue",
                                                    "count1","count2","total_count1","total_count2",
                                                    "rank1.avg", "rank2.avg","method","var1","var2"))))
  nr = 1
  for (n in 1:nspecies){
    cat("\r",n, " / ", nspecies)
    temp_dt = dt[n,]
    for(c in 1:nrow(com)){
      g1 = com[c,1]
      g2 = com[c,2]
      
      g1s = sample_map[which(sample_map[,group] == g1), ID] # group
      g2s = sample_map[which(sample_map[,group] == g2), ID]
      
      dt1 = as.matrix(temp_dt[,g1s]) # data
      dt2 = as.matrix(temp_dt[,g2s])
      
      c1 = sum(dt1 != 0, na.rm=T)  # count !0
      c2 = sum(dt2 != 0, na.rm=T)
      
      tc1 = length(dt1) # total count
      tc2 = length(dt2)
      
      m1 = mean(dt1, na.rm=T) # mean data
      m2 = mean(dt2, na.rm=T)
      
      var_1 = var(as.numeric(dt1), na.rm=T) # var data
      var_2 = var(as.numeric(dt2), na.rm=T)
      
      am = mean(c(dt1,dt2), na.rm=T) # total mean
      a_var=var(c(dt1,dt2), na.rm=T) # total var
      p = my_test(dt1,dt2)$p.value # pvalue
      fold_change = ifelse(m1>m2, m1/m2, m2/m1) # fold change
      enriched = ifelse(m1>m2, g1,g2) # enriched
      
      m = sample_map[which(sample_map[,group] %in% c(g1,g2)), ID]
      all_rank = rank(temp_dt[,m])
      
      rank1 = all_rank[colnames(dt1)] # rank
      rank2 = all_rank[colnames(dt2)]
      
      rank1.avg = mean(rank1) # mean rank
      rank2.avg = mean(rank2)
      
      result[nr,] = c(names_[n], g1, g2, m1, m2, 
                      fold_change, enriched, am, a_var,
                      p, c1, c2, tc1, tc2,rank1.avg, rank2.avg,
                      method=p.method, var1=var_1, var2=var_2)
      nr = nr+1
    }
  }
  result[,c(4:6,8:16,18,19)] = lapply(result[,c(4:6,8:16,18,19)], as.numeric)
  result
}

zy_qvalue = function(dt=NA, sample_map=NULL, pvalue_dt=NULL, group="Group", ID="Sample",p.method="wilcox.test",
                     adj.method="BH", min_count=0, min_prevalence=0, min_avg=0,min_fd=0){
  if(min_prevalence < 0 || min_prevalence > 1){
    stop("min_prevalence should be in 0~1")
  }
  result = rbind()
  if (typeof(pvalue_dt) != "NULL" && typeof(sample_map) != "NULL"){
    result.pval = pvalue_dt
  }else{
    result.pval <- as.data.frame(zy_pvalue(dt, sample_map, group=group, ID=ID, p.method=p.method))
  }
  comp <- combn(unique(sample_map[,group]),2,list)
  for ( i in 1:length(comp)){
    result.tmp1 <- result.pval %>%
      dplyr::filter(g1 == comp[[i]][1] & g2 == comp[[i]][2])
    
    result.tmp1$qvalue = NA
    
    result1 <- result.tmp1 %>%
      dplyr::filter((count1 >= min_count | count2 >= min_count)
                    & (count1/total_count1 > min_prevalence | count2/total_count2 > min_prevalence)
                    & fold_change >= min_fd & all.avg >= min_avg)
    
    result1$qvalue = p.adjust(result1$pvalue, method=adj.method)
    
    result2 <- result.tmp1 %>%
      dplyr::filter((count1 < min_count & count2 < min_count)
                    & (count1/total_count1 < min_prevalence & count2/total_count2 < min_prevalence)
                    | fold_change < min_fd | all.avg < min_avg) %>%
      rbind(result1)
    result = rbind(result,result2)
  }
  result
}

unique_prjs <- unique(sample_map$PRJ)

process_prj <- function(prj) {
  message("Processing PRJ: ", prj)
  
  sample_map_sub <- sample_map %>% filter(PRJ == prj)
  dt_sub <- dt[, colnames(dt) %in% sample_map_sub$Sample]
  
  if (ncol(dt_sub) == 0) {
    message("No samples matched for PRJ = ", prj)
    return(NULL)
  }
  
  result_sub <- zy_pvalue(dt = dt_sub, sample_map = sample_map_sub,
                          group = "Group", ID = "Sample", p.method = "wilcox.test")
  
  if (is.null(result_sub) || nrow(result_sub) == 0) {
    message("No results from zy_pvalue for PRJ = ", prj)
    return(NULL)
  }
  
  result_sub <- zy_qvalue(
    sample_map = sample_map_sub,
    pvalue_dt = result_sub,
    adj.method = "BH",
    min_count = 0,
    min_prevalence = 0,
    min_avg = 0,
    min_fd = 0 
  )
  
  result_sub$PRJ <- prj
  
  return(result_sub)
}

num_cores <- detectCores() - 1 
all_results_list <- mclapply(unique_prjs, process_prj, mc.cores = num_cores)
#save(all_results_list,file = "wilcox_all_results_list.Rdata")

all_results <- do.call(rbind, all_results_list[!sapply(all_results_list, is.null)])
missing_prjs <- setdiff(unique_prjs, unique(all_results$PRJ))
if(length(missing_prjs) > 0){
  message("Missing PRJs after processing and q-value adjustment: ", paste(missing_prjs, collapse = ", "))
}

output_dir = "./"
lapply(unique_prjs, function(prj) {
  sub_result <- all_results[all_results$PRJ == prj, , drop = FALSE]
  
  if(nrow(sub_result) > 0){
    output_file <- file.path(output_dir, paste0("wilcox_", prj, ".txt"))
    write.table(sub_result, 
                file = output_file,
                sep = "\t", row.names = FALSE, quote = FALSE)
    message("Saved: ", output_file)
  } else {
    message("No data to save for PRJ: ", prj)
  }
})
#write.table(all_results, "wilcox/merged_all_results.txt",sep = "\t", row.names = F, quote = FALSE)

#######################maaslin
dt = dt[, colSums(dt)!=0]
intersect_id = intersect(colnames(dt), sample_map$Sample)

sampf = unique(subset(sample_map, Sample %in% intersect_id, c("PRJ","Sample","Group")))
dtf = t(dt[,sampf$Sample])
rownames(dtf) = rownames(sampf)

library(Maaslin2)
main_output_dir <- "./"

dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)

projects <- unique(sampf$PRJ)

run_maaslin2_project <- function(proj) {
  message("Processing project: ", proj)
  proj_output_dir <- file.path(main_output_dir, proj)
  proj_output_dir <- normalizePath(proj_output_dir, winslash = "/")
  message("Output directory: ", proj_output_dir)
  
  dir.create(proj_output_dir, showWarnings = TRUE, recursive = TRUE)
  
  proj_samples <- sampf[sampf$PRJ == proj, ]
  proj_dtf <- dtf[rownames(proj_samples), , drop = FALSE]
  proj_dtf <- as.data.frame(proj_dtf) 
  
  rownames(proj_samples) <- rownames(proj_dtf)
  
  if (ncol(proj_dtf) == 0 || nrow(proj_dtf) == 0) {
    stop("Input data is empty for project: ", proj)
  }
  
  fit_data <- Maaslin2(
    input_data = proj_dtf,
    input_metadata = proj_samples,
    output = proj_output_dir,
    min_prevalence = 0,
    random_effects = c(),
    fixed_effects = c("Group"),
    normalization = "NONE",
    reference = c("Group", "Donor"),
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  save_path <- file.path(proj_output_dir, paste0(proj, "_results.RData"))
  message("Saving results to: ", save_path)
  save(fit_data, file = save_path)
  
  return(fit_data)
}

num_cores <- detectCores() - 1 
results <- mclapply(projects, run_maaslin2_project, mc.cores = num_cores)

projects <- unique(sampf$PRJ)
main_output_dir <- "./"
all_results_list <- list()
for (proj in projects) {
  file_path <- file.path(main_output_dir, proj, "all_results.tsv")
  
  if (file.exists(file_path)) {
    df <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df$Project <- proj 
    all_results_list[[proj]] <- df
  } else {
    message("Warning: File not found for project ", proj)
  }
}
combined_df <- do.call(rbind, all_results_list)
combined_df <- combined_df %>%
  filter((value == "Donor" ) | (value == "pre-FMT" )) 
write.table(
  combined_df,
  file = file.path(main_output_dir, "combined_all_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

##################丰度气泡图
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
library(dplyr)
library(stringr)
filtered_data <- read.table("meta.analysis_0.0002_0.25_filtered.tsv", header = T,  sep = "\t")

unique_features <- filtered_data %>%
  distinct(feature) %>%
  pull(feature)

feature <- as.data.frame(unique_features)
names(feature) <- "feature"

tax <- read.table("taxonomy.txt", header = T,  sep = "\t")
merge_dt <- merge(tax, feature, by.x = "name",by.y="feature")

merge_dt <- merge_dt$species
merge_dt <- as.data.frame(merge_dt)
names(merge_dt) <-"species"

filtered_data2 <- read.table("../05.keystone/meta.analysis/meta.analysis_0.0002_0.25.tsv", header = T,  sep = "\t")

filtered_data3 <- filtered_data2 %>%
  filter(I2 <= 50, qval <= 0.05) %>%
  group_by(feature) %>%
  filter((ci_lb > 0 & ci_ub > 0) | (ci_lb < 0 & ci_ub < 0))

unique_features <- filtered_data2 %>%
  distinct(feature) %>%
  pull(feature)

feature2 <- as.data.frame(unique_features)
names(feature2) <- "feature"

###########
tax <- read.table("../02.taxon/taxonomy.txt", header = T,  sep = "\t")
all_data <- merge(tax, filtered_data2, by.x = "name",by.y="feature")

qval <- all_data %>% dplyr::select(species,I2, qval)
qval <- qval %>%
  mutate(across(-(1), as.numeric))

qval <- qval %>%
  mutate(group = "Others") %>% 
  mutate(group = ifelse(species %in% merge_dt$species, "Keystone", group)) 

qval$qval <--log10(qval$qval)

unique_qval <- qval %>%
  distinct(species,.keep_all = T) 

library(ggplot2)
library(scales)  

ggplot(qval, aes(x = qval, y = I2)) +
  geom_point(aes(color = group), size = 2) +
  scale_color_manual(values = c("Keystone" = "red", "Others" = "gray")) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey") +  
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20)) +  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(0.25, "cm"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = rel(1.5)),  
    axis.title = element_text(size = rel(1.25)),           
    axis.text = element_text(size = rel(1.1)),       
    legend.title = element_text(size = rel(1.2)),    
    legend.text = element_text(size = rel(1.1)),   
    aspect.ratio = 0.7
  ) +
  labs(title = "",
       x = "-log10(qvalue)",
       y = "I2")

##################富集方向热图
library(readxl)
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
library(dplyr)
library(stringr)
library(scales)

taxo_sample_to_grps <- function(sample_map=NA, group=NA,ID=NA){
  grps = unique(sample_map[,group])
  grps_list = list()
  for(grp in grps){
    gg <- sample_map[which(sample_map[,group] == grp), ID]
    grps_list[[grp]] = gg
  }
  return(grps_list)
}

grc_profile_to_test = function(Dt=NA, sample_map=NA, group=NA, ID=NA){
  grps <- taxo_sample_to_grps(sample_map = sample_map, group = group, ID = ID)
  matched_columns <- intersect(colnames(Dt), sample_map[, ID])
  if (length(matched_columns) == 0) {
    stop("No matching columns found between species data and sample map.")
  }
  Dt_filtered <- Dt[, matched_columns]
  
  com = t(combn(names(grps), 2))
  nspecies = nrow(Dt_filtered)
  names = rownames(Dt_filtered)
  result = matrix(NA, nrow = nrow(com)*nspecies, ncol = 8,
                  dimnames = list(NULL, c("name","g1","g2","m1","m2","enriched","fold_change","pvalue")))
  nr = 1
  
  for (n in 1:nspecies){
    temp_dt = Dt_filtered[n,]
    for(c in 1:nrow(com)){
      g1 = com[c,1]
      g2 = com[c,2]
      
      if (!all(grps[[g1]] %in% colnames(temp_dt))) {
        stop(paste("Group", g1, "contains samples not found in the data."))
      }
      if (!all(grps[[g2]] %in% colnames(temp_dt))) {
        stop(paste("Group", g2, "contains samples not found in the data."))
      }
      
      dt1 = as.matrix(temp_dt[,grps[[g1]]])
      dt2 = as.matrix(temp_dt[,grps[[g2]]])
      m1 = mean(dt1)
      m2 = mean(dt2)
      enrich = ifelse(m1 > m2, g1, ifelse(m1 < m2, g2, NA))
      fold = max(m1/m2, m2/m1)
      p = wilcox.test(dt1, dt2, exact = FALSE, correct = TRUE)$p.value#p = wilcox.test(dt1, dt2)$p.value#这是最初的版本
      result[nr,] = c(names[n], g1, g2, m1, m2, enrich, fold, p)
      nr = nr + 1
    }
  }
  result = as.data.frame(result)
  
  for(x in c("m1","m2","fold_change","pvalue")){
    result[,x] = as.numeric(result[,x])
  }
  
  result$qvalue = NA
  for (c in 1:nrow(com)){
    ff = which(result$g1 == com[c,1] & result$g2 == com[c,2] & !is.na(result$fold_change))
    result$qvalue[ff] = p.adjust(result$pvalue[ff], method = 'BH')
  }
  
  return(result)
}


grc_test_to_volcano = function(Test,g1='CON',g2='Case',p_or_q='pvalue',p_or_q_cutoff=0.05,fold_change_cutoff=0){
  
  dat = Test[which(Test$g1==g1 & Test$g2==g2 & !is.na(Test$pvalue)),]
  
  dat$point_size = dat$m1
  dat$point_size[dat$enriched==g2] = dat$m2[dat$enriched==g2]
  
  dat$point_color = dat$enriched
  dat$point_color[dat$fold_change < abs(fold_change_cutoff) | dat[,p_or_q] > p_or_q_cutoff] = 'PASS'
  
  dat$point_axis_x = dat$m1/dat$m2
  dat$point_axis_x = log2(dat$point_axis_x)
  dat$point_axis_x[dat$point_axis_x > 10] = 10
  dat$point_axis_x[dat$point_axis_x < -10] = -10
  
  dat$point_axis_y = -log10(dat[,p_or_q])
  dat$point_axis_y[dat$point_axis_y>10] = 10
  
  dat$point_color = factor(dat$point_color,c(g1,g2,'PASS'))
  
  fold_change_cutoff = ifelse(fold_change_cutoff==0,0,c(-log2(fold_change_cutoff),log2(fold_change_cutoff)))
  
  Fig = ggplot(dat,aes(x=point_axis_x,y=point_axis_y,color=point_color))+ 
    geom_point(aes(size=point_size),shape=1)+
    geom_vline(xintercept=unique(c(as.numeric(fold_change_cutoff),-as.numeric(fold_change_cutoff))),linetype="dashed")+
    geom_hline(yintercept=-log10(p_or_q_cutoff),linetype="dashed")+
    theme_bw()+
    theme(
    )+
    scale_color_manual(values = c('red','blue','grey'))+
    xlab('log2(Fold Change)')+
    ylab(paste('-log10(',p_or_q,')',sep=''))
  
  return(Fig)
}

Dt <- read.table("../05.keystone/meta.analysis_0.0002_0.25_filtered.tsv", header = T,  sep = "\t")

unique_features <- Dt %>%
  distinct(feature) %>%
  pull(feature)

feature <- as.data.frame(unique_features)
names(feature) <- "feature"

dt <- read.table("../01.abun/re_filter_Donor_pre-FMT.txt", header = T,  sep = "\t")
dt <- rownames_to_column(dt,var = "feature")

sig_dt <- merge(feature,dt,by = "feature")
sig_dt <- column_to_rownames(sig_dt,var="feature")
#write.table(sig_dt,file = "significance_abundance.txt", sep = "\t", row.names = TRUE, quote = FALSE)

sample_map = read.table("Donor_pre-FMT_sample.info.txt", sep="\t", header=T, check.names=F)

batch_process_by_PRJ <- function(abundance_data, sample_map) {
  colnames(abundance_data) <- as.character(colnames(abundance_data))
  prjs <- unique(sample_map$PRJ)
  all_results <- list()
  for (prj in prjs) {
    message("Processing project: ", prj)
    current_sample_map <- sample_map %>%
      filter(PRJ == prj,
             Group %in% c("Donor", "pre-FMT"))
    if (!all(c("Donor", "pre-FMT") %in% unique(current_sample_map$Group))) {
      message("Skipping PRJ: ", prj, " due to missing groups.")
      next
    }
    result <- grc_profile_to_test(
      Dt = abundance_data,
      sample_map = current_sample_map,
      group = "Group",
      ID = "Sample"
    )
    result$PRJ <- prj
    all_results[[prj]] <- result
  }
  final_result <- do.call(rbind, all_results)
  return(final_result)
}

final_result <- batch_process_by_PRJ(abundance_data = sig_dt,
                                     sample_map = sample_map)
#write.table(final_result,file = "fold_change.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

data <- final_result
data$log10_fold_change <- log10(data$fold_change)

data2 <- data %>%
  mutate(
    fold_change_display = case_when(
      log10_fold_change == Inf & m1 == 0 & m2 == 0 ~ 0,
      log10_fold_change == Inf & m1 == 0 ~ m2,
      log10_fold_change == Inf & m2 == 0 ~ m1,
      log10_fold_change == Inf ~ Inf,
      TRUE ~ log10_fold_change
    )
  )

data3 <- data2 %>%
  mutate(
    log10_fold_change_display = case_when(
      enriched == "Donor" ~ fold_change_display,
      enriched == "pre-FMT" & fold_change_display != 0 ~ -fold_change_display,
      TRUE ~ 0
    )
  )

data4 <- data3 %>%
  dplyr::select(name, log10_fold_change_display,PRJ) %>%
  mutate(name = as.character(name)) 

library(data.table)
data_long <- data4 %>%
  melt(id.vars = c("name", "PRJ")) 

data_long_sorted <- data_long %>%
  arrange(value)

order_data <- unique(as.data.frame(data_long_sorted$name))
names(order_data) <- "name"
#write.table(order_data,file = "species_order.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

data_long_sorted$name <- factor(data_long_sorted$name, levels = unique(order_data$name))

# 计算用于颜色渐变的数据范围
max_abs_value <- max(abs(data_long_sorted$value))

p <- ggplot(data_long_sorted, aes(x = name, y = variable, fill = value)) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(~ PRJ, scales = "free_x", ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Heatmap of log10 Fold Change and Significance (q-value)") +
  scale_fill_gradientn(
    name = "Fold Change",
    colours = c(
      "#DC143C",       # 极端负值（可选）
      "#FFE6E6",       # 粉色（负值接近0）
      "white",         # 中心点（0）
      "#E6F5FF",       # 淡蓝（正值接近0）
      "#0077B5"        # 极端正值（可选）
    ),
    values = rescale(c(-max_abs_value, -1e-08,0, 1e-08, max_abs_value)),  # 控制各颜色对应的位置
    limits = c(-max_abs_value, max_abs_value),
    guide = guide_colourbar(tickreverse = TRUE)
  )

print(p)

####################外部验证数据metap分析
#两个数据分别进行maaslin分析，获得p值进行metap分析
pvalue <- read.table("pvalue.txt", sep="\t", header=T, check.names=F)
pvalue[is.na(pvalue)] <- 0

chisq_values <- numeric(nrow(pvalue))
df_values <- integer(nrow(pvalue))
p_combined <- numeric(nrow(pvalue))

for (i in 1:nrow(pvalue)) {
  p1 <- pvalue$Jieun_2023[i]
  p2 <- pvalue$Bushman_2020[i]
  
  p_vec <- c(p1, p2)
  # 执行 Fisher's sumlog 方法
  result <- try(sumlog(p_vec), silent = TRUE)
  
  if (inherits(result, "try-error")) {
    chisq_values[i] <- NA
    p_combined[i] <- NA
  } else {
    chisq_values[i] <- result$chisq
    p_combined[i] <- result$p
  }
  df_values[i] <- 4 
}

p_values_meta <- data.frame(
  pvalue,
  meta_chisq = chisq_values,
  meta_df = df_values,
  meta_p = p_combined
)

pvalue_adj <- p_values_meta %>%
  mutate(adj_meta_p = p.adjust(meta_p, method = "fdr"))
#write.table(pvalue_adj, "metap.txt",sep = "\t", quote = FALSE, row.names = FALSE)

##########data.txt为不同样本的数据
data <- read.table("data.txt", sep="\t", header=T, check.names=F)
meta_metap <- data[,c(1,12:15)]

Inconsistent <- meta_metap[meta_metap$enriched == "-",]
Inconsistent_sig <- Inconsistent[Inconsistent$significance=="Significant",]
nrow(Inconsistent_sig)

Control <- meta_metap[meta_metap$enriched=="Control",]
Control_sig <- Control[Control$significance=="Significant",]
nrow(Control_sig)

Case <- meta_metap[meta_metap$enriched=="Case",]
Case_sig <- Case[Case$significance=="Significant",]
nrow(Case_sig)

bar_data <- data.frame(
  Group = c("Control", "Case", "Inconsistent"),
  Significant = c(117, 104, 9),
  Not_Significant = c(129 - 117, 112 - 104, 42 - 9)
)

data_long <- reshape2::melt(bar_data, id.vars = "Group", variable.name = "Status", value.name = "Count")

data_long$Status <- factor(data_long$Status, levels = c("Not_Significant","Significant"))
data_long$Group <- factor(data_long$Group, levels = c("Control","Case","Inconsistent"))

colors <- c("#0077B5", "#ADD8E6", "#DC143C", "#ffcccc", "grey30", "grey")
names(colors) <- c("Control.Significant", "Control.Not_Significant",
                   "Case.Significant", "Case.Not_Significant",
                   "Inconsistent.Significant", "Inconsistent.Not_Significant")

data_long$fill_var <- interaction(data_long$Group, data_long$Status)

ggplot(data_long, aes(x = Group, y = Count, fill = fill_var)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = colors, 
                    labels = c("Not_Significant", "Not_Significant", 
                               "Not_Significant", "Significant", 
                               "Significant", "Significant")) +
  labs(x = "Groups", y = "Counts", fill = "Status") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    text = element_text(size = 17),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    plot.title = element_text(size = 16, face = "bold"),
    aspect.ratio = 2,
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, 150))

################丰度气泡图
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

dt <- read.table("/share/data1/limin/FMT/rCDI/Analysis/00.data/CDI/abundance_CDI.txt", header = T,  sep = "\t")
dt <- column_to_rownames(dt,"Key")
keystone <- read.table("/share/data1/limin/FMT/rCDI/Analysis/05.keystone/meta.analysis_0.0002_0.25_unique.tsv", header = T,  sep = "\t")
keystone <- keystone$feature
dt_filtered <- dt[row.names(dt) %in% keystone, ]

sample_map <-  read.table("/share/data1/limin/FMT/rCDI/Analysis/00.data/metadata_CDI.txt", header = T,  sep = "\t")
sample_map <- sample_map[,c(2,4,6)]
names(sample_map) <- c("Project","Sample","Group")
sample_map1 <- sample_map[sample_map$Project == "Jieun_2023",]
sample_map2 <- sample_map[sample_map$Project == "Bushman_2020",]

result <- align_dt_sample(dt_filtered, sample_map2, ID = "Sample")
aligned_dt <- result$dt

aligned_dt$average_re_abun <- rowMeans(aligned_dt, na.rm = TRUE)
dt_filtered1 <- rownames_to_column(aligned_dt,var="feature")
dt_filtered2 <- rownames_to_column(aligned_dt,var="feature")

meta_metap1 <- read.table("Jieun_2023_maaslin.txt", header = T,  sep = "\t")
meta_metap2 <- read.table("Bushman_2020_maaslin.txt", header = T,  sep = "\t")

keystone <- as.data.frame(meta_metap1[,1])
names(keystone) <- "feature"
meta_metap2 <- left_join(keystone,meta_metap2,by="feature")
meta_metap2$coef[is.na(meta_metap2$coef)] <- "0"
meta_metap2$qval[is.na(meta_metap2$qval)] <- "1"

#############
maaslin <- meta_metap1 %>%
  mutate(enriched = case_when(
    coef < 0 ~ "Donor",
    coef > 0 ~ "pre-FMT",
    TRUE ~ "-"
  ))

masslin2 <- merge(maaslin,dt_filtered1[,c(1,89)],by="feature")

masslin3 <- masslin2 %>%
  select(feature,coef,qval,average_re_abun)%>%
  left_join(meta_metap[, c("feature", "enriched","significance")], by = "feature")

masslin3$qval <- as.numeric(masslin3$qval)
masslin3$log_qvalue <- -log10(masslin3$qval)

df1 <- masslin3 %>%
  mutate(
    border_color = case_when(
      enriched == "Control" & significance == "Significant" ~ "#0077B5",  
      enriched == "Control" & significance == "Not significant" ~ "#ADD8E6",                  
      enriched == "Case" & significance == "Significant" ~ "#DC143C",     
      enriched == "Case"  & significance == "Not significant" ~"#ffcccc",                      
      enriched == "-" & significance == "Significant"  ~ "grey30",     
      enriched == "-" & significance == "Not significant" ~ "grey",
      TRUE ~ NA_character_
    ),
    
    fill_color = NA_character_  
  )

#############
maaslin <- meta_metap2 %>%
  mutate(enriched = case_when(
    coef < 0 ~ "Donor",
    coef > 0 ~ "pre-FMT",
    TRUE ~ "-"
  ))

masslin2 <- left_join(maaslin,dt_filtered2[,c(1,20)],by="feature")
masslin2$average_re_abun[is.na(masslin2$average_re_abun)] <- "0"

masslin3 <- masslin2 %>%
  select(feature,coef,qval,average_re_abun)%>%
  left_join(meta_metap[, c("feature", "enriched","significance")], by = "feature")

masslin3$qval <- as.numeric(masslin3$qval)
masslin3$log_qvalue <- -log10(masslin3$qval)
masslin3$average_re_abun <- as.numeric(masslin3$average_re_abun)
masslin3$coef <- as.numeric(masslin3$coef)

df2 <- masslin3 %>%
  mutate(
    border_color = case_when(
      enriched == "Control" & significance == "Significant" ~ "#0077B5",  
      enriched == "Control" & significance == "Not significant" ~ "#ADD8E6",                  
      enriched == "Case" & significance == "Significant" ~ "#DC143C",     
      enriched == "Case"  & significance == "Not significant" ~"#ffcccc",                      
      enriched == "-" & significance == "Significant"  ~ "grey30",     
      enriched == "-" & significance == "Not significant" ~ "grey",
      TRUE ~ NA_character_
    ),
    
    fill_color = NA_character_ 
  )

#################
global_min <- min(min(df1$average_re_abun, na.rm = TRUE), min(df2$average_re_abun, na.rm = TRUE))
global_max <- max(max(df1$average_re_abun, na.rm = TRUE), max(df2$average_re_abun, na.rm = TRUE))

p1 <- ggplot(df1, aes(x = coef, y = log_qvalue)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_point(aes(size = average_re_abun, fill = fill_color, color = border_color), shape = 21, stroke = 1.2) +  
  scale_fill_identity() +  
  scale_color_identity() + 
  scale_size_continuous(range = c(1, 5), limits = c(global_min, global_max), name = "Average\nAbundance") +
  labs(
    x = "Coefficient",
    y = "-log₁₀(q-value)",
    title = "Bubble Plot of Coef vs. Significance (Dataset 1)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.background = element_rect(color = "white") 
  ) +
  
  ylim(0, max(df1$log_qvalue, na.rm = TRUE) * 1.1)

p2 <- ggplot(df2, aes(x = coef, y = log_qvalue)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_point(aes(size = average_re_abun, fill = fill_color, color = border_color), shape = 21, stroke = 1.2) +  
  scale_fill_identity() +  
  scale_color_identity() + 
  scale_size_continuous(range = c(1, 5), limits = c(global_min, global_max), name = "Average\nAbundance") +
  labs(
    x = "Coefficient",
    y = "-log₁₀(q-value)",
    title = "Bubble Plot of Coef vs. Significance (Dataset 2)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 17),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    plot.margin = margin(10, 10, 10, 10),
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.background = element_rect(color = "white") 
  ) +
  
  ylim(0, max(df2$log_qvalue, na.rm = TRUE) * 1.1)

library(patchwork)
combined_plot <- (p1 / p2 ) + plot_layout(heights = c(1, 1))
print(combined_plot)