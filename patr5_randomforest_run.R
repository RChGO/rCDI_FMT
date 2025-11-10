############重要性排序

#!/bin/bash

ABD_DIR="../01.abundance"
GROUP_DIR="../02.group"
OUTPUT_DIR="./"

mkdir -p "$OUTPUT_DIR"

for abd_file in "$ABD_DIR"/*_abundance.txt; do
base=$(basename "$abd_file" _abundance.txt)

group_file="$GROUP_DIR/${base}_group_modified.txt"

if [ ! -f "$group_file" ]; then
echo "⚠️ 缺少对应的 group 文件: $group_file"
continue
fi

output_prefix="$OUTPUT_DIR/${base}"

echo "Running: rf.py importance -i $abd_file -g $group_file -o $output_prefix --threads 30 -K 3 -R 50 --stratified --fast"
rf.py importance -i "$abd_file" -g "$group_file" -o "$output_prefix" --threads 30 -K 3 -R 50 --stratified --fast
done

echo "✅ 所有任务完成。"

#############AUC
abundance_dir="./"
group_dir="../02.group"
output_dir="auc"
script="rf.py"
seed_list="seed.list"

topNs=( 1 $(seq 5 5 100) )

mkdir -p "$output_dir"

projects=("Fungi" "Virus" )

for project in "${projects[@]}"; do
echo "Processing $project..."

abundance_file="${abundance_dir}/${project}.sorted"
group_file="${group_dir}/${project}_group_modified.txt"
out_subdir="${output_dir}/${project}"

if [[ ! -f "$abundance_file" ]]; then
echo "Error: Abundance file not found: $abundance_file"
continue
fi
if [[ ! -f "$group_file" ]]; then
echo "Error: Group file not found: $group_file"
continue
fi

mkdir -p "$out_subdir"

parallel -j 20 \
python "$script" KF \
-i "$abundance_file" \
-g "$group_file" \
--threads 1 \
-o "$out_subdir/predict_top{1}_seed{2}.tsv" \
-K 3 \
--topN {1} \
--seed {2} \
>> "$out_subdir/top{1}_seed{2}.log" 2>&1 \
::: $(printf "%s " "${topNs[@]}") \
:::: "$seed_list"
done


#############计算AUC
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

my_join <- function(dt, cols) {
  dt$zy_tmp_group <- do.call(paste, c(dt[, ..cols, with = FALSE], sep = "_"))
  return(dt)
}

my_split <- function(df, tmp_col, old_cols) {
  df[, (old_cols) := tstrsplit(zy_tmp_group, "_", fixed = FALSE)]
  df[, zy_tmp_group := NULL]
  return(df)
}

analyze_file <- function(file_path, output_path) {
  library(pROC)
  
  dt <- read.table(file_path, sep = '\t', header = TRUE, fill = TRUE)
  dt$X <- paste(dt$seed, dt$nspecies, sep = "_")
  
  otu <- dt
  result <- calc_auc(otu, pred = "Multiple_1", true = "Actual", group = "X")
  auc_table <- result$table
  
  rownames(auc_table) <- auc_table$X
  auc_table$X <- sub("^(.*?)_", "", rownames(auc_table))
  
  write.table(auc_table, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Saved: ", output_path)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript run_auc_analysis.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

analyze_file(input_file, output_file)
