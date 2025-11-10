#!/share/data1/software/miniconda3/envs/bio38/bin/python
# -*- encoding: utf-8 -*-
##########################################################
# Creater       :  夜下凝月
# Created  date :  2024-07-23, 17:03
# Modiffed date :  2025-02-12, 16:20
##########################################################

'''
Tips:
    如果分组中，有一列全是数字，脚本会将这一列字数字前面的0替换掉，这个需要注意
    
    v20250212:
        新增：
            1、如果两个分组中，样本数量差距大于0.8， 则会在抽样是采用欠采样的方式
            2、计算重要性排序时，因为permutation_importance太慢了，所以可以启用参数`--fast`，使用构建模型时的特征重要性(速度超快)。但相应的也会高估部分feature的重要性
'''

import argparse
import sys, os

def get_args():
  
  parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

## 公共参数 
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument("-i", metavar="", required=True,  type=str, help="input Matrix")
parent_parser.add_argument("-g", metavar="", required=True,  type=str, help="group file")
parent_parser.add_argument("-o", metavar="", required=True,  type=str, help="Output")
# parent_parser.add_argument("-f", metavar="", required=False, type=str, default="r", choices=['r','l', 'randomforest', 'lasso'], help="预测方式,(l[asso],r[andomforest])")
# parent_parser.add_argument("-c", metavar="", required=False, type=tuple, help="控制变量,分组中的列名. 如('Age', 'Sex','BMI')")
parent_parser.add_argument("-gid", metavar="", required=False, type=str, default="Group", help="Group ID. [default: Group]")
parent_parser.add_argument("-sid", metavar="", required=False, type=str, default="Sample", help="Sample ID. [default: Sample]")
parent_parser.add_argument("-s", metavar="", required=False, type=int, default=1, choices=[0,1,2,3], help="Split Str fo matrix and group.[default: 1] \033[31mThey must have similary split str\033[0m\n0 -> \\s+\n1 -> \\t\n2 -> |\n3 -> ,")
parent_parser.add_argument("--topN", metavar="", type=int, default=0, help="使用前N个特征[0:使用所有的; 如果使用的大于所有的，那也使用所有的特征. default: 0]")
parent_parser.add_argument("--seed", metavar="", type=int, default=2024, help="随机种子. [default: 2024]")
parent_parser.add_argument("--nest", metavar="", type=int, default=1000, help="n_estimators,对应R语言中的ntree.[default: 1000]")
parent_parser.add_argument("--threads", metavar="", type=int, default=10, help="线程数.[default: 10]")

## 子命令参数
subparsers = parser.add_subparsers(dest='command', required=True)

### 1、获取重要性排序，以及输出排序后的matrix
import_desc = '''输出特征重要性，以及排序后的丰度表. 采用多折交叉检验的方式计算重要性, 在计算特征的重要性时，会进行置换检验，这边使用-R实现 '''
parser_1 = subparsers.add_parser('importance', parents=[parent_parser], help= import_desc, formatter_class=argparse.RawTextHelpFormatter)
parser_1.add_argument("-R", default=10, type=int, help="重复次数，.[default: 10]")
parser_1.add_argument("-K", default=10, type=int, help="利用多折检验来计算重要性")
parser_1.add_argument("--stratified", action="store_true", help="是否分层划分数据，每一折的数据标签保持比例和总的相同.")
parser_1.add_argument("--fast", action="store_true", help="因为默认的计算重要性需要置换，所以太慢了，如果特征数太多了，就启用这个参数，使用模型自带的重要性")

### 2、K-Fold
parser_2 = subparsers.add_parser("KF", parents=[parent_parser], help="K-Fold检验", formatter_class=argparse.RawTextHelpFormatter)
parser_2.add_argument("-K", default=10, type=int, help="K折检验.[default: 10]")
parser_2.add_argument("--stratified", action="store_true", help="是否分层划分数据，每一折的数据标签保持比例和总的相同.")
parser_2.add_argument("--down", dest="downrate", default=0.5, help="如果分组比例低于0.5，就会使用用欠采样的方式，也就是在较多的分组中，抽取同等数量的样本.")

### 3、留一法LOOCV
parser_3 = subparsers.add_parser("LOOCV", parents=[parent_parser], help="留一法", formatter_class=argparse.RawTextHelpFormatter)

### 4、按比例分割后预测
parser_4 = subparsers.add_parser("easy", parents=[parent_parser], help="按照比例分割数据集，一部分进行建模型，一部分作为测试集", formatter_class=argparse.RawTextHelpFormatter)
parser_4.add_argument("-S", default=0.3, type=float, help="测试集大小0.01 ~ 0.99。[default: 0.3]")
parser_4.add_argument("-om", type=str, help="输出建立的随机森林模型")

### 5、输入模型，然后预测别的数据
#### 按理说，作为预测，应当只有一个预测值，不用输入-G, -S, -g，后面有时间改一下
parser_5 = subparsers.add_parser("predict", parents=[parent_parser], help="输入已有的模型，预测新的数据", formatter_class=argparse.RawTextHelpFormatter)
parser_5.add_argument("-M", required=True, type=str, help="已有的模型")
parser_5.add_argument("--fillna", action="store_true", help="自动填补缺失的feature为0")

### 6、将输入的数据作为整体建立模型
parser_6 = subparsers.add_parser("build", parents=[parent_parser], help="根据给定的表建立随机森林模型", formatter_class=argparse.RawTextHelpFormatter)

### 7、通过调整参数来获取最优解
parser_7 = subparsers.add_parser("best", parents=[parent_parser], help="给定某些参数范围，脚本自行获取最优解", formatter_class=argparse.RawTextHelpFormatter)
parser_7.add_argument("-est", default=(0,9999), type=tuple, help="将生成随机")

args = parser.parse_args()
return(args)


def rf_predict(model, X, y, threads=10, n_est=999, seed=42, fillna=False):
  
  mf = model.myfeatures
intersect_id = mf.intersection(X.columns)

if (len(mf) != len(intersect_id)):
  miss_features = set(mf).difference(intersect_id)
if fillna:
  ##  填充缺失的features
  logging.warning("输入的特诊缺少了模型中的feature，将自动补充为0")
logging.warning(miss_features)
miss_dt = pd.DataFrame(0, index=X.index, columns=list(miss_features))
X = pd.concat([X, miss_dt], axis=1)
else:
  logging.error("输入数据中缺少模型所需要的features")
logging.error(miss_features)
logging.error("Failure!!!")
exit(127)

X = X.loc[:, mf]
res = pd.DataFrame(model.predict_proba(X), columns=model.classes_, index=X.index)
res.loc[:,'Group'] = y.values
# auc = roc_auc_score(y.values, res.iloc[:,1])
# logging.info(f"AUC:\t{auc}")
return(res)


def easy_predict():
  logging.error("deprecated")
exit(127)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train) # 训练模型
res = pd.DataFrame(model.predict_proba(X_test), columns=model.classes_, index=X_test.index)
res.loc[:,'Group'] = y_test.values
return(res)


def KF_scores(X, y):
  model = RandomForestClassifier(n_estimators=999) # 初始化模型
kf = KFold(n_splits=10, shuffle=True, random_state=1) # 设置K-fold交叉验证
scores = cross_val_score(model, X, y, cv=kf) # 这个得到的是每一折的准确率
return(scores)



def KF_predctions(X, y, Sample='Sample', k=10, n_est=1000, n_job=10, seed=2024, Stratified=False, downrate=0.5, **kwargs):
  '''
        X 2d-matrix
        y np.array
    '''
# 欠采样函数
def undersample(X, y):
  # 统计每个类别的样本数
  unique, counts = np.unique(y, return_counts=True)
min_class = unique[np.argmin(counts)]  # 找到最少类别
min_class_count = np.min(counts)

# 找到多数类的索引
majority_class = unique[np.argmax(counts)]
majority_class_indices = np.where(y == majority_class)[0]

# 随机选择与少数类相同数量的多数类样本
np.random.seed(42)  # 为了可复现性
random_majority_class_indices = np.random.choice(majority_class_indices, size=min_class_count, replace=False)

# 组合少数类和随机选择的多数类样本
minority_class_indices = np.where(y == min_class)[0]
selected_indices = np.concatenate([minority_class_indices, random_majority_class_indices])

return X.iloc[selected_indices,:], y[selected_indices]

# 欠采样
a,b = y.value_counts()
if min(a/b, b/a) < downrate:
  logging.info(f"Sampling: RandomUndersampler. [{downrate}]")
X, y = undersample(X, y)

# 设置K-fold交叉验证
if Stratified:
  # 分层划分，就是保持每个类别中，标签的比例和原始的一样
  logging.info(f"Method: {k}-Fold (Stratified)")
kf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
kf_sp = kf.split(X, y)
else:
  logging.info(f"Method: {k}-Fold (Random)")
# 随机划分
kf = KFold(n_splits=k, shuffle=True, random_state=seed)
kf_sp = kf.split(X)

logging.info(f"Build model")
model = RandomForestClassifier(n_estimators=n_est, n_jobs=n_job)
model.fit(X,y) # 只有建模以后，才能获取model.classes_

names = X.index
n = 1
sample_split = []
for _, test_index in kf_sp:
  tmp_split = pd.DataFrame({Sample: names[test_index], "split": n})
sample_split.append(tmp_split)
n += 1
X_split = pd.concat(sample_split, axis=0)


logging.info("Running... ")
predictions = cross_val_predict(model, X, y, cv=kf, method="predict_proba")

res = pd.concat([pd.DataFrame(predictions,columns=model.classes_), pd.DataFrame(y.values,columns=['Actual'])], axis=1)
res.index = X.index
auc = roc_auc_score(res.loc[:,'Actual'], res.iloc[:,1])
topn = kwargs.get('topn')
logging.info(f"auc-seed-nspecies:\t{auc}\t{seed}\t{topn}")

res = pd.merge(X_split, res, left_on=Sample, right_index=True) # 将预测结果和数据的分割情况分开

return(model, res)


def get_importance_rf_perm(X, y,repeat=10, seed=2024, nest=1000, threads=-1, k=10, Stratified=False):
  '''
        这边计算重要性时，可以使用多折检验的方式
        我估计还有其他的方式，可以看看，再重新想想整个脚本怎么构建
            是利用类函数呢，还是传递函数名的方式来构建
    '''
logging.info(f"Get features' RandomFoest importance.")
logging.info(f"\t\\_Repeat:\t{repeat}")
logging.info(f"\t\\_K-Fold:\t{k}")
logging.info(f"\t\\_random_state:\t{seed}")
logging.info(f"\t\\_n_estimators:\t{nest}")
logging.info(f"\t\\_Stratified:\t{Stratified}")
logging.info(f"\t\\_Get features' importance by `permutation_importance`")
logging.warning(f"\033[31m不建议使用这个方式排序，推荐使用参数--fast\033[0m")
print(f"\033[31m不建议使用这个方式排序，推荐使用参数--fast\033[0m")

if Stratified:
  kf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
kf_sp = kf.split(X,y)
else:
  kf = KFold(n_splits=k, shuffle=True, random_state=seed)
kf_sp = kf.split(X)


print("Calc importance")
all_importances = pd.DataFrame(index=X.columns)
for index, (train_index, test_index) in enumerate(kf_sp):
  start_time = time.time()
print(f"Part: [ {index+1} / {k} ]", end="")
X_train, X_test = X.iloc[train_index], X.iloc[test_index]
y_train, y_test = y[train_index], y[test_index]

model = RandomForestClassifier(n_estimators=nest, oob_score=False, n_jobs=threads, random_state=seed)
model.fit(X_train, y_train)

results = permutation_importance(model, X_test, y_test, n_repeats=repeat, random_state=seed, n_jobs=threads)

fold_importances = pd.DataFrame(results.importances, index=X.columns)
all_importances = pd.concat([all_importances, fold_importances], axis=1)
print(f"\t{time.time()-start_time} Sec.")

mean_importances = all_importances.mean(axis=1).sort_values(ascending=False)
features_ord = mean_importances.index
all_importances = all_importances.reindex(index=features_ord)
return(all_importances)


def get_importance_rf_fast(X, y,repeat=10, seed=2024, nest=1000, threads=-1, k=10, Stratified=False):
  '''
        这边计算重要性时，可以使用多折检验的方式
        我估计还有其他的方式，可以看看，再重新想想整个脚本怎么构建
            是利用类函数呢，还是传递函数名的方式来构建
    '''
logging.info(f"Get features' RandomFoest importance.")
logging.info(f"\t\\_Repeat:\t{repeat}")
logging.info(f"\t\\_K-Fold:\t{k}")
logging.info(f"\t\\_random_state:\t{seed}")
logging.info(f"\t\\_n_estimators:\t{nest}")
logging.info(f"\t\\_Stratified:\t{Stratified}")
logging.info(f"\t\\_Get features' importance by `model.feature_importances_`")


## 生成n次随机种子
random.seed(seed)
seed = random.randint(0, 1000000)
seeds = [random.randint(0, 1000000) for x in range(0,repeat)]

all_importances = list()
print("Calc importance")
for i in range(0, repeat):
  seed = seeds[i]
if Stratified:
  kf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
kf_sp = kf.split(X,y)
else:
  kf = KFold(n_splits=k, shuffle=True, random_state=seed)
kf_sp = kf.split(X)

for index, (train_index, test_index) in enumerate(kf_sp):
  start_time = time.time()
print(f"Repeat: [ {i+1}/{repeat} ] - Part: [ {index+1} / {k} ]", end="")
X_train, X_test = X.iloc[train_index], X.iloc[test_index]
y_train, y_test = y[train_index], y[test_index]

model = RandomForestClassifier(n_estimators=nest, oob_score=False, n_jobs=threads, random_state=seed)
model.fit(X_train, y_train)

all_importances.append(pd.DataFrame({
  'Feature': X.columns,
  'Importance': model.feature_importances_
}) )
print(f"\t{time.time()-start_time} Sec.")
merged_importances = pd.concat(all_importances, axis=0, ignore_index=True)
importance_stats = ( merged_importances.groupby('Feature')['Importance'].agg(['mean', 'std']).rename(columns={'mean': 'Mean_Importance', 'std': 'Std_Importance'}).reset_index()
)
mean_importances = importance_stats.sort_values(by=['Mean_Importance', 'Std_Importance'], ascending=[False, True]).reset_index(drop=True)
mean_importances.index = mean_importances['Feature']
return(mean_importances)



def get_importance_lasso(X, y, repeat=10, seed=2024, k=10, threaeds=-1):
  '''
        有待完成
    '''

logging.info(f"Get features' Lasso importance. (N={repeat})")
model = LassoCV(cv=k, random_state=seed)
model.fit(X,y) # 只有建模以后，才能获取重要度
best_alpha = model.alpha_ # 获取最好的alpha

results = permutation_importance(model, X, y, n_repeats=repeat, random_state=seed, n_jobs=threaeds)

features_imp = pd.DataFrame(results.importances, X.columns)
# features_imp.to_csv("xxx.tsv", sep="\t")
features_ord = features_imp.mean(axis=1).sort_values(ascending=False).index
features_imp = features_imp.reindex(index=features_ord)
return(features_imp)


def format_data(features, metadata, Sample, Group, sstr, topn):
  
  features = os.path.abspath(features)

logging.info(f"Read file: {features}")
if topn == 0:
  topn = None
df = pd.read_csv(features, sep=sstr, header=0, index_col=0, nrows=topn, engine='python')
rawN_nf, rawN_df = df.shape # 获取原始矩阵的样本数量和特征数量

# dff = df.loc[:, ~(df==0).all(axis=0)] # 删除全为0的列
dff = df
cleanN_nf = dff.shape[0]

metadata = os.path.abspath(metadata)
logging.info(f"Read file: {metadata}")
sample_map = pd.read_csv(metadata, sep=sstr, header=0, engine='python')
rawN_sample = sample_map.shape[0] # 原始的行数

sampf = sample_map.loc[ :, (Sample, Group) ].drop_duplicates()
cleanN_sample = sampf.shape[0] # 去重以后的样本数量

print("\n\n--------------------------------")
print(f"Info\tFeatures\tSamples")
print("--------------------------------")
print(f"raw\t{rawN_nf}\t\t{rawN_df}")
print(f"format\t{cleanN_nf}\t\t{cleanN_sample}")

## 对齐分组和数据
intersect_id = dff.columns.intersection( sampf.loc[:, Sample] )
sampf = sampf.loc[sampf.loc[:,Sample].isin(intersect_id),:]
dff = dff.reindex(columns=sampf.loc[:, Sample]).T
alignN_sample, alignN_features = dff.shape

print(f"align\t{alignN_features}\t\t{alignN_sample}")
print("--------------------------------\n\n")

X = dff
y = sampf.loc[:,Group].reset_index(drop=True) # 如果不重置索引，后续筛选的时候会出错
return (X, y, alignN_features)



def main_lasso(args):
  '''
        Lasso模型
    '''
sstr = {0:r'\s+', 1:r'\t', 2:r'\|', 3:","}[args.s]
X, y, topn = format_data( features = args.i, metadata = args.g, Sample = args.sid, Group = args.gid, sstr=sstr, topn=args.topN)

outf = args.o
n_est = args.nest # 生成决策树的数量
threads = args.threads # 线程数
seed = args.seed

command = args.command

if command == "importance":
  X_importance = get_importance_lasso(X, y, args.R, seed, args.K, threads)

def main(args):
  '''
        因为所有的输入都需要对数据进行处理，所以这边先对数据和分组排序对齐
    '''
sstr = {0:r'\s+', 1:r'\t', 2:r'\|', 3:","}[args.s]

X, y, topn = format_data( features = args.i, metadata = args.g, Sample = args.sid, Group = args.gid, sstr=sstr, topn=args.topN)

outf = args.o
n_est = args.nest # 生成决策树的数量
threads = args.threads # 线程数
seed = args.seed

command = args.command


if command == "importance":
  if args.fast:
  X_importance = get_importance_rf_fast(X, y, args.R, seed, n_est, threads, args.K, args.stratified)
else:
  X_importance = get_importance_rf_perm(X, y, args.R, seed, n_est, threads, args.K, args.stratified) ## 获取特征的重要性
X = X.reindex(columns=X_importance.index)
X.T.to_csv(f"{outf}.sorted", sep="\t")
X_importance.to_csv(f"{outf}", sep="\t")
elif command == "KF":
  model, res = KF_predctions(X, y, args.sid, args.K, n_est, threads, seed, args.stratified, topn=topn)
res['nspecies'] = topn
res['seed'] = seed
res.to_csv(f"{outf}", sep="\t", index=False)
elif command == "LOOCV":
  pass
elif command == "easy":
  outm = args.om
with open (f"{outm}.pkl", 'wb') as file:
  pickle.dump(model, file)
pass
elif command == "predict":
  with open(f"{args.M}", 'rb') as file:
  model = pickle.load(file)
res = rf_predict(model, X, y, threads=threads, n_est=n_est, seed=seed, fillna=args.fillna)
res.to_csv(f"{outf}", sep="\t")
elif command == "build":
  outm = args.o
model = RandomForestClassifier(n_estimators=n_est,n_jobs=threads)
model.fit(X,y) # 只有建模以后，才能获取重要度
model.myfeatures = X.columns
with open (f"{outf}", 'wb') as file:
  pickle.dump(model, file)
else:
  exit(127)
logging.info("Done.")


if __name__ == "__main__":
  args = get_args()

import numpy as np
# import joblib # 用于保存随机森林模型
import pickle # 用于保存随机森林模型,只不过性能没有joblib好
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, KFold, StratifiedKFold, cross_val_score, cross_val_predict, RandomizedSearchCV
from sklearn.linear_model import LassoCV, Lasso # lasso线性模型
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.inspection import permutation_importance # 重要性置换检验
import pandas as pd
import dask.dataframe as dd
import random
import time
import logging

logging.basicConfig(level=logging.INFO,
                    filename=f"{args.o}.{args.command}.log", filemode="w",
                    format='%(asctime)s - %(levelname)s - %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)
cmd = ' '.join(sys.argv)

main(args)
# main_lasso(args)

