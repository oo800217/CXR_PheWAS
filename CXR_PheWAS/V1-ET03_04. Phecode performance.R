
library(magrittr)
library(data.table)
library(ggrepel)
library(glmnet)
#library(doMC)
library(ggplot2)
library(xgboost)
library(scales) 
library(survival)
library(cowplot)
library(gridExtra)
library(pROC)

# 0. settings

DATASET_NAME <- c('Internal hold-out set' = 'internal', 'External dataset' = 'external', 'MIMIC dataset' = 'MIMIC_png')

# 0. pathway

data_path <- 'data/summary.RData'

table_path.1 <-  'result/ETable 03 (binary).csv'
table_path.2 <-  'result/ETable 04 (survival).csv'

# 1. Load data

load(data_path)

## 2. Preprocessing

binary_model_info <- binary_model_info[,c('phenotype', 'PhecodeString', 'group', 'internal.auc_txt', 'internal.auc.comb_txt', 'internal.auc.comb_pval', 'external.auc_txt', 'external.auc.comb_txt', 'external.auc.comb_pval', 'MIMIC_png.auc_txt', 'MIMIC_png.auc.comb_txt', 'MIMIC_png.auc.comb_pval')]

## 3. Preprocessing

survival_model_info <- survival_model_info[,c('phenotype', 'PhecodeString', 'group', 'internal.auc_txt', 'internal.auc.comb_txt', 'internal.auc.comb_pval', 'external.auc_txt', 'external.auc.comb_txt', 'external.auc.comb_pval', 'MIMIC_png.auc_txt', 'MIMIC_png.auc.comb_txt', 'MIMIC_png.auc.comb_pval')]

## 4. Write out

write.csv(binary_model_info, table_path.1, na = '', row.names = FALSE, quote = TRUE, fileEncoding = 'CP950')
write.csv(survival_model_info, table_path.2, na = '', row.names = FALSE, quote = TRUE, fileEncoding = 'CP950')
