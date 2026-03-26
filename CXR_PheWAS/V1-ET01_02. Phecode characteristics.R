
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

table_path.1 <-  'result/ETable 01 (binary).csv'
table_path.2 <-  'result/ETable 02 (survival).csv'

# 1. Load data

load(data_path)

## 2. Preprocessing

binary_model_info <- binary_model_info[,c('phenotype', 'PhecodeString', 'group', 'training.n', 'training.prev', 'internal.n', 'internal.prev', 'external.n', 'external.prev', 'MIMIC_png.n', 'MIMIC_png.prev')]

binary_model_info[,'training.case'] <- round(binary_model_info[,'training.n'] * binary_model_info[,'training.prev'])
binary_model_info[,'training.ctrl'] <- round(binary_model_info[,'training.n'] * (1 - binary_model_info[,'training.prev']))

binary_model_info[,'internal.case'] <- round(binary_model_info[,'internal.n'] * binary_model_info[,'internal.prev'])
binary_model_info[,'internal.ctrl'] <- round(binary_model_info[,'internal.n'] * (1 - binary_model_info[,'internal.prev']))

binary_model_info[,'external.case'] <- round(binary_model_info[,'external.n'] * binary_model_info[,'external.prev'])
binary_model_info[,'external.ctrl'] <- round(binary_model_info[,'external.n'] * (1 - binary_model_info[,'external.prev']))

binary_model_info[,'MIMIC_png.case'] <- round(binary_model_info[,'MIMIC_png.n'] * binary_model_info[,'MIMIC_png.prev'])
binary_model_info[,'MIMIC_png.ctrl'] <- round(binary_model_info[,'MIMIC_png.n'] * (1 - binary_model_info[,'MIMIC_png.prev']))

binary_model_info <- binary_model_info[,c('phenotype', 'PhecodeString', 'group', 'training.case', 'training.ctrl', 'internal.case', 'internal.ctrl', 'external.case', 'external.ctrl', 'MIMIC_png.case', 'MIMIC_png.ctrl')]

## 3. Preprocessing

survival_model_info <- survival_model_info[,c('phenotype', 'PhecodeString', 'group', 'training.n', 'training.prev', 'training.preid', 'internal.n', 'internal.prev', 'internal.preid', 'external.n', 'external.prev', 'external.preid', 'MIMIC_png.n', 'MIMIC_png.prev', 'MIMIC_png.preid')]

survival_model_info[,'training.event'] <- round(survival_model_info[,'training.n'] * survival_model_info[,'training.prev'])
survival_model_info[,'internal.event'] <- round(survival_model_info[,'internal.n'] * survival_model_info[,'internal.prev'])
survival_model_info[,'external.event'] <- round(survival_model_info[,'external.n'] * survival_model_info[,'external.prev'])
survival_model_info[,'MIMIC_png.event'] <- round(survival_model_info[,'MIMIC_png.n'] * survival_model_info[,'MIMIC_png.prev'])

survival_model_info <- survival_model_info[,c('phenotype', 'PhecodeString', 'group', 'training.n', 'training.event', 'training.preid', 'internal.n', 'internal.event', 'internal.preid', 'external.n', 'external.event', 'external.preid', 'MIMIC_png.n', 'MIMIC_png.event', 'MIMIC_png.preid')]

## 4. Write out

write.csv(binary_model_info, table_path.1, na = '', row.names = FALSE, quote = TRUE, fileEncoding = 'CP950')
write.csv(survival_model_info, table_path.2, na = '', row.names = FALSE, quote = TRUE, fileEncoding = 'CP950')
