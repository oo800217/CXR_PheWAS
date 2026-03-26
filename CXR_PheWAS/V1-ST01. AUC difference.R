
library(magrittr)
library(data.table)
library(ggrepel)
library(glmnet)
library(ggplot2)
library(xgboost)
library(scales) 
library(survival)
library(cowplot)
library(gridExtra)
library(pROC)

# 0. settings

DATASET_NAME <- c('Hold-out dataset', 'External dataset', 'MIMIC dataset')

# 0. pathway

data_path <- 'data/explainable_aucs (NNLS-pathology).RData'

phecode_path <- 'data/summary.RData'

table_path <-  'result/STable 01.csv'

# 1. Load data

load(phecode_path)
load(data_path)

# 2. Processing

ref_phe <- binary_model_info[,c('phenotype', 'description', 'group')]
colnames(ref_phe) <- c('Phecode', 'Disease', 'Group')

AUC_data[,'AUC_diff'] <- AUC_data[,'auc.radio_less'] - AUC_data[,'auc.ori']
AUC_data[,'AUC_diff'] <- formatC(AUC_data[,'AUC_diff'], digits = 3, format = 'f')

AUC_data[,'Phecode'] <- gsub('\\-.*', '', AUC_data[,'y']) %>% gsub(' ', '', .)
AUC_data[,'Type'] <- ifelse(gsub('.*\\-', '', AUC_data[,'y']) %in% 'survival', 'Incident', 'Prevalent')
AUC_data[,'Order'] <- 1:nrow(AUC_data)

AUC_data <- merge(AUC_data, ref_phe, all.x = TRUE, by = 'Phecode')
AUC_data <- AUC_data[order(AUC_data[,'Order']),]

Write_AUC_data <- AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[1],c('y', 'cluster', 'Phecode', 'Type', 'Disease', 'Group')]

for (u in 1:length(DATASET_NAME)) {
  
  Write_AUC_data[,paste0('AUC(', DATASET_NAME[u], ')')] <- AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[u],'auc_txt.radio_less']
  Write_AUC_data[,paste0('Diff_AUC(', DATASET_NAME[u], ')')] <- AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[u],'AUC_diff']
  
}

## 3. Write out

write.csv(Write_AUC_data, table_path, na = '', row.names = FALSE, quote = TRUE, fileEncoding = 'CP950')
