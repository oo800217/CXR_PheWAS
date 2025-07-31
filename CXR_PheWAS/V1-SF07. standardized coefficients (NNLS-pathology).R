
library(scales) 
library(ggrepel)
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(Rtsne)
library(dplyr)

# 0. settings

ignore_cluster <- 12

cluster_names <- c('1. Future cardiovascular disease',
                   '2. Future critical illness',
                   '3. Cachexia and advanced malignancies',
                   '4. Advanced heart failure spectrum',
                   '5. Hemorrhage and coma',
                   '6. Sepsis and infarction',
                   '7. Cardiac arrest and respiratory failure',
                   '8. Dialysis-dependent disease',
                   '9. Advanced chronic kidney disease',
                   '10. Rheumatic lung disease',
                   '11. Osteoporosis',
                   '12. Prostate disorder',
                   '13. Pacemaker',
                   '14. Placental complications',
                   '15. Cardiogenic shock',
                   '16. Hypernatremic hyperosmolality',
                   '17. Valve replacement',
                   '18. Prevalent morbid obesity',
                   '19. Dementias',
                   '20. Empyema and pneumothorax',
                   '21. Pleurisy and effusion',
                   '22. Bronchial carcinoma',
                   '23. Ascites',
                   '24. Hydrocephalus',
                   '25. Fibrosing alveolitis',
                   '26. Laryngeal carcinoma',
                   '27. Future morbid obesity',
                   '28. Shock')

DATASET_NAME <- c('Hold-out dataset', 'External dataset', 'MIMIC dataset')

# 0. pathway

radio_order_path <-  'data/radio_order.RData'


data_path <- 'data/explainable_aucs (NNLS-pathology).RData'

plot_path <- 'result/Fig S07.pdf'

# 1. Load data

load(data_path)
load(radio_order_path)

radio_list <- radio_list.less
rsq_name <- 'rsq_less'



# 2. Preprocessing

radio_coef <- do.call('cbind', radio_list)
colnames(radio_coef) <- names(radio_list)

radio_coef <- radio_coef * matrix(sd_list.radio, nrow = nrow(radio_coef), ncol = ncol(radio_coef), byrow = FALSE)
radio_coef <- radio_coef / t(matrix(sd_list.phe, nrow = ncol(radio_coef), ncol = nrow(radio_coef), byrow = FALSE))

radio_order_info <- radio_order_info[radio_order_info[,'phenotype'] %in% rownames(radio_coef),]

radio_coef[radio_coef > 1] <- 1
radio_coef[radio_coef <= 0] <- 0

col_range <- c(0.005, 1)

my_colors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))[6:11] %>% paste0(., 'D0')

my_colors <- c(
  # "#08306B",  # 深藍
  # "#2171B5",
  # "#4292C6",
  # "#6BAED6",
  # "#9ECAE1",
  # "#FFFFFF",  # 淡藍接近白
  "#FFFFFF",  # 淡藍接近白
  "#FEE0D2",  # 淡紅接近白
  "#FC9272",
  "#FB6A4A",
  "#EF3B2C",
  "#99000D"   # 深紅
) 

# 3. Clustering for radio


features_order <- radio_order_info[,'phenotype']

# 3. Plottings (For all)

colour_list <- hue_pal()(ceiling(ignore_cluster / 3) * 3)
colour_list <- matrix(colour_list, nrow = 3, byrow = TRUE) %>% as.character()

cluster_df.radio <- NULL

for (m in 1:max(radio_order_info[,'new_cluster.radio'])) {
  
  if (m %in% radio_order_info[,'new_cluster.radio']) {
    
    cluster_df.radio <- rbind(cluster_df.radio, data.frame(cluster = m,
                                                           start = min(which(radio_order_info[,'new_cluster.radio'] %in% m)),
                                                           end = max(which(radio_order_info[,'new_cluster.radio'] %in% m))))
    
  }
  
}

auc_df <- AUC_data[,c('dataset', 'y', rsq_name)]
colnames(auc_df) <- c("Var1", "Var2", "value")

cor_df <- melt(radio_coef)
colnames(cor_df) <- c("Var1", "Var2", "value")

cor_df <- rbind(cor_df, auc_df)

cor_df[,'Var1'] <- factor(cor_df[,'Var1'], levels = c(features_order, 'Hold-out dataset', 'External dataset', 'MIMIC dataset'))
levels(cor_df[,'Var1'])[!levels(cor_df[,'Var1']) %in% features_order] <- levels(cor_df[,'Var1'])[!levels(cor_df[,'Var1']) %in% features_order] %>% gsub(' dataset', '', .) %>% paste0('R-quare (', ., ')')

cor_df[,'Var2'] <- factor(cor_df[,'Var2'], levels = unique(AUC_data[,'y']))
levels(cor_df[,'Var2']) <- levels(cor_df[,'Var2']) %>% gsub('\\-', ' (', .) %>% gsub('survival', 'incident) [cluster ', .) %>% gsub('binary', 'prevalent) [cluster ', .) %>% paste0(., AUC_data[!duplicated(AUC_data[,'y']),'cluster'], ']')

coef_p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value))

coef_p <- coef_p + geom_tile(color = "white")
coef_p <- coef_p + scale_fill_gradientn(colors = my_colors, limits = col_range, name = "Standarized regression coefficients\nor R-square", na.value = "grey50")

coef_p <- coef_p + labs(title = '', x = '', y = '')
coef_p <- coef_p + theme_minimal(base_size = 8)

coef_p <- coef_p + theme(plot.title = element_blank(),
                             axis.text.x = element_text(angle = 45, size = 7, vjust = 1, hjust = 1, color = 'black'),
                             axis.text.y = element_text(size = 7, hjust = 1, color = 'black'),
                             legend.position = "bottom",
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             panel.grid = element_blank())

for (i in 1:nrow(cluster_df.radio)) {
  
  if (i <= ignore_cluster & (cluster_df.radio[i,'end'] - cluster_df.radio[i,'start']) > 0) {current_col <- colour_list[i]} else {current_col <- 'gray50'}
  
  coef_p <- coef_p + annotate(geom = "rect",
                                  xmin = cluster_df.radio[i,'start'] - 0.49, xmax = cluster_df.radio[i,'end'] + 0.49,
                                  ymin = dim(radio_coef)[2] + 1, ymax = dim(radio_coef)[2] + 2.5,
                                  color = "black",
                                  fill = current_col, size = 0.5, alpha = 0.5)
  
  coef_p <- coef_p + annotate(geom = "rect",
                                  xmin = cluster_df.radio[i,'start'] - 0.49, xmax = cluster_df.radio[i,'end'] + 0.49,
                                  ymin = 0, ymax = -1.5,
                                  color = "black",
                                  fill = current_col, size = 0.5, alpha = 0.5)
  
}

## Cluster 

sub_AUC_data <- AUC_data[AUC_data[,'dataset'] %in% 'Hold-out dataset',]
cluster_df <- NULL

for (m in unique(sub_AUC_data[,'cluster'])) {
  
  cluster_df <- rbind(cluster_df, data.frame(cluster = m,
                                             start = min(which(sub_AUC_data[,'cluster'] %in% m)),
                                             end = max(which(sub_AUC_data[,'cluster'] %in% m))))
  
}

for (i in 1:nrow(cluster_df)) {
  
  if (i <= ignore_cluster) {current_col <- colour_list[i]} else {current_col <- 'gray50'}
  
  coef_p <- coef_p + annotate(geom = "rect",
                                  ymin = cluster_df[i,'start'] - 0.49, ymax = cluster_df[i,'end'] + 0.49,
                                  xmin = length(levels(cor_df[,'Var1'])) + 1, xmax = length(levels(cor_df[,'Var1'])) + 2,
                                  color = "black",
                                  fill = current_col, size = 0.5, alpha = 0.5)
  
  coef_p <- coef_p + annotate(geom = "rect",
                                  ymin = cluster_df[i,'start'] - 0.49, ymax = cluster_df[i,'end'] + 0.49,
                                  xmin = 0, xmax = -1,
                                  color = "black",
                                  fill = current_col, size = 0.5, alpha = 0.5)
}

coef_p <- coef_p + annotate(geom = "rect",
                                ymin = 0.5, ymax = length(levels(cor_df[,'Var2'])) + 0.5,
                                xmin = length(levels(cor_df[,'Var1'])) - 2.5, xmax = length(levels(cor_df[,'Var1'])) + 0.5,
                                color = "black",
                                fill = NA, size = 0.5, alpha = 0.5)


# 3. Merge plot


pdf(plot_path, width = 12, height = 12)

print(coef_p)

dev.off()
