
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

data_path <- 'data/explainable_aucs (NNLS-all).RData'

plot_path <- 'result/Fig S08.pdf'

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

dist_mat.radio <- dist(radio_coef, method = "euclidean")

hc.radio <- hclust(dist_mat.radio, method = "average")
features_order <- hc.radio[['labels']][hc.radio[['order']]]

features_order <- radio_order_info[,'phenotype']
features_order


# 4. Plottings (For cluster)

new_radio_coef <- matrix(0, nrow = dim(radio_coef)[1], ncol = length(unique(AUC_data[,'cluster'])))
rownames(new_radio_coef) <- rownames(radio_coef)
colnames(new_radio_coef) <- cluster_names

for (i in 1:dim(radio_coef)[1]) {
  
  new_radio_coef[i,] <- tapply(radio_coef[i,], AUC_data[!duplicated(AUC_data[,'y']),'cluster'], mean)
  
}

# new_radio_coef <- new_radio_coef[which(apply(new_radio_coef, 1, max) > 0.005),]

for (i in 1:ncol(new_radio_coef)) {
  
  pos <- which(new_radio_coef[,i] >= 0.1)
  if (length(pos) > 3) {pos <- pos[1:3]}
  
  message('cluster ', i, ' (', gsub('^[0-9]{1,}\\. ', '', tolower(cluster_names[i])), ') by ', paste(paste0(tolower(rownames(new_radio_coef)[pos]), ' (', formatC(new_radio_coef[pos,i], digits = 2, format = 'f'), ')')[order(new_radio_coef[pos,i], decreasing = TRUE)], collapse = ', '))
  
}

mean_auc <- tapply(AUC_data[,rsq_name], AUC_data[,c('cluster', 'dataset')], mean) %>% t()

auc_df <- melt(mean_auc)
colnames(auc_df) <- c("Var1", "Var2", "value")
auc_df[,'Var2'] <- cluster_names[auc_df[,'Var2']]

cor_df <- melt(new_radio_coef)
colnames(cor_df) <- c("Var1", "Var2", "value")

cor_df <- rbind(cor_df, auc_df)

cor_df[,'Var1'] <- factor(cor_df[,'Var1'], levels = c(features_order, 'Hold-out dataset', 'External dataset', 'MIMIC dataset'))
levels(cor_df[,'Var1'])[!levels(cor_df[,'Var1']) %in% features_order] <- levels(cor_df[,'Var1'])[!levels(cor_df[,'Var1']) %in% features_order] %>% gsub(' dataset', '', .) %>% paste0('Mean R-quare (', ., ')')

cor_df[,'Var2'] <- factor(cor_df[,'Var2'], levels = cluster_names)

coef_p <- ggplot(cor_df, aes(x = Var2, y = Var1, fill = value))

coef_p <- coef_p + geom_tile(color = "white")
coef_p <- coef_p + scale_fill_gradientn(colors = my_colors, limits = col_range, name = "Standarized regression coefficients\nor R-square", na.value = "grey50")
coef_p <- coef_p + geom_text(aes(label = ifelse(value < 0.005, '', sprintf("%.2f", value))), size = 2, color = "black")

coef_p <- coef_p + labs(title = '', x = '', y = '')
coef_p <- coef_p + theme_minimal(base_size = 8)

coef_p <- coef_p + theme(plot.title = element_blank(),
                             axis.text.x = element_text(angle = 45, size = 7, vjust = 1, hjust = 1, color = 'black'),
                             axis.text.y = element_text(size = 7, hjust = 1, color = 'black'),
                             legend.position = "bottom",
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             panel.grid = element_blank())

coef_p <- coef_p + annotate(geom = "rect",
                                xmin = 0.5, xmax = length(levels(cor_df[,'Var2'])) + 0.5,
                                ymin = length(levels(cor_df[,'Var1'])) - 2.5, ymax = length(levels(cor_df[,'Var1'])) + 0.5,
                                color = "black",
                                fill = NA, size = 0.5, alpha = 0.5)

for (i in 1:nrow(cluster_df.radio)) {
  
  if (i <= ignore_cluster & (cluster_df.radio[i,'end'] - cluster_df.radio[i,'start']) > 0) {current_col <- colour_list[i]} else {current_col <- 'gray50'}
  
  coef_p <- coef_p + annotate(geom = "rect",
                                      ymin = cluster_df.radio[i,'start'] - 0.49, ymax = cluster_df.radio[i,'end'] + 0.49,
                                      xmin = dim(new_radio_coef)[2] + 0.7, xmax = dim(new_radio_coef)[2] + 1.2,
                                      color = "black",
                                      fill = current_col, size = 0.5, alpha = 0.5)
  
  coef_p <- coef_p + annotate(geom = "rect",
                                      ymin = cluster_df.radio[i,'start'] - 0.49, ymax = cluster_df.radio[i,'end'] + 0.49,
                                      xmin = 0.3, xmax = -0.2,
                                      color = "black",
                                      fill = current_col, size = 0.5, alpha = 0.5)
  
}

# 3. Merge plot

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

pdf(plot_path, width = 10, height = 10)

print(coef_p)

dev.off()

