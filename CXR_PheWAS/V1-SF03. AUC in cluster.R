
library(scales) 
library(ggrepel)
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(Rtsne)
library(dplyr)

# 0. settings

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

data_path <- 'data/cluster_aucs.RData'

plot_path <- 'result/Fig S03.pdf'

# 1. Load data

load(data_path)

# 2. Plottings

my_colors <- rev(RColorBrewer::brewer.pal(11, "Spectral")) %>% paste0(., 'D0')

delta_p_list <- list()

for (i in 1:length(DATASET_NAME)) {
  
  sub_AUC_data <- AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[i],]
  
  table_cluster <- table(sub_AUC_data[,'cluster'])
  table_cluster <- table_cluster[table_cluster > 1]
  
  for (j in 1:length(table_cluster)) {
    
    current_cluster <- names(table_cluster)[j]
    
    if (i == 1) {delta_p_list[[paste0('Cluster ', current_cluster)]] <- list()}
    
    cor_df <- sub_AUC_data[sub_AUC_data[,'cluster'] %in% current_cluster,]
    
    same_df <- cor_df[cor_df[,'y'] == cor_df[,'x'],c('y', 'auc')]
    colnames(same_df)[2] <- 'base_auc'
    
    cor_df <- merge(cor_df, same_df, by = 'y', all.x = TRUE)
    cor_df[,'value'] <- cor_df[,'auc'] - cor_df[,'base_auc']
    
    disease_order <- AUC_data[,'y'] %>% unique()
    disease_order <- disease_order[disease_order %in% cor_df[,'y']]
    
    cor_df[,'Var1'] <- factor(cor_df[,'x'], levels = disease_order)
    cor_df[,'Var2'] <- factor(cor_df[,'y'], levels = disease_order)
    
    levels(cor_df[,'Var1']) <- levels(cor_df[,'Var1']) %>% gsub('\\-', ' (', .) %>% gsub('survival', 'incident)', .) %>% gsub('binary', 'prevalent)', .)
    levels(cor_df[,'Var2']) <- levels(cor_df[,'Var2']) %>% gsub('\\-', ' (', .) %>% gsub('survival', 'incident)', .) %>% gsub('binary', 'prevalent)', .)
    
    cor_df[cor_df[,'value'] < -0.05,'value'] <- -0.05
    cor_df[cor_df[,'value'] > 0.05,'value'] <- 0.05
    
    ## Plottings
    
    cosine_p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value))
    
    cosine_p <- cosine_p + geom_tile(color = "white")
    cosine_p <- cosine_p + scale_fill_gradientn(colors = my_colors, limits = c(-0.05, 0.05), name = "Delta AUC/C-index")
    
    cosine_p <- cosine_p + labs(title = paste0('Cluster ', current_cluster, '\n(', gsub(' dataset', '', DATASET_NAME[i]), ')'), x = 'Model', y = 'Outcome')
    cosine_p <- cosine_p + theme_minimal(base_size = 8)
    
    cosine_p <- cosine_p + theme(plot.title = element_text(size = 10, color = 'black'),
                                 axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5, hjust = 1, color = 'black'),
                                 axis.text.y = element_text(size = 7, hjust = 1, color = 'black'),
                                 legend.position = "none",
                                 axis.title.x = element_text(size = 8, color = 'black'),
                                 axis.title.y = element_text(size = 8, color = 'black'),
                                 panel.grid = element_blank())
    
    cosine_p <- cosine_p + coord_equal()
    
    delta_p_list[[paste0('Cluster ', current_cluster)]][[DATASET_NAME[i]]] <- cosine_p
    
  }
  
}

cor_df[,'value'] <- -5

legend_p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value))

legend_p <- legend_p + geom_tile(color = "white")
legend_p <- legend_p + scale_fill_gradientn(colors = my_colors, limits = c(-0.05, 0.05), name = "Delta AUC/C-index", na.value = "white")

legend_p <- legend_p + labs(title = paste0('Cluster ', current_cluster, ' (', DATASET_NAME[i], ')'), x = 'Model', y = 'Outcome')
legend_p <- legend_p + theme_minimal(base_size = 8)

legend_p <- legend_p + theme(plot.title = element_blank(),
                             axis.text.x = element_blank(),
                             axis.text.y = element_blank(),
                             legend.position = "bottom",
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             panel.grid = element_blank())

legend_p <- legend_p + coord_equal()

# 3. Merge plot

merge_plot.list <- list()

for (u in 1:length(delta_p_list)) {
  
  if (u <= 2 | u %in% 10:11) {
    
    merge_plot.list[[u]] <- arrangeGrob(delta_p_list[[u]][[1]], delta_p_list[[u]][[2]], delta_p_list[[u]][[3]], 
                                        ncol = 3)  
    
  } else if (u %in% 12) {
    
    merge_plot.list[[u]] <- arrangeGrob(delta_p_list[[u]][[1]], delta_p_list[[u]][[2]], delta_p_list[[u]][[3]], 
                                        ncol = 2)  
    
  } else {
    
    merge_plot.list[[u]] <- arrangeGrob(delta_p_list[[u]][[1]], delta_p_list[[u]][[2]], delta_p_list[[u]][[3]], 
                                        ncol = 1)  
    
  }
  
}

final_p <- ggdraw()
final_p <- final_p + draw_plot(legend_p, x = 0.4, y = 0.05, width = 1, height = 1)
final_p <- final_p + draw_plot(merge_plot.list[[1]], x = 0, y = 0.7, width = 1, height = 0.3)
final_p <- final_p + draw_plot(merge_plot.list[[2]], x = 0, y = 0.5, width = 0.65, height = 0.2)
final_p <- final_p + draw_plot(merge_plot.list[[3]], x = 0.65, y = 0.25, width = 0.175, height = 0.45)
final_p <- final_p + draw_plot(merge_plot.list[[4]], x = 0.825, y = 0.25, width = 0.175, height = 0.45)
final_p <- final_p + draw_plot(merge_plot.list[[5]], x = 0, y = 0.12, width = 0.13, height = 0.37)
final_p <- final_p + draw_plot(merge_plot.list[[6]], x = 0.13, y = 0.12, width = 0.13, height = 0.37)
final_p <- final_p + draw_plot(merge_plot.list[[7]], x = 0.26, y = 0.12, width = 0.13, height = 0.37)
final_p <- final_p + draw_plot(merge_plot.list[[8]], x = 0.39, y = 0.12, width = 0.13, height = 0.37)
final_p <- final_p + draw_plot(merge_plot.list[[9]], x = 0.52, y = 0.12, width = 0.13, height = 0.37)
final_p <- final_p + draw_plot(merge_plot.list[[10]], x = 0, y = 0, width = 0.333, height = 0.11)
final_p <- final_p + draw_plot(merge_plot.list[[11]], x = 0.333, y = 0, width = 0.333, height = 0.11)
final_p <- final_p + draw_plot(merge_plot.list[[12]], x = 0.70, y = 0.02, width = 0.222, height = 0.22)


if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

pdf(plot_path, width = 15, height = 16)

print(final_p)

dev.off()
