
library(scales) 
library(ggrepel)
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(Rtsne)
library(dplyr)

# 0. settings

DATASET_NAME <- c('internal' = 'Hold-out dataset', 'external' = 'External dataset', 'MIMIC_png' = 'MIMIC dataset')
TIME_INT_NAME <- c('(30, 90]', '(90, 365]', '(365, 1095]', '(1095, Inf]')

# 0. pathway

data_path <- 'data/aucs_by_time.RData'
phecode_path <- 'data/summary.RData'

plot_path <- 'result/Fig S04.pdf'

# 1. Load data

load(phecode_path)
load(data_path)

# 2. Preprocessing

survival_model_info[,'y'] <- paste0(survival_model_info[,'label'], '-', survival_model_info[,'label_type'])
baseline_info <- survival_model_info[,c('y', 'internal.auc', 'external.auc', 'MIMIC_png.auc')]

AUC_data[AUC_data[,'time_interval'] %in% '(-Inf, 90]','time_interval'] <- '(30, 90]'
AUC_data <- merge(AUC_data, baseline_info, by = 'y', all.x = TRUE)

for (u in 1:length(DATASET_NAME)) {
  
  AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[u],'auc_diff'] <- AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[u],'auc'] - AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[u],paste0(names(DATASET_NAME)[u], '.auc')]  
    
}

AUC_data[,'time_interval'] <- factor(AUC_data[,'time_interval'], levels = TIME_INT_NAME)

# 3. Define color

unique_disease <- c()

for (u in 1:length(DATASET_NAME)) {
  
  sub_AUC_data <- AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[u],] 
  diff_range <- tapply(sub_AUC_data[,'auc_diff'], sub_AUC_data[,'y'], range, na.rm = TRUE) %>% lapply(., diff) %>% sapply(., abs)
  
  current_names <- ifelse(sub_AUC_data[,'y'] %in% names(diff_range)[diff_range > 0.1], paste0(gsub('-survival', '', sub_AUC_data[,'y']), ' [cluster ', sub_AUC_data[,'cluster'], ']'), 'Other')
  
  unique_disease <- c(unique_disease, current_names[!current_names %in% 'Other'])
  
}

unique_disease <- unique_disease %>% unique()

my_colour_list <- hue_pal()(length(unique_disease))
names(my_colour_list) <- unique_disease

# 4. Plottings

delta_p_list <- list()

for (u in 1:length(DATASET_NAME)) {
  
  sub_AUC_data <- AUC_data[AUC_data[,'dataset'] %in% DATASET_NAME[u],] 
  diff_range <- tapply(sub_AUC_data[,'auc_diff'], sub_AUC_data[,'y'], range, na.rm = TRUE) %>% lapply(., diff) %>% sapply(., abs)
  
  sub_AUC_data[,'with_col'] <- ifelse(sub_AUC_data[,'y'] %in% names(diff_range)[diff_range > 0.1], paste0(gsub('-survival', '', sub_AUC_data[,'y']), ' [cluster ', sub_AUC_data[,'cluster'], ']'), 'Other')

  colour_list <- c(my_colour_list[names(my_colour_list) %in% sub_AUC_data[,'with_col']], 'Other' = '#A0A0A0')
  
  delta_p <- ggplot(sub_AUC_data, aes(x = time_interval, y = auc_diff, group = y))
  delta_p <- delta_p + geom_abline(slope = 0, intercept = 0, colour = 'blue', linetype = 3)
  
  delta_p <- delta_p + geom_line(data = subset(sub_AUC_data, with_col == "Other"), aes(color = with_col), alpha = 0.3, size = 0.5)
  delta_p <- delta_p + geom_line(data = subset(sub_AUC_data, with_col != "Other"), aes(color = with_col), alpha = 0.7, size = 1)
  
  delta_p <- delta_p + theme_classic()
  
  delta_p <- delta_p + labs(title = DATASET_NAME[u], x = 'Follow-up time interval', y = 'Delta C-index', color = "")
  
  delta_p <- delta_p + scale_y_continuous(expand = c(0, 0), limits = c(-0.2, 0.2),
                                          breaks = seq(-0.2, 0.2, by = 0.05))
  
  delta_p <- delta_p + theme(plot.title = element_text(color = "#000000", size = 10, face = 1),
                             legend.position = "bottom",
                             legend.title = element_blank(),
                             legend.text = element_text(size = 6),
                             axis.text.x = element_text(color = "#000000", size = 6),
                             axis.text.y = element_text(color = "#000000", size = 6),
                             axis.title.x = element_text(color = "#000000", size = 8),
                             axis.title.y = element_text(color = "#000000", size = 8),
                             legend.text.align = 0) + labs(fill = '')
  
  delta_p <- delta_p + scale_color_manual(values = colour_list)
  delta_p <- delta_p + guides(color = guide_legend(ncol = c(2, 1, 3)[u], byrow = FALSE))
  
  delta_p_list[[u]] <- delta_p
  
}

# 3. Merge plot

merge_plot <- arrangeGrob(delta_p_list[[1]], delta_p_list[[2]], delta_p_list[[3]], 
                          ncol = 3)  


final_p <- ggdraw()
final_p <- final_p + draw_plot(merge_plot, x = 0, y = 0, width = 1, height = 1)

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

pdf(plot_path, width = 10, height = 6)

print(final_p)

dev.off()
