
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

cut_internal <- c(0.85, 0.8)
cut_external <- c(0.85, 0.8)
cut_MIMIC <- c(0.8, 0.75)

DATASET_NAME <- c('Hold-out dataset' = 'internal', 'External dataset' = 'external', 'MIMIC dataset' = 'MIMIC_png')

# 0. pathway

data_path <- 'data/summary.RData'

plot_path <- 'result/Fig 02.pdf'

# 1. Load data

load(data_path)

model_info <- rbind(binary_model_info[,colnames(survival_model_info)], survival_model_info)

# 2. Plot

all_p_list <- list()

for (q in 1:2) {
  
  if (q == 1) {
    
    current_type <- 'binary'
    y_title <- 'AUC'
    
  } else {
    
    current_type <- 'survival'
    y_title <- 'C-index'
    
  }
  
  sub_model_info <- model_info[model_info[,'label_type'] %in% current_type,c('group', 'phenotype', 'PhecodeString', 'internal.auc', 'internal.auc.comb_adj.pval', 'external.auc', 'external.auc.comb_adj.pval', 'MIMIC_png.auc', 'MIMIC_png.auc.comb_adj.pval')]
  sub_model_info[,'group'] <- factor(sub_model_info[,'group'], levels = unique(sub_model_info[,'group']))
  
  n_group <- length(levels(sub_model_info[,'group']))
  
  colour_list <- hue_pal()(n_group)
  colour_list <- hue_pal()(ceiling(n_group / 3) * 3)
  colour_list <- matrix(colour_list, nrow = 3, byrow = TRUE) %>% as.character()
  
  sub_model_info[,'colour'] <- sub_model_info[,'group']
  levels(sub_model_info[,'colour']) <- colour_list[1:n_group]
  sub_model_info[,'colour'] <- as.character(sub_model_info[,'colour'])
  
  sub_model_info <- sub_model_info[order(sub_model_info[,'internal.auc']),]
  sub_model_info <- sub_model_info[order(sub_model_info[,'group']),]
  
  sub_model_info[,'highlight.1'] <- (sub_model_info[,'internal.auc'] >= cut_internal[q] & sub_model_info[,'external.auc'] >= cut_external[q] & sub_model_info[,'MIMIC_png.auc'] >= cut_MIMIC[q]) + 0L
  sub_model_info[,'highlight.2'] <- (sub_model_info[,'internal.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info) & sub_model_info[,'external.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info) & sub_model_info[,'MIMIC_png.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info)) + 0L
  sub_model_info[,'highlight'] <- sub_model_info[,'highlight.1'] * sub_model_info[,'highlight.2']
  sub_model_info[,'x'] <- 1:nrow(sub_model_info)
  
  # print(table(sub_model_info[,c('group', 'highlight')]))
  
  gg_p_list <- list()
  
  for (u in 1:length(DATASET_NAME)) {
    
    main_txt <- names(DATASET_NAME)[u]
    
    new_df <- sub_model_info
    new_df[,'auc'] <- new_df[,paste0(DATASET_NAME[u], '.auc')]
    new_df[new_df[,'auc'] < 0.5,'auc'] <- 0.5
    
    
    gg_p <- ggplot(new_df, aes(x = x, y = auc, color = group))
    gg_p <- gg_p + geom_abline(slope = 0, intercept = rbind(cut_internal, cut_external, cut_MIMIC)[u,q], colour = 'gray70', linetype = 3, size = 0.4)
    gg_p <- gg_p + geom_point(size = 0.5, alpha = 0.5)
    gg_p <- gg_p + theme_classic()
    gg_p <- gg_p + labs(title = main_txt, x = '', y = y_title)
    gg_p <- gg_p + theme(plot.title = element_text(color = "#000000", size = 12, face = 1),
                         legend.position = "none",
                         legend.title = element_blank(),
                         legend.text = element_text(size = 10),
                         axis.ticks.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_text(color = "#000000", size = 7),
                         axis.title.x = element_blank(),
                         axis.title.y = element_text(color = "#000000", size = 9)) + labs(fill = '')
    
    gg_p <- gg_p + scale_color_manual(values = unique(new_df[,'colour']))
    gg_p <- gg_p + scale_y_continuous(limits = c(0.5, 1.1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = expression(paste('' <= 0.5), 0.6, 0.7, 0.8, 0.9, '1.0'))
    gg_p <- gg_p + scale_x_continuous(limits = c(0.5, max(new_df[,'x']) + 0.5), breaks = mean(new_df[,'x']))
    
    gg_p <- gg_p + annotate(geom = "point",
                            x = new_df[new_df[,'highlight'] %in% 1,'x'],
                            y = new_df[new_df[,'highlight'] %in% 1,'auc'],
                            color = "#00000080",
                            fill = paste0(new_df[new_df[,'highlight'] %in% 1,'colour'], 'A0'),
                            size = 1.5, shape = 21)
    
    gg_p <- gg_p +  geom_text_repel(data = subset(new_df, highlight == 1),
                                    aes(x = x, y = auc, label = phenotype),   # label 欄位請換成你的文字欄位名稱
                                    nudge_y      = 0.01,                      # 垂直方向微調
                                    size         = 2,                        # 文字大小
                                    box.padding  = 0.3,                   # 文字框周圍保留空間
                                    point.padding = 0.5,                    # 點與文字框之間的距離
                                    segment.color = "#80808080",             # 連線顏色
                                    segment.size  = 0.2                      # 連線粗細
    )
    
    gg_p_list[[u]] <- gg_p
    
  }
  
  all_p_list[[q]] <- gg_p_list
  
}
  
## Legend plot

new_df[,'auc'] <- -1

legend_p <- ggplot(new_df, aes(x = x, y = auc, color = group))
legend_p <- legend_p + geom_point(size = 1, alpha = 1)
legend_p <- legend_p + theme_classic()
legend_p <- legend_p + labs(title = main_txt, x = '', y = y_title)
legend_p <- legend_p + theme(plot.title = element_blank(),
                             legend.position = "bottom",
                             legend.title = element_blank(),
                             legend.text = element_text(size = 8),
                             axis.text.x = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank()) + labs(fill = '')

legend_p <- legend_p + scale_y_continuous(limits = c(0.5, 1.05), breaks = NULL)
legend_p <- legend_p + scale_x_continuous(breaks = NULL)
legend_p <- legend_p + scale_color_manual(values = unique(new_df[,'colour']))

legend_p <- legend_p + guides(color = guide_legend(ncol = 9, byrow = FALSE))

# 3. Merge plot

merge_plot <- arrangeGrob(all_p_list[[1]][[1]], all_p_list[[2]][[1]],
                          all_p_list[[1]][[2]], all_p_list[[2]][[2]], 
                          all_p_list[[1]][[3]],  all_p_list[[2]][[3]], 
                          ncol = 2)  

final_p <- ggdraw()
final_p <- final_p + draw_plot(legend_p, x = 0, y = 0, width = 1, height = 1)
final_p <- final_p + draw_plot(merge_plot, x = 0, y = 0.06, width = 1, height = 0.89)
final_p <- final_p + draw_plot_label(c('Prevalent', 'Incident'), c(0.25, 0.75), c(0.98, 0.98), size = 14, hjust = 0.5, fontface = 1)
final_p <- final_p + draw_plot_label(letters[1:2], c(0.005, 0.505), c(1, 1), size = 15, hjust = 0)

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

pdf(plot_path, width = 12, height = 10)

print(final_p)

dev.off()

