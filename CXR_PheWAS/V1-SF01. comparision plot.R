
library(reshape2)
library(ggplot2)
library(gridExtra)
library(data.table)
library(cowplot)
library(dplyr)

# 0. settings

my_colors <- c('internal' = "#F8766D", 'external' = "#00BA38", 'MIMIC_png' = "#619CFF")

cut_internal <- c(0.85, 0.80) 
cut_external <- c(0.85, 0.80) 
cut_MIMIC <- c(0.80, 0.75) 

strat_var_list <- list()

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'CXR_VIEW_G', 
                                                      title_name = 'Imaging type',
                                                      lvl = c('PA', 'AP'),
                                                      lvl_name = expression('PA view', 'AP view'), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'CXR_POS_G', 
                                                      title_name = 'Position',
                                                      lvl = c('OPD/PEC', 'ER/IPD'),
                                                      lvl_name = expression('Outpatient', 'Emergency/Inpatient'), 
                                                      dataset = c('internal', 'external'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'RACE_G', 
                                                      title_name = 'Race',
                                                      lvl = c('ASIAN', 'WHITE'),
                                                      lvl_name = expression('Asian', 'White'), 
                                                      dataset = c('MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'RACE_G', 
                                                      title_name = 'Race',
                                                      lvl = c('ASIAN', 'BLACK'),
                                                      lvl_name = expression('Asian', 'Black'), 
                                                      dataset = c('MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'RACE_G', 
                                                      title_name = 'Race',
                                                      lvl = c('ASIAN', 'OTHER/UNKNOWN'),
                                                      lvl_name = expression('Asian', 'Other'), 
                                                      dataset = c('MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'GENDER_G', 
                                                      title_name = 'Gender',
                                                      lvl = c('male', 'female'),
                                                      lvl_name = expression('Male', 'Female'), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'AGE_G', 
                                                      title_name = 'Age',
                                                      lvl = c('<50', '50-64'),
                                                      lvl_name = expression(paste('' < 50, ' y/o'), '50-64 y/o'), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'AGE_G', 
                                                      title_name = 'Age',
                                                      lvl = c('50-64', '>=65'),
                                                      lvl_name = expression('50-64 y/o', paste('' >= 65, ' y/o')), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'N_DIS_HIST_G', 
                                                      title_name = 'Number of disease histories',
                                                      lvl = c('0-2', '3-10'),
                                                      lvl_name = expression(paste('' <= 2, ' disease histories'), '3-10 disease histories'), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'N_DIS_HIST_G', 
                                                      title_name = 'Number of disease histories',
                                                      lvl = c('3-10', '>10'),
                                                      lvl_name = expression('3-10 disease histories', paste('' >= 11, ' disease histories')), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'N_DIS_PRES_G', 
                                                      title_name = 'Number of comorbidities',
                                                      lvl = c('0-2', '3-5'),
                                                      lvl_name = expression(paste('' <= 2, ' comorbidities'), '3-5 comorbidities'), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

strat_var_list[[length(strat_var_list) + 1L]] <- list(strat_var = 'N_DIS_PRES_G', 
                                                      title_name = 'Number of comorbidities',
                                                      lvl = c('3-5', '>5'),
                                                      lvl_name = expression('3-5 comorbidities', paste('' >= 6, ' comorbidities')), 
                                                      dataset = c('internal', 'external', 'MIMIC_png'))

DATASET_NAME <- c('Hold-out dataset' = 'internal', 'External dataset' = 'external', 'MIMIC dataset' = 'MIMIC_png')

# 0. pathway

data_path <- 'data/summary.RData'

plot_path <- 'result/Fig S01.pdf'

# 1. Load data

load(data_path)

# 2. Plotting

all_p_list <- list()

sub_data_list <- list()

for (q in 1:3) {
  
  if (q == 3) {
    
    sub_model_info <- do.call('rbind', sub_data_list)
    
   } else if (q == 1) {
  
     sub_model_info <- binary_model_info
  
   } else if (q == 2)  {
  
     sub_model_info <- survival_model_info
  
   }

  # Selected diseases

   if (q != 3) {
  
     sub_model_info[,'highlight_1.1'] <- (sub_model_info[,'internal.auc'] >= cut_internal[q] & sub_model_info[,'external.auc'] >= cut_external[q] & sub_model_info[,'MIMIC_png.auc'] >= cut_MIMIC[q]) + 0L
     sub_model_info[,'highlight_1.2'] <- (sub_model_info[,'internal.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info) & sub_model_info[,'external.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info) & sub_model_info[,'MIMIC_png.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info)) + 0L
     sub_model_info[,'highlight_1'] <- sub_model_info[,'highlight_1.1'] * sub_model_info[,'highlight_1.2']
  
     sub_model_info <- sub_model_info[sub_model_info[,'highlight_1'] %in% 1L,colnames(survival_model_info)]
  
     sub_data_list[[q]] <- sub_model_info
  
   }

  ## Selected

  gg_p_list <- list()
  
  for (i in 1:length(strat_var_list)) {
    
    scatter_df <- NULL
    
    for (j in 1:length(strat_var_list[[i]][['dataset']])) {
      
      scatter_df <- rbind(scatter_df, data.frame(dataset = strat_var_list[[i]][['dataset']][j],
                                                 x = sub_model_info[,paste0(strat_var_list[[i]][['dataset']][j], '.', strat_var_list[[i]][['strat_var']], '[', strat_var_list[[i]][['lvl']][1], '].auc')],
                                                 y = sub_model_info[,paste0(strat_var_list[[i]][['dataset']][j], '.', strat_var_list[[i]][['strat_var']], '[', strat_var_list[[i]][['lvl']][2], '].auc')],
                                                 stringsAsFactors = FALSE))
      
    }
    
    scatter_df[,'dataset'] <- factor(scatter_df[,'dataset'], levels = names(my_colors))
    scatter_df[,'diff'] <- scatter_df[,'x'] - scatter_df[,'y']
    
    gg_p <- ggplot(scatter_df, aes(x = x, y = y, color = dataset))
    
    gg_p <- gg_p + scale_x_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, by = 0.1))
    gg_p <- gg_p + scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, by = 0.1))
    gg_p <- gg_p + coord_equal()
    
    # gg_p <- gg_p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", size = 0.5)
    gg_p <- gg_p + geom_abline(slope = 1, intercept = 0.05, linetype = "dashed", color = "gray40", size = 0.5)
    gg_p <- gg_p + geom_abline(slope = 1, intercept = -0.05, linetype = "dashed", color = "gray40", size = 0.5)
    gg_p <- gg_p + geom_abline(slope = 1, intercept = -0.15, linetype = "dotted", color = "gray40", size = 0.5)
    gg_p <- gg_p + geom_abline(slope = 1, intercept = 0.15, linetype = "dotted", color = "gray40", size = 0.5)
    
    gg_p <- gg_p + annotate(geom = 'text',
                            x = 0.55,
                            y = 0.55,
                            label = paste0('', round(100 * mean(scatter_df[,'diff'] < 0.05 & scatter_df[,'diff'] >= -0.05, na.rm = TRUE)), '%'),
                            size = 3)
    
    gg_p <- gg_p + annotate(geom = 'text',
                            x = 0.65,
                            y = 0.55,
                            label = paste0('', round(100 * mean(scatter_df[,'diff'] < 0.15 & scatter_df[,'diff'] >= 0.05, na.rm = TRUE)), '%'),
                            size = 3)
    
    gg_p <- gg_p + annotate(geom = 'text',
                            x = 0.75,
                            y = 0.55,
                            label = paste0('', round(100 * mean( scatter_df[,'diff'] >= 0.15, na.rm = TRUE)), '%'),
                            size = 3)
    
    gg_p <- gg_p + annotate(geom = 'text',
                            x = 0.55,
                            y = 0.65,
                            label = paste0('', round(100 * mean(scatter_df[,'diff'] < -0.05 & scatter_df[,'diff'] >= -0.15, na.rm = TRUE)), '%'),
                            size = 3)
    
    gg_p <- gg_p + annotate(geom = 'text',
                            x = 0.55,
                            y = 0.75,
                            label = paste0('', round(100 * mean(scatter_df[,'diff'] < -0.15, na.rm = TRUE)), '%'),
                            size = 3)
    
    gg_p <- gg_p + geom_point(size = 1, alpha = 0.5)
    gg_p <- gg_p + scale_color_manual(values = my_colors)
    gg_p <- gg_p + theme_classic(base_size = 12)
    
    gg_p <- gg_p + labs(title = strat_var_list[[i]][['title_name']], x = strat_var_list[[i]][['lvl_name']][1], y = strat_var_list[[i]][['lvl_name']][2]) 
    
    gg_p <- gg_p + theme(plot.title = element_text(size = 12, color = 'black'),
                         legend.position = "none",
                         axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, color = 'black'),
                         axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, color = 'black'),
                         axis.title.x = element_text(color = "#000000", size = 10, face = 1),
                         axis.title.y = element_text(color = "#000000", size = 10, face = 1),
                         panel.grid = element_blank())
    
    gg_p_list[[i]] <- gg_p
    
  }
  
  all_p_list[[q]] <- gg_p_list
  
}

## Legend p

scatter_df[,'x'] <- -1
scatter_df[,'y'] <- -1
levels(scatter_df[,'dataset']) <- names(DATASET_NAME)
names(my_colors) <- names(DATASET_NAME)

legend_p <- ggplot(scatter_df, aes(x = x, y = y, color = dataset))

legend_p <- legend_p + geom_point(size = 2, alpha = 1)
legend_p <- legend_p + scale_color_manual(values = my_colors)

legend_p <- legend_p + theme(plot.title = element_blank(),
                             legend.position = "bottom",
                             legend.title = element_blank(),
                             legend.text = element_text(size = 10),
                             axis.text.x = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank()) + labs(fill = '')

legend_p <- legend_p + scale_y_continuous(limits = c(0.5, 1.05), breaks = NULL)
legend_p <- legend_p + scale_x_continuous(breaks = NULL)

# 3. Merge plot

merge_plot.1 <- arrangeGrob(all_p_list[[1]][[1]], all_p_list[[1]][[2]], all_p_list[[1]][[3]], all_p_list[[1]][[4]], 
                            all_p_list[[1]][[5]], all_p_list[[1]][[6]], all_p_list[[1]][[7]], all_p_list[[1]][[8]], 
                            all_p_list[[1]][[9]], all_p_list[[1]][[10]], all_p_list[[1]][[11]], all_p_list[[1]][[12]], 
                            ncol = 4)  

merge_plot.2 <- arrangeGrob(all_p_list[[2]][[1]], all_p_list[[2]][[2]], all_p_list[[2]][[3]], all_p_list[[2]][[4]], 
                            all_p_list[[2]][[5]], all_p_list[[2]][[6]], all_p_list[[2]][[7]], all_p_list[[2]][[8]], 
                            all_p_list[[2]][[9]], all_p_list[[2]][[10]], all_p_list[[2]][[11]], all_p_list[[2]][[12]], 
                            ncol = 4)  

merge_plot.3 <- arrangeGrob(all_p_list[[3]][[1]], all_p_list[[3]][[2]], all_p_list[[3]][[3]], all_p_list[[3]][[4]], 
                            all_p_list[[3]][[5]], all_p_list[[3]][[6]], all_p_list[[3]][[7]], all_p_list[[3]][[8]], 
                            all_p_list[[3]][[9]], all_p_list[[3]][[10]], all_p_list[[3]][[11]], all_p_list[[3]][[12]], 
                            ncol = 4)  

final_p.1 <- ggdraw()
final_p.1 <- final_p.1 + draw_plot(legend_p, x = 0.35, y = 0, width = 1, height = 1)
final_p.1 <- final_p.1 + draw_plot(merge_plot.1, x = 0, y = 0.05, width = 1, height = 0.95)

final_p.2 <- ggdraw()
final_p.2 <- final_p.2 + draw_plot(legend_p, x = 0.35, y = 0, width = 1, height = 1)
final_p.2 <- final_p.2 + draw_plot(merge_plot.2, x = 0, y = 0.05, width = 1, height = 0.95)

final_p.3 <- ggdraw()
final_p.3 <- final_p.3 + draw_plot(legend_p, x = 0.35, y = 0, width = 1, height = 1)
final_p.3 <- final_p.3 + draw_plot(merge_plot.3, x = 0, y = 0.05, width = 1, height = 0.95)

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

pdf(plot_path, width = 12, height = 10)

#print(final_p.1)
#print(final_p.2)
print(final_p.3)

dev.off()