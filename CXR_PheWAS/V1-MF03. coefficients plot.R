
library(scales) 
library(ggrepel)
library(reshape2)
library(ggplot2)
library(Rtsne)
library(dplyr)

# 0. settings

n_clusters <- 28
ignore_cluster <- 12
n_clusters.radio <- 27 

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

Pathology_interested <- c(1:20, 44:57)

cut_internal <- c(0.85, 0.80) 
cut_external <- c(0.85, 0.80) 
cut_MIMIC <- c(0.80, 0.75) 

DATASET_NAME <- c('Hold-out dataset' = 'internal', 'External dataset' = 'external', 'MIMIC dataset' = 'MIMIC_png')

# 0. pathway

radio_path <- 'data/radio_feat.RData'
data_path <- 'data/summary.RData'
model_path <- 'data/model_info.RData'

plot_path <- 'result/Fig 03.pdf'

# 1. Load data

load(data_path)


# 2. Radio data

load(radio_path)

model_list <- model_list %>% data.frame()

radio_info <- model_list[,c('label_name'),drop = FALSE] %>% data.frame()

radio_coef <- model_list[,-1:-2] %>% as.matrix()
rownames(radio_coef) <- radio_info[,'label_name'] %>% gsub('[RADIO] ', '', ., fixed = TRUE) %>% gsub('.', ' ', ., fixed = TRUE)
rownames(radio_coef) <- paste0(toupper(substr(rownames(radio_coef), 1, 1)), tolower(substr(rownames(radio_coef), 2, nchar(rownames(radio_coef)))))
rownames(radio_coef)[rownames(radio_coef) %in% 'Icd implantation'] <- 'ICD implantation'
rownames(radio_coef)[rownames(radio_coef) %in% 'Picc insertion'] <- 'PICC insertion'

radio_info[,'phenotype'] <- rownames(radio_coef) 

Pathology_list <- radio_info[Pathology_interested,'phenotype']

radio_coef <- radio_coef * matrix(sd_list, nrow = nrow(radio_coef), ncol = ncol(radio_coef), byrow = TRUE)
  


# 1. Load data

load(model_path)

# 2. Merge data

sub_model_info_list <- list()

for (q in 1:2) {
  
  if (q == 1) {
    
    sub_model_info <- binary_model_info

  } else {
    
    sub_model_info <- survival_model_info

  }
  
  ## Selected diseases
  
  sub_model_info[,'highlight_1.1'] <- (sub_model_info[,'internal.auc'] >= cut_internal[q] & sub_model_info[,'external.auc'] >= cut_external[q] & sub_model_info[,'MIMIC_png.auc'] >= cut_MIMIC[q]) + 0L
  sub_model_info[,'highlight_1.2'] <- (sub_model_info[,'internal.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info) & sub_model_info[,'external.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info) & sub_model_info[,'MIMIC_png.auc.comb_adj.pval'] < 0.05 / nrow(survival_model_info)) + 0L
  sub_model_info[,'highlight_1'] <- sub_model_info[,'highlight_1.1'] * sub_model_info[,'highlight_1.2']
  
  sub_model_info_list[[q]] <- sub_model_info[sub_model_info[,'highlight_1'] %in% 1,c('phenotype', 'label', 'PhecodeString', 'group', 'label_type', 'internal.auc', 'internal.auc_txt', 'external.auc_txt', 'MIMIC_png.auc_txt')]
   
}

selected_model_info <- do.call('rbind', sub_model_info_list)

# 3. Get model_coef

selected_phecodes <- paste0(selected_model_info[,'label'], '-', selected_model_info[,'label_type'])
selected_model_list <- model_list[selected_phecodes]

selected_model_coef <- do.call('cbind', selected_model_list) %>% t()
selected_model_coef <- selected_model_coef[,-1]
rownames(selected_model_coef) <- selected_phecodes

selected_model_coef <- selected_model_coef * matrix(sd_list, nrow = nrow(selected_model_coef), ncol = ncol(selected_model_coef), byrow = TRUE)
  



## Cosine

normalize <- function(x) {x / sqrt(sum(x^2))}
mat_norm <- t(apply(selected_model_coef, 1, normalize))
cosine_sim <- mat_norm %*% t(mat_norm)
diag(cosine_sim) <- 1L

mat_norm.radio <- t(apply(radio_coef, 1, normalize))
cosine_sim.radio <- mat_norm.radio %*% t(mat_norm.radio)
diag(cosine_sim.radio) <- 1L

## T-sne

set.seed(0) # Sets seed for reproducibility

tsne_out <- Rtsne(rbind(mat_norm, mat_norm.radio), dims = 2, perplexity = 5, max_iter = 5000, eta = 200) # Run TSNE

selected_model_info[,'x'] <- tsne_out[['Y']][1:nrow(mat_norm),1]
selected_model_info[,'y'] <- tsne_out[['Y']][1:nrow(mat_norm),2]

radio_info[,'x'] <- tsne_out[['Y']][-1:-nrow(mat_norm),1]
radio_info[,'y'] <- tsne_out[['Y']][-1:-nrow(mat_norm),2]

## Clustering for features

dist_mat.radio <- dist(mat_norm.radio, method = "euclidean")

hc.radio <- hclust(dist_mat.radio, method = "average")
features_order <- hc.radio[['labels']][hc.radio[['order']]]

clusters.radio <- cutree(hc.radio, k = n_clusters.radio)
clusters_data.radio <- data.frame(phenotype = names(clusters.radio), cluster = as.integer(clusters.radio), stringsAsFactors = FALSE)

radio_info <- radio_info[hc.radio[['order']],]
radio_coef <- radio_coef[hc.radio[['order']],]

radio_info[,'order'] <- 1:nrow(radio_info)
radio_info <- merge(radio_info, clusters_data.radio, by = 'phenotype', all.x = TRUE)

radio_info <- radio_info[order(radio_info[,'order']),]
representatives_radio_info <- radio_info[!duplicated(radio_info[,'cluster']),]

table_cluster <- table(radio_info[,'cluster'])
table_cluster.data <- data.frame(cluster = names(table_cluster), Freq = as.numeric(table_cluster), stringsAsFactors = FALSE)
table_cluster.data <- merge(table_cluster.data, representatives_radio_info[,c('cluster', 'order')], by = 'cluster', all.x = TRUE)
table_cluster.data[,'Freq'] <- table_cluster.data[,'Freq'] - table_cluster.data[,'order'] * 1e-3
table_cluster.data <- table_cluster.data[order(table_cluster.data[,'Freq'], decreasing = TRUE),]
table_cluster.data[,'new_cluster.radio'] <- 1:nrow(table_cluster.data)

radio_info <- merge(radio_info, table_cluster.data[,c('cluster', 'new_cluster.radio')], by = 'cluster', all.x = TRUE)
radio_info <- radio_info[order(radio_info[,'order']),]

## Clustering for diseases

dist_mat <- dist(mat_norm, method = "euclidean")

hc <- hclust(dist_mat, method = "average")
disease_order <- hc[['labels']][hc[['order']]]

clusters <- cutree(hc, k = n_clusters)
clusters_data <- data.frame(label = names(clusters), cluster = as.integer(clusters), stringsAsFactors = FALSE)
clusters_data[,'label_type'] <- clusters_data[,'label'] %>% gsub('.*\\-', '', .)
clusters_data[,'label'] <- clusters_data[,'label'] %>% gsub('\\-.*', '', .)

selected_model_info <- selected_model_info[hc[['order']],]
selected_model_info[,'order'] <- 1:nrow(selected_model_info)
selected_model_info <- merge(selected_model_info, clusters_data, by = c('label', 'label_type'), all.x = TRUE)
selected_model_info <- selected_model_info[!selected_model_info[,'cluster'] %in% NA,]
selected_model_info <- selected_model_info[order(selected_model_info[,'order']),]

selected_model_info <- selected_model_info[order(selected_model_info[,'internal.auc'], decreasing = TRUE),]
representatives_info <- selected_model_info[!duplicated(selected_model_info[,'cluster']),]
phecode_info <- selected_model_info[!duplicated(selected_model_info[,'phenotype']),]
phecode_info[,'high_auc'] <- phecode_info[,'internal.auc']

table_cluster <- table(selected_model_info[,'cluster'])
table_cluster.data <- data.frame(cluster = names(table_cluster), Freq = as.numeric(table_cluster), stringsAsFactors = FALSE)
table_cluster.data <- merge(table_cluster.data, representatives_info[,c('cluster', 'internal.auc')], by = 'cluster', all.x = TRUE)
table_cluster.data[,'Freq'] <- table_cluster.data[,'Freq'] + table_cluster.data[,'internal.auc']
table_cluster.data <- table_cluster.data[order(table_cluster.data[,'Freq'], decreasing = TRUE),]
table_cluster.data[,'new_cluster'] <- 1:nrow(table_cluster.data)
table_cluster.data[3:4,'new_cluster'] <- 4:3
table_cluster.data[5:6,'new_cluster'] <- 6:5
table_cluster.data[9:11,'new_cluster'] <- c(10, 11, 9)

selected_model_info <- merge(selected_model_info, table_cluster.data[,c('cluster', 'new_cluster')], by = 'cluster', all.x = TRUE)
selected_model_info <- merge(selected_model_info, phecode_info[,c('phenotype', 'high_auc')], by = 'phenotype', all.x = TRUE)

selected_model_info <- selected_model_info[order(selected_model_info[,'order']),]
selected_model_info[,'label_type'] <- factor(selected_model_info[,'label_type'])
levels(selected_model_info[,'label_type']) <- c('Prevalent', 'Incident')
selected_model_info[,'label_type'] <- as.character(selected_model_info[,'label_type'])
  
#selected_model_info[,c('group', 'label_type', 'new_cluster')] %>% table %>% print()
#selected_model_info[,c('label_type', 'new_cluster')] %>% table %>% print()



# 5. T-sne plot

radio_info[,'new_cluster'] <- 0L
radio_info[radio_info[,c('phenotype')] %in% Pathology_list,'new_cluster'] <- -1L

radio_info[,'label_type'] <- 'radio'

rad_df <- radio_info[,c('phenotype', 'label_type', 'x', 'y', 'new_cluster')]
phe_df <- selected_model_info[,c('phenotype', 'label_type', 'x', 'y', 'new_cluster')]

df <- rbind(rad_df, phe_df)


## 建立顏色與分群標籤

colour_list <- hue_pal()((ceiling(ignore_cluster / 3) + 0) * 3)

colour_list <- matrix(colour_list, nrow = 3, byrow = TRUE) %>% as.character()

df[,'cluster_group'] <- ifelse(df[,'new_cluster'] <= ignore_cluster, df[,'new_cluster'], ignore_cluster + 1) %>% factor(levels = -1:(ignore_cluster + 1))
levels(df[,'cluster_group'])[1] <- 'Pathological changes or anatomical/morphological signs'
levels(df[,'cluster_group'])[2] <- 'Implants/devices or surgical scars/traces'
levels(df[,'cluster_group'])[ignore_cluster + 3] <- paste0('Other clusters (', ignore_cluster + 1, '-', n_clusters, ')')
df[,'color'] <- df[,'cluster_group']
levels(df[,'color']) <- c('#202020', '#808080', colour_list[1:ignore_cluster], '#000000')

df[,'label_type'] <- factor(df[,'label_type'])
levels(df[,'label_type']) <- c('Incident', 'Prevalent', 'Radiological')

levels(df[,'cluster_group'])[3:14] <- cluster_names[1:ignore_cluster]

## 凸包

sub_df <- df[df[,'new_cluster'] <= ignore_cluster & df[,'new_cluster'] >= 1,]

cluster_list <- split(sub_df, sub_df[,'new_cluster'])

# 使用 base R 建立每個 cluster 的 convex hull

hull_list <- lapply(cluster_list, function (cluster_data, scale = 1.2, min_val = 1) {
  
  cluster_data <- cluster_data[,c('x', 'y', 'new_cluster', 'color')]
  cluster_data[,'color'] <- as.character(cluster_data[,'color'])
  
  if (nrow(cluster_data) >= 3) {
    
    hull_idx <- chull(cluster_data[,'x'], cluster_data[,'y'])
    cluster_data <- cluster_data[hull_idx, ]
    
  }
  
  if (nrow(cluster_data) >= 4 & cluster_data[1,'new_cluster'] != 4) {
    
    center_x <- mean(cluster_data[,'x'])
    center_y <- mean(cluster_data[,'y'])
    
    cluster_data[,'x_diff'] <- cluster_data[,'x'] - center_x
    cluster_data[,'y_diff'] <- cluster_data[,'y'] - center_y
    cluster_data[,'x'] <- center_x + cluster_data[,'x_diff'] + 0.03 * abs(diff(range(df[,'x']))) * sign(cluster_data[,'x_diff'])
    cluster_data[,'y'] <- center_y + cluster_data[,'y_diff'] + 0.03 * abs(diff(range(df[,'y']))) * sign(cluster_data[,'y_diff'])
    
    return(cluster_data[,c('x', 'y', 'color')])
    
  } else {
    
    center_x <- mean(cluster_data[,'x'])
    center_y <- mean(cluster_data[,'y'])
    
    cluster_data.1 <- cluster_data[1:2,]
    cluster_data.2 <- cluster_data[1:2,]
    
    cluster_data.1[,'x'] <- center_x + c(0.02, 0.02) * (abs(diff(range(cluster_data[,'x']))) * 25 + abs(diff(range(df[,'x']))))
    cluster_data.1[,'y'] <- center_y + c(0.02, -0.02) * (abs(diff(range(cluster_data[,'y']))) * 25 + abs(diff(range(df[,'y']))))
    
    cluster_data.2[,'x'] <- center_x + c(-0.02, -0.02) * (abs(diff(range(cluster_data[,'x']))) * 25 + abs(diff(range(df[,'x']))))
    cluster_data.2[,'y'] <- center_y + c(-0.02, 0.02) * (abs(diff(range(cluster_data[,'y']))) * 25 + abs(diff(range(df[,'y']))))
    
    cluster_data <- rbind(cluster_data.1, cluster_data.2)
    
    return(cluster_data[,c('x', 'y', 'color')])
    
  }
  
})

## 畫圖

dr_p <- ggplot(df, aes(x = x, y = y))

for (i in 1:length(hull_list)) {
  
  dr_p <- dr_p + annotate(geom = "polygon",
                          x = hull_list[[i]][,'x'],
                          y = hull_list[[i]][,'y'],
                          fill = hull_list[[i]][1,'color'],
                          alpha = 0.2, color = NA)
  
  dr_p <- dr_p + annotate(geom = "text",
                          x = mean(hull_list[[i]][,'x']),
                          y = mean(hull_list[[i]][,'y']),
                          label = i,
                          color = hull_list[[i]][1,'color'],
                          alpha = 0.8, size = 4.5, fontface = 2)
  
}

dr_p <- dr_p + geom_point(aes(shape = label_type, color = cluster_group), size = 2, alpha = 0.5)
dr_p <- dr_p + theme_minimal()

dr_p <- dr_p + scale_color_manual(values = levels(df[,'color']))
dr_p <- dr_p + labs(color = "Cluster", shape = "Label type", fill = "Cluster")

should_be_highlighted <- subset(df, new_cluster <= 0 | new_cluster > ignore_cluster)

name_revise_pos <- which(should_be_highlighted[,'new_cluster'] > ignore_cluster)
should_be_highlighted[name_revise_pos,'phenotype'] <- cluster_names[should_be_highlighted[name_revise_pos,'new_cluster']]

dr_p <- dr_p +  geom_text_repel(data = should_be_highlighted,
                                aes(x = x, y = y, label = phenotype), 
                                nudge_y = 0.01,
                                nudge_x = 0.01,
                                size = 2.5,
                                box.padding = 0.3,
                                point.padding = 0.5,
                                colour = should_be_highlighted[,'color'],
                                segment.color = paste0(should_be_highlighted[,'color'], '60'),
                                segment.size  = 0.2,
                                min.segment.length = 0,
                                arrow= arrow(length = unit(0.005, "npc")))

dr_p <- dr_p + theme(plot.title = element_blank(),
                     legend.position = "right",
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank()) 


# 6. Merge plot

if (!dir.exists(dirname(plot_path))) {dir.create(dirname(plot_path), recursive = TRUE)}

pdf(plot_path, width = 12, height = 8)

print(dr_p)

dev.off()


