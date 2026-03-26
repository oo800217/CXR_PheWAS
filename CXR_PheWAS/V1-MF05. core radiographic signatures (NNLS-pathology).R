library(scales) 
library(ggrepel)
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(Rtsne)
library(dplyr)
library(RColorBrewer)

# =========================
# 0. settings
# =========================

selected_clusters <- c(1,2,3,4,6,7,8,9,11,12)

cluster_names_full <- c(
  '1. Future cardiovascular disease',
  '2. Future critical illness',
  '3. Cachexia and advanced malignancies',
  '4. Advanced heart failure spectrum',
  '6. Sepsis and infarction',
  '7. Cardiac arrest and respiratory failure',
  '8. Dialysis-dependent disease',
  '9. Advanced chronic kidney disease',
  '11. Osteoporosis',
  '12. Prostate disorder'
)

selected_features <- c(
  "Atherosclerosis",
  "Ground glass opacity",
  "Emphysematous change",
  "Atelectasis",
  "Cardiomegaly",
  "Nodule",
  "Osteophyte formation"
)

radio_order_path <- 'data/radio_order.RData'
data_path <- 'data/explainable_aucs (NNLS-pathology).RData'
plot_path <- 'result/Fig 05.pdf'

rsq_name <- 'rsq_less'

# =========================
# 1. Load data
# =========================

load(data_path)
load(radio_order_path)

radio_list <- radio_list.less

# =========================
# 2. Build coefficient matrix
# =========================

radio_coef <- do.call('cbind', radio_list)
colnames(radio_coef) <- names(radio_list)

radio_coef <- radio_coef * matrix(
  sd_list.radio,
  nrow = nrow(radio_coef),
  ncol = ncol(radio_coef),
  byrow = FALSE
)

radio_coef <- radio_coef / t(matrix(
  sd_list.phe,
  nrow = ncol(radio_coef),
  ncol = nrow(radio_coef),
  byrow = FALSE
))

radio_coef[radio_coef > 1] <- 1
radio_coef[radio_coef <= 0] <- 0

# =========================
# 3. Aggregate to cluster level
# =========================

new_radio_coef <- matrix(
  NA,
  nrow = nrow(radio_coef),
  ncol = length(selected_clusters)
)

rownames(new_radio_coef) <- rownames(radio_coef)
colnames(new_radio_coef) <- cluster_names_full

for (i in seq_along(selected_clusters)) {
  
  cl <- selected_clusters[i]
  
  phe <- AUC_data[!duplicated(AUC_data$y) & AUC_data$cluster == cl, 'y']
  phe <- phe[phe %in% colnames(radio_coef)]
  
  if (length(phe) > 0) {
    new_radio_coef[, i] <- rowMeans(radio_coef[, phe, drop = FALSE], na.rm = TRUE)
  }
}

# =========================
# 4. Feature filtering
# =========================

feature_df <- data.frame(
  raw = rownames(new_radio_coef),
  lower = tolower(rownames(new_radio_coef)),
  stringsAsFactors = FALSE
)

feature_df$clean <- feature_df$lower

feature_df$clean[feature_df$lower == "atherosclerosis"] <- "Atherosclerosis"
feature_df$clean[feature_df$lower == "ground glass opacity"] <- "Ground glass opacity"
feature_df$clean[feature_df$lower == "emphysematous change"] <- "Emphysematous change"
feature_df$clean[feature_df$lower == "atelectasis"] <- "Atelectasis"
feature_df$clean[feature_df$lower == "cardiomegaly"] <- "Cardiomegaly"
feature_df$clean[feature_df$lower == "nodule"] <- "Nodule"
feature_df$clean[feature_df$lower == "osteophyte formation"] <- "Osteophyte formation"
feature_df$clean[feature_df$lower == "aneurysm dilatation"] <- "Aneurysmal dilatation"

keep <- feature_df$clean %in% selected_features

new_radio_coef <- new_radio_coef[keep, , drop = FALSE]
feature_df <- feature_df[keep, , drop = FALSE]

rownames(new_radio_coef) <- feature_df$clean

# 依指定順序排列
new_radio_coef <- new_radio_coef[selected_features, , drop = FALSE]

# =========================
# 5. Mean R-square
# =========================

sub_auc <- AUC_data[AUC_data$cluster %in% selected_clusters, ]

mean_auc <- tapply(
  sub_auc[, rsq_name],
  sub_auc[, c('cluster', 'dataset')],
  mean
) %>% t()

auc_df <- melt(mean_auc)
colnames(auc_df) <- c("Var1", "Var2", "value")

cluster_map <- data.frame(
  cluster = selected_clusters,
  name = cluster_names_full,
  stringsAsFactors = FALSE
)

auc_df$Var2 <- cluster_map$name[match(as.numeric(as.character(auc_df$Var2)), cluster_map$cluster)]

auc_df$Var1 <- gsub(" dataset", "", auc_df$Var1)
auc_df$Var1 <- paste0("Mean R-square (", auc_df$Var1, ")")

# =========================
# 6. Combine data
# =========================

coef_df <- melt(new_radio_coef)
colnames(coef_df) <- c("Var1", "Var2", "value")

cor_df <- rbind(coef_df, auc_df)

y_levels <- c(
  "Mean R-square (MIMIC)",
  "Mean R-square (External)",
  "Mean R-square (Hold-out)",
  selected_features
)

cor_df$Var1 <- factor(cor_df$Var1, levels = rev(y_levels))
cor_df$Var2 <- factor(cor_df$Var2, levels = cluster_names_full)

# =========================
# 7. Highlight cells
# =========================

highlight_pairs <- data.frame(
  Var1 = c(
    # C2
    "Atherosclerosis","Ground glass opacity","Emphysematous change","Atelectasis",
    # C1
    "Cardiomegaly","Atherosclerosis",
    # C3
    "Ground glass opacity","Atelectasis","Nodule",
    # C11
    "Atherosclerosis","Cardiomegaly",
    # C12
    "Osteophyte formation"
  ),
  Var2 = c(
    cluster_names_full[2],cluster_names_full[2],cluster_names_full[2],cluster_names_full[2],
    cluster_names_full[1],cluster_names_full[1],
    cluster_names_full[3],cluster_names_full[3],cluster_names_full[3],
    cluster_names_full[9],cluster_names_full[9],
    cluster_names_full[10]
  ),
  stringsAsFactors = FALSE
)

highlight_pairs$flag <- TRUE

cor_df <- merge(cor_df, highlight_pairs, all.x = TRUE)
cor_df$flag[is.na(cor_df$flag)] <- FALSE

# =========================
# 8. Plot
# =========================

my_colors <- c(
  "#FFFFFF",
  "#FEE0D2",
  "#FC9272",
  "#FB6A4A",
  "#EF3B2C",
  "#99000D"
)

coef_p <- ggplot(cor_df, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0.005, 1),
    name = "Standardized regression coefficients\nor R-square",
    na.value = "grey55"
  ) +
  geom_text(
    aes(label = ifelse(is.na(value) | value < 0.005, "", sprintf("%.2f", value))),
    size = 2.8,
    color = "black"
  ) +
  geom_tile(
    data = cor_df[cor_df$flag == TRUE, ],
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, size = 7, vjust = 1, hjust = 1, color = 'black'),
    axis.text.y = element_text(size = 7, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")

# =========================
# 9. Top R-square border
# =========================
# 修正白色突出：只框住最上面3列，不超出 panel

n_y <- length(levels(cor_df$Var1))

coef_p <- coef_p + annotate(
  "rect",
  xmin = 0.5,
  xmax = length(cluster_names_full) + 0.5,
  ymin = n_y - 2.5,
  ymax = n_y + 0.5 - 1e-6,
  fill = NA,
  color = "black",
  linewidth = 0.6
)

# =========================
# 10. Save
# =========================

if (!dir.exists(dirname(plot_path))) {
  dir.create(dirname(plot_path), recursive = TRUE)
}

pdf(plot_path, width = 10, height = 8)
print(coef_p)
dev.off()

print(coef_p)