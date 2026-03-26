library(dplyr)
library(ggplot2)
library(reshape2)

load('baseline.RData')
selected_model_info = read.csv('selected_model_info.csv')




# 1. 102 phecodes（binary + survival）

selected_labels <- selected_model_info$label
selected_labels2 <- as.character(selected_labels)


mat_102 <- CXR_label[, selected_labels2, drop = FALSE]


mat_102 <- data.frame(lapply(mat_102, function(x){
  if(is.factor(x)) x <- as.character(x)
  as.numeric(x)
}))


keep <- sapply(mat_102, function(x){
  ux <- unique(x[!is.na(x)])
  length(ux) > 1
})

mat_102 <- mat_102[, keep, drop = FALSE]

cat("Remaining phecodes:", ncol(mat_102), "\n")

#========================
# 4. correlation（EHR）
#========================
cor_mat <- cor(mat_102, use = "pairwise.complete.obs")


cor_mat[is.na(cor_mat)] <- 0
diag(cor_mat) <- 1


hc <- hclust(as.dist(1 - cor_mat), method = "average")
ord <- hc$order

cor_mat <- cor_mat[ord, ord]


info_ord <- selected_model_info[ord, ]

display_names <- info_ord$PhecodeString
display_names[is.na(display_names) | display_names == ""] <- info_ord$label

display_names <- paste0(
  display_names,
  " (",
  ifelse(info_ord$label_type == "binary", "prevalent", "incident"),
  ")"
)

rownames(cor_mat) <- display_names
colnames(cor_mat) <- display_names


cor_df <- melt(cor_mat)
colnames(cor_df) <- c("Var1", "Var2", "Correlation")

p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Phi\ncorrelation"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) 

ggsave("Fig S03.png", p, width = 10, height = 9, dpi = 300)
ggsave("Fig S03.pdf", p, width = 10, height = 9)

p