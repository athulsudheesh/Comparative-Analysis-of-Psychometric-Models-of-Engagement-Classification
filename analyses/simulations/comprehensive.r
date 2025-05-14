library(tinytable)
C10matrix <- readRDS("analyses/simulations/C_10/confusion_matrix.rds")
numeric_cols <- names(C10matrix)[names(C10matrix) != "Model"]
C10matrix[numeric_cols] <- C10matrix[numeric_cols] * 100
rownames(C10matrix) <-  NULL
C10matrix <- C10matrix[c('Model','Accuracy', 'Sensitivity', 'Specificity')]
tt(C10matrix, caption = "Engagement classification performance for simulation condition where probability of engagement is 10%")|>
            format_tt(j = 2:4, digits = 2,num_fmt = "decimal")  |>
            print("latex")


C30matrix <- readRDS("analyses/simulations/C_30/confusion_matrix.rds")
C30matrix <- C30matrix[c('Model','Accuracy', 'Sensitivity', 'Specificity')]
numeric_cols <- names(C30matrix)[names(C30matrix) != "Model"]
C30matrix[numeric_cols] <- C30matrix[numeric_cols] * 100
rownames(C30matrix) <-  NULL

tt(C30matrix, caption = "Engagement classification performance for simulation condition where probability of engagement is 30%")|>
            format_tt(j = 2:4, digits = 2,num_fmt = "decimal")  |>
            print("latex")


C60matrix <- readRDS("analyses/simulations/C_60/confusion_matrix.rds")
C60matrix <- C60matrix[c('Model','Accuracy', 'Sensitivity', 'Specificity')]
numeric_cols <- names(C60matrix)[names(C60matrix) != "Model"]
C60matrix[numeric_cols] <- C60matrix[numeric_cols] * 100
rownames(C60matrix) <-  NULL

tt(C60matrix, caption = "Engagement classification performance for simulation condition where probability of engagement is 60%")|>
            format_tt(j = 2:4, digits = 2,num_fmt = "decimal")  |>
            print("latex")


C90matrix <- readRDS("analyses/simulations/C_90/confusion_matrix.rds")
C90matrix <- C60matrix[c('Model','Accuracy', 'Sensitivity', 'Specificity')]
numeric_cols <- names(C90matrix)[names(C90matrix) != "Model"]
C90matrix[numeric_cols] <- C90matrix[numeric_cols] * 100
rownames(C90matrix) <-  NULL

tt(C90matrix, caption = "Engagement classification performance for simulation condition where probability of engagement is 90%")|>
            format_tt(j = 2:4, digits = 2,num_fmt = "decimal")  |>
            print("latex")

# Combine the dataframes and add a condition column
C10matrix$Condition <- "C10"
C30matrix$Condition <- "C30"
C60matrix$Condition <- "C60"
C90matrix$Condition <- "C90"

# Combine all matrices into one dataframe
combined_data <- rbind(C10matrix, C30matrix, C60matrix, C90matrix)

# Convert Condition to a factor with proper ordering
combined_data$Condition <- factor(combined_data$Condition, 
                                 levels = c("C10", "C30", "C60", "C90"))


library(ggplot2)
library(viridis)  # For a better color palette

# Create the plot
a<- ggplot(combined_data, aes(x = Condition, y = Accuracy, 
                         group = Model, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Accuracy (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Common-k" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),

    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
  )


b<- ggplot(combined_data, aes(x = Condition, y = Specificity, 
                         group = Model, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Specificity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Common-k" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
  )


c <- ggplot(combined_data, aes(x = Condition, y = Sensitivity, 
                         group = Model, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Sensitivity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Common-k" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
  )


d<- ggplot(combined_data, aes(x = Condition, y = Precision, 
                         group = Model, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Precision (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Common-k" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),

    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
  )
library(ggpubr)
p1 <- a + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))+  coord_fixed(0.015)
p2 <- b + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))+  coord_fixed(0.015)
p3 <- c + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))+ coord_fixed(0.0166)
p4 <- d + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))+ coord_fixed(0.0165)
ggarrange(p1,p2,p3,nrow=3, common.legend=TRUE,labels=c("A","B", "C","D"), hjust=-1, legend="right")
ggsave("NewConfusionMatrix.pdf")
library(knitr)
plot_crop("ConfusionMatrix.pdf")


model_fit_results <- readRDS("analyses/simulations/C_10/model_fit_results.rds")

library(tinytable)
model_fit_results$SE <- formatC(model_fit_results$SE,format = "e", digits = 2)
tt(model_fit_results,
caption = "Model fit results from MCMC analysis for all models applied to simulated data where engagement percentage is 10") |>
            group_tt(j = list(
                "3-fold CV Error" = 5:6
            )) |>
            format_tt(j = 2:6, digits = 2,num_fmt = "decimal") |>
            style_tt(i = 7, bold=TRUE) |>
            print("latex")


model_fit_results <- readRDS("analyses/simulations/C_30/model_fit_results.rds")

library(tinytable)
model_fit_results$SE <- formatC(model_fit_results$SE,format = "e", digits = 2)
tt(model_fit_results,
caption = "Model fit results from MCMC analysis for all models applied to simulated data where engagement percentage is 30") |>
            group_tt(j = list(
                "3-fold CV Error" = 5:6
            )) |>
            format_tt(j = 2:6, digits = 2,num_fmt = "decimal") |>
            style_tt(i = 7, bold=TRUE) |>
            print("latex")

model_fit_results <- readRDS("analyses/simulations/C_60/model_fit_results.rds")

library(tinytable)
model_fit_results$SE <- formatC(model_fit_results$SE,format = "e", digits = 2)
tt(model_fit_results,
caption = "Model fit results from MCMC analysis for all models applied to simulated data where engagement percentage is 60") |>
            group_tt(j = list(
                "3-fold CV Error" = 5:6
            )) |>
            format_tt(j = 2:6, digits = 2,num_fmt = "decimal") |>
            style_tt(i = 7, bold=TRUE) |>
            print("latex")



model_fit_results <- readRDS("analyses/simulations/C_90/model_fit_results.rds")

library(tinytable)
model_fit_results$SE <- formatC(model_fit_results$SE,format = "e", digits = 2)
tt(model_fit_results,
caption = "Model fit results from MCMC analysis for all models applied to simulated data where engagement percentage is 90") |>
            group_tt(j = list(
                "3-fold CV Error" = 5:6
            )) |>
            format_tt(j = 2:6, digits = 2,num_fmt = "decimal") |>
            style_tt(i = 7, bold=TRUE) |>
            print("latex")

### HIGHLISGHTS PREDICTIVE

# Create the plot
a <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Accuracy, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Accuracy, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Accuracy, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Accuracy, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Accuracy (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Common-k" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )


b <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Specificity, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Specificity, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Specificity, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Specificity, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Specificity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Common-k" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )


c <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Sensitivity, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Sensitivity, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Sensitivity, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("DLCSL", "DLCTL", "GoMRT")),
            aes(x = Condition, y = Sensitivity, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Sensitivity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Common-k" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )


library(ggpubr)
p1 <- a + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))+  coord_fixed(0.015)
p2 <- b + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))+  coord_fixed(0.015)
p3 <- c + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))+ coord_fixed(0.0166)
p4 <- d + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))+ coord_fixed(0.0165)
predictive <- ggarrange(p1,p2,p3,nrow=3, common.legend=TRUE, hjust=-1, legend="bottom")
ggsave("NewPredictiveConfusionMatrix.pdf")
library(knitr)
plot_crop("NewPredictiveConfusionMatrix.pdf")


### HIGHLISGHTS CAUSAL

# Create the plot
a <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Accuracy, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Accuracy, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Accuracy, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Accuracy, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Accuracy (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )


b <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Specificity, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Specificity, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Specificity, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Specificity, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Specificity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )


c <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Sensitivity, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Sensitivity, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Sensitivity, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("MHM", "ILCRE", "ILCRI")),
            aes(x = Condition, y = Sensitivity, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Sensitivity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )
p1 <- a + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))+  coord_fixed(0.015)
p2 <- b + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))+  coord_fixed(0.015)
p3 <- c + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))+ coord_fixed(0.0166)
p4 <- d + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))+ coord_fixed(0.0165)
causal <- ggarrange(p1,p2,p3,nrow=3, common.legend=TRUE, hjust=-1, legend="bottom")
ggsave("CausalConfusionMatrix.pdf")
library(knitr)
plot_crop("CausalConfusionMatrix.pdf")


#### Threshold Highlight

# Create the plot
a <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Accuracy, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Accuracy, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Accuracy, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Accuracy, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Accuracy (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Normative" = "#cab2d6",
    "Common-k" = "#6a3d9a",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )


b <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Specificity, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Specificity, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Specificity, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Specificity, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Specificity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Normative" = "#cab2d6",
    "Common-k" = "#6a3d9a",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )


c <- ggplot() +
  # Add gray lines for non-highlighted models
  geom_line(data = subset(combined_data, !Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Sensitivity, group = Model), 
            color = "#eeeeee", linewidth = 1) +
  geom_point(data = subset(combined_data, !Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Sensitivity, group = Model), 
            color = "#eeeeee", size = 3) +
  # Add colored lines for the highlighted models
  geom_line(data = subset(combined_data, Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Sensitivity, group = Model, color = Model), linewidth = 1) +
  geom_point(data = subset(combined_data, Model %in% c("Normative", "Common-k", "VICS", "VII")),
            aes(x = Condition, y = Sensitivity, group = Model, color = Model), size = 3) +
  theme_minimal() +
  labs(
       x = "Condition",
       y = "Sensitivity (%)",
       color = "Model") +
  scale_color_manual(values = c(
    "Normative" = "#cab2d6",
    "Common-k" = "#6a3d9a",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF"
  )) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() + 
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.border = element_rect(size = 1.5),
    axis.title = element_text(face = "bold")
  )
p1 <- a + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))+ coord_fixed()
p2 <- b + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 15,  # Bottom margin
                             l = 5))
p3 <- c + theme(plot.margin = margin(t = 5,  # Top margin
                             r = 5,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5))
threshold <- ggarrange(p1,p2,p3,nrow=3, common.legend=TRUE, hjust=-1, legend="bottom")
ggarrange(threshold, predictive, causal, ncol=3)

ggsave("ThresholdConfusionMatrix.pdf")
library(knitr)
plot_crop("ThresholdConfusionMatrix.pdf")