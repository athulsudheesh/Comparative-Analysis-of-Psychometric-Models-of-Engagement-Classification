source("scripts/model-fits.r") 
source("scripts/data_handling_utils.r")
data <- load_and_clean_data(680810)
library(MCMCvis)
header_names <- c("Model", "lppd", "pWAIC", "WAIC", "MSE", "SE")
model_fit_results <- data.frame(matrix(ncol = length(header_names), nrow=0))
colnames(model_fit_results) <- header_names


model <- CommonK.fit(data)
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "Common-k",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/commonk.rds")
#model_fit_results <- model_fit_results[-nrow(model_fit_results), ]
model <- readRDS("analyses/model-fit_results/680810/commonk.rds")

model <- Normative.fit(data)
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "Normative",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/normative.rds")
model <- readRDS("analyses/model-fit_results/680810/normative.rds")

model <- VICS.fit(data)
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "VICS",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/vics.rds")
model <- readRDS("analyses/model-fit_results/680810/vics.rds")

model <- VII.fit(data)
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "VII",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/vii.rds")
model <- readRDS("analyses/model-fit_results/680810/vii.rds")

model <- GoMRT.fit(data, 
                   params.to.monitor = c("correct_resp_prob",
                                          "class_p","C"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "GoMRT",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/gomrt.rds")
model <- readRDS("analyses/model-fit_results/680810/gomrt.rds")

model <- DLCSL.fit(data, 
                   params.to.monitor = c("correct_resp_prob", 
                                          "class_p","C", "tau"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "DLCSL",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/dlcsl.rds")
model <- readRDS("analyses/model-fit_results/680810/dlcsl.rds")

model <- DLCTL.fit(data, 
                   params.to.monitor = c("correct_resp_prob",
                                          "class_p","C","theta","b","zeta","tau"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "DLCTL",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/dlctl.rds")
model <- readRDS("analyses/model-fit_results/680810/dlctl.rds")
# model <- nimbleModel(code = MHM, constants = const, 
#                      data = data, inits = inits(),
#                      check = TRUE, calculate = FALSE)
model <- MHM.fit(data, 
                   params.to.monitor = c("correct_resp_prob", "theta","b",
                                          "class_p","C", "eta"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "MHM",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/mhm.rds")
model <- readRDS("analyses/model-fit_results/680810/mhm.rds")

model <- ILCRE.fit(data, 
                   params.to.monitor = c("correct_resp_prob","theta","b",
                                          "class_p","C", "eta","xi","kappa"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "ILCRE",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/ilcre.rds")
model <- readRDS("analyses/model-fit_results/680810/ilcre.rds")

model <- ILCRI.fit(data, 
                   params.to.monitor = c("correct_resp_prob", "theta","b",
                                          "class_p","C", "eta","xi","kappa"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "ILCRI",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/680810/ilcri.rds")
model <- readRDS("analyses/model-fit_results/680810/ilcri.rds")



saveRDS(model_fit_results,"analyses/model-fit_results/680810/model_fit_results.rds")

#===========================================================================
#                            Model Fit results table
#===========================================================================
model_fit_results <- readRDS("analyses/model-fit_results/680810/model_fit_results.rds")
library(tinytable)
model_fit_results <- model_fit_results[,c("Model","lppd", "pWAIC", "WAIC")]
model_fit_results$SE <- formatC(model_fit_results$SE,format = "e", digits = 2)
tt(model_fit_results, width=1,
caption = "Model fit results from MCMC analysis for all models applied to Assessment 747119 data") |>
            format_tt(j = 2:4, digits = 2,num_fmt = "decimal") |>
            style_tt(i = 7, bold=TRUE) |>
            print("latex")
library(tidyverse)

# Assuming model_fit_results is your data frame
library(tidyverse)
model_fit_results <- model_fit_results |>
  mutate(
    # For lppd, higher is better
    lppd_rank = rank(-lppd),  # Negative to reverse the ordering
    
    # For the other metrics, lower is better
    pWAIC_rank = rank(pWAIC),
    WAIC_rank = rank(WAIC),
    MSE_rank = rank(MSE)
  )
library(tidyverse)
library(ggplot2)

library(tidyverse)
library(ggplot2)

# First, reshape the data to long format but only include WAIC and lppd ranks
model_rankings_long <- model_fit_results %>%
  select(Model, lppd_rank, WAIC_rank) %>%
  pivot_longer(
    cols = c(lppd_rank, WAIC_rank),
    names_to = "Metric",
    values_to = "Rank"
  ) %>%
  # Clean up metric names by removing "_rank" suffix
  mutate(Metric = str_remove(Metric, "_rank"),
         # Rename the metrics for the legend
         Metric = case_when(
           Metric == "WAIC" ~ "out-of-sample performance",
           Metric == "lppd" ~ "in-sample performance",
           TRUE ~ Metric
         ))

# Get the WAIC rank order to use for sorting models
waic_order <- model_fit_results %>%
  arrange(WAIC_rank) %>%
  pull(Model)

# Convert Model to a factor with levels ordered by WAIC rank
model_rankings_long <- model_rankings_long %>%
  mutate(Model = factor(Model, levels = waic_order))

# Create the rank order plot with the specified colors
rank_order <- ggplot(model_rankings_long, aes(x = Model, y = Rank, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal()+
  scale_y_continuous(breaks = 1:10) +
  coord_flip() +  # Flip coordinates for horizontal bars
  labs(title = "680810", 
    x = "Model",
    y = "Rank",
    fill = "Performance Metric"  # This renames the legend title
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  ) +
  # Use the two specific colors you provided
  scale_fill_manual(values = c("in-sample performance" = "#0F425CFF", 
                               "out-of-sample performance" = "#800000FF"))

saveRDS(rank_order, "680810rankorder.rds")
ggsave("rankorderplot.pdf")

#===========================================================================
#                            MCMC Coverage Table 
#===========================================================================
compute_coverage <- function(model,params){
    sum(MCMCsummary(model$posterior_samples$samples, params=params)$Rhat < 1.1, na.rm = TRUE)*100/length(MCMCsummary(model$posterior_samples$samples,params=params)$Rhat)
}

header_names <- c("Model","{\\theta}", "{d}", "{C}", "{\\eta}","{\\kappa}","{\\xi}")
coverage_results <- data.frame(matrix(ncol = length(header_names), nrow=0))
colnames(coverage_results) <- header_names


library(MCMCvis)
model <- readRDS("analyses/model-fit_results/680810/commonk.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "Common-k", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))


model <- readRDS("analyses/model-fit_results/680810/normative.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "Normative", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/680810/vics.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "VICS", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/680810/vii.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "VII", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/680810/gomrt.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "GoMRT", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/680810/dlcsl.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "DLCSL", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/680810/dlctl.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "DLCTL", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/680810/mhm.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "MHM", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/680810/ilcre.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
etacoverage <- compute_coverage(model, "eta")
kappacoverage <- compute_coverage(model, "kappa")
xicoverage <- compute_coverage(model, "xi")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "ILCRE", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = etacoverage, 
                        "kappa" = kappacoverage, "xi"= xicoverage))

model <- readRDS("analyses/model-fit_results/680810/ilcri.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
etacoverage <- compute_coverage(model, "eta")
kappacoverage <- compute_coverage(model, "kappa")
xicoverage <- compute_coverage(model, "xi")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "ILCRI", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = etacoverage, 
                        "kappa" = kappacoverage, "xi"= xicoverage))

saveRDS(coverage_results,"analyses/model-fit_results/680810/coverage_results.rds")
library(tinytable)
coverage_results<- coverage_results |> 
  mutate_at(vars("C", "eta","kappa","xi"), as.numeric)
tt(coverage_results,
    caption = "MCMC parameter coverage based on Rhat statistics for Assessment 680810") |>
     format_tt(j = 3:7, digits = 2,num_fmt = "decimal", replace = TRUE) |>
     style_tt(i =5:7, j = 4, color="red") |>
     style_tt(i =9, j = c(4,7), color="red") |>
     print("latex")

#===========================================================================
#                            Engagement Classification Plots and Tables
#===========================================================================
model <- readRDS("analyses/model-fit_results/680810/vics.rds")
library(coda)
library(codatools)

I = nrow(data$X)
J = ncol(data$X)
smode <- function(x) {
        return(as.numeric(names(which.max(table(x)))))
    }
#C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
#C <- matrix(apply(C ,2,smode),I,J)
RT <- data$RT
RT[is.na(RT)] <- 0
# Common-k Threshold
C_commonk <- RT >5
C_commonk <- +C_commonk
# Normative Threshold 
k = 0.1
thresholds = k* colMeans(RT <- RT, na.rm = TRUE)
C_normative <- t(t(RT) > thresholds)
C_normative <- +C_normative
# VICS Threshold 
model <- readRDS("analyses/model-fit_results/680810/vics.rds")
thresholds = model$thresholds
C_vics <- t(t(RT) > thresholds)
C_vics <- +C_vics

# VII Threshold
model <- readRDS("analyses/model-fit_results/680810/vii.rds")
thresholds = model$thresholds
C_vii <- t(t(RT) > thresholds)
C_vii <- +C_vii
# GoMRT 
library(codatools)
model <- readRDS("analyses/model-fit_results/680810/gomrt.rds")
C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
C_gomrt <- matrix(apply(C ,2,smode),I,J)

model <- readRDS("analyses/model-fit_results/680810/dlcsl.rds")
C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
C_dlcsl <- matrix(apply(C ,2,smode),I,J)

model <- readRDS("analyses/model-fit_results/680810/dlctl.rds")
C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
C_dlctl <- matrix(apply(C ,2,smode),I,J)


model <- readRDS("analyses/model-fit_results/680810/mhm.rds")
C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
C_mhm <- matrix(apply(C ,2,smode),I,J)

model <- readRDS("analyses/model-fit_results/680810/ilcre.rds")
C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
C_ilcre <- matrix(apply(C ,2,smode),I,J)


model <- readRDS("analyses/model-fit_results/680810/ilcri.rds")
C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
C_ilcri <- matrix(apply(C ,2,smode),I,J)

C_raw_disagreement <- C_commonk + C_normative + C_vics + C_vii + C_gomrt + C_dlcsl + C_dlctl + C_mhm + C_ilcre + C_ilcri

library(tidyverse)
library(reshape2)


plot_engagement_disagreement <- function(C_raw_disagreement, color_fill="red", element = element_text()){
    
    I <- nrow(C_raw_disagreement)
    J <- ncol(C_raw_disagreement)
    rownames(C_raw_disagreement) <- paste0("S", 1:I)
    colnames(C_raw_disagreement) <- colnames(C_raw_disagreement)
    melted_mat <- melt(C_raw_disagreement)
    names(melted_mat) <- c("row", "column", "value")
    melted_mat$column <- factor(melted_mat$column, levels = unique(melted_mat$column))

    ggplot(melted_mat, aes(x = column, y = row, fill = value)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_tile(color = "black") +  # add black borders around tiles [[3]]
    scale_fill_gradient(low = "white", high =color_fill) +  # color gradient [[2]]
    theme_bw() +
    #coord_fixed() +  # ensure tiles are square [[2]]
    labs(
        x = "Item ID", 
        y = "Students", 
        fill = "Disagreement")  +
                theme(legend.position = "none",
                        plot.title = element_text(hjust=0.5),
                        axis.text.x = element_text(angle = 90, size = 8),
                        axis.text.y = element_text(),
                        panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold"),
                        axis.title.y = element
                    )
}

## Updated plotfunction 
plot_engagement_disagreement <- function(C_raw_disagreement, X_matrix =data$X, color_fill="red", element = element_text()){
    
    I <- nrow(C_raw_disagreement)
    J <- ncol(C_raw_disagreement)
    
    # Calculate student scores and item difficulties if X_matrix is provided
    if(!is.null(X_matrix)) {
        student_scores <- rowMeans(X_matrix, na.rm = TRUE)
        item_difficulties <- colMeans(X_matrix, na.rm = TRUE)
        
        # Create student labels with scores
        student_labels <- paste0("S", 1:I, " (", round(student_scores, 2), ")")
        
        # Order the rows of C_raw_disagreement by student scores (descending)
        score_order <- order(student_scores, decreasing = TRUE)
        C_raw_disagreement <- C_raw_disagreement[score_order, ]
        rownames(C_raw_disagreement) <- student_labels[score_order]
        
        # Order the columns of C_raw_disagreement by item difficulty (ascending)
        difficulty_order <- order(item_difficulties)
        C_raw_disagreement <- C_raw_disagreement[, difficulty_order]
        colnames(C_raw_disagreement) <- colnames(X_matrix)[difficulty_order]
    } else {
        rownames(C_raw_disagreement) <- paste0("S", 1:I)
        colnames(C_raw_disagreement) <- colnames(C_raw_disagreement)
    }
    
    melted_mat <- melt(C_raw_disagreement)
    names(melted_mat) <- c("row", "column", "value")
    melted_mat$column <- factor(melted_mat$column, levels = unique(melted_mat$column))
    melted_mat$row <- factor(melted_mat$row, levels = rev(unique(melted_mat$row)))  # Reverse to have highest scores at top

    ggplot(melted_mat, aes(x = column, y = row, fill = value)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_tile(color = "black") +  # add black borders around tiles
    scale_fill_gradient(low = "white", high = color_fill) +  # color gradient
    theme_bw() +
    labs(
        x = "Item ID (ordered by difficulty)", 
        y = "Students (ordered by score)", 
        fill = "Disagreement") +
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(angle = 90, size = 8),
          axis.text.y = element_text(),
          panel.border = element_rect(size = 1.5),
          axis.title = element_text(face = "bold"),
          axis.title.y = element_text(size=8),
          axis.title.x = element_text(size=8)
    )
}

plot_engagement_disagreement(C_dlctl, color="#155F83FF", element = element_blank()) + theme(aspect.ratio = 2)
ggsave("680810dlctl_engagement_matrix.pdf")
full<- plot_engagement_disagreement(C_raw_disagreement, color="#155F83FF", element = element_blank()) + theme(aspect.ratio = 2)
C_threshold_disagreement <- C_commonk + C_normative + C_vics + C_vii
threshold_plots <- plot_engagement_disagreement(C_threshold_disagreement,color_fill = "#155F83FF") + theme(aspect.ratio =2.1)

C_latent_diagreement <- C_gomrt+ C_dlcsl + C_dlctl + C_mhm + C_ilcre + C_ilcri
latent_plots<- plot_engagement_disagreement(C_latent_diagreement,color_fill = "#155F83FF",
                        element = element_blank()) + 
    scale_x_discrete(labels = colnames(C_raw_disagreement)) + theme(aspect.ratio = 2)
library(knitr)
library(ggpubr)

ggsave("disagreement_ordered.pdf")
plot_crop("disagreement.pdf")

ggarrange(threshold_plots,latent_plots, full, ncol = 3, labels = c("A","B","C"),
 font.label = list(size = 12, face = "bold"), hjust = -1.5, # Move labels right
vjust = 4)
ggplot()


full <- plot_engagement_disagreement(C_raw_disagreement, color="#155F83FF", element = element_blank()) + 
  theme(
    aspect.ratio = 2,
    axis.text.x = element_text(angle = 90, size = 6, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6)
  )

C_threshold_disagreement <- C_commonk + C_normative + C_vics + C_vii
threshold_plots <- plot_engagement_disagreement(C_threshold_disagreement, color_fill = "#155F83FF") + 
  theme(
    aspect.ratio = 2.1,
    axis.text.x = element_text(angle = 90, size = 6, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6)
  )

C_latent_diagreement <- C_gomrt + C_dlcsl + C_dlctl + C_mhm + C_ilcre + C_ilcri
latent_plots <- plot_engagement_disagreement(C_latent_diagreement, color_fill = "#155F83FF",
                        element = element_blank()) + 
    scale_x_discrete(labels = colnames(C_raw_disagreement)) + 
    theme(
      aspect.ratio = 2,
      axis.text.x = element_text(angle = 90, size = 6, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 6)
    )

# Add more spacing between plots and adjust label position
ggarrange(threshold_plots, latent_plots, full, 
          ncol = 3, 
          labels = c("A", "B", "C"),
          font.label = list(size = 12, face = "bold"), 
          hjust = -0.5,  # Less negative hjust value to move labels right
          vjust = 2,     # Reduced vjust to move labels down
          widths = c(1, 1, 1),  # Equal widths for all plots
          align = "h"    # Horizontal alignment
)

full<- plot_engagement_disagreement(C_raw_disagreement, color="#155F83FF", element = element_blank()) + theme(aspect.ratio = 2.1)
C_threshold_disagreement <- C_commonk + C_normative + C_vics + C_vii
threshold_plots <- plot_engagement_disagreement(C_threshold_disagreement,color_fill = "#155F83FF") + theme(aspect.ratio =2.1)

C_latent_diagreement <- C_gomrt+ C_dlcsl + C_dlctl + C_mhm + C_ilcre + C_ilcri
latent_plots<- plot_engagement_disagreement(C_latent_diagreement,color_fill = "#155F83FF",
                        element = element_blank()) + 
    scale_x_discrete(labels = colnames(C_raw_disagreement)) + theme(aspect.ratio = 2.1)
library(knitr)
library(ggpubr)

ggsave("680810disagreement.pdf")
plot_crop("680810disagreement.pdf")

ggarrange(threshold_plots,latent_plots, full, ncol = 3, labels = c("A","B","C"),
 font.label = list(size = 12, face = "bold"), hjust = -1.5, vjust = 12)
ggplot()

commonk_long <- as.numeric(C_commonk)
normative_long <- as.numeric(C_normative)
vics_long <- as.numeric(C_vics)
vii_long <- as.numeric(C_vii)
gomrt_long <- as.numeric(C_gomrt)
dlcsl_long <- as.numeric(C_dlcsl)
dlctl_long <- as.numeric(C_dlctl)
mhm_long <- as.numeric(C_mhm)
ilcre_long <- as.numeric(C_ilcre)
ilcri_long <- as.numeric(C_ilcri)

ratings_data <- cbind(commonk_long,normative_long, vics_long, vii_long, gomrt_long, dlcsl_long, dlctl_long, mhm_long, ilcre_long,ilcre_long)
colnames(ratings_data) <- c("Common K","Normative", "VICS", "VII", "GoMRT", "DLCSL","DLCTL","MHM","ILCRE", "ILCRI")
n_raters <- ncol(ratings_data)

kappa_matrix <- matrix(NA, nrow = n_raters, ncol = n_raters)
rownames(kappa_matrix) <- colnames(ratings_data)
colnames(kappa_matrix) <- colnames(ratings_data)
library(irr)
for(i in 1:(n_raters-1)) {
    for(j in (i+1):n_raters) {
      # Extract data for the pair of raters
      pair_data <- ratings_data[, c(i, j)]
      
      # Remove rows with NA values
      pair_data <- na.omit(pair_data)
      
      # Calculate kappa if there are enough observations
      if(nrow(pair_data) > 1) {
        k <- kappa2(pair_data)
        kappa_matrix[i, j] <- k$value
        kappa_matrix[j, i] <- k$value  # Matrix is symmetric
      }
    }
  }
diag(kappa_matrix) <- 1
library(tinytable)
as.data.frame(kappa_matrix)
heatmap(kappa_matrix, 
        Rowv = NA, Colv = NA, 
        col = colorRampPalette(c("white", "steelblue"))(100),
        main = "Cohen's Kappa Between Rater Pairs",
        symm = TRUE)+theme_bw()

melted_mat <- melt(kappa_matrix)
    names(melted_mat) <- c("row", "column", "value")
    melted_mat$column <- factor(melted_mat$column, levels = unique(melted_mat$column))

    ggplot(melted_mat, aes(x = column, y = row, fill = value)) + geom_tile() + theme_minimal()+
  scale_fill_gradient(low="white", high="#C16622FF") + geom_text(aes(label = sprintf("%.2f",value)),  # Format to 2 decimal places
            size = 5.5,                           # Adjust text size as needed
            color = "black")+ theme(aspect.ratio = 1,
                    axis.text.x = element_text(size=18, angle =90),
                    axis.text.y = element_text(size=18)) +
            labs(x="", y = "")
ggsave("680810kappa_agreement.pdf")
plot_crop("680810kappa_agreement.pdf")
kappa_matrix
ratings_data
kappa2(ratings_data[,c(1,10)])

X <- data$X

#===========================================================================
#                            NEW KAPPA MATRIX CODE:WITH SIGNIFICANCE
#===========================================================================

# Create a p-value matrix
p_value_matrix <- matrix(NA, nrow = n_raters, ncol = n_raters)
rownames(p_value_matrix) <- colnames(ratings_data)
colnames(p_value_matrix) <- colnames(ratings_data)
diag(p_value_matrix) <- 0  # Diagonal p-values are 0 (self-comparison)

# Recalculate the p-values using the same loop structure as your kappa calculations
for(i in 1:(n_raters-1)) {
  for(j in (i+1):n_raters) {
    # Extract data for the pair of raters
    pair_data <- ratings_data[, c(i, j)]
    
    # Remove rows with NA values
    pair_data <- na.omit(pair_data)
    
    # Calculate kappa if there are enough observations
    if(nrow(pair_data) > 1) {
      k <- kappa2(pair_data)
      p_value_matrix[i, j] <- k$p.value
      p_value_matrix[j, i] <- k$p.value  # Matrix is symmetric
    }
  }
}

# Create significance matrix (TRUE if significant at alpha = 0.05)
significance_matrix <- p_value_matrix < 0.05
diag(significance_matrix) <- TRUE  # Diagonal is always significant (self-comparison)

# Load necessary libraries
library(reshape2)
library(ggplot2)

# Convert matrices to long format for ggplot
kappa_df <- melt(kappa_matrix)
colnames(kappa_df) <- c("Rater1", "Rater2", "Kappa")

# Add significance information
significance_df <- melt(significance_matrix)
colnames(significance_df) <- c("Rater1", "Rater2", "Significant")

# Add p-values
pvalue_df <- melt(p_value_matrix)
colnames(pvalue_df) <- c("Rater1", "Rater2", "Pvalue")

# Combine all information
kappa_df$Significant <- significance_df$Significant
kappa_df$Pvalue <- pvalue_df$Pvalue

# Add significance symbol
kappa_df$SignificanceSymbol <- ifelse(kappa_df$Significant, "*", "")

# Handle NA values for display
kappa_df$KappaText <- ifelse(is.na(kappa_df$Kappa), 
                            "NA", 
                            sprintf("%.2f", kappa_df$Kappa))

# Create a custom color variable based on your logic
kappa_df$CustomColor <- ifelse(kappa_df$Kappa == 1.0, 
                              "#006d2c",  # Dark orange for perfect agreement
                              kappa_df$CustomColor)

kappa_df$CustomColor <- ifelse(!is.na(kappa_df$Kappa) & 
                              kappa_df$Kappa > 0.4 & 
                              kappa_df$Kappa < 1.0 & 
                              kappa_df$Significant == TRUE,
                              "#31a354",  # Lighter orange for kappa between 0.4 and 1.0
                              kappa_df$CustomColor)

kappa_df$CustomColor <- ifelse(!is.na(kappa_df$Kappa) & 
                              kappa_df$Kappa > 0.2 & 
                              kappa_df$Kappa < 0.4 & 
                              kappa_df$Significant == TRUE,
                              "#bae4b3",  # Lighter orange for kappa between 0.4 and 1.0
                              kappa_df$CustomColor)

kappa_df$Category <- "Not significant or low agreement"

# Assign categories based on your existing color logic
kappa_df$Category <- ifelse(!is.na(kappa_df$Kappa) & 
                           kappa_df$Kappa > 0.2 & 
                           kappa_df$Kappa < 0.4 & 
                           kappa_df$Significant == TRUE,
                           "Fair significant agreement (0.2 < k < 0.4)", 
                           kappa_df$Category)

kappa_df$Category <- ifelse(!is.na(kappa_df$Kappa) & 
                           kappa_df$Kappa > 0.4 & 
                           kappa_df$Kappa < 1.0 & 
                           kappa_df$Significant == TRUE,
                           "Moderate-to-high significant agreement (0.4 < k < 1.0)", 
                           kappa_df$Category)

kappa_df$Category <- ifelse(kappa_df$Kappa == 1.0, 
                          "Perfect agreement (k = 1.0)", 
                          kappa_df$Category)

# Create a factor with levels in the desired order
kappa_df$Category <- factor(kappa_df$Category, 
                           levels = c("Perfect agreement (k = 1.0)",
                                     "Moderate-to-high significant agreement (0.4 < k < 1.0)",
                                     "Fair significant agreement (0.2 < k < 0.4)",
                                     "Not significant or low agreement"))

# Create a data frame for the legend
legend_df <- data.frame(
  Category = levels(kappa_df$Category),
  stringsAsFactors = FALSE
)
library(latex2exp)
TeX(r"(Perfect agreement (\kappa = 1.0))")
# Create the plot with the manual legend
kappa_heatmap <- ggplot(kappa_df, aes(x = Rater2, y = Rater1)) +
  geom_tile(aes(fill = Category)) +
  geom_text(aes(label = paste0(KappaText, SignificanceSymbol)), 
            color = "black", size = 3) +
  scale_fill_manual(values = c("Perfect agreement (k = 1.0)" = "#006d2c",
                              "Moderate-to-high significant agreement (0.4 < k < 1.0)" = "#31a354",
                              "Fair significant agreement (0.2 < k < 0.4)" = "#bae4b3",
                              "Not significant or low agreement" = "white"),
                    name = "Agreement Levels")  +
  theme_minimal() +
  theme(aspect.ratio = 1,
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 0, vjust = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8)
  ) + labs(x="",y="")

library(knitr)
# Optional: Save the plot
 ggsave("Updated_kappa_heatmap680810.pdf", kappa_heatmap, width = 15, height = 8)
 plot_crop("Updated_kappa_heatmap680810.pdf")


plot_proportion_correct <- function(X,C_state){
  X[is.na(X)] <- 0
  P <- rowMeans(X)
  E<- rowMeans(X * C_state)
  D <- rowMeans(X* (1-C_state))
  df<- data.frame("E"=E,"D"=D, "P" = P)


  df_long <- df %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(cols = c(E, D), names_to = "Type", values_to = "Value")
  # Create the stacked bar plot with bars ordered by P
  ggplot(df_long, aes(x = factor(row_id), y = Value, fill = Type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("E" = "#CC8214FF", "D" = "#800000FF")) +
    labs(
      x = "Students",
      y = "Proportion Correct",
      fill = "Proportion Correct Given the Engagement Status"
    ) +
    theme_minimal() +
    theme(aspect.ratio = 0.5,
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "top", 
      legend.title =  element_text(size = 15),
      legend.text = element_text(size=15)
    ) 
}


p_commonk<- plot_proportion_correct(X, C_commonk)+ labs(title="Common K") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                       axis.title.x = element_blank()) 
p_normative<- plot_proportion_correct(X, C_normative)+ labs(title="Normative") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                      axis.title.x = element_blank(),  )
p_vics<- plot_proportion_correct(X, C_vics)+ labs(title="VICS") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                      axis.title.x = element_blank()) 
p_vii<- plot_proportion_correct(X, C_vii)+ labs(title="VII") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                      axis.title.x = element_blank()) 
p_gomrt<- plot_proportion_correct(X, C_gomrt)+ labs(title="GoMRT") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                      axis.title.x = element_blank()) 
p_dlcsl<- plot_proportion_correct(X, C_dlcsl)+ labs(title="DLCSL") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                      axis.title.x = element_blank()) 
p_dlctl<- plot_proportion_correct(X, C_dlctl)+ labs(title="DLCTL") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                      axis.title.x = element_blank()) 
p_mhm<- plot_proportion_correct(X, C_mhm)+ labs(title="MHM") + theme(plot.title = element_text(hjust=0.5, size = 16),
                                                                                      axis.title.x = element_blank())
p_ilcre<- plot_proportion_correct(X, C_ilcre)+ labs(title="ILCRE") + theme(plot.title = element_text(hjust=0.5, size = 16)) 
p_ilcri<- plot_proportion_correct(X, C_ilcri) + labs(title="ILCRI") + theme(plot.title = element_text(hjust=0.5, size = 16)) 

ggarrange(p_commonk, p_normative, p_vics, p_vii,
          p_gomrt, p_dlcsl, p_dlctl,
          p_mhm, p_ilcre, p_ilcri, 
          common.legend = TRUE, ncol=2, nrow=5, label.x = 1,heights = c(1, 1, 1, 1, 1.1))
ggsave("680810proportion_correct.pdf")
plot_crop("680810proportion_correct.pdf")

#===========================================================================
#                            Theta comparisons
#===========================================================================
library(MCMCvis)
model <- readRDS("analyses/model-fit_results/680810/commonk.rds")
theta_commonk <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/normative.rds")
theta_normative <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/vics.rds")
theta_vics <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/vii.rds")
theta_vii <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/gomrt.rds")
theta_gomrt <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/dlcsl.rds")
theta_dlcsl <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/dlctl.rds")
theta_dlctl <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/mhm.rds")
theta_mhm <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/ilcre.rds")
theta_ilcre <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

model <- readRDS("analyses/model-fit_results/680810/ilcri.rds")
theta_ilcri <- MCMCsummary(model$posterior_samples$samples,"theta")$mean

theta_compare <- data.frame("Common K" = theta_commonk,
            "Normative" = theta_normative,
            "VICS" = theta_vics,
            "VII" = theta_vii,
            "GoMRT" = theta_gomrt,
            "DLCSL" = theta_dlcsl,
            "DLCTL" = theta_dlctl,
            "MHM" = theta_mhm,
            "ILCRE" = theta_ilcre,
            "ILCRI" = theta_ilcri)
df_long <- theta_compare %>%
  mutate(ID = row_number()) %>%
  pivot_longer(cols = -ID, names_to = "Model", values_to = "Value")
theta_compare_plot<- ggplot(df_long, aes(x = ID, y = Value, color = Model, group = Model)) +
  geom_point(size = 2.5, alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Students",
    y = expression(theta)
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0),
    axis.title.y = element_text(size=15),
    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
  ) + coord_fixed(5)  +
  scale_color_manual(values = c(
    "Common.K" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  ))+ scale_x_continuous(
    breaks = 1:nrow(X),  # Show every student  # Label with student numbers
  )
theta_compare_plot
ggsave("theta_compare.pdf")
plot_crop("theta_compare.pdf")
ggplot(df_long, aes(x = Model, y = Estimate, fill = Model)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Distribution of Parameter Estimates Across Models",
    x = "Model",
    y = "Parameter Estimate"
  )
library(ggpubr)

df_long <- theta_compare %>%
  # Add row identifier
  mutate(Observation = row_number()) %>%
  # Reshape to long format
  pivot_longer(
    cols = -Observation,
    names_to = "Model",
    values_to = "Estimate"
  )
aov_result <- aov(Estimate ~ Model + Error(Observation), data = df_long)
summary(aov_result)

# Error: Observation
#           Df Sum Sq Mean Sq F value Pr(>F)
# Residuals  1   1.75    1.75               

# Error: Within
#            Df Sum Sq Mean Sq F value Pr(>F)
# Model       9    0.8  0.0890    0.13  0.999
# Residuals 659  451.6  0.6854 

lm_model <- lm(Estimate ~ Model, data = df_long)

# Normality of residuals
shapiro.test(residuals(lm_model))

# Homogeneity of variances
bartlett.test(Estimate ~ Model, data = df_long)

# Diagnostic plots
par(mfrow = c(2, 2))
plot(lm_model)
model_summary <- df_long %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(Estimate, na.rm = TRUE),
    SD = sd(Estimate, na.rm = TRUE),
    SE = SD / sqrt(n()),
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE,
    .groups = "drop"
  )

# Create enhanced plot with points, means, and error bars
ability_estimates_compare_plot  <- ggpubr::ggline(model_summary, 
       x = "Model", 
       y = "Mean", 
       group = 1,
       color = "black",
       size = 1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, 
                color = "darkblue") +
  # Add individual data points
  geom_jitter(data = df_long, 
              aes(x = Model, y = Estimate, color = Model),
              width = 0.2, 
              alpha = 0.6, 
              size = 2) +
  labs(
    x = "Model",
    y = expression(theta)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
    axis.text.y = element_text(size=15),
    legend.position = "none",
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=11),
    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
  )+
  scale_color_manual(values = c(
    "Common.K" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) + coord_fixed(0.75)
saveRDS(ability_estimates_compare_plot, "680810ability.rds")

posthoc <- pairwise.t.test(
  df_long$Estimate, 
  df_long$Model,
  paired = TRUE
)
convert_to_stars <- function(p_values) {
  # Replace NA values with empty strings
  p_matrix <- p_values$p.value
  p_matrix[is.na(p_matrix)] <- ""
  
  # Convert p-values to stars
  star_matrix <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
  rownames(star_matrix) <- rownames(p_matrix)
  colnames(star_matrix) <- colnames(p_matrix)
  
  for(i in 1:nrow(p_matrix)) {
    for(j in 1:ncol(p_matrix)) {
      if(p_matrix[i,j] == "") {
        star_matrix[i,j] <- ""
      } else if(p_matrix[i,j] < 0.001) {
        star_matrix[i,j] <- "***"
      } else if(p_matrix[i,j] < 0.01) {
        star_matrix[i,j] <- "**"
      } else if(p_matrix[i,j] < 0.05) {
        star_matrix[i,j] <- "*"
      } else {
        star_matrix[i,j] <- ""
      }
    }
  }
  
  return(star_matrix)
}

# Apply the function to get significance stars
star_table <- convert_to_stars(posthoc)

# Print the star table
print(star_table)
print(posthoc)
ggplot(df_long, aes(x = Model, y = Estimate, fill = Model)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Distribution of Parameter Estimates Across Models",
    x = "Model",
    y = "Parameter Estimate"
  ) + coord_flip()

data.frame(posthoc)
library(tinytable)
ttest_res <- data.frame(posthoc$p.value)
ttest_res$Model <- c("DLCSL", "DLCTL", "GoMRT", "ILCRE", "ILCRI", "MHM", "Normative", "VICS", "VII")
ttest_res <- ttest_res %>%
  select(Model, everything())
ttest_res <-ttest_res %>%
  mutate(across(where(is.numeric), ~formatC(., 
                                           format = "g", 
                                           digits = 2)),
        across(everything(), ~replace_na(., "")),
        across(everything(), 
                ~if(is.character(.)) str_replace_all(., c(" NA" = "", "^NA$" = "")) else .))

ttest_matrix <- as.matrix(ttest_res)
library(tinytable)
tt <- tt(ttest_res,caption = "Bonferroni-Corrected Pairwise T-Tests Comparing Ability Estimates Between Engagement Models for Assessment 747119.")
for(i in 1:nrow(ttest_matrix)) {
  for(j in 2:ncol(ttest_matrix)) {
    # Extract the value as text
    val_text <- ttest_matrix[i,j]
    
    # Skip empty cells
    if(val_text == "") next
    
    # Convert to numeric (handle scientific notation)
    val_num <- as.numeric(val_text)
    
    # Bold if significant
    if(!is.na(val_num) && val_num < 0.05) {
      tt <- style_tt(
        tt, 
        i = i,
        j = j,
        bold = TRUE
      )
    }
  }
}
tt  |> theme_tt("rotate")|> print("latex")

#===========================================================================
#                            Difficulty Comparisons
#===========================================================================

#===========================================================================
#                            Theta comparisons
#===========================================================================
library(MCMCvis)
model <- readRDS("analyses/model-fit_results/680810/commonk.rds")
d_commonk <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/normative.rds")
b_normative <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/vics.rds")
b_vics <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/vii.rds")
b_vii <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/gomrt.rds")
b_gomrt <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/dlcsl.rds")
b_dlcsl <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/dlctl.rds")
b_dlctl <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/mhm.rds")
b_mhm <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/ilcre.rds")
b_ilcre <- MCMCsummary(model$posterior_samples$samples,"b")$mean

model <- readRDS("analyses/model-fit_results/680810/ilcri.rds")
b_ilcri <- MCMCsummary(model$posterior_samples$samples,"b")$mean

b_compare <- data.frame("Common K" = d_commonk,
            "Normative" = b_normative,
            "VICS" = b_vics,
            "VII" = b_vii,
            "GoMRT" = b_gomrt,
            "DLCSL" = b_dlcsl,
            "DLCTL" = b_dlctl,
            "MHM" = b_mhm,
            "ILCRE" = b_ilcre,
            "ILCRI" = b_ilcri)
df_long <- b_compare %>%
  mutate(ID = row_number()) %>%
  pivot_longer(cols = -ID, names_to = "Model", values_to = "Value")



df_long <- b_compare %>%
  # Add row identifier
  mutate(Observation = row_number()) %>%
  # Reshape to long format
  pivot_longer(
    cols = -Observation,
    names_to = "Model",
    values_to = "Estimate"
  )
model_summary <- df_long %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(Estimate, na.rm = TRUE),
    SD = sd(Estimate, na.rm = TRUE),
    SE = SD / sqrt(n()),
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE,
    .groups = "drop"
  )
df_long <-  b_compare %>%
  # Add row identifier
  mutate(Observation = row_number()) %>%
  # Reshape to long format
  pivot_longer(
    cols = -Observation,
    names_to = "Model",
    values_to = "Estimate"
  )
aov_result <- aov(Estimate ~ Model + Error(Observation), data = df_long)
summary(aov_result)

# Error: Observation
#           Df Sum Sq Mean Sq F value Pr(>F)
# Residuals  1  3.641   3.641               

# Error: Within
#           Df Sum Sq Mean Sq F value   Pr(>F)    
# Model      9  5.405  0.6006   5.525 5.04e-06 ***
# Residuals 89  9.675  0.1087                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Create enhanced plot with points, means, and error bars
d_estimates_compare_plot <- ggpubr::ggline(model_summary, 
       x = "Model", 
       y = "Mean", 
       group = 1,
       color = "black",
       size = 1) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.2, 
                color = "darkblue") +
  # Add individual data points
  geom_jitter(data = df_long, 
              aes(x = Model, y = Estimate, color = Model),
              width = 0.2, 
              alpha = 0.6, 
              size = 2) +
  labs(
    x = "Model",
    y = expression(d)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
    axis.text.y = element_text(size=15),
    legend.position = "none",
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=11),
    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
  )+
  scale_color_manual(values = c(
    "Common.K" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  )) + coord_fixed(1.25)
saveRDS(d_estimates_compare_plot, "680810difficulty.rds")

ggarrange(ability_estimates_compare_plot, d_estimates_compare_plot, nrow = 2, 
            labels = c("A","B"), widths = c(1,1), heights = c(1,1), vjust = 7)
ggsave("680810estimates_compare.pdf")
library(knitr)
plot_crop("680810estimates_compare.pdf")
d_compare_plot<- ggplot(df_long, aes(x = ID, y = Value, color = Model, group = Model)) +
  geom_point(size = 2.5, alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Items",
    y = expression(d)
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90),
    axis.title.y = element_text(size=15),
    panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold"),
  ) + coord_fixed(5)  +
  scale_color_manual(values = c(
    "Common.K" = "#6a3d9a",
    "Normative" = "#cab2d6",
    "VICS" = "#5B8FA8FF",
    "VII" = "#0F425CFF",
    "GoMRT" = "#C16622FF",
    "DLCSL" = "#B1746FFF",
    "DLCTL" = "#800000FF",
    "MHM" = "#ADB17DFF",
    "ILCRE" = "#616530FF",
    "ILCRI" = "#3E3E23FF"
  ))+ scale_x_continuous(
    breaks = 1:ncol(X),  # Show every student  # Label with student numbers
    label =  colnames(X)
  ) + coord_fixed(2.75)
ggarrange(theta_compare_plot,d_compare_plot,nrow = 2, common.legend = TRUE, labels=c("A", "B"), hjust = -6)
ggsave("estimates_compare.pdf")
plot_crop("estimates_compare.pdf")


df_long <- b_compare %>%
  # Add row identifier
  mutate(Observation = row_number()) %>%
  # Reshape to long format
  pivot_longer(
    cols = -Observation,
    names_to = "Model",
    values_to = "Estimate"
  )



posthoc <- pairwise.t.test(
  df_long$Estimate, 
  df_long$Model,
  paired = TRUE,
  p.adjust.method = "bonferroni"
)
star_table <- convert_to_stars(posthoc)

print(posthoc)

library(tinytable)
ttest_res <- data.frame(star_table)
ttest_res$Model <- c("DLCSL", "DLCTL", "GoMRT", "ILCRE", "ILCRI", "MHM", "Normative", "VICS", "VII")
ttest_res <- ttest_res %>%
  select(Model, everything())
ttest_res <-ttest_res %>%
  mutate(across(where(is.numeric), ~formatC(., 
                                           format = "g", 
                                           digits = 2)),
        across(everything(), ~replace_na(., "")),
        across(everything(), 
                ~if(is.character(.)) str_replace_all(., c(" NA" = "", "^NA$" = "")) else .))

ttest_matrix <- as.matrix(ttest_res)
tt <- tt(ttest_res,caption = "Bonferroni-Corrected Pairwise T-Tests Comparing Item Difficulty Estimates Between Engagement Models for Assessment 747119.")
for(i in 1:nrow(ttest_matrix)) {
  for(j in 2:ncol(ttest_matrix)) {
    # Extract the value as text
    val_text <- ttest_matrix[i,j]
    
    # Skip empty cells
    if(val_text == "") next
    
    # Convert to numeric (handle scientific notation)
    val_num <- as.numeric(val_text)
    
    # Bold if significant
    if(!is.na(val_num) && val_num < 0.05) {
      tt <- style_tt(
        tt, 
        i = i,
        j = j,
        bold = TRUE
      )
    }
  }
}
tt|> theme_tt("rotate")|> print("latex")

#===========================================================================
#                            PATH ANALYSIS
#===========================================================================
library(MCMCvis)
library(codatools)
X <- data$X
RT <- data$RT
X[is.na(X)] <- 0
RT[is.na(RT)] <- 0
I = nrow(X)
J = ncol(X)
model <- readRDS("analyses/model-fit_results/680810/ilcri.rds")
Prob_X <- MCMCsummary(model$posterior_samples$samples,"correct_resp_prob")$mean
Prob_X <-matrix(Prob_X,I,J)
Prob_C <- MCMCsummary(model$posterior_samples$samples,"class_p")$mean
Prob_C <-matrix(Prob_C,I,J)
theta_ilcri <- MCMCsummary(model$posterior_samples$samples,"theta")$mean
b_ilcri <- MCMCsummary(model$posterior_samples$samples,"b")$mean
xi_ilcri <- MCMCsummary(model$posterior_samples$samples,"xi")$mean
kappa_ilcri <- MCMCsummary(model$posterior_samples$samples,"kappa")$mean
eta_ilcri <- MCMCsummary(model$posterior_samples$samples,"eta")$mean
# C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
# C_ilcri <- matrix(apply(C ,2,smode),I,J)

library(tidyverse)

# Assuming all matrices have the same dimensions
# n_students = number of rows in matrices
# n_items = number of columns in matrices

# First, create indices for all student-item combinations
student_indices <- 1:nrow(Prob_X)
item_indices <- 1:ncol(Prob_X)

# Create a grid of all student-item combinations
student_item_grid <- expand.grid(
  student_id = student_indices,
  item_id = item_indices
)

# Convert matrices to data frames
prob_correct_df <- as.data.frame(Prob_X)
prob_engaged_df <- as.data.frame(Prob_C)  # You mentioned both are similar
rt_df <- as.data.frame(as.matrix(RT))  # Convert to matrix then dataframe to ensure proper conversion

# Add student IDs as row names
rownames(prob_correct_df) <- 1:nrow(prob_correct_df)
rownames(prob_engaged_df) <- 1:nrow(prob_engaged_df)
rownames(rt_df) <- 1:nrow(rt_df)
colnames(rt_df) <- 1:ncol(rt_df)

# Convert each matrix to long format
prob_correct_long <- prob_correct_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_correct") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

prob_engaged_long <- prob_engaged_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_engaged") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

rt_long <- rt_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "response_time") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

# Join the long dataframes
combined_df <- prob_correct_long %>%
  select(student_id, item_id, prob_correct) %>%
  left_join(
    prob_engaged_long %>% select(student_id, item_id, prob_engaged),
    by = c("student_id", "item_id")
  )

combined_df <- combined_df |> left_join(rt_long|>  select(student_id, item_id,response_time) , by = c("student_id", "item_id"))

semdata<- combined_df %>%
  mutate(
    theta = theta_ilcri[student_id],
    b = b_ilcri[item_id],
    kappa = kappa_ilcri[item_id],
    eta = eta_ilcri[student_id],
    xi = xi_ilcri[student_id]
  ) |> select(-c(student_id,item_id))

library(lavaan)
library(lavaanPlot)
library(piecewiseSEM)
library(tidyverse)
library(matrixStats)

model1 <- 'prob_correct ~ theta + prob_engaged + b
           prob_engaged ~ eta + kappa
           response_time ~ prob_engaged + xi
           xi ~ eta
           eta ~ theta
           theta ~ xi'
semdata$response_time <- scale(semdata$response_time)
fit1 <- sem(model1, data = semdata) 
lavaanPlot(name = "MODEL1", fit1,coefs = TRUE)

model <- readRDS("analyses/model-fit_results/680810/ilcre.rds")
Prob_X <- MCMCsummary(model$posterior_samples$samples,"correct_resp_prob")$mean
Prob_X <-matrix(Prob_X,I,J)
Prob_C <- MCMCsummary(model$posterior_samples$samples,"class_p")$mean
Prob_C <-matrix(Prob_C,I,J)
theta_ilcri <- MCMCsummary(model$posterior_samples$samples,"theta")$mean
b_ilcri <- MCMCsummary(model$posterior_samples$samples,"b")$mean
xi_ilcri <- MCMCsummary(model$posterior_samples$samples,"xi")$mean
kappa_ilcri <- MCMCsummary(model$posterior_samples$samples,"kappa")$mean
eta_ilcri <- MCMCsummary(model$posterior_samples$samples,"eta")$mean
# C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
# C_ilcri <- matrix(apply(C ,2,smode),I,J)

library(tidyverse)

# Assuming all matrices have the same dimensions
# n_students = number of rows in matrices
# n_items = number of columns in matrices

# First, create indices for all student-item combinations
student_indices <- 1:nrow(Prob_X)
item_indices <- 1:ncol(Prob_X)

# Create a grid of all student-item combinations
student_item_grid <- expand.grid(
  student_id = student_indices,
  item_id = item_indices
)

# Convert matrices to data frames
prob_correct_df <- as.data.frame(Prob_X)
prob_engaged_df <- as.data.frame(Prob_C)  # You mentioned both are similar
rt_df <- as.data.frame(as.matrix(RT))  # Convert to matrix then dataframe to ensure proper conversion

# Add student IDs as row names
rownames(prob_correct_df) <- 1:nrow(prob_correct_df)
rownames(prob_engaged_df) <- 1:nrow(prob_engaged_df)
rownames(rt_df) <- 1:nrow(rt_df)
colnames(rt_df) <- 1:ncol(rt_df)

# Convert each matrix to long format
prob_correct_long <- prob_correct_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_correct") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

prob_engaged_long <- prob_engaged_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_engaged") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

rt_long <- rt_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "response_time") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

# Join the long dataframes
combined_df <- prob_correct_long %>%
  select(student_id, item_id, prob_correct) %>%
  left_join(
    prob_engaged_long %>% select(student_id, item_id, prob_engaged),
    by = c("student_id", "item_id")
  )

combined_df <- combined_df |> left_join(rt_long|>  select(student_id, item_id,response_time) , by = c("student_id", "item_id"))

semdata<- combined_df %>%
  mutate(
    theta = theta_ilcri[student_id],
    b = b_ilcri[item_id],
    kappa = kappa_ilcri[item_id],
    eta = eta_ilcri[student_id],
    xi = xi_ilcri[student_id]
  ) |> select(-c(student_id,item_id))
semdata$response_time <- scale(semdata$response_time)
library(lavaan)
library(lavaanPlot)
library(piecewiseSEM)
library(tidyverse)
library(matrixStats)

model1 <- 'prob_correct ~ theta + prob_engaged + b
           prob_engaged ~ eta + kappa
           response_time ~ prob_engaged + xi
           xi ~ eta
           eta ~ theta
           theta ~ xi'
fit1 <- sem(model1, data = semdata) 
lavaanPlot(name = "MODEL1", fit1,coefs = TRUE)


model <- readRDS("analyses/model-fit_results/680810/mhm.rds")
Prob_X <- MCMCsummary(model$posterior_samples$samples,"correct_resp_prob")$mean
Prob_X <-matrix(Prob_X,I,J)
Prob_C <- MCMCsummary(model$posterior_samples$samples,"class_p")$mean
Prob_C <-matrix(Prob_C,I,J)
theta_ilcri <- MCMCsummary(model$posterior_samples$samples,"theta")$mean
b_ilcri <- MCMCsummary(model$posterior_samples$samples,"b")$mean
#xi_ilcri <- MCMCsummary(model$posterior_samples$samples,"xi")$mean
kappa_ilcri <- MCMCsummary(model$posterior_samples$samples,"kappa")$mean
eta_ilcri <- MCMCsummary(model$posterior_samples$samples,"eta")$mean
# C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
# C_ilcri <- matrix(apply(C ,2,smode),I,J)

library(tidyverse)

# Assuming all matrices have the same dimensions
# n_students = number of rows in matrices
# n_items = number of columns in matrices

# First, create indices for all student-item combinations
student_indices <- 1:nrow(Prob_X)
item_indices <- 1:ncol(Prob_X)

# Create a grid of all student-item combinations
student_item_grid <- expand.grid(
  student_id = student_indices,
  item_id = item_indices
)

# Convert matrices to data frames
prob_correct_df <- as.data.frame(Prob_X)
prob_engaged_df <- as.data.frame(Prob_C)  # You mentioned both are similar
rt_df <- as.data.frame(as.matrix(RT))  # Convert to matrix then dataframe to ensure proper conversion

# Add student IDs as row names
rownames(prob_correct_df) <- 1:nrow(prob_correct_df)
rownames(prob_engaged_df) <- 1:nrow(prob_engaged_df)
rownames(rt_df) <- 1:nrow(rt_df)
colnames(rt_df) <- 1:ncol(rt_df)

# Convert each matrix to long format
prob_correct_long <- prob_correct_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_correct") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

prob_engaged_long <- prob_engaged_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_engaged") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

rt_long <- rt_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "response_time") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

# Join the long dataframes
combined_df <- prob_correct_long %>%
  select(student_id, item_id, prob_correct) %>%
  left_join(
    prob_engaged_long %>% select(student_id, item_id, prob_engaged),
    by = c("student_id", "item_id")
  )

combined_df <- combined_df |> left_join(rt_long|>  select(student_id, item_id,response_time) , by = c("student_id", "item_id"))

semdata<- combined_df %>%
  mutate(
    theta = theta_ilcri[student_id],
    b = b_ilcri[item_id],
    kappa = kappa_ilcri[item_id],
    eta = eta_ilcri[student_id],
  ) |> select(-c(student_id,item_id))
semdata$response_time <- scale(semdata$response_time)
library(lavaan)
library(lavaanPlot)
library(piecewiseSEM)
library(tidyverse)
library(matrixStats)

model1 <- 'prob_correct ~ theta + prob_engaged + b
           prob_engaged ~ eta + kappa
           response_time ~ prob_engaged
      
           eta ~ theta'
fit1 <- sem(model1, data = semdata) 
lavaanPlot(name = "MODEL1", fit1,coefs = TRUE)


model <- readRDS("analyses/model-fit_results/680810/dlctl.rds")
Prob_X <- MCMCsummary(model$posterior_samples$samples,"correct_resp_prob")$mean
Prob_X <-matrix(Prob_X,I,J)
Prob_C <- MCMCsummary(model$posterior_samples$samples,"class_p")$mean
Prob_C <-matrix(Prob_C,I,J)
theta_ilcri <- MCMCsummary(model$posterior_samples$samples,"theta")$mean
b_ilcri <- MCMCsummary(model$posterior_samples$samples,"b")$mean
#xi_ilcri <- MCMCsummary(model$posterior_samples$samples,"xi")$mean
tau_ilcri <- MCMCsummary(model$posterior_samples$samples,"tau")$mean
zeta_ilcri <- MCMCsummary(model$posterior_samples$samples,"zeta")$mean
# C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
# C_ilcri <- matrix(apply(C ,2,smode),I,J)

library(tidyverse)

# Assuming all matrices have the same dimensions
# n_students = number of rows in matrices
# n_items = number of columns in matrices

# First, create indices for all student-item combinations
student_indices <- 1:nrow(Prob_X)
item_indices <- 1:ncol(Prob_X)

# Create a grid of all student-item combinations
student_item_grid <- expand.grid(
  student_id = student_indices,
  item_id = item_indices
)

# Convert matrices to data frames
prob_correct_df <- as.data.frame(Prob_X)
prob_engaged_df <- as.data.frame(Prob_C)  # You mentioned both are similar
rt_df <- as.data.frame(as.matrix(RT))  # Convert to matrix then dataframe to ensure proper conversion

# Add student IDs as row names
rownames(prob_correct_df) <- 1:nrow(prob_correct_df)
rownames(prob_engaged_df) <- 1:nrow(prob_engaged_df)
rownames(rt_df) <- 1:nrow(rt_df)
colnames(rt_df) <- 1:ncol(rt_df)

# Convert each matrix to long format
prob_correct_long <- prob_correct_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_correct") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

prob_engaged_long <- prob_engaged_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_engaged") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

rt_long <- rt_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "response_time") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

# Join the long dataframes
combined_df <- prob_correct_long %>%
  select(student_id, item_id, prob_correct) %>%
  left_join(
    prob_engaged_long %>% select(student_id, item_id, prob_engaged),
    by = c("student_id", "item_id")
  )

combined_df <- combined_df |> left_join(rt_long|>  select(student_id, item_id,response_time) , by = c("student_id", "item_id"))

semdata<- combined_df %>%
  mutate(
    theta = theta_ilcri[student_id],
    b = b_ilcri[item_id],
    tau = tau_ilcri[item_id],
    zeta = eta_ilcri[student_id],
  ) |> select(-c(student_id,item_id))
semdata$response_time <- scale(semdata$response_time)
library(lavaan)
library(lavaanPlot)
library(piecewiseSEM)
library(tidyverse)
library(matrixStats)

model1 <- 'prob_correct ~ theta + prob_engaged + b
           prob_engaged ~ zeta + tau
           response_time ~ prob_engaged
      
           zeta ~ theta'
fit1 <- sem(model1, data = semdata) 
lavaanPlot(name = "MODEL1", fit1,coefs = TRUE)



model <- readRDS("analyses/model-fit_results/680810/dlcsl.rds")
Prob_X <- MCMCsummary(model$posterior_samples$samples,"correct_resp_prob")$mean
Prob_X <-matrix(Prob_X,I,J)
Prob_C <- MCMCsummary(model$posterior_samples$samples,"class_p")$mean
Prob_C <-matrix(Prob_C,I,J)
theta_ilcri <- MCMCsummary(model$posterior_samples$samples,"theta")$mean
b_ilcri <- MCMCsummary(model$posterior_samples$samples,"b")$mean
#xi_ilcri <- MCMCsummary(model$posterior_samples$samples,"xi")$mean
tau_ilcri <- MCMCsummary(model$posterior_samples$samples,"tau")$mean
#zeta_ilcri <- MCMCsummary(model$posterior_samples$samples,"zeta")$mean
# C <- coda_grep(model$posterior_samples$samples,"C", return.matrix = TRUE)
# C_ilcri <- matrix(apply(C ,2,smode),I,J)

library(tidyverse)

# Assuming all matrices have the same dimensions
# n_students = number of rows in matrices
# n_items = number of columns in matrices

# First, create indices for all student-item combinations
student_indices <- 1:nrow(Prob_X)
item_indices <- 1:ncol(Prob_X)

# Create a grid of all student-item combinations
student_item_grid <- expand.grid(
  student_id = student_indices,
  item_id = item_indices
)

# Convert matrices to data frames
prob_correct_df <- as.data.frame(Prob_X)
prob_engaged_df <- as.data.frame(Prob_C)  # You mentioned both are similar
rt_df <- as.data.frame(as.matrix(RT))  # Convert to matrix then dataframe to ensure proper conversion

# Add student IDs as row names
rownames(prob_correct_df) <- 1:nrow(prob_correct_df)
rownames(prob_engaged_df) <- 1:nrow(prob_engaged_df)
rownames(rt_df) <- 1:nrow(rt_df)
colnames(rt_df) <- 1:ncol(rt_df)

# Convert each matrix to long format
prob_correct_long <- prob_correct_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_correct") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

prob_engaged_long <- prob_engaged_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "prob_engaged") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

rt_long <- rt_df %>%
  rownames_to_column("student_id") %>%
  pivot_longer(-student_id, 
               names_to = "item_name", 
               values_to = "response_time") %>%
  mutate(student_id = as.integer(student_id),
         item_id = as.integer(str_replace(item_name, "V", "")))

# Join the long dataframes
combined_df <- prob_correct_long %>%
  select(student_id, item_id, prob_correct) %>%
  left_join(
    prob_engaged_long %>% select(student_id, item_id, prob_engaged),
    by = c("student_id", "item_id")
  )

combined_df <- combined_df |> left_join(rt_long|>  select(student_id, item_id,response_time) , by = c("student_id", "item_id"))

semdata<- combined_df %>%
  mutate(
    theta = theta_ilcri[student_id],
    b = b_ilcri[item_id],
    tau = tau_ilcri[item_id]
  ) |> select(-c(student_id,item_id))
semdata$response_time <- scale(semdata$response_time)
library(lavaan)
library(lavaanPlot)
library(piecewiseSEM)
library(tidyverse)
library(matrixStats)

model1 <- 'prob_correct ~ theta + prob_engaged + b
           prob_engaged ~ tau
           response_time ~ prob_engaged'
fit1 <- sem(model1, data = semdata) 
lavaanPlot(name = "MODEL1", fit1,coefs = TRUE)