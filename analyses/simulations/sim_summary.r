C_90 <- readRDS("analyses/simulations/C_90/sim_data.rds")
C_60 <- readRDS("analyses/simulations/C_60/sim_data.rds")
C_30 <- readRDS("analyses/simulations/C_30/sim_data.rds")
C_10 <- readRDS("analyses/simulations/C_10/sim_data.rds")


get_summarystats_sim <- function(sim_data){

    scores <- rowSums(sim_data$Y,na.rm=TRUE)
    meanRT <- rowMeans(exp(sim_data$T_mat),na.rm=TRUE)
    alpha <- suppressMessages({ psych::alpha(sim_data$Y)$total[[1]]})
    return(list(
                scores_m = mean(scores),
                scores_sd = sd(scores),
                RT_m= mean(meanRT),
                RT_sd = sd(meanRT),
            
                Calpha = alpha))
}
summ <- list()
sim_summ <- get_summarystats_sim(C_10)
summ <- rbind(summ,sim_summ)

sim_summ <- get_summarystats_sim(C_30)
summ <- rbind(summ,sim_summ)

sim_summ <- get_summarystats_sim(C_60)
summ <- rbind(summ,sim_summ)

sim_summ <- get_summarystats_sim(C_90)
summ <- rbind(summ,sim_summ)

summ <- cbind(p_engagement = c(10, 30, 60, 90), summ)
summ_df <- as.data.frame(summ)
summ_df <- sapply(summ_df, as.numeric)

library(tinytable)
library(scales)
tt(as.data.frame(summ_df),
digits = 2,
    caption = "Summary Statistics of Simulated Datasets") |>
    setNames(c("Simulation Condition \\\\ Prob. of Engagement", "M", "SD", "M", "SD", "Cronbach's alpha")) |>
    group_tt(j = list(
                "{\\\\ Score}" = 2:3,
                "{Avg. \\\\ Response \\\\ Time (s)}" = 4:5
            ))|>
    style_tt(align="c")

plot_score_distributions <- function(dat, cond){
    #dat <- load_and_clean_data(assessment_id)
    longRT <- data.frame(dat$T_mat) |> rownames_to_column("student") |> pivot_longer(
            cols = -student,
            names_to = "item",
            values_to = "response_time"
            )
    scores <- dat$Y |> rowSums(na.rm = TRUE) |> as.data.frame() |> rownames_to_column("student") |> rename(ts = "rowSums(dat$Y, na.rm = TRUE)")
    ggdensity(scores,x="ts", color = "#155F83FF")+theme_bw() +
                theme(legend.position="none",
                        plot.title = element_text(hjust=0.5),
                        axis.text.x = element_text(angle = 0, size = 10),
                        axis.text.y = element_text(angle = 0, size = 10),
                        axis.title.x = element_text(size = 10),
                        panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
                    ) +
                    labs(title = cond,x = "Sum Score", y = "Density") + xlim(0,15)
}


plot_RT_distributions <- function(dat, cond){
    dat$T_mat <- exp(dat$T_mat)
    longRT <- data.frame(dat$T_mat) |> rownames_to_column("student") |> pivot_longer(
            cols = -student,
            names_to = "item",
            values_to = "response_time"
            )
    RTs <- dat$T_mat |> rowMeans(na.rm = TRUE) |> as.data.frame() |> rownames_to_column("student") |> rename(ts = "rowMeans(dat$T_mat, na.rm = TRUE)")
    ggdensity(RTs,x="ts", color = "#800000FF")+ theme_bw() +
                theme(legend.position="none",
                        plot.title = element_text(hjust=0.5),
                        axis.text.x = element_text(angle = 0, size = 10),
                        axis.text.y = element_text(angle = 0, size = 10),
                        axis.title.x = element_text(size = 10),
                        panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
                    ) +
                    labs(title = cond,x = "Avg. Response Time (s)", y = "Density") + xlim(0,50)
}

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)
a <- plot_score_distributions(C_10, "Engagement Probability: 10%")
b <- plot_score_distributions(C_30, "Engagement Probability: 30%")
c <- plot_score_distributions(C_60, "Engagement Probability: 60%")
d <- plot_score_distributions(C_90, "Engagement Probability: 90%")
x<- ggarrange(a,b,c,d, nrow=4)

l <- plot_RT_distributions(C_10, "Engagement Probability: 10%")
m <- plot_RT_distributions(C_30, "Engagement Probability: 30%")
n <- plot_RT_distributions(C_60, "Engagement Probability: 60%")
o <- plot_RT_distributions(C_90, "Engagement Probability: 90%")
y <- ggarrange(l,m,n,o,nrow = 4)
ggarrange(x,y,ncol = 2,labels = c("A","B"))
ggsave("Newsim_distrbutions.pdf")
library(knitr)
plot_crop("Newsim_distrbutions.pdf")
