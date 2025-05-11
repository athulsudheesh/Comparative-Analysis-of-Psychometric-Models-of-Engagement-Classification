library(tinytable)
library(tidyverse)
library(psych)
library(ggplot2)
library(ggridges)

library(dplyr)
#===========================================================================
#                            FUNCTION DEFINITIONS
#===========================================================================

source("scripts/data_handling_utils.r")
get_summarystats <- function(assessment_id){
    data <- load_and_clean_data(assessment_id)
    
    NStudents <- dim(data$X)[1]
    NItems <- dim(data$X)[2]
    scores <- rowSums(data$X,na.rm=TRUE)
    meanRT <- rowMeans(data$RT,na.rm=TRUE)
    percentNA <-  (sum(is.na(data$X)) / prod(dim(data$X))) * 100
    alpha <- suppressMessages({ psych::alpha(data$X)$total[[1]]})
    return(list(ID = assessment_id,
                NStudents = NStudents,
                NItems = NItems,
                scores_m = mean(scores),
                scores_sd = sd(scores),
                RT_m= mean(meanRT),
                RT_sd = sd(meanRT),
                percentNA = percentNA,
                Calpha = alpha))
}
#===========================================================================
#                            Ploting: Item Response Time
#===========================================================================
plot_item_response_times <- function(assessment_id){
     dat <- load_and_clean_data(assessment_id)
    longRT <- dat$RT |> rownames_to_column("student") |> pivot_longer(
        cols = -student,
        names_to = "item",
        values_to = "response_time"
        )

    D <- dat$X |> colMeans(na.rm = TRUE) |> as.data.frame() |> rownames_to_column() |> 
            rename(item = "rowname" , d="colMeans(dat$X, na.rm = TRUE)")|>
            mutate(difficulty = ifelse(d > 0.5, "H", "L")) |>
            select(c(item, difficulty))

    longRT <- left_join(longRT,D)
    p <- ggplot(longRT, aes(x = response_time, y = item, fill = difficulty)) +
        geom_density_ridges(alpha = 0.4, scale = 10,rel_min_height = 0.01,bandwidth = 1.5, color=NA) +
            theme_bw() +
            theme(legend.position="none",
                    plot.title = element_text(hjust=0.5, size =),
                    axis.text.x = element_text(angle = 0, size = 10),
                    axis.text.y = element_text(angle = 0, size = 7.5),
                    axis.title.x = element_text(size = 10),
                    panel.border = element_rect(size = 1.5),
                    axis.title = element_text(face = "bold")
                ) +
            labs(
                title = as.character(assessment_id),
                x = "Response Time (s)",
                y = "Items"
            )+
    scale_fill_manual(
        values = c("H" = "#1b9e77", "L" = "#7570b3"),
        name = "Difficulty",
        labels = c("High",  "Low")
    )
    
}
#===========================================================================
#                            Ploting: Student Response Time
#===========================================================================
plot_student_response_times <- function(assessment_id){
    dat <- load_and_clean_data(assessment_id)
    longRT <- dat$RT |> rownames_to_column("student") |> pivot_longer(
            cols = -student,
            names_to = "item",
            values_to = "response_time"
            )

    scores <- dat$X |> rowSums() |> as.data.frame() |> rownames_to_column("student")

    score_d <- dat$X |> 
    as.data.frame() |>
    rownames_to_column("student") |>
    mutate(
        sum = rowSums(across(-student), na.rm=TRUE)
    ) |>
    mutate(
        # Create three equally sized groups
        score = case_when(
        sum > quantile(sum, 2/3) ~ "H",
        sum > quantile(sum, 1/3) ~ "M",
        TRUE ~ "L"
        )
    ) |> select(student, score)

    longRT<- left_join(longRT, score_d)
    p<- ggplot(longRT, aes(x = response_time, y = student, fill = score)) +
        geom_density_ridges(alpha = 0.4, scale = 10,rel_min_height = 0.01,bandwidth = 1.5, color=NA) +
            theme_bw() +
            theme(legend.position="none",
                    plot.title = element_text(hjust=0.5),
                    axis.text.x = element_text(angle = 0, size = 10),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.x = element_text(size = 10),
                    panel.border = element_rect(size = 1.5),
                    axis.title = element_text(face = "bold")
                ) +
            labs(
                title = as.character(assessment_id),
                x = "Response Time (s)",
                y = "Students"
            )+
    scale_fill_manual(
        values = c("H" = "#1b9e77", "M" = "#d95f02", "L" = "#7570b3"),
        name = "Score",
        labels = c("High", "Medium", "Low")
    )
}

plot_score_distributions <- function(assessment_id){
    dat <- load_and_clean_data(assessment_id)
    longRT <- dat$RT |> rownames_to_column("student") |> pivot_longer(
            cols = -student,
            names_to = "item",
            values_to = "response_time"
            )
    scores <- dat$X |> rowSums(na.rm = TRUE) |> as.data.frame() |> rownames_to_column("student") |> rename(ts = "rowSums(dat$X, na.rm = TRUE)")
    ggdensity(scores,x="ts", color = "#155F83FF")+theme_bw() +
                theme(legend.position="none",
                        plot.title = element_text(hjust=0.5),
                        axis.text.x = element_text(angle = 0, size = 10),
                        axis.text.y = element_text(angle = 0, size = 10),
                        axis.title.x = element_text(size = 10),
                        panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
                    ) +
                    labs(title = as.character(assessment_id),x = "Sum Score", y = "Density")
}

plot_RT_distributions <- function(assessment_id){
    dat <- load_and_clean_data(assessment_id)
    longRT <- dat$RT |> rownames_to_column("student") |> pivot_longer(
            cols = -student,
            names_to = "item",
            values_to = "response_time"
            )
    RTs <- dat$RT |> rowMeans(na.rm = TRUE) |> as.data.frame() |> rownames_to_column("student") |> rename(ts = "rowMeans(dat$RT, na.rm = TRUE)")
    ggdensity(RTs,x="ts", color = "#800000FF")+ theme_bw() +
                theme(legend.position="none",
                        plot.title = element_text(hjust=0.5),
                        axis.text.x = element_text(angle = 0, size = 10),
                        axis.text.y = element_text(angle = 0, size = 10),
                        axis.title.x = element_text(size = 10),
                        panel.border = element_rect(size = 1.5),
                        axis.title = element_text(face = "bold")
                    ) +
                    labs(title = as.character(assessment_id),x = "Avg. Response Time (s)", y = "Density")
}
#===========================================================================
#                            SCRIPTS
#===========================================================================
#------------------- Summary of datasets -------------------


dat <- data.frame(
    school_id = c(5759,7007,5450,11998),
    class_id = c(23078,23428,23017,24902),
    teacher_id = c(156390,98463,59965,191175),
    assessment_id = c(747119,680810,735540,763628),
    tutor_mode = c("No Feedback", "With Feedback","With Feedback", "With Feedback")
)
tt(dat,
    caption = "Overview of Datasets considered in this thesis.") |> setNames(c("School ID", "Class ID", "Teacher ID", "Assessment ID", "Tutor Mode")) |>
    style_tt(align="c")|>
    print("latex")
#-------------------------------------------------------------------------------------

#--------------------------------- Summary Statistics --------------------------------
dat <- get_summarystats(821218)
assessments <- c("584955", "585165","591509","616900","633168","678696","680810",
                "735540","736883","747119","763628","821218")
assessments <- c(747119,680810,735540,763628)
summary_table <- list()
for(assessment in assessments){
    summary_table <- rbind(summary_table, get_summarystats(assessment))
    }

summary_stats <- data.frame(summary_table) 
summary_stats[] <- sapply(summary_stats, as.numeric) 
library(tinytable)
tt(summary_stats, digits = 2,
    caption = "Summary Statistics of Assessments") |>
    setNames(c("{Assessment \\\\ ID}", "{No. of \\\\ Students}", "{No. of \\\\ Items}", "M", "SD",
            "M", "SD", "{\\% Missing}", "{Cronbach's \\\\ alpha}")) |>
    group_tt(j = list(
                "{\\\\ Score}" = 4:5,
                "{Avg. \\\\ Response \\\\ Time (s)}" = 6:7
            )) |>
    style_tt(align="c") |>
    print("latex")
#-------------------------------------------------------------------------------------

library(ggplot2)
library(ggridges)
library(viridis)
library(hrbrthemes)


library(ggpubr)
a<- plot_student_response_times(747119)
b<- plot_student_response_times(680810)
c<- plot_student_response_times(735540)
d<- plot_student_response_times(763628)
ggarrange(a,b,c,d,legend="bottom", common.legend=TRUE)
ggsave("RTStudents.pdf")
plot_crop("RTStudents.pdf")

ggsave("test.pdf")

j <- plot_item_response_times(747119)
k <- plot_item_response_times(680810)
l <- plot_item_response_times(735540)
m <- plot_item_response_times(763628)
ggarrange(j,k,l,m,legend="bottom", common.legend=TRUE) 
ggsave("RTItems.pdf")
plot_crop("RTItems.pdf")


library(knitr)
a<- plot_score_distributions(747119)
b<- plot_score_distributions(680810)
c<- plot_score_distributions(735540)
d<- plot_score_distributions(763628)
z<- ggarrange(a,b,c,d,nrow=4)
plot_crop("ScoreDistributions.pdf")
ggsave('ScoreDistributions.pdf')
system2(command = "pdfcrop", 
        args    = c("ScoreDistributions.pdf", 
                    "ScoreDistributions.pdf") 
        )

i<- plot_RT_distributions(747119)
j<- plot_RT_distributions(680810)
k<- plot_RT_distributions(735540)
l<- plot_RT_distributions(763628)
x<- ggarrange(i,j,k,l,nrow=4) 
ggarrange(z,x, labels = c("A", "B"))
ggsave('ScoresAndRTDistributions.pdf')
plot_crop('ScoresAndRTDistributions.pdf')

#==========================================================================================================================
#                                                     Simulation Summary
#==========================================================================================================================