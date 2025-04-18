library(tinytable)
library(tidyverse)
library(psych)
library(ggplot2)
library(ggridges)

library(dplyr)
#===========================================================================
#                            FUNCTION DEFINITIONS
#===========================================================================

load_data <- function(assessment_id, directory = "data"){
    all_files <- list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)
    matching_files <- all_files[grepl(as.character(assessment_id), basename(all_files))]
    print(matching_files)
   RT <- read_csv(matching_files[1]) |> select(-1)
   X <- read_csv(matching_files[2]) |> select(-1)
      return(list(X=X,RT=RT))
}


load_and_clean_data <- function(assessment_id, directory = "data"){
    suppressMessages({ 
    data <- load_data(assessment_id)
    })
    #-------- removing rows where the no of missing is more than 60% ---------
    data$X <- data$X |> 
        filter(rowSums(is.na(across(everything()))) <= 0.60 * ncol(data$X))
    data$RT <- data$RT |> 
        filter(rowSums(is.na(across(everything()))) <= 0.60 * ncol(data$X))
    #-------- removing COLS where the no of missing is more than 60% ---------
    data$X <- data$X |> select(where(~ sum(is.na(.)) <= 0.75 * nrow(data$X)))
    data$RT <- data$RT |> select(where(~ sum(is.na(.)) <= 0.75 * nrow(data$RT)))
    return(list(X = data$X, RT = data$RT))
}
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
    p <- ggplot(longRT, aes(x = response_time, y = item, fill = item)) +
    geom_density_ridges(alpha = 0.4, scale = 5,rel_min_height = 0.01,bandwidth = 1.5,fill = "#800000FF", color=NA) +
        theme_bw() +
        theme(
                plot.title = element_text(hjust=0.5),
                axis.text.x = element_text(angle = 0, size = 10),
                axis.text.y = element_text(angle = 0, size = 7.5),
                axis.title.x = element_text(size = 10),
                panel.border = element_rect(size = 1.5),
                axis.title = element_text(face = "bold"),
                legend.position = "none"
            ) +
        labs(
            title = as.character(assessment_id),
            x = "Response Time (s)",
            y = "Items"
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

 |> mutate(student = factor(student, levels = sort(unique(student))))
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

b
a
ggsave("test.pdf")

j <- plot_item_response_times(747119)
k <- plot_item_response_times(680810)
l <- plot_item_response_times(735540)
m <- plot_item_response_times(763628)
plot_grid(j,k,l,m,align="hv")
ggsave("RTItems.pdf")


#==========================================================================================================================
#                                                     WORK IN PROGRESS SCRIPTS
#==========================================================================================================================
ggplot()
dat <- load_and_clean_data(763628)
longRT <- dat$RT |> rownames_to_column("student") |> pivot_longer(
        cols = -student,
        names_to = "item",
        values_to = "response_time"
        )

scores <- dat$X |> rowSums() |> as.data.frame() |> rownames_to_column("student")

left_join(longRT)



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
ggplot(longRT, aes(x = response_time, y = student, fill = score)) +
    geom_density_ridges(alpha = 0.4, scale = 10,rel_min_height = 0.01,bandwidth = 1.5, color=NA) +
        theme_bw() +
        theme(
                plot.title = element_text(hjust=0.5),
                axis.text.x = element_text(angle = 0, size = 10),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.x = element_text(size = 10),
                panel.border = element_rect(size = 1.5),
                axis.title = element_text(face = "bold")
            ) +
        labs(
            title = as.character(1),
            x = "Response Time (s)",
            y = "Students"
        )+
  scale_fill_manual(
    values = c("H" = "#8A9045FF", "M" = "#155F83FF", "L" = "#800000FF"),
    name = "Score",
    labels = c("High", "Medium", "Low")
  )
ggsave('test.pdf')

