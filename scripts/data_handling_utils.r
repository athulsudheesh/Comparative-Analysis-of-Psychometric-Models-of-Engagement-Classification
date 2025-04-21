library(tidyverse)
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