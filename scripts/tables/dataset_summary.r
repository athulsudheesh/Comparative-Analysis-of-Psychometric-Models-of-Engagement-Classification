library(tinytable)

dat <- data.frame(
    school_id = c(5759,5759,5450,5450),
    class_id = c(23078, 23078, 23017, 22755),
    teacher_id = c(156390,156390,59965,68444),
    assessment_id = c(747119,591509,735540,633168),
    tutor_mode = c("No Feedback", "With Feedback","With Feedback", "With Feedback")
)
tt(dat,
    caption = "Overview of Datasets considered in this thesis.") |> setNames(c("School ID", "Class ID", "Teacher ID", "Assessment ID", "Tutor Mode")) |>
    style_tt(align="c")|>
    print("latex")
