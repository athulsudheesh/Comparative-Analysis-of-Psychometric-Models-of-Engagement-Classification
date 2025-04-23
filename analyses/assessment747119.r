source("scripts/model-fits.r") 
source("scripts/data_handling_utils.r")
data <- load_and_clean_data(747119)
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
saveRDS(model, "analyses/model-fit_results/commonk.rds")
#model_fit_results <- model_fit_results[-nrow(model_fit_results), ]

model <- Normative.fit(data)
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "Normative",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/normative.rds")


model <- VICS.fit(data)
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "VICS",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/vics.rds")

model <- VII.fit(data)
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "VII",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/vii.rds")


model <- GoMRT.fit(data, 
                   params.to.monitor = c("correct_resp_prob",
                                          "class_p","C"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "GoMRT",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/gomrt.rds")


model <- DLCSL.fit(data, 
                   params.to.monitor = c("correct_resp_prob", 
                                          "class_p","C"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "DLCSL",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/dlcsl.rds")


model <- DLCTL.fit(data, 
                   params.to.monitor = c("correct_resp_prob",
                                          "class_p","C","theta","b"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "DLCTL",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/dlctl.rds")

model <- nimbleModel(code = MHM, constants = const, 
                     data = data, inits = inits(),
                     check = TRUE, calculate = FALSE)
model <- MHM.fit(data, 
                   params.to.monitor = c("correct_resp_prob", "theta","b",
                                          "class_p","C"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "MHM",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/mhm.rds")


model <- ILCRE.fit(data, 
                   params.to.monitor = c("correct_resp_prob","theta","b",
                                          "class_p","C", "eta","xi","kappa"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "ILCRE",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/ilcre.rds")

model <- ILCRI.fit(data, 
                   params.to.monitor = c("correct_resp_prob", "theta","b",
                                          "class_p","C", "eta","xi","kappa"))
model_fit_results <- rbind(model_fit_results, data.frame("Model" = "ILCRI",
                        "lppd"=model$posterior_samples$WAIC$lppd,
                        "pWAIC"=model$posterior_samples$WAIC$pWAIC,
                        "WAIC"=model$posterior_samples$WAIC$WAIC,
                        "MSE"= model$cvvals$CVvalue,"SE"= model$cvvals$CVstandardError))
saveRDS(model, "analyses/model-fit_results/ilcri.rds")




saveRDS(model_fit_results_reordered,"analyses/model_fit_results.rds")

#===========================================================================
#                            Model Fit results table
#===========================================================================
model_fit_results <- readRDS("analyses/model_fit_results.rds")
library(tinytable)
model_fit_results$SE <- formatC(model_fit_results$SE,format = "e", digits = 2)
tt(model_fit_results,
caption = "Model fit results from MCMC analysis for all models applied to Assessment 747119 data") |>
            group_tt(j = list(
                "3-fold CV Error" = 5:6
            )) |>
            format_tt(j = 2:6, digits = 2,num_fmt = "decimal") |>
            style_tt(i = 7, bold=TRUE) |>
            print("latex")

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
model <- readRDS("analyses/model-fit_results/commonk.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "Common-k", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))


model <- readRDS("analyses/model-fit_results/normative.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "Normative", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/vics.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "VICS", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/vii.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "VII", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= "", "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/gomrt.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "GoMRT", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/dlcsl.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "DLCSL", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/dlctl.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "DLCTL", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/mhm.rds")
Thetacoverage <-  compute_coverage(model,"theta")
Bcoverage <- compute_coverage(model, "b")
Ccoverage <- compute_coverage(model,"class_p")
coverage_results <- rbind(coverage_results, 
                        data.frame("Model" = "MHM", "Theta" = Thetacoverage,
                        "d" = Bcoverage,"C"= Ccoverage, "eta" = "", "kappa" = "", "xi"= ""))

model <- readRDS("analyses/model-fit_results/ilcre.rds")
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

model <- readRDS("analyses/model-fit_results/ilcri.rds")
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

saveRDS(coverage_results,"analyses/model-fit_results/coverage_results.rds")
library(tinytable)
coverage_results<- coverage_results |> 
  mutate_at(vars("C", "eta","kappa","xi"), as.numeric)
tt(coverage_results,
    caption = "MCMC parameter coverage based on Rhat statistics for Assessment 747119") |>
     format_tt(j = 3:7, digits = 2,num_fmt = "decimal", replace = TRUE) |>
     style_tt(i =5:8, j = 4, color="red") |>
     style_tt(i =9, j = 7, color="red") |>
     print("latex")

#===========================================================================
#                            debugging
#===========================================================================


coverage_results <- head(coverage_results, -1)
