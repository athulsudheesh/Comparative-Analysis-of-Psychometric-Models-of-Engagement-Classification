files.sources = list.files("models", full.names=TRUE)
sapply(files.sources, source)
source("scripts/threshold_functions.r")

model <- nimbleModel(code=ILCRE, constants = list(I=1,J=1))
model$plotGraph()
