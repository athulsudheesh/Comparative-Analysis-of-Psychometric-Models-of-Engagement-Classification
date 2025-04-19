# The below code is from Rios paper...need lots of tweaking to make it
# work for my thesis

# R codes for calculating the CUMP threshold value
Chance = 1/4  # the random guessing chance correct for an item with 4 options
ULimit = 100  # the upper limit for threshold; user choice
LLimit = 5    # the lower limit for threshold; user choice
# x[,j]       # the response time matrix; the response time of item j
# y[,j]       # the response (0/1) matrix; the response (0/1) of item j

dat <- load_and_clean_data(763628)
dat$RT[is.na(dat$RT)] <- 0
dat$X[is.na(dat$X)] <- 0
x <- as.matrix(dat$RT)
y <- as.matrix(dat$X)
CUMP <- function(j) {
  item <- table(x[, j], y[, j])
  item <- as.matrix(item)
  S = as.numeric(row.names(item))
  Pplus = item[, 2]/apply(item, 1, sum)  # P+ at each RT
  cumP = Pplus
  K = length(item[, 2])
  
  # calculate cumulative P+ for each RT
  for (k in 2:K) {
    cumP[k] = sum(item[1:k, 2])/sum(item[1:k,])
  }
  
  ind = (cumP < Chance)
  value.p = cumP[ind]  # cumP less than chance
  LL = which(ind, arr.ind = TRUE)
  LL = as.vector(LL)
  value.T = LL[length(LL)]  # threshold value
  
  if (value.T > ULimit) {
    value.T = ULimit
  } else if (value.T < LLimit) {
    value.T = LLimit
  }
  
  drop(value.T)
}
threshold.value(1)