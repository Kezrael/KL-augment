
#   This simulation does the following:

#   ** Find the LRT distribution via bootstrap estimation under the null.

# ** Simulate $n=75$ observations from the homoskedastic and heteroskedatic models, using the two chosen designs. So we have:

#   ** Design1 + homoskedastic

# ** Design2 + homoskedastic

# ** Design1 + heteroskedastic

# ** Design2 + heteroskedastic

# ** For each of the 4 cases, test the hypothesis:
#   $$H_0: \mbox{ homoskedastic} \quad vs \quad  H_1: \mbox{ heteroskedastic}$$

#   The p-value is computed as $P(Reject \mid H_0)$, so, using the distribution under the homoskedastic case.


## Initialization: useful functions and libraries, initialization of the parameters and definition of the designs



rm(list=ls())

a <- 8.07131
b <- 1730.63
c <- 233.426

n <- 75

library(stats4)
library(ggplot2)

simulate_hom <- function(x, sigma){
  return(10^(a-b/(c+x)) + rnorm(1, mean = 0, sd = sigma))
}

simulate_het <- function(x, sigma){
  return(10^(a-b/(c+x)) + rnorm(1, mean = 0, sd = sigma*10^(a-b/(c+x))))
}

design1 <- data.frame("Point" = c(1, 1.63, 41.88, 100), "Weight" = c(0.2, 0.4, 0.2, 0.2))
design2 <- data.frame("Point" = c(1, 41.88, 100), "Weight" = c(0.76, 0.08, 0.16))

exact_design1 <- rep(design1$Point, times = n * design1$Weight)
exact_design2 <- rep(design2$Point, times = n * design2$Weight)

likelihood_hom <- function(x, y, a, b, c, sigma){
  1/sqrt(2*pi*sigma^2)*exp(-(y-10^(a-b/(c+x))^2/(2*sigma^2)))
}

loglikelihood_hom <- function(x, y, a, b, c, sigma){
  -log(sqrt(2*pi*sigma^2))-(y-10^(a-b/(c+x))^2/(2*sigma^2))
}

likelihood_het <- function(x, a, b, c, sigma){
  1/sqrt(2*pi*(sigma*10^(a-b/(c+x)))^2)*exp(-(y-10^(a-b/(c+x))^2/(2*(sigma*10^(a-b/(c+x)))^2)))
}

loglikelihood_het <- function(x, a, b, c, sigma){
  -log(sqrt(2*pi*(sigma*10^(a-b/(c+x)))^2)) -(y-10^(a-b/(c+x))^2/(2*(sigma*10^(a-b/(c+x)))^2))
}



# The output of each of the bootstrap distributions give the 'N' bootstrap values of the LRT generated under each of the 2 scenarios (designs) for the null hypothesis, homoskedasticity.

# The test should reject for high values of $logLik(result1) - logLik(result0)$, where 'result0' indicates the loglikelihood under the null (homoskedasticity) and 'result1' the loglikelihood under the alternative (heterosk).

# The bootstrap resamples are fixed here: 

N=5000 # change to a bigger number



## 1. Distribution of the LRT, case 1: design1 + H0


#### Design 1
# Simulate H0 hom
set.seed(42)

sigma <- 2

r <- vector(mode = "double", length = N)

start_time <- Sys.time()
# Progress bar
i <- 1
cli::cli_progress_bar("Computing R distribution", total = N)
while(i <= N){
  sim_design <- purrr::map_dbl(exact_design1, ~simulate_hom(., sigma))
  log_likehom <- function(a, b, c){
    -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*sigma^2))-(sim_design[.]-10^(a-b/(c+exact_design1[.])))^2/(2*sigma^2)))
  }
  log_likehet <- function(a, b, c){
    -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*(sigma*10^(a-b/(c+exact_design1[.])))^2))-(sim_design[.]-10^(a-b/(c+exact_design1[.])))^2/(2*(sigma*10^(a-b/(c+exact_design1[.])))^2)))
  }
  try({
    result0 <- mle(minuslogl = log_likehom, start = list(a = 8.07131, b = 1730.63, c = 233.426))
    result1 <- mle(minuslogl = log_likehet, start = list(a = 8.07131, b = 1730.63, c = 233.426))
   
    r[i] <- as.numeric(logLik(result1) - logLik(result0))
    i <- i+1
    cli::cli_progress_update()
  })
}
# Calculate total execution time
end_time <- Sys.time()
end_time - start_time

quant95 <- quantile(r, probs = 0.95)

ggplot(mapping = aes(x=r)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("R statistic distribution homoscedastic") +
  theme_bw() +
  geom_vline(xintercept = quant95, color = "firebrick4")



## 2. Distribution of the LRT, case 2: design2 + H0

set.seed(84)
sigma <- 2

rr <- vector(mode = "double", length = N)

start_time <- Sys.time()
# Progress bar
cli::cli_progress_bar("Computing R distribution", total = N)
for(i in 1:N){
  sim_design <- purrr::map_dbl(exact_design2, ~simulate_hom(., sigma))
  log_likehom <- function(a, b, c){
    -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*sigma^2))-(sim_design[.]-10^(a-b/(c+exact_design2[.])))^2/(2*sigma^2)))
  }
  log_likehet <- function(a, b, c){
    -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*(sigma*10^(a-b/(c+exact_design2[.])))^2))-(sim_design[.]-10^(a-b/(c+exact_design2[.])))^2/(2*(sigma*10^(a-b/(c+exact_design2[.])))^2)))
  }

  try({
    result0 <- mle(minuslogl = log_likehom, start = list(a = 8.07131, b = 1730.63, c = 233.426))
    result1 <- mle(minuslogl = log_likehet, start = list(a = 8.07131, b = 1730.63, c = 233.426))
    
    rr[i] <- as.numeric(logLik(result1) - logLik(result0))
    i <- i+1
    cli::cli_progress_update()
  })
}
# Calculate total execution time
end_time <- Sys.time()
end_time - start_time

quant95_2 <- quantile(rr, probs = 0.95)

ggplot(mapping = aes(x=rr)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("R statistic distribution homoscedastic") +
  theme_bw() +
  geom_vline(xintercept = quant95_2, color = "firebrick4")




# Save in a table the 2 distributions


distributions <- data.frame("des1_hom" = r, "des2_hom" = rr)

readr::write_csv(distributions, "distributions_2.csv")

# distributions <- readr::read_csv("C:/Users/drake/Documents/R/KL-sim/distributions.csv")

## Simulation of data under H0 from design 1

# In this case, the p-value **should be higher than $\alpha$ **, because H0 is true


set.seed(1)
sigma <- 2
sim_design1_truehom <- purrr::map_dbl(exact_design1, ~simulate_hom(., sigma))

log_likehom <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*sigma^2))-(sim_design1_truehom[.]-10^(a-b/(c+exact_design1[.])))^2/(2*sigma^2)))
}
log_likehet <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*(sigma*10^(a-b/(c+exact_design1[.])))^2))-(sim_design1_truehom[.]-10^(a-b/(c+exact_design1[.])))^2/(2*(sigma*10^(a-b/(c+exact_design1[.])))^2)))
}
result0 <- mle(minuslogl = log_likehom, start = list(a = 8.07131, b = 1730.63, c = 233.426))
result1 <- mle(minuslogl = log_likehet, start = list(a = 8.07131, b = 1730.63, c = 233.426))

p_val <- sum(distributions$des1_hom > as.numeric(logLik(result1) - logLik(result0)))/N
print(p_val)

if(p_val > 0.05){
  print(attr(result0, "coef"))
  print(vcov(result0))
  print(sqrt(diag(vcov(result0))))
} else{
  print(attr(result1, "coef"))
  print(vcov(result1))
  print(sqrt(diag(vcov(result1))))
}



## Simulation of data under H1 from design 1

# In this case, the p-value **should be low**, because H0 is false.



set.seed(2)
# Design1 het sim
sigma <- 0.1
sim_design1_truehet <- purrr::map_dbl(exact_design1, ~simulate_het(., sigma))

log_likehom <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*sigma^2))-(sim_design1_truehet[.]-10^(a-b/(c+exact_design1[.])))^2/(2*sigma^2)))
}
log_likehet <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*(sigma*10^(a-b/(c+exact_design1[.])))^2))-(sim_design1_truehet[.]-10^(a-b/(c+exact_design1[.])))^2/(2*(sigma*10^(a-b/(c+exact_design1[.])))^2)))
}
result0 <- mle(minuslogl = log_likehom, start = list(a = 8.07131, b = 1730.63, c = 233.426))
result1 <- mle(minuslogl = log_likehet, start = list(a = 8.07131, b = 1730.63, c = 233.426))

p_val <- sum(distributions$des1_hom > as.numeric(logLik(result1) - logLik(result0)))/N
print(p_val)

if(p_val > 0.05){
  print(attr(result0, "coef"))
  print(vcov(result0))
  print(sqrt(diag(vcov(result0))))
} else{
  print(attr(result1, "coef"))
  print(vcov(result1))
  print(sqrt(diag(vcov(result1))))
}



## Simulation of data under H0 from design 2

# In this case, the p-value **should be higher than $\alpha$ **, because H0 is true


set.seed(3)
# Design2 hom sim
sigma <- 2
sim_design2_truehom <- purrr::map_dbl(exact_design2, ~simulate_hom(., sigma))

log_likehom <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*sigma^2))-(sim_design2_truehom[.]-10^(a-b/(c+exact_design2[.])))^2/(2*sigma^2)))
}
log_likehet <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*(sigma*10^(a-b/(c+exact_design2[.])))^2))-(sim_design2_truehom[.]-10^(a-b/(c+exact_design2[.])))^2/(2*(sigma*10^(a-b/(c+exact_design2[.])))^2)))
}
result0 <- mle(minuslogl = log_likehom, start = list(a = 8.07131, b = 1730.63, c = 233.426))
result1 <- mle(minuslogl = log_likehet, start = list(a = 8.07131, b = 1730.63, c = 233.426))

p_val <- sum(distributions$des2_hom > as.numeric(logLik(result1) - logLik(result0)))/N
print(p_val)

if(p_val > 0.05){
  print(attr(result0, "coef"))
  print(vcov(result0))
  print(sqrt(diag(vcov(result0))))
} else{
  print(attr(result1, "coef"))
  print(vcov(result1))
  print(sqrt(diag(vcov(result1))))
}




## Simulation of data under H1 from design 2

# In this case, the p-value **should be low**, because H0 is false.


set.seed(4)
# Design1 het sim
sigma <- 0.1
sim_design2_truehet <- purrr::map_dbl(exact_design2, ~simulate_het(., sigma))

log_likehom <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*sigma^2))-(sim_design2_truehet[.]-10^(a-b/(c+exact_design2[.])))^2/(2*sigma^2)))
}
log_likehet <- function(a, b, c){
  -sum(purrr::map_dbl(1:n, ~ -log(sqrt(2*pi*(sigma*10^(a-b/(c+exact_design2[.])))^2))-(sim_design2_truehet[.]-10^(a-b/(c+exact_design2[.])))^2/(2*(sigma*10^(a-b/(c+exact_design2[.])))^2)))
}
result0 <- mle(minuslogl = log_likehom, start = list(a = 8.07131, b = 1730.63, c = 233.426))
result1 <- mle(minuslogl = log_likehet, start = list(a = 8.07131, b = 1730.63, c = 233.426))

p_val <- sum(distributions$des2_hom > as.numeric(logLik(result1) - logLik(result0)))/N
print(p_val)

if(p_val > 0.05){
  print(attr(result0, "coef"))
  print(vcov(result0))
  print(sqrt(diag(vcov(result0))))
} else{
  print(attr(result1, "coef"))
  print(vcov(result1))
  print(sqrt(diag(vcov(result1))))
}






