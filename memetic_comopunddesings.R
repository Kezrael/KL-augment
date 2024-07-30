# Algorithm parameters
total_pop_size <- 100
mutation_probability <- 0.2
crossover_probability <- 0.8
number_iterations <- 100
tour_size <- 2
ls_iterations <- 20

join <- 0.5
delete <- 0.02

# Design Parameters 
design_space <- c(1, 100)
min_points <- 3
max_points <- 7




# Antoine's Equation D- KL- optimality
# Antoine's Equation Water
# Initial Parameters
A <- 8.07131
B <- 1730.63
C <- 233.426

# Model Information
mean_antoine <- function(x) (10^(A - B/(C + x)))
model <- y ~ 10^(A-B/(C+x))
parameters <- c("A","B","C")
par_values <- c(A, B, C)
weight_fun <- function(x) (10^(A - B/(C + x)))^(-1)
var_fun_antoine <- function(x) (10^(A - B/(C + x)))^2

# Gradient vector for the information matrix
gradient <- function(model, char_vars, values, weight_fun = function(x) 1) {
  ext_char_vars <- c(char_vars, "x")
  arglist <- lapply(ext_char_vars, function(x) NULL)
  f <- as.function(append(stats::setNames(arglist, ext_char_vars), quote({})))
  f1 <- stats::deriv(model, char_vars, f)
  f2 <- function(x_val) {
    attr(do.call(f1, as.list(c(values, x_val))), "gradient")
  }
  f3 <- function(x){
    return(f2(x)*weight_fun(x))
  }
  return(f3)
}

# Information matrix function
inf_mat <- function(grad, points, weights) {
  matrix_ret <- 0 * diag(length(grad(points[[1]])))
  for (i in seq_along(weights)) {
    f_col <- as.matrix(grad(points[[i]]), nrow = 1, byrow = TRUE, dimnames = NULL)
    matrix_ret <- matrix_ret + (t(f_col) %*% f_col) * weights[[i]]
  }
  return(matrix_ret)
}

# Calculation of the gradient for this set of parameters and model
grad <- gradient(model, parameters, par_values, weight_fun)

# D-optimality function
D_opt_antoine_water <- function(points, weights){
  inf_mat <- inf_mat(grad, points, weights)
  det <- tryCatch({
    det(inf_mat)^(-1/3)
  }, error = function(e) {
    Inf
  })
  if(is.nan(det))
    det <- -Inf
  return(det)
}

#  Log Arithmetic - Log Geometric
arith_geom_means <- function(var_func, points, weights){
  k <- length(points)
  arith <- sum(map_dbl(1:k, function(i) do.call(var_func, list(points[i])) * weights[i]))
  geom <- prod(map_dbl(1:k, function(i) do.call(var_func, list(points[i])) ^ weights[i]))
  return(log(arith) - log(geom))
}

# KL- divergence
KL_antoine <- function(var_func, points, weights){
  return(1/2 * arith_geom_means(var_func, points, weights))
}

# KL-divergence for Antoine's Equation (eta^2)
KL_opt_antoine_water <- function(points, weights){
  KL_val <- KL_antoine(var_fun_antoine, points, weights)
  if(is.nan(KL_val))
    return(0)
  else
    return(KL_val)
}

# Define functions to minimise
function1 <- D_opt_antoine_water
function2 <- KL_opt_antoine_water

D_opt_val <- 6165.957
KL_opt_val <- 3.391391




lambda <- 0.50

  
# Define the function that calculates a linear combination of both metrics
fitness_comb <- function(points, weights){
  fitness <- lambda * D_opt_val / function1(points, weights) + (1 - lambda) * function2(points, weights) / KL_opt_val
  if(is.nan(fitness))
    return(-Inf)
  else 
    return(fitness)
}





### GENETIC ALGORITHM ####
memetic_alg <- function(lambda){
  
  fitness_comb <- function(points, weights){
    fitness <- lambda * D_opt_val / function1(points, weights) + (1 - lambda) * function2(points, weights) / KL_opt_val
    if(is.nan(fitness))
      return(-Inf)
    else 
      return(fitness)
  }
  
  # Generate initial population
  population <- generate_design_population(total_pop_size, min_points, max_points, design_space[1], design_space[2])
  
  # Calculate fitness values (KL- and D-criteria)
  population$Metrics <- fitness(population, fitness_comb)
  
  # # Calculate pareto front (first set of non-dominated solutions, second set, etc.)
  # population$Layers <- pareto_front(population)
  # 
  # # Calculate Crowding distance (manhattan distance to neighbours in same front)
  # population$CrowdingDistance <- crowding_distance(population) 
  
  # Start time of the algorithm
  start_time <- Sys.time()
  # Progress bar
  cli::cli_progress_bar("Evolving population", total = number_iterations)
  # Loop of generations: each iterations, one new generations and evolutionary selection is made
  for(i in 1:number_iterations){
    # Tournament selection (extract 100 parents from population after confronting a sample of 2 at a time)
    parents <- tournament_selection(population, total_pop_size, tour_size)
    parents_extract <- map(population , ~.[parents])
    
    # Generate the offsprings from the pairs of parents with probability 0.8
    sons <- crossover(parents_extract, crossover_probability)
    
    # Join offsprings and parents
    candidate_population <- map2(sons, population, c)
    
    # Generate random mutations to the population with probability 0.2
    candidate_population <- mutation(candidate_population, design_space[1], design_space[2], mutation_probability)
    
    # Calculate fitness values (KL- and D-criteria)
    candidate_population$Metrics <- fitness(candidate_population, fitness_comb)
    
    # Perform local search: mutate n times and keep mutated individual if fitness is dominant to the previous solution
    candidate_population <- local_search(candidate_population, design_space[1], design_space[2], ls_iterations, fitness_comb)
    
    # Heuristics to delete points with little weight and merge together similar points
    candidate_population <- heuristics(candidate_population, join, delete)
    
    # Calculate fitness values (KL- and D-criteria)
    candidate_population$Metrics <- fitness(candidate_population, fitness_comb)
    
    # Eliminate individuals with duplicated fitness
    candidate_population <- map(candidate_population, ~.[remove_duplicates(candidate_population)])
    
    # # Calculate pareto front (first set of non-dominated solutions, second set, etc.)
    # candidate_population$Layers <- pareto_front(candidate_population)
    # 
    # # Calculate Crowding distance (manhattan distance to neighbours in same front)
    # candidate_population$CrowdingDistance <- crowding_distance(candidate_population)
    
    # Select the best individuals for next generation
    population <- map(candidate_population, ~.[order(-candidate_population$Metrics)[1:total_pop_size]])
    
    # Progress bar update
    cli::cli_progress_update()
  }
  # Calculate total execution time
  end_time <- Sys.time()
  end_time - start_time
  
  # Final population sorted by metrics
  best_index <- which.max(population$Metrics)
  
  D_eff <- D_opt_val / D_opt_antoine_water(population$Points[[best_index]], population$Weights[[best_index]])
  
  KL_eff <- KL_opt_antoine_water(population$Points[[best_index]], population$Weights[[best_index]]) / KL_opt_val
  
  
  return(list("Point" = population$Points[[best_index]], "Weight" = population$Weights[[best_index]], "Deff" = D_eff, "KLeff" = KL_eff, "lambda" = lambda))
}

library(furrr)
library(purrr)

plan(multisession, workers = 3)

# compound_results <- map(seq(0, 1, 0.01), memetic_alg)

compound_results <- list()
for (i in 0:100){
  print(i)
  compound_results[[i+1]] <- memetic_alg(0.01 * i)
}

save(compound_results, file="../Documents/R/KL-opt/compound_results.RData")


df_eff_comp <- bind_rows(map(compound_results, ~.[c(3,4,5)]))

# Create a sorted dataframe with KL and D efficiencies
# df_eff <- data.frame(KL_eff = (1 / sapply(population$Metrics,"[[",2)) / kl_opt_val, D_eff = d_opt_val / sapply(population$Metrics,"[[",1), klcrit = sapply(population$Metrics,"[[",2), dcrit = sapply(population$Metrics,"[[",1)) %>% arrange(KL_eff)

# Add index to the solutions
# df_eff$Index <- 1:1000

# Plot colors
KL_color <- "#69b3a2"
D_color <- rgb(0.2, 0.6, 0.9, 1)

library(hrbrthemes)

df_eff_comp <- df_eff_comp %>% 
  arrange(desc(lambda))

which(df_eff_comp$KLeff > 0.8)
compound_results[47]

print(df_eff_comp, n = 60)

which.min(abs(df_eff_comp$KLeff - df_eff_comp$Deff))

compound_results[41]

library(latex2exp)

# Plot double axis efficiencies
ggplot(df_eff_comp, aes(x=lambda)) +
  geom_line( aes(y=KLeff), size=1, color=KL_color) + 
  geom_line( aes(y=Deff), size=1, color=D_color) +
  scale_y_continuous(
    # Features of the first axis
    name = "KL-efficiency",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~., name="D-efficiency")
  ) + 
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = KL_color, size=18),
    axis.title.y.right = element_text(color = D_color, size=18),
    axis.title.x = element_text(size=20),
    plot.title = element_text(size=22),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  ) +
  ggtitle("KL- vs D-efficiency of DKL-optimal designs") +
  xlab(TeX("$\\lambda$"))

df_eff$Index2 <- 1000:1


# Plot double axis efficiencies
ggplot(df_eff, aes(x=Index2)) +
  geom_line( aes(y=KL_eff), size=1, color=KL_color) + 
  geom_line( aes(y=D_eff), size=1, color=D_color) +
  scale_y_continuous(
    # Features of the first axis
    name = "KL-efficiency",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~., name="D-efficiency")
  ) + 
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = KL_color, size=13),
    axis.title.y.right = element_text(color = D_color, size=13)
  ) +
  ggtitle("KL- vs D-efficiency of Pareto Front")


df_eff %>% mutate(dif = -abs(KL_eff - D_eff)) %>% top_n(1, dif)
# df_eff %>% filter(KL_eff > 0.8) %>% top_n(1, D_eff)

out_lgl <- map_lgl(population$Metrics, ~all(near(., c(7366.899, 0.3514389), 0.01)))
which(out_lgl)

map(population, ~.[[210]])

df_eff %>% filter(KL_eff > 0.825) %>% top_n(2, D_eff)

out_lgl <- map_lgl(population$Metrics, ~all(near(., c(7225.980, 0.356726), 0.01)))
which(out_lgl)

map(population, ~.[[363]])

out_lgl <- map_lgl(population$Metrics, ~all(near(., c(6958.296, 0.3678842), 0.01)))
which(out_lgl)

map(population, ~.[[382]])



#################################################################
########### Auxiliary functions to run the algorithm ############
#################################################################

# Generate a random populations of "size" number of designs with a number of points in a range, and points within the design space
generate_design_population <- function(size, min_point_number, max_point_number, min_design_space, max_design_space){
  n_points <- sample(min_point_number:max_point_number, size, replace = TRUE)
  points <- lapply(n_points, function(x) sort(runif(x, min_design_space, max_design_space)))
  weights <- lapply(n_points, function(x) runif(x, 0, 1))
  weights <- lapply(weights, function(x) x / sum(x))
  return(list("Points" = points, "Weights" = weights))
}

# Calculate fitness for a population, given the two criterion function fun1 and fun2
fitness <- function(population, fun){
  # Calculations
  metrics <- map_dbl(1:length(population$Points), function(index) fun(population$Points[[index]], population$Weights[[index]]))
  
  return(metrics)
}

# Calculate all pareto fronts with package rPref to sort the solutions, and return the vector of front indexes
pareto_front <- function(population) {
  metrics_df <- as.data.frame(do.call(rbind, population$Metrics)) 
  pareto_front_levels <- metrics_df %>% rPref::psel(low(V1) | low(V2), top = nrow(metrics_df)) 
  
  return(merge(x = metrics_df, y = pareto_front_levels, by = c("V1", "V2"), all.x = TRUE, sort = FALSE)$.level)
}

# Calculate the crowding distance for the solutions in the different pareto fronts
crowding_distance <- function(population){
  population_size <- length(population$Weights)
  n_layers <- max(population$Layers)
  cd <- matrix(Inf, nrow = population_size, ncol = length(population$Metrics[[1]]))
  range <- vector(mode = "numeric", length(population$Metrics[[1]]))
  for (j in 1:length(population$Metrics[[1]])) {
    j_metric <- sapply(population$Metrics,"[[",j)
    range[[j]] <- max(j_metric[!is.infinite(j_metric)]) - min(j_metric[!is.infinite(j_metric)])
  }
  for (i in 1:n_layers) {
    indexes <- which(population$Layers == i)
    len <- length(indexes)
    if (len > 2) {
      for (j in 1:length(population$Metrics[[1]])) {
        # i <- 3
        j_metric <- sapply(population$Metrics,"[[",j)
        # range <- max(j_metric) - min(j_metric)
        j_metrics_layer <- j_metric[indexes]
        original_indexes <- indexes[order(j_metrics_layer)]
        cd[original_indexes[2:(len - 1)], j] <- abs(j_metric[original_indexes[3:len]] - j_metric[original_indexes[1:(len - 2)]]) / range[[j]]
      }
    }
  }
  return(rowSums(cd))
}

# Sample "pool_size" individuals, by sampling "tournament_size" individuals at a time and keeping the one with the best pareto front and crowding distance
tournament_selection <- function(population, pool_size, tournament_size){
  population_size <- length(population$Weights)
  winners <- purrr::map_int(1:pool_size, function(i) {
    candidates <- sample(population_size, tournament_size)
    # Revisar si positivo o negativo
    return(candidates[order(-population$Metrics[candidates])[1]])
  })
  return(winners)
}

# Given a population, randomise it, and generate offsprings by swaping between 1 to min_number_of_points among the designs
crossover <- function(population, cross_probability){
  population_size <- length(population$Weights)
  crossover_order <- sample(1:population_size)
  for(index in 1:(population_size/2)){
    if(runif(1) < cross_probability){
      max_point_swap <- min(length(population$Weights[[index]]), length(population$Weights[[index + population_size/2]]))
      n_point_swap <- sample(1:max_point_swap, 1)
      points_index_1 <- sample(1:length(population$Points[[index]]), n_point_swap)
      points_index_2 <- sample(1:length(population$Points[[index + population_size/2]]), n_point_swap)
      temp <- population$Points[[index]][points_index_1]
      population$Points[[index]][points_index_1] <- population$Points[[index + population_size/2]][points_index_2]
      population$Points[[index + population_size/2]][points_index_2] <- temp
    }
  }
  return(population)
}

# Mutate, aka produce some small change, in one weight or points of each individual with probability "mut_probability"
mutation <- function(population, min_design_space, max_design_space, mut_probability){
  indexes <- purrr::map_int(1:length(population$Weights), function(i) sample(1:length(population$Weights[[i]]), 1))
  mutate <- (1:length(population$Weights))[sample(c(FALSE, TRUE), size = length(population$Weights), replace = TRUE, prob = c(1 - mut_probability, mut_probability))]
  for(index in mutate){
    if(runif(1) < 0.5){
      new_weight <- (population$Weights[[index]][indexes[index]] + sample(c(-1, 1), 1) * 0.1 * sum(sample(0:1, 16, replace = TRUE, prob = c(15/16, 1/16)) * 2^seq(0, -15, -1)))
      if(new_weight < 0){
        new_weight <- 0
      } else if(new_weight > 1){
        new_weight <- 1
      }
      population$Weights[[index]][indexes[index]] <- new_weight
      population$Weights[[index]] <- population$Weights[[index]] / sum(population$Weights[[index]])
    } else{
      new_point <- (population$Points[[index]][indexes[index]] + sample(c(-1, 1), 1) * 0.1 * (max_design_space - min_design_space) * sum(sample(0:1, 16, replace = TRUE, prob = c(15/16, 1/16)) * 2^seq(0, -15, -1)))
      if(new_point < min_design_space){
        new_point <- min_design_space
      } else if(new_point > max_design_space){
        new_point <- max_design_space
      }
      population$Points[[index]][indexes[index]] <- new_point
    }
  }
  return(population)
}

# Produce "iterations" mutation that are kept only if the produce dominant solutions compared to the initial solution
local_search <- function(population, min_design_space, max_design_space, iterations, fitness){
  for(index in 1:length(population$Weights)){
    indexes <- purrr::map_int(1:length(population$Weights), function(i) sample(1:length(population$Weights[[index]]), 1))
    for(i in 1:iterations){
      if(runif(1) < 0.5){
        temp_weight <- population$Weights[[index]]
        temp_weight[indexes[index]] <- (population$Weights[[index]][indexes[index]] + sample(c(-1, 1), 1) * 0.1 * sum(sample(0:1, 16, replace = TRUE, prob = c(15/16, 1/16)) * 2^seq(0, -15, -1)))
        if(temp_weight[indexes[index]] < 0){
          temp_weight[indexes[index]] <- 0
        } else if(temp_weight[indexes[index]] > 1){
          temp_weight[indexes[index]] <- 1
        }
        temp_weight <- temp_weight / sum(temp_weight)
        if(population$Metrics[index] < fitness(population$Points[[index]], temp_weight)){
          population$Weights[[index]] <- temp_weight
          population$Metrics[index] <- fitness(population$Points[[index]], temp_weight)
        }
      } else{
        temp_point <- population$Points[[index]]
        temp_point[indexes[index]] <- (population$Points[[index]][indexes[index]] + sample(c(-1, 1), 1) * 0.1 * (max_design_space - min_design_space) * sum(sample(0:1, 16, replace = TRUE, prob = c(15/16, 1/16)) * 2^seq(0, -15, -1)))
        if(temp_point[indexes[index]] < min_design_space){
          temp_point[indexes[index]] <- min_design_space
        } else if(temp_point[indexes[index]] > max_design_space){
          temp_point[indexes[index]] <- max_design_space
        }
        if(fitness(population$Points[[index]], population$Weights[[index]]) < fitness(temp_point, population$Weights[[index]])){
          population$Points[[index]] <- temp_point
          population$Metrics[index] <- fitness(temp_point, population$Weights[[index]])
        }
      }
    }
  }
  return(population)
}

# Eliminate points with less than "delete_thresh" weight and join points closer than "join_thresh"
heuristics <- function(population, join_thresh, delete_thresh){
  for(index in 1:length(population$Points)){
    # Eliminate points with little weight
    to_delete <- population$Weights[[index]] <  delete_thresh
    population$Weights[[index]] <- population$Weights[[index]][!to_delete] / sum(population$Weights[[index]][!to_delete])
    population$Points[[index]] <- population$Points[[index]][!to_delete]
    
    # Order
    order <- order(population$Points[[index]])
    population$Points[[index]] <- population$Points[[index]][order]
    population$Weights[[index]] <- population$Weights[[index]][order]
    
    # Join similar points
    i <- 1
    while(i <= length(population$Points[[index]]) - 1) {
      if (population$Points[[index]][i + 1] - population$Points[[index]][i] < join_thresh) {
        added_point <- (population$Weights[[index]][i + 1] * population$Points[[index]][i + 1] +  population$Weights[[index]][i] * population$Points[[index]][i]) / (population$Weights[[index]][i + 1] + population$Weights[[index]][i])
        added_weight <- population$Weights[[index]][i + 1] + population$Weights[[index]][i]
        if(i == 1 && i == (length(population$Points[[index]]) - 1)){
          population$Points[[index]] <- c(added_point)
          population$Weights[[index]] <- c(added_weight)
          break
        }
        else if(i == 1){
          population$Points[[index]] <- c(added_point, population$Points[[index]][(i+2):length(population$Points[[index]])])
          population$Weights[[index]] <- c(added_weight, population$Weights[[index]][(i+2):length(population$Weights[[index]])])
        } else if (i == (length(population$Points[[index]]) - 1)) {
          population$Points[[index]] <- c(population$Points[[index]][1:(i-1)], added_point)
          population$Weights[[index]] <- c(population$Weights[[index]][1:(i-1)], added_weight)
        } else{
          population$Points[[index]] <- c(population$Points[[index]][1:(i-1)], added_point, population$Points[[index]][(i+2):length(population$Points[[index]])])
          population$Weights[[index]] <- c(population$Weights[[index]][1:(i-1)], added_weight, population$Weights[[index]][(i+2):length(population$Weights[[index]])])
        }
      } else {
        i <- i + 1
      }
    }
  }
  return(population)
}

# Remove points with the same fitness for both criterion to improve diversity 
remove_duplicates <- function(population){
  metrics_data <- data.frame(V1 = population$Metrics) %>% rowid_to_column("index") %>% arrange(V1)
  i <- 1
  while(i < nrow(metrics_data)){
    if(any(is.infinite(as.numeric(metrics_data[i,][c("V1")])))){
      metrics_data <- metrics_data[-i,]
    }
    else if(any(is.infinite(as.numeric(metrics_data[i+1,][c("V1")])))){
      metrics_data <- metrics_data[-(i+1),]
    } else if(all(near(metrics_data[i,][c("V1")], metrics_data[i+1,][c("V1")], tol = 0.000001))){
      metrics_data <- metrics_data[-i,]
    } else{
      i <- i + 1
    }
  }
  return(metrics_data$index)
}
