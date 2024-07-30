library(tidyverse)

kl_opt_des <- function(var_func, des_space, grid.length = 10000){
  x_val <- seq(des_space[1], des_space[2], length.out = grid.length)
  y_val <- map_dbl(x_val, var_func)
  h_max <- max(y_val)
  x_max <- x_val[which.max(y_val)]
  h_min <- min(y_val)
  x_min <- x_val[which.min(y_val)]
  omega <- h_max/(h_max - h_min) - 1 / (log(h_max) - log(h_min))
  return(data.frame("Point" = c(x_min, x_max), "Weight" = c(omega, 1 - omega) ))
}

kl_opt_des_multivar <- function(mean_func, des_space){
  
}

arith_geom_means <- function(var_func, design){
  k <- dim(design)[1]
  arith <- sum(map_dbl(1:k, function(i) do.call(var_func, list(design$Point[i])) * design$Weight[i]))
  geom <- prod(map_dbl(1:k, function(i) do.call(var_func, list(design$Point[i])) ^ design$Weight[i]))
  return(log(arith) - log(geom))
}
  
KL_antoine <- function(var_func, design){
  return(1/2 * arith_geom_means(var_func, design))
}

sum_design_point <- function(design, x, alpha){
  design$Weight <- design$Weight * (1 - alpha)
  design[nrow(design)+1, ] <- c(x, alpha)
  return(design)
}

get_KL_augment_data <- function(KL_divergence, des_space, init_design, alpha, grid.length = 10000){
  x_val <- seq(des_space[1], des_space[2], length.out = grid.length)
  y_val <- map_dbl(x_val, function(x) do.call(KL_divergence, list(sum_design_point(init_design, x, alpha))))
  return(data.frame("x" = x_val, "KL_value" = y_val))
}

get_KL_augment_val <- function(KL_divergence, init_design, alpha, x){
  return(do.call(KL_divergence, list(sum_design_point(init_design, x, alpha))))
}

# Antoine's Equationb
a <- 8.07131
b <- 1730.63
c <- 233.426
design_space <- c(1, 100)
mean_antoine <- function(x) (10^(a - b/(c + x)))
var_antoine <- function(x) (10^(a - b/(c + x)))^2
weight_antoine <- function(x) (10^(a - b/(c + x)))^(-1)

kl_optdes <- kl_opt_des(var_antoine, c(1, 100))

arith_geom_means(var_antoine, kl_optdes)
KL_opt_val <- KL_antoine(var_antoine, kl_optdes)

KL_antoine_water <- function(design){
  KL_antoine(var_antoine, design)
}

library(optedr)

resAnt.D <-opt_des("D-Optimality", y ~ 10^(a-b/(c+x)),
                   c("a","b","c"), c(8.07131,  1730.63, 233.426),
                   c(1, 100), weight_fun = weight_antoine)
resAnt.D$crit_value

KL_antoine_water(resAnt.D$optdes)
KL_antoine_water(kl_optdes)

alpha <- 0.25

eff_plot_data <- get_KL_augment_data(KL_antoine_water, design_space, kl_optdes, alpha)

eff_plot_data$eff <- eff_plot_data$KL_value / KL_opt_val

library(latex2exp)

ggplot() +
  geom_line(data = eff_plot_data, mapping = aes(x = x, y = eff), color = "steelblue") +
  ggplot2::xlim(design_space[[1]], design_space[[2]]) +
  ggplot2::labs(x = "x", y = "Efficiency") +
  ggplot2::theme_bw() +
  theme(
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  ) +
  ylab(TeX("eff$_{KL}(\\xi_{k+1}^{(x)})$"))


kl_near_opt <- sum_design_point(kl_optdes, 50, 0.01)

KL_nearopt_val <- KL_antoine(var_antoine, kl_near_opt)

eff_plot_data_n <- get_KL_augment_data(KL_antoine_water, design_space, kl_near_opt, alpha)

eff_plot_data_n$eff <- eff_plot_data_n$KL_value / KL_nearopt_val

ggplot() +
  geom_line(data = eff_plot_data_n, mapping = aes(x = x, y = eff), color = "steelblue") +
  ggplot2::xlim(design_space[[1]], design_space[[2]]) +
  ggplot2::labs(x = "x", y = "Efficiency") +
  ggplot2::theme_bw()


update_sequence <- function(points, tol){
  i <- 1
  imax <- (length(points)-1)
  while(i < imax) {
    absdiff <- c(rep(T, i), abs(points[-seq(1, i)] - points[i]) > tol)
    points <- points[absdiff]
    i <- i+1
    imax <- (length(points)-1)
  }
  return(points)
}


crosspoints <- function(val, sens, gridlength, tol, xmin, xmax){
  
  sensfix <- function(x){
    return(sens(x) - val)
  }
  
  sols <- vector(mode = "numeric", length = gridlength)
  cli::cli_progress_bar("Calculating regions", total = gridlength)
  startsx <- seq(xmin, xmax, length.out = gridlength)
  for(i in 1:gridlength){
    sols[i] <- nleqslv::nleqslv(startsx[i], fn = sensfix)$x
    cli::cli_progress_update()
  }
  
  # sols <- unlist(purrr::map(purrr::map(cli::cli_progress_along(seq(xmin, xmax, length.out = gridlength), name = "Calculating regions"), nleqslv::nleqslv, fn = sensfix), function(x) x$x))
  
  # Eliminar duplicados
  sols_upd <- update_sequence(sols, tol)
  
  # Quedarnos con puntos dentro del espacio de diseño
  sols_upd <- sols_upd[sols_upd >= xmin & sols_upd <= xmax]
  
  # Eliminar los que no son soluciones
  sols_upd <- sols_upd[abs(purrr::map_dbl(sols_upd, sens) - val) < tol]
  
  
  return(sort(sols_upd))
}

getStart <- function(cross, min, max, val, sens_opt){
  if(length(cross)==1){
    if(round(cross[1]-min, 6)==0){
      if((val < sens_opt((max+cross[1])/2))){
        start <- F
      }
      else{
        start <- T
      }
    }
    else{
      if((val < sens_opt((min+cross[1])/2))){
        start <- T
      }
      else{
        start <- F
      }
    }
  }
  else if(val < sens_opt((cross[2]+cross[1])/2)){
    start <- F
  }
  else {
    start <- T
  }
  return(start)
}

getPar <- function(cross){
  return(length(cross)%%2 == 0)
}

getCross2 <- function(cross, min, max, start, par){
  if(par & start){
    cross2 <- c(min, cross, max)
  }
  else if(par & !start){
    cross2 <- cross
  }
  else if(!par & start){
    cross2 <- c(min, cross)
  }
  else{
    cross2 <- c(cross, max)
  }
  return(cross2)
}


crossregion <- function(kl_div_func, init_design, alpha, eff, design_space, init_KL){
  kl_eff_func <- function(x){
    get_KL_augment_val(kl_div_func, init_design, alpha, x) / init_KL
  }
    
  cross <- sort(crosspoints(eff, kl_eff_func, 10000, 10^(-3), design_space[[1]], design_space[[2]]))
  start <- getStart(cross, design_space[[1]], design_space[[2]], eff, kl_eff_func)
  par <- getPar(cross)
  cand_points_reg <- getCross2(cross, design_space[[1]], design_space[[2]], start, par)
  
  x_val <- seq(design_space[[1]], design_space[[2]], length.out = 10000)
  # eff <- function(x){
  #   return((1 - alpha) * (1 + alpha * sens_1(x)/(1 - alpha))^(1 / length(parameters)))
  # }
  y_val <- purrr::map_dbl(x_val, kl_eff_func)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(mapping = ggplot2::aes(x = x_val, y = y_val), color = "steelblue3") +
    ggplot2::geom_hline(yintercept =  kl_eff_func(cross[1]), color = "goldenrod3") +
    ggplot2::xlim(design_space[[1]], design_space[[2]]) +
    ggplot2::labs(x = "x") +
    ylab(TeX("eff$_{KL}(\\xi_{k+1}^{(x)})$")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
  
  efficiency <- kl_eff_func(cross[1])
  for(i in 1:(length(cand_points_reg)/2)){
    loop_input = paste("ggplot2::geom_segment(ggplot2::aes(x=",cand_points_reg[2*i-1],",xend=",cand_points_reg[2*i],",y=efficiency,yend=efficiency), size = 1.5, color = 'green3')", sep="")
    p <- p + eval(parse(text=loop_input))
  }
  
  values <- purrr::map_dbl(cross, kl_eff_func)
  cutoffpoints <- data.frame("x_value" = cross, "eff" = values)
  
  p <- p + ggplot2::geom_point(data = cutoffpoints, ggplot2::aes(x = x_value, y = eff), shape = 16, size = 2, color = "firebrick3")
  
  # p <- p + ggplot2::geom_segment(ggplot2::aes(x = design_space[[1]], xend = design_space[[1]], y = delta_range[1], yend = delta_range[2]), col = "mediumpurple2", size = 1.5)
  
  
  return(list("region" = cand_points_reg, "plot" = p))
}

efficiency <- 0.8

# crossregion(KL_antoine_water, kl_optdes, alpha, efficiency, design_space, KL_opt_val)

antoine_regions <- crossregion(KL_antoine_water, kl_optdes, alpha, efficiency, design_space, KL_opt_val)


antoine_regions$region

antoine_regions$plot +
  theme(
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  ) 


















#####################################
######## KL Wynn-Fedorov ############
#####################################

arith_geom_means <- function(var_func, design){
  k <- dim(design)[1]
  arith <- sum(map_dbl(1:k, function(i) do.call(var_func, list(design$Point[i])) * design$Weight[i]))
  geom <- prod(map_dbl(1:k, function(i) do.call(var_func, list(design$Point[i])) ^ design$Weight[i]))
  return(log(arith) - log(geom))
}

KL_antoine <- function(var_func, design){
  return(1/2 * arith_geom_means(var_func, design))
}

sum_design_point <- function(design, x, alpha){
  design$Weight <- design$Weight * (1 - alpha)
  design[nrow(design)+1, ] <- c(x, alpha)
  return(design)
}

get_KL_augment_data <- function(KL_divergence, des_space, init_design, alpha, grid.length = 10000){
  x_val <- seq(des_space[1], des_space[2], length.out = grid.length)
  y_val <- map_dbl(x_val, function(x) do.call(KL_divergence, list(sum_design_point(init_design, x, alpha))))
  return(data.frame("x" = x_val, "KL_value" = y_val))
}

get_KL_augment_val <- function(KL_divergence, init_design, alpha, x){
  return(do.call(KL_divergence, list(sum_design_point(init_design, x, alpha))))
}

# Antoine's Equation
a <- 8.07131
b <- 1730.63
c <- 233.426
design_space <- c(1, 100)
mean_antoine <- function(x) (10^(a - b/(c + x)))
var_antoine <- function(x) (10^(a - b/(c + x)))^2


delete_points <- function(design, delta) {
  updatedDesign <- design[design$Weight > delta, ]
  updatedDesign$Weight <- updatedDesign$Weight / sum(updatedDesign$Weight)
  return(updatedDesign)
}

update_design <- function(design, xmax, delta, new_weight) {
  absdiff <- abs(design$Point - xmax) < delta
  design$Weight <- design$Weight * (1 - new_weight)
  if (any(absdiff)) {
    pos <- min(which(absdiff == TRUE))
    design$Point[[pos]] <- (design$Point[[pos]] + xmax) / 2
    design$Weight[[pos]] <- design$Weight[[pos]] + new_weight
  }
  else {
    design[nrow(design) + 1, ] <- c(xmax, new_weight)
  }
  # design$Weight <- rep(1 / nrow(design), nrow(design))
  return(design)
}

KL_opt_des <- function(KL_div_func, des_space, init_design, n_pars, delete_thresh, join_thresh, tol = 10^(-4), init_pars = NULL, grid.length = 10000, max_iter = 200){
  if(is.null(init_pars)){
    init_pars <- rep(0, n_pars)
  }
  x_grid <- seq(des_space[1], des_space[2], length.out = grid.length)
  KL_values <- vector(mode = "numeric")
  
  for(i in 1:max_iter){
    alpha <- 1/(3 + i)
    opt_theta2 <- function(pars){
      return(sum(map_dbl(1:nrow(init_design), function(i) KL_div_func(pars, init_design$Point[i]) * init_design$Weight[i])))
    }
    opt_pars <- optim(opt_theta2, init_pars)$par
    
    KL_values <- c(KL_values, opt_theta2(opt_pars))
    
    sens_val <- map_dbl(x_grid, function(x) KL_div_func(opt_pars, x))
    x_add <- x_grid[which.max(sens_val)]
    
    init_design <- update_design(init_design, x_add, join_thresh, alpha) 
    
    if(i %% 10 == 0){
      init_design <- delete_points(init_design, delete_thresh)
    }
  }
  
  return(list("optdes" = init_design, "convergence" = KL_values))
}





























#########################################################
############# D-optimality split des_space ##############
#########################################################

DWFMult <- function(init_design, grad, min, max, x.grid, join_thresh, delete_thresh, k, delta_weights, tol, tol2) {
  Point <- NULL
  crit_val <- numeric(2122)
  index <- 1
  # Maximum iterations for the optimize weights loop
  maxiter <- 100
  cli::cli_progress_bar("Calculating optimal design")
  for (i in 1:21) {
    cli::cli_progress_update()
    M <- inf_mat(grad, init_design)
    crit_val[index] <- dcrit(M, k)
    index <- index + 1
    sensM <- dsens(grad, M)
    xmax <- findmax(sensM, x.grid)
    if ((sensM(xmax) - k) / k < tol2) {
      message("\n", crayon::blue(cli::symbol$info), " Stop condition reached: difference between sensitivity and criterion < ", tol2)
      break
    }
    init_design <- update_design(init_design, xmax, join_thresh, 1/(index + 2))
    iter <- 1
    stopw <- FALSE
    while (!stopw) {
      weightsInit <- init_design$Weight
      M <- inf_mat(grad, init_design)
      crit_val[index] <- dcrit(M, k)
      index <- index + 1
      sensM <- dsens(grad, M)
      init_design$Weight <- update_weights(init_design, sensM, k, delta_weights)
      stopw <- any(max(abs(weightsInit - init_design$Weight) < tol)) || iter >= maxiter
      iter <- iter + 1
    }
    init_design <- delete_points(init_design, delete_thresh)
    if (i %% 5 == 0) {
      init_design <- update_design_total(init_design, join_thresh)
    }
    if (i == 21) {
      message("\n", crayon::blue(cli::symbol$info), " Stop condition not reached, max iterations performed")
    }
  }
  cli::cli_progress_update(force = TRUE)
  base::cat("")
  crit_val[index] <- dcrit(M, k)
  crit_val <- crit_val[1:(length(crit_val) - sum(crit_val == 0))]
  conv <- data.frame("criteria" = crit_val, "step" = seq(1, length(crit_val), 1))
  conv_plot <- plot_convergence(conv)
  
  init_design <- dplyr::arrange(init_design, Point)
  rownames(init_design) <- NULL
  
  M <- inf_mat(grad, init_design)
  sensM <- dsens(grad, M)
  xmax <- findmax(sensM, x.grid)
  
  atwood <- k / sensM(xmax) * 100
  
  message(crayon::blue(cli::symbol$info), " The lower bound for efficiency is ", atwood, "%")
  
  plot_opt <- plot_sens(min, max, x.grid, sensM, k)
  l_return <- list(
    "optdes" = init_design, "convergence" = conv_plot,
    "sens" = plot_opt, "criterion" = "D-Optimality", "crit_value" = crit_val[length(crit_val)]
  )
  attr(l_return, "hidden_value") <- k
  attr(l_return, "gradient") <- grad
  attr(l_return, "atwood") <- atwood
  class(l_return) <- "optdes"
  l_return
}

a <- 8.07131
b <- 1730.63
c <- 233.426
weight_antoine <- function(x) (10^(a - b/(c + x)))^(-1)
resAnt.D <-opt_des("D-Optimality", y ~ 10^(a-b/(c+x)),
                   c("a","b","c"), c(8.07131,  1730.63, 233.426),
                   c(1, 100))


x_grid <- split_grid(antoine_regions$region)
model <- y ~ 10^(a-b/(c+x))
parameters <- c("a","b","c")
par_values <- c(a, b, c)
min_val <- min(x_grid)
max_val <- max(x_grid)
init_design <- NULL
join_thresh <- (max_val - min_val) / 20 
delete_thresh <- 0.02
delta <- 1 / 2
tol <- 0.00001
tol2 <- 0.00001
weight_fun <- function(x) (10^(a - b/(c + x)))^(-1)
k <- length(parameters)
grad <- gradient(model, parameters, par_values, weight_fun)

points_index <- vector(mode = "numeric")
for(i in 1:(length(x_grid) %/% 1000)){
  points_index <- c(points_index, round(seq(1, 1000, length.out = k * (k + 1) / 2 + 1)) + (i - 1) * 1000)
}
init_design <- data.frame("Point" = x_grid[points_index], "Weight" = rep(1 / length(points_index), times = length(points_index)))

test <- DWFMult(init_design, grad, min_val, max_val, x_grid, join_thresh, delete_thresh, k, delta, tol, tol2) 

test$sens +
  theme(
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  ) 
test$optdes

ggplot2::ggplot(data = test$data) +
  ggplot2::theme_bw() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y, group = grp), color = "steelblue3") +
  ggplot2::stat_function(fun = function(x) 3, col = "goldenrod3") +
  ggplot2::xlim(1, 100) +
  ggplot2::labs(x = "X", y = "Y")
sensibility


#########################################################
############ D-optimality restriced design ##############
#########################################################
DWFMultRegularised <- function(fixed_design, init_design, lambda, grad, min, max, x.grid, join_thresh, delete_thresh, k, delta_weights, tol, tol2) {
  Point <- NULL
  crit_val <- numeric(2122)
  index <- 1
  # Maximum iterations for the optimize weights loop
  maxiter <- 100
  cli::cli_progress_bar("Calculating optimal design")
  for (i in 1:401) {
    cli::cli_progress_update()
    design <- add_design(fixed_design, init_design, lambda)
    M <- inf_mat(grad, design)
    crit_val[index] <- dcrit(M, k)
    index <- index + 1
    sensM <- dsens(grad, M)
    xmax <- findmax(sensM, x_grid)
    if ((sensM(xmax) - k) / k < tol2) {
      message("\n", crayon::blue(cli::symbol$info), " Stop condition reached: difference between sensitivity and criterion < ", tol2)
      break
    }
    init_design <- update_design(init_design, xmax, join_thresh, 1/(index + 2))
    # iter <- 1
    # stopw <- FALSE
    # while (!stopw) {
    #   design <- add_design(fixed_design, init_design, lambda)
    #   weightsInit <- init_design$Weight
    #   M <- inf_mat(grad, design)
    #   crit_val[index] <- dcrit(M, k)
    #   index <- index + 1
    #   sensM <- dsens(grad, M)
    #   init_design$Weight <- update_weights(init_design, sensM, tr(solve(M) %*% inf_mat(grad, init_design)), delta_weights)
    #   stopw <- any(max(abs(weightsInit - init_design$Weight) < tol)) || iter >= maxiter
    #   iter <- iter + 1
    # }
    init_design <- delete_points(init_design, delete_thresh)
    if (i %% 5 == 0) {
      init_design <- update_design_total(init_design, join_thresh)
    }
    if (i == 21) {
      message("\n", crayon::blue(cli::symbol$info), " Stop condition not reached, max iterations performed")
    }
  }
  cli::cli_progress_update(force = TRUE)
  base::cat("")
  crit_val[index] <- dcrit(M, k)
  crit_val <- crit_val[1:(length(crit_val) - sum(crit_val == 0))]
  conv <- data.frame("criteria" = crit_val, "step" = seq(1, length(crit_val), 1))
  conv_plot <- plot_convergence(conv)
  
  design <- add_design(fixed_design, init_design, lambda)
  design <- update_design_total(design, join_thresh)
  design <- dplyr::arrange(design, Point)
  rownames(design) <- NULL
  
  M <- inf_mat(grad, design)
  sensM <- dsens(grad, M)
  xmax <- findmax(sensM, x.grid)
  
  atwood <-  tr(solve(M) %*% inf_mat(grad, init_design)) / sensM(xmax) * 100
  
  message(crayon::blue(cli::symbol$info), " The lower bound for efficiency is ", atwood, "%")
  
  plot_opt <- plot_sens(min, max, x.grid, sensM,  tr(solve(M) %*% inf_mat(grad, init_design)))
  plot_opt
  l_return <- list(
    "optdes" = design, "convergence" = conv_plot,
    "sens" = plot_opt, "criterion" = "D-Optimality", "crit_value" = crit_val[length(crit_val)]
  )
  attr(l_return, "hidden_value") <- tr(solve(M) %*% inf_mat(grad, init_design))
  attr(l_return, "gradient") <- grad
  attr(l_return, "atwood") <- atwood
  class(l_return) <- "optdes"
  l_return
}

a <- 8.07131
b <- 1730.63
c <- 233.426
weight_fun <- function(x) (10^(a - b/(c + x)))^(-1)
resAnt.D <-opt_des("D-Optimality", y ~ 10^(a-b/(c+x)),
                   c("a","b","c"), c(8.07131,  1730.63, 233.426),
                   c(1, 100))


x_grid <- split_grid(antoine_regions$region)
model <- y ~ 10^(a-b/(c+x))
parameters <- c("a","b","c")
par_values <- c(a, b, c)
min_val <- min(x_grid)
max_val <- max(x_grid)
init_design <- NULL
join_thresh <- (max_val - min_val) / 20 
delete_thresh <- 0.02
delta <- 1 / 2
tol <- 0.00001
tol2 <- 0.00001
weight_fun <- function(x) (10^(a - b/(c + x)))^(-1)
k <- length(parameters)
grad <- gradient(model, parameters, par_values, weight_fun)

points_index <- vector(mode = "numeric")
for(i in 1:(length(x_grid) %/% 1000)){
  points_index <- c(points_index, round(seq(1, 1000, length.out = k * (k + 1) / 2 + 1)) + (i - 1) * 1000)
}
init_design <- data.frame("Point" = x_grid[points_index], "Weight" = rep(1 / length(points_index), times = length(points_index)))

start_time <- Sys.time()
test_reg <- DWFMultRegularised(kl_optdes, init_design, 0.25, grad, min_val, max_val, x_grid, join_thresh, delete_thresh, k, delta, tol, tol2) 
end_time <- Sys.time()
end_time - start_time

test_reg$sens +
  theme(
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  ) 
test_reg$optdes

opt_plus_kl <- add_design(kl_optdes, test$optdes, 0.25)
opt_plus_kl <- update_design_total(opt_plus_kl, join_thresh)

dcrit(inf_mat(grad, test_reg$optdes), 3)

dcrit(inf_mat(grad, opt_plus_kl), 3)

eff("D-Optimality", inf_mat(grad, opt_plus_kl), inf_mat(grad, test_reg$optdes), k = 3)

d_opt <- optedr::opt_des("D-Optimality", model, parameters, par_values, design_space = c(1, 100), init_design = init_design,
                  join_thresh = join_thresh, delete_thresh = delete_thresh, delta = 1/2, weight_fun = weight_fun)

eff("D-Optimality", inf_mat(grad, opt_plus_kl), inf_mat(grad, d_opt$optdes), k = 3)

eff("D-Optimality", inf_mat(grad, test_reg$optdes), inf_mat(grad, d_opt$optdes), k = 3)

KL_antoine_water(opt_plus_kl) / KL_opt_val

KL_antoine_water(test_reg$optdes) / KL_opt_val

KL_antoine_water(d_opt$optdes) / KL_opt_val


# get_D_augment_val <- function(grad, init_design, alpha, x){
#   function(x)
#   return(do.call(dcrit, list((inf_mat(grad, sum_design_point(init_design, x, alpha)))),k))
# }


get_D_augment_data <- function(des_space, init_design, alpha, grad, k, grid.length = 10000){
  x_val <- seq(des_space[1], des_space[2], length.out = grid.length)
  y_val <- map_dbl(x_val, function(x) do.call(dcrit, list(inf_mat(grad, sum_design_point(init_design, x, alpha)), k)))
  return(data.frame("x" = x_val, "D_crit" = y_val))
}




alpha <- 0.25

Deff_plot_data <- get_D_augment_data(c(1, 100), kl_optdes, alpha, grad, 3)

Deff_plot_data$eff <- dcrit(inf_mat(grad, d_opt$optdes), 3) / Deff_plot_data$D_crit 

which.max(Deff_plot_data$eff)

Deff_plot_data[4129,]

new_data1 <- Deff_plot_data %>%
  filter(x <= 45.22589)

new_data2 <- Deff_plot_data %>%
  filter(92.12225 <= x)



ggplot() +
  geom_line(data = Deff_plot_data, mapping = aes(x = x, y = eff), color = "steelblue") +
  ggplot2::xlim(design_space[[1]], design_space[[2]]) +
  ggplot2::labs(x = "x", y = "D-Efficiency") +
  ggplot2::theme_bw() +
  ggplot2::geom_line(mapping = ggplot2::aes(x=new_data1$x,y=new_data1$eff), color = 'green3', linewidth = 0.7) +
  ggplot2::geom_line(mapping = ggplot2::aes(x=new_data2$x,y=new_data2$eff), color = 'green3', linewidth = 0.7) +
  theme(
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  ) 


antoine_regions$region[2]

naive_d_aug <- sum_design_point(kl_optdes, 41.86136, alpha)

KL_antoine_water(naive_d_aug) / KL_opt_val

eff("D-Optimality", inf_mat(grad, naive_d_aug), inf_mat(grad, d_opt$optdes), k = 3)

#########################################################
############ AUGMENT COMPOUND D-OPTIMALITY ##############
#########################################################

get_D_augment_data_compound <- function(des_space, init_design, alpha, lambda, grad1, grad2, k, grid.length = 10000){
  x_val <- seq(des_space[1], des_space[2], length.out = grid.length)
  y_val <- map_dbl(x_val, function(x) (1 - lambda) * do.call(dcrit, list(inf_mat(grad1, sum_design_point(init_design, x, alpha)), k)) 
                   + lambda * do.call(dcrit, list(inf_mat(grad2, sum_design_point(init_design, x, alpha)), k)))
  return(data.frame("x" = x_val, "D_crit" = y_val))
}


grad2 <- gradient(model, parameters, par_values)



alpha <- 0.25
lambda <- 0.54

Deff_comp_plot_data <- get_D_augment_data_compound(c(1, 100), kl_optdes, alpha, lambda, grad, grad2, 3)

comp_optdes <- data.frame("Point" = c(1.00000, 43.31377, 82.95844, 100.00000), "Weight" = c(0.1673876, 0.3053031, 0.2266398, 0.3006696))   

Deff_comp_plot_data$eff <- ((1 - lambda) * dcrit(inf_mat(grad, comp_optdes), 3) + lambda * dcrit(inf_mat(grad2, comp_optdes), 3)) / Deff_comp_plot_data$D_crit 

which.max(Deff_comp_plot_data$eff)

Deff_comp_plot_data[4131,]

new_data_comp1 <- Deff_comp_plot_data %>%
  filter(x <= 45.22589)

new_data_comp2 <- Deff_comp_plot_data %>%
  filter(92.12225 <= x)


left_join(Deff_comp_plot_data, Deff_plot_data, by = join_by(x == x))

ggplot() +
  geom_line(data = Deff_comp_plot_data, mapping = aes(x = x, y = eff), color = "steelblue") +
  ggplot2::xlim(design_space[[1]], design_space[[2]]) +
  ggplot2::labs(x = "x", y = "D-Efficiency") +
  ggplot2::theme_bw() +
  ggplot2::geom_line(mapping = ggplot2::aes(x=new_data_comp1$x,y=new_data_comp1$eff), color = 'green3', linewidth = 0.7) +
  ggplot2::geom_line(mapping = ggplot2::aes(x=new_data_comp2$x,y=new_data_comp2$eff), color = 'green3', linewidth = 0.7)


D_hom_antoine <- opt_des("D-Optimality", y ~ 10^(a-b/(c+x)), c("a","b","c"), c(8.07131,  1730.63, 233.426),  
                         c(1, 100))


D_het_antoine <- opt_des("D-Optimality", y ~ 10^(a-b/(c+x)), c("a","b","c"), c(8.07131,  1730.63, 233.426),  
                         c(1, 100), weight_fun = weight_antoine)



naive_d_aug_comp <- sum_design_point(kl_optdes, 41.89109, alpha)

design_efficiency(D_hom_antoine, naive_d_aug_comp)
design_efficiency(D_het_antoine, naive_d_aug_comp)

design_efficiency(D_hom_antoine, comp_optdes)
design_efficiency(D_het_antoine, comp_optdes)

KL_antoine_water(naive_d_aug_comp) / KL_opt_val





#################################################################
############ AUGMENT COMPOUND D-OPTIMALITY DIV EFF ##############
#################################################################

D_hom_antoine$crit_value 
D_het_antoine$crit_value


get_D_augment_data_compound2 <- function(des_space, init_design, alpha, lambda, grad1, grad2, k, grid.length = 10000){
  x_val <- seq(des_space[1], des_space[2], length.out = grid.length)
  y_val <- map_dbl(x_val, function(x) (1 - lambda) * D_het_antoine$crit_value / do.call(dcrit, list(inf_mat(grad, sum_design_point(init_design, x, alpha)), k)) 
                   + lambda * D_hom_antoine$crit_value / do.call(dcrit, list(inf_mat(grad2, sum_design_point(init_design, x, alpha)), k)))
  return(data.frame("x" = x_val, "D_crit" = y_val))
}


grad2 <- gradient(model, parameters, par_values)


alpha <- 0.25
lambda <- 0.50

Deff_comp_plot_data2 <- get_D_augment_data_compound2(c(1, 100), kl_optdes, alpha, lambda, grad, grad2, 3)

comp_optdes <- data.frame("Point" = c(1.00000, 43.31377, 82.95844, 100.00000), "Weight" = c(0.1673876, 0.3053031, 0.2266398, 0.3006696))   

Deff_comp_plot_data2$eff <- Deff_comp_plot_data2$D_crit / ((1 - lambda) * D_het_antoine$crit_value / dcrit(inf_mat(grad, comp_optdes), 3) + lambda * D_hom_antoine$crit_value / dcrit(inf_mat(grad2, comp_optdes), 3)) 

which.max(Deff_comp_plot_data2$eff)

Deff_comp_plot_data2[5921,]


new_data_comp1 <- Deff_comp_plot_data2 %>%
  filter(x <= 45.22589)

new_data_comp2 <- Deff_comp_plot_data2 %>%
  filter(92.12225 <= x)


left_join(Deff_comp_plot_data2, Deff_plot_data, by = join_by(x == x))

ggplot() +
  geom_line(data = Deff_comp_plot_data2, mapping = aes(x = x, y = eff), color = "steelblue") +
  ggplot2::xlim(design_space[[1]], design_space[[2]]) +
  ggplot2::labs(x = "x", y = "D-Efficiency") +
  ggplot2::theme_bw() +
  ggplot2::geom_line(mapping = ggplot2::aes(x=new_data_comp1$x,y=new_data_comp1$eff), color = 'green3', linewidth = 0.7) +
  ggplot2::geom_line(mapping = ggplot2::aes(x=new_data_comp2$x,y=new_data_comp2$eff), color = 'green3', linewidth = 0.7)


D_hom_antoine <- opt_des("D-Optimality", y ~ 10^(a-b/(c+x)), c("a","b","c"), c(8.07131,  1730.63, 233.426),  
                         c(1, 100))


D_het_antoine <- opt_des("D-Optimality", y ~ 10^(a-b/(c+x)), c("a","b","c"), c(8.07131,  1730.63, 233.426),  
                         c(1, 100), weight_fun = weight_antoine)



naive_d_aug_comp2 <- sum_design_point(kl_optdes, 59.61386, alpha)

design_efficiency(D_hom_antoine, naive_d_aug_comp)
design_efficiency(D_het_antoine, naive_d_aug_comp)

design_efficiency(D_hom_antoine, naive_d_aug_comp2)
design_efficiency(D_het_antoine, naive_d_aug_comp2)

KL_antoine_water(naive_d_aug_comp2) / KL_opt_val
attributes(D_hom_antoine$convergence)



D_hom_antoine$crit_value

#########################################################
############ INVERSE: AUGMENT D-OPTIMAL DESIGN ##########
#########################################################
weight_fun <- function(x) (10^(a - b/(c + x)))^(-1)
resAnt.D <-opt_des("D-Optimality", y ~ 10^(a-b/(c+x)),
                   c("a","b","c"), c(8.07131,  1730.63, 233.426),
                   c(1, 100), weight_fun = weight_fun)

eff_plot_data_D <- get_KL_augment_data(KL_antoine_water, design_space, resAnt.D$optdes, alpha)

eff_plot_data_D$eff <- eff_plot_data_D$KL_value / KL_opt_val

init_eff <- KL_antoine_water(resAnt.D$optdes) / KL_opt_val

kl_daug1 <- sum_design_point(resAnt.D$optdes, 1, alpha)

KL_antoine_water(kl_daug1) / KL_opt_val

design_efficiency(resAnt.D, kl_daug1)

plot25 <- ggplot() +
  geom_line(data = eff_plot_data_D, mapping = aes(x = x, y = eff), color = "steelblue") +
  ggplot2::xlim(design_space[[1]], design_space[[2]]) +
  ggplot2::labs(x = "x", y = "Efficiency") +
  ggplot2::theme_bw() +
  ylab(TeX("eff$_{KL}(\\xi_{k+1}^{(x)})$")) +
  theme(
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  ) # +
  geom_hline(yintercept = init_eff, color = "goldenrod3")

test_Desgin <- data.frame( Point= c(1.00000,  41.86306, 100.00000), Weight = c(0.3335231,  0.3331976,  0.3332792))



eff_plot_data_D2 <- get_KL_augment_data(KL_antoine_water, design_space, resAnt.D$optdes, 0.4)

eff_plot_data_D2$eff <- eff_plot_data_D2$KL_value / KL_opt_val

init_eff <- KL_antoine_water(resAnt.D$optdes) / KL_opt_val

kl_daug2 <- sum_design_point(resAnt.D$optdes, 1.633663, 0.4)

KL_antoine_water(kl_daug2) / KL_opt_val

design_efficiency(resAnt.D, kl_daug2)


plot40 <- ggplot() +
  geom_line(data = eff_plot_data_D2, mapping = aes(x = x, y = eff), color = "steelblue") +
  ggplot2::xlim(design_space[[1]], design_space[[2]]) +
  ggplot2::labs(x = "x", y = "Efficiency") +
  ggplot2::theme_bw() +
  ylab(TeX("eff$_{KL}(\\xi_{k+1}^{(x)})$")) +
  theme(
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )# +
  geom_hline(yintercept = init_eff, color = "goldenrod3")

test_Desgin <- data.frame( Point= c(1.00000,  41.86306, 100.00000), Weight = c(0.3335231,  0.3331976,  0.3332792))

eff_plot_data_D2[1,]

library(gridExtra)

grid.arrange(plot25, plot40, ncol=2)



#########################################################
############ KL-D Compound design #######################
#########################################################
design_space <- c(1, 100)
comp_design_KLD <- CompoundDKL(0.5, init_design, grad, design_space[[1]], design_space[[2]], 1000, join_thresh, delete_thresh, k, tol, tol2)

comp_design_KLD$sens


CompoundDKL <- function(alpha, init_design, grad, min, max, grid.length, join_thresh, delete_thresh, k, tol, tol2) {
  Point <- NULL
  crit_val <- numeric(2122)
  index <- 1
  # Maximum iterations for the optimize weights loop
  maxiter <- 300
  cli::cli_progress_bar("Calculating optimal design")
  for (i in 1:300) {
    cli::cli_progress_update()
    M <- inf_mat(grad, init_design)
    crit_val[index] <- alpha * dcrit(M, k) + (1 - alpha) * KL_antoine_water(init_design)
    index <- index + 1
    sensM <- function(xval) {
      dsensval <- dsens(grad, M)
      klsensval <- KL_antoine_water(init_design)
      alpha * (dsensval(xval) - k) - (1 - alpha) * klsensval
    }
    xmax <- findmax_extremes(sensM, min, max, grid.length)
    # if ((sensM(xmax)) / (alpha * k - (1 - alpha) * KL_antoine_water(init_design)) < tol2) {
    #   message("\n", crayon::blue(cli::symbol$info), " Stop condition reached: difference between sensitivity and criterion < ", tol2)
    #   break
    # }
    init_design <- update_design(init_design, xmax, join_thresh, 1/(index + 2))
    # iter <- 1
    # stopw <- FALSE
    # while (!stopw) {
    #   weightsInit <- init_design$Weight
    #   M <- inf_mat(grad, init_design)
    #   crit_val[index] <- alpha * dcrit(M, k) + (1 - alpha) * icrit(M, matB)
    #   index <- index + 1
    #   sensM <- dsens(grad, M)
    #   init_design$Weight <- update_weights_compound(init_design, sensM, k, icrit(M, matB), alpha, delta_weights)
    #   stopw <- any(max(abs(weightsInit - init_design$Weight) < tol)) || iter >= maxiter
    #   iter <- iter + 1
    # }
    init_design <- delete_points(init_design, delete_thresh)
    if (i %% 5 == 0) {
      init_design <- update_design_total(init_design, join_thresh)
    }
    if (i == 300) {
      message("\n", crayon::blue(cli::symbol$info), " Stop condition not reached, max iterations performed")
    }
  }
  cli::cli_progress_update(force = TRUE)
  base::cat("")
  crit_val[index] <- alpha * dcrit(M, k) + (1 - alpha) * KL_antoine_water(init_design)
  crit_val <- crit_val[1:(length(crit_val) - sum(crit_val == 0))]
  conv <- data.frame("criteria" = crit_val, "step" = seq(1, length(crit_val), 1))
  conv_plot <- plot_convergence(conv)
  
  init_design <- dplyr::arrange(init_design, Point)
  rownames(init_design) <- NULL
  
  M <- inf_mat(grad, init_design)
  sensM <- function(xval) {
    dsensval <- dsens(grad, M)
    klsensval <- KL_antoine_water(init_design)
    alpha * (dsensval(xval) - k) - (1 - alpha) * klsensval
  }
  xmax <- findmax_extremes(sensM, min, max, grid.length * 10)
  dsenstest <- dsens(grad, M)
  
  atwood <- (alpha * k - (1 - alpha) * KL_antoine_water(init_design)) / dsenstest(xmax) * 100
  
  message(crayon::blue(cli::symbol$info), " The lower bound for efficiency is ", atwood, "%")
  
  plot_opt <- plot_sens_orig(min, max, sensM, 0)
  l_return <- list(
    "optdes" = init_design, "convergence" = conv_plot,
    "sens" = plot_opt, "criterion" = "D-Optimality", "crit_value" = crit_val[length(crit_val)]
  )
  attr(l_return, "hidden_value") <- k
  attr(l_return, "gradient") <- grad
  attr(l_return, "atwood") <- atwood
  class(l_return) <- "optdes"
  l_return
}

