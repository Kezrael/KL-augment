# Auxiliar function for algorithms --------------------------------------



#' Weight function per distribution
#'
#' @param model formula describing the model to use. Must use x as the variable.
#' @param char_vars character vector with the parameters of the models, as written in the \code{formula}
#' @param values numeric vector with the parameters nominal values, in the same order as given in \code{parameters}.
#' @param distribution character variable specifying the probability distribution of the response. Can be one of the following:
#'   * 'Homoscedasticity'
#'   * 'Gamma', which can be used for exponential or normal heteroscedastic with constant relative error
#'   * 'Poisson'
#'   * 'Logistic'
#'   * 'Log-Normal' (work in progress)
#'
#' @return one variable function that represents the square of the structure of variance, in case of heteroscedastic variance of the response.
#'
weight_function <- function(model, char_vars, values, distribution = "Homoscedasticity") {
  # vars <- as.list(match.call())[-(1:2)]
  # char_vars <- sapply(vars, as.character)
  if(!(distribution %in% c("Poisson", "Gamma", "Logistic", #"log-normal",
                           "Homoscedasticity"))){
    warning(crayon::yellow(cli::symbol$warning), " Not a valid distribution specified, using a normal homoscedastic")
    return(function(x) 1)
  }
  else if(distribution == "Homoscedasticity"){
    return(function(x) 1)
  }
  else if(distribution == "Poisson"){
    cmd <- utils::tail(as.character(model),1)
    expres <- parse(text=cmd)
    lista <- values
    names(lista) <- char_vars
    f <- function(x_val) {
      exp(eval(expres, c(lista, list("x" = x_val)))/2)
    }
  }
  else if(distribution == "Gamma"){
    cmd <- utils::tail(as.character(model),1)
    expres <- parse(text=cmd)
    lista <- values
    names(lista) <- char_vars
    f <- function(x_val) {
      (eval(expres, c(lista, list("x" = x_val))))^(-1)
    }
  }
  else if(distribution == "Logistic"){
    cmd <- utils::tail(as.character(model),1)
    expres <- parse(text=cmd)
    lista <- values
    names(lista) <- char_vars
    f <- function(x_val) {
      sqrt(exp(-eval(expres, c(lista, list("x" = x_val))))/(1 + exp(-eval(expres, c(lista, list("x" = x_val)))))^2)
    }
  }
  else if(distribution == "Log-normal"){
    return(function(x) 1)
  }
  return(f)
}


#' Find Minimum Value
#'
#' @description
#' Searches the maximum of a function over a grid on a given grid.
#'
#' @param sens a single variable numeric function to evaluate.
#' @param min minimum value of the search grid.
#' @param max maximum value of the search grid.
#' @param grid.length length of the search grid.
#'
#' @return The value of the minimum
findminval <- function(sens, min, max, grid.length) {
  if(min <= max){
    grid <- seq(min, max, length.out = grid.length)
  }
  else {
    grid <- seq(max, min, length.out = grid.length)
  }
  minval <- min(purrr::map_dbl(grid, sens))
  return(minval)
}

#' Find Maximum Value
#'
#' @description
#' Searches the maximum of a function over a grid on a given interval.
#'
#' @param sens A single variable numeric function to evaluate.
#' @param min Minimum value of the search interval.
#' @param max Maximum value of the search interval.
#' @param grid.length Length of the search interval.
#'
#' @return The value of the maximum
findmaxval <- function(sens, min, max, grid.length) {
  if(min <= max){
    grid <- seq(min, max, length.out = grid.length)
  }
  else {
    grid <- seq(max, min, length.out = grid.length)
  }
  maxval <- max(purrr::map_dbl(grid, sens))
  return(maxval)
}



#' Deletes duplicates points
#'
#' @description
#' Within a vector of points, deletes points that are close enough (less than the tol parameter).
#' Returns the points without the "duplicates"
#'
#' @param points Points to be updated
#' @param tol Tolerance for which two points are considered the same
#'
#' @return The points without duplicates
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


#' Find where the candidate points region starts
#'
#' @description
#' Given the crosspoints and the sensitivity function, this function
#' finds where the candidate points region starts, either on the extreme of
#' the space of the design or the first crosspoints
#'
#' @param cross Vector of crosspoints in the sensitivity function given an efficiency and weight
#' @param min Minimum of the space of the design
#' @param max Maximum of the space of the design
#' @param val Value of the sensitivity function at the crosspoints
#' @param sens_opt Sensitivity function
#'
#' @return True if the candidate points region starts on the minimum, False otherwise
#'
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


#' Parity of the crosspoints
#'
#' @description
#' Determines if the number of crosspoints is even or odd given the vector
#' of crosspoints
#'
#' @param cross Vector of crosspoints in the sensitivity function given an efficiency and weight
#'
#' @return True if the number of crosspoints is even, false otherwise
getPar <- function(cross){
  return(length(cross)%%2 == 0)
}

#' Give effective limits to candidate points region
#'
#' @description
#' Given the start of the candidates points region, the parity of the crosspoints
#' and the boundaries of the space of the design returns the effective limits of
#' the candidate points region. Those points, taken in pairs from the first to
#' the last delimit the region.
#'
#' @param cross Vector of crosspoints in the sensitivity function given an efficiency and weight
#' @param min Minimum of the space of the design
#' @param max Maximum of the space of the design
#' @param start Boolean that gives the effective start of the candidate points region
#' @param par Boolean with the parity of the region
#'
#' @return Vector of effective limits of the candidate points region. Taken
#' in pairs from the beginning delimit the region.
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




#' Update weight D-Optimality
#'
#' @description
#' Implementation of the weight update formula for D-Optimality used to optimize the weights of a design,
#' which is to be applied iteratively until no sizable changes happen.
#'
#' @param design Design to optimize the weights from. It's a dataframe with two columns:
#'   * \code{Point} contains the support points of the design.
#'   * \code{Weight} contains the corresponding weights of the \code{Point}s.
#' @param sens Sensibility function for the design and model.
#' @param k Number of parameters of the model.
#' @param delta A parameter of the algorithm that can be tuned. Must be \eqn{0< delta < 1}.
#'
#' @return returns the new weights of the design after one iteration.
update_weights <- function(design, sens, k, delta) {
  weights <- design$Weight * (purrr::map_dbl(design$Point, sens) / k)^delta
  weights[is.nan(weights)] <- 0
  return(weights / sum(weights))
}


#' Update weight Ds-Optimality
#'
#' @description
#' Implementation of the weight update formula for Ds-Optimality used to optimize the weights of a design,
#' which is to be applied iteratively until no sizable changes happen.
#'
#'
#' @param s number of parameters of interest of the model
#' @inheritParams update_weights
#'
#' @return returns the new weights of the design after one iteration.
update_weightsDS <- function(design, sens, s, delta) {
  weights <- design$Weight * (purrr::map_dbl(design$Point, sens) / s)^delta
  weights[is.nan(weights)] <- 0
  return(weights)
}

#' Update weight I-Optimality
#'
#' @description
#' Implementation of the weight update formula for I-Optimality used to optimize the weights of a design,
#' which is to be applied iteratively until no sizable changes happen. A-Optimality if instead of the
#' integral matrix the identity function is used.
#'
#'
#' @param crit Value of the criterion function for I-Optimality.
#' @inheritParams update_weights
#'
#' @return returns the new weights of the design after one iteration.
update_weightsI <- function(design, sens, crit, delta) {
  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  weights <- design$Weight * exponent((purrr::map_dbl(design$Point, sens) / crit), delta)
  weights[is.nan(weights)] <- 0
  # weights <- weights/sum(weights)
  return(weights/sum(weights))
}


#' Update Design with new point
#'
#' @description
#' Updates a design adding a new point to it. If the added point is closer than \code{delta} to an existing
#' point of the design, the two points are merged together as their arithmetic average. Then updates the weights
#' to be equal to all points of the design.
#'
#' @param design Design to update. It's a dataframe with two columns:
#'   * \code{Point} contains the support points of the design.
#'   * \code{Weight} contains the corresponding weights of the \code{Point}s.
#' @param xmax The point to add as a numeric value.
#' @param delta Threshold which defines how close the new point has to be to any of the existing ones in order to
#'   merge with them.
#' @param new_weight Number with the weight for the new point.
#'
#' @return The updated design.
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


#' Merge close points of a design
#'
#' @description
#' Takes a design and merge together all points that are closer between them than a certain threshold \code{delta}.
#'
#' @param design The design to update. It's a dataframe with two columns:
#'   * \code{Point} contains the support points of the design.
#'   * \code{Weight} contains the corresponding weights of the \code{Point}s.
#' @param delta Threshold which defines how close two points have to be to any of the existing ones in order to
#'   merge with them.
#'
#' @return The updated design.
update_design_total <- function(design, delta) {
  updated <- FALSE
  finished <- FALSE
  while (!finished) {
    for (i in 1:(length(design$Point) - 1)) {
      absdiff <- abs(design$Point[-seq(1, i)] - design$Point[i]) < delta
      if (any(absdiff)) {
        updated <- TRUE
        design$Weight[which(absdiff)[1] + i] <- design$Weight[which(absdiff)[1] + i] + design$Weight[i]
        design <- design[-i, ]
        break
      }
    }
    finished <- !updated
    updated <- FALSE
  }
  return(design)
}


#' Remove low weight points
#'
#' @description
#' Removes the points of a design with a weight lower than a threshold, \code{delta}, and distributes that weights
#' proportionally to the rest of the points.
#'
#' @param design The design from which to remove points as a dataframe with two columns:
#'   * \code{Point} contains the support points of the design.
#'   * \code{Weight} contains the corresponding weights of the \code{Point}s.
#' @param delta The threshold from which the points with such a weight or lower will be removed.
#'
#' @return The design without the removed points.
delete_points <- function(design, delta) {
  updatedDesign <- design[design$Weight > delta, ]
  updatedDesign$Weight <- updatedDesign$Weight / sum(updatedDesign$Weight)
  return(updatedDesign)
}

# General auxiliar functions --------------------------
# Cálculo de la traza de una matriz
#' Trace
#'
#' @description
#' Return the mathematical trace of a matrix, the sum of its diagonal elements.
#'
#'
#' @param M The matrix from which to calculate the trace.
#'
#' @return The trace of the matrix.
tr <- function(M) {
  return(sum(diag(M)))
}



plot_sens <- function(min, max, x.grid, sens_function, criterion_value) {
  x <- y <- NULL
  sens_grid <- purrr::map_dbl(x.grid, sens_function)
  
  data <- data.frame(x = x.grid, y = sens_grid, grp = rep(1:(length(x_grid) %/% 1000), each = 1000))
  
  sensibility <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y, group = grp), color = "steelblue3") +
    ggplot2::stat_function(fun = function(x) criterion_value, col = "goldenrod3") +
    ggplot2::xlim(min(x.grid), max(x.grid)) +
    ggplot2::labs(x = "X", y = "Y")
  return(sensibility)
}

plot_sens_orig <- function(min, max, sens_function, criterion_value) {
  x <- y <- NULL
  grid <- seq(min, max, length.out = 10000)
  sens_grid <- purrr::map_dbl(grid, sens_function)
  
  sensibility <- ggplot2::ggplot(data = data.frame(x = grid, y = sens_grid), mapping = ggplot2::aes(x = x)) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y), color = "steelblue3") +
    ggplot2::stat_function(fun = function(x) criterion_value, col = "goldenrod3") +
    ggplot2::xlim(min, max) +
    ggplot2::labs(x = "X", y = "Y")
}

#' Plot Convergence of the algorithm
#'
#' @description
#' Plots the criterion value on each of the steps of the algorithm, both for optimizing weights and points,
#' against the total step number.
#'
#' @param convergence A dataframe with two columns:
#'   * \code{criteria} contains value of the criterion on each step.
#'   * \code{step} contains number of the step.
#'
#' @return A ggplot object with the \code{criteria} in the \code{y} axis and \code{step} in the \code{x} axis.
plot_convergence <- function(convergence) {
  step <- criteria <- NULL
  ggplot2::ggplot(data = convergence, ggplot2::aes(x = step, y = criteria)) +
    ggplot2::geom_line(color = "coral1") +
    ggplot2::theme_bw()
}


#' Integrate IM
#'
#' @description
#' Integrates the information matrix over the region of interest to calculate matrix B to be used in I-Optimality
#' calculation.
#'
#' @param grad function of partial derivatives of the model.
#' @param k number of unknown parameters of the model.
#' @param reg_int optional numeric vector of two components with the bounds of the interest region for I-Optimality.
#'
#' @return The integrated information matrix.
integrate_reg_int <- function(grad, k, reg_int) {
  matrix_int <- 0 * diag(k)
  for (i in 1:k) {
    for (j in 1:k) {
      if (j >= i) {
        int_part <- function(x_value) {
          purrr::map_dbl(x_value, function(x_value) grad(x_value)[i] * grad(x_value)[j] / (reg_int[2] - reg_int[1]))
        }
        matrix_int[i, j] <- stats::integrate(int_part, lower = reg_int[1], upper = reg_int[2])$value
      }
      else {
        matrix_int[i, j] <- matrix_int[j, i]
      }
    }
  }
  return(matrix_int)
}


#' Summary function for optdes
#'
#' @param object An object of class \code{optdes}.
#' @param ... Possible extra arguments for the summary
#'
#' @export
#'
#' @examples
#' rri <- opt_des(Criterion = "I-Optimality", model = y ~ a * exp(-b / x),
#'   parameters = c("a", "b"), par_values = c(1, 1500), design_space = c(212, 422),
#'   reg_int = c(380, 422))
#' summary(rri)
summary.optdes <- function(object, ...) {
  cat("Model: \n")
  print(attr(object, "model"))
  cat("and weight function: \n")
  print(attr(attr(object, "weight_fun"), "srcref"))
  cat("Optimal design for ", object$criterion, ":\n")
  print.data.frame(object$optdes, ...)
  cat("\n Minimum efficiency (Atwood): ", paste0(attr(object, "atwood"), "%"))
  cat("\n Criterion value: ", object$crit_value)
}


#' Print function for optdes
#'
#' @param x An object of class \code{optdes}.
#' @param ... Possible extra arguments for printing dataframes
#'
#' @export
#'
#' @examples
#' rri <- opt_des(Criterion = "I-Optimality", model = y ~ a * exp(-b / x),
#'   parameters = c("a", "b"), par_values = c(1, 1500), design_space = c(212, 422),
#'   reg_int = c(380, 422))
#' print(rri)
print.optdes <- function(x, ...) {
  print.data.frame(x$optdes, ...)
}



#' Plot function for optdes
#'
#' @param x An object of class \code{optdes}.
#' @param ... Possible extra arguments for plotting dataframes
#'
#' @export
#'
#' @examples
#' rri <- opt_des(Criterion = "I-Optimality", model = y ~ a * exp(-b / x),
#'   parameters = c("a", "b"), par_values = c(1, 1500), design_space = c(212, 422),
#'   reg_int = c(380, 422))
#' plot(rri)
plot.optdes <- function(x, ...) {
  Point <- Value <- Weight <- NULL
  x$optdes[["Value"]] <- rep(0, nrow(x$optdes))
  x$optdes[["Weight"]] <- round(x$optdes[["Weight"]], 2)
  p <- x$sens + ggplot2::geom_point(data = x$optdes, ggplot2::aes(x = Point, y = Value)
                                    , size = 4, shape = 16, color = "darkgreen") +
    ggplot2::geom_text(data = x$optdes, ggplot2::aes(x = Point, y = Value, label = Weight),
                       hjust=1.5, vjust=1.5) +
    ggplot2::labs(x = "Design Space", y = "Sensitivity Function")
  p
}

# Gradient Vector ------------------------------------


#' Gradient function
#'
#' @description
#' Calculates the gradient function of a \code{model} with respect to the parameters, \code{char_vars}, evaluates
#' it at the provided \code{values} and returns the result as a function of the variable \code{x}.
#'
#' @param model formula describing the model, which must contain only \code{x}, the parameters defined in
#'   \code{char_vars} and the numerical operators.
#' @param char_vars character vector of the parameters of the model.
#' @param values numeric vector with the nominal values of the parameters in \code{char_vars}.
#' @param weight_fun optional function variable that represents the square of the structure of variance, in case of heteroscedastic variance of the response
#'
#' @return A function depending on \code{x} that's the gradient of the \code{model} with respect to \code{char_vars}
gradient <- function(model, char_vars, values, weight_fun = function(x) 1) {
  # vars <- as.list(match.call())[-(1:2)]
  # char_vars <- sapply(vars, as.character)
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


# Information Matrix ------------------------------------


#' Information Matrix
#'
#' @description
#' Given the gradient vector of a model in a single variable model and a design, calculates the information matrix.
#'
#' @param grad A function in a single variable that returns the partial derivatives vector of the model.
#' @param design A dataframe that represents the design. Must have two columns:
#'   * \code{Point} contains the support points of the design.
#'   * \code{Weight} contains the corresponding weights of the \code{Point}s.
#'
#' @return The information matrix of the design, a \eqn{k\times k} matrix where k is the length of the gradient.
inf_mat <- function(grad, design) {
  matrix_ret <- 0 * diag(length(grad(design$Point[[1]])))
  for (i in seq_along(design$Weight)) {
    f_col <- as.matrix(grad(design$Point[[i]]), nrow = 1, byrow = TRUE, dimnames = NULL)
    matrix_ret <- matrix_ret + (t(f_col) %*% f_col) * design$Weight[[i]]
  }
  return(matrix_ret)
}

# Sensibility Function ------------------------------------

#' Master function to calculate the sensitivity function
#'
#' @description
#' Calculates the sensitivity function given the desired \code{Criterion}, an information matrix and other
#' necessary values depending on the chosen criterion.
#'
#' @param Criterion character variable with the chosen optimality criterion. Can be one of the following:
#'   * 'D-Optimality'
#'   * 'Ds-Optimality'
#'   * 'A-Optimality'
#'   * 'I-Optimality'
#' @param grad A function in a single variable that returns the partial derivatives vector of the model.
#' @param M Information Matrix for the sensitivity function.
#' @param par_int Numeric vector of the indexes of the parameters of interest for Ds-Optimality.
#' @param matB Matrix resulting from the integration of the one-point Information Matrix along the interest
#'   region.
#'
#' @return The sensitivity function as a matrix of single variable.
sens <- function(Criterion, grad, M, par_int = c(1), matB = NA) {
  if (identical(Criterion, "D-Optimality")) {
    return(dsens(grad, M))
  }
  else if (identical(Criterion, "Ds-Optimality")) {
    return(dssens(grad, M, par_int))
  }
  else if (identical(Criterion, "A-Optimality")) {
    return(isens(grad, M, diag(nrow(M))))
  }
  else if (identical(Criterion, "I-Optimality")) {
    return(isens(grad, M, matB))
  }
}


#' Sensitivity function for D-Optimality
#'
#' @description
#' Calculates the sensitivity function from the gradient vector and the Identity Matrix.
#'
#' @inherit sens params return
#'
dsens <- function(grad, M) {
  invMat <- solve(M)
  sens_ret <- function(xval) {
    f_col <- as.matrix(grad(xval), nrow = 1, ncol = 3, byrow = TRUE, dimnames = NULL)
    return(f_col %*% invMat %*% t(f_col))
  }
}


#' Sensitivity function for Ds-Optimality
#'
#' @description
#' Calculates the sensitivity function from the gradient vector, the Identity Matrix and the parameters of
#' interest.
#'
#' @inherit sens params return
#'
dssens <- function(grad, M, par_int) {
  invMat <- solve(M)
  if (length(M[-par_int, -par_int]) == 1) {
    invMat22 <- 1 / M[-par_int, -par_int]
  } else {
    invMat22 <- solve(M[-par_int, -par_int])
  }
  sens_ret <- function(xval) {
    f_col <- as.matrix(grad(xval), nrow = 1, byrow = TRUE, dimnames = NULL)
    return(f_col %*% invMat %*% t(f_col) - f_col[-par_int] %*% invMat22 %*% as.matrix(f_col[-par_int], ncol = 1))
  }
}


#' Sensitivity function for I-Optimality
#'
#' @description
#' Calculates the sensitivity function from the gradient vector, the Information Matrix and the integral of the
#' one-point Identity Matrix over the interest region. If instead the identity matrix is used, it can be used
#' for A-Optimality.
#'
#' @inherit sens params return
#'
isens <- function(grad, M, matB) {
  invMat <- solve(M)
  sens_ret <- function(xval) {
    f_col <- as.matrix(grad(xval), nrow = 1, ncol = 3, byrow = TRUE, dimnames = NULL)
    return(f_col %*% invMat %*% matB %*% invMat %*% t(f_col))
  }
}

# Optimality Criteria -----------------------

#' Master function for the criterion function
#'
#' @description
#' Depending on the Criterion input, the function returns the output of the corresponding criterion function given
#' the information matrix.
#'
#' @param Criterion character variable with the chosen optimality criterion. Can be one of the following:
#'   * 'D-Optimality'
#'   * 'Ds-Optimality'
#'   * 'A-Optimality'
#'   * 'I-Optimality'
#' @param M information matrix for which the criterion value wants to be calculated.
#' @param k numeric variable with the number of parameters of the model. Taken from the number of rows of the matrix if omitted.
#' @param par_int numeric vector with the index of the parameters of interest of the model. Only for "Ds-Optimality".
#' @param matB matrix of the integral of the information matrix over the interest region. Only for "I-Optimality".
#'
#' @return Numeric value of the optimality criterion for the information matrix.
crit <- function(Criterion, M, k = 0, par_int = c(1), matB = NA) {
  if (identical(Criterion, "D-Optimality")) {
    return(dcrit(M, k))
  }
  else if (identical(Criterion, "Ds-Optimality")) {
    return(dscrit(M, par_int))
  }
  else if (identical(Criterion, "A-Optimality")) {
    return(icrit(M, diag(k)))
  }
  else if (identical(Criterion, "I-Optimality")) {
    return(icrit(M, matB))
  }
  else if (identical(Criterion, "L-Optimality")) {
    return(icrit(M, matB))
  }
}


#' Criterion function for D-Optimality
#'
#' @description
#' Calculates the value of the D-Optimality criterion function, which follows the expression:
#' \deqn{\phi_D = \left(\frac{1}{|M|}\right)^{1/k}}
#'
#'
#' @param M information matrix for which the criterion value wants to be calculated.
#' @param k numeric variable with the number of parameters of the model. Taken from the number of rows of the matrix if omitted.
#'
#'
#' @return numeric value of the D-optimality criterion for the information matrix.
dcrit <- function(M, k) {
  if (k == 0) k <- nrow(M)
  return((1 / det(M))^(1 / k))
}


#' Criterion function for Ds-Optimality
#'
#' @description
#' Calculates the value of the Ds-Optimality criterion function, which follows the expression:
#' \deqn{\phi_{Ds} = \left(\frac{|M_{22}|}{|M|}\right)^{1/s}}
#'
#'
#' @param M information matrix for which the criterion value wants to be calculated.
#' @param par_int numeric vector with the index of the parameters of interest of the model.
#'
#'
#' @return Numeric value of the Ds-optimality criterion for the information matrix.
dscrit <- function(M, par_int) {
  if (length(M[-par_int, -par_int]) == 1) {
    return((M[-par_int, -par_int] / det(M))^(1 / length(par_int)))
  }
  else {
    return((det(M[-par_int, -par_int]) / det(M))^(1 / length(par_int)))
  }
}

#' Criterion function for I-Optimality
#'
#' @description
#' Calculates the value of the I-Optimality criterion function, which follows the expression:
#' \deqn{\phi_I = Tr(M^{-1} \cdot B)}
#'
#'
#' @param M information matrix for which the criterion value wants to be calculated.
#' @param matB matrix of the integral of the information matrix over the interest region. Identity matrix for
#'   A-Optimality.
#'
#'
#' @return Numeric value of the I-optimality criterion for the information matrix.
icrit <- function(M, matB) {
  return(tr(matB %*% solve(M)))
}



# Efficiency-------- -----------------------


#' Efficiency between two Information Matrices
#'
#' @param Criterion character variable with the chosen optimality criterion. Can be one of the following:
#'   * 'D-Optimality'
#'   * 'Ds-Optimality'
#'   * 'A-Optimality'
#'   * 'I-Optimality'
#' @param mat1 first information matrix, for the numerator.
#' @param mat2 second information matrix, for the denominator.
#' @param k number of parameters of the model. Taken from the number of rows of the matrix if omitted.
#' @param intPars numeric vector with the index of the parameters of interest of the model. Only for "Ds-Optimality".
#' @param matB matrix of the integral of the information matrix over the interest region. Only for "I-Optimality".
#'
#' @return Efficiency of first design with respect to the second design, as a decimal number.
eff <- function(Criterion, mat1, mat2, k = 0, intPars = c(1), matB = NA) {
  if (identical(Criterion, "D-Optimality")) {
    if (k == 0) k <- nrow(mat1)
    return((det(mat1) / det(mat2))^(1 / k))
  }
  else if (identical(Criterion, "Ds-Optimality")) {
    if (length(mat1[-intPars, -intPars]) == 1) {
      return((mat2[-intPars, -intPars] / det(mat2) / (mat1[-intPars, -intPars] / det(mat1)))^(1 / length(intPars)))
    }
    else {
      return((det(mat2[-intPars, -intPars]) / det(mat2) / (det(mat1[-intPars, -intPars]) / det(mat1)))^(1 / length(intPars)))
    }
  }
  else if (identical(Criterion, "A-Optimality")) {
    if (k == 0) k <- nrow(mat1)
    return(tr(diag(k) %*% solve(mat2)) / tr(diag(k) %*% solve(mat1)))
  }
  else if (identical(Criterion, "I-Optimality")) {
    return(tr(matB %*% solve(mat2)) / tr(matB %*% solve(mat1)))
  }
}


#' Efficiency between optimal design and a user given design
#'
#' @description
#' Takes an optimal design provided from the function \code{opt_des} and a user given design and compares their
#' efficiency
#'
#' @seealso opt_des
#'
#' @param opt_des_obj an object given by the function \code{opt_des}.
#' @param design dataframe that represents the design. Must have two columns:
#'   * \code{Point} contains the support points of the design.
#'   * \code{Weight} contains the corresponding weights of the \code{Point}s.
#'
#' @return The efficiency as a value between 0 and 1
#' @export
#'
#' @examples
#' result <- opt_des("D-Optimality", y ~ a * exp(-b / x), c("a", "b"), c(1, 1500), c(212, 422))
#' design <- data.frame("Point" = c(220, 240, 400), "Weight" = c(1 / 3, 1 / 3, 1 / 3))
#' design_efficiency(result, design)
design_efficiency <- function(opt_des_obj, design) {
  # check_efficiency_input(opt_des_obj, design) COMPROBAR QUE EL NUMERO DE POUNTOS ES >= LENGTH GRAD/NROW MAT
  mat1 <- inf_mat(attr(opt_des_obj, "gradient"), design)
  mat2 <- inf_mat(attr(opt_des_obj, "gradient"), opt_des_obj$optdes)
  if (identical(opt_des_obj$criterion, "D-Optimality")) {
    eff <- (det(mat1) / det(mat2))^(1 / attr(opt_des_obj, "hidden_value"))
    message(crayon::blue(cli::symbol$info), " The efficiency of the design is ", eff * 100, "%")
    return(eff)
  }
  else if (identical(opt_des_obj$criterion, "Ds-Optimality")) {
    int_pars <- attr(opt_des_obj, "hidden_value")
    if (length(mat1[-int_pars, -int_pars]) == 1) {
      eff <- (mat2[-int_pars, -int_pars] / det(mat2) / (mat1[-int_pars, -int_pars] / det(mat1)))^(1 / length(int_pars))
      message(crayon::blue(cli::symbol$info), " The efficiency of the design is ", eff * 100, "%")
      return(eff)
    }
    else {
      eff <- (det(mat2[-int_pars, -int_pars]) / det(mat2) / (det(mat1[-int_pars, -int_pars]) / det(mat1)))^(1 / length(int_pars))
      message(crayon::blue(cli::symbol$info), " The efficiency of the design is ", eff * 100, "%")
      return(eff)
    }
  }
  else if (identical(opt_des_obj$criterion, "I-Optimality")) {
    eff <- tr(attr(opt_des_obj, "hidden_value") %*% solve(mat2)) / tr(attr(opt_des_obj, "hidden_value") %*% solve(mat1))
    message(crayon::blue(cli::symbol$info), " The efficiency of the design is ", eff * 100, "%")
    return(eff)
  }
}


findmax <- function(sensM, x.grid){
  xmax <- x.grid[which.max(purrr::map_dbl(x.grid, sensM))]
  return(xmax)
}

findmax_extremes <- function(sens, min, max, grid.length) {
  if (min <= max) {
    grid <- seq(min, max, length.out = grid.length)
  }
  else {
    grid <- seq(max, min, length.out = grid.length)
  }
  xmax <- grid[which.max(purrr::map(grid, sens))]
  return(xmax)
}

split_grid <- function(region, grid.length = 1000){
  x.grid <- vector(mode = "numeric")
  for(i in 1:(length(region)/2)){
    x.grid <- c(x.grid, seq(region[2 * i - 1], region[ 2 * i], length.out = grid.length))
  }
  return(x.grid)
}


add_design <- function(design_1, design_2, alpha){
  design_2[, c("Weight")] <- design_2[, c("Weight")] / sum(design_2[, c("Weight")]) * alpha
  design_1[, c("Weight")] <- design_1[, c("Weight")] * (1 - alpha)
  return(rbind(design_1, design_2))
}


sum_design_point <- function(design, x, alpha){
  design$Weight <- design$Weight * (1 - alpha)
  design[nrow(design)+1, ] <- c(x, alpha)
  return(design)
}


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

get_augment_val <- function(KL_divergence, init_design, alpha, x){
  return(do.call(KL_divergence, list(sum_design_point(init_design, x, alpha))))
}


crossregion <- function(kl_div_func, init_design, alpha, eff, design_space, init_KL){
  kl_eff_func <- function(x){
    init_KL / get_augment_val(kl_div_func, init_design, alpha, x)
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
