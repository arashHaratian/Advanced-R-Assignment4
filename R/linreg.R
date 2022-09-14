#' Ordinary Linear Regression
#'
#' This function calculates all the interesting parameters for the linear regression.
#'
#' @param formula is a `formula` object representing the model
#' @param data is a `data.frame` object containing the dataset that is being used
#'
#' @return an `RC` class object containing the following items:
#' \enumerate{
#'  \item Regressions Coefficients
#'  \item Fitted Values
#'  \item Residuals
#'  \item Degrees of Freedom
#'  \item Residual Variance
#'  \item Variance of the Regression Coefficients
#'  \item t-values
#'  \item p-values
#' }
#' @export
#'
#' @examples

linreg <- function(formula,data){
  X <- model.matrix(formula,data)
  dependent <- all.vars(formula)[1]

  n <- length(data[[dependent]])
  p <- length(colnames(X))

  beta <- solve(t(X) %*% X) %*% t(X) %*% data[[dependent]]
  fitted_val <- X %*% beta
  residual_val <- data[[dependent]] - fitted_val
  df <- n-p
  res_var <- as.vector((t(residual_val) %*% residual_val)/df)
  var_reg <- diag(res_var * solve(t(X) %*% X))
  t_values <- beta / sqrt(var_reg)

  p_values <- pt(-abs(t_values),df=df)

  beta_vector <- as.vector(beta)
  names(beta_vector) <- rownames(beta)

  formula_string <- as.character(formula)
  call <- paste(
    "linreg(formula =",
    formula_string[[2]],
    "~",
    formula_string[[3]],
    ", data = ",
    deparse(substitute(data)),
    ")"
  )

  result <- linreg_RC$new(
    "call" = match.call(),
    "beta" = beta_vector,
    "fitted_val" = as.vector(fitted_val),
    "residual_val" = as.vector(residual_val),
    "df" = df,
    "res_var" = res_var,
    "var_reg" = var_reg,
    "t_values" = t_values,
    "p_values" = p_values
  )
  return(result)
}

