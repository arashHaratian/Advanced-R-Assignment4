#' Ordinary Linear Regression
#'
#' This function calculates all the interesting parameters for the linear regression.
#'
#' @param formula is a `formula` object representing the model
#' @param data is a `data.frame` object containing the dataset that is being used
#'
#' @return an `RC` class object containing the following items:
#'\enumerate{
#'  \item Regressions Coefficients
#'  \item Fitted Values
#'  \item Residuals
#'  \item Degrees of Freedom
#'  \item Residual Variance
#'  \item Variance of the Regression Coefficients
#'  \item t-values
#'  \item p-values
#'}
#' @export
#'
#' @examples
#'
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
  return(list(
    beta,
    as.vector(fitted_val),
    as.vector(residual_val),
    df,
    res_var,
    var_reg,
    t_values,
    p_values
  ))
}

