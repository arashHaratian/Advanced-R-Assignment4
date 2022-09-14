linreg <- setRefClass(
  "linreg",
  fields = list(
    beta = "numeric",
    fitted_val = "numeric",
    residual_val = "numeric",
    df = "numeric",
    res_var = "numeric",
    var_reg = "numeric",
    t_values = "numeric",
    p_values = "numeric",
    call = "character"
  ),
  methods = list(
    initialize = function(formula, data){

      stopifnot(
        is.data.frame(data),
        inherits(formula, "formula")
      )

      X <- model.matrix(formula,data)
      dependent <- all.vars(formula)[1]

      n <- length(data[[dependent]])
      p <- length(colnames(X))

      beta_local <- solve(t(X) %*% X) %*% t(X) %*% data[[dependent]]
      .self$fitted_val <- as.vector(X %*% beta_local)
      residual_val_local <- data[[dependent]] - .self$fitted_val
      .self$residual_val <- as.vector(residual_val_local)
      .self$df <- n-p
      .self$res_var <- as.vector((t(residual_val_local) %*% residual_val_local)/.self$df)
      .self$var_reg <- diag(.self$res_var * solve(t(X) %*% X))
      .self$t_values <- as.vector(beta_local / sqrt(.self$var_reg))

      .self$p_values <- as.vector(pt(-abs(t_values),df = .self$df))

      beta_vector <- as.vector(beta_local)
      names(beta_vector) <- rownames(beta_local)
      .self$beta <- beta_vector

      formula_string <- as.character(formula)
      .self$call <- paste(
        "linreg(formula =",
        formula_string[[2]],
        "~",
        formula_string[[3]],
        ", data = ",
        deparse(substitute(data)),
        ")"
      )
    },

    print = function() {
      cat("Call:\n",
          format(.self$call),
          "\n\nCoefficients:\n")
      .self$beta
    },

    plot = function(){
      fig_1_data_frame <- data.frame(
        x_axis = .self$fitted_val,
        y_axis = .self$residual_val
      )
      fig_1 <- ggplot2::ggplot(data = fig_1_data_frame) +
        ggplot2::geom_point(ggplot2::aes(x = x_axis, y = y_axis)) +
        ggplot2::stat_summary(ggplot2::aes(x = x_axis, y = y_axis), fun = median, geom = "line") +
        ggplot2::labs(
          title = "Residuals vs Fitted",
          x ="Fitted Values",
          y = "Residuals"
        )



      fig_2_data_frame <- data.frame(
        x_axis = .self$fitted_val,
        y_axis = sqrt(
          abs(
            (.self$residual_val - mean(.self$residual_val)) / sd(.self$residual_val)
          )
        )
      )
      fig_2 <- ggplot2::ggplot(data = fig_2_data_frame) +
        ggplot2::geom_point(ggplot2::aes(x = x_axis, y = y_axis)) +
        ggplot2::stat_summary(ggplot2::aes(x = x_axis, y = y_axis), fun = mean, geom = "line") +
        ggplot2::labs(
          title = "Scale-Location",
          x = "Fitted Values",
          y = expression(sqrt(abs("Standardized Residuals")))
        )

      result <- list(
        "fig1" <- fig_1,
        "fig2" <- fig_2
      )
      return(result)
    },

    resid = function(){
      return(.self$residual_val)
    },

    pred = function(){
      return(.self$fitted_val)
    },

    coef = function(){
      return(.self$beta)
    },

    summary = function(){
      result <- cbind(
        "Estimate" = .self$beta,
        "Std. Error" = .self$var_reg,
        "t value" = .self$t_values,
        "Pr(>|t|)" = .self$p_values
      )

      cat("Residual standard error: ", .self$res_var," on ", .self$df, " degrees of freedom\n\n",
          "Coefficients:\n", sep = "")
      print.default(result)

      # print.default(result)
      # cat("Residual standard error: ", .self$res_var," on ", .self$df, " degrees of freedom\n\n",
      #     "Coefficients:\n", sep = "")


    }

  )


)
