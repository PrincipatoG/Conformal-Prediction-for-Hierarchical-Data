### Generate the structural matrix for the synthetical experiment ###
generate_S <- function(n = n, l = l, balanced=TRUE){
  exponent <- log(n,base = 12)  # Warning, this requires n to be a power of 12
  S <- rep(1,n)
  if(l==3){
    for (i in 1:(3**exponent)){
      S <- S %>% rbind(c(rep(0,(i-1)*(4**exponent)),rep(1,4**exponent),rep(0,(3**exponent-i)*(4**exponent)) ))
    }
    S <- S %>% rbind(diag(1,n))
  }
  if(l==4){
    for (i in 1:(2**exponent)){
      S <- S %>% rbind(c(rep(0,(i-1)*(6**exponent)),rep(1,6**exponent),rep(0,(2**exponent-i)*(6**exponent)) ))
    }
    for (i in 1:(4**exponent)){
      S <- S %>% rbind(c(rep(0,(i-1)*(3**exponent)),rep(1,3**exponent),rep(0,(4**exponent-i)*(3**exponent)) ))
    }
    S <- S %>% rbind(diag(1,n))
  }
  return(S)
}
### Generate the observation (X,Y) from the most disaggregated level and reconstruct the entire hierarchy ###
generate_multi <- function(S,n_obs){
  n <- ncol(S)
  # Generate the features
  X1 <- fl(rnorm(n_obs, mean = 10, sd = 2))
  X2 <- fl(rnorm(n_obs, mean = -5, sd = 2))
  X3 <- fl(rnorm(n_obs, mean = 5, sd = 1))
  
  # Define the base explainatory functions
  base_functions <- list(
    function(X1, X2, X3) { X1 },
    function(X1, X2, X3) { X2 },
    function(X1, X2, X3) { X3 },
    function(X1, X2, X3) { X1^2 },
    function(X1, X2, X3) { X2^2 },
    function(X1, X2, X3) { X3^2 },
    function(X1, X2, X3) { sin(X1) },
    function(X1, X2, X3) { cos(X2) },
    function(X1, X2, X3) { exp(X3) },
    function(X1, X2, X3) { log(abs(X1) + 1) },
    function(X1, X2, X3) { sqrt(abs(X2)) }
  )
  num_base_functions <- length(base_functions)
  # Generate the observations Y_i
  Y <- fl(matrix(0, n_obs, n))
  for (i in 1:n) {
    k <- sample(1:num_base_functions, 1)
    selected_functions <- sample(base_functions, k, replace = TRUE)
    weights <- 2*(rbernoulli(k, p = 0.5)-0.5)
    combined_function <- function(X1, X2, X3) {
      result <- 0
      for (j in 1:k) {
        result <- result + weights[j] * selected_functions[[j]](X1, X2, X3)
      }
      return(result)
    }
    Y[, i] <- fl(combined_function(X1, X2, X3)) 
  }
  # Add a multivariate noise to make the forecasting task more realistic 
  if (noise_type == "B"){
    characeristics <- readRDS("...")
    noise <- fl(rmvt(n_obs, mu = characeristics[[1]], sigma = characeristics[[2]], df = 3))
  }else{
    # Use Cholesky to generate a random covariance matrix
    generate_random_corr_matrix <- function(n) {
      A <- fl(matrix(rnorm(n^2), nrow=n))
      cov_matrix <- A %*% t(A)
      D <- diag(1 / sqrt(diag(cov_matrix)))
      corr_matrix <- D %*% cov_matrix %*% D
      return(corr_matrix)
    }
    Sigma <- generate_random_corr_matrix(n)
    noise <- mvrnorm(n_obs, mu = rep(10,n), Sigma = 100*Sigma)
  }
  Y_bruite <- fl(Y + noise)
  # Reconstruct the multivariate data series
  b <- as.matrix(t(Y_bruite))
  Y_level <- as.matrix(S %*% b) 
  rownames(Y_level) <- paste0("V", seq_len(nrow(Y_level)))
  Features <- cbind(X1,X2,X3)
  return(list(fl(Features),fl(Y_level)))
}  
### Compute the Coverage and Length from the simultaneous individual SCP ###
conformal <- function(data_calib, data_test, n, m, method, alpha = 0.1){
  info_interval <- data.frame(metric = c("Coverage", "Length", "Quantile"))
  for (i in 1:m){
    y_i_calib <- data_calib[[paste0("V", i)]]
    y_i_test <- data_test[[paste0("V", i)]]
    if (method == "Direct") {
      y_tilde_i_calib <- data_calib[[paste0("V", i, "_hat")]]
      y_tilde_i_test <- data_test[[paste0("V", i, "_hat")]]
    }else {
      y_tilde_i_calib <- data_calib[[paste0("V", i, "_",method)]]
      y_tilde_i_test <- data_test[[paste0("V", i, "_",method)]]
    }  
    score <- y_i_calib-y_tilde_i_calib
    q_Lower <- -quantile(-score,  (length(score)+1)*(1-alpha/2)/(length(score) ), type = 1)
    q_Upper <- quantile(score,  (length(score)+1)*(1-alpha/2)/(length(score) ), type = 1) 
    # The Length is the squared length (to correspond to the objective of the article)
    Length <- (q_Upper - q_Lower)^2
    C <- cbind(y_tilde_i_test + q_Lower, y_tilde_i_test + q_Upper)
    Coverage <- mean( (C[,1] <= y_i_test) & (y_i_test <= C[,2]) )  
    new_col <- data.frame(V = c(Coverage, Length, 1-alpha))
    names(new_col) <- c(paste0("V", i))
    info_interval <- info_interval %>% cbind(new_col)
  }
  return(info_interval)
}
### Compute the 2m-th root of the determinant of A in a robust fashion  ###
root_det <- function(A) {
  m <- nrow(A)
  d <- determinant(A, logarithm = TRUE) 
  sign <- d$sign
  modulus <- d$modulus
  root_modulus <- exp(as.numeric(modulus) / (2*m))
  return(sign * root_modulus)
}
### Compute the Coverage and Volume from the Joint SCP ###
joint_conformal <- function(data_calib, Sigma, data_test, n, m, method, alpha = 0.1){
  info_interval <- data.frame(metric = c("Coverage", "Volume", "Quantile"))
    y_calib <- data_calib[,paste0("V", 1:m)]
    y_test <- data_test[,paste0("V", 1:m)]
    if (method == "Direct") {
      y_tilde_calib <- data_calib[,paste0("V", 1:m, "_hat")]
      y_tilde_test <- data_test[,paste0("V", 1:m, "_hat")]
    }else {
      y_tilde_calib <- data_calib[,paste0("V", 1:m, "_",method)]
      y_tilde_test <- data_test[,paste0("V", 1:m, "_",method)]
    }  
    print("Computing the calibration scores")
    score <- rowSums( (as.matrix(y_calib-y_tilde_calib) %*% ginv(Sigma) ) * as.matrix(y_calib-y_tilde_calib) )
    q_joint <- quantile(score,  (length(score)+1)*(1-alpha)/(length(score) ), type = 1) 
    determinant_sigma <- root_det(Sigma)
    Volume <- sqrt(q_joint)*determinant_sigma
    print("Computing the test scores")
    score_test <- rowSums( (as.matrix(y_test-y_tilde_test) %*% ginv(Sigma) ) * as.matrix(y_test-y_tilde_test) )
    Coverage <- mean( (score_test <= q_joint) )  
    new_col <- data.frame(V = c(Coverage, Volume, 1-alpha))
    info_interval <- info_interval %>% cbind(new_col)
  return(info_interval)
}
### Synthesize the metric from the functions conformal and joint_conformal ###
metrics_proba_multi <- function(data_calib, Sigma, data_test, n, m,  methods, alpha = 0.1) {
  data_list <- list()  
  joint_results <- list()
  for (method in methods) {
    print("First the component-wise coverage")
    info_interval <- conformal(data_calib, data_test, n, m, method, alpha = alpha)
    info_interval$method <- method
    data_list[[method]] <- info_interval
    print("Then the joint coverage")
    joint_results[[method]] <- joint_conformal(data_calib, Sigma, data_test, n, m, method, alpha = alpha)
  }
  combined_data <- do.call(rbind, data_list)
  long_data <- combined_data %>%
    pivot_longer(cols = starts_with("V"), names_to = "variable", values_to = "value")
  return(list(long_data,joint_results))
}
### Perform the simulation ###
simulation_multi <- function(S, n_obs, model){
  alpha <- 0.1
  # Split the available data
  train <- seq(1,floor((0.4)*n_obs))
  valid <- seq(floor((0.4)*n_obs)+1, floor((0.6)*n_obs))
  calib <- seq(floor((0.6)*n_obs)+1, floor((0.8)*n_obs))
  test  <- seq(floor((0.8)*n_obs)+1, n_obs)
  # Generate the data
  data  <- generate_multi(S, n_obs = n_obs)
  Features <- data[[1]]
  Y <- data[[2]]
  print("Data Generated Successfully")
  # Training of the base forecasts
  m <- nrow(Y)
  data <- as.data.frame(cbind(t(Y),Features)) 
  colnames(data) <- c(paste0("V",1:m),"X1","X2","X3")
  prediction_fast <- function(data, model = "gam", train, m){
  forecasts <- fl(matrix(0, nrow(data)-length(train), m))
  for (i in 1:m){
      print(paste0("i = ",i))
      # Start with the top levels
      if(i < m-n+1){
        if (model == "gam"){
          formula_i <- as.formula(paste0("V",i,"~s(X1)+s(X2)+s(X3)"))            
          # model_i   <- bam(formula_i,data = training_data, sp = c(1,1,1), discrete = TRUE)
          model_i   <- bam(formula_i, data = data[train,c(paste0("V",i),"X1","X2","X3")], sp = c(1,1,1), discrete = TRUE)
          model_i$model <- NULL
        }
        if (model =="lm"){
          formula_i <- as.formula(paste0("V",i,"~X1+X2+X3"))
          model_i <- lm(formula_i, data = data[train,c(paste0("V",i),"X1","X2","X3")])
        }
        forecasts[,i] <- fl(predict(model_i,data[-train,]))
      }else{
        # For the leaves, we force some of the prediction tasks to be harder by hidding the explanatory variable.
        if (model == "gam"){
          if(rbernoulli(1,0.8)){ # The number of 'bad prediction' has to be tuned to illustrate different behaviors
            formula_i <- as.formula(paste0("V",i,"~s(X1)+s(X2)+s(X3)"))  
            model_i   <- bam(formula = formula_i,data = data[train,c(paste0("V",i),"X1","X2","X3")],sp = c(1,1,1), discrete = TRUE)
          }else{
            formula_i <- as.formula(paste0("V",i,"~s(X1)+s(X2)")) 
            model_i   <- bam(formula = formula_i,data = data[train,c(paste0("V",i),"X1","X2")],sp = c(1,1), discrete = TRUE)
          }
        }
        if (model =="lm"){
          # For the leaves, we force some of the prediction tasks to be harder by hidding the explanatory variable.
          if(rbernoulli(1,0.8)){
            formula_i <- as.formula(paste0("V",i,"~X1+X2+X3"))
          }else{
            formula_i <- as.formula(paste0("V",i,"~X1+X2"))
          }
          model_i <- lm(formula_i, data = data[train,c(paste0("V",i),"X1","X2","X3")])
        }
        forecasts[,i] <- fl(predict(model_i,data[-train,]))
      }
    }
    return(forecasts)
  }
  forecasts <- prediction_fast(data, model, train, m)
  colnames(forecasts) <- paste0(paste0("V",1:m,"_hat"))
  # Reconciliation Step 
  names   <- c(paste0("V",1:m))
  W_OLS <- diag(m)
  print("Compute the reconciliation matrices")
  error <- as.matrix(fl(data[valid,1:m])-(forecasts[valid-length(train),]))
  W_MinT <- cov(as.numeric(error))
  W_WLS <- diag(diag(W_MinT))
  G_MinT  <- ginv(t(S)%*%ginv(W_MinT)%*%S)%*%t(S)%*%ginv(W_MinT)
  P_MinT  <- S %*%G_MinT
  G_OLS <- ginv(t(S)%*%ginv(W_OLS)%*%S)%*%t(S)%*%ginv(W_OLS)
  P_OLS <- S %*% G_OLS
  ## Uncomment if you want to consider Bottom-Up Reconciliation
  # G_BT    <- matrix(0,ncol = m-n, nrow = n) %>% cbind(diag(n))
  # P_BT    <- S %*% G_BT
  G_WLS <- ginv(t(S)%*%ginv(W_WLS)%*%S)%*%t(S)%*%ginv(W_WLS)
  P_WLS <- S %*% G_WLS
  G_Combi <- (1/3)*(G_MinT+G_OLS+G_WLS)
  P_Combi <- S %*% G_Combi
  projections <- list(MinT = P_MinT, WLS = P_WLS, OLS = P_OLS, Combi = P_Combi)
  methods <- names(projections)
  print("Reconcile the forecasts")
  data_calib <- cbind(data[calib,], as.numeric(as.matrix(forecasts[calib-length(train),])))
  data_test <- cbind(data[test,], as.numeric(as.matrix(forecasts[test-length(train),])))
  data <- NULL
  for (method in methods) {
    projection_matrix <- projections[[method]]
    predicteurs_calib <- t(as.numeric(as.matrix(forecasts[calib-length(train),])))
    predicteurs_test <- t(as.numeric(as.matrix(forecasts[test-length(train),])))
    transformed_data_calib  <- as.data.frame(t(projection_matrix %*% predicteurs_calib))
    transformed_data_test <- as.data.frame(t(projection_matrix %*% predicteurs_test))
    colnames(transformed_data_calib) <- paste0(names, '_', method)    
    colnames(transformed_data_test) <- paste0(names, '_', method)
    data_calib <- cbind(data_calib, transformed_data_calib)
    data_test <- cbind(data_test, transformed_data_test)
  }
  methods <- c("Direct","MinT","WLS","OLS","Combi")
  # Compute the probabilistic metrics
  print("Run the Conformal step")
  result_conformal <- metrics_proba_multi(data_calib, Sigma = W_MinT, data_test, n, m, methods)
  return(result_conformal)
}
### Summarize the result in a LaTeX table ###
get_table <- function(data, uncertainty, bold = TRUE, digits = 3) {
  latex_rows <- sapply(seq_len(nrow(data)),  function(i) {
    row <- data[i, ]
    row_unc <- uncertainty[i, ]
    # Look for the minimum for each line
    min_val <- min(row, na.rm = TRUE)
    formatted_row <- mapply(function(x, u) {
      value_str <- paste0(format(x, digits = digits), " $\\pm$ ", format(u, digits = 2))
      if (x == min_val && bold){
        return(paste0("\\textbf{", value_str, "}"))
      } else {
        return(value_str)
      }
    }, row, row_unc)
    # Gather everything and Add a configuration column
    line_with_config <- paste(c(paste0(i), formatted_row), collapse = " & ")
    return(line_with_config)
  })
  col_names <- paste(c("Configuration", colnames(data)), collapse = " & ")
  body <- paste(latex_rows, collapse = " \\\\\n")
  table_latex <- paste0(
    "\\begin{adjustbox}{width=\\textwidth}\n",
    "\\begin{tabular}{", paste(rep("c", ncol(data) + 1), collapse = ""), "}\n",
    col_names, " \\\\\n\\hline\n",
    body, " \\\\\n",
    "\\end{tabular}\n",
    "\\end{adjustbox} "
  )
  return(table_latex)
  # The LaTeX table can be displayed with the cat function
}
### Summarize the Occurence of a method being the best in a LaTeX table ###
get_table_occurence <- function(data) {
  n_total <- nrow(data)
  # The number of layers is either 3 or 4
  config_sizes <- rep(c(3, 4), length.out = ceiling(n_total / 4))
  config_sizes <- config_sizes[cumsum(config_sizes) <= n_total]
  if (sum(config_sizes) < n_total) {
    config_sizes <- c(config_sizes, n_total - sum(config_sizes))
  }
  latex_rows <- c()
  row_index <- 1
  config_number <- 1
  # Look for the maximum for each line
  for (size in config_sizes) {
    for (j in 1:size) {
      i <- row_index
      row <- data[i, ]
      max_val <- max(row, na.rm = TRUE)
      formatted_row <- mapply(function(x) {
        value_str <- paste0(format(x, digits = 3))
        if (!is.na(x) && x == max_val) {
          return(paste0("\\textbf{", value_str, "}"))
        } else {
          return(value_str)
        }
      }, row)
      config_cell <- if (j == 1) {
        paste0("\\multirow{", size, "}{*}{Config ", config_number, "}")
      } else {
        ""
      }
      layer_cell <- paste0("L", j)
      # Gather everything and Add a configuration column
      line <- paste(c(config_cell, layer_cell, formatted_row), collapse = " & ")
      latex_rows <- c(latex_rows, line)
      row_index <- row_index + 1
    }
    if (config_number < length(config_sizes)) {
      latex_rows <- c(latex_rows, "\\addlinespace")
    }
    config_number <- config_number + 1
  }
  col_names <- paste(c("Configuration", "Level", colnames(data)), collapse = " & ")
  body <- paste(latex_rows, collapse = " \\\\\n")
  table_latex <- paste0(
    "\\begin{adjustbox}{width=\\textwidth}\n",
    "\\begin{tabular}{c|", paste(rep("c", ncol(data) + 1), collapse = " "), "}\n",
    "\\toprule\n",
    col_names, " \\\\\n\\midrule\n",
    body, " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{adjustbox} "
  )
  return(table_latex)
  # The LaTeX table can be displayed with the cat function
}
