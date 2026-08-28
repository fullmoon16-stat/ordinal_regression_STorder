
## baseline survival functions ##
make_surv_func <- function(m){
  
  list(
    
    ###### Simple baseline Sfs ######
    
    S1 = list(
      func  = function(t) punif(t, 0, m, lower.tail = FALSE),
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    S2 = list(
      func  = function(t) (1 - (t / m))^2,
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    S3 = list(
      func  = function(t) sqrt(1 - (t / m)),
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    S4 = list(
      func  = function(t) plogis(-2 * qlogis(t / m)),
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    S5 = list(
      func  = function(t) plogis(-0.5 * qlogis(t / m)),
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    ##### S^{*}(t) #####
    
    H1 = list(
      func  = function(t, log_k) {
        k <- exp(log_k)
        1 - (t / m)^k
      },
      npar  = 1,
      init  = c(log_k = log(0.5)),
      lower = c(log_k = -3),
      upper = c(log_k = 3)
    ),
    
    H2 = list(
      func  = function(t, log_k) {
        k <- exp(log_k)
        (1 - (t / m))^k
      },
      npar  = 1,
      init  = c(log_k = log(2)),
      lower = c(log_k = -3),
      upper = c(log_k = 3)
    ),
    
    H3 = list(
      func  = function(t, k) plogis(k - qlogis(t / m)),
      npar  = 1,
      init  = c(k = -1),
      lower = c(k = -3),
      upper = c(k = 3)
    ),
    
    H4 = list(
      func  = function(t, k, log_s) {
        s <- exp(log_s)
        plogis(k - s * qlogis(t / m))
      },
      npar  = 2,
      init  = c(k = -1, log_s = log(1)),
      lower = c(k = -3, log_s = -3),
      upper = c(k = 3,  log_s = 3)
    ),
    
    ###### Dist ######
    
    S_beta = list(
      func  = function(t, log_a, log_b) {
        a <- exp(log_a)
        b <- exp(log_b)
        Rbeta(1 - (t / m), b, a, lower = TRUE)
      },
      npar  = 2,
      init  = c(log_a = log(0.5), log_b = log(1.7)),
      lower = c(log_a = -3, log_b = -3),
      upper = c(log_a = 3,  log_b = 3)
    )
    
  )
}


## transformations of the linear predictor g(eta) ##
geta_candidates <- list(
  inc_1 = list(func = function(eta) exp(eta)),
  zero_1 = list(func = function(eta) eta)
)

## covariate-dependent transformation functions rho(x;j) ##
make_rhoj_candidates <- function() {
  list(
    
    #### rho_1(x;j) ####
    rho_11 = list(
      fun = function(eta_val, j, m, ...) {
        m * (j / m)^eta_val
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_12 = list(
      fun = function(eta_val, j, m, ...) {
        m * (1 - (1 - j / m)^eta_val)
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_13 = list(
      fun = function(eta_val, j, m, ...) {
        m * plogis(qlogis(j / m) - eta_val)
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_14 = list(
      fun = function(eta_val, j, m, log_c, ...) {
        c <- exp(log_c)
        m * plogis(-eta_val + c * qlogis(j / m))
      },
      npar  = 1,
      init  = c(log_c = log(1)),
      lower = c(log_c = -3),
      upper = c(log_c = 3)
    ),
    
    #### rho_2(x;j) ####
    rho_21 = list(
      fun = function(eta_val, j, m, ...) {
        m * (1 - sqrt(1 - (j / m)^eta_val))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_22 = list(
      fun = function(eta_val, j, m, ...) {
        m * (1 - sqrt((1 - j / m)^eta_val))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_23 = list(
      fun = function(eta_val, j, m, ...) {
        m * (1 - sqrt(plogis(eta_val - qlogis(j / m))))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_24 = list(
      fun = function(eta_val, j, m, log_c, ...) {
        c <- exp(log_c)
        m * (1 - sqrt(plogis(eta_val - c * qlogis(j / m))))
      },
      npar  = 1,
      init  = c(log_c = log(2.5)),
      lower = c(log_c = -3),
      upper = c(log_c = 3)
    ),
    
    #### rho_3(x;j) ####
    rho_31 = list(
      fun = function(eta_val, j, m, ...) {
        m * (1 - (1 - (j / m)^eta_val)^2)
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_32 = list(
      fun = function(eta_val, j, m, ...) {
        m * (1 - ((1 - j / m)^eta_val)^2)
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_33 = list(
      fun = function(eta_val, j, m, ...) {
        m * (1 - plogis(eta_val - qlogis(j / m))^2)
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_34 = list(
      fun = function(eta_val, j, m, log_c, ...) {
        c <- exp(log_c)
        m * (1 - plogis(eta_val - c * qlogis(j / m))^2)
      },
      npar  = 1,
      init  = c(log_c = log(1)),
      lower = c(log_c = -3),
      upper = c(log_c = 3)
    ),
    
    #### rho_4(x;j) ####
    rho_41 = list(
      fun = function(eta_val, j, m, ...) {
        m * plogis(-0.5 * qlogis(1 - (j / m)^eta_val))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_42 = list(
      fun = function(eta_val, j, m, ...) {
        m * plogis(-0.5 * qlogis((1 - j / m)^eta_val))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_43 = list(
      fun = function(eta_val, j, m, ...) {
        m * plogis(-0.5 * (eta_val - qlogis(j / m)))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_44 = list(
      fun = function(eta_val, j, m, log_c, ...) {
        c <- exp(log_c)
        m * plogis(-0.5 * (eta_val - c * qlogis(j / m)))
      },
      npar  = 1,
      init  = c(log_c = log(1)),
      lower = c(log_c = -3),
      upper = c(log_c = 3)
    ),
    
    #### rho_5(x;j) ####
    rho_51 = list(
      fun = function(eta_val, j, m, ...) {
        m * plogis(-2 * qlogis(1 - (j / m)^eta_val))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_52 = list(
      fun = function(eta_val, j, m, ...) {
        m * plogis(-2 * qlogis((1 - j / m)^eta_val))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_53 = list(
      fun = function(eta_val, j, m, ...) {
        m * plogis(-2 * (eta_val - qlogis(j / m)))
      },
      npar  = 0,
      init  = numeric(0),
      lower = numeric(0),
      upper = numeric(0)
    ),
    
    rho_54 = list(
      fun = function(eta_val, j, m, log_c, ...) {
        c <- exp(log_c)
        m * plogis(-2 * (eta_val - c * qlogis(j / m)))
      },
      npar  = 1,
      init  = c(log_c = log(1)),
      lower = c(log_c = -3),
      upper = c(log_c = 3)
    )
    
  )
}


## split the optimization parameter vector ##
split_par <- function(par, p, qS, qJ, pT = 0, m = NULL) {
  idx <- 1
  
  ## common beta
  beta <- par[idx:(idx + p - 1)]
  idx <- idx + p
  
  ## gamma_2, ..., gamma_{m-1}
  if (pT > 0) {
    if (is.null(m)) stop("m must be supplied when pT > 0.")
    
    qG <- pT * (m - 2)
    
    gamma_free <- par[idx:(idx + qG - 1)]
    idx <- idx + qG
    
    Gamma <- matrix(0, nrow = pT, ncol = m - 1)
    Gamma[, 2:(m - 1)] <- matrix(gamma_free, nrow = pT, ncol = m - 2)
    # Gamma 는 (0,gamma1,gamma2,...) 미리 만들어놓는 용도임
    
  } else {
    gamma_free <- numeric(0)
    Gamma <- matrix(numeric(0), nrow = 0, ncol = ifelse(is.null(m), 0, m - 1))
  }
  
  ## baseline survival parameters
  surv_par <- if (qS > 0) par[idx:(idx + qS - 1)] else numeric(0)
  idx <- idx + qS
  
  ## rhoj parameters
  rhoj_par <- if (qJ > 0) par[idx:(idx + qJ - 1)] else numeric(0)
  
  list(
    beta       = beta,
    gamma_free = gamma_free,
    Gamma      = Gamma,
    surv_par   = surv_par,
    rhoj_par   = rhoj_par
  )
}


## conditional pmf induced by S(rho(eta_j;j)) ##
pmf_from_rhoj <- function(S_base, eta_val, y, m, rhoj_fun,
                           surv_par = numeric(0),
                           rhoj_par = numeric(0)) {
  
  y <- as.numeric(y)
  n <- length(y)
  j_grid <- 1:(m - 1)
  
  ## eta_val: vector이면 proportional version
  ## eta_val: matrix이면 PPOM-like version
  if (is.null(dim(eta_val))) {
    if (length(eta_val) == 1 && n > 1) {
      eta_val <- rep(eta_val, n)
    }
    eta_mat <- matrix(eta_val, nrow = n, ncol = m - 1)
  } else {
    eta_mat <- eta_val
    if (nrow(eta_mat) != n) stop("nrow(eta_val) must equal length(y).")
    if (ncol(eta_mat) != m - 1) stop("ncol(eta_val) must be m - 1.")
  }
  
  ## Q_j = P(Y > j | x), j = 1,...,m-1
  Q_mid <- sapply(j_grid, function(jval) {
    t_j <- do.call(
      rhoj_fun,
      c(list(eta_val = eta_mat[, jval], j = jval, m = m),
        as.list(rhoj_par))
    )
    
    do.call(
      S_base,
      c(list(t = t_j), as.list(surv_par))
    )
  })
  
  ## Q_0 = 1, Q_m = 0
  Q <- cbind(1, Q_mid, 0)
  
  ## p_y = Q_{y-1} - Q_y
  p <- Q[cbind(seq_len(n), y)] - Q[cbind(seq_len(n), y + 1)]
  
  p
}


## fit one baseline-rho combination by maximum likelihood ##
fit_one_combo <- function(X, y, m, n,
                           surv_info, eta_info, rhoj_info,
                           T_np = NULL,
                           method = "BFGS") {
  
  p  <- ncol(X)
  pT <- if (is.null(T_np)) 0 else ncol(T_np)
  
  qS <- surv_info$npar
  qJ <- rhoj_info$npar
  qG <- pT * (m - 2)
  
  init_par <- c(
    rep(0, p),
    rep(0, qG),
    surv_info$init,
    rhoj_info$init
  )
  
  nerho_loglik <- function(par) {
    sp <- split_par(par, p = p, qS = qS, qJ = qJ, pT = pT, m = m)
    
    eta_common <- drop(X %*% sp$beta)
    if (any(!is.finite(eta_common))) return(1e12)
    
    if (pT > 0) {
      eta <- matrix(eta_common, nrow = n, ncol = m - 1) + T_np %*% sp$Gamma
    } else {
      eta <- eta_common
    }
    
    if (any(!is.finite(eta))) return(1e12)
    
    eta_val <- eta_info$func(eta)
    if (any(!is.finite(eta_val))) return(1e12)
    
    probs <- pmf_from_rhoj(
      S_base   = surv_info$func,
      eta_val  = eta_val,
      y        = y,
      m        = m,
      rhoj_fun = rhoj_info$fun,
      surv_par = sp$surv_par,
      rhoj_par = sp$rhoj_par
    )
    
    if (any(!is.finite(probs))) return(1e12)
    if (any(probs <= 0)) return(1e12)
    
    probs <- pmax(pmin(probs, 1 - 1e-12), 1e-12)
    -sum(log(probs))
  }
  
  opt <- optim(par = init_par, fn = nerho_loglik,
               method = method, hessian = TRUE)
  
  opt <- optim(par = opt$par, fn = nerho_loglik,
               method = method, hessian = TRUE)
  
  opt <- optim(par = opt$par, fn = nerho_loglik,
               method = method, hessian = TRUE)
  
  H <- opt$hessian
  H <- 0.5 * (H + t(H))
  
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V) || any(diag(V) <= 0)) {
    se_raw <- rep(NA_real_, length(opt$par))
  } else {
    se_raw <- sqrt(diag(V))
  }
  
  sp_hat <- split_par(opt$par, p = p, qS = qS, qJ = qJ, pT = pT, m = m)
  
  idx_beta  <- seq_len(p)
  idx_gamma <- if (qG > 0) (p + 1):(p + qG) else integer(0)
  
  list(
    par        = opt$par,
    beta       = sp_hat$beta,
    gamma_free = sp_hat$gamma_free,
    Gamma      = sp_hat$Gamma,
    surv_par   = sp_hat$surv_par,
    rhoj_par   = sp_hat$rhoj_par,
    se_raw     = se_raw,
    se_beta    = se_raw[idx_beta],
    se_gamma   = if (qG > 0) se_raw[idx_gamma] else numeric(0),
    V          = V,
    H          = H,
    nll        = opt$value,
    logLik     = -opt$value,
    aic        = 2 * opt$value + 2 * length(opt$par),
    bic        = 2 * opt$value + log(n) * length(opt$par),
    conv       = opt$convergence
  )
}


## initial values based on the empirical survival function ##
get_data_driven_init <- function(y, m, surv_info, eta_info, rhoj_info) {
  
  qS <- surv_info$npar
  qJ <- rhoj_info$npar
  
  # empirical survival
  empS <- sapply(0:m, function(j) mean(y > j))
  
  # beta = 0 -> eta = 0
  eta0 <- 0
  j_grid <- 1:(m - 1)
  
  # g(0)
  rho0 <- eta_info$func(eta0)
  
  # 최적화할 추가모수가 전혀 없으면 그대로 반환
  if (qS == 0 && qJ == 0) {
    t_j <- do.call(rhoj_info$fun, c(list(eta_val = rho0, j = j_grid, m = m)))
    S_j <- do.call(surv_info$func, c(list(t = t_j)))
    
    return(list(
      surv_init = numeric(0),
      rhoj_init = numeric(0),
      obj_value = sum((S_j - empS[j_grid + 1])^2),
      conv = 0
    ))
  }
  
  init_par <- c(surv_info$init, rhoj_info$init)
  lower    <- c(surv_info$lower, rhoj_info$lower)
  upper    <- c(surv_info$upper, rhoj_info$upper)
  
  obj_fun <- function(par) {
    idx <- 1
    
    surv_par <- if (qS > 0) par[idx:(idx + qS - 1)] else numeric(0)
    idx <- idx + qS
    
    rhoj_par <- if (qJ > 0) par[idx:(idx + qJ - 1)] else numeric(0)
    
    t_j <- tryCatch(
      do.call(
        rhoj_info$fun,
        c(list(eta_val = rho0, j = j_grid, m = m), as.list(rhoj_par))
      ),
      error = function(e) rep(NaN, length(j_grid))
    )
    
    if (any(!is.finite(t_j))) return(1e12)
    
    S_j <- tryCatch(
      do.call(
        surv_info$func,
        c(list(t = t_j), as.list(surv_par))
      ),
      error = function(e) rep(NaN, length(j_grid))
    )
    
    if (any(!is.finite(S_j))) return(1e12)
    
    # survival 값의 기본 제약 확인
    if (any(S_j < 0 | S_j > 1)) return(1e12)
    
    sum((S_j - empS[j_grid + 1])^2)
  }
  
  # box constraint 있으므로 L-BFGS-B 사용
  opt <- optim(
    par    = init_par,
    fn     = obj_fun,
    method = "L-BFGS-B",
    lower  = lower,
    upper  = upper
  )
  
  par_hat <- opt$par
  idx <- 1
  
  surv_init <- if (qS > 0) par_hat[idx:(idx + qS - 1)] else numeric(0)
  idx <- idx + qS
  
  rhoj_init <- if (qJ > 0) par_hat[idx:(idx + qJ - 1)] else numeric(0)
  
  list(
    surv_init = surv_init,
    rhoj_init = rhoj_init,
    obj_par = opt$par,
    conv = opt$convergence
  )
}

## Select g(eta) according to the rho function ##
get_eta_name_for_rhoj <- function(rhoj_name) {
  
  # rhoj_name examples: "rho_11", "rho_12", ..., "rho_54"
  rho_type <- as.integer(substr(rhoj_name, nchar(rhoj_name), nchar(rhoj_name)))
  
  if (rho_type %in% c(1, 2)) {
    return("inc_1")   # g(eta) = exp(eta)
  }
  
  if (rho_type %in% c(3, 4)) {
    return("zero_1")  # g(eta) = eta
  }
  
  stop("Unknown rhoj type in rhoj_name: ", rhoj_name)
}

## Utility: return y when x is NULL ##
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


## Parse rhoj name such as "rho_11", "rho_24", ..., "rho_54" ##
parse_rhoj_name <- function(rhoj_name) {
  
  rhoj_name <- as.character(rhoj_name)
  
  out <- regexec("^rho_([0-9]+)([0-9]+)$", rhoj_name)
  parsed <- regmatches(rhoj_name, out)[[1]]
  
  list(
    rho_l = as.integer(parsed[2]),
    rho_a = as.integer(parsed[3])
  )
}


## Convert numeric vector to a compact character string ##
compact_num <- function(x, digits = 4) {
  if (is.null(x) || length(x) == 0) return("")
  paste(round(as.numeric(x), digits = digits), collapse = ", ")
}


## Convert all_results list to a data frame ##
all_results_to_df <- function(all_results, digits = 4) {
  
  rows <- lapply(names(all_results), function(key) {
    
    res <- all_results[[key]]
    
    rhoj_name <- res$rhoj
    parsed <- parse_rhoj_name(rhoj_name)
    
    data.frame(
      key      = key,
      surv     = as.character(res$surv),
      geta     = as.character(res$geta),
      rhoj     = as.character(res$rhoj),
      rho_l    = parsed$rho_l,
      rho_a    = parsed$rho_a,
      logLik   = -as.numeric(res$nll),
      AIC      = as.numeric(res$aic),
      BIC      = as.numeric(res$bic),
      conv     = as.integer(res$conv %||% NA_integer_),
      beta     = compact_num(res$beta, digits = digits),
      gamma    = compact_num(res$gamma_free, digits = digits),
      surv_par = compact_num(res$surv_par, digits = digits),
      rhoj_par = compact_num(res$rhoj_par, digits = digits),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}



## Select the best fitted model for each rho_a based on AIC or BIC ##
select_best_by_rhoa <- function(all_results,
                                 criterion = "BIC",
                                 require_convergence = TRUE,
                                 digits = 4,
                                 tie_tol = 1e-6) {
  
  criterion <- toupper(criterion)
  
  result_df <- all_results_to_df(all_results, digits = digits)
  
  if (require_convergence) {
    result_df <- result_df[result_df$conv == 0, , drop = FALSE]
  }
  
  best_list <- lapply(split(result_df, result_df$rho_a), function(df_a) {
    
    df_a <- df_a[is.finite(df_a[[criterion]]), , drop = FALSE]
    
    if (nrow(df_a) == 0) return(NULL)
    
    min_val <- min(df_a[[criterion]], na.rm = TRUE)
    
    df_tie <- df_a[
      abs(df_a[[criterion]] - min_val) <= tie_tol,
      ,
      drop = FALSE
    ]
    
    df_tie <- df_tie[
      order(df_tie$rho_l, df_tie$surv, df_tie$rhoj),
      ,
      drop = FALSE
    ]
    
    df_tie[1, , drop = FALSE]
  })
  
  best_df <- do.call(rbind, best_list)
  rownames(best_df) <- NULL
  
  best_df <- best_df[order(best_df$rho_a), , drop = FALSE]
  
  best_df[, c(
    "rho_a", "surv", "geta", "rhoj", "rho_l",
    "logLik", "AIC", "BIC", "conv",
    "beta", "gamma", "surv_par", "rhoj_par", "key"
  )]
}


}
