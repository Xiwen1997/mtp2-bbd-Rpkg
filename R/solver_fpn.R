#' @title Fast projected Newton-like algorithm for precision matrix estimation under nonnegative partial correlations
#'
#' @description Fast projected Newton-like algorithm for precision matrix estimation under nonnegative partial correlations
#' \preformatted{
#'   minimize     -log det(Theta) + <Theta, S> + sum_{i neq j} Lambda_{ij} |Theta_{ij}|
#'   subject to   Theta > 0, Theta_{ij} <=0 forall i neq j
#' }
#'
#' @author Jiaxi Ying, Xiwen Wang, and Daniel P. Palomar
#'
#' @references J.-F. Cai, J. V. de Miranda Cardoso, D. P. Palomar, and J. Ying, "Fast Projected Newton-like Method for Precision Matrix Estimation under Total Positivity", accepted in Neural Information Processing Systems (NeurIPS), New Orleans, LA, USA, Dec. 2023.
#'
#' @param S Sample covariacne matrix
#' @param lmbd Regularization coefficient or regularization matrix
#' @param opts stopping criterion options, passed as a list. Possible field name values:
#' (1) 'max_iter'  : maximum number of iterations
#' (2) 'max_time'  : maximum time
#' (3) 'tol'       : tolerance, the algorithm stops if ||X_k+1 - X_k||_F / ||X_k||_F < tol
#' (4) 'beta'      : parameter used in backtracking line search
#' (5) 'X_opt'     : exists if the optimal solution is available for computing the objective function error
#' (6) 'edge'      : exists if a disconnectivity set is imposed
#' (7) 'display'   : "1" if display the iteration results, "0" otherwise.
#'
#' @return A list containing the following elements:
#' \item{\code{time}}{Run time cost before algorithm stopped.}
#' \item{\code{X_est}}{The estimated precision matrix.}
#' \item{\code{objective}}{The objective function value when algorithm stopped}
#' \item{\code{obj_itr}}{Store the objective function value for each iteration}
#' \item{\code{time_itr}}{Store the cpu time for each iteration.}
#' \item{\code{iterate}}{The number of iterations cost before algorithm stopped.}
#' \item{\code{relobj_itr}}{Exists if 'X_opt' exists, store the relative error of the objective function value for each iteration, i.e., |fk - fopt|/|fopt|}
#' \item{\code{relerr_iter}}{Exists if 'X_opt' exists, store the relative error of each iteration, i.e., ||X_k - X_opt||_F / ||X_opt||_F}
#' \item{\code{converge}}{ "1": algorithm converges; "0" algorithm does not converge.}
#'
#' @import MASS
#' @examples
#' library(mtp2bbd)
#' library(igraph)
#' p <- 100 # problem dimension
#' BA_graph <- barabasi.game(p,  directed = FALSE)
#' adjacency_matrix <- as_adjacency_matrix(BA_graph,  type = c("both"))
#' max_eig          <- eigen(adjacency_matrix)$values[1]
#' A                <- 1.05*max_eig*diag(p) - adjacency_matrix
#' inv_A            <- solve(A)
#' D                <- diag(sqrt(diag(inv_A)))
#' Mtrue            <- D %*% A %*% D
#' X                <- MASS::mvrnorm(5 * p , mu = rep(0, p), Sigma = solve(Mtrue))
#' S                <- cov(X)
#' fpn_res          <- solver_fpn(S, 0.2)
#'
#' @export
solver_fpn <- function(S, lmbd, opts=NULL) {
  check_stopping_fpn <- function(opts, x, y, iter, t0){
    if(!'max_time' %in% names(opts)){
      opts$max_time <- 1e4
    }
    if(!'max_iter' %in% names(opts)){
      opts$max_iter <- 1e3
    }
    if(!'tol' %in% names(opts)){
      opts$tol <- 1e-6
    }

    stop <- 0
    converge <- 0

    if (Sys.time() - t0 > opts$max_time){
      stop <- 1  # Maximum CPU time exceeded
    }
    if (iter > opts$max_iter){
      stop <- 1  # Maximum iterations exceeded
    }
    if (norm(x - y, "F") / norm(y, "F") < opts$tol){
      stop <- 1  # Condition on successive iterations holds
      converge <- 1
    }

    return(list(stop = stop, converge = converge))
  }

  chol_decomposition <- function(X){
    if(all(eigen(X)$values > 0)){
      L <- chol(X)
      flag_pd <- 0
    } else {
      svd_decomp <- svd(X)
      rank <- sum(svd_decomp$d > 1e-10)
      flag_pd <- rank + 1
      L <- NULL
    }
    return(list(L = L, flag_pd = flag_pd))
  }

  objective_function <- function(T, S){
    tryCatch({
      A <- chol(T)
      FLAG <- 0
      lg_det <- 2 * sum(log(diag(A)))
      fun <- -1 * lg_det + sum(T * S)
      return(list(value = fun, flag = FLAG))
    }, error = function(e) {
      return(list(value = Inf, flag = 1))
    })
  }

  t0 <- Sys.time()

  p <- ncol(S)
  iter <- 0
  delta <- 1e-15
  alpha <- 0.5

  if(length(dim(lmbd)) > 1) {
    Lamb <- lmbd
  } else {
    Lamb <- lmbd * (matrix(1, nrow=p, ncol=p) - diag(p))
  }

  if(is.null(opts)) {
    opts <- list()
  }
  if(!"max_iter" %in% names(opts)) opts$max_iter <- 1e4
  if(!"max_time" %in% names(opts)) opts$max_time <- 1e4
  if(!"tol" %in% names(opts)) opts$tol <- 1e-12
  beta <- ifelse("beta" %in% names(opts), opts$beta, 0.5)
  display <- ifelse("display" %in% names(opts), opts$display, 1)

  flag_edge <- "edge" %in% names(opts)
  edgeset <- which(!as.logical(opts$edge))

  flag_opt <- "X_opt" %in% names(opts)
  if(flag_opt) {
    X_opt <- opts$X_opt
    f_opt <- objective_function(X_opt, S - Lamb)$value
    relobj_iter <- vector()
    relerr_iter <- vector()
  }

  proj <- function(T) {
    pmin(0, T - diag(diag(T))) + diag(diag(T))
  }

  X <- diag(1/diag(S))
  Fcur <- objective_function(X, S - Lamb)
  if(Fcur$flag) {
    X <- diag(1 / (diag(S) + rep(1, p) * 1e-3))
    Fcur <- objective_function(X, S - Lamb)
  }
  objcur <- Fcur$value
  obj_iter <- c(objcur)
  time_iter <- c(Sys.time() - t0)

  if(flag_opt) {
    rel_object <- abs(f_opt - objcur) / abs(f_opt)
    relobj_iter <- c(relobj_iter, rel_object)
    rel_err <- norm(X_opt - X, type="F") / norm(X_opt, type="F")
    relerr_iter <- c(relerr_iter, rel_err)
  }

  grad <- function(T) {
    -solve(T) + S - Lamb
  }
  gradf <- grad(X)

  check <- check_stopping_fpn(opts, X, X + matrix(1e16, nrow=p, ncol=p), iter, t0)

  if(check$stop) {
    X_new <- X
    objnew <- objcur
  }

  while(!check$stop) {
    iter <- iter + 1

    rstset <- (X - diag(diag(X)) - diag(rep(1e3, p)) > -delta) & (gradf < 0)

    if(flag_edge) {
      rstset <- union(rstset, edgeset)
    }

    X_up <- X
    X_up[rstset] <- 0

    grady <- gradf
    grady[rstset] <- 0

    descent <- X_up %*% grady %*% X_up
    descent[rstset] <- 0

    step_size <- 1

    Theta_f <- function(gamma) {
      proj(X_up - gamma * descent)
    }
    X_new <- Theta_f(step_size)

    list_A_flag_pd <- chol_decomposition(X_new)
    A <- list_A_flag_pd[[1]]
    flag_pd <- list_A_flag_pd[[2]]
    if(!flag_pd) {
      lg_det <- 2 * sum(log(diag(A)))
      objnew <- -1 * lg_det + sum(X_new * (S - Lamb))
    } else {
      objnew <- 1e8
    }

    gd <- abs(sum(gradf * descent))
    gdI <- abs(sum(gradf[rstset] * X[rstset]))

    while ((objnew > objcur - alpha * step_size * gd - alpha * gdI) && step_size > .Machine$double.eps) {
      step_size <- step_size * beta
      X_new <- Theta_f(step_size)

      list_A_flag_pd <- chol_decomposition(X_new)
      A <- list_A_flag_pd[[1]]
      flag_pd <- list_A_flag_pd[[2]]
      if(!flag_pd) {
        lg_det <- 2 * sum(log(diag(A)))
        objnew <- -1 * lg_det + sum(X_new * (S - Lamb))
      } else {
        objnew <- 1e8
      }
    }

    check <- check_stopping_fpn(opts, X_new, X, iter, t0)

    obj_iter <- c(obj_iter, objnew)
    time_iter <- c(time_iter, Sys.time() - t0)


    if(flag_opt) {
      rel_object <- abs(f_opt - objnew) / abs(f_opt)
      relobj_iter <- c(relobj_iter, rel_object)
      rel_err <- norm(X_opt - X_new, type="F") / norm(X_opt, type="F")
      relerr_iter <- c(relerr_iter, rel_err)
    }

    if(display && iter %% 5 == 0) {
      cat(paste("iter:", iter, "objective:", sprintf("%.11f", obj_iter[length(obj_iter)]), "cpu time:", sprintf("%.11f", time_iter[length(time_iter)])), "\n")
    }

    gradf <- grad(X_new)
    X <- X_new
    objcur <- objnew
  }

  run_time <- Sys.time() - t0

  out <- list()
  if(flag_opt) {
    out$relobj_itr <- relobj_iter
    out$relerr_itr <- relerr_iter
  }

  out$time <- run_time
  out$X_est <- X_new
  out$objective <- objnew
  out$obj_itr <- obj_iter
  out$time_itr <- time_iter
  out$iterate <- iter
  out$converge <- check$converge

  return(out)
}
