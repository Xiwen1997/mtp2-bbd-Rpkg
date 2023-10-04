#' @title Bridge-block decomposition approach for precision matrix estimation under nonnegative partial correlations
#'
#' @description Bridge-block decomposition approach for precision matrix estimation under nonnegative partial correlations
#' \preformatted{
#'   minimize     -log det(Theta) + <Theta, S> + sum_{i neq j} Lambda_{ij} |Theta_{ij}|
#'   subject to   Theta > 0, Theta_{ij} <=0 forall i neq j
#' }
#'
#' @author Xiwen Wang, Jiaxi Ying, and Daniel P. Palomar
#'
#' @references X. Wang, J. Ying, and D. P. Palomar, ‘Learning Large-Scale MTP2 Gaussian Graphical Models via Bridge-Block Decomposition,’ accepted in Neural Information Processing Systems (NeurIPS), New Orleans, LA, USA, Dec. 2023.
#'
#' @param S Sample covariacne matrix
#' @param Lambda Regularization coefficient or regularization matrix
#'
#' @return A list containing the following elements:
#' \item{\code{time}}{Run time cost before algorithm stopped.}
#' \item{\code{X_est}}{The estimated precision matrix.}
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
#' bbd_res          <- solver_bbd(S, 0.2)
#'
#' @import igraph
#' @import MASS
#' @export
solver_bbd <- function(S, Lambda) {
  t0 <- Sys.time()
  p <- dim(S)[1]
  if(length(dim(Lambda)) <= 1) {
    Lambda <- Lambda * (matrix(1, nrow=p, ncol=p) - diag(p))
  }

  S_res               <- S - Lambda
  S_res[S_res < 0] <- 0
  diag(S_res)         <- 0
  S_supp              <- as.matrix(S_res > 0)
  G_thresholded       <- graph_from_adjacency_matrix(S_supp, mode = "undirected")

  ## Compute set of bridges in the thresholded graph

  comps <- biconnected_components(G_thresholded)
  bridges <- comps$component_edges[ lapply(comps$component_edges, length) == 1 ]

  ## Removing bridges

  G_reduced <- G_thresholded
  G_reduced <- delete_edges(G_reduced, unlist(bridges))

  ## Get the indexes for sub-graphs containing more than one node

  Final_components <- clusters(G_reduced)$membership
  Final_components_csize <- table(Final_components)

  Subgraphs_index <- which(Final_components_csize > 1)
  Subgraph_list <- vector("list", length(Subgraphs_index))

  kit <- 1
  for (k_index in Subgraphs_index) {
    nodes_indexes <- which(Final_components == k_index)
    Subgraph_list[[kit]] <- nodes_indexes
    kit <- kit + 1
  }

  out_FPN_sub_opt <- vector("list", length(Subgraphs_index))
  for (i in 1:length(Subgraphs_index)) {
    sub_index <- Subgraph_list[[i]]
    S_sub <- S[sub_index, sub_index]
    Lambda_sub <- Lambda[sub_index, sub_index]
    out_FPN_sub_opt[[i]] <- solver_fpn(S_sub, Lambda_sub)
  }
  # Solve with FPN with bridge-block decomposition

  Theta_hat <- matrix(0, p, p)

  Edge_array <- as_edgelist(G_thresholded)

  for (e in 1:nrow(Edge_array)) {
    i <- Edge_array[e, 1]
    j <- Edge_array[e, 2]
    if (Final_components[i] != Final_components[j]) {
      Theta_hat[i, j] <- -(S_res[i, j]) / (S[i, i] * S[j, j] - (S_res[i, j] * S_res[i, j]))
      Theta_hat[j, i] <- Theta_hat[i, j]
    }
  }

  for (i in 1:p) {
    if (Final_components_csize[Final_components[i]] == 1) {
      Theta_ii <- 1
      for (j in neighbors(G_thresholded, i)) {
        Theta_ii <- Theta_ii + (S_res[i, j] * S_res[i, j]) / (S[i, i] * S[j, j] - S_res[i, j] * S_res[i, j])
      }
      Theta_ii <- Theta_ii / S[i, i]
      Theta_hat[i, i] <- Theta_ii
    } else {
      Theta_hat[i, i] <- 1 / S[i, i]
    }
  }

  for (k in as.numeric(which(Final_components_csize > 1))) {
    k_id <- as.numeric(which(Subgraphs_index == k))
    sub_i <- Subgraph_list[[k_id]]
    Theta_sub <- out_FPN_sub_opt[[k_id]][['X_est']]
    Theta_hat[sub_i, sub_i] <- Theta_sub

    for (i in sub_i) {
      for (j in neighbors(G_thresholded, i)) {
        if (Final_components[i] != Final_components[j]) {
          Theta_hat[i, i] <- Theta_hat[i, i] + (1 / S[i, i]) * (S_res[i, j] * S_res[i, j]) / (S[i, i] * S[j, j] - S_res[i, j] * S_res[i, j])
        }
      }
    }
  }
  run_time <- Sys.time() - t0
  out <- list()

  out$time <- run_time
  out$X_est <- Theta_hat
  return(out)
}
