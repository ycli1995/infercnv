.get_simulated_cell_matrix_using_meanvar_trend <- function(
    infercnv_obj,
    gene_means,
    num_cells,
    include.dropout = FALSE
) {
  # should be working on the total sum count normalized data.
  # model the mean variance relationship
  mean_var_table <- .get_mean_var_table(infercnv_obj)

  dropout_logistic_params <- NULL
  if (include.dropout) {
    mean_p0_table <- .get_mean_vs_p0_table(infercnv_obj)
    dropout_logistic_params <- .get_logistic_params(mean_p0_table)
  }
  .get_simulated_cell_matrix_using_meanvar_trend_helper(
    gene_means = gene_means,
    mean_var_table = mean_var_table,
    num_cells = num_cells,
    dropout_logistic_params = dropout_logistic_params
  )
}

.get_simulated_cell_matrix_using_meanvar_trend_helper <- function(
    gene_means,
    mean_var_table,
    num_cells,
    dropout_logistic_params = NULL
) {

  ngenes <- length(gene_means)

  logm <- log1p(mean_var_table$m)
  logv <- log1p(mean_var_table$v)

  mean_var_spline <- smooth.spline(logv ~ logm)

  spike_cell_names <- paste0('sim_cell_', seq_len(num_cells))
  sim_cell_matrix <- matrix(0, nrow = ngenes, ncol = num_cells)
  rownames(sim_cell_matrix) <- names(gene_means)
  colnames(sim_cell_matrix) <- spike_cell_names

  sim_expr_vals <- function(gene_idx) {
    m <- gene_means[gene_idx]
    .sim_expr_val_mean_var_no_dropout(m = m, mean_var_spline = mean_var_spline)
  }
  for (i in seq_len(num_cells)) {
    newvals <- sapply(seq_len(ngenes), FUN = sim_expr_vals)
    sim_cell_matrix[, i] <- newvals
  }
  ## apply dropout
  if (!is.null(dropout_logistic_params)) {
    sim_cell_matrix <- .apply_dropout(
      counts.matrix = sim_cell_matrix,
      dropout_logistic_params = dropout_logistic_params
    )
  }
  return(sim_cell_matrix)
}

.get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix <- function(gene_means, normal_counts_matrix, num_cells, include.dropout=TRUE, cell_groupings=NULL) {

  mean_var_table <- .get_mean_var_given_matrix(normal_counts_matrix, cell_groupings)

  dropout_logistic_params <- NULL
  if (include.dropout) {
    mean_vs_p0_table <- .get_mean_vs_p0_table_from_matrix(normal_counts_matrix, cell_groupings)
    dropout_logistic_params <- .get_logistic_params(mean_vs_p0_table)
  }

  sim_matrix <- .get_simulated_cell_matrix_using_meanvar_trend_helper(gene_means, mean_var_table, num_cells, dropout_logistic_params)

  return(sim_matrix)
}


##' @keywords internal
##' @noRd
##'
.sim_expr_val_mean_var <- function(
    m,
    mean_var_spline,
    dropout_logistic_params
) {
  # include drop-out prediction
  val <- 0
  if (m <= 0) {
    return(val)
  }
  logm <- log1p(m)
  pred_log_var <- predict(mean_var_spline, logm)$y
  var <- max(expm1(pred_log_var), 0)
  val <- round(max(rnorm(n = 1, mean = m, sd = sqrt(var)), 0))

  if ((!is.null(dropout_logistic_params)) & val > 0) {
    dropout_prob <- predict(dropout_logistic_params$spline, log(val))$y[1]
    if (runif(1) <= dropout_prob) {
      ## a drop-out
      val <- 0
    }
  }
  return(val)
}

.sim_expr_val_mean_var_no_dropout <- function(m, mean_var_spline) {
  val <- 0
  if (m <= 0) {
    return(val)
  }
  logm <- log1p(m)
  pred_log_var <- predict(mean_var_spline, logm)$y
  var <- max(expm1(pred_log_var), 0)
  val <- round(max(rnorm(n = 1, mean = m, sd = sqrt(var)), 0))
  return(val)
}

.apply_dropout <- function(counts.matrix, dropout_logistic_params) {
  counts.matrix <- apply(counts.matrix, 1, function(x) {
    mean.val <- mean(x)
    dropout_prob <- predict(dropout_logistic_params$spline, log(mean.val))$y[1]
    nzeros <- sum(x == 0)
    ntotal <- length(x)

    # padj = ( (pzero*total) - (current_nzero) ) / remaining

    padj <- max(0, ( (dropout_prob * ntotal) - (nzeros) ) / ntotal - nzeros)
    flog.debug(sprintf(
      "mean.val: %g, dropout_prob: %g, adj_dropout_prob: %g",
      mean.val,
      dropout_prob,
      padj
    ))
    x.adj <- vapply(x, function(y) {
      if (runif(1) <= padj) {
        return(0)
      }
      return(y)
    }, FUN.VALUE = numeric(1L))
    x.adj
  })
  return(Matrix::t(counts.matrix))
}

##' .get_mean_var_table()
##'
##' Computes the gene mean/variance table based on all defined cell groupings (reference and observations)
##'
##' @param infercnv_obj An infercnv object populated with raw count data
##'
##' @return data.frame with 3 columns: group_name, mean, variance
##'
##'
##' @keywords internal
##' @noRd
##'
.get_mean_var_table <- function(infercnv_obj) {
  group_indices <- c(
    infercnv_obj@observation_grouped_cell_indices,
    infercnv_obj@reference_grouped_cell_indices
  )
  mean_variance_table <- .get_mean_var_given_matrix(
    expr.matrix = infercnv_obj@expr.data,
    cell_cluster_groupings = group_indices
  )
  return(mean_variance_table)
}

#' @importFrom Matrix rowMeans
#' @importFrom MatrixGenerics rowVars
.get_mean_var_given_matrix <- function(
    expr.matrix,
    cell_cluster_groupings = NULL
) {
  if (is.null(cell_cluster_groupings)) {
    ## use all cells
    cell_cluster_groupings <- list(allcells = seq(ncol(expr.matrix)))
  }

  mean_var_table <- NULL
  for (group_name in names(cell_cluster_groupings)) {
    idx <- cell_cluster_groupings[[group_name]]
    expr.data <- expr.matrix[, idx, drop = FALSE]
    m <- Matrix::rowMeans(expr.data)
    v <- MatrixGenerics::rowVars(expr.data)
    if (is.null(mean_var_table)) {
      mean_var_table <- data.frame(g = group_name, m = m, v = v)
    } else {
      mean_var_table <- rbind(
        mean_var_table,
        data.frame(g = group_name, m = m, v = v)
      )
    }
  }
  return(mean_var_table)
}

##' .get_spike_in_average_bounds()
##'
##' return mean bounds for expression of all cells in the spike-in
##'
##' @param infercnv_obj An infercnv object populated with raw count data
##'
##' @return c(left_bound, right_bound)
##'
##' @keywords internal
##' @noRd
##'

