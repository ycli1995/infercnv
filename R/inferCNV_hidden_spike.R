
#' @importFrom Matrix colSums rowMeans
.build_and_add_hspike <- function(
    infercnv_obj,
    sim_method = c('meanvar', 'simple', 'splatter'),
    aggregate_normals = FALSE
) {
  sim_method <- match.arg(sim_method)
  flog.info("Adding h-spike")

  if (has_reference_cells(infercnv_obj)) {
    ## the usual method to use
    normal_cells_idx_lists <- infercnv_obj@reference_grouped_cell_indices
    if (aggregate_normals) {
      ## for experimental use / data exploration
      normal_cells_idx_lists <- list(
        normalsToUse = unlist(infercnv_obj@reference_grouped_cell_indices)
      )
    }
  } else {
    ## the reference-less case:
    idx <- unlist(infercnv_obj@observation_grouped_cell_indices)
    normal_cells_idx_lists <- list(normalsToUse = idx)
    flog.info("-no normals defined, using all observation cells as proxy")
    infercnv_obj_tmp <- infercnv_obj
    infercnv_obj_tmp@observation_grouped_cell_indices <- normal_cells_idx_lists
  }

  params <- list()

  ## build a fake genome with fake chromosomes,
  ## alternate between 'normal' and 'variable' regions.
  num_cells <- 100
  num_genes_per_chr <- 400
  num_total_genes <- nrow(infercnv_obj@expr.data)

  chr_info <- .get_hspike_chr_info(
    num_genes_each = num_genes_per_chr,
    num_total = num_total_genes
  )
  gene_order <- do.call(rbind, lapply(X = chr_info, FUN = function(x) {
    data.frame(
      chr = x$name,
      start = seq_len(x$ngenes),
      stop = seq_len(x$ngenes)
    )
  }))
  num_genes <- nrow(gene_order)
  rownames(gene_order) <- paste0("gene_", seq_len(num_genes))
  genes_means_use_idx <- sample(
    seq_len(nrow(infercnv_obj@expr.data)),
    size = num_genes,
    replace = TRUE
  )

  ## do for each group of normal cells
  sim.counts.matrix <- NULL
  reference_grouped_cell_indices <- observation_grouped_cell_indices <- list()
  cell_counter <- 0
  for (normal_type in names(normal_cells_idx_lists)) {
    flog.info(sprintf("-hspike modeling of %s", normal_type))
    normal_cells_idx <- normal_cells_idx_lists[[normal_type]]

    ## sample gene info from the normal data
    normal_cells_expr <- infercnv_obj@expr.data[, normal_cells_idx, drop = FALSE]
    gene_means_orig <- Matrix::rowMeans(normal_cells_expr)
    gene_means <- gene_means_orig[genes_means_use_idx]

    gene_means[gene_means == 0] <- 1e-3 # just make small nonzero values
    names(gene_means) <- rownames(gene_order)
    mean_p0_table <- NULL

    ## simulate normals:
    if (sim_method == 'splatter') {
      params = .estimateSingleCellParamsSplatterScrape(counts=infercnv_obj@count.data[,normal_cells_idx])
      params[["nGenes"]] <- num_genes
      params[["nCells"]] <- num_cells
      sim.scExpObj = .simulateSingleCellCountsMatrixSplatterScrape(params, use.genes.means=gene_means)

      sim_normal_matrix = counts(sim.scExpObj)
    } else if (sim_method == 'simple') {
      ## using simple
      mean_p0_table <- .get_mean_vs_p0_from_matrix(normal_cells_expr)
      sim_normal_matrix <- .get_simulated_cell_matrix(
        gene_means = gene_means,
        mean_p0_table=mean_p0_table,
        num_cells=num_cells,
        common_dispersion=0.1
      )
    } else if (sim_method == 'meanvar') {
      ## using mean,var trend
      if (exists("infercnv_obj_tmp")) {
        sim_normal_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(
          infercnv_obj = infercnv_obj_tmp,
          gene_means = gene_means,
          num_cells = num_cells,
          include.dropout = TRUE
        )
      } else {
        sim_normal_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(
          infercnv_obj = infercnv_obj,
          gene_means = gene_means,
          num_cells = num_cells,
          include.dropout = TRUE
        )
      }
    }

    spike_norm_name <- sprintf("simnorm_cell_%s", normal_type)
    colnames(sim_normal_matrix) <- paste0(spike_norm_name, seq_len(num_cells))
    rownames(sim_normal_matrix) <- rownames(gene_order)

    ## apply spike-in multiplier vec
    hspike_gene_means <- gene_means
    for (info in chr_info) {
      chr_name <- info$name
      cnv <- info$cnv
      if (cnv != 1) {
        gene_idx <- which(gene_order$chr == chr_name)
        hspike_gene_means[gene_idx] <- hspike_gene_means[gene_idx] * cnv
      }
    }
    if (sim_method == 'splatter') {
      sim_spiked_cnv.scExpObj = .simulateSingleCellCountsMatrixSplatterScrape(
        params = params,
        use.genes.means = hspike_gene_means
      )
      sim_spiked_cnv_matrix <- counts(sim_spiked_cnv.scExpObj)
    } else if (sim_method == 'simple') {
      ## using simple
      sim_spiked_cnv_matrix <- .get_simulated_cell_matrix(
        gene_means = hspike_gene_means,
        mean_p0_table = mean_p0_table,
        num_cells = num_cells,
        common_dispersion = 0.1
      )
    } else if (sim_method == 'meanvar') {
      if (exists("infercnv_obj_tmp")) {
        sim_spiked_cnv_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(
          infercnv_obj = infercnv_obj_tmp,
          gene_means = hspike_gene_means,
          num_cells = num_cells,
          include.dropout = TRUE
        )
      } else {
        sim_spiked_cnv_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(
          infercnv_obj = infercnv_obj,
          gene_means = hspike_gene_means,
          num_cells = num_cells,
          include.dropout = TRUE
        )
      }
    }

    spike_tumor_name <- sprintf("spike_tumor_cell_%s", normal_type)
    colnames(sim_spiked_cnv_matrix) <- paste0(
      spike_tumor_name,
      seq_len(num_cells)
    )
    rownames(sim_spiked_cnv_matrix) <- rownames(gene_order)

    if (is.null(sim.counts.matrix)) {
      sim.counts.matrix <- cbind(sim_normal_matrix, sim_spiked_cnv_matrix)
    } else {
      sim.counts.matrix <- cbind(
        sim.counts.matrix,
        sim_normal_matrix,
        sim_spiked_cnv_matrix
      )
    }

    reference_grouped_cell_indices[[spike_norm_name]] <- (cell_counter + 1):(cell_counter + num_cells)
    cell_counter <- cell_counter + num_cells

    observation_grouped_cell_indices[[spike_tumor_name]] <- (cell_counter + 1):(cell_counter + num_cells)
    cell_counter <- cell_counter + num_cells
  }

  .hspike <- new(
    Class = "infercnv",
    expr.data = sim.counts.matrix,
    count.data = sim.counts.matrix,
    gene_order = gene_order,
    reference_grouped_cell_indices = reference_grouped_cell_indices,
    observation_grouped_cell_indices = observation_grouped_cell_indices
  )

  validate_infercnv_obj(.hspike)

  .hspike <- normalize_counts_by_seq_depth(
    infercnv_obj = .hspike,
    normalize_factor = median(Matrix::colSums(normal_cells_expr))
  ) # make same target counts/cell as normals

  infercnv_obj@.hspike <- .hspike
  return(infercnv_obj)
}


.get_hspike_chr_info <- function(num_genes_each, num_total) {
  num_remaining <- num_total - 10*num_genes_each
  if (num_remaining < num_genes_each) {
    num_remaining <- num_genes_each
  }

  ## design for fake chr
  chr_info <- list(
    list(name = 'chrA', cnv = 1, ngenes = num_genes_each),
    list(name = 'chr_0', cnv = 0.01, ngenes = num_genes_each),
    list(name = 'chr_B', cnv = 1, ngenes = num_genes_each),
    list(name = 'chr_0pt5', cnv = 0.5, ngenes = num_genes_each),
    list(name = 'chr_C', cnv = 1, ngenes = num_genes_each),
    list(name = 'chr_1pt5', cnv = 1.5, ngenes = num_genes_each),
    list(name = 'chr_D', cnv = 1, ngenes = num_genes_each),
    list(name = 'chr_2pt0', cnv = 2.0, ngenes = num_genes_each),
    list(name = 'chr_E', cnv = 1, ngenes = num_genes_each),
    list(name = 'chr_3pt0', cnv = 3, ngenes = num_genes_each),
    list(name = 'chr_F', cnv = 1, ngenes = num_remaining)
  )
  return(chr_info)
}

##
.sim_foreground <- function(infercnv_obj, sim_method) {

  flog.info("## simulating foreground")

  expr.matrix <- infercnv_obj@expr.data

  samples <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)


  normal_data <- expr.matrix[, unlist(infercnv_obj@reference_grouped_cell_indices)]
  target_total_counts <- median(colSums(normal_data))

  params <- NULL
  if (sim_method == 'splatter') {
    params <- .estimateSingleCellParamsSplatterScrape(infercnv_obj@count.data[, unlist(infercnv_obj@reference_grouped_cell_indices)])
  }


  mean_p0_table <- NULL
  if (sim_method == 'simple') {
    mean_p0_table <- .get_mean_vs_p0_from_matrix(normal_data)
  }

  for (sample_name in names(samples)) {

    cell_idx = samples[[ sample_name ]]

    sample_expr_data = expr.matrix[, cell_idx]
    gene_means = rowMeans(sample_expr_data)
    gene_means[gene_means==0] <- 1e-3 #avoid zeros, breaks splatter sim

    num_cells = ncol(sample_expr_data)

    ## sim the tumor matrix
    if (sim_method == 'simple') {

      sim_matrix <- .get_simulated_cell_matrix(gene_means,
                                               mean_p0_table=mean_p0_table,
                                               num_cells=num_cells,
                                               common_dispersion=0.1)
    } else if (sim_method == 'splatter') {
      params[['nCells']] <- num_cells
      sim_matrix <- .simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
      sim_matrix <- counts(sim_matrix)
    } else if (sim_method == 'meanvar') {
      sim_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj,
                                                                   gene_means,
                                                                   num_cells,
                                                                   include.dropout=TRUE)
    } else {
      stop(sprintf("not recognizing --sim_method: %s", args$sim_method))
    }

    expr.matrix[, cell_idx] <- sim_matrix

  }

  infercnv_obj@expr.data <- expr.matrix

  infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj, target_total_counts)

  return(infercnv_obj)
}

