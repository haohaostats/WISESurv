.wise_new_rows <- function(n_columns) {
  e <- new.env(parent = emptyenv())
  e$n_columns <- n_columns
  e$i <- integer()
  e$j <- integer()
  e$x <- numeric()
  e$rhs <- numeric()
  e$group <- character()
  e$scale <- numeric()
  e
}

.wise_add_row <- function(rows, columns, values, rhs,
                          group = "survival", scale = 1) {
  row <- length(rows$rhs) + 1L
  columns <- as.integer(columns)
  values <- as.numeric(values)
  keep <- values != 0
  rows$i <- c(rows$i, rep.int(row, sum(keep)))
  rows$j <- c(rows$j, columns[keep])
  rows$x <- c(rows$x, values[keep])
  rows$rhs <- c(rows$rhs, as.numeric(rhs))
  rows$group <- c(rows$group, group)
  rows$scale <- c(rows$scale, scale)
  invisible(rows)
}

.wise_rows_matrix <- function(rows) {
  Matrix::sparseMatrix(i = rows$i, j = rows$j, x = rows$x,
                       dims = c(length(rows$rhs), rows$n_columns))
}

.wise_build_problem <- function(grid, cutoff_survival, specification,
                                external = list()) {
  n <- length(grid)
  i0 <- function(j) j
  i1 <- function(j) n + j
  rows <- .wise_new_rows(2L * n)

  for (offset in c(0L, n)) {
    for (j in seq.int(2L, n)) {
      .wise_add_row(rows, c(offset + j, offset + j - 1L), c(1, -1), 0,
                    group = "survival")
    }
  }

  if (specification$non_harm) {
    for (j in seq_len(n)) {
      .wise_add_row(rows, c(i0(j), i1(j)), c(1, -1), 0,
                    group = "non-harm")
    }
  }

  upper <- .wise_eval_parameter(specification$effect_upper, grid,
                                "effect_upper")
  if (any(upper < 0 | upper > 1)) {
    stop("The effect upper envelope must lie in [0, 1].", call. = FALSE)
  }
  for (j in seq_len(n)) {
    .wise_add_row(rows, c(i1(j), i0(j)), c(1, -1), upper[j],
                  group = "envelope")
  }

  .wise_add_row(
    rows,
    c(i1(n), i0(n), i1(1L), i0(1L)),
    c(1, -1, -specification$residual_upper,
      specification$residual_upper),
    0, group = "residual"
  )
  .wise_add_row(
    rows,
    c(i1(n), i0(n), i1(1L), i0(1L)),
    c(-1, 1, specification$residual_lower,
      -specification$residual_lower),
    0, group = "residual"
  )

  widths <- diff(grid)
  for (j in seq.int(2L, n)) {
    columns <- c(i1(j), i0(j), i1(j - 1L), i0(j - 1L))
    .wise_add_row(
      rows, columns, c(1, -1, -1, 1),
      specification$rebound_rate * widths[j - 1L],
      group = "rate", scale = widths[j - 1L]
    )
    .wise_add_row(
      rows, columns, c(-1, 1, 1, -1),
      specification$decline_rate * widths[j - 1L],
      group = "rate", scale = widths[j - 1L]
    )
  }

  external_intervals <- lapply(external, .wise_external_interval, grid = grid)
  for (source in external_intervals) {
    offset <- source$arm * n
    active <- which(is.finite(source$lower) & is.finite(source$upper))
    for (j in active) {
      group <- paste0("external:", source$name)
      .wise_add_row(rows, offset + j, 1, source$upper[j], group = group)
      .wise_add_row(rows, offset + j, -1, -source$lower[j], group = group)
    }
  }

  lower <- rep(0, 2L * n)
  upper_bound <- rep(1, 2L * n)
  lower[i0(1L)] <- upper_bound[i0(1L)] <- cutoff_survival[1L]
  lower[i1(1L)] <- upper_bound[i1(1L)] <- cutoff_survival[2L]

  list(
    A = .wise_rows_matrix(rows),
    rhs = rows$rhs,
    lower = lower,
    upper = upper_bound,
    group = rows$group,
    scale = rows$scale,
    grid = grid,
    external = external_intervals
  )
}
