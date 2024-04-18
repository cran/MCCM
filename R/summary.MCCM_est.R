#' Summary a MCCM Estimation Result
#' @import stats
#' @description
#' Display the estimated correlation matrix and std matrix for a MCCM_est list.
#'
#' @param out_MCCM output of function MCCM_est.
#'
#' @return The summary of estimation.
#'
#' @seealso
#' [MCCM_est],
#' [draw_correlation_matrix]
#'
#' @export
#'
summary_MCCM_est <- function(out_MCCM) {
  pair_est = out_MCCM$pair_est
  if(pair_est){
    return(list(
      Rmatrix = round(out_MCCM$Rmatrix,4),
      std_matrix = round(out_MCCM$std_matrix,4)
    ))
  }else{
    return(list(
      Rmatrix = round(out_MCCM$Rmatrix,4),
      std_matrix = round(out_MCCM$std_matrix,4),
      COV = out_MCCM$COV
    ))
  }
}
