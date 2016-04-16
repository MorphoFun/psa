################### OPPORTUNITY FOR SELECTION #################

#' @title Total opportunity for selection
#'
#' @name I.total
#'
#' @description Calculation of the total opportunity for selection based on the variance in relative fitness: \code{var(W/mean(W))} (Crow 1958, Arnold and Wade 1984, Brodie et al. 1995), where \code{W} = absolute fitness.
#'
#' @usage I.total(fitness, type = c("W", "w"))
#'
#' @param \code{fitness} Vector of numeric or integer values that represent the fitness metric or proxy.
#' @param \code{type} Type of fitness metric or proxy, either absolute fitness (\code{"W"}) or relative fitness (\code{"w"}).
#'
#' @return \code{I.total} returns a single numeric value.
#' @references Arnold SJ, Wade MJ. 1984. On the measurement of natural and sexual selectionL applications. \emph{Evolution} 38(4): 720-734.
#' @references Brodie ED III, Moore AJ, Janzen FJ. 1995. Visualizing and quantifying natural selection. \emph{TREE} 10(8): 313-318.
#' @references Crow JF. 1958. Some possibilities for measuring selection intensities in man. \emph{Human Biology} 30: 1-13. \url{http://www.jstor.org/stable/41449168?seq=1#page_scan_tab_contents}
#' @example
#' # load the dataset
#' data(BumpusMales)
#'

# I.total <- function(fitness, type = c("W", "w")) {
# 	I.tot <- ifelse(type=="w", var(fitness), var(fitness/mean(fitness)))
# 	return(I.tot)
# }

#' # Calculate the total opportunity for selection
#' I.total(BumpusMales$W, type = "W")
