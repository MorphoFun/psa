################### OPPORTUNITY FOR SELECTION #################

#' @title Total opportunity for selection
#'
#' @name I_total
#'
#' @description Calculation of the total opportunity for selection based on the variance in relative fitness: \code{var(W/mean(W))} (Crow 1958, Arnold and Wade 1984, Brodie et al. 1995), where \code{W} = absolute fitness.
#'
#' @usage I_total(fitness, type = c("W", "w"))
#'
#' @param \code{fitness} Vector of numeric or integer values that represent the fitness metric or proxy.
#' @param \code{type} Type of fitness metric or proxy, either absolute fitness (\code{"W"}) or relative fitness (\code{"w"}).
#'
#' @return \code{I_total} returns a single numeric value.
#' @references Arnold SJ, Wade MJ. 1984. On the measurement of natural and sexual selection applications. \emph{Evolution} 38(4): 720-734. \url{http://www.jstor.org/stable/2408384?seq=1#page_scan_tab_contents}
#' @references Brodie ED III, Moore AJ, Janzen FJ. 1995. Visualizing and quantifying natural selection. \emph{TREE} 10(8): 313-318. \url{http://www.sciencedirect.com/science/article/pii/S016953470089117X}
#' @references Crow JF. 1958. Some possibilities for measuring selection intensities in man. \emph{Human Biology} 30: 1-13. \url{http://www.jstor.org/stable/41449168?seq=1#page_scan_tab_contents}
#'
#' @export
I_total <- function(fitness, type = c("W", "w")) {
	ifelse(type=="w", w <- fitness, w <- fitness/mean(fitness))
	var(w)
}

#' @examples
#' # load the dataset
#' data(BumpusMales)

#' # Calculate the total opportunity for selection
#' I_total(BumpusMales$W, type = "W")



#' @title Partitioning the opportunity for selection due to individual phenotypic traits 
#'
#' @name I_traits
#'
#' @description Calculation of the opportunity for selection based on the method provided by Moorad and Wade (2013), in which the variation in fitness due to a given trait can be calculated as the product of its selection gradient and selection differential. The I_traits function was written using the code provided in the supplementary material of Moorad and Wade (2013), with slight modifications including a slight correction to  s = b * q * unimputed.var rather than s = B * q * unimputed.var. In addition, the original method by Moorad and Wade (2013) estimated the selection differentials using simple linear regressions that could incorporate weights, but \code{I_traits} estimates the selection differentials using the covariance method (Phillips and Arnold 1983) that does not incorporate weights. The \code{I_traits} function also automatically mean-centers the phenotypic traits. 
#' 
#' @usage I_traits(z, fitness, fitType = c("rel", "abs"))
#'
#' @param \code{z} Phenotypic traits. 
#' @param \code{fitness} Vector of numeric or integer values that represent the fitness metric or proxy.
#' @param \code{fitType} Type of fitness metric or proxy, either absolute fitness (\code{"W"}) or relative fitness (\code{"w"}).
#'
#' @return \code{I_traits} returns a matrix of numeric values. 
#' 
#' \itemize{
#'    \item \code{B:} A vector of selection gradients (normalized to mean of zero and unit variance)
#'    \item \code{b:} slope of a simple linear regression of relative fitness against an individual phenotypic trait.
#'    \item \code{q:} Weights to be applied to the phenotypic traits (defaults to 1).
#'    \item \code{var:} Imputed variance due to missing data, divided by the weight.
#'    \item \code{s:} A vector of selection differentials (normalized to mean of zero and unit variance)
#'    \item \code{i:} Additive opportunity for selection for an individual phenotypic trait. 
#'    \item \code{I.model:} Opportunity of selection, based on the model of phenotypic traits given. 
#'    \item \code{I.total:} The total opportunity for selection, based on weighted variance of fitness.
#' }
#' @references Moorad JA, Wade MJ. 2013. Selection gradients, the opportunity for selection, and the coefficient of determination. \emph{The American Naturalist} 181(3): 291-300. \url{http://www.journals.uchicago.edu/doi/abs/10.1086/669158}
#' @references Phillips PC, Arnold SJ. 1989. Visualizing multivariate selection. \emph{Evolution} 43(6): 1209-1222. \url{http://www.jstor.org/stable/2409357}
#' @import matrixStats
#' @import corpcor
#' @export

I_traits<-function(z, fitness, fitType = c("rel", "abs"), wt = NULL){
  ifelse(fitType == "rel", w <- fitness, w <- fitness/mean(fitness))
  ifelse(is.numeric(wt), wt <- wt, wt <- rep(1, nrow(z)))

  nonlinz <- data.frame(scale(multiprod(z)), stringsAsFactors = FALSE)

  isScale <- function(x) {
    ifelse(class(x) == "data.frame", x <- x, x <- data.frame(x, stringsAsFactors = FALSE))
    if (is.vector(x) == TRUE) {
      ifelse(round(mean(x),6)==0 && sd(x)==1, "TRUE", "FALSE")
    } else {
      ifelse(round(colMeans(x),6)==0 && sapply(x, FUN=function(x) sd(x))==1, "TRUE", "FALSE")
    }
  }

  ifelse(isScale(z) == TRUE, z <- z, z <- data.frame(scale(z), stringsAsFactors = FALSE))

  totalz <- cbind(z, nonlinz)

  q = colWeightedMeans(!is.na(totalz), w=wt)
  means = colWeightedMeans(as.matrix(totalz), w=wt, na.rm = TRUE)
  n = ncol(z)
  for (i in 1:n) {z[is.na(z[,i]),i] = means[i]}
  model.B.lin = lm(w ~ ., data = z, weights = wt)
  qq <- list()
  for (i in 1:length(z)) {
    qq[i] <- paste("I(0.5*", names(z)[i], "^2)", sep="")
    QT <- paste(unlist(qq), collapse=" + ", sep=" + ")
    NLM <- paste("w ~ ", paste(names(z), collapse=" + "), " + ", QT, " + ", ".:.")
  }
  model.B.nonlin <- lm(as.formula(NLM), data = z, weights = wt)

  B = c(model.B.lin$coef[-1], model.B.nonlin$coefficients[-c(1:(n+1))])

  X = data.frame(B=B, p = rep(NA,length(B)), name = names(B))
  p = c(summary(model.B.lin)$coefficients[-1,4], summary(model.B.nonlin)$coefficients[-1,4])
  for (i in 1:length(p)) {X$p[which(X$name == names(p[i]))] = p[i] }

  b <- dCov(w, z)

  imputed.var = NULL

  for (i in 1:ncol(totalz)) {imputed.var = c(imputed.var,wt.moments(totalz[,i], w = wt)$var)}
  unimputed.var = imputed.var/q
  data.frame(B = X$B, b = b, q = q, var = unimputed.var, s = b * q * unimputed.var, i = B * q * unimputed.var * b,
             I.model = sum(B * q * unimputed.var * b, na.rm = TRUE), I.total = wt.moments(x = w,w = wt)$var)}

#' @examples
#'

#' # load the dataset
#' data(BumpusMales)

#' # Calculate the opportunity for selection for the total variation in fitness and variation from phenotypic traits
#' I_traits(BumpusMales[,3:11], BumpusMales$w, fitType = "rel")
