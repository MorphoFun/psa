###### CALCULATING INDIRECT SELECTION - COMPARISON OF METHODS #######
### Outputs a comparison of the selection differentials based on different methods

## From Falconer 1989:
# selection differential = i*sd(P), which i = intesity of selection and sd(P) is the phenotypic standard deviation
# the standardized selection differential would then be S/sd(P)
# See notes about i, though, b/c it depends on normal data

# See Arnold and Wade (1984). On Table 2, they say: "When using Arnold and Wade's (1984) expression (7) to calculate selection differentials, multiply the result by n/(n-1). Expression (7) is for a parameter rather than its statistical estimate."  The equation shown was for both sMean (7a) and sCov (7b). Thus, has the selection differential been miscalculated in the past?

#' @title Comparing estimates of selection differentials
#'
#' @name dCompare
#'
#' @usage dCompare(w, z)
#'
#' @param \code{w} Relative fitness.
#' @param \code{z} Phenotypic trait(s). Character values are not accepted.
#'
#' @section Value: \code{dCov} is the selection differential calculated from the covariance between the relative fitness, \code{w}, and each phenotypic trait \code{z} (Lande and Arnold 1983).
#' @section Value: \code{dMean} is the selection differential calculated as mean(z*) - mean(z), where z* is the phenotypic trait before selection and z is the phenotypic trait after selection (Lande and Arnold 1983).
#' @section Value: \code{dReg} is the selection differential calculated as the partial regression coefficient of relative fitness, \code{w}, against each individual phenotypic trait through univariate linear regressions.
#' @section Value: \code{dCovScale} calculates the selection differential using the equation in \code{dCovScale}, but uses z data that are standardized to a mean of zero and unit variance.
#' @section Value: \code{dMeanScale} calculates the selection differential using the equation in \code{dMeanScale}, but uses z and z* data that are standardized to a mean of zero and unit variance.
#' @section Value: \code{dRegScale} calculates the selection differential using the equation in \code{dRegScale}, but uses z data that are standardized to a mean of zero and unit variance.
#'
#' @return \code{dCompare} returns a matrix of numeric values with 6 rows and X columns, where X = number of phenotpyic traits.
#'
#' @references Lande R, Arnold SJ. 1983. The measurement of selection on correlated characters. \emph{Evolution} 37(6): 1210-1226. \url{http://www.jstor.org/stable/2408842}
#'
#' @examples
#' data(BumpusMales)
#'
#' dCompare(BumpusMales$w, BumpusMales[,3:11])
#' @export


dCompare <- function(w, z, wType = c("gaussian", "binomial")) {
	zScale <- data.frame(scale(z), stringsAsFactors = FALSE)
	df <- cbind(w,z)
	dfScale <- cbind(w, zScale)

		NL <- function(z) {
		z <- z
		z2 <- sapply(z, function(x) x^2)
		ifelse(ncol(z) > 1, colnames(z2) <- paste(names(z), ".Sq"), colnames(z2) <- paste(names(z), ".Sq", sep=""))
		
		if(ncol(z) > 1) {
		# Calculating cross-products for correlational selection
		cp <- array(1, c(ncol(z), ncol(z),nrow(z)))
		for (i in 1:nrow(z)) {cp[,,i] <- t(tcrossprod(as.numeric(z[i,]), as.numeric(z[i,])))}

		## Creating empty matrix with number of individuals as rows, and # trait pair combinations as columns
		# Just using matrix 1 as an example in the first line to determine how many pair combinations there are
		cp.empty.matrix <- matrix(1, nrow=nrow(z), ncol=length(cp[,,1][lower.tri(cp[,,1])]))
		# Pullin the lower triangle of the cp matrix
		for (i in 1:nrow(z)) {cp.empty.matrix[i,] <- cp[,,i][lower.tri(cp[,,i])]}
		cp.paircombos <- data.frame(cp.empty.matrix, stringsAsFactors = F)

		# Generating all non-repeating combinations of trait pairs (e.g., 1vs2, 1vs3, 3vs5, etc.)
		# row 1 = column in data matrix from above
		# row 2 = row in data matrix from above
		PairCombos <- combn(colnames(z[,1:ncol(z)]), 2)

		# Generating names of pair combinations
		PairComboNames <- matrix( "NA", nrow=1, ncol=ncol(PairCombos))
		for (i in 1:ncol(PairCombos)) {PairComboNames[,i] <- paste(PairCombos[1,i], "x", PairCombos[2,i])}

		# Adding those pair combination names as the column names to the data
		colnames(cp.paircombos) <- PairComboNames
		dd <- cbind(z2, cp.paircombos)
		return(dd)
		} else {
      dd <- z2
      return(dd)
      }
		}
		
	zNL <- NL(z)
	zNLScale <- data.frame(scale(NL(z)), stringsAsFactors = FALSE)

	# differential based on Covariance (Price equation)
	dCov <- as.numeric(cov(w, z))
	dCovScale <- as.numeric(cov(w, zScale))
  
	if (wType == "binomial") {
	# differential based on Diffs; this is primarily useful for when you have 'before' and 'after' selection categorizations
	zt <- colMeans(z)
	zs <- colMeans(df[which(df[,1]>0),-1])
	dMean <- zs-zt
	ztScale <- colMeans(zScale)
	zsScale <- colMeans(scale(df[which(df$w>0),-1]))
	dMeanScale <- zsScale-ztScale
	}

	# nonlinear selection differentials need to be corrected for directional selection; double-check that it is s^2 for quadratic selection, and crossproduct of s1*s2 for correlational selection;
	# C = cov(w, (z-zbar)*t(z-zbar))
	# covt <- cov(z)
	# covs <- cov(df[which(df$w>0),-1])

	# differential based on regression
	dReg <- sapply(z, function(x) lm(w ~ x, data = z)$coefficients[2])
	names(dReg) <- names(z)
	dRegScale <- sapply(zScale, function(x) lm(w ~ x, data = zScale)$coefficients[2])
	names(dRegScale) <- names(zScale)

	ifelse(wType == "binomial", dAll <- rbind(dCov, dMean, dReg, dCovScale, dMeanScale, dRegScale), dAll <- rbind(dCov, dReg, dCovScale, dRegScale))
  return(dAll)
	}

