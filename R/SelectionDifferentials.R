###### CALCULATING INDIRECT SELECTION - COMPARISON OF METHODS #######
### Outputs a comparison of the selection differentials based on different methods

#' @title Multiple pair-wise products of data
#' 
#' @name multiprod
#' @usage multiprod(x)
#' 
#' @description Generates the product of all pairwise combinations of variables. xi * xi products are labeled with VariableName.Sq, whereas xi * xj are labeled with VariableNamei x VariableNamej.
#' @param \code{x} A data frame of numeric values. 
#' @return \code{multiprod} returns a data frame of numeric values with X number of columns, where X = number of phenotpyic traits. 
#'
#' @examples
#' data(BumpusMales)
#'
#' multiprod(BumpusMales[,3:11])
#' @export

multiprod <- function(x) {
  ifelse(class(x)== "data.frame", x <- x, x <- data.frame(x, stringsAsFactors = FALSE))
  x2 <- sapply(x, function(x) x^2)
  ifelse(ncol(x) > 1, colnames(x2) <- paste(names(x), ".Sq", sep=""), colnames(x2) <- paste(names(x), ".Sq", sep=""))
  
  if (ncol(x) > 1) {
    cp <- array(1, c(ncol(x), ncol(x), nrow(x)))
    for (i in 1:nrow(x)) {cp[,,i] <- t(tcrossprod(as.numeric(x[i,]), as.numeric(x[i,])))}
    
    ## Creating empty matrix with number of individuals as rows, and # trait pair combinations as columns
    # Just using matrix 1 as an example in the first line to determine how many pair combinations there are
    cp.empty.matrix <- matrix(1, nrow=nrow(x), ncol=length(cp[,,1][lower.tri(cp[,,1])]))
    # Pullin the lower triangle of the cp matrix
    for (i in 1:nrow(x)) {cp.empty.matrix[i,] <- cp[,,i][lower.tri(cp[,,i])]}
    cp.paircombos <- data.frame(cp.empty.matrix, stringsAsFactors = F)
    
    # Generating all non-repeating combinations of trait pairs (e.g., 1vs2, 1vs3, 3vs5, etc.)
    # row 1 = column in data matrix from above
    # row 2 = row in data matrix from above
    PairCombos <- combn(colnames(x[,1:ncol(x)]), 2)
    
    # Generating names of pair combinations
    PairComboNames <- matrix( "NA", nrow=1, ncol=ncol(PairCombos))
    for (i in 1:ncol(PairCombos)) {PairComboNames[,i] <- paste(PairCombos[1,i], " x ", PairCombos[2,i], sep = "")}
    
    # Adding those pair combination names as the column names to the data
    colnames(cp.paircombos) <- PairComboNames
    dd <- cbind(x2, cp.paircombos)
  } else {
    dd <- x2
  }
    return(dd)
}



## From Falconer 1989:
# selection differential = i*sd(P), which i = intesity of selection and sd(P) is the phenotypic standard deviation
# the standardized selection differential would then be S/sd(P)
# See notes about i, though, b/c it depends on normal data

# See Arnold and Wade (1984). On Table 2, they say: "When using Arnold and Wade's (1984) expression (7) to calculate selection differentials, multiply the result by n/(n-1). Expression (7) is for a parameter rather than its statistical estimate."  The equation shown was for both sMean (7a) and sCov (7b). Thus, has the selection differential been miscalculated in the past?

#' @title Comparing estimates of selection differentials
#'
#' @name dCompare
#' 
#' @description \code{dCompare} allows the user to evaluate how different methods for calculating phenotypic section differentials influence the output and, thus, the interpretation of indirect selection on phenotypic traits. 
#'
#' @usage dCompare(w, z, method = c("cov", "mean", "reg", "all"), normalize = TRUE)
#'
#' @param \code{w} Relative fitness.
#' @param \code{z} Phenotypic trait(s). Character values are not accepted.
#' @param \code{method} Method to calculate the selection differentials.
#'
#' @section Output:  \code{dCov} is the selection differential calculated from the covariance between the relative fitness, \code{w}, and each phenotypic trait \code{z} (Lande and Arnold 1983).
#' @section Output:  \code{dMean} is the selection differential calculated as mean(z*) - mean(z), where z* is the phenotypic trait before selection and z is the phenotypic trait after selection (Lande and Arnold 1983). Not output for Gaussian fitness measures.
#' @section Output:  \code{dReg} is the selection differential calculated as the partial regression coefficient of relative fitness, \code{w}, against each individual phenotypic trait through univariate linear regressions.
#' @section Output:  \code{dCovScale} calculates the selection differential using the equation in \code{dCovScale}, but uses z data that are standardized to a mean of zero and unit variance.
#' @section Output:  \code{dMeanScale} calculates the selection differential using the equation in \code{dMeanScale}, but uses z and z* data that are standardized to a mean of zero and unit variance. Not output for Gaussian fitness measures.
#' @section Output:  \code{dRegScale} calculates the selection differential using the equation in \code{dRegScale}, but uses z data that are standardized to a mean of zero and unit variance.
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

dCompare <- function(w, z) {
  z <- data.frame(z, stringsAsFactors = FALSE)
  z2 <- sapply(z, function(x) x^2)
  ifelse(ncol(z) > 1, colnames(z2) <- paste(names(z), ".Sq", sep=""), colnames(z2) <- paste(names(z), ".Sq", sep=""))
	df <- cbind(w,z)
	
	zScale <- data.frame(scale(z), stringsAsFactors = FALSE)
	zScales <- data.frame(scale(df[which(df[,1]>0),-1]), stringsAsFactors = FALSE)
	
	dfScale <- cbind(w, zScale)

		if(ncol(z) > 1) {
		  sds <- c(sapply(z, function(x) sd(x)), sapply(data.frame(multiprod(z)), function(x) sd(x)))  
		  dd <- multiprod(z)
		  dd_scale <- scale(multiprod(z))
		  sds_scale <- c(sapply(zScale, function(x) sd(x)), sapply(data.frame(dd_scale), function(x) sd(x))) 
		} else {
      dd <- z2
      }
		

	# differential based on Covariance (Price equation)
	if(ncol(z) > 1) {
	  dCov_linear <- as.numeric(cov(w, z))
	  dCov_linear_scale <- as.numeric(cov(w, zScale))
	  z_dev <- sapply(z, function(x) x - mean(x))
	  dCov_quad <- cov(w, z_dev^2)
	  dCov_quad_scale <- cov(w, (sapply(data.frame(scale(z)), function(x) x - mean(x)))^2)

	  cpdev <- array(1, c(ncol(z_dev), ncol(z_dev),nrow(z_dev)))
	  for (i in 1:nrow(z_dev)) {cpdev[,,i] <- t(tcrossprod(as.numeric(z_dev[i,]), as.numeric(z_dev[i,])))}
	  cpdev.empty.matrix <- matrix(1, nrow=nrow(z_dev), ncol=length(cpdev[,,1][lower.tri(cpdev[,,1])]))
	  # Pullin the lower triangle of the cp matrix
	  for (i in 1:nrow(z_dev)) {cpdev.empty.matrix[i,] <- cpdev[,,i][lower.tri(cpdev[,,i])]}
	  cpdev.paircombos <- data.frame(cpdev.empty.matrix, stringsAsFactors = F)
	  PairCombos <- combn(colnames(z[,1:ncol(z)]), 2)
	  PairComboNames <- matrix( "NA", nrow=1, ncol=ncol(PairCombos))
	  for (i in 1:ncol(PairCombos)) {PairComboNames[,i] <- paste(PairCombos[1,i], "x", PairCombos[2,i])}
	  colnames(cpdev.paircombos) <- PairComboNames
	  dCov_corr <- cov(w, cpdev.paircombos)
	  dCov_corr_scale <- cov(w, data.frame(scale(cpdev.paircombos)))
	  dCov <- c(dCov_linear, dCov_quad, dCov_corr)
	  dCovScale <- c(dCov_linear_scale, dCov_quad_scale, dCov_corr_scale)
	  fullNames <- c(colnames(z), colnames(z2), PairComboNames)
	  names(dCov) <- fullNames; names(dCovScale) <- fullNames
  } else {
	  dCov <- as.numeric(cov(w, z))
	  dCovScale <- as.numeric(cov(w, zScale))
	}

	  
	if (length(levels(as.factor(w))) == 2) {
	# differential based on taking differences between 'before' and 'after' selection events
	  if(ncol(z) > 1) {
	    zt <- colMeans(z)
	    zs <- colMeans(df[which(df[,1]>0),-1])
	    dMean_linear <- zs-zt
	    
	    ztScale <- colMeans(zScale)
	    zsScale <- colMeans(zScales)
	    dMean_linear_scale <- zsScale-ztScale
	    
	    varzt <- sapply(z, function(x) var(x))
	    varzs <- sapply(df[which(df[,1]>0),-1], function(x) var(x))
	    dMean_quad <- varzs - varzt + dMean_linear^2 
	    
	    varzt_scale <- sapply(zScale, function(x) var(x))
	    varzs_scale <- sapply(zScales, function(x) var(x))
	    dMean_quad_scale <- varzs_scale - varzt_scale + dMean_linear_scale^2 
	    
	    covzt <- cov(z)
	    covzs <- cov(df[which(df[,1]>0),-1])
	    dMean_corr <- covzs[lower.tri(covzs)] - covzt[lower.tri(covzt)] + tcrossprod(dMean_linear)[lower.tri(tcrossprod(dMean_linear))]
	    
	    covzt_scale <- cov(zScale)
	    covzs_scale <- cov(zScales)
	    dMean_corr_scale <- covzs_scale[lower.tri(covzs_scale)] - covzt_scale[lower.tri(covzt_scale)] + tcrossprod(dMean_linear_scale)[lower.tri(tcrossprod(dMean_linear_scale))]
	    
	    dMean <- c(dMean_linear, dMean_quad, dMean_corr)
	    names(dMean) <- fullNames
	    dMean_stdsd <- dMean/sds
	    dMean_stdmean <- dMean/(c(colMeans(z), colMeans(multiprod(z))))
	    
	    dMeanScale <- c(dMean_linear_scale, dMean_quad_scale, dMean_corr_scale)
	    names(dMeanScale) <- fullNames
	    
	  } else {
	    zt <- colMeans(z)
	    zs <- mean(df[which(df[,1]>0),-1])
	    dMean <- zs-zt
	    
	    dMean_stdsd <- dMean/sapply(z, function(x) sd(x))
	    
	    dMean_stdmean <- dMean/sapply(z, function(x) mean(x))
	    
	    ztScale <- colMeans(zScale)
	    zsScale <- colMeans(scale(df[which(df$w>0),-1]))
	    dMeanScale <- zsScale-ztScale
	  }
	}

	# nonlinear selection differentials need to be corrected for directional selection; double-check that it is s^2 for quadratic selection, and crossproduct of s1*s2 for correlational selection;

	# differential based on regression
	if(ncol(z) > 1) {
	  dReg_linear <- sapply(z, function(x) lm(w ~ x, data = z)$coefficients[2])
	  names(dReg_linear) <- names(z)
	  dRegScale_linear <- sapply(zScale, function(x) lm(w ~ x, data = zScale)$coefficients[2])
	  names(dRegScale_linear) <- names(zScale)
	  
	  dReg_quad <- sapply(z, function(x) lm(w ~ x + I(0.5*x^2), data = z)$coefficients[3])
	  names(dReg_quad) <- colnames(z2)
	  dRegScale_quad <- sapply(zScale, function(x) lm(w ~ x + I(0.5*x^2), data = zScale)$coefficients[3])
	  names(dRegScale_quad) <- colnames(z2)
	  
	  PairCombos <- combn(colnames(z[,1:ncol(z)]), 2)
	  corrmod <- list()
	  for (i in 1:ncol(PairCombos)) {
	    corrmod[i] <- paste("w ~", PairCombos[1,i], "*", PairCombos[2,i])
	  }
	  
	  dReg_corrmods <- list()
	  for (j in 1:length(corrmod)) {
	    dReg_corrmods[[j]] <-  lm(as.formula(corrmod[[j]]), data = df) 
	  }
	  dReg_corr <- sapply(dReg_corrmods, function(x) c(x$coefficients[4]))
	  dReg <- c(dReg_linear, dReg_quad, dReg_corr) 
	  names(dReg) <- fullNames
	  
	  dRegScale_corrmods <- list()
	  for (j in 1:length(corrmod)) {
	    dRegScale_corrmods[[j]] <-  lm(as.formula(corrmod[[j]]), data = dfScale) 
	  }
	  dRegScale_corr <- sapply(dRegScale_corrmods, function(x) c(x$coefficients[4]))
	  dRegScale <- c(dRegScale_linear, dRegScale_quad, dRegScale_corr) 
	  names(dRegScale) <- fullNames
	  
	} else {
	  dReg <- sapply(z, function(x) lm(w ~ x, data = z)$coefficients[2])
	  names(dReg) <- names(z)
	  dRegScale <- sapply(zScale, function(x) lm(w ~ x, data = zScale)$coefficients[2])
	  names(dRegScale) <- names(zScale)
	}


	ifelse(length(levels(as.factor(w))) == 2, dAll <- rbind(dCov, dMean, dReg, dMean_stdsd, dMean_stdmean, dCovScale, dMeanScale, dRegScale), dAll <- rbind(dCov, dReg, dCovScale, dRegScale))
  return(dAll)
	}

