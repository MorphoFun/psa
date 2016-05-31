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
#' @usage dCompare(w, z)
#'
#' @param \code{w} Relative fitness.
#' @param \code{z} Phenotypic trait(s). Character values are not accepted.
#'
#' @section Output:  \code{dCov} is the selection differential calculated from the covariance between the relative fitness, \code{w}, and each phenotypic trait \code{z} (Lande and Arnold 1983).
#' @section Output:  \code{dMean} is the selection differential calculated as  as described in Brodie et al. (1995). Not output for Gaussian fitness measures.
#' @section Output:  \code{dMean_stdsd} is the selection differential calculated as described in Gvoždík and Smolinský (2015). Not output for Gaussian fitness measures.
#' @section Output:  \code{dReg} is the selection differential calculated as the partial regression coefficient of relative fitness, \code{w}, against each individual phenotypic trait through univariate linear regressions.
#' @section Output:  \code{dCovScale} calculates the selection differential using the equation in \code{dCovScale}, but uses z data that are standardized to a mean of zero and unit variance.
#' @section Output:  \code{dMeanScale} calculates the selection differential using the equation in \code{dMeanScale}, but uses z and z* data that are standardized to a mean of zero and unit variance. Not output for Gaussian fitness measures.
#' @section Output:  \code{dRegScale} calculates the selection differential using the equation in \code{dRegScale}, but uses z data that are standardized to a mean of zero and unit variance.
#'
#' @return \code{dCompare} returns a matrix of numeric values with 7 rows and X columns, where X = number of phenotpyic traits.
#'
#' @references Brodie III ED, Moore AJ, Janzen FJ. 1995. Visualizing and quantifying natural selection. \emph{Trends in Ecology and Evolution} 10(8): 313-318. \url{http://www.sciencedirect.com/science/article/pii/S016953470089117X}
#' @references Gvoždík L, Smolinský R. 2015. Body size, swimming speed, or thermal sensitivity? Predator-imposed selection on amphibian larvae. \emph{BMC Evolutionary Biology} 15: 238. \url{http://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-015-0522-y}
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
	  dCov_linear <- as.numeric(cov(w, z))
	  dCov_linear_scale <- as.numeric(cov(w, zScale))
	  
	  z_dev <- sapply(z, function(x) x - mean(x))
	  dCov_quad <- cov(w, z_dev^2)
	  dCov_quad_scale <- cov(w, (sapply(data.frame(scale(z)), function(x) x - mean(x)))^2)
	  
	  dCov <- c(dCov_linear, dCov_quad)
	  names(dCov) <- c(colnames(z), paste(colnames(z), "^2", sep = ""))
	 
	  dCovScale <- c(dCov_linear_scale, dCov_quad_scale)
	  names(dCovScale) <- c(colnames(zScale), paste(colnames(zScale), "^2", sep = ""))             
	                   
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
	    dMean_stdsd_linear <- (zs-zt)/zt
	    dMean_stdsd_quad <- ((varzs - varzt)/varzt) + dMean_stdsd_linear^2
	    dMean_stdsd_corr <- ((covzs[lower.tri(covzs)] - covzt[lower.tri(covzt)])/covzt[lower.tri(covzt)]) + tcrossprod(dMean_stdsd_linear)[lower.tri(tcrossprod(dMean_stdsd_linear))]
	    
	    dMean_stdsd <- c(dMean_stdsd_linear, dMean_stdsd_quad, dMean_stdsd_corr)
	    names(dMean_stdsd) <- fullNames
	    
	    dMeanScale <- c(dMean_linear_scale, dMean_quad_scale, dMean_corr_scale)
	    names(dMeanScale) <- fullNames
	    
	  } else {
	    zt <- colMeans(z)
	    zs <- mean(df[which(df[,1]>0),-1])
	    dMean_linear <- zs-zt
	    
	    ztScale <- colMeans(zScale)
	    zsScale <- colMeans(zScales)
	    dMean_linear_scale <- zsScale-ztScale
	    
	    varzt <- var(z)
	    varzs <- var(df[which(df[,1]>0),-1])
	    dMean_quad <- varzs - varzt + dMean_linear^2 
	    
	    varzt_scale <- var(zScale)
	    varzs_scale <- var(zScales)
	    dMean_quad_scale <- varzs_scale - varzt_scale + dMean_linear_scale^2 
	    
	    dMean <- c(dMean_linear, dMean_quad)
	    names(dMean) <- c(colnames(z), paste(colnames(z), "^2", sep = ""))
	    
	    dMeanScale <- c(dMean_linear_scale, dMean_quad_scale)
	    names(dMeanScale) <- c(colnames(z), paste(colnames(z), "^2", sep = ""))
	   
	    dMean_stdsd_linear <- (zs-zt)/zt
	    dMean_stdsd_quad <- ((varzs - varzt)/varzt) + dMean_stdsd_linear^2
	    dMean_stdsd <- c(dMean_stdsd_linear, dMean_stdsd_quad)
	    names(dMean_stdsd) <- names(dMeanScale)
	                           
	  }
	}

	# nonlinear selection differentials need to be corrected for directional selection; double-check that it is s^2 for quadratic selection, and crossproduct of s1*s2 for correlational selection;

	# differential based on regression
	if(ncol(z) > 1) {
	  dReg_linear <- sapply(z, function(x) lm(w ~ x, data = z)$coefficients[2])
	  names(dReg_linear) <- names(z)
	  dRegScale_linear <- sapply(zScale, function(x) lm(w ~ x, data = zScale)$coefficients[2])
	  names(dRegScale_linear) <- names(zScale)
	  
	  dReg_quad <- sapply(z, function(x) lm(w ~ x + I(x^2), data = z)$coefficients[3])
	  names(dReg_quad) <- colnames(z2)
	  dRegScale_quad <- sapply(zScale, function(x) lm(w ~ x + I(x^2), data = zScale)$coefficients[3])
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
	  dReg_linear <- sapply(z, function(x) lm(w ~ x, data = z)$coefficients[2])
	  names(dReg_linear) <- names(z)
	  dRegScale_linear <- sapply(zScale, function(x) lm(w ~ x, data = zScale)$coefficients[2])
	  names(dRegScale_linear) <- names(zScale)
	  
	  dReg_quad <- sapply(z, function(x) lm(w ~ x + I(x^2), data = z)$coefficients[3])
	  names(dReg_quad) <- colnames(z2)
	  dRegScale_quad <- sapply(zScale, function(x) lm(w ~ x + I(x^2), data = zScale)$coefficients[3])
	  names(dRegScale_quad) <- colnames(z2)
	  
	  dReg <- c(dReg_linear, dReg_quad)
	  names(dReg) <- c(colnames(z), paste(colnames(z), "^2", sep = ""))
	  
	  dRegScale <- c(dRegScale_linear, dRegScale_quad)
	  names(dRegScale) <- c(colnames(z), paste(colnames(z), "^2", sep = ""))
	}


	ifelse(length(levels(as.factor(w))) == 2, dAll <- rbind(dCov, dMean, dReg, dMean_stdsd, dCovScale, dMeanScale, dRegScale), dAll <- rbind(dCov, dReg, dCovScale, dRegScale))
  return(dAll)
	}


#### Calculating the selection differentials, using the covariance method
#' @title Estimating linear and nonlinear selection differentials
#'
#' @name dCov
#' 
#' @description \code{dCov} estimates the linear and nonlinear selection differentials for one or more phenotypic traits, based on descriptions in Phillips and Arnold (1989). Estimations are based on phenotypic traits that have been normalized to a mean of zero and unit variance. If one phenotypic trait is input, linear and quadratic selection differentials will be output. If more than one phenotypic trait is input, correlational selection differentials will also be output.
#'
#' @usage dCov(w, z)
#'
#' @param \code{w} Relative fitness.
#' @param \code{z} Phenotypic trait(s). Character values are not accepted.
#'
#' @section Output:  \code{dCov} is the selection differential calculated from the covariance between the relative fitness, \code{w}, and each phenotypic trait \code{z}, using the approach described in Table 1 of Phillips and Arnold (1989).
#'
#' @return \code{dCov} returns a vector of numeric values. 
#'
#' @references Phillips PC, Arnold SJ. 1989. Visualizing multivariate selection. \emph{Evolution} 43(6): 1209-1222.  \url{http://www.jstor.org/stable/2409357}
#' @references Lande R, Arnold SJ. 1983. The measurement of selection on correlated characters. \emph{Evolution} 37(6): 1210-1226. \url{http://www.jstor.org/stable/2408842}
#'
#' @examples
#' data(BumpusMales)
#'
#' dCov(BumpusMales$w, BumpusMales[,3:11])
#' @export


# differential based on Covariance (Price equation)

dCov <- function(w, z) {
  
  isScale <- function(x) {
    ifelse(class(x) == "data.frame", x <- x, x <- data.frame(x, stringsAsFactors = FALSE))
    if (is.vector(x) == TRUE) {
      ifelse(round(mean(x),6)==0 && sd(x)==1, "TRUE", "FALSE")
    } else {
      ifelse(round(colMeans(x),6)==0 && sapply(x, FUN=function(x) sd(x))==1, "TRUE", "FALSE")
    }
  }
  
  ifelse(isScale(z) == TRUE, z <- z, z <- data.frame(scale(z), stringsAsFactors = FALSE))
  
  d <- cbind(w, z)  

  if(ncol(d) > 2) {
    dCov_linear <- as.numeric(cov(w, z))
    
    z_dev <- sapply(z, function(x) x - mean(x))
    dCov_nonlinear <- cov(w, multiprod(z_dev))
    
    differentials <- c(dCov_linear, dCov_nonlinear)
    fullNames <- c(colnames(z), colnames(dCov_nonlinear))
    names(differentials) <- fullNames
    
  } else {
    dCov_linear <- as.numeric(cov(w, z))
    
    z_dev <- sapply(z, function(x) x - mean(x))
    dCov_quad <- cov(w, z_dev^2)

    differentials <- c(dCov_linear, dCov_quad)
    names(differentials) <- c("z", "z^2")
  }
  return(differentials)
}


#' @param \code{method}: method to estimate the selection differential. 1 = covariance of relative fitness to the trait; 2 = differences in mean, variance, and covariance before and after selection; 3 = matrix algebra approach of phenotypic distributions before and after selection; 4 = ordinary least-squares regression of relative fitness against the trait; "all" = use all of the methods to produce multiple estimates. 
#'  

differentials <- function(w, z, method = c(1,2,3, "all"), normalize = TRUE, fitType = "gaussian", ...) {
  
  isScale <- function(x) {
    ifelse(class(x) == "data.frame", x <- x, x <- data.frame(x, stringsAsFactors = FALSE))
    if (is.vector(x) == TRUE) {
      ifelse(round(mean(x),6)==0 && sd(x)==1, "TRUE", "FALSE")
    } else {
      ifelse(round(colMeans(x),6)==0 && sapply(x, FUN=function(x) sd(x))==1, "TRUE", "FALSE")
    }
  }
  
  ifelse(isTRUE(normalize) && isScale(z) == FALSE, z <- data.frame(scale(z), stringsAsFactors = FALSE), z <- data.frame(z, stringsAsFactors = FALSE))
  
      ## method 1: Based on covariance equations from Table 1 of Brodie et al. 1995
      dCov <- function(w, z) {
        d <- cbind(w,z)
        dCov_linear <- as.numeric(cov(w, z))
        if (ncol(d) > 2) {
          z_dev <- sapply(z, function(x) x - mean(x))
          dCov_nonlinear <- cov(w, multiprod(z_dev))
          diffs <- c(dCov_linear, dCov_nonlinear)
          fullNames <- c(colnames(z), colnames(dCov_nonlinear))
          names(diffs) <- fullNames
        } else {
          z_dev <- z - colMeans(z)
          dCov_quad <- cov(w, z_dev^2)
          diffs <- c(dCov_linear, dCov_quad)
          names(diffs) <- c("z", "z^2")
        }
        diffs <- c(Method = "dCov", diffs)
        return(diffs)
      }
      
      ## method 2: Based on before-after equations from Table 1 in Brodie et al. 1995
      dBeforeAfter <- function(w,z) { 
        d <- cbind(w,z)
        z_before <- data.frame(z, stringsAsFactors = FALSE) 
        
        if (ncol(d) >2) {
          z_after <- subset(d, w > 0, select = c(names(z)))
          dBA_linear <- sapply(z_after, function(x) mean(x)) - sapply(z, function(x) mean(x))
          dBA_quad <- diag(var(z_after)) - diag(var(z_before)) + dBA_linear^2
          dBA_corr <- cov(z_after)[lower.tri(cov(z_after))] - cov(z_before)[lower.tri(cov(z_before))] + tcrossprod(dBA_linear)[lower.tri(tcrossprod(dBA_linear))]
          diffs <- c(dBA_linear, dBA_quad, dBA_corr)
          fullNames <- c(names(z), names(multiprod(z)))
          names(diffs) <- fullNames
        } else {
          z_after <- data.frame(subset(d, w > 0)[,2], stringsAsFactors = FALSE)
          dBA_linear <- colMeans(z_after) - colMeans(z_before)
          dBA_quad <- var(z_after) - var(z_before) + dBA_linear^2
          diffs <- c(dBA_linear, dBA_quad)
          names(diffs) <- c("z", "z^2")
        }
        diffs <- c(Method = "dBeforeAfter", diffs)
        return(diffs)
      }
      
      ## method 3: matrix algebra approach from Lande and Arnold (1983)
      dMatrix <- function(w,z) {
          d <- cbind(w,z)
          s <- cov(w,z)
          ssT <- as.vector(s) * as.vector(t(s))
          P <- cov(as.matrix(z))
          if (ncol(d) > 2) {
            P_star <- cov(subset(d, w > 0, select = c(names(z)))) 
            C = P_star - P + ssT
            diffs <- c(s, diag(C), C[lower.tri(C, diag = FALSE)])
            fullNames <- c(names(z), names(multiprod(z)))
            names(diffs) <- fullNames
          } else {
            P_star <- cov(as.matrix(subset(d, w > 0)[,-c(1)]))
            C = P_star - P + ssT
            diffs <- c(s,C)
            names(diffs) <- c("z", "z^2")
          }
          diffs <- c(Method = "dMatrix", diffs)
          return(diffs)
      }
      
      dReg <- function(w,z) {
        # double-check whether correlational models should have quadratic terms (maybe not)
        d <- cbind(w, z)
        fitType == fitType
        dReg_linear <- lm(w ~ ., data = d)
        if (ncol(d) > 2) {
          dReg_nonlinearmod <- glam(w, z, fitType = fitType, prep = FALSE)
          dReg_nonlinear <- dReg_nonlinearmod$GNL$coefficients[-c(1:(length(z)+1))]
        } else {
          names(d) <- c("w", "z")
          dReg_nonlinearmod <- lm(w ~ z + I(0.5*z^2), data = d)
          dReg_nonlinear <- dReg_nonlinearmod$coefficients[-c(1:2)]
        }
        diffs <- c("dReg", dReg_linear$coefficients[-1], dReg_nonlinear)
        names(diffs)[1] <- "Method"
        return(diffs)
      }
      
      if (method == 1) {
        output <- data.frame(dCov(w,z), check.names = FALSE)
      } else if(method == 2) {
        output <- data.frame(dBeforeAfter(w,z), check.names = FALSE)
      } else if(method == 3) {
        output <- data.frame(dMatrix(w,z), check.names = FALSE)
      } else if(method == 4) {
        output <- data.frame(dReg(w,z), check.names = FALSE)
      } else if(method == "all") {
        output <- data.frame(dCov = dCov(w,z), dBeforeAfter = dBeforeAfter(w,z), dMatrix = dMatrix(w,z), dReg = dReg(w,z))
      }
      return(output)
}
  
