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



#### Calculating the selection differentials, using the covariance method
#' @title Estimating linear and nonlinear selection differentials
#'
#' @name dCov
#' 
#' @description \code{dCov} estimates the linear and nonlinear selection differentials for one or more phenotypic traits, based on descriptions in Phillips and Arnold (1989). Estimations are based on phenotypic traits that have been standardized to a mean of zero and unit variance. If one phenotypic trait is input, linear and quadratic selection differentials will be output. If more than one phenotypic trait is input, correlational selection differentials will also be output.
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

## From Falconer 1989:
# selection differential = i*sd(P), which i = intesity of selection and sd(P) is the phenotypic standard deviation
# the standardized selection differential would then be S/sd(P)
# See notes about i, though, b/c it depends on normal data

# See Arnold and Wade (1984). On Table 2, they say: "When using Arnold and Wade's (1984) expression (7) to calculate selection differentials, multiply the result by n/(n-1). Expression (7) is for a parameter rather than its statistical estimate."  The equation shown was for both sMean (7a) and sCov (7b). Thus, has the selection differential been miscalculated in the past?

#' @title Estimate phenotypic selection differentials
#'
#' @name differentials
#' 
#' @description \code{differentials} allows the user to evaluate how different methods for calculating phenotypic section differentials influence the output and, thus, the interpretation of indirect selection on phenotypic traits. 
#'
#' @usage differentials(w, z, method = c(1,2,3,4, "all"), standardize = TRUE, ...)
#'
#' @param \code{w} Relative fitness.
#' @param \code{z} Phenotypic trait(s). Character values are not accepted.
#' @param \code{method}: method to estimate the selection differential. 1 = covariance of relative fitness to the trait; 2 = differences in mean, variance, and covariance before and after selection; 3 = matrix algebra approach of phenotypic distributions before and after selection; 4 = ordinary least-squares regression of relative fitness against the trait; "all" = use all of the methods to produce multiple estimates. 
#' @param \code{standardize} Indicate whether phenotypic trait data should be standardized to a mean of zero and unit variance.
#'
#' @section Output:  \code{dCov} is the selection differential calculated from the covariance between the relative fitness, \code{w}, and each phenotypic trait \code{z} (Lande and Arnold 1983).
#' @section Output:  \code{dBeforeAfter} is the selection differential calculated as  as described in Brodie et al. (1995). 
#' @section Output: \code{dMatrix} is the selection differential calculated using matrix algebra (e.g., equation 13b in Lande and Arnold (1983))
#' @section Output:  \code{dReg} is the selection differential calculated as the partial regression coefficient of relative fitness, \code{w}, against each individual phenotypic trait through univariate linear regressions.
#' @section Output:  \code{dCovScale} calculates the selection differential using the equation in \code{dCovScale}, but uses z data that are standardized to a mean of zero and unit variance.

#'
#' @return \code{differentials} returns a data.frame
#'
#' @references Brodie III ED, Moore AJ, Janzen FJ. 1995. Visualizing and quantifying natural selection. \emph{Trends in Ecology and Evolution} 10(8): 313-318. \url{http://www.sciencedirect.com/science/article/pii/S016953470089117X}
#' @references Lande R, Arnold SJ. 1983. The measurement of selection on correlated characters. \emph{Evolution} 37(6): 1210-1226. \url{http://www.jstor.org/stable/2408842}
#'
#' @examples
#' data(BumpusMales)
#'
#' differentials(BumpusMales$w, BumpusMales[,3:11], "all")
#' @export



differentials <- function(w, z, method = c(1,2,3,4, "all"), standardize = TRUE, ...) {
  
  isScale <- function(x) {
    ifelse(class(x) == "data.frame", x <- x, x <- data.frame(x, stringsAsFactors = FALSE))
    if (is.vector(x) == TRUE) {
      ifelse(round(mean(x),6)==0 && sd(x)==1, "TRUE", "FALSE")
    } else {
      ifelse(round(colMeans(x),6)==0 && sapply(x, FUN=function(x) sd(x))==1, "TRUE", "FALSE")
    }
  }
  
  ifelse(isTRUE(standardize) && isScale(z) == FALSE, z <- data.frame(scale(z), stringsAsFactors = FALSE), z <- data.frame(z, stringsAsFactors = FALSE))
  
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
        diffs <- data.frame(t(diffs), stringsAsFactors = FALSE, check.names = FALSE)
        diffs$Method <- "dCov"
        diffsfinal <- diffs[,c(ncol(diffs), 1:(ncol(diffs)-1))]
        return(diffsfinal[-1])
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
        diffs <- data.frame(t(diffs), stringsAsFactors = FALSE, check.names = FALSE)
        diffs$Method <- "dBeforeAfter"
        diffsfinal <- diffs[,c(ncol(diffs), 1:(ncol(diffs)-1))]
        return(diffsfinal[-1])
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
          diffs <- data.frame(t(diffs), stringsAsFactors = FALSE, check.names = FALSE)
          diffs$Method <- "dMatrix"
          diffsfinal <- diffs[,c(ncol(diffs), 1:(ncol(diffs)-1))]
          return(diffsfinal[-1])
      }
      
      dReg <- function(w,z) {
        # double-check whether correlational models should have quadratic terms (maybe not)
        d <- data.frame(w, z)
  
        if (ncol(d) > 2) {
          ll <- list()
          dReg_linearcoeffs <- list()
          qq <- list()
          dReg_quadcoeffs <- list()
          for (i in 1:length(z)) {
            ll[i] <- paste("w ~ ", names(z)[i])
            dReg_linearmod <- lapply(ll, function(x) lm(as.formula(x), data = d))
            dReg_linearcoeffs[[i]] <- dReg_linearmod[[i]]$coefficients[-1]
            dReg_linear <- unlist(dReg_linearcoeffs)
            qq[i] <- paste("w ~ ", names(z)[i], " + I(0.5*", names(z)[i], "^2)", sep="")
            dReg_quadmod <- lapply(qq, function(x) lm(as.formula(x), data = d))
            dReg_quadcoeffs[[i]] <- dReg_quadmod[[i]]$coefficients[-c(1:2)]
            dReg_quad <- unlist(dReg_quadcoeffs)
          }
          idn <- combn(1:length(z), 2)
          modx <- list()
          dReg_corrmod <- list()
          dReg_corrcoeffs <- list()
          for (i in 1:ncol(idn)) {
            modx[[i]] <- paste("w ~ ", paste(names(z)[idn[,i]], collapse = " + "), " + ", paste(names(z)[idn[,i]], collapse = ":"))
            dReg_corrmod <- lapply(modx, function(x) lm(as.formula(x), data = d))
            dReg_corrcoeffs[[i]] <- c(dReg_corrmod[[i]]$coefficients[-c(1:3)])
          }
            dReg_corr <- unlist(dReg_corrcoeffs)
            dReg_nonlinear <- c(dReg_quad, dReg_corr)  
            names(dReg_nonlinear) <- names(multiprod(z))
        } else {
          names(d) <- c("w", "z")
          dReg_linear <- lm(w ~ ., data = d)$coefficients[-1]
          dReg_nonlinearmod <- lm(w ~ z + I(0.5*z^2), data = d)
          dReg_nonlinear <- data.frame(t(dReg_nonlinearmod$coefficients[-c(1:2)]))
          names(dReg_nonlinear) <- "z.Sq"
        }
        diffs <- data.frame(t(c(dReg_linear, dReg_nonlinear)), check.names = FALSE)
        diffs$Method <- "dReg"
        diffsfinal <- diffs[,c(ncol(diffs), 1:(ncol(diffs)-1))]
        return(diffsfinal[-1])
      }
      
      if (method == 1) {
        output <- data.frame((dCov(w,z)), check.names = FALSE, stringsAsFactors = FALSE)
      } else if(method == 2) {
        output <- data.frame((dBeforeAfter(w,z)), check.names = FALSE, stringsAsFactors = FALSE)
      } else if(method == 3) {
        output <- data.frame((dMatrix(w,z)), check.names = FALSE, stringsAsFactors = FALSE)
      } else if(method == 4) {
        output <- data.frame((dReg(w,z)), check.names = FALSE, stringsAsFactors = FALSE)
      } else if(method == "all") {
        output <- data.frame(dCov = t(dCov(w,z)), dBeforeAfter = t(dBeforeAfter(w,z)), dMatrix = t(dMatrix(w,z)), dReg = t(dReg(w,z)), row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
      }
      return(output)
}
  

### BOOTSTRAPPING STANDARD ERRORS AND CONFIDENCE INTERVALS FOR DIFFERENTIALS ###
#' @title Use bootstrapping to estimate standard errors and confidence intervals for selection differentials
#'
#' @name differentials_bootstats
#' 
#' @description \code{differentials_bootstats} allows the user to calculate the standard deviations and confidence intervals for phenotypic selection differientials that are estimated using the \code{differentials} function in \code{psa}.
#'
#' @usage differentials_bootstats(w, z, conf = 0.95, R = 2000, method = c(1,2,3,4))
#'
#' @param \code{w} Relative fitness.
#' @param \code{z} Phenotypic trait(s). Character values are not accepted.
#' @param \code{conf} Confidence interval. 95% confidence interval is set as a default. See ?boot.ci for more details.
#' @param \code{R} Number of bootstrap replicates. 2000 is set as the default. See ?boot for more details.
#' @param \code{method} Method to estimate the selection differential. 1 = covariance of relative fitness to the trait; 2 = differences in mean, variance, and covariance before and after selection; 3 = matrix algebra approach of phenotypic distributions before and after selection; 4 = ordinary least-squares regression of relative fitness against the trait.
#'
#' @section Output: \code{bootoutput} contains the estimates for the phenotypic selection differentials, bias, and standard errors using an "ordinary" resampling method (see the "sim" option in boot::boot for more details).
#' @section Output: \code{se} contains the bootstrapped standard errors.
#' @section Output: \code{ci} contains the confidence intervals for four bootstrapping methods (basic, student, percent, and bca). See boot::boot.ci for more details.
#'
#' @return \code{differentials_bootstats} returns a list of three objects (boot output, standard errors, and confidence intervals).
#'
#' @examples
#' data(BumpusMales)
#'
#' differentials_bootstats(BumpusMales$w, scale(BumpusMales[,3:11]), method = 1)
#' @import boot
#' @export

differentials_bootstats <- function(w, z, conf = 0.95, R = 2000, method = c(1,2,3,4)) {
  df <- data.frame(w, z, stringsAsFactors = FALSE)
    
    differentialsFunc <- function(df, i) {
      d <- df[i,]
      mod <- t(differentials(d[i,1], d[i,-1], method = method))
      return(mod)
    }
    
    boot.out <- boot(data = df, differentialsFunc, R = R)
    
    booted_se <- NULL
    for (i in 1:length(boot.out$t0)) {
      booted_se[i] <- sd(boot.out$t[,i])
    }
    
    trait.names = row.names(boot.out$t0)
    n.traits = length(trait.names)
    ci = numeric(n.traits * 8); dim(ci)<-c(n.traits,8)
    ci.types = c("norm","basic", "perc", "bca")
    for (i in 1:n.traits) {
      y = boot.ci(boot.out, conf = conf, type = ci.types ,index = i)
      ci[i,] = c(y$norm[2:3],y$basic[4:5],y$perc[4:5],y$bca[4:5])
    }
    ci = data.frame(ci)
    rownames(ci) = trait.names
    int = c((1-conf)/2,1-(1-conf)/2)
    v = NULL
    for (i in 1:length(ci.types)) {
      for (j in 1:2) {
        v = c(v,paste(ci.types[i],int[j]))
      }}
    colnames(ci) = v
    output <- list(
      bootoutput = boot.out,
      se = booted_se,
      ci = ci
    )
  return(output)
}
