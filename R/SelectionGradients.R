########################## SELECTION GRADIENTS ###############################
######## ESTIMATING BETA COEFFICIENTS USING THE LANDE-ARNOLD METHOD ##########
#### Automatically conducts analyses of linear, quadratic, and correlational selection
## w is the fitness measure
## z is the phenotypic trait

## Warning message:
## glm.fit: fitted probabilities numerically 0 or 1 occurred
## arises due to linear separability: https://stat.ethz.ch/pipermail/r-help/2012-March/307352.html
## which can be due to too many predictors and/or small sample size
## This gets output from the nonlinear model (not surprising)
## It still runs, but it means that the glm is forcing separation between the observations
## that have a w=0 and w=1; that could potentially exaggerate differences and
## inflate the important predictors

## Currently spits out warning if the nonlinear model is overparameterized.

## If the data are standardized to unit variance beforehand, standardizing using the var doesn't
# do anything b/c the sd's of all of the traits==1.


#' @title Calculating gradients based on the Lande-Arnold Method
#'
#' @name glam
#' @description \code{glam} is used to fit generalized linear models, specified by error distributions and link functions as denoted by \code{family}, using approaches based on the quantitative framework established by Lande and Arnold (1983). Statistical methods are based on \code{glm} for linear models and \code{glmer} for linear mixed effects models. An option to employ the Janzen and Stern (1998) correction factor for logistic regression models is available via \code{JS = TRUE}. Model formulae are constructed such that regression coefficients and standard errors for quadratic terms do NOT need to be doubled (Stinchcombe et al. 2008)
#'
#' @usage glam(w, z, wType=c("gaussian", "binomial"), JS = FALSE, prep = TRUE, st= NULL, RE = NULL)
#'
#' @param \code{fitness} Fitness measure. Gaussian fitness types should be use relative which, which is calculated as the absolute fitness for each individual \code{W(z)} divided by the mean absolute fitness \code{W}. Binomial fitness types should use the absolute fitness measures (e.g., 0 = failed, 1 = survived)
#' @param \code{z} Phenotypic traits.
#' @param \code{fitType} Type of distribution for relative fitness \code{w}. Option to either "gaussian" or "binomial".
#' @param \code{JS} Janzen and Stern (1998) correction factor for make logistic regression coefficients congruent with linear regression coefficients for estimating multivariate selection. Default is \code{FALSE}. Setting \code{JS = TRUE} applies the Janzen and Stern correction factor to the logistic regression coefficients.
#' @param \code{prep} Option to scale the phenotypic trait data (\code{z}) to a mean of zero and unit variance in preparation for running regression models. Default is set to \code{TRUE}, with the assumption that \code{z} are the raw data for the morphological traits.
#' @param \code{st} Option to standardize the regression coefficients by either the mean or the standard deviation of \code{z}.
#' @param \code{RE} Random effects of the model, if applicable.
#'
#' @details The Lande-Arnold Method is based on the 1983 paper by Russell Lande and Stevan Arnold, entitled "The measurement of selection on correlated characters". Their method involves applying ordinary least-squares (OLS) regression to estimate selection gradients. 
#'
#' @return The function returns an object of classes "\code{glam}", "\code{lm}", and "\code{glm}."
#'
#' @section Value: \code{GL}
#'
#' @section Warning: \strong{These analyses are currently only available for longitudinal data}. Selection gradients for cross-sectional data must be evalated using matrix algebra rather than OLS regressions (Lande and Arnold 1983).
#'
#'
#' @references Janzen FJ, Stern HL. 1998. Logistic regression for empirical studies of multivariate selection. \emph{Evolution} 52(6): 1564-1571. \url{http://www.jstor.org/stable/2411330?seq=1#page_scan_tab_contents}
#' @references Lande R, Arnold SJ. 1983. The measurement of selection on correlated characters. \emph{Evolution} 37(6): 1210-1226. \url{http://www.jstor.org/stable/2408842}
#' @references Stinchcombe JR, Agrawal AF, Hohenlohe PA, Arnold SJ, Blows MW. 2008. Estimating nonlinear selection gradients using quadratic regression coefficients: double or nothing? \emph{Evolution} 62(9): 2435-2440. \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2008.00449.x/abstract}
#'
#' @seealso \code{glamx}, \code{glm}, \code{lm}, \code{summary.blam}, \code{glmer}
#' @examples
#' # use the BumpusMales data set from the tremors package
#' data(BumpusMales)
#' 
#' # Define the input
#' w <- BumpusMales$w
#' z <- BumpusMales[,3:11]
#'
#' # Calculate the selection gradients using glam
#' mod1 <- glam(w, z, "gaussian")
#'
#' # Review the summary statistics
#' summary.glam(mod1)
#' @export

# Need to add an option where the selection gradients are from a LM and the test statistics are from the Log Reg
# Also, consider whether to log-transform data as a default? 

glam <- function(fitness, z, fitType=c("gaussian", "binomial"), JS = FALSE, prep = TRUE, st= NULL, RE = NULL,...) {
	if (!is.null(fitness) && any(fitness < 0))
		stop("negative fitness is not allowed")

	if (fitType=="gaussian" && isTRUE(JS))
		stop("Janzen and Stern (JS) correction only applicable to binomial response types")

	if (round(colMeans(z),6)==0 && sapply(z, FUN=function(x) sd(x))==1 && isTRUE(prep))
		warning("z data seem to be standardized to mean of zero and unit variance already. Consider setting 'prep' to FALSE and re-run glam model.")

	if (!round(colMeans(z),6)==0 && !sapply(z, FUN=function(x) sd(x))==1 && prep==FALSE)
		warning("z data are not standardized to mean of zero and unit variance. Results may not be reliable.")

	# Creating function that will generate data for quadratic and correlational selection using
	# squared terms and cross-products, in addition to linear terms
	prepr <- function(z) {
		z <- z
		z2 <- sapply(z, function(x) x^2)
		colnames(z2) <- paste(names(z), ".Sq", sep="")
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

		dd <- cbind(z, z2, cp.paircombos)
		return(dd)
	}

	all <- prepr(z)
	ifelse(isTRUE(prep), all <- data.frame(scale(all), check.names = FALSE, stringsAsFactors = FALSE) , all <- all)
	ifelse(isTRUE(prep), z <- data.frame(scale(z), check.names = FALSE, stringsAsFactors = FALSE) , z <- z)
	# If prep = TRUE, sds for all traits will equal 1 b/c the data will be scaled tp unit variance
	sds <- sapply(all, function(x) sd(x))

	dL <- data.frame(z, fitness = fitness, stringsAsFactors = FALSE)
	qq <- list()
	if (levels(as.factor(fitness))[1]==0 & levels(as.factor(fitness))[2]==1 & length(levels(as.factor(fitness)))==2 & fitType=="gaussian")
		stop(gettextf("gaussian selected for fitType but w is binomial"))

	# Creating the linear and nonlinear model parameters
	for (i in 1:length(z)) {
		qq[i] <- paste("I(0.5*", names(z)[i], "^2)", sep="")
		QT <- paste(unlist(qq), collapse=" + ", sep=" + ")
		NLM <- paste("fitness ~ ", paste(names(z), collapse=" + "), " + ", QT, " + ", ".:.")
	}
	# Incorporating random effects (RE) into a LMM, if applicable
	if (is.null(RE)==TRUE) {
		GL <- glm(fitness ~ ., data = dL, family = fitType)
		GNL <- glm(as.formula(NLM), data = dL, family = fitType)
	}
	else {
		if (!require("lme4")) {
			install.packages("lme4", dependencies = TRUE)
			library(lme4)
		}
		for (i in 1:length(RE)) {
			REt <- paste(RE, collapse=" + ")
			LM <- paste("fitness ~ ", paste(names(z), collapse=" + "), "+", REt)
			NLM <- paste("fitness ~ ", paste(names(z), collapse=" + "), " + ", QT, " + ", ".:. + ", REt)
			GL <- glmer(as.formula(LM), data=dL, family = fitType, ...)
			GNL <- glmer(as.formula(NLM), data = dL, family = fitType, ...)
		}
	}

	# Standard logistic
	if (fitType=="binomial" && !isTRUE(JS)) {
		GL <- glm(fitness ~ ., data = dL, family = "binomial")
		GNL <- glm(as.formula(NLM), data = dL, family = "binomial")
	}

	# Applying the Janzen and Stern (1998) correction factor for logistic regressions to coefficients
	# Based on the SAS code at: http://www.public.iastate.edu/~fjanzen/software/regress.htm
	if (fitType=="binomial" && isTRUE(JS)) {
		GL <- glm(fitness ~ ., data = dL, family = "binomial")
		GNL <- glm(as.formula(NLM), data = dL, family = "binomial")
		JS.Correct.L <- mean(GL$fitted.values*(1-GL$fitted.values))
		JS.Correct.NL <- mean(GNL$fitted.values*(1-GNL$fitted.values))
		GL$coefficientsorig <- GL$coefficients
		GNL$coefficientsorig <- GNL$coefficients
		GL$coefficients[-1] <- GL$coefficients[-1]*sds[1:length(z)]*JS.Correct.L/mean(fitness)
		GNL$coefficients[-1] <- GNL$coefficients[-1]*JS.Correct.NL*sds*JS.Correct.NL/mean(fitness)
	}

	# Standardizing the beta coefficients
	if (is.null(st)==FALSE && st=="mean") {
		for (i in 2:length(z)) { # starting at 2 b/c 1 is the intercept
			means <- sapply(z, function(x) mean(x))
			GL$coefficients[i] <- GL$coefficients[i]*means[i]
		}
		for (j in 2:length(all)) {
			meansNL <- colMeans(all)
			GNL$coefficients[j] <- GNL$coefficients[j]*meansNL[j]
		}
	}
		if (is.null(st)==FALSE && st=="var") {
			for (i in 2:length(z)) {
				vars <- sapply(z, function(x) var(x))
				GL$coefficients[i] <- GL$coefficients[i]*vars[i]
			}
			for (j in 2:length(all)) {
				varsNL <- sapply(all, function(x) var(x))
				GNL$coefficients[j] <- GNL$coefficients[j]*varsNL[j]
			}
		}

	GL$JS = JS
	GNL$JS = JS
	dd <- list(GL = GL, GNL = GNL)
	#class(dd) <- "glam"
	return(dd)
}

# Could do it so that I query the blam and determine whether it's binomial or continuous;
# if continuous, use summary.glm, but if binomial, need to use summary.blam
# Should also verify that summary.glm foundation is okay to use for calculating parameters (e.g., se)
# for nonlinear models
# Is it confusing that the linear terms in nonlinear Comparisons are from the linear models?

## Need to update/double-check glamx, so that it is congruent with glam and usable with summary.glam


glamx <- function(w, z, method = c("linear", "nonlinear", "both")) {
	# Creating the different lm models
	method = method
	dLx <- cbind(z, w = w)
	n <- length(z)
	idn <- unlist(lapply(1:n, function(x) combn(1:n, x, simplify = F)), recursive = F, use.names = T)
	qq <- list()
	for (i in 1:length(z)) {
		qq[i] <- paste("I(0.5*", names(z)[i], "^2)", sep="")
	}

	# Linear selection
	modL <-lapply(idn, function(x) paste("w ~", paste(names(z)[x], collapse=" + ")))
	# Running the different lm models
	BLx <- lapply(modL, function(x) lm(as.formula(x),data=z))
	names(BLx) <- modL
	# creating table that compares the test statistic output from the models
	BcompareL <- list()
	for (i in 1:length(BLx)) {
		BcompareL[[i]] <- data.frame(
			Type = "Linear",
			Variable = row.names(summary(BLx[[i]])$coefficients),
			Model = names(BLx)[i],
			summary(BLx[[i]])$coefficients,
			Fstat = summary(BLx[[i]])$fstatistic[1],
			numdf = summary(BLx[[i]])$fstatistic[2],
			dendf = summary(BLx[[i]])$fstatistic[3],
			R2 = summary(BLx[[i]])$r.squared,
			aR2 = summary(BLx[[i]])$adj.r.squared,
			stringsAsFactors = FALSE,
			check.names = FALSE,
			row.names = NULL
			)
	}
	TableL <- do.call("rbind", BcompareL)

	# Nonlinear selection
	modx <-lapply(idn, function(x) paste("w ~", paste( paste(names(z)[x], collapse = " + "), paste((qq)[x], collapse=" + "), sep = " + "), sep=" "))
	modNL <- list()
	cc <- list()
	ll <- list()
	mm <- list()
	for (i in 1:length(idn)) {
		if (length(idn[[i]])==1) { cc[i] <- NA; modNL[i] <- modx[i] }
		if (length(idn[[i]])==2) { cc[i] <- lapply(idn[i], function(x) paste(names(z)[x], collapse = ":")); modNL[i] <- paste(modx[i], cc[i], sep=" + ") }
		if (length(idn[[i]])>2) {
			ll[[i]] <- combn(idn[[i]], 2, simplify = F)
			mm[[i]] <- ll[[i]]
			for (j in 1:length(ll[[i]])) {
				mm[[i]][j] <- lapply(ll[[i]][j], function(x) paste(names(z)[x], collapse=":"))
			}
			cc[i] <- lapply(mm[i], function(x) paste(x, collapse = " + "))
			modNL[i] <- paste(modx[i], cc[i], sep=" + ")
		}
	}
	# Running the different lm models
	BNLx <- lapply(modNL, function(x) lm(as.formula(x),data=z))
	names(BNLx) <- modNL # Note names(BNLx) that are too long are abbreviated as ellipses (...)
	# creating table that compares the test statistic output from the models
	BcompareNL <- list()
	for (i in 1:length(BNLx)) {
		BcompareNL[[i]] <- data.frame(
			Type = "Nonlinear",
			Variable = row.names(summary(BNLx[[i]])$coefficients),
			Model = names(BNLx)[i],
			summary(BNLx[[i]])$coefficients,
			Fstat = summary(BNLx[[i]])$fstatistic[1],
			numdf = summary(BNLx[[i]])$fstatistic[2],
			dendf = summary(BNLx[[i]])$fstatistic[3],
			R2 = summary(BNLx[[i]])$r.squared,
			aR2 = summary(BNLx[[i]])$adj.r.squared,
			row.names = NULL,
			check.names= FALSE,
			stringsAsFactors = FALSE
			)
		# Removing linear terms b/c they are correlated with the nonlinear terms
		BcompareNL[[i]] <- BcompareNL[[i]][-c(2:(1+length(idn[[i]]))),]
	}
	TableNL <- do.call("rbind", BcompareNL)

	# Linear output
	if (method == "linear") {
		Boutput <- list(
			Models = BLx,
			Comparisons = TableL
		)
	}

	# Nonlinear output
	if (method == "nonlinear") {
# 		if (!require("lme4")) {
# 			install.packages("lme4", dependencies = TRUE)
# 			library(lme4)
# 		}
		Boutput <- list(
			Models = BNLx,
			Comparisons = TableNL
		)
	}

	# Both linear and nonlinear output
	if (method == "both") {
		Boutput <- list(
			Models = c(BLx,BNLx),
			Comparisons = rbind(TableL, TableNL)
		)
	}
	Boutput
	#class(Boutput) <- c("glam", "lm")
}


### Save the output as a .txt file
saveR <- function (object, file) {
	sink(file)
	print(object)
	sink(NULL)
}



#' @title Summarizing regression models for class 'glam'
#'
#' @name summary.glam
#'
#' @usage summary.glam(object, JS = FALSE, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...)
#' @description \code{summary.glam} provides the summary information for class \code{glam, glm} or \code{lm}. Output is similar to \code{summary}, except that an option to apply the Janzen and Stern (1998) correction for logistic regressions is available in \code{summary.glam}.
#' @param \code{object}
#'
#' @details Ordinary least-squares (OLS) regressions are used to estimate selection gradients, based on methods outlined by Lande and Arnold (1983). Output for logistic regressions have the option of applying the Janzen and Stern (1998) correction factor to make estimates from logistic regressions congruent with those of the OLS method from Lande and Arnold (1983). Quadratic terms have already been coded, so that their regression estimates and standard errors do NOT need to be doubled (Stinchcombe et al. 2008).
#'
#' @references Janzen FJ, Stern HL. 1998. Logistic regression for empirical studies of multivariate selection. \emph{Evolution} 52(6): 1564-1571. \url{http://www.jstor.org/stable/2411330?seq=1#page_scan_tab_contents}
#' @references Lande R, Arnold SJ. 1983. The measurement of selection on correlated characters. \emph{Evolution} 37(6): 1210-1226. \url{http://www.jstor.org/stable/2408842}
#' @references Stinchcombe JR, Agrawal AF, Hohenlohe PA, Arnold SJ, Blows MW. 2008. Estimating nonlinear selection gradients using quadratic regression coefficients: double or nothing? \emph{Evolution} 62(9): 2435-2440. \url{http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2008.00449.x/abstract}
#' @examples
#' 
#' data(BumpusMales)
#' bm <- glam(BumpusMales[,1], BumpusMales[,3:11])
#' summary.glam(bm$GL)
#' @export 


summary.glam <- function (object, JS = FALSE, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE,
	...) {
	if (is.null(object$JS)==FALSE && identical(object$JS, JS)==FALSE)
		stop("mismatch between JS in glm object and summary.blam's optional JS designation")

	if (is.null(object$BNL)==FALSE) {
		prepr <- function(z) {
			z <- z
			z2 <- sapply(z, function(x) x^2)
			colnames(z2) <- paste(names(z), ".Sq", sep="")
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

			dd <- cbind(z, z2, cp.paircombos)
			return(dd)
		}
		all <- prepr(object$data)
		sds <- c(1, sapply(all, function(x) sd(x)))
	}
	else {
		sds <- c(1, sapply(object$data, function(x) sd(x)))
	}


	est.disp <- FALSE
	df.r <- object$df.residual
	if (is.null(dispersion))
		dispersion <- if (object$family$family %in% c("poisson",
			"binomial"))
			1
	else if (df.r > 0) {
		est.disp <- TRUE
		if (any(object$weights == 0))
			warning("observations with zero weight not used for calculating dispersion")
		sum((object$weights * object$residuals^2)[object$weights >
				0])/df.r
	}
	else {
		est.disp <- TRUE
		NaN
	}
	aliased <- is.na(coef(object))
	p <- object$rank
	if (p > 0) {
		p1 <- 1L:p

		# from lm()
		qr.lm <- function(x, ...) {
			if(is.null(r <- x$qr))
				stop("lm object does not have a proper 'qr' component.
					Rank zero or should not have used lm(.., qr=FALSE).")
			r
		}

		Qr <- qr.lm(object)
		coef.p <- object$coefficients[Qr$pivot[p1]]
    #ifelse(isTRUE(JS), coef.p <- object$coefficientsorig[Qr$pivot[p1]], coef.p <- object$coefficients[Qr$pivot[p1]])
    
		covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
		dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
		covmat <- dispersion * covmat.unscaled
		var.cf <- diag(covmat)
		s.err <- sqrt(var.cf)

		#tvalue <- coef.p/s.err
    ifelse(isTRUE(JS), tvalue <- object$coefficientsorig[Qr$pivot[p1]]/s.err, tvalue <- coef.p/s.err)
		JS.Correct.L <- mean(object$fitted.values*(1-object$fitted.values))
		ifelse(JS==FALSE,   s.err <- sqrt(var.cf), s.err <- sqrt(var.cf)*JS.Correct.L)
		
		dn <- c("Estimate", "Std. Error")
		if (!est.disp) {
			pvalue <- 2 * pnorm(-abs(tvalue))
			coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
			dimnames(coef.table) <- list(names(coef.p), c(dn,
				"z value", "Pr(>|z|)"))
		}
		else if (df.r > 0) {
			pvalue <- 2 * pt(-abs(tvalue), df.r)
			coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
			dimnames(coef.table) <- list(names(coef.p), c(dn,
				"t value", "Pr(>|t|)"))
		}
		else {
			coef.table <- cbind(coef.p, NaN, NaN, NaN)
			dimnames(coef.table) <- list(names(coef.p), c(dn,
				"t value", "Pr(>|t|)"))
		}
		df.f <- NCOL(Qr$qr)
	}
	else {
		coef.table <- matrix(, 0L, 4L)
		dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
			"t value", "Pr(>|t|)"))
		covmat.unscaled <- covmat <- matrix(, 0L, 0L)
		df.f <- length(aliased)
	}
	keep <- match(c("call", "terms", "family", "deviance", "aic",
		"contrasts", "df.residual", "null.deviance", "df.null",
		"iter", "na.action"), names(object), 0L)
	ans <- c(object[keep], list(deviance.resid = residuals(object,
		type = "deviance"), coefficients = coef.table, aliased = aliased,
		dispersion = dispersion, df = c(object$rank, df.r, df.f),
		cov.unscaled = covmat.unscaled, cov.scaled = covmat, JS = JS))
	if (correlation && p > 0) {
		dd <- sqrt(diag(covmat.unscaled))
		ans$correlation <- covmat.unscaled/outer(dd, dd)
		ans$symbolic.cor <- symbolic.cor
	}

	print.summary.glam <- function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...)
	{
		cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
			"\n\n", sep = "")
		cat("Deviance Residuals: \n")
		if (x$df.residual > 5) {
			x$deviance.resid <- setNames(quantile(x$deviance.resid,
				na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
		}
		xx <- zapsmall(x$deviance.resid, digits + 1L)
		print.default(xx, digits = digits, na.print = "", print.gap = 2L)
		if (length(x$aliased) == 0L) {
			cat("\nNo Coefficients\n")
		}
		else {
			df <- if ("df" %in% names(x))
				x[["df"]]
			else NULL
			if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
				cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
					sep = "")
			else cat("\nCoefficients:\n")
			coefs <- x$coefficients

			if (!is.null(aliased <- x$aliased) && any(aliased)) {
				cn <- names(aliased)
				coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
					colnames(coefs)))
				coefs[!aliased, ] <- x$coefficients
			}
			printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
				na.print = "NA", ...)
		}
		cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
			format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
				"Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
					"deviance")]), digits = max(5L, digits + 1L)), " on",
				format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
				1L, paste, collapse = " "), sep = "")
		if (nzchar(mess <- naprint(x$na.action)))
			cat("  (", mess, ")\n", sep = "")
		cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),
			"\n\n", "Number of Fisher Scoring iterations: ", x$iter,
			"\n", sep = "")
		if (JS==TRUE)
			cat("\nLogistic regression coefficients were corrected with the Janzen and Stern method\n")
		correl <- x$correlation
		if (!is.null(correl)) {
			p <- NCOL(correl)
			if (p > 1) {
				cat("\nCorrelation of Coefficients:\n")
				if (is.logical(symbolic.cor) && symbolic.cor) {
					print(symnum(correl, abbr.colnames = NULL))
				}
				else {
					correl <- format(round(correl, 2L), nsmall = 2L,
						digits = digits)
					correl[!lower.tri(correl)] <- ""
					print(correl[-1, -p, drop = FALSE], quote = FALSE)
				}
			}
		}
		cat("\n")
		invisible(x)
	}

	output <- print.summary.glam(ans)
	class(output) <- c("summary.glm", "summary.glam")
	print(ans$coefficients[,1])
}



