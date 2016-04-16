################### INTENSITY OF SELECTION AKA STANDARDIZED SELECTION DIFFERENTIAL #################

#' @title Intensity of selecion using the selection differential and phenotypic standard deviation AKA standardized selection differentials
#'
#' @name i.AW
#'
#' @description The total intensity of selection represents the maximum amount of response to selection possible for that sample (Crow 1958). According to Arnold and Wade (1984), selection intensities can be calculated as "selection differentials divided by phenotypic standard deviations before selection", which is equivalent to the standardized selection differential. This metric can be applied to quantify how much directional selection shifts the mean of the phenotypic trait, in units of standard deviations, which is "useful for comparing selective pressures in different populations" (Arnold and Wade 1984). However, there is an assumption of normality if one calculates selection intensities in this fashion. According to Falconer (1989),"The intensity of selection, i, depends only on the proportion of the population included in the selected group and, provided the distribution of phenotypic values is normal, it can be determined from tables of the properties of the normal distribution. If p is the proportion selected, i.e., the proportion of the population falling beyond the point of trunction, and z is the height of the ordinate at the point of truncation, then it follows from mathematical properties of the normal distribution that": (selection differential/phenotypic standard deviation) = standardized selection differential = intensity of selection = z/p, where p is the proportion selected from the sample and z is the phenotypic trait.
#'
#' @usage i.AW
#'
#' @references Arnold SJ, Wade MJ. 1984. On the measurement of natural and sexual selection: applications. \emph{Evolution} 38(4): 720-734.
#' @references Crow JF. 1958. Some possibilities for measuring selection intensities in man. \emph{Human Biology} 30: 1-13. \url{http://www.jstor.org/stable/41449168?seq=1#page_scan_tab_contents}
#' @references Falconer DS. 1989. Introduction to Quantitative Genetics. Third Edition. Longman Scientific and Technical. Copublished in the United States with John Wiley and Sons, Inc., New York. p. 192.


i.AW <- function(w, z, scalez = TRUE) {
	if (apply(z, 2, function(x) sd(x)) == 1 && round(colMeans(z), 6)==0) {
		stop("Trait data are already standardized to mean of zero and unit variance. Change to scalez = FALSE")
	}
	if (!apply(z, 2, function(x) sd(x)) == 1 && !round(colMeans(z), 6)==0 && scalez == FALSE) {
		warning("Trait data are not standardized to mean of zero and unit variance. Change to scalez = TRUE")
	}
	if (isTRUE(scalez)) {
		z <- data.frame(scale(z), stringsAsFactors = FALSE)
	}

	dd <- cbind(w,z)
	meanz <- colMeans(z)
	sdz <- sapply(z, function(x) sd(x))

	# calculating selection differentials in different ways
	dCov <- cov(w, z)
	dReg <- apply(z, 2, function(x) lm(w ~ x, data = z)$coefficients[2])

	VCov <- dCov/sdz; row.names(VCov) <- "VCov"
	VReg <- dReg/sdz

	is.dichotomous <- function(object) {
		ff <- factor(object)
		ifelse(length(levels(ff))==2, "TRUE", "FALSE")
	}

	if (is.dichotomous(w)==TRUE) {
		zs <- colMeans(dd[which(dd$w>0),-1])
		dMean <- zs-meanz
		VMean <- dMean/sdz
	}

	if (exists("VMean")) {
		rbind(VCov, VReg, VMean)
	} else {
		rbind(VCov, VReg)
	}

}


#' @references Crow JF. 1958. Some possibilities for measuring selection intensities in man. \emph{Hum. Biol.} 30: 1-13. \url{http://www.jstor.org/stable/41449168?seq=1#page_scan_tab_contents}

I.Crow.m <- function(fitness) {
	is.dichotomous <- function(object) {
		ff <- factor(object)
		ifelse(length(levels(ff))==2, "TRUE", "FALSE")
	}
	if (is.dichotomous(fitness)==FALSE)
		stop(gettextf("Fitness values do not fit criterion of a dichotomous measure of survival and death. Verify that these data are for mortality selection. Crow's (1958) calculation for the intensity of fertility selection is given by I.Crow.f."))
	pd <- length(fitness[which(fitness==0)])/length(fitness)
	ps <- length(fitness[which(fitness>0)])/length(fitness)
	Im <- pd/ps
	return(Im)
}

I.Crow.f <- function() {

}

#' @references Crow JF. 1958. Some possibilities for measuring selection intensities in man. \emph{Human Biology} 30: 1-13. \url{http://www.jstor.org/stable/41449168?seq=1#page_scan_tab_contents}
#' @references Haldane JBS. 1954. The measurement of natural selection. \emph{Caryologia} 6, suppl. 1: 480.


I.Haldane <- function(opt, fitness) {
	is.dichotomous <- function(object) {
		ff <- factor(object)
		ifelse(length(levels(ff))==2, "TRUE", "FALSE")
	}
	if (is.dichotomous(fitness)==TRUE)
		warning("Haldane's equation was originally based on continuous fitness measures; may not be applicable to dichotomous fitness measures, such as survival/death.")
	sbar <- mean(fitness)
}
