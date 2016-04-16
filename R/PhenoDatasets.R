################ Phenotypic selection datasets ################

#' @title Read raw data files that are uploaded on GitHub
#'
#' @name readGitHub
#'
#' @description \code{readGitHub} allows one to read raw data from GitHub repositories
#'
#' @usage readGitHub(url, header = TRUE, fill = TRUE, stringAsFactors = FALSE, ...)
#'
#' @param \code{url} Character string containing the HTML url for the raw data. Supports the same file formats as \code{read.table}.
#' @param \code{header} Logical. Indicate whether header lines should be included. Default is \code{TRUE}.
#' @param \code{fill} Logical. Indicate whether blank fields in rows should be filled with "NA" values to insure that rows have equal length. Default is \code{TRUE}.
#' @param \code{stringsAsFactors}. Logical. Indicate whether character vectors should be automatically converted to factors. Default is \code{FALSE}.
#' @param \code{...} Further arguments that can be passed to the read.table family
#'
#' @details See documentation under ?read.table for specific information regarding data import.
#'
#' @seealso \code{read.table}, \code{connections}
#'
#' @examples
#' # Load the Bumpus.csv data file from GitHub
#' FlowersData <- readGitHub("https://raw.githubusercontent.com/MorphoFun/psa/master/dataraw/Baranzelli_etal_Flowers.csv", sep = ",")
#'
#' # Review the first 6 rows
#' head(FlowersData)
#' @export

readGitHub <- function(url, header = TRUE, fill = TRUE, stringAsFactors = FALSE, ...) {
	store <- tempfile()
	download.file(url, destfile=store, method="curl")
	dd <- read.table(store, header = header, ...)
	return(dd)
}


####### BUMPUS ########

#' @title Bumpus
#'
#' @name Bumpus
#' @description Raw data from the \code{Bumpus} (1899) study on morphological selection in sparrows. Humerus length, femur length, tibiotarsus length, skull width, and sternum keel length were originally in inches, but are converted to millimeters in this version in order to make the units of scale congruent amongst the traits. \code{Bumpus} contains the following variables:
#'	\describe{
#'			\item{\code{w}}{relative fitness, based on survival for all males and females pooled together.}
#'			\item{\code{W}}{absolute fitness, based on survival for all males and females pooled together.}
#'			\item{\code{TotalLength.mm}}{length from the tip of the beak to the tip of the tail in millimeters.}
#'			\item{\code{AlarExtent.mm}}{length of the wing span, from the tip of the left wing to the tip of the right wing in millimeters.}
#'			\item{\code{Weight.g}}{total weight in grams.}
#'			\item{\code{SkullL.mm}}{length of the skull, from the tip of the beak to the occiput in millimeters.}
#'			\item{\code{HumerusL.mm}}{length of the humerus (upper arm bone) in millimeters.}
#'			\item{\code{FemurL.mm}}{length of the femur (thigh bone) in millimeters.}
#'			\item{\code{TibioTarL.mm}}{length of the tibiotarsus (lower leg bone) in millimeters.}
#'			\item{\code{SkullW.mm}}{width of the skull, from the left postorbital bone to the right postorbital bone, in millimeters.}
#'			\item{\code{SternumL.mm}}{length of the keel of the sternum in millimeters.}
#'			\item{\code{Sex}}{sex of the individual (1 = male, 0 = female).}
#'			\item{\code{Age.Group}}{ontogenetic age group (1 = adult, 0 = young, NA = not available).}
#'			}
#' @source A spreadsheet of the morphological data was obtained from the Field Museum: \url{http://www.fieldmuseum.org/science/blog/hermon-bumpus-and-house-sparrows}, with measurements in inches converted to millimeters.
#' @references Bumpus, H.C. 1899. The elimination of the unfit as illustrated by the introduced sparrow, \emph{Passer domesticus}. \emph{Biol. Lectures, Woods Hole Marine Biol. Station}: 209-226.
#' @usage Bumpus
#' @seealso \code{\link{BumpusFemales}, \link{BumpusMales}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(Bumpus)
#'
#' # Look at the structure of the data.frame
#' str(Bumpus)
#'
#' # Run a linear regression with w as the response and the morphological traits as the predictors
#' lm(w ~ ., data=Bumpus[,c(1,3:11)])
#'


####### BUMPUS FEMALES ########

#' @title Bumpus - Females only
#'
#' @name BumpusFemales
#' @description Raw data from the \code{Bumpus} (1899) study on morphological selection in sparrows. Humerus length, femur length, tibiotarsus length, skull width, and sternum keel length were originally in inches, but are converted to millimeters in this version in order to make the units of scale congruent amongst the traits. Bumpus (1899) analyzed the females and males separately, so this subset includes only the data from the females. \code{BumpusFemales} contains the following variables:
#'
#'	\describe{
#'			\item{\code{w}}{relative fitness, based on survival for only females pooled together.}
#'			\item{\code{W}}{absolute fitness, based on survival for only females pooled together.}
#'			\item{\code{TotalLength.mm}}{length from the tip of the beak to the tip of the tail in millimeters.}
#'			\item{\code{AlarExtent.mm}}{length of the wing span, from the tip of the left wing to the tip of the right wing in millimeters.}
#'			\item{\code{Weight.g}}{total weight in grams.}
#'			\item{\code{SkullL.mm}}{length of the skull, from the tip of the beak to the occiput in millimeters.}
#'			\item{\code{HumerusL.mm}}{length of the humerus (upper arm bone) in millimeters.}
#'			\item{\code{FemurL.mm}}{length of the femur (thigh bone) in millimeters.}
#'			\item{\code{TibioTarL.mm}}{length of the tibiotarsus (lower leg bone) in millimeters.}
#'			\item{\code{SkullW.mm}}{width of the skull, from the left postorbital bone to the right postorbital bone, in millimeters.}
#'			\item{\code{SternumL.mm}}{length of the keel of the sternum in millimeters.}
#'			\item{\code{Sex}}{sex of the individual (1 = male, 0 = female).}
#'			\item{\code{Age.Group}}{ontogenetic age group (1 = adult, 0 = young, NA = not available).}
#'			}
#' @source A spreadsheet of the morphological data was obtained from the Field Museum: \url{http://www.fieldmuseum.org/science/blog/hermon-bumpus-and-house-sparrows}, with measurements in inches converted to millimeters.
#' @format data.frame with 49 observations and 13 variables
#' @references Bumpus, H.C. 1899. The elimination of the unfit as illustrated by the introduced sparrow, \emph{Passer domesticus}. \emph{Biol. Lectures, Woods Hole Marine Biol. Station}: 209-226.
#' @usage BumpusFemales
#' @seealso \code{\link{Bumpus}, \link{BumpusMales}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(BumpusFemales)
#'
#' # Look at the structure of the data.frame
#' str(BumpusFemales)
#'
#' # Run a linear regression with w as the response and the morphological traits as the predictors
#' lm(w ~ ., data=BumpusFemales[,c(1,3:11)])




####### BUMPUS MALES ########

#' @title Bumpus - Males only
#'
#' @name BumpusMales
#' @description Raw data from the \code{Bumpus} (1899) study on morphological selection in sparrows. Humerus length, femur length, tibiotarsus length, skull width, and sternum keel length were originally in inches, but are converted to millimeters in this version in order to make the units of scale congruent amongst the traits.  Bumpus (1899) analyzed the females and males separately, so this subset includes only the data from the males. \code{BumpusMales} contains the following variables:
#'
#'	\describe{
#'			\item{\code{w}}{relative fitness, based on survival for only males pooled together.}
#'			\item{\code{W}}{absolute fitness, based on survival for only males pooled together.}
#'			\item{\code{TotalLength.mm}}{length from the tip of the beak to the tip of the tail in millimeters.}
#'			\item{\code{AlarExtent.mm}}{length of the wing span, from the tip of the left wing to the tip of the right wing in millimeters.}
#'			\item{\code{Weight.g}}{total weight in grams.}
#'			\item{\code{SkullL.mm}}{length of the skull, from the tip of the beak to the occiput in millimeters.}
#'			\item{\code{HumerusL.mm}}{length of the humerus (upper arm bone) in millimeters.}
#'			\item{\code{FemurL.mm}}{length of the femur (thigh bone) in millimeters.}
#'			\item{\code{TibioTarL.mm}}{length of the tibiotarsus (lower leg bone) in millimeters.}
#'			\item{\code{SkullW.mm}}{width of the skull, from the left postorbital bone to the right postorbital bone, in millimeters.}
#'			\item{\code{SternumL.mm}}{length of the keel of the sternum in millimeters.}
#'			\item{\code{Sex}}{sex of the individual (1 = male, 0 = female).}
#'			\item{\code{Age.Group}}{ontogenetic age group (1 = adult, 0 = young, NA = not available).}
#'			}
#' @source A spreadsheet of the morphological data was obtained from the Field Museum: \url{http://www.fieldmuseum.org/science/blog/hermon-bumpus-and-house-sparrows}, with measurements in inches converted to millimeters.
#' @format data.frame with 87 observations and 13 variables
#' @references Bumpus, H.C. 1899. The elimination of the unfit as illustrated by the introduced sparrow, \emph{Passer domesticus}. \emph{Biol. Lectures, Woods Hole Marine Biol. Station}: 209-226.
#' @usage BumpusMales
#' @seealso \code{\link{BumpusFemales}, \link{Bumpus}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(BumpusMales)
#'
#' # Look at the structure of the data.frame
#' str(BumpusMales)
#'
#' # Run a linear regression with w as the response and the morphological traits as the predictors
#' lm(w ~ ., data=BumpusMales[,c(1,3:11)])


####### CameroonCichlids ########

#' @title Cameroon crater lake cichlids from the Barombi Mbo lake
#'
#' @name CichlidBarombi
#' @description Raw data from the Martin (2012) study on morphological selection in Cameroon crater lake cichlids.  These data were downloaded from Dryad and then formatted for use in this R package. \code{CichlidBarombi} is a subset of the Martin (2012) data that includes the fish from the Barombi Mbo lake, and contains the following variables:
#'
#'\describe{
#'	\item{\code{w}}{relative fitness, based on relative growth rates averaged from multiple measurements of the spacing between scale circuli of fishes from Barombi.}
#'	\item{\code{W}}{absolute fitness, based on relative growth rates averaged from multiple measurements of the spacing between scale circuli of fishes from Barombi.}
#'	\item{\code{ascending.process.cm}}{distance "from the dorsal tip of the ascending process of the premaxilla to the tip of the most anterior tooth on the dentigerous arm of the premaxilla" (in centimeters).}
#'	\item{\code{body.depth.cm}}{the distance "from the insertion of the first dorsal fin ray to the ventral surface, perpendicular to the major axis of the fish" (in centimeters).}
#'	\item{\code{orbit.diameter.cm}}{average of the diameters of the orbit (eye socket), taken from the horizontal and vertical axes of the fish (in centimeters).}
#'	\item{\code{head.depth.cm}}{distance "from the insertion of the epaxial muscle on the dorsal surface of the neurocranium to the ventral surface, perpendicular to the major axis of the fish" (in centimeters).}
#'	\item{\code{jaw.length.cm}}{distance "from the center of the protruding quadrate-articular joint on the external jaw line to the tip of the most anterior tooth on the mandible" (in centimeters).}
#'	\item{\code{row.ID}}{integer to identify the row number for each observation.}
#'	\item{\code{fishID}}{individual identification number for each fish.}
#'	\item{\code{lake}}{lake that the fish were collected from.}
#'	\item{\code{species.diagnostic.coloration}}{assignment of species based on species-diagnostic-melanin-coloration. Individuals that could not be accurately identified to species were labeled as "ambiguous".}
#'}
#' @source These data were obtained from the "fitness dataset" tab of the "Complete Dataset" Excel spreadsheet that is available on Dryad: \url{http://datadryad.org/resource/doi:10.5061/dryad.rn30d}, and are associated with Martin (2012).
#' @format data.frame with 573 observations and 11 variables
#' @references Martin CH. 2012. Weak disruptive selection and incomplete phenotypic divergence in two classic examples of sympatric speciation: Cameroon crater lake cichlids. \emph{American Naturalist} 180:E90-109. \url{http://dx.doi.org/10.1086/667586}
#' @usage CichlidBarombi
#' @seealso \code{\link{CichlidEjagham}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(CichlidBarombi)
#'
#' # Look at the structure of the data.frame
#' str(CichlidBarombi)
#'
#' # Run a linear regression with w as the response and the morphological traits as the predictors
#' lm(w ~ ., data=CichlidBarombi[,1:5])



#' @title Cameroon crater lake cichlids from the Ejagham lake
#'
#' @name CichlidEjagham
#' @description Raw data from the Martin (2012) study on morphological selection in Cameroon crater lake cichlids.  These data were downloaded from Dryad and then formatted for use in this R package. \code{CichlidEjagham} is a subset of the Martin (2012) data that includes the fish from the Ejagham lake, and contains the following variables:
#'
#'\describe{
#'	\item{\code{w}}{relative fitness, based on relative growth rates averaged from multiple measurements of the spacing between scale circuli of fishes from Ejagham.}
#'	\item{\code{W}}{absolute fitness, based on relative growth rates averaged from multiple measurements of the spacing between scale circuli of fishes from Ejagham.}
#'	\item{\code{ascending.process.cm}}{distance "from the dorsal tip of the ascending process of the premaxilla to the tip of the most anterior tooth on the dentigerous arm of the premaxilla" (in centimeters).}
#'	\item{\code{body.depth.cm}}{the distance "from the insertion of the first dorsal fin ray to the ventral surface, perpendicular to the major axis of the fish" (in centimeters).}
#'	\item{\code{orbit.diameter.cm}}{average of the diameters of the orbit (eye socket), taken from the horizontal and vertical axes of the fish (in centimeters).}
#'	\item{\code{head.depth.cm}}{distance "from the insertion of the epaxial muscle on the dorsal surface of the neurocranium to the ventral surface, perpendicular to the major axis of the fish" (in centimeters).}
#'	\item{\code{jaw.length.cm}}{distance "from the center of the protruding quadrate-articular joint on the external jaw line to the tip of the most anterior tooth on the mandible" (in centimeters).}
#'	\item{\code{row.ID}}{integer to identify the row number for each observation.}
#'	\item{\code{fishID}}{individual identification number for each fish.}
#'	\item{\code{lake}}{lake that the fish were collected from.}
#'	\item{\code{species.diagnostic.coloration}}{assignment of species based on species-diagnostic-melanin-coloration. Individuals that could not be accurately identified to species were labeled as "ambiguous".}
#'}
#' @source These data were obtained from the "fitness dataset" tab of the "Complete Dataset" Excel spreadsheet that is available on Dryad: \url{http://datadryad.org/resource/doi:10.5061/dryad.rn30d}, and are associated with Martin (2012).
#' @format data.frame with 523 observations and 11 variables
#' @references Martin CH. 2012. Weak disruptive selection and incomplete phenotypic divergence in two classic examples of sympatric speciation: Cameroon crater lake cichlids. \emph{American Naturalist} 180:E90-109. \url{http://dx.doi.org/10.1086/667586}
#' @usage CichlidEjagham
#' @seealso \code{\link{CichlidBarombi}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(CichlidEjagham)
#'
#' # Look at the structure of the data.frame
#' str(CichlidEjagham)
#'
#' # Run a linear regression with w as the response and the morphological traits as the predictors
#' lm(w ~ ., data=CichlidEjagham[,1:5])


###### LarvalSquirts #######

#' @title Larvae of sea squirts (\emph{Styela plicata})
#'
#' @name LarvalSquirts
#' @description Raw data from the Crean et al. (2011) study on phenotypic selection across the metamorphic boundary in \emph{Styela plicata} larvae.  These data were downloaded from Dryad and then formatted for use in this R package. \code{LarvalSquirts} includes the data for both the high density and low density treatments, and contains the following variables:
#'
#'\describe{
#'	\item{\code{wr}}{relative fitness based on postmetamorphic survival to reproduction after 4 months in the field (\code{survivalToReprod}) for all larvae.}
#'	\item{\code{ws}}{relative fitness based on postmetamorphic growth, using wet weight (in grams) after 4 months in the field as a proxy (\code{wetWeight4months.g}) for all larvae.}
#'	\item{\code{Wr}}{absolute fitness based on postmetamorphic survival to reproduction after 4 months in the field (\code{survivalToReprod}) for all larvae. 1 = survived, 0 = died.}
#'	\item{\code{Ws}}{absolute fitness based on postmetamorphic growth, using wet weight (in grams) after 4 months in the field as a proxy (\code{wetWeight4months.g}) for all larvae.}
#'	\item{\code{larvalArea.um2}}{area of the larva, "estimated from the average of at least three measurements", in units of squared micrometers.}
#'	\item{\code{hatchTime.mins}}{time to hatching, a premetamorphic trait used as a proxy for development time, in units of minutes.}
#'	\item{\code{settleTime.mins}}{time to settlement, a premetamorphic trait measured every 4 hours, in units of minutes. "Larvae that had not settled within 25 hours posthatch were excluded from their analyses (n = 3)".}
#'	\item{\code{trial}}{one of five replicate runs ("trials") that were run betwen May and October 2009.}
#'	\item{\code{individual}}{individual identification code.}
#'	\item{\code{treatment}}{experimental treatment group (LD = low density, HD = high density).}
#'}
#' @source These data were obtained from the Excel spreadsheet that is available on Dryad: \url{http://datadryad.org/resource/doi:10.5061/dryad.574r7}, and are associated with Crean et al. (2011).
#' @format data.frame with 165 observations and 10 variables
#' @references Crean AJ, Monro K, Marshall DJ. 2011. Fitness consequences of larval traits persist across the metamorphic boundary. \emph{Evolution} 65(11): 3079-3089. \url{doi:10.1111/j.1558-5646.2011.01372.x}
#' @usage LarvalSquirts
#' @seealso \code{\link{LarvalSquirtsHD}, \link{LarvalSquirtsLD}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(LarvalSquirts)
#'
#' # Look at the structure of the data.frame
#' str(LarvalSquirts)
#'
#' # Run a linear regression with wr as the response and the phenotypic traits as the predictors
#' lm(wr ~ ., data = LarvalSquirts[,5:7])



#' @title Larvae of sea squirts (\emph{Styela plicata}) in high density
#'
#' @name LarvalSquirtsHD
#' @description Raw data from the Crean et al. (2011) study on phenotypic selection across the metamorphic boundary in \emph{Styela plicata} larvae.  These data were downloaded from Dryad and then formatted for use in this R package. \code{LarvalSquirtsHD} includes the data from the high density treatment, and consists of the following variables:
#'
#'\describe{
#'	\item{\code{wr}}{relative fitness based on postmetamorphic survival to reproduction after 4 months in the field (\code{survivalToReprod}) for larvae in the high density treatment.}
#'	\item{\code{ws}}{relative fitness based on postmetamorphic growth, using wet weight (in grams) after 4 months in the field as a proxy (\code{wetWeight4months.g}) for larvae in the high density treatment.}
#'	\item{\code{Wr}}{absolute fitness based on postmetamorphic survival to reproduction after 4 months in the field (\code{survivalToReprod}) for larvae in the high density treatment. 1 = survived, 0 = died.}
#'	\item{\code{Ws}}{absolute fitness based on postmetamorphic growth, using wet weight (in grams) after 4 months in the field as a proxy (\code{wetWeight4months.g}) for larvae in the high density treatment.}
#'	\item{\code{larvalArea.um2}}{area of the larva, "estimated from the average of at least three measurements", in units of squared micrometers.}
#'	\item{\code{hatchTime.mins}}{time to hatching, a premetamorphic trait used as a proxy for development time, in units of minutes.}
#'	\item{\code{settleTime.mins}}{time to settlement, a premetamorphic trait measured every 4 hours, in units of minutes. "Larvae that had not settled within 25 hours posthatch were excluded from their analyses (n = 3)".}
#'	\item{\code{trial}}{one of five replicate runs ("trials") that were run betwen May and October 2009.}
#'	\item{\code{individual}}{individual identification code.}
#'	\item{\code{treatment}}{experimental treatment group (LD = low density, HD = high density).}
#'}
#' @source These data were obtained from the Excel spreadsheet that is available on Dryad: \url{http://datadryad.org/resource/doi:10.5061/dryad.574r7}, and are associated with Crean et al. (2011).
#' @format data.frame with 83 observations and 10 variables
#' @references Crean AJ, Monro K, Marshall DJ. 2011. Fitness consequences of larval traits persist across the metamorphic boundary. \emph{Evolution} 65(11): 3079-3089. \url{doi:10.1111/j.1558-5646.2011.01372.x}
#' @usage LarvalSquirtsHD
#' @seealso \code{\link{LarvalSquirts}, \link{LarvalSquirtsLD}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(LarvalSquirtsHD)
#'
#' # Look at the structure of the data.frame
#' str(LarvalSquirtsHD)
#'
#' # Run a linear regression with wr as the response and the phenotypic traits as the predictors
#' lm(wr ~ ., data = LarvalSquirtsHD[,5:7])



#' @title Larvae of sea squirts (\emph{Styela plicata}) in low density
#'
#' @name LarvalSquirtsLD
#' @description Raw data from the Crean et al. (2011) study on phenotypic selection across the metamorphic boundary in \emph{Styela plicata} larvae.  These data were downloaded from Dryad and then formatted for use in this R package. \code{LarvalSquirtsLD} includes the data from the high density treatment, and consists of the following variables:
#'
#'\describe{
#'	\item{\code{wr}}{relative fitness based on postmetamorphic survival to reproduction after 4 months in the field (\code{survivalToReprod}) for larvae in the low density treatment.}
#'	\item{\code{ws}}{relative fitness based on postmetamorphic growth, using wet weight (in grams) after 4 months in the field as a proxy (\code{wetWeight4months.g}) for larvae in the low density treatment.}
#'	\item{\code{Wr}}{absolute fitness based on postmetamorphic survival to reproduction after 4 months in the field (\code{survivalToReprod}) for larvae in the low density treatment. 1 = survived, 0 = died.}
#'	\item{\code{Ws}}{absolute fitness based on postmetamorphic growth, using wet weight (in grams) after 4 months in the field as a proxy (\code{wetWeight4months.g}) for larvae in the low density treatment.}
#'	\item{\code{larvalArea.um2}}{area of the larva, "estimated from the average of at least three measurements", in units of squared micrometers.}
#'	\item{\code{hatchTime.mins}}{time to hatching, a premetamorphic trait used as a proxy for development time, in units of minutes.}
#'	\item{\code{settleTime.mins}}{time to settlement, a premetamorphic trait measured every 4 hours, in units of minutes. "Larvae that had not settled within 25 hours posthatch were excluded from their analyses (n = 3)".}
#'	\item{\code{trial}}{one of five replicate runs ("trials") that were run betwen May and October 2009.}
#'	\item{\code{individual}}{individual identification code.}
#'	\item{\code{treatment}}{experimental treatment group (LD = low density, HD = high density).}
#'}
#' @source These data were obtained from the Excel spreadsheet that is available on Dryad: \url{http://datadryad.org/resource/doi:10.5061/dryad.574r7}, and are associated with Crean et al. (2011).
#' @format data.frame with 82 observations and 10 variables
#' @references Crean AJ, Monro K, Marshall DJ. 2011. Fitness consequences of larval traits persist across the metamorphic boundary. \emph{Evolution} 65(11): 3079-3089. \url{doi:10.1111/j.1558-5646.2011.01372.x}
#' @usage LarvalSquirtsLD
#' @seealso \code{\link{LarvalSquirts}, \link{LarvalSquirtsHD}}
#' @keywords datasets
#' @examples
#' # Load the data
#' data(LarvalSquirtsLD)
#'
#' # Look at the structure of the data.frame
#' str(LarvalSquirtsLD)
#'
#' # Run a linear regression with wr as the response and the phenotypic traits as the predictors
#' lm(wr ~ ., data = LarvalSquirtsLD[,5:7])



###### Galls #######

#' @title Gall size
#'
#' @name Galls
#' @description Raw data from the Egan et al. (2011) study on phenotypic selection on gall size produced by \emph{Belonocnema treatae} on a host tree (\emph{Quercus fusiformis}).  These data were downloaded from Dryad and then formatted for use in this R package. \code{Galls} consists of the following variables:
#'
#'\describe{
#'	\item{\code{w}}{relative fitness, based on the emergence/survivorship of the host-specific gall-former \emph{Belonocnema treatae} (\code{BtreataeEmergence}) on it host plant (the plateau live oak \emph{Quercus fusiformis}).}
#'	\item{\code{W}}{absolute fitness, based on the emergence/survivorship of the host-specific gall-former \emph{Belonocnema treatae} (\code{BtreataeEmergence}) on it host plant (the plateau live oak \emph{Quercus fusiformis}).  1 = emerged/survived, 0 = no emergence.}
#'	\item{\code{gallSizeDiameter.mm}}{diameter of the galls produced by \emph{B. treatae} on its host plant, in units of millimeters and measured to the nearest 0.01 mm.}
#'	\item{\code{year}}{year that the data were collected (either 2004 or 2008)}
#'	\item{\code{tree}}{individual identification of the tree that was sampled. Each tree was treated as a separate population.}
#'	\item{\code{bagPerReplicate}}{identication number of the bag/enclosure that was treated as a blocking factor within each tree.}
#'}
#' @source These data were obtained from the Excel spreadsheet that is available on Dryad: \url{http://datadryad.org/resource/doi:10.5061/dryad.1js1n}, and are associated with Egan et al. (2011).
#' @format data.frame with 14283 observations and 6 variables
#' @references Egan SP, Hood GR, Ott JR. 2011. Natural selection on gall size: variable contributions of individual host plants to population-wide patterns. \emph{Evolution} 65(12): 3543-3557. \url{doi:10.1111/j.1558-5646.2011.01396.x}
#' @usage Galls
#' @keywords datasets
#' @examples
#' # Load the data
#' data(Galls)
#'
#' # Look at the structure of the data.frame
#' str(Galls)
#'
#' # Run a linear regression with wr as the response and the phenotypic traits as the predictors
#' lm(wr ~ ., data = Galls[,3])



###### FlyCHC #######

#' @title Phenotypic selection on \emph{Drosophila serrata} CHC traits
#'
#' @name FlyCHC
#' @description Raw data from the Gosden et al. (2014) study on phenotypic selection on CHC traits in \emph{Drosophila serrata}.  These data were downloaded from Dryad and then formatted for use in this R package. \code{FlyCHC} consists of the following variables:
#'
#' \describe{
#' 		\item{\code{w}}{relative fitness, based on mating success.}
#' 		\item{\code{W}}{absolute fitness, based on mating success (1 = chosen in mating trial, 0 = rejected in mating trial).}
#' 		\item{\code{trait_1}}{logcontrast value of the cuticular hydrocarbon pheromone 1.}
#' 		\item{\code{trait_2}}{logcontrast value of the cuticular hydrocarbon pheromone 2.}
#' 		\item{\code{trait_3}}{logcontrast value of the cuticular hydrocarbon pheromone 3.}
#' 		\item{\code{trait_4}}{logcontrast value of the cuticular hydrocarbon pheromone 4.}
#' 		\item{\code{trait_5}}{logcontrast value of the cuticular hydrocarbon pheromone 5.}
#' 		\item{\code{trait_6}}{logcontrast value of the cuticular hydrocarbon pheromone 6.}
#' 		\item{\code{trait_7}}{logcontrast value of the cuticular hydrocarbon pheromone 7.}
#' 		\item{\code{sex}}{sex of the individual.}
#' }
#' @source These data were obtained from the Excel spreadsheet called "phenotypic_selection" that is available on Dryad: \url{http://datadryad.org/resource/doi:10.5061/dryad.t0g84}, and are associated with Gosden et al. (2014).
#' @format data.frame with 3689 observations and 10 variables
#' @references Gosden TP, Rundle HD, Chenoweth SF. 2014. Testing the correlated response hypothesis for the evolution and maintenance of male mating preferences in \emph{Drosophila serrata. Journal of Evolutionary Biology} 27(10): 2106-2112. \url{doi: 10.1111/jeb.12461}
#' @usage FlyCHC
#' @keywords datasets
#' @examples
#' # Load the data
#' data(FlyCHC)
#'
#' # Look at the structure of the data.frame
#' str(FlyCHC)
#'
#' # Run a linear regerssion with w as the response and the morphological traits as the predictors
#' lm(w ~ ., data = FlyCHC[,3:9])



####### Flowers ########

#' @title Flower morphology Morrenia brachtephana
#'
#' @name Flowers
#' @description Raw data from the Baranzelli et al. (2014) study on morphological selection in the flowers of \emph{Morrenia brachtephana}.  These data were downloaded from Dryad and then formatted for use in this R package. \code{Flowers} consists of the following variables:
#'
#' \describe{
#' 		\item{\code{wm}}{relative male fitness, based on the average number of pollinaria exported per flower from all the flowers sampled for each individual (\code{Pollinarium.export}).}
#' 		\item{\code{wf}}{relative female fitness, based on the average number of pollinia received per flower from the total number collected per individual (\code{pollinia.import}).}
#' 		\item{\code{Wm}}{absolute male fitness, based on the average number of pollinaria exported per flower from all the flowers sampled for each individual (\code{Pollinarium.export}).}
#' 		\item{\code{Wf}}{absolute female fitness, based on the average number of pollinia received per flower from the total number collected per individual (\code{pollinia.import}).}
#' 		\item{\code{At3}}{width of the corona, a trait related to the visual attraction of the flower.}
#' 		\item{\code{Gy2}}{width of the inner guide rail, a trait related to pollinaria transfer.}
#' 		\item{\code{Gy3}}{width of the inner guide rail at the level of the upper edges, a trait related to pollinaria transfer.}
#' 		\item{\code{Gy4}}{length of the inner guide rail, a trait related to pollinaria transfer.}
#' 		\item{\code{Gy5}}{length of the outer guide rail, a trait related to pollinaria transfer.}
#' 		\item{\code{Co3}}{length of the corpusculum, a measure of the pollinarium.}
#' 		\item{\code{Ca3}}{length of the caudicles, a measure of the pollinarium.}
#' 		\item{\code{individual}}{individual identification number.}
#' }
#' @source These data were obtained from "Phenotypic Selection Analysis" data file available on Dryad: \url{http://dx.doi.org/10.5061/dryad.kq72j}, and are associated with Baranzelli et al. (2014).
#' @format data.frame with 125 observations and 12 variables
#' @references Baranzelli MC, Sersic AN, Cocucci AA. 2014. The search for Pleiades in trait constellations: functional integration and phenotypic selection in the complex flowers of \emph{Morrenia brachystephana} (Apocynaceae). \emph{Journal of Evolutionary Biology} 27(4): 724-736. \url{http://dx.doi.org/10.1111/jeb.12341}
#' @usage Flowers
#' @keywords datasets
#' @examples
#' # Load the data
#' data(Flowers)
#'
#' # Look at the structure of the data.frame
#' str(Flowers)
#'
#' # Run a linear regression with wm as the response and the morphological traits as the predictors
#' lm(wm ~ ., data = Flowers[,5:11])



####### DesertPlants ########

#' @title Phenotypic selection in winter annual plants from the Sonoran Desert
#'
#' @name DesertPlants
#' @description Standardized data from the Kimball et al. (2013) study on phenotypic selection in four winter annual plant species from the Sonoran Desert that identified a trade-off between relative growth rate and water use efficiency. Traits were standardized by "substracting the mean and dividing by the standard deviation". These data were downloaded from Dryad and then formatted for use in this R package. \code{DesertPlants} consists of the following variables:
#'
#' \describe{
#' 		\item{\code{w}}{relative fitness, based on total plant biomass (\code{RelFitness}), which was used as a proxy for seed production.}
#' 		\item{\code{stdRGR}}{standardized relative growth rate.}
#' 		\item{\code{stdSLA}}{standardized specific leaf area, which was calculated as leaf area/dry mass of the leaf.}
#' 		\item{\code{stdRMR}}{standardized root mass ratio, which was calculated as root dry mass/total dry mass.}
#' 		\item{\code{stdN}}{standardized leaf Nitrogen content.}
#' 		\item{\code{stdDELTA}}{standardized integrated water-use efficiency over the lifetime of the leaf, based on carbon isotope ratios.}
#' 		\item{\code{Site}}{collection site within the Sonoran Desert. TH = University of Arizona's Desert Laboratory at Tumamoc Hill; cooler and wetter site. OPNM = Organ Pipe National Monument; warmer and drier site.}
#' 		\item{\code{Species}}{species of plant. STMI = \emph{Stylocline micropoides}, ERLA = \emph{Eriophyllum lanosum}, PERE = \emph{Pectocarya recurvata}, and ERTE = \emph{Erodium texanum}.}
#' 		\item{\code{ID}}{individual identification number within each treatment/species combination.}
#' }
#' @source These data were obtained from the "KimballetalData.xlsx" data file available on Dryad: \url{http://dx.doi.org/10.5061/dryad.c8c58}, and are associated with Kimball et al. (2013)
#' @format data.frame with 913 observations and 9 variables
#' @references Kimball S, Gremer JR, Huxman TE, Venable DL, Angert AL (2013) Phenotypic selection favors missing trait combinations in coexisting annual plants. \emph{The American Naturalist} 182(2): 191-207.  \url{http://dx.doi.org/10.1086/671058}
#' @usage DesertPlants
#' @keywords datasets
#' @examples
#' # Load the data
#' data(DesertPlants)
#'
#' # Look at the structure of the data.frame
#' str(DesertPlants)
#'
#' # Run a linear regression with w as the response and the morphological traits as the predictors
#' lm(w ~ ., data = DesertPlants[,1:6])



####### Bullfrogs ########

#' @title Reproductive success of male bullfrogs
#'
#' @name Bullfrogs
#' @description Raw data from the Howard (1979) study on sexual selection on body size in male \emph{Rana catesbeiana} bullfrogs during 1976.  These data were transcribed from Table 1 in Howard (1979) and then formatted for use in this R package. \code{Bullfrogs} contains the following variables:
#'
#'\describe{
#'	\item{\code{BodySize.mm}}{Body size of the individual, based on snout-ischium length in millimeters}
#'	\item{\code{w.TotalCops}}{Relative fitness for mating success (sexual selection), based on the total number of copulations.}
#'	\item{\code{w.TotalZygotes}}{Relative fitness for fertility per mate (natural selection), based on the total number of zygotes produced.}
#'	\item{\code{w.TotalHatchlings}}{Relative fitness for offspring survivorship (natural selection), based on the total number of hatchlings surviving.}
#'	\item{\code{W.TotalCops}}{Absolute fitness for mating success (sexual selection), based on the total number of copulations.}
#'	\item{\code{W.TotalZygotes}}{Absolute fitness for fertility per mate (natural selection), based on the total number of zygotes produced.}
#'	\item{\code{W.TotalHatchlings}}{Absolute fitness for offspring survivorship (natural selection), based on the total number of hatchlings surviving.}
#'	\item{\code{Ind}}{Individual identification number. Note: one individual had a multi-number ID (i.e., 3-7) and another was not assigned a specific ID and was designated as "NA".}
#'}
#' @source These data were obtained from Table 1 in Howard (1979): \url{http://www.jstor.org/stable/2460219}.
#' @format data.frame with 38 observations and 8 variables
#' @references Howard RD. 1979. Estimating reproductive success in natural populations. \emph{American Naturalist} 114(2):221-231. \url{http://dx.doi.org/10.1086/667586}
#' @usage Bullfrogs
#' @keywords datasets
#' @examples
#' # Load the data
#' data(Bullfrogs)
#'
#' # Look at the structure of the data.frame
#' str(Bullfrogs)
#'
#' # Run a linear regression with w.TotalCops as the response and the morphological traits as the predictors
#' lm(wTotalZygotes ~ BodySize.mm, data=Bullfrogs)

