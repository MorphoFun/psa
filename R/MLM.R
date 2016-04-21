########### Manipulated Linear Model ############
###### Function to simulate data that will generate specific output in a full linear model
###### Based on: https://stat.ethz.ch/pipermail/r-sig-teaching/2008q4/000076.html

## Use this to generate summary statistics based on the full model
## and then see how that changes as I change the #/types of parameters entered

# b = vector of values for the beta coefficients
# means = means of your predictors
# sds = std devs of your predictors


simz <- function(n, means, sds) {
	set.seed(123)
	x <- list()
	xx <- list()
	for (i in 1:length(means)) {
		x[[i]] <- scale(rnorm(n,means[i], sds[i]), center=F, scale=F)
	}
	xx <- do.call("cbind", x)
	xx <- data.frame(xx, stringsAsFactors=F)
	return(xx)
}

# Need to figure out how to apply this to a binomial case... or whether this is the proper way to
# address the issue overall.
# info for the binomial data were taken from this: http://www.r-bloggers.com/example-7-2-simulate-data-from-a-logistic-regression/

MLM <- function(x, b, ydistn = c("normal", "binomial"), intercept=1, R2 = NULL, prob = NULL) {
	bx <- list()

	for (i in 1:length(b)) {
		bx[[i]] <- (b[i]*x[,i])
		yhat <- intercept + rowSums(do.call("cbind", bx))
	}
	# compute ssr
	ssr <- sum( (yhat-mean(yhat))^2 )

	# generate errors
	set.seed(123)
	e <- list()
	if (ydistn=="normal") {
		e <- rnorm(nrow(x))
		e <- resid( glm( e ~ ., data=x, family = gaussian) )
		
		# to get R^2 of 0.8, ssr/(ssr+sse)=0.8 so sse=0.2/0.8*ssr
		e <- e* sqrt((1-R2)/R2*ssr/(sum(e^2)))
		
		# now for y
		y <- yhat + e
	}
	if (ydistn=="binomial") {
    prob <- exp(yhat)/(1 + exp(yhat))
    runis <- runif(nrow(x), 0, 1)
    y <- ifelse(runis < prob, 1, 0)
	}

	# put into a data frame and test
	mydata <- data.frame( y=y, x)
	return(mydata)

}

# # verify with
# fit <- lm(y ~ x1 + x2 + x3, data=mydata )
# summary(fit)
