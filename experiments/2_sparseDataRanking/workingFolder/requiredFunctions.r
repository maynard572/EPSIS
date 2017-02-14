####################################################################
# Required Functions


#------------------------------------------------------------------
getPdf.gaussian <- function(mu, sigma, numSig = 4, len = 1024){
	x <- seq(mu - numSig * sigma, mu + numSig * sigma, length.out = len )
	y <- dnorm(x, mean = mu, sd=sigma)
	return(list(x=x,  y=y))
	
}




#------------------------------------------------------------------
# afunction which takes a verification distribution and spits out verification variables.


getVerification.gaussian <- function(N, verificationDists){
	
	verification <- rnorm(N, mean=verificationDists$param$mu, sd = verificationDists$param$sigma)
	return(verification)

	
} 


#------------------------------------------------------------------
calcIgnorance <- function(forecastPDF, verification){
	
	f <- approxfun(forecastPDF$x, forecastPDF$y, yleft=0, yright=0)
	ig <- -log2(f(verification))
	return(ig)
}


#------------------------------------------------------------------
calcNaiveLinear <- function(forecastPDF, verification){
	
	f <- approxfun(forecastPDF$x, forecastPDF$y, yleft=0, yright=0)
	nl <- -f(verification)
	return(nl)
}

#------------------------------------------------------------------
calcProperLinear <- function(forecastPDF, verification, simpleIntegrate=FALSE){
	
	x <- forecastPDF$x
	y <- forecastPDF$y
	f <- approxfun(x,y, yleft=0, yright=0)	
	
	if(simpleIntegrate){
		IntFSquared <- simpleIntegrate(list(x=x, y=y^2))
	} else {
		fSquared <- approxfun(x, y^2 )
		IntFSquared <- integrate(fSquared, range(x)[1] , range(x)[2])$value
	}
	
	plr <- IntFSquared - 2 * f(verification)
	return(plr)
	
}

#------------------------------------------------------------------
calcPowerRule <- function(forecastPDF, verification, alpha, simpleIntegrate=FALSE){
	
	x <- forecastPDF$x
	y <- forecastPDF$y
	f <- approxfun(x,y, yleft=0, yright=0)
	
	if(simpleIntegrate){
		IntFToAlpha <- simpleIntegrate(list(x=x, y=y^alpha))
	} else {
		fToAlpha <- approxfun(x, y^alpha )
		IntFToAlpha <- integrate(fToAlpha, range(x)[1] , range(x)[2])$value
	}
	
	plr <- (alpha - 1)*IntFToAlpha - alpha * f(verification) ^(alpha - 1)
	return(plr)
	
}


#------------------------------------------------------------------
calcSpherical <- function(forecastPDF, verification, simpleIntegrate=FALSE){
	
	x <- forecastPDF$x
	y <- forecastPDF$y
	f <- approxfun(x,y, yleft=0, yright=0)
	
	if(simpleIntegrate){
		IntFSquared <- simpleIntegrate(list(x=x, y=y^2))
		
	} else {
		fSquared <- approxfun(x, y^2 )
		IntFSquared <- integrate(fSquared, range(x)[1] , range(x)[2])$value
	}
	
	plr <-  -f(verification)/ sqrt(IntFSquared)
	return(plr)
	
}

#------------------------------------------------------------------
calcContinuousRankedProbability <- function(forecastPDF, verification, simpleIntegrate=FALSE){
	
	# The Heaviside function is defined as H(x)=0, x<0; H(x)=1, x>0.  Its value at zero is left to the user depending on how they intend to use it - as we  use it here within an integral and p(y) is right continuous we'll set H(0)=1 so it is right continuous as well.
	
	#S(p(X), x) = &int (p(y) -H(x-y))² dy  - I think this definition is wrong  p(y) should be P(y) (i.e) its integral from -Inf to y.  
	
	x <- forecastPDF$x
	y <- forecastPDF$y
	p <- approxfun(x,y, yleft=0, yright=0)
	
	H <- function(x){
		H <- 0 + (x>=0)
		return(H)
	}
	

	if(simpleIntegrate){
		
		P.sgl <- function(t,p, lwr,xVals, yVals){
			xTot <- x[x<=t]
			yTot <- y[x<=t]
			P.sgl <- simpleIntegrate(list(x=xTot,y=yTot))
			return(P.sgl)
		}
		
	} else {
	
		P.sgl <- function(t,p, lwr,xVals=NULL, yVals=NULL){
			P.sgl <- integrate(f=p, lower = lwr , upper = t)$value
			return(P.sgl)
		}
	
	}
	
	
	P.vec <- Vectorize(FUN=P.sgl, vectorize.args=list("t"))
	
	integrand <- function(y, v, xVals=NULL, yVals=NULL){
		integrand <- ( P.vec(y, p, lwr=range(x)[1],xVals,yVals) - H(y-v) )^2
		return(integrand)	
	}
	
	
	if(simpleIntegrate){
		
		intg <- integrand(y=x ,v=verification, xVals=x, yVals=y)
		crps <- simpleIntegrate(list(x=x,y=intg))
		
	} else {
		crps <- integrate(f=integrand, lower = range(x)[1] , upper = range(x)[2],  v=verification)$value
	}
	
	return(crps)
	
}

#------------------------------------------------------------------
calcMeanSquaredError <- function(forecastPDF, verification, simpleIntegrate=FALSE){
	
	#S(p(X),x) = &int (x-z)² p(z)dz, also equals
	#          = (x-m)² + v

	xvals <- forecastPDF$x
	yvals <- forecastPDF$y
	p <- approxfun(xvals,yvals, yleft=0, yright=0)
		
	integrand <- function(z, x){
		integrand <- (x-z)^2 * p(z)
		return(integrand)	
	}
	
	if(simpleIntegrate){
		mse <- simpleIntegrate(list(x=xvals, y=integrand(z=xvals, x=verification)))
	} else {	
		mse <- integrate(f=integrand, lower = range(xvals)[1] , upper = range(xvals)[2],  x=verification)$value
	}
	
	return(mse)
}


#------------------------------------------------------------------
calcScore <- function(scoreType, forecastPDF, verification, scoreParams = NULL){
	
	#scoreParams a list of parameters to determine the particular score to use within a parametrised family - or NULL if not used.
	
	if(tolower(scoreType)=="ignorance") return(calcIgnorance(forecastPDF, verification))
	
	if(tolower(scoreType) == "properlinear") return(calcProperLinear(forecastPDF, verification,simpleIntegrate=scoreParams$simpleIntegrate))
	
	if(tolower(scoreType) == "continuousrankedprobability") return(calcContinuousRankedProbability(forecastPDF, verification, simpleIntegrate=scoreParams$simpleIntegrate))
	
	if(tolower(scoreType) == "meansquarederror") return(calcMeanSquaredError(forecastPDF, verification,simpleIntegrate=scoreParams$simpleIntegrate))
	
	if(tolower(scoreType) == "naivelinear") return(calcNaiveLinear(forecastPDF, verification))
	
	if(tolower(scoreType) == "powerrule") return(calcPowerRule(forecastPDF, verification, alpha = scoreParams$alpha, simpleIntegrate=scoreParams$simpleIntegrate))
	
	if(tolower(scoreType) == "spherical") return(calcSpherical(forecastPDF, verification, simpleIntegrate=scoreParams$simpleIntegrate))

	
 # add new scores here	
	
}


#------------------------------------------------------------------
getForecasts.pdf <- function(N, forecastDists, numSig, len){

	forecasts.pdf <- list()
	for (i in 1:N){
		forecasts.pdf[[i]] <- getPdf.gaussian(mu=forecastDists$param$mu[i], sigma = forecastDists$param$sigma[i], numSig = numSig, len = len ) 
		
	}
	return(forecasts.pdf)
}


#------------------------------------------------------------------
getScores <- function(N, scoreNames, forecasts.pdf, verifications, chnl, seedVal){
	
	# first check whether a seed has been used before if it has then you want this experiment to match the others - give a warning if not and user has to delete prior files...
	
	# assumes odbc connection is already open
	tables <- sqlTables(channel=chnl, catalog="scoreData")  #seems to need the catalog parameter otherwise only sees the first table
	
	seedExists <- "seedVal" %in% tables$TABLE_NAME
	if(seedExists){
		# check the seed is equal to the current seed - else stop
		oldSeedVal <- sqlFetch(chnl, sqtable="seedVal")[1,1]
		if(seedVal != oldSeedVal) stop("seeds dont match! remove old experiments?")
		
	} else {
		sqlSave(chnl, dat=data.frame(seedVal=seedVal), tablename="seedVal")
		
	}
	
	
	
	experimentName <- gsub("forecasts.pdf.", "", deparse(substitute(forecasts.pdf)))

	scoreTabName <- paste("scores_", experimentName, sep="")
	scoreValuesExist <- scoreTabName %in% tables$TABLE_NAME
	
	if(scoreValuesExist){
		
		scores <- t(sqlFetch(chnl, sqtable = scoreTabName ))
		scores[scores == 999999] <- Inf  # SQL doesnt cope with Inf
		
		
	} else {
	
		scores <- getScoresArray(N=N, scoreNames=scoreNames, forecasts.pdf=forecasts.pdf, verifications=verifications)
	
		rownames(scores) <- scoreNames
		
		scoresDat <- scores
		scoresDat[scoresDat==Inf] <- 999999
		
		sqlSave(chnl, dat=as.data.frame(t(scoresDat)), scoreTabName)
		
	
	}
	
	
	return(scores)
}


#----------------------------------------------------------------- 
getIntegerFromParmX <- function(parm, level=1, levelMax = 10){
	if (level == levelMax) stop("max number of levels reached: result not provided. increase levelMax?")
	
	if(abs(floor(parm + 0.000000001) - parm) < 0.00000001) {
		return(parm)
	} else {
		parm <- parm * 10
		level <- level + 1
		intParm <- getIntegerFromParmX(parm, level)
	}
	
	return(intParm)
	
}


getIntegerFromParmY <- function(integerVal, decimalPointPlace){
	N <- floor(log10(integerVal))
	
	D<- 10^(N - decimalPointPlace + 1)
	leftOfPoint <- floor(integerVal / (D))
	rightOfPoint <- integerVal - leftOfPoint * D
	
	return(paste(leftOfPoint, rightOfPoint, sep="o"))
}

getIntegerFromParm <- function(parm, level=1, levelMax = 10){
	parmVal <- parm
	
	integerVal <- getIntegerFromParmX(parm, level=1, levelMax = 10)
	
	decimalPointPlace <- floor(log10(parmVal))+1
	
	IntFromParm <- getIntegerFromParmY(integerVal, decimalPointPlace)
	
	return(IntFromParm)
	
}

#----------------------------------------------------------------- 
getCombinedScoreName.sgl <- function(scoreType, scoreParams){

	if(tolower(scoreType)=="ignorance") return(scoreType)
	
	if(tolower(scoreType) == "properlinear") return(scoreType)
	
	if(tolower(scoreType) == "continuousrankedprobability") return(scoreType)
	
	if(tolower(scoreType) == "meansquarederror") return(scoreType)
	
	if(tolower(scoreType) == "naivelinear") return(scoreType)
	
	if(tolower(scoreType) == "powerrule"){
		parm <- scoreParams$alpha
		parmI <- getIntegerFromParm(parm)
		scoreType <- paste(scoreType, parmI, sep="_")
		return(scoreType)
		}
	
	if(tolower(scoreType) == "spherical") return(scoreType)
	
	# Add other scores here
		
}



#------------------------------------------------------------------

getScoresArray <- function(N, scoreNames, scoreParams, forecasts.pdf, verifications){
	

	# it is assumed that scoreParams and scoreNames are in matched order; i.e. the params vector matches the scores 
	

	
	scoreNameAug <- NULL
	scores <- array(dim=c(length(scoreNames),N))
		for (i in 1:N){

			for (j in 1:length(scoreNames)){
				scores[j,i] <- calcScore(scoreType = scoreNames[j], forecastPDF = forecasts.pdf[[i]], verification = verifications[i], scoreParams=scoreParams[[j]])
				scoreNameAug <- c(scoreNameAug,  getCombinedScoreName.sgl(scoreType = scoreNames[j], scoreParams=scoreParams[[j]]))
			}
		
		}
	
		rownames(scores) <- scoreNameAug
		return(scores)	
	
}

#------------------------------------------------------------------

getScoresAllOneWayArray <- function(meanMoves, sigmaMoves, scoreNames, verificationDists, verifications, numSig, len){

	d1 <- length(meanMoves)
	d2 <- length(sigmaMoves)
	d3 <- length(scoreNames)
	d4 <- length(verifications)

	scores <- array(dim=c( d1,d2,d3,d4))
	
	for( i1 in 1:d1){
		for (i2 in 1:d2){
			forecastDists <- getForecastDists(verificationDists, type="allOneWay", typeParams = list(meanMove= meanMoves[i1] , sigmaMove=sigmaMoves[i2]))
			forecasts.pdf <- getForecasts.pdf(N=N, forecastDists = forecastDists, numSig=numSig, len=len)
		scores[i1,i2,,] <- getScoresArray(N=d4, scoreNames=scoreNames, forecasts.pdf=forecasts.pdf, verifications=verifications)	
			
		}
	}

	dimnames(scores) <- list(meanMoves=meanMoves, sigmaMoves=sigmaMoves, scoreNames=scoreNames, verifications=verifications)
	return(scores)

}

#------------------------------------------------------------------
getForecastDists.peturb <- function(verificationDists, typeParams){

	# typeParams = list(seed=?,  peturbWidth=?)
	set.seed(typeParams$seed)
	N <- length(verificationDists$param$mu)
	meanPeturb <- 1 + rnorm(N)/typeParams$peturbWidth
	sigmaPeturb <- 1 + rnorm(N)/typeParams$peturbWidth
	forecastDists <- verificationDists
	forecastDists$param$mu <- verificationDists$param$mu * meanPeturb
	forecastDists$param$sigma <- verificationDists$param$sigma * sigmaPeturb
	
	return(forecastDists)
}



#------------------------------------------------------------------
getForecastDists.allOneWay <- function(verificationDists, typeParams){
	
	
	# typeParams = list(meanMove=? , sigmaMove=?)
	
	N <- length(verificationDists$param$mu)
	meanPeturb <- 1 + typeParams$meanMove
	sigmaPeturb <- 1 + typeParams$sigmaMove
	forecastDists <- verificationDists
	forecastDists$param$mu <- verificationDists$param$mu * meanPeturb
	forecastDists$param$sigma <- verificationDists$param$sigma * sigmaPeturb
	
	return(forecastDists)
}



#------------------------------------------------------------------
getForecastDists <- function(verificationDists, type, typeParams=NULL){
	
	if(type == "equalsVerification") return(verificationDists)
	
	if(type == "peturb") return(getForecastDists.peturb(verificationDists, typeParams=typeParams))
	
	if(type == "allOneWay") return(getForecastDists.allOneWay(verificationDists, typeParams=typeParams))
	
	# other types here ….
	
}



#------------------------------------------------------------------
createLegendInMargin <- function(desc, col, cex=0.5, startAtLine=2){
	# assumes that the xlab is in line 1 so that lines 2,3,4 are available
	# desc a list of descriptions (strings)
	# col a list of colours for the descriptions
	
	N <- length(desc)
	if (length(col)!=N) stop("col vector not same length as desc")
	
	maxSlots <- (4-startAtLine)*2 + 1
	numBlocks <- ceiling( N / maxSlots)
	usr <- par("usr")
	blockSize <- (usr[2] - usr[1])/numBlocks
	
	block <- 0
	line <- 2
	for (d in seq(along=desc)){
		at = usr[1] + block * blockSize
		mtext(at=at, side=1, text=paste("---", desc[[d]], sep="") , line=line, adj=0, col=col[d], cex=0.5)
		
		line <- line + 0.5
		if (line == 4.5){
			line <- 2
			block <- block + 1			
			}	
		
		}
	return(NULL)
	}
	
	
#------------------------------------------------------------------	
	
plotNormalisedScores <- function(scores, verifications, mu, muhat=mu, desc=scoreNames.prty){
	
	scores.normalise <- scores / abs(apply(scores,1,median))
#	mx <- max(apply(scores.normalise,1,max))
#	mn <- min(apply(scores.normalise,1,min))

	mx <- max(scores.normalise[scores.normalise!=Inf])
	mn <- min(scores.normalise[scores.normalise!=Inf])
	x <- range(verifications)


	plot(x, c(mn,mx) , type="n", xlab="verification", ylab="score(normalised)")
	for (i in 1:5){
		points(verifications, scores.normalise[i,] , pch=".", col=i)
	}
	createLegendInMargin(desc=desc, col=1:5, cex=0.5, startAtLine=2)
	mtext(text=formatC(mu, digits=2, format="fg"), side=1, at=mu, col="orange")
	abline(v=mu, col="orange")

	mtext(text=formatC(muhat, digits=2, format="fg"), side=1, at=muhat, col="olivedrab")
	abline(v=muhat, col="olivedrab")
	
}
	
	
	
#------------------------------------------------------------------	
	
getPairwiseScorePlot <- function(scores, scoreNames = scoreNames.prty){


	panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
	{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, method="spearman"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * max(r, 0.5))
	}


	panel.xy <- function (x, y){ 
    points(x, y, pch = ".")
	}

	scr <- t(scores)
	colnames(scr) <- scoreNames.prty 
	pairs(scr, lower.panel=panel.xy, upper.panel=panel.cor)

}


#------------------------------------------------------------------

mean1 <- function(x){x[x==Inf]<-NA;return(mean(x, na.rm=TRUE))}

#------------------------------------------------------------------
mean999 <- function(x){x[x==Inf]<-999;return(mean(x, na.rm=TRUE))}


#------------------------------------------------------------------
kernelPDF.sgl <- function(t, sigmaK, S){
	# S = set of data points
	N <- length(S)
	f <- sum(dnorm((t-S)/sigmaK))/(N*sigmaK)
	return(f)
	
}

#------------------------------------------------------------------
kernelPDF <- Vectorize(FUN=kernelPDF.sgl, vectorize.args=c("t"))


#------------------------------------------------------------------
getVerificationFromElementX <- function(x, sigmaK){
	N <- length(x)
	
	eta <- rnorm(N,mean=0, sd=rep(sigmaK,N))
	return(x+eta)
}

#------------------------------------------------------------------
getForecastPDFM <- function(S, sigmaM, len, numSig=NULL, xRange=NULL){
		if(is.null(xRange)) xRange <- range(S) + c(-numSig*sigmaM, numSig*sigmaM)
		xVals  <- seq(xRange[1], xRange[2], length.out=len)
		yVals  <- kernelPDF(xVals, sigmaK=sigmaM, S=S)
		forecastPDFM <- list(x=xVals, y=yVals) 
		return(forecastPDFM)
}

#------------------------------------------------------------------

getMin <- function(x,y){
	# assumes a functional relationship between x and y
	# returns the x value corresponding to the total minimum of y (i.e. local minima ignored)
	# if there are more than one minima (where the y value is equal) then all x values are returned, corresponding to the one y value
	
	g <- data.frame(x,y)
	g <- g[order(g[,2]),]
	gMin <- g[g[,2]==g[1,2],]
	yMin <- gMin[1,2]
	xMins <- gMin[,1]
	
	return(list(y=yMin, x=xMins))
	
}


#------------------------------------------------------------------
# data frame version

getKernelSmoothedScores.dataframe <- function(
	S_V,
	S_E,
	sigmaU,
	M,
	len,
	sN,
	sP,
	seed.Sm,
	seed.v,
	sigmaMs){
	
	scoreVals <- data.frame()
	#dimnames(scoreVals) <- list(sigmaMs,NULL, sN)





	# get verifications
	#  - the same verifications will be used in each of the forecasts
	set.seed(seed.Sm)
	Sm <- sample(S_V, M, replace=TRUE)
	set.seed(seed.v)
	v <- getVerificationFromElementX(x=Sm , sigmaK= sigmaU)


	t<- 1 # initialise counter
	for (sigmaM in sigmaMs){
	
		# get forecast distribution p(sigmaM)
		forecastPDFM <- getForecastPDFM(S=S_E, sigmaM = sigmaM, len=len, xRange = range(v))

		for (i in 1:M){
		
			# calc scores
			scAry <- getScoresArray(N=1, scoreNames = sN, scoreParams = sP, forecasts.pdf=list(forecastPDFM), verifications = c(v[i]))
			scAry.df <- data.frame(sigmaM = sigmaM, verification = i, scoreType=rownames(scAry), scoreValue = scAry[,1])
			
			scoreVals <- rbind(scoreVals,scAry.df )
	
			print(paste("i:", i, "  t:", t, sep=""))
	
		}  # next i
	
		t <- t+1
	} #next sigmaM
	rownames(scoreVals) <- NULL
	return(scoreVals)


}



#------------------------------------------------------------------

getMinVals <- function(meanVals){

	splitFac <- as.factor(paste(meanVals$seedIndex,meanVals$fcstSize,meanVals$scoreType, sep="_"))


	meanVals.split <- split(meanVals, f=splitFac)

	getMinScore <- function(mvb){
		minxy <- getMin(x=mvb$sigmaM,  y=mvb$scoreValue)
		rslt <- data.frame(
				seedIndex = mvb$seedIndex[1],
				fcstSize = mvb$fcstSize,
				scoreType = mvb$scoreType,
				xmin = minxy$x,
				ymin = minxy$y )
		return(rslt)
	}


	minVals <- do.call(rbind, lapply(meanVals.split, getMinScore))
	rownames(minVals) <- NULL
	return(minVals)
}

#------------------------------------------------------------------

getMinVals12 <- function(meanVals){

	splitFac <- as.factor(paste(meanVals$seedIndex,meanVals$fcstSize,meanVals$fcstPos,meanVals$scoreType, sep="_"))


	meanVals.split <- split(meanVals, f=splitFac)

	getMinScore <- function(mvb){
		minxy <- getMin(x=mvb$sigmaM,  y=mvb$scoreValue)
		rslt <- data.frame(
				seedIndex = mvb$seedIndex[1],
				fcstSize = mvb$fcstSize,
				fcstPos = mvb$fcstPos,
				scoreType = mvb$scoreType,
				xmin = minxy$x,
				ymin = minxy$y )
		return(rslt)
	}


	minVals <- do.call(rbind, lapply(meanVals.split, getMinScore))
	rownames(minVals) <- NULL
	return(minVals)
}


#------------------------------------------------------------------

getMinVals13 <- function(meanVals){

	splitFac <- as.factor(paste(meanVals$seedIndex,meanVals$fcstSize,meanVals$fcstPos,meanVals$vPower,meanVals$scoreType, sep="_"))


	meanVals.split <- split(meanVals, f=splitFac)

	getMinScore <- function(mvb){
		minxy <- getMin(x=mvb$sigmaM,  y=mvb$scoreValue)
		rslt <- data.frame(
				seedIndex = mvb$seedIndex[1],
				fcstSize = mvb$fcstSize,
				fcstPos = mvb$fcstPos,
				vPower = mvb$vPower,
				scoreType = mvb$scoreType,
				xmin = minxy$x,
				ymin = minxy$y )
		return(rslt)
	}


	minVals <- do.call(rbind, lapply(meanVals.split, getMinScore))
	rownames(minVals) <- NULL
	return(minVals)
}


#-----------------------------------------------------------------
getMinVals19 <- function(meanVals){

	splitFac <- as.factor(paste(meanVals$initialisation,meanVals$seedIndex,meanVals$fcstSize,meanVals$fcstPos,meanVals$vPower,meanVals$scoreType, sep="_"))


	meanVals.split <- split(meanVals, f=splitFac)

	getMinScore <- function(mvb){
		minxy <- getMin(x=mvb$sigmaM,  y=mvb$scoreValue)
		rslt <- data.frame(
				initialisation = mvb$initialisation[1],
				seedIndex = mvb$seedIndex[1],
				fcstSize = mvb$fcstSize[1],
				fcstPos = mvb$fcstPos[1],
				vPower = mvb$vPower[1],
				scoreType = mvb$scoreType[1],
				xmin = minxy$x,
				ymin = minxy$y )
		return(rslt)
	}


	minVals <- do.call(rbind, lapply(meanVals.split, getMinScore))
	rownames(minVals) <- NULL
	return(minVals)
}


#------------------------------------------------------------------

compareScores <- function(s1,s2,s1Name, s2Name, trueParm, pch=19, forecastQuality=NULL, cex=0.5, rad1=NULL, rad2=NULL, radVals=NULL, box=NULL){
	
	#forecastQuality is a number between 0 (very poor) and 1 (very good)
	
	if(is.null(forecastQuality)) forecastQuality = rep(1, length(s1))
	f <- colorRamp(colors = c("red", "blue"))
	cols <- rgb(f(forecastQuality), maxColorValue=255)
	
	
	axRange <- range(s1,s2)
	plot(axRange,axRange, type="n", xlab=s1Name, ylab=s2Name)
	
	
	polyLimit <- 5 * max(abs(axRange[1]), abs(axRange[2]))  #make this more clever
	
	polygon(x=c(-polyLimit,trueParm,polyLimit + 2*trueParm,-polyLimit), y=c(-polyLimit,trueParm,-polyLimit,-polyLimit), col="palegoldenrod", border=NA)
	
	polygon(x=c(-polyLimit,trueParm,polyLimit + 2*trueParm,-polyLimit), y=c(polyLimit + 2 * trueParm,trueParm,polyLimit + 2 * trueParm,polyLimit + 2 * trueParm), col="palegoldenrod", border=NA)

	
	abline(v=trueParm, col="grey")
	abline(h=trueParm, col="grey")
	
	if(length(pch)==1){
		points(s1,s2,pch=pch, col=cols)
	} else {
		text(s1,s2,labels=pch, col=cols, cex=cex)
		
	}
	
	
	if(!is.null(rad1)){
		# want to draw elipses on the gaphics
		
		if(is.null(rad2)) stop("you must also specify rad2 if you want to use rad 1")
		if(is.null(radVals)) stop("you must also specify radVals if you want to use rad 1")
		
		N <- length(rad1)
		for (k in 1:N) {
			elps <- getEllipse(a=rad1[k], b=rad2[k],c=c(trueParm, trueParm), box=box)
			lines(elps$x, elps$y, col="grey")
			text(elps$x[1], elps$y[1], labels=radVals[k], col="grey")
			
		}
		
		
		
	}
	
	
	
}


#------------------------------------------------------------------

distanceFromTruth <- function( sigmaU, sigmaM, S_V, S_E, L=2){
	
	if (L==1){
		f <- function(x, sigmaU, sigmaM, S_V, S_E, L=2){
			val <- abs(kernelPDF(x, sigmaU ,S=S_V) - kernelPDF(x, sigmaM ,S=S_E) )
			return(val)
			}		
		
	}
	
	if (L==Inf) stop("sup norm not implemented - cant use L = Inf")

	if (L>1 & L<Inf){
		f <- function(x, sigmaU, sigmaM, S_V, S_E, L=2){
			val <- ( kernelPDF(x, sigmaU ,S=S_V) - kernelPDF(x, sigmaM ,S=S_E) )^L
			return(val)
			}
	}

	
	fv <- Vectorize(FUN=f, vectorize.args=c("x"))
	
	allSamples <- c(S_V,S_E)
	lower <- min(allSamples) - sd(allSamples)
	upper <- max(allSamples) + sd(allSamples)

	dist <- (integrate(f=fv , lower=lower, upper=upper, sigmaU=sigmaU, sigmaM=sigmaM, S_V=S_V, S_E=S_E, L=L)$value)^(1/L)
	
	return(dist)
	
}


#-----------------------------------------------------------------

getEllipse <- function(a,b,theta=0, c=c(0,0), len=100, box=FALSE){
	
	
	# if box is FALSE
	# produces xy points on an ellipse with axes of lengths a and b
	# rotated through an angle theta
	# centred on the point c (expressed as c=c(x,y))
	
	# if box is TRUE 
	# produces a rectangle with lengths a and b

	
	
	if(!box){
	
	getY <- function(x,a,b){
		y <- (  (1/1)^2 - (x/a)^2 )^(1/2) * b	
		return(y)
	}
	
	getX <- function(y,a,b){
		x <- (  (1/1)^2 - (y/b)^2 )^(1/2) * a	
		return(x)
		
	}
	
	x1Vals <- seq( -a/(2) , a/(2) , length.out = len)
	y1Vals <- getY(x1Vals,a,b)
	
	y2Vals <- c(
		seq(getY(a/(2),a,b),0, length.out=len/2 ), 
		seq(0,-getY(a/(2),a,b), length.out=len/2 ))
	x2Vals <- getX(y2Vals,a,b)
	
	x3Vals <- seq(a/(2), -a/(2), length.out = len)
	y3Vals <- - getY(x3Vals,a,b)	
	
	y4Vals <- c(
		seq(-getY(-a/(2),a,b),0, length.out=len/2 ), 
		seq(0,getY(-a/(2),a,b), length.out=len/2 ))
	x4Vals <- - getX(y4Vals,a,b)
	

	
	} else {
		
		x1Vals <- c(-a,a)
		x2Vals <- c(a,a)
		x3Vals <- c(a,-a)
		x4Vals <- c(-a,-a)
		
		y1Vals <- c(b,b)
		y2Vals <- c(b,-b)
		y3Vals <- c(-b,-b)
		y4Vals <- c(-b,b)
		
		
	}
	
	
	
	xVals <- c(x1Vals, x2Vals,x3Vals,x4Vals )
	yVals <- c(y1Vals, y2Vals,y3Vals,y4Vals )	
	
	g <- matrix(c(xVals, yVals), nrow=2, byrow=TRUE)
	k <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2, byrow=TRUE)
	
	g2 <- k %*% g
	
	
	
	return(list(x=g2[1,]+c[1], y=g2[2,]+c[2]))
	
	
}



#-----------------------------------------------------------------

getRadii <- function(
	radCol,
	s1Name,
	s2Name,
	trueParm,
	FUN = "mean"
){

	minMean <- aggregate(minVals$xmin, by=list(minVals[,radCol], scoreType=minVals$scoreType), FUN=FUN)
	names(minMean)[1] <- radCol

	radVals <- unique(minVals[, radCol])
	s1Tab <- minMean[minMean$scoreType == s1Name,]
	s2Tab <- minMean[minMean$scoreType == s2Name,]

	rads1 <- abs(s1Tab[match(s1Tab[, radCol], radVals), "x"] - trueParm)
	rads2 <- abs(s2Tab[match(s2Tab[, radCol], radVals), "x"] - trueParm)

	return(list(radVals=radVals, rads1=rads1, rads2=rads2))
}


#-----------------------------------------------------------------


getForecastPDFPlot <- function(
	xRange,
	yRange,
	underlying,
	forecast,
	v,
	S_E,
	d=d, 
	L=L
){
	plot(xRange, yRange, type="n", xlab="", ylab="Density")
	lines(forecast$x, forecast$y, col="darkgreen", lwd=1)
	lines(underlying$x, underlying$y, col="red", lwd=2)	
	rug(S_E[1:2^8])	
	rug(v, side=3, ticksize = 0.02, col="olivedrab", lwd=0.2)
	title(sub=paste("L",L, " distance = ", d, sep=""))
}

#-----------------------------------------------------------------


getFcstPDFFileName <- function(SType, typeParms,mk, k2,r,p, sigmaM ){
	
	expString <- getExpString(SType, typeParms,mk, k2,r,p, sigmaM )
	
	fileName <- paste(
		"fcstPDF_",
		expString,
		".pdf", 
		sep="")
	
	return(fileName)
	
}


#-----------------------------------------------------------------
getExpString <- function(SType, typeParms,mk, k2,r,p , sigmaM){
	
	if (SType == "Gaussian"){
		root1 <- paste("_seedS_",typeParms$seedS, "_lenS_", typeParms$lenS, sep="" )
		
	} else {
		
		root1 <- paste("_nStep_", typeParms$nStep, "_dufCol_", typeParms$dufCol, sep="")
		
	}
	
	expString <- paste(
		
		"SType_", SType,
		root1,
		"_mk_", mk,
		"_k2_", k2,
		"_r_", r,
		"_p_", p,
		"_sigmaM_", sigmaM*1000,
		sep="")
	
	return(expString)
	
}


#-----------------------------------------------------------------
getVforPlot <- function(seed.Sm, seed.v, sigmaU, S_V, M){
	set.seed(seed.Sm)
	Sm <- sample(S_V, M, replace=TRUE)
	set.seed(seed.v)
	v <- getVerificationFromElementX(x=Sm , sigmaK= sigmaU)
	return(v)
}

#-----------------------------------------------------------------


getPlotFromSamplingLoops <- function(SType, typeParms,runInfo, S){
	S_V <- S[(1:2^12)]
	sigmaMs <- typeParms$sigmaMs
	for (mk in runInfo$mks){
		M <- 2^mk   # the number of verifications

		for (k2 in runInfo$k2s){

			v <- getVforPlot(
						seed.Sm= runInfo$seed.Sm + k2, 
						seed.v = runInfo$seed.v, 
						sigmaU = typeParms$sigmaU, 
						S_V = S_V, 
						M = M)
			
			for (r in runInfo$rs){

				for (p in runInfo$ps){

					for (sigmaM in sigmaMs){

						S_E = S[seq(2^12+1+2^r * p, 2^12 + 2^r +2^r * p)]
				
						d <- distanceFromTruth(
							sigmaU = typeParms$sigmaU,
							sigmaM = sigmaM,
							S_V = S_V,
							S_E = S_E,
							L = runInfo$L
						)
				
						underlying <- getForecastPDFM(
							S=S_V, 
							sigmaM = typeParms$sigmaU, 
							len=runInfo$len, 
							xRange=c(min(S_V)*1.2, max(S_V)*1.2))
						
						forecast <- getForecastPDFM(
							S=S_E, 
							sigmaM = sigmaM, 
							len=runInfo$len, 
							xRange=c(min(S_V)*1.2, max(S_V)*1.2))
	
						yRange <- c(
							min(underlying$y,  forecast$y),
							max(underlying$y,  forecast$y))

						xRange <- c(
							min(underlying$x, forecast$x),
							max(underlying$x, forecast$x))

					
						fileName <- getFcstPDFFileName(SType=SType, typeParms=typeParms,mk=mk, k2=k2,r=r,p=p, sigmaM=sigmaM )		
				
						print(fileName)
				
						pdf(fileName)
							getForecastPDFPlot(
								xRange = xRange,
								yRange = yRange,
								underlying = underlying,
								forecast = forecast,
								v = v,
								S_E = S_E, 
								d=d, 
								L=runInfo$L)
						dev.off()	
			
					} # next sigmaM				
				} # next p							
			} # next r 							
		} # next k2							
	} # next mk	
}



#-----------------------------------------------------------------



getForecastPDFGraphic <- function(
	runInfo,
	mainFolder = "/Users/trevormaynard/Documents/Phd/robustIndex/SkillScores/"){

	require(grDevices)
	setwd( paste(mainFolder, "requiredFunctions", sep=""))
	source("requiredFunctions.r")
	setwd( paste(mainFolder, "ProduceForecastPDFs/pdfs", sep=""))
	require(RODBC)
	chnl <- odbcConnect("scoreData")

	for (SType in runInfo$STypes){

		if (SType == "Gaussian"){
			set.seed(runInfo$seed.S)
			S <- rnorm(2^runInfo$lenS, 0,1) 

				getPlotFromSamplingLoops(
					SType=SType, 
					typeParms=list(seedS = runInfo$seed.S, lenS = runInfo$lenS, sigmaMs = runInfo$sigmaMsG, sigmaU = runInfo$sigmaU$Gaussian),
					runInfo=runInfo,
					S=S)
	
		} else {
		
			for (nStep in runInfo$nSteps){ # either 8 or 32
				
				dufDat <- sqlFetch(channel=chnl, sqtable=paste("dufDat_", nStep,sep=""))

				for (dufCol in runInfo$dufCols){ #any number from 1 to 32
	
					S <- dufDat[,dufCol][1:2^runInfo$lenS]
				
					getPlotFromSamplingLoops(
						SType=SType,
						typeParms = list(nStep=nStep, dufCol=dufCol, sigmaMs = runInfo$sigmaMsD, sigmaU = runInfo$sigmaU$Duffing),
						runInfo=runInfo,
						S=S)	
						
				} # next nStep
			} # next dufCol
	
		} # endif 

	} # next SType


} # end function

#-----------------------------------------------------------------

getID <- function(x){
	return(paste(x, collapse="_"))
}

#-----------------------------------------------------------------

getScoreComparisonResult <- function(
	S1,
	S2,
	S1Name = S1,
	S2Name = S2,
	sigmaU,
	scoreResults,
	wantLatex = TRUE,
	rowDescription = "Row description"
){

e <- scoreResults

e1 <- e[e$scoreType == S1,]
e2 <- e[e$scoreType == S2,]

e1 <- cbind(e1, d=abs(e1$xmin - sigmaU))
e2 <- cbind(e2, d=abs(e2$xmin - sigmaU))

N1 <- dim(e1)[2]-5
ID1 <- apply(e1[,1:N1],1,getID)
e1 <- cbind(e1, ID=ID1)

N2 <- dim(e2)[2]-5
ID2 <- apply(e2[,1:N1],1,getID)
e2 <- cbind(e2, ID=ID2)


g <- data.frame(e1[1:N1], d1=e1$d, d2=e2[match(e1$ID, e2$ID),"d"])

s1BeatS2   <- sum(g$d1 < g$d2)
s1EqualS2  <- sum(g$d1 == g$d2)
s1LoseToS2 <- sum(g$d1 > g$d2)

if(wantLatex){
	rslt <- paste(rowDescription, " & ", S1Name, " vs ", S2Name,  " & " , s1BeatS2 , " & ",s1EqualS2 , " & " ,s1LoseToS2, " \\\\", sep="")
} else {
	rslt <- list(s1BeatS2 = s1BeatS2, s1EqualS2=s1EqualS2, s1LoseToS2=s1LoseToS2 )
	
}


return(rslt)
}

#-----------------------------------------------------------------

getScoreComparisonResult2 <- function(
	S1,
	S2,
	S1Name = S1,
	S2Name = S2,
	sigmaU,
	scoreResults,
	wantLatex = TRUE,
	rowDescription = "Row description"
){

e <- scoreResults

e1 <- e[e$scoreType == S1,]
e2 <- e[e$scoreType == S2,]

e1 <- cbind(e1, d=abs(e1$xmin - sigmaU))
e2 <- cbind(e2, d=abs(e2$xmin - sigmaU))

N1 <- dim(e1)[2]-5
ID1 <- apply(e1[,1:N1],1,getID)
e1 <- cbind(e1, ID=ID1)

N2 <- dim(e2)[2]-5
ID2 <- apply(e2[,1:N1],1,getID)
e2 <- cbind(e2, ID=ID2)


g <- data.frame(e1[1:N1], d1=e1$d, d2=e2[match(e1$ID, e2$ID),"d"])

s1BeatS2   <- sum(g$d1 < g$d2)
s1EqualS2  <- sum(g$d1 == g$d2)
s1LoseToS2 <- sum(g$d1 > g$d2)

s1Ratio <- (s1BeatS2 + s1EqualS2/2) / (s1BeatS2 + s1EqualS2 + s1LoseToS2  ) 
s1Ratio <- formatC(s1Ratio, digits=3, format="f")

if(wantLatex){
	rslt <- paste(rowDescription, " & ", S1Name, " vs ", S2Name,  " & " , s1BeatS2 , " & ",s1EqualS2 , " & " ,s1LoseToS2, " & ", s1Ratio,  " \\\\", sep="")
} else {
	rslt <- list(s1BeatS2 = s1BeatS2, s1EqualS2=s1EqualS2, s1LoseToS2=s1LoseToS2 )
	
}


return(rslt)
}


#-----------------------------------------------------------------

readInExperiment <- function(root, chnl){
	print("this may take a couple of minutes")

	
	rm(scoreVals)
	if (!exists(x="scoreVals")) {
		scoreVals <- sqlFetch(channel = chnl, sqtable = paste(root, "scoreVals",sep=""))
	
	}


#	expectedFields <- c("sigmaU" ,"nStep", "initialisation", "seedIndex", "fcstSize", "fcstPos", "vPower")
	
	
	
	#create NA values if expected fields arent available
#	for (eF in expectedFields){
		
#		if( !(eF %in% names(scoreVals)) ) {
#			print(eF)
#			scoreVals <- cbind(x=999, scoreVals)	
#			names(scoreVals)[1]<-eF
#		} 
#	}


	meanVals <- aggregate(
		scoreVals$scoreValue, 
		by=list(
			initialisation= scoreVals$initialisation,
			nStep = scoreVals$nStep,
			sigmaU = scoreVals$sigmaU,
			seedIndex = scoreVals$seedIndex, 
			fcstSize = scoreVals$fcstSize,
			fcstPos = scoreVals$fcstPos, 
			vPower = scoreVals$vPower, 
			sigmaM = scoreVals$sigmaM, 
			scoreType = scoreVals$scoreType),
		FUN = mean999)

	colnames(meanVals)[10] <- "scoreValue"

	minVals <- getMinVals13(meanVals)	
	minVals <- cbind(
		minVals, 
		fcstID = paste(minVals$fcstSize, minVals$xmin, sep="_"))
	
	minVals <- minVals[!duplicated(minVals),]

	return(
		list(
			minVals=minVals, 
			meanVals=meanVals, 
			scoreVals=scoreVals))
}


#-----------------------------------------------------------------


getTestTex <- function(
	rowDescription, 
	sigmaU, scoreResults, 
	wantLatex=TRUE, 
	file="testTextOutput", 
	sN=c("ignorance", "continuousrankedprobability", "properlinear"), 
	sNPretty = c("ignorance", "CRPS", "proper linear"))
	{

 

scoreCombn <- combn(sN,2)
scoreCombnPretty <- combn(sNPretty,2)
numComb <- dim(scoreCombn)[2]

cat("   \n", file=file, append=TRUE)
cat(paste("%", rowDescription, "\n", sep=""), file=file, append=TRUE)
cat("%---------------------------------------------\n", file=file, append=TRUE)


for(i in 1:numComb){

	sR <- getScoreComparisonResult(
		S1=scoreCombn[1,i],
		S2=scoreCombn[2,i],
		S1Name = scoreCombnPretty[1,i],
		S2Name = scoreCombnPretty[2,i],
		sigmaU=sigmaU,
		scoreResults=scoreResults,
		wantLatex = TRUE,
		rowDescription = rowDescription
	)

print(sR)
#cat(paste(rowDescription, "\n", sep=""), file=file, append=TRUE)

cat(sR, file=file, append=TRUE)
cat("\n", file=file, append=TRUE)


}

}



#-----------------------------------------------------------------


getTestTex2 <- function(
	rowDescription, 
	sigmaU, scoreResults, 
	wantLatex=TRUE, 
	file="testTextOutput", 
	sN=c("ignorance", "continuousrankedprobability", "properlinear"), 
	sNPretty = c("ignorance", "CRPS", "proper linear"))
	{

 

scoreCombn <- combn(sN,2)
scoreCombnPretty <- combn(sNPretty,2)
numComb <- dim(scoreCombn)[2]

cat("   \n", file=file, append=TRUE)
cat(paste("%", rowDescription, "\n", sep=""), file=file, append=TRUE)
cat("%---------------------------------------------\n", file=file, append=TRUE)


for(i in 1:numComb){

	sR <- getScoreComparisonResult2(
		S1=scoreCombn[1,i],
		S2=scoreCombn[2,i],
		S1Name = scoreCombnPretty[1,i],
		S2Name = scoreCombnPretty[2,i],
		sigmaU=sigmaU,
		scoreResults=scoreResults,
		wantLatex = TRUE,
		rowDescription = rowDescription
	)

print(sR)
#cat(paste(rowDescription, "\n", sep=""), file=file, append=TRUE)

cat(sR, file=file, append=TRUE)
cat("\n", file=file, append=TRUE)


}

}




#-----------------------------------------------------------------

getBimodal <- function(mu1, sigma1, mu2, sigma2, w1=0.5, numSig = 4, len = 1024, makeSparse=FALSE ){

	# always make mu1 the smallest
	if (mu1 > mu2){
		temp <- mu2
		mu2 <- mu1
		mu1 <- temp
	}
	
	
	
	x <- seq(mu1 - numSig * sigma1, mu2 + numSig * sigma2, length.out = len )
	y <- w1 * dnorm(x, mean = mu1, sd=sigma1) + (1-w1)* dnorm(x, mean = mu2, sd=sigma2) 
	
	if(makeSparse){
		tol <- 0.3
		sparcePropn <- 0.1
		slope <- diff(y)/diff(x)
		slope2 <- diff(slope)/diff(x)[-1]
		xVarying <- x[abs(slope2) > tol]
		yVarying <- y[abs(slope2) > tol]
		xStable <-  x[abs(slope2) <= tol]
		yStable <-  y[abs(slope2) <= tol]
		sparseSet <- ceiling(seq(1,length(xStable), length.out = floor(length(xStable)*sparcePropn)))
		xStable <- xStable[sparseSet]
		yStable <- yStable[sparseSet]
		
		g <- data.frame(x=c(xVarying, xStable), y=c(yVarying, yStable))
		g <- g[order(g$x),]
		
		x <- g$x
		y <- g$y
		
	}
	
	
	
	return(list(x=x,  y=y))
	
}


#-----------------------------------------------------------------

calcExpectedScore <- function(scoreType, forecastPDF, underlyingPDF, scoreParams){


	sT <- scoreType
	fP <- forecastPDF
	uX <- underlyingPDF
	sP <- scoreParams


	integrand <- function(z, sT, fP, sP, uX){
	
		uXpdf <- approxfun(uX$x, uX$y)
	
		intg <- calcScore(scoreType = sT, forecastPDF = fP, verification = z,  scoreParams = sP)  * uXpdf(z)
		return(intg)
	}

	expScore <- integrate( integrand, lower = min(uX$x), upper = max(uX$x), sT=scoreType, fP=fP, sP=NULL, uX=uX )

	return(expScore)

}


#-----------------------------------------------------------------



#-----------------------------------------------------------------

calcExpectedScore2 <- function(scoreType, forecastPDF, underlyingPDF, scoreParams, simpleIntegrate=TRUE){


	sT <- scoreType
	fP <- forecastPDF
	uX <- underlyingPDF
	sP <- scoreParams


	integrand <- function(z, sT, fP, sP, uX){
	
		uXpdf <- approxfun(uX$x, uX$y)
	
		score_z <- calcScore(scoreType = sT, forecastPDF = fP, verification = z,  scoreParams = sP)

		score_z[score_z == Inf] <- 9999 # New code 20/3/13

	
		intg <-  score_z * uXpdf(z)
		#if(is.nan(intg)) intg <- 0  # New code 20/3/13
		return(intg)
	}
	
	if(sT=="continuousrankedprobability"){
		int.sgl <- integrand
		integrand <- Vectorize(int.sgl, vectorize.args = "z")
		
	}
	
	if(simpleIntegrate){
		
		x <- sort(uX$x)
		y <- integrand(z=x, sT=scoreType, fP=fP, sP=sP, uX=uX  )
		
		expScore <- simpleIntegrate(list(x=x,y=y)) 
		
	} else {
		expScore <- integrate( integrand, lower = min(uX$x)*0.99, upper = max(uX$x)*0.99, sT=scoreType, fP=fP, sP=sP, uX=uX , subdivisions=10000)$value #new code - add * 0.99 20/3/13
	}


	return(expScore)

}


#--------------------------------------------------------------------------------------

getTruthsRanges <- function(truths){

	getRange <- function(truth){
		return(range(truth$truth$x))
	}

	ranges <- lapply(truths, getRange)
	rngTruths <- range(do.call(c, ranges))
	
	return(rngTruths)
}

#---------------------------------------------------------------------------------------
getFcstNumSig <- function(truths, fcst_m1, fcst_s1, fcst_m2, fcst_s2){
	rngTruths <- getTruthsRanges(truths)
	
	numSig1 <- ceiling((fcst_m1 - rngTruths[1])/fcst_s1) + 1
	numSig2 <- ceiling((rngTruths[2]-fcst_m2)/fcst_s2) + 1
	
	return(max(numSig1, numSig2))
	
}

#---------------------------------------------------------------------------------------
getSepSpreadFromTruth <- function(truths){
	
	getSepSprd <- function(truth){
		return(data.frame(sep=truth$sep, spread=truth$spread))
	}
	
	return(do.call(rbind, lapply(truths, getSepSprd)))
}

#---------------------------------------------------------------------------------------
getExpectedScoreResultGrid <- function(experimentName, scoreType, scoreParam, truths, expectedScores, vsBestScore){
	
	sepSpreads <- getSepSpreadFromTruth(truths)
	eScores <- data.frame(expectedScores = do.call(rbind, expectedScores))
	vsBestScore <- data.frame(vsBestScore = do.call(rbind, vsBestScore))
	scoreType <- data.frame(scoreType=rep(scoreType, length(vsBestScore)))
	experimentName <- data.frame(experimentName=rep(experimentName, length(vsBestScore)))
	
	if(is.null(scoreParam$alpha)) {
		scoreParam <- data.frame(scoreParam=rep(NA, length(vsBestScore)))
	} else {
		scoreParam <- data.frame(scoreParam=rep(scoreParam$alpha, length(vsBestScore)))
	}	

	
	return(cbind(experimentName, scoreType, scoreParam,sepSpreads, eScores, vsBestScore))
	
	
}

#---------------------------------------------------------------------------------------
simpleIntegrate <- function(func){
	# assumes a function y=f(x) has been discretised into x and y values
	# assumes the function is actually continuous and the y values are equal to
	#  the function f evaluated at the x values.
	
	x <- func$x
	y <- func$y
	lX <- length(x)
	
	if(lX==1) {
		int <- 0
	}
	
	if(lX==2){
		int1 <-  diff(x)[1] * y[1]
		int2 <-  diff(x)[1] * y[2]
		int <- (int1+int2)/2	
	}
	
	if(lX>2){
		int1 <- sum(diff(x) * y[1:(lX-1)])
		int2 <- sum(diff(x) * y[2:lX])
		int <- (int1+int2)/2
	}
	
	return(int)

}


#--------------------------------------------------------------------------------------
getMedian <- function(fcst){

	intToN <- function(N,fcst){
		return(simpleIntegrate(list(x=fcst$x[1:N], y=fcst$y[1:N])))
	}

	dfFromHalf <- function(N, fcst){
		return((intToN(N,fcst)-0.5)^2)
	}

	dfHlf <- Vectorize(dfFromHalf, vectorize.args=list("N"))

	NVals <- 1:length(fcst$x)
	dfHlf_At_NVals <- dfHlf(N=NVals, fcst)

	Nmin <- NVals[dfHlf_At_NVals == min(dfHlf_At_NVals)]

	if(dfHlf_At_NVals[Nmin]>=0.5){
		Ni <- Nmin - 1
	} else {
		Ni <- Nmin + 1
	}

	x1 <- fcst$x[Nmin]
	x2 <- fcst$x[Ni]
	y1 <- intToN(Nmin, fcst)
	y2 <- intToN(Ni, fcst)
	
		median <- (x2-x1)/(y2-y1) * (0.5 - y2) + x1 
	
	return(median)
	
}



#---------------------------------------------------------------------------------------
latexSection <- function(sctn){
	latexTex <- paste(
		"\\clearpage \n",
		"\\section{", sctn, "} \n \n",
		 sep="")
	return(latexTex)
}


#---------------------------------------------------------------------------------------
latexFigure <- function(figName, cptn=NULL,label=NULL, scale=0.5){
	latexTex <- paste(
		"\\begin{figure}[ht] \n",
		"\\centering \n", 
		"\\includegraphics[scale=" , scale, " ,draft=FALSE]{", figName, "} \n",  
		"\\caption{" , cptn, "} \n",
		"\\label{", label, "} \n", 
		"\\end{figure} \n", 
 		sep="") 

	 return(latexTex)
}


#--------------------------------------------------------------------------------------
getMinRow <- function(thingToMin){
	xx <-thingToMin
	minRw <- seq(1, length(xx))[xx == min(xx)]
	return(minRw)
}
