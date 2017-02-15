
#----------------------------------------------------------------------------------------

weightedPDF <- function(x, weight.lnorm, weight.gamma, parms.lnorm, parms.gamma, parms.pareto){
	
	weight.pareto <- 1- weight.lnorm - weight.gamma
	
	f.gamma_x <- dgamma(x,shape= parms.gamma$shape, scale=parms.gamma$scale )
	f.lnorm_x <- dlnorm(x, meanlog=parms.lnorm$meanlog, sdlog=parms.lnorm$sdlog)
	f.pareto_x <- dpareto(x, alpha=parms.pareto$alpha, beta=parms.pareto$beta)
	f.pareto_x[x < parms.pareto$beta ]<-0

	f.weighted_x <- weight.lnorm * f.lnorm_x +
					weight.gamma * f.gamma_x + 
					weight.pareto * f.pareto_x 

	return(f.weighted_x)

	
}


#----------------------------------------------------------------------------------------

getTriangleArray <- function(N){
	weights <- array(dim=c((N+1)*(N+2)/2,3),dimnames=c("w1","w2", "w3"))
	row <- 0
	for(i in 0:N){
		for(j in 0:i){
			row <- row + 1
			weights[row,1] <- w1 <- i/N - j/N
			weights[row,2] <- w2 <- j/N
			weights[row,3] <- 1 - w1 - w2
		
		}	# next j
	}  # next i

	return(weights)

}


#----------------------------------------------------------------------------------------

createForecastFromDistn2.pareto <- function(x, parms, info){
	
	
	x1 <- x[x<parms$beta]
	x2 <- x[x>= parms$beta]
	
	y1 <- rep(0, length(x1))
	y2 <- dpareto(x2, alpha = parms$alpha, beta=parms$beta)
	return(list(x=x, y=c(y1,y2)))
	
}

#----------------------------------------------------------------------------------------

createForecastFromDistn2.paretoE <- function(x, parms, info){
	# this extends the pareto distribution to zero by adding a lower right tail below beta - and then renormalising..
	
	
	pwr <- info$pwr
	fracOfPeak <- info$fracOfPeak
	
	x1 <- x[x<=parms$beta]
	x2 <- x[x> parms$beta]
	
	
	y2 <- dpareto(x2, alpha = parms$alpha, beta=parms$beta)
	
	f <- approxfun(y2,x2)

	maxY2 <- max(y2)

	k <- (maxY2/ fracOfPeak)/(f(maxY2)^pwr )

	y1 <- k * x1^pwr

	x <- c(x1,x2)
	y <- c(y1,y2)
	
	Ival <- simpleIntegrate(list(x=x, y=y))

	y <- y / Ival


	return(list(x=x, y=y))
	
}

#----------------------------------------------------------------------------------------
createForecastFromDistn2.lnorm <- function(x,parms,  info){
	
	y <- dlnorm(x= x, meanlog = parms$meanlog, sdlog = parms$sdlog)
	return(list(x=x, y=y))
	
}

#----------------------------------------------------------------------------------------

createForecastFromDistn2.gamma <- function(x, parms,  info){
	
	y <- dgamma(x, shape = parms$shape, scale = parms$scale)
	return(list(x=x, y=y))
	
}

#----------------------------------------------------------------------------------------
createForecastFromDistn2 <- function(x, parms, info){
	
	# verification is supplies so that the forecast width is always large enough to give a non-zero probability of the verification in that year
	distn <- info$distn
	if(distn == "lnorm")  forecast <- createForecastFromDistn2.lnorm(x,parms, info)
	if(distn == "gamma")  forecast <- createForecastFromDistn2.gamma(x,parms, info)
	if(distn == "pareto")  forecast <- createForecastFromDistn2.pareto(x,parms, info)
	if(distn == "paretoE")  forecast <- createForecastFromDistn2.paretoE(x,parms, info)
	return(forecast)
	
}


#----------------------------------------------------------------------------------------
plotTriangleWithDot <- function(distn, scoreType, scoreParam, varVal, w1,w2){	
	distCol <- rep("grey",3)
	if(distn == "gamma") distCol[3] <- "blue"
	if(distn == "lognormal") distCol[2] <- "blue"
	if(distn == "hybridPareto") distCol[1] <- "blue"

	vtx.pareto <- c(0,0)
	vtx.hybridPareto <- c(0.025, 0.025)
	vtx.lnorm <- c(1,0)
	vtx.gamma <- c(0,1)

	pts <- rbind(vtx.pareto,vtx.lnorm,vtx.gamma, vtx.pareto)
	pts2 <-rbind(vtx.hybridPareto,vtx.lnorm,vtx.gamma) 
	par(mar=c(1,0,0,0))

	plot(pts[,1], pts[,2], type="l", axes=FALSE, ann=FALSE)
	
	points(pts2[,1], pts2[,2], pch=19, col="black", cex=0.5)
	points(w1, w2, pch=19, col="red")
	text(c(-0.01,1,-.02), c(-0.025,-.05,0.95), labels=c("hP", "L", "G"), cex=1, pos=c(4,3,4), col=distCol)
	text(x=0.75,y=0.75, labels=scoreType, cex=0.5)
	if(scoreParam < 9999) text(x=0.75,y=0.65, labels=scoreParam, cex=0.5)
	text(x=0.75,y=0.55, labels=paste("var = ", varVal, sep=""), cex=0.5)

}	


#----------------------------------------------------------------------------------------

simpleCDF <- function(func){
	# assumes a function y=f(x) has been discretised into x and y values
	# assumes the function is actually continuous and the y values are equal to
	#  the function f evaluated at the x values.
	
	x <- func$x
	y <- func$y
	lX <- length(x)
	
	#if(lX==1) {
#		int <- 0
#	}
	
	#if(lX==2){
#		int1 <-  diff(x)[1] * y[1]
#		int2 <-  diff(x)[1] * y[2]
#		int <- (int1+int2)/2	
#	}
	
	if(lX>2){
		int1 <- cumsum(diff(x) * y[1:(lX-1)])
		int2 <- cumsum(diff(x) * y[2:lX])
		int <- (int1+int2)/2
	}
	
	return(int)

}



#----------------------------------------------------------------------------------------




getForecastPDF_weighted <- function(xVals, mean, var, weight.lnorm, weight.gamma){
	parms.lnorm <- getParamsFromMeanAndVar(distn="lnorm", mean=mean, var=var)
	parms.gamma <- getParamsFromMeanAndVar(distn="gamma", mean=mean, var=var)
	parms.pareto <- getParamsFromMeanAndVar(distn="pareto", mean=mean, var=var)

	fcst.wghtd <- list(
		x=xVals, 
		y=weightedPDF(x=xVals,
			weight.lnorm = weight.lnorm,
			weight.gamma = weight.gamma,
			parms.lnorm = parms.lnorm,
			parms.gamma = parms.gamma,
			parms.pareto =parms.pareto ))

	return(fcst.wghtd)
	
}



#----------------------------------------------------------------------------------------


getGlidePath <- function(truth, forecast, scoreType, scoreParam, glideResults){
	
	gR <- glideResults
	gR <- gR[gR$trueMean == truth$mean &
				gR$trueVar == truth$var &
				gR$trueWeightL == truth$weightL &
				gR$trueWeightG == truth$weightG &
				gR$fcstMean == forecast$mean &
				gR$fcstVar == forecast$var &
				gR$fcstWeightL == forecast$weightL &
				gR$fcstWeightG == forecast$weightG &
				gR$scoreType == scoreType &
				gR$scoreParam == scoreParam  , ]

	return(gR)
		
}




#----------------------------------------------------------------------------------------

plotGlidePathOld <- function(truth,forecast, scoreType, scoreParam, glideResults, quantiles = c(0.1,0.25,0.50, 0.75, 0.9), keyQuantile=quantiles(length(quantiles)), plotKeyQuantile=FALSE, diffTolerance=0.000000001, ... ){
	
	if(!(keyQuantile %in% quantiles) )  quantiles <- sort(c(quantiles, keyQuantile))
	if(!(round(1-keyQuantile,digits=3) %in% quantiles) )  quantiles <- sort(c(quantiles, round(1- keyQuantile,digits=3)))
	
	
	glide <- getGlidePathPair(truth=truth, forecast=forecast, scoreType=scoreType, scoreParam = scoreParam, glideResults = glideResults)
	
	glide.ref <- getGlidePathPair(truth=truth, forecast=truth, scoreType=scoreType, scoreParam = scoreParam, glideResults = glideResults)	
	
	timeQls <- NULL
	timeQls.ref <- NULL
	for (fT in sort(unique(glide$fcstTime))){
		
		glide.time <- glide[glide$fcstTime == fT, ]
		glide.ref.time <- glide.ref[glide.ref$fcstTime == fT, ]
		
		qls <- quantile(glide.time$informationDeficit, quantiles)
		qls.ref <- quantile(glide.ref.time$informationDeficit, quantiles)
		
		qlsRow <- c(time=fT, qls)
		qlsRow.ref <- c(time=fT, qls.ref)
		
		timeQls <- rbind(timeQls, qlsRow )
		timeQls.ref <- rbind(timeQls.ref, qlsRow.ref )  
		
	}
	
	yRange <- range(rbind(timeQls[,-1], timeQls.ref[,-1] ))
	xRange <- range(rbind(timeQls[,1],timeQls.ref[,1] ))
	plot(x=xRange, y=yRange, type="n",...)
	abline(h=0, col="grey")
	
	for (j in seq(along=quantiles)){
		labPos <- 0.9
		lX <- length(timeQls[,1])
		posX <- floor(labPos * lX)
		lines(x=timeQls[,1], y=timeQls[,1+j])
		text(timeQls[posX,1], timeQls[posX,1+j], labels=quantiles[j]*100, col="black", cex=0.5 )

		labPos <- 0.1
		lX <- length(timeQls.ref[,1])
		posX <- floor(labPos * lX)
		lines(x=timeQls.ref[,1], y=timeQls.ref[,1+j], col="darkgrey")
		text(timeQls.ref[posX,1], timeQls.ref[posX,1+j], labels=quantiles[j]*100, col="darkgrey", cex=0.5 )
	}
	
	# now find when the key quantile curves cross
	if(plotKeyQuantile){
		keyQl.ref <- timeQls.ref[, c(FALSE, (quantiles == keyQuantile))]
		keyQl <- timeQls[, c(FALSE, (
			abs(1-keyQuantile - quantiles)< diffTolerance))]
		
		diffQl <- keyQl - keyQl.ref
		
		lnQ <- length(diffQl)
		i <- 1
		
		repeat {
			if((all(diffQl[i:lnQ]>0) | i==lnQ)) break
			i <- i + 1
			
		}
		# need the time when afterwards the difference is always negative (or positive?)
	
		timeCrossing <- i  
		if(timeCrossing != lnQ ){
			abline(v=timeCrossing, col="red")
			points(x=timeCrossing, y=keyQl[timeCrossing], pch=19, col="red", cex=0.75)
			mtext(text=timeCrossing, at=timeCrossing, side=1, line=0, col="red", cex=0.75)
		} else {
			timeCrossing <- 999
		}
	}
	
	title(main=paste("score:", scoreType, ifelse(scoreParam==9999,"", scoreParam), sep=" "))


	col1 <- 1
	col2 <- floor(45/256 * lnQ)
	col3 <- floor(75/256 * lnQ)

	mtext(side=1, line = 2, at= col1, "mean:", cex=0.5)
	mtext(side=1, line = 2.5, at= col1, "var:", cex=0.5)
	mtext(side=1, line = 3, at= col1, "weightL:", cex=0.5)	
	mtext(side=1, line = 3.5, at= col1, "weightG:", cex=0.5)
	mtext(side=1, line = 1.5, at= col2, "truth", cex=0.5)
	mtext(side=1, line = 1.5, at= col3, "forecast", cex=0.5)

	mtext(side=1, line = 2, at= col2, truth$mean, cex=0.5)
	mtext(side=1, line = 2.5, at= col2, truth$var, cex=0.5)
	mtext(side=1, line = 3, at= col2, truth$weightL, cex=0.5)
	mtext(side=1, line = 3.5, at= col2, truth$weightG, cex=0.5)

	mtext(side=1, line = 2, at= col3, forecast$mean, cex=0.5)
	mtext(side=1, line = 2.5, at= col3, forecast$var, cex=0.5)
	mtext(side=1, line = 3, at= col3, forecast$weightL, cex=0.5)
	mtext(side=1, line = 3.5, at= col3, forecast$weightG, cex=0.5)


		
	return(list(timeQls=timeQls, timeQls.ref=timeQls.ref, keyQuantile=keyQuantile, timeCrossing=timeCrossing))
}





#----------------------------------------------------------------------------------------

getGlideResultsOLD <- function(compParms, trueParms, scoreType,scoreParam, xVals, numGlides = 2^7, numVerifications = 2^7 ,start.seed=123456, len=2^13){
	
	names(compParms)<- c("fcstMean","fcstVar","fcstWeightL", "fcstWeightG")
	names(trueParms)<- c("trueMean","trueVar","trueWeightL", "trueWeightG")

	N <- numVerifications
	yrsTol <- NULL
	glideResults <- NULL
	for (i in 1:numGlides){
		print(i)
		seedVal <- start.seed + i * N
		companyForecast <- getForecastPDF_weighted(
			xVals, mean=compParms[1], 
			var=compParms[2], 
			weight.lnorm=compParms[3], weight.gamma=compParms[4])

		trueDist <- getForecastPDF_weighted(
			xVals, 
			mean=trueParms[1], 
			var=trueParms[2], weight.lnorm=trueParms[3], weight.gamma=trueParms[4])

		invCDF_true <- approxfun( 
			x=c(0,simpleCDF(list(x=trueDist$x, y=trueDist$y)),1), 
			y=c(trueDist$x[1],trueDist$x[1:(length(trueDist$x)-1)],max(xVals)))

		set.seed(seedVal)
		U<-runif(N)
		verifications <- invCDF_true(U)

		companyScores <- calcScore(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			verification = verifications, 
			scoreParams= scoreParam)

		scoreStack <- list()
		for (y in 1:N){
			scoreStack[[y]] <- companyScores[1:y]	
		}

	
		meanStack <- lapply(scoreStack, mean)
		meanStack.vec <- do.call(rbind, meanStack)

		eScore <- calcExpectedScore2(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			underlyingPDF = companyForecast, 
			scoreParams=scoreParam)


		informationDeficit <- meanStack.vec - eScore
	
		glides <- data.frame(
			cbind(t(trueParms), t(compParms)),
			seedVal = seedVal,
			scoreType = scoreType,
			scoreParam = ifelse(!is.null(scoreParams$alpha),scoreParams$alpha, 9999),
			len=len,
			rangeX0 = range(xVals)[1],
			rangeX1 = range(xVals)[2], 
			fcstTime = 1:N,
			informationDeficit=informationDeficit
			)
		
		glideResults <- rbind(glideResults, glides)	
	
	} # next seed

	return(glideResults)
}

#----------------------------------------------------------------------------------------

extractGlideQuantiles <- function(truth,forecast, scoreType, scoreParam, glideResults, quantiles = c(0.1,0.25,0.50, 0.75, 0.9)){
		
	glide <- getGlidePath(truth=truth, forecast=forecast, scoreType=scoreType, scoreParam = scoreParam, glideResults = glideResults)
	
	timeQls <- NULL
	for (fT in sort(unique(glide$fcstTime))){
		glide.time <- glide[glide$fcstTime == fT, ]
		qls <- quantile(glide.time$informationDeficit, quantiles)
		qlsRow <- c(time=fT, qls)
		timeQls <- rbind(timeQls, qlsRow )
	}
	
	colnames(timeQls) <- c("time", paste("Q", quantiles*100, sep=""))	
	return(list(timeQls=timeQls, truth=truth, forecast=forecast, scoreType = scoreType, scoreParam=scoreParam))
}


#----------------------------------------------------------------------------------------
getResIndices <- function(xVals, len, scoreType, scoreParam, trueParms, compParms,numVerifications, seedVal){
		
		
		resIndices <- data.frame(
			cbind(t(trueParms), t(compParms)),
			seedVal = seedVal,
			numVerifications = numVerifications,
			scoreType = scoreType,
			scoreParam = ifelse(!is.null(scoreParams$alpha),scoreParams$alpha, 9999),
			len=len,
			rangeX0 = range(xVals)[1],
			rangeX1 = range(xVals)[2]
			)
	
		return(resIndices)
}

#----------------------------------------------------------------------------------------

getResBlock <- function(resIndices,gQ){
	
	resBlock <- data.frame(resIndices, gQ$timeQls)
	return(resBlock)
		
}


#----------------------------------------------------------------------------------------

checkIfExperimentDoneBefore <- function(indexVals , resTabName = "triangleInfDefctQuantiles", dsn){
	
	iV <- indexVals
	require(RODBC)
	chnl <- odbcConnect(dsn=dsn)
    tables <- sqlTables(channel=chnl, catalog=dsn)  #seems to need the catalog parameter otherwise only sees the first table
	
	resTableExists <- resTabName %in% tables$TABLE_NAME	
	
	if (resTableExists){

			tableExists <- TRUE
			
			# run the sql query from R - dont pull the data in (wasteful)
			sqlQry <- paste(
				"select sum(trueMean) from ", resTabName, 
				" where ", 
					" trueMean = ", iV$trueMean," and ",
					" trueVar = ", iV$trueVar," and ",
					" trueWeightL = ", iV$trueWeightL," and ",					
					" trueWeightG = ", iV$trueWeightG," and ",					
					" fcstMean = ", iV$fcstMean," and ",
					" fcstVar = ", iV$fcstVar," and ",
					" fcstWeightL = ", iV$fcstWeightL," and ",					
					" fcstWeightG = ", iV$fcstWeightG," and ",
					" numVerifications = ", iV$numVerifications, " and ",
					" seedVal = ", iV$seedVal," and ",					
					" scoreType = '", iV$scoreType,"' and ",
					" scoreParam = ", iV$scoreParam," and ",
					" rangeX0 = ", iV$rangeX0, " and ",
					" rangeX1 = ", iV$rangeX1,											
				" ;", sep="")
			qryResult <- sqlQuery(channel = chnl, query = sqlQry)
			
			
			# if result is NA then the experiment has NOT been done before
			doneBefore <- as.logical(!is.na(qryResult))   # to do
			

		} else {

			tableExists <- FALSE
			doneBefore <- FALSE

		}
	
	odbcClose(chnl)
	return(list(tableExists=tableExists, doneBefore = doneBefore))
	
}


#----------------------------------------------------------------------------------------


triangle_saveResults <- function(resData , resTabName = "triangleInfDefctQuantiles",experimentDoneBefore, dsn="modOffScores"){
	
	require(RODBC)
	chnl <- odbcConnect(dsn=dsn)

	if(experimentDoneBefore$doneBefore){
		print("warning: experiment done before this instance not saved")
		} else {
			if (experimentDoneBefore$tableExists){
				# append to table
				sqlSave(channel=chnl, dat=resData, 
					tablename=resTabName, rownames=FALSE, append=TRUE)    
			} else {
				sqlSave(channel=chnl, dat=resData, 
					tablename=resTabName, rownames=FALSE)    
			}
	}
	odbcClose(chnl)
}

#----------------------------------------------------------------------------------------

plotGlidePath <- function(truth,forecast, scoreInfo, otherInfo, keyQuantile, plotKeyQuantile=TRUE, dsn="modOffScores", resTabName = "triangleInfDefctQuantiles"){
	
	require(Hmisc)
	scoreParam <- ifelse(is.null(scoreInfo$scoreParam$alpha), 9999,scoreInfo$scoreParam$alpha) 
	
	
	qDatForecast <- getQuantileData(
		truth=list(mean=truth$mean, var=truth$var, weightL=truth$weightL,
			weightG=truth$weightG), 
		forecast=list(mean=forecast$mean, var=forecast$var, weightL=forecast$weightL, 				weightG=forecast$weightG), 
		scoreInfo = scoreInfo,
		otherInfo = otherInfo,
		dsn=dsn,
		resTabName=resTabName
		)
	
	
	qDatTruth <- getQuantileData(
		truth=list(mean=truth$mean, var=truth$var, weightL=truth$weightL,
			weightG=truth$weightG), 
		forecast=list(mean=truth$mean, var=truth$var, weightL=truth$weightL,
			weightG=truth$weightG), 
		scoreInfo = scoreInfo,
		otherInfo = otherInfo,
		dsn=dsn,
		resTabName=resTabName
		)
	
	timeCrossing <- getTimeCrossing(qFcstData = qDatForecast, qTruthData = qDatTruth, keyQuantile=keyQuantile )	
	
	quantiles <- getQuantileList(qDatTruth)/100
	
	if(!(keyQuantile %in% quantiles) )  stop("upper key quantile not available")
	if(!(round(1-keyQuantile,digits=3) %in% quantiles) )  stop("low key quantile not available")
	
	timeQls <- qDatForecast[, c("time", paste("Q", quantiles*100, sep=""))]
	timeQls.ref <- qDatTruth[, c("time", paste("Q", quantiles*100, sep=""))]
		
	yRange <- range(rbind(timeQls[,-1], timeQls.ref[,-1] ))
	xRange <- range(rbind(timeQls[,1],timeQls.ref[,1] ))
	plot( x=xRange, y=yRange, type="n", xlab="time", ylab="informationDeficit" )  
	abline(h=0, col="grey")
	
	for (j in seq(along=quantiles)){
		labPos <- 0.9
		lX <- length(timeQls[,1])
		posX <- floor(labPos * lX)
		lines(x=timeQls[,1], y=timeQls[,1+j])
		text(timeQls[posX,1], timeQls[posX,1+j], labels=quantiles[j]*100, col="black", cex=0.5 )

		labPos <- 0.1
		lX <- length(timeQls.ref[,1])
		posX <- floor(labPos * lX)
		lines(x=timeQls.ref[,1], y=timeQls.ref[,1+j], col="darkgrey")
		text(timeQls.ref[posX,1], timeQls.ref[posX,1+j], labels=quantiles[j]*100, col="darkgrey", cex=0.5 )
	}
	
	# now find when the key quantile curves cross
	if(plotKeyQuantile){
		keyQl.ref <- timeQls.ref[, c(FALSE, (quantiles == keyQuantile))]
		lnQ <- length(timeQls$time)
		if(timeCrossing != lnQ ){
			abline(v=timeCrossing, col="red")
			points(x=timeCrossing, y=keyQl.ref[timeCrossing], pch=19, col="red", cex=0.75)
			mtext(text=timeCrossing, at=timeCrossing, side=1, line=0, col="red", cex=0.75)
		} 
	}
	
	title(main=paste("score:", scoreInfo$scoreType, ifelse(scoreParam==9999,"", scoreParam), sep=" "))


	col1 <- 1
	col2 <- floor(45/256 * lnQ)
	col3 <- floor(75/256 * lnQ)

	mtext(side=1, line = 2, at= col1, "mean:", cex=0.5)
	mtext(side=1, line = 2.5, at= col1, "var:", cex=0.5)
	mtext(side=1, line = 3, at= col1, "weightL:", cex=0.5)	
	mtext(side=1, line = 3.5, at= col1, "weightG:", cex=0.5)
	mtext(side=1, line = 1.5, at= col2, "truth", cex=0.5)
	mtext(side=1, line = 1.5, at= col3, "forecast", cex=0.5)

	mtext(side=1, line = 2, at= col2, formatC(truth$mean, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 2.5, at= col2, formatC(truth$var, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3, at= col2, formatC(truth$weightL, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3.5, at= col2, formatC(truth$weightG, format="f", digits=3), cex=0.5)

	mtext(side=1, line = 2, at= col3, formatC(forecast$mean, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 2.5, at= col3, formatC(forecast$var, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3, at= col3, formatC(forecast$weightL, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3.5, at= col3, formatC(forecast$weightG, format="f", digits=3), cex=0.5)
	
	gC <- par("usr")  #for grid Coordinates
	xL <- (gC[2]-gC[1]) * 0.8 + gC[1]
	xR <- (gC[2]-gC[1]) * 0.95 + gC[1]
	yB <- (gC[4]-gC[3]) * 0.8 + gC[3]
	yT <- (gC[4]-gC[3]) * 0.95 + gC[3]

	subplot(plotTriangleWithTwoDots(forecast$weightL, forecast$weightG, truth$weightL, truth$weightG), x=c(xL,xR), y=c(yB,yT), par=list(mar=c(1,0,0,0)))

	
	
}


#----------------------------------------------------------------------------------------

getTimeCrossing <- function(qFcstData, qTruthData,  keyQuantile){
	
	# this assumes that qData has already been subsetted into two tables one containing the truth,truth pair quantiles and one containing the truth,forecast pair quantiles 
	
	quantiles <- getQuantileList(qData = qFcstData)/100

		
	if(!(keyQuantile %in% quantiles) )  stop("upper key quantile not available")
	if(!(round(1-keyQuantile,digits=3) %in% quantiles) )  stop("low key quantile not available")
	
		keyQl.ref <- qTruthData[, paste("Q", keyQuantile*100, sep="")]
		keyQl <- qFcstData[, paste("Q",(1-keyQuantile)*100, sep="")]
		
		diffQl <- keyQl - keyQl.ref
		
		lnQ <- length(diffQl)
		i <- 1
		
		repeat {
			if((all(diffQl[i:lnQ]>0) | i==lnQ)) break
			i <- i + 1
			
		}	
		
		if (i ==  lnQ) i <- 10^(ceiling(log(lnQ, base=10))+1) -1
		
		return(i)
	
	
}


#----------------------------------------------------------------------------------------

getQuantileList <- function(qData){

	findQ <- apply(as.array(names(qData)), 1, strsplit, "Q")

	getC2 <- function(x){
		return(x[[1]][2])
	
	}

	quantiles <- do.call(rbind, lapply(findQ, getC2))
	quantiles <- as.numeric(quantiles[!is.na(quantiles)])

	return(quantiles)

}

#----------------------------------------------------------------------------------------

getQuantileDataOLD <- function(truth, forecast,scoreInfo, otherInfo, dsn="modOffScores", resTabName ="triangleInfDefctQuantiles" ){
	
	require(RODBC)
	chnl <- odbcConnect(dsn)
	
	
	scoreParam <- ifelse(is.null(scoreInfo$scoreParam$alpha), 9999,scoreInfo$scoreParam$alpha) 
	sqlQry <- paste(
				"select * from ", resTabName, 
				" where ", 
					" trueMean = ", truth$mean," and ",
					" trueVar = ", truth$var," and ",
					" trueWeightL = ", truth$weightL," and ",					
					" trueWeightG = ", truth$weightG," and ",					
					" fcstMean = ", forecast$mean," and ",
					" fcstVar = ", forecast$var," and ",
					" fcstWeightL = ", forecast$weightL," and ",					
					" fcstWeightG = ", forecast$weightG," and ",
					" numVerifications = ", otherInfo$numVerifications, " and ",
					" len = ", otherInfo$len, " and ",
					" seedVal = ", otherInfo$seedVal," and ",					
					" scoreType = '", scoreInfo$scoreType,"' and ",
					" scoreParam = ", scoreParam," and ",
					" rangeX0 = ", otherInfo$rangeX0, " and ",
					" rangeX1 = ", otherInfo$rangeX1,											
				" ;", sep="")
	qryResult <- sqlQuery(channel = chnl, query = sqlQry)
		
	odbcClose(chnl)	
	return(qryResult)
	
}


#----------------------------------------------------------------------------------------

plotTriangleWithTwoDots <- function( w1,w2, wt1,wt2){	

	vtx.pareto <- c(0,0)
	vtx.lnorm <- c(1,0)
	vtx.gamma <- c(0,1)

	pts <- rbind(vtx.pareto,vtx.lnorm,vtx.gamma, vtx.pareto)
	

	plot(pts[,1], pts[,2], type="l", axes=FALSE, ann=FALSE)
	points(w1, w2, pch=19, col="blue")
	points(wt1, wt2, pch=19, col="red")

}


#----------------------------------------------------------------------------------------

getGlideResultsOLD2 <- function(compParms, trueParms, scoreType,scoreParam, xVals, numGlides = 2^7, numVerifications = 2^7 ,start.seed=123456, len=2^13){
	
	names(compParms)<- c("fcstMean","fcstVar","fcstWeightL", "fcstWeightG")
	names(trueParms)<- c("trueMean","trueVar","trueWeightL", "trueWeightG")

	companyForecast <- getForecastPDF_weighted(
			xVals, mean=compParms[1], 
			var=compParms[2], 
			weight.lnorm=compParms[3], weight.gamma=compParms[4])

	trueDist <- getForecastPDF_weighted(
			xVals, 
			mean=trueParms[1], 
			var=trueParms[2], weight.lnorm=trueParms[3], weight.gamma=trueParms[4])

	invCDF_true <- approxfun( 
			x=c(0,simpleCDF(list(x=trueDist$x, y=trueDist$y)),1), 
			y=c(trueDist$x[1],trueDist$x[1:(length(trueDist$x)-1)],max(xVals)))



	N <- numVerifications
	yrsTol <- NULL
	glideResults <- NULL
	for (i in 1:numGlides){
		print(i)
		seedVal <- start.seed + i * N
		
		set.seed(seedVal)
		U<-runif(N)
		verifications <- invCDF_true(U)

		companyScores <- calcScore(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			verification = verifications, 
			scoreParams= scoreParam)

		scoreStack <- list()
		for (y in 1:N){
			scoreStack[[y]] <- companyScores[1:y]	
		}

	
		meanStack <- lapply(scoreStack, mean)
		meanStack.vec <- do.call(rbind, meanStack)

		eScore <- calcExpectedScore2(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			underlyingPDF = companyForecast, 
			scoreParams=scoreParam)


		informationDeficit <- meanStack.vec - eScore
	
		glides <- data.frame(
			cbind(t(trueParms), t(compParms)),
			seedVal = seedVal,
			scoreType = scoreType,
			scoreParam = ifelse(!is.null(scoreParams$alpha),scoreParams$alpha, 9999),
			len=len,
			rangeX0 = range(xVals)[1],
			rangeX1 = range(xVals)[2], 
			fcstTime = 1:N,
			informationDeficit=informationDeficit
			)
		
		glideResults <- rbind(glideResults, glides)	
	
	rm(scoreStack)
	} # next seed

	return(glideResults)
}


#------------------------------------------------------------------------------------

getGlideResultsOLD3 <- function(compParms, trueParms, scoreType,scoreParam, xVals, numGlides = 2^7, numVerifications = 2^7 ,start.seed=123456, len=2^13){
	
	names(compParms)<- c("fcstMean","fcstVar","fcstWeightL", "fcstWeightG")
	names(trueParms)<- c("trueMean","trueVar","trueWeightL", "trueWeightG")

	
	# This is the distribution that will be used to calculate the expected score
	companyForecast <- getForecastPDF_weighted(
			xVals, mean=compParms[1], 
			var=compParms[2], 
			weight.lnorm=compParms[3], weight.gamma=compParms[4])

	# This is the distribution that will be used to generate verification points
	trueDist <- getForecastPDF_weighted(
			xVals, 
			mean=trueParms[1], 
			var=trueParms[2], weight.lnorm=trueParms[3], weight.gamma=trueParms[4])

	invCDF_true <- approxfun( 
			x=c(0,simpleCDF(list(x=trueDist$x, y=trueDist$y)),1), 
			y=c(trueDist$x[1],trueDist$x[1:(length(trueDist$x)-1)],max(xVals)))



	N <- numVerifications
	yrsTol <- NULL
	glideResults <- list()
	for (i in 1:numGlides){
		print(i)
		seedVal <- start.seed + i * N
		
		set.seed(seedVal)
		U<-runif(N)
		verifications <- invCDF_true(U)

		companyScores <- calcScore(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			verification = verifications, 
			scoreParams= scoreParam)

		scoreStack <- list()
		for (y in 1:N){
			scoreStack[[y]] <- companyScores[1:y]	
		}

	
		meanStack <- lapply(scoreStack, mean)
		meanStack.vec <- do.call(rbind, meanStack)

		eScore <- calcExpectedScore2(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			underlyingPDF = companyForecast, 
			scoreParams=scoreParam)


		informationDeficit <- meanStack.vec - eScore
	
		glides <- data.frame(
			cbind(t(trueParms), t(compParms)),
			seedVal = seedVal,
			scoreType = scoreType,
			scoreParam = ifelse(!is.null(scoreParams$alpha),scoreParams$alpha, 9999),
			len=len,
			rangeX0 = range(xVals)[1],
			rangeX1 = range(xVals)[2], 
			fcstTime = 1:N,
			informationDeficit=informationDeficit
			)
		
		glideResults[[i]] <- glides
	
	rm(scoreStack)
	} # next seed


	glideResults <- do.call(rbind, glideResults)
	return(glideResults)
}


#----------------------------------------------------------------------------------------

getQuantileDataOLD2 <- function(truth, forecast,scoreInfo, otherInfo, dsn="modOffScores", resTabName ="triangleInfDefctQuantiles", tolerancePC = 0.00001 ){
	
	require(RODBC)
	chnl <- odbcConnect(dsn)
	
	
	scoreParam <- ifelse(is.null(scoreInfo$scoreParam$alpha), 9999,scoreInfo$scoreParam$alpha) 
	sqlQry <- paste(
				"select * from ", resTabName, 
				" where ", 
					" trueMean < ", truth$mean * (1+tolerancePC)," and ",
					" trueMean > ", truth$mean / (1+tolerancePC)," and ",
					" trueVar < ", truth$var * (1+tolerancePC)," and ",
					" trueVar > ", truth$var / (1+tolerancePC)," and ",
					" trueWeightL < ", truth$weightL * (1+tolerancePC)," and ",								" trueWeightL > ", truth$weightL / (1+tolerancePC)," and ",			
					" trueWeightG < ", truth$weightG * (1+tolerancePC)," and ",								" trueWeightG > ", truth$weightG / (1+tolerancePC)," and ",			
					" fcstMean < ", forecast$mean * (1+tolerancePC)," and ",
					" fcstMean > ", forecast$mean / (1+tolerancePC)," and ",
					" fcstVar < ", forecast$var * (1+tolerancePC)," and ",
					" fcstVar > ", forecast$var / (1+tolerancePC)," and ",
					" fcstWeightL < ", forecast$weightL * (1+tolerancePC)," and ",							" fcstWeightL > ", forecast$weightL / (1+tolerancePC)," and ",		
					" fcstWeightG < ", forecast$weightG * (1+tolerancePC)," and ",
					" fcstWeightG > ", forecast$weightG / (1+tolerancePC)," and ",
					" numVerifications = ", otherInfo$numVerifications, " and ",
					" len = ", otherInfo$len, " and ",
					" seedVal = ", otherInfo$seedVal," and ",					
					" scoreType = '", scoreInfo$scoreType,"' and ",
					" scoreParam = ", scoreParam," and ",
					" rangeX0 = ", otherInfo$rangeX0, " and ",
					" rangeX1 = ", otherInfo$rangeX1,											
				" ;", sep="")
	qryResult <- sqlQuery(channel = chnl, query = sqlQry)
		
	odbcClose(chnl)	
	return(qryResult)
	
}

#----------------------------------------------------------------------------------------

getSqlTol <- function(param, tolerancePC, paramValue){
	if (paramValue == 0){
		sqlString <- paste(param, " = 0", sep="")
	} else {
		sqlString <- paste(
			param, " < ", paramValue * (1 + tolerancePC), " and ", 
			param, " > ", paramValue / (1 + tolerancePC), 
		sep="")
	}
	return(sqlString)	
}

#----------------------------------------------------------------------------------------

getQuantileData <- function(truth, forecast,scoreInfo, otherInfo, dsn="modOffScores", resTabName ="triangleInfDefctQuantiles", tolerancePC = 0.00001 ){
	
	require(RODBC)
	chnl <- odbcConnect(dsn)
	
	
	scoreParam <- ifelse(is.null(scoreInfo$scoreParam$alpha), 9999,scoreInfo$scoreParam$alpha) 
	sqlQry <- paste(
				"select * from ", resTabName, 
				" where ", 
					getSqlTol("trueMean", tolerancePC, truth$mean ) ," and ",
					getSqlTol("trueVar", tolerancePC, truth$var ) ," and ",
					getSqlTol("trueWeightL", tolerancePC, truth$weightL ) ," and ",
					getSqlTol("trueWeightG", tolerancePC, truth$weightG ) ," and ",
					getSqlTol("fcstMean", tolerancePC, forecast$mean) ," and ",
					getSqlTol("fcstVar", tolerancePC, forecast$var) ," and ",
					getSqlTol("fcstWeightL", tolerancePC, forecast$weightL) ," and ",
					getSqlTol("fcstWeightG", tolerancePC, forecast$weightG) ," and ",
					" numVerifications = ", otherInfo$numVerifications, " and ",
					" len = ", otherInfo$len, " and ",
					" seedVal = ", otherInfo$seedVal," and ",					
					" scoreType = '", scoreInfo$scoreType,"' and ",
					" scoreParam = ", scoreParam," and ",
					" rangeX0 = ", otherInfo$rangeX0, " and ",
					" rangeX1 = ", otherInfo$rangeX1,											
				" ;", sep="")
	qryResult <- sqlQuery(channel = chnl, query = sqlQry)
		
	odbcClose(chnl)	
	return(qryResult)
	
}



#----------------------------------------------------------------------------------------

plotGlidePathFAB <- function(truth,forecast, scoreInfo, otherInfo, probObsOsExp, quantileExp, plotKeyQuantile=TRUE, dsn="modOffScores", resTabName = "triangleInfDefctQuantiles"){
	
	require(Hmisc)
	scoreParam <- ifelse(is.null(scoreInfo$scoreParam$alpha), 9999,scoreInfo$scoreParam$alpha) 

#qDatTruth = qObservedGlides, 
#qDatForecast = qExpectedGlides
	
	
	
	# Here the quantile is the expected glide path assuming that the forecast is correct
	qExpectedGlides <- getQuantileData(
		truth=list(mean=forecast$mean, var=forecast$var, weightL=forecast$weightL, 				weightG=forecast$weightG),  
		forecast=list(mean=forecast$mean, var=forecast$var, weightL=forecast$weightL, 				weightG=forecast$weightG), 
		scoreInfo = scoreInfo,
		otherInfo = otherInfo,
		dsn=dsn,
		resTabName=resTabName
		)
	
	# Here the quantile is the glide path when the truth is possibly different from the forecast but the forecast is as above
	 qObservedGlides <- getQuantileData(
		truth=list(mean=truth$mean, var=truth$var, weightL=truth$weightL,
			weightG=truth$weightG), 
		forecast=list(mean=forecast$mean, var=forecast$var, weightL=forecast$weightL, 				weightG=forecast$weightG),  
		scoreInfo = scoreInfo,
		otherInfo = otherInfo,
		dsn=dsn,
		resTabName=resTabName
		)
	
	timeCrossing <- getTimeCrossingFAB2(
		qObservedGlides = qObservedGlides , 
		qExpectedGlides = qExpectedGlides,  
		probObsOsExp = probObsOsExp, 
		quantileExp = quantileExp)	
	
	quantiles <- getQuantileList(qObservedGlides )/100
	
	if(!(round(probObsOsExp,5) %in% quantiles) )  stop("probObsOsExp not available")
	if(!(round(quantileExp,5) %in% quantiles) )  stop("quantileExp not available")
	
	
	timeQls.exp <- qExpectedGlides[, c("time", paste("Q", quantiles*100, sep=""))]
	timeQls.obs <- qObservedGlides[, c("time", paste("Q", quantiles*100, sep=""))]
		
	yRange <- range(rbind(timeQls.exp[,-1], timeQls.obs[,-1] ))
	xRange <- range(rbind(timeQls.exp[,1],timeQls.obs[,1] ))
	plot( x=xRange, y=yRange, type="n", xlab="time", ylab="informationDeficit" )  
	abline(h=0, col="grey")
	
	for (j in seq(along=quantiles)){
		labPos <- 0.9
		lX <- length(timeQls.exp[,1])
		posX <- floor(labPos * lX)
		lines(x=timeQls.exp[,1], y=timeQls.exp[,1+j])
		text(timeQls.exp[posX,1], timeQls.exp[posX,1+j], labels=quantiles[j]*100, col="black", cex=0.5 )

		labPos <- 0.1
		lX <- length(timeQls.obs[,1])
		posX <- floor(labPos * lX)
		lines(x=timeQls.obs[,1], y=timeQls.obs[,1+j], col="darkgrey")
		text(timeQls.obs[posX,1], timeQls.obs[posX,1+j], labels=quantiles[j]*100, col="darkgrey", cex=0.5 )
	}
	
	# now find when the key quantile curves cross
	if(plotKeyQuantile){
		
		qPlot <- quantileExp
		
		lDat <- length(qObservedGlides$Q50)
		if(qObservedGlides$Q50[lDat] < qExpectedGlides$Q50[lDat]){
			qPlot <- round(1 - qPlot, 5)
		}
		
		
		keyQl.exp <- timeQls.exp[, c(FALSE, (quantiles == qPlot))]
		lnQ <- length(timeQls.exp$time)
		
		if(timeCrossing != lnQ ){
			abline(v=timeCrossing, col="red")
			points(x=timeCrossing, y=keyQl.exp[timeCrossing], pch=4, col="red", cex=0.75)
			mtext(text=timeCrossing, at=timeCrossing, side=1, line=0, col="red", cex=0.75)
		} 
	}
	
	title(main=paste("score:", scoreInfo$scoreType, ifelse(scoreParam==9999,"", scoreParam), sep=" "))


	col1 <- 1
	col2 <- floor(45/256 * lnQ)
	col3 <- floor(75/256 * lnQ)

	mtext(side=1, line = 2, at= col1, "mean:", cex=0.5)
	mtext(side=1, line = 2.5, at= col1, "var:", cex=0.5)
	mtext(side=1, line = 3, at= col1, "weightL:", cex=0.5)	
	mtext(side=1, line = 3.5, at= col1, "weightG:", cex=0.5)
	mtext(side=1, line = 1.5, at= col2, "truth", cex=0.5)
	mtext(side=1, line = 1.5, at= col3, "forecast", cex=0.5)

	mtext(side=1, line = 2, at= col2, formatC(truth$mean, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 2.5, at= col2, formatC(truth$var, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3, at= col2, formatC(truth$weightL, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3.5, at= col2, formatC(truth$weightG, format="f", digits=3), cex=0.5)

	mtext(side=1, line = 2, at= col3, formatC(forecast$mean, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 2.5, at= col3, formatC(forecast$var, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3, at= col3, formatC(forecast$weightL, format="f", digits=3), cex=0.5)
	mtext(side=1, line = 3.5, at= col3, formatC(forecast$weightG, format="f", digits=3), cex=0.5)
	
	gC <- par("usr")  #for grid Coordinates
	xL <- (gC[2]-gC[1]) * 0.8 + gC[1]
	xR <- (gC[2]-gC[1]) * 0.95 + gC[1]
	yB <- (gC[4]-gC[3]) * 0.8 + gC[3]
	yT <- (gC[4]-gC[3]) * 0.95 + gC[3]

	subplot(plotTriangleWithTwoDots(forecast$weightL, forecast$weightG, truth$weightL, truth$weightG), x=c(xL,xR), y=c(yB,yT), par=list(mar=c(1,0,0,0)))

	
	
}





#----------------------------------------------------------------------------------------

getTimeCrossingFAB <- function(qFcstData, qTruthData,  keyQuantile){
	
	# this assumes that qData has already been subsetted into two tables one containing the truth,truth pair quantiles and one containing the truth,forecast pair quantiles 
	
	quantiles <- getQuantileList(qData = qFcstData)/100

		
	if(!(keyQuantile %in% quantiles) )  stop("upper key quantile not available")
	if(!(round(1-keyQuantile,digits=3) %in% quantiles) )  stop("low key quantile not available")
	
		keyQl.ref <- qTruthData[, paste("Q", (1-keyQuantile)*100, sep="")]
		keyQl <- qFcstData[, paste("Q",keyQuantile*100, sep="")]
		
		diffQl <- -(keyQl - keyQl.ref)
		
		lnQ <- length(diffQl)
		i <- 1
		
		repeat {
			if((all(diffQl[i:lnQ]>0) | i==lnQ)) break
			i <- i + 1
			
		}	
		
		if (i ==  lnQ) i <- 10^(ceiling(log(lnQ, base=10))+1) -1
		
		return(i)
	
	
}




#------------------------------------------------------------------------------------

getGlideResults <- function(compParms, trueParms, scoreType,scoreParam, xVals, numGlides = 2^7, numVerifications = 2^7 ,start.seed=123456, len=2^13){
	
	names(compParms)<- c("fcstMean","fcstVar","fcstWeightL", "fcstWeightG")
	names(trueParms)<- c("trueMean","trueVar","trueWeightL", "trueWeightG")

	
	# This is the distribution that will be used to calculate the expected score
	companyForecast <- getForecastPDF_weighted(
			xVals, mean=compParms[1], 
			var=compParms[2], 
			weight.lnorm=compParms[3], weight.gamma=compParms[4])

	# This is the distribution that will be used to generate verification points
	trueDist <- getForecastPDF_weighted(
			xVals, 
			mean=trueParms[1], 
			var=trueParms[2], weight.lnorm=trueParms[3], weight.gamma=trueParms[4])

	invCDF_true <- approxfun( 
			x=c(0,simpleCDF(list(x=trueDist$x, y=trueDist$y)),1), 
			y=c(trueDist$x[1],trueDist$x[1:(length(trueDist$x)-1)],max(xVals)))



	N <- numVerifications
	yrsTol <- NULL
	glideResults <- list()
	for (i in 1:numGlides){
		print(i)
		seedVal <- start.seed + i * N
		
		set.seed(seedVal)
		U<-runif(N)
		verifications <- invCDF_true(U)

		companyScores <- calcScore(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			verification = verifications, 
			scoreParams= scoreParam)

		#scoreStack <- list()
		#for (y in 1:N){
		#	scoreStack[[y]] <- companyScores[1:y]	
		#}

	
		#meanStack <- lapply(scoreStack, mean)
		#meanStack.vec <- do.call(rbind, meanStack)
		meanStack.vec <- cumsum(companyScores)/1:length(companyScores)
		
		eScore <- calcExpectedScore2(
			scoreType=scoreType, 
			forecastPDF = companyForecast, 
			underlyingPDF = companyForecast, 
			scoreParams=scoreParam)


		informationDeficit <- meanStack.vec - eScore
	
		glides <- data.frame(
			cbind(t(trueParms), t(compParms)),
			seedVal = seedVal,
			scoreType = scoreType,
			scoreParam = ifelse(!is.null(scoreParams$alpha),scoreParams$alpha, 9999),
			len=len,
			rangeX0 = range(xVals)[1],
			rangeX1 = range(xVals)[2], 
			fcstTime = 1:N,
			informationDeficit=informationDeficit
			)
		
		glideResults[[i]] <- glides
	
	rm(scoreStack)
	} # next seed


	glideResults <- do.call(rbind, glideResults)
	return(glideResults)
}


#----------------------------------------------------------------------------------------

getTimeCrossingFAB2 <- function(qObservedGlides, qExpectedGlides,  probObsOsExp, quantileExp ){

	# qFcst  qObservedGlides, 
	# qTruth qExpectedGlides
	# quantileExp - the glide quantile of the expected information deficit which determines the level of deficit D(t) to be tested against at each time t.   
	# probObsOsExp  - this is the probability that the observed information deficit is outside of the specified confidence interval (quantileExp) around the expected information deficit
	
	# The function finds the time t* when the observed information deficit has the specified probabilty (probObsLtExp) of being outside of the specified confidenc interval for all t greater than t*
	
	# NB observations arise from the true distribution and are used to calculate the score of the forecasted distribution.   Then the expected score from the forecast is deducted to calculate the information deficit.

	# this assumes that qData has already been subsetted into two tables one containing the truth,truth pair quantiles and one containing the truth,forecast pair quantiles 
	
	quantiles <- getQuantileList(qData = qObservedGlides)/100
	diffScalar <- 1 
	# qTruth is black  (it is expected information deficit assuming the foreacst to be true)
	# qTruth is big quantile for low qFcst quantile

# if keyQuantile = Y (grey)  and keyQuantile.ref = X  (black)
# define the time crossing as the time at which the observed information deficit (grey) has probability Y of being outside of the X quantile lines


	lDat <- length(qObservedGlides$Q50)
	if(qObservedGlides$Q50[lDat] < qExpectedGlides$Q50[lDat]){
		quantileExp <- 1 - quantileExp
		probObsOsExp <- 1 - probObsOsExp
		diffScalar <- -1  #This reverses the test below 
	}

		
	if(!(round(probObsOsExp,5) %in% quantiles) )  stop("probObsOsExp not available")
	if(!(round(quantileExp,5) %in% quantiles) )  stop("quantileExp not available")
	

		keyQl.exp <- qExpectedGlides[, paste("Q", quantileExp*100, sep="")]
		keyQl.obs <- qObservedGlides[, paste("Q",  (1-probObsOsExp)*100, sep="")]
		
		diffQl <- (keyQl.obs - keyQl.exp)  #if t crossing < lnQ this will be positive for diffScalar =1 and negative otherwise
		
		lnQ <- length(diffQl)
		i <- 1
		
		repeat {
			if((all(diffScalar * diffQl[i:lnQ] > 0) | i==lnQ)) break
			i <- i + 1
			
		}	
		
		if (i ==  lnQ) i <- 10^(ceiling(log(lnQ, base=10))+1) -1
		
		return(i)
	
	
}


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

createPcInBinsOneSeed <- function(experiment){

	
	trueDist <- getForecastPDF_weighted(
			xVals = experiment$xVals, 
			mean = experiment$truth$mean, 
			var = experiment$truth$var, 
			weight.lnorm = experiment$truth$weightL, 
			weight.gamma = experiment$truth$weightG)

	invCDF_true <- approxfun( 
			x=c(0,simpleCDF(list(x=trueDist$x, y=trueDist$y)),1), 
			y=c(trueDist$x[1],trueDist$x[1:(length(trueDist$x)-1)],max(experiment$xVals)))

	set.seed(experiment$seed)
	U<-runif(experiment$numObs)
	observations <- invCDF_true(U)

	
	pcInBins <- list()
	for (i in 1:length(observations)){
	pcInBins[[i]] <- hist(observations[1:i], breaks=experiment$breaks, plot=FALSE)$counts /  i
	#print(i)
	}
	pcInBins <- do.call(rbind, pcInBins)

	return(pcInBins)	
}	

#----------------------------------------------------------------------------

createPcInBins <- function(experiment, numBlocks){
	seed <- experiment$seed
	numCols <- length(experiment$breaks)+2
	numObs <- experiment$numObs
	resTab <- array(dim=c(numObs * numBlocks, numCols))
	for (i in 1:numBlocks){
		exp_i <- experiment
		exp_i$seed <- seed + i * numObs
		pcRes <- createPcInBinsOneSeed(exp_i)	
		block <- cbind(ID=experiment$ID, seed=exp_i$seed,obNum = 1:numObs, pcRes)
		resTab[seq(1+(i-1)*numObs, i*numObs),] <- block
		if(floor(i/2^4)*2^4 == i) print(i)
	}
	
	binNames <- paste("Bin", 1:(numCols-3), sep="")
	colnames(resTab) <- c("ID", "seed", "obNum", binNames)
	return(as.data.frame(resTab))
	
}


#----------------------------------------------------------------------------

checkIfIDExists <- function(ID, tableName, dsn){
	require(RODBC)
	chnl <- odbcConnect(dsn)
	tableExists <- TRUE  #assume this is true to start with
	IDExists <- TRUE # assume true to start with (so wont create duplicates)
	qry <- sqlQuery(
		channel=chnl, 
		query=paste("select * from ", tableName, " where ID = ", ID," ;", sep="")
		)

	if (is.null(dim(qry))){
		tableExists <- FALSE
		IDExists <- FALSE
	} else {
		if(dim(qry)[1] == 0){
			IDExists <- FALSE  
		}
	}

	odbcClose(chnl)
	return(list(tableExists = tableExists, IDExists = IDExists))
} 	


#----------------------------------------------------------------------------


saveBinExp.exp <-function(experiment, dsn="modOffScores", blockTabName="pcInRange_ExpDesign",  breaksTabName = "pcInRange_Breaks"){
	
	require(RODBC)
	chnl <- odbcConnect("modOffScores")
	tabNames <- sqlTables(channel=chnl, catalog=dsn)
	
		
	rX <- range(experiment$xVals)

	dBlock <- data.frame(experiment$ID, experiment$truth, t(rX), length(xVals), experiment$numObs, experiment$seed, experiment$binQuantilesTabName)
	names(dBlock) <- c(
		"ID",
		paste("trth", names(experiment$truth), sep=""),
		"xMin",
		"xMax",
		"xLen",
		"numObs",
		"seed",
		"binQuantilesTabName"
		)
	
	binBlock <- data.frame(ID=experiment$ID, experiment$breaks)
	
	# decide whether ID is new - if not then do not save
	IDisNew <- FALSE
	
	blockTableExists <-FALSE
	if(blockTabName %in% tabNames[,"TABLE_NAME"]) blockTableExists <- TRUE

	breakTableExists <- FALSE
	if(breaksTabName %in% tabNames[,"TABLE_NAME"]) breakTableExists <- TRUE

	if (blockTableExists & !breakTableExists) stop("something has gone wrong some experiments exist but there are no break points")

	if (!blockTableExists & breakTableExists) stop("something has gone wrong a break points table exists but there are no experiments")

	bothTabsExist <- blockTableExists & breakTableExists
	
	if(bothTabsExist){

		# check if ID used before
		rtnIDs <- sqlQuery(channel=chnl, query=paste("select * from  ",blockTabName, "  where ID =  ", experiment$ID, " ;", sep="" ))
		if(!is.null(rtnIDs)){
			if(dim(rtnIDs)[1]>0) {
				stop("ID has been used before -  not saved")
			} else {
				# the results tables exist - but the ID is unique so can append the new results to the table
				
				sqlSave(channel = chnl, dat=dBlock, tablename = blockTabName , append=TRUE, rownames=FALSE)


		sqlSave(channel = chnl, dat=binBlock, tablename = breaksTabName  , append=TRUE, rownames=FALSE)
				
			}
		}
				
		
		
	} else {
		# given previous checks - in this case neither tables exist
		# therefore it is obvious that the ID is new!:  need to create both tables
		sqlSave(channel = chnl, dat=dBlock, tablename = blockTabName , append=FALSE, rownames=FALSE)


		sqlSave(channel = chnl, dat=binBlock, tablename = breaksTabName  , append=FALSE, rownames=FALSE)

		

	}
	
	odbcClose(chnl)
	return(TRUE)
	
}


#----------------------------------------------------------------------------

saveBinExp.results <- function(experiment, binPCResults, binResultsTabName = "pcInRange_Res", dsn="modOffScores"){
		
	require(RODBC)
	chnl <- odbcConnect(dsn)
	sqlInfo <- checkIfIDExists(ID = experiment$ID, tableName = binResultsTabName, dsn=dsn)
	
	if(sqlInfo$IDExists) stop("ID already exists, either change ID number or experiment already run")
		
	binPCResults <- as.data.frame(binPCResults)	
	if(sqlInfo$tableExists){
		# append
		sqlSave(channel, chnl, dat=binPCResults, tableName=binResultsTabName, append=TRUE, rownames=FALSE)
	} else {
		#create table
		sqlSave(channel= chnl, dat=binPCResults, tablename=binResultsTabName, append=FALSE, rownames=FALSE)
	}
	odbcClose(chnl)
	
}


#----------------------------------------------------------------------------


saveBinExp <- function(experiment, binPCQuantiles,binQuantilesTabName = "pcInRange_Qntls",blockTabName="pcInRange_ExpDesign",  breaksTabName = "pcInRange_Breaks", dsn="modOffScores"){
	
	
	saveSuccessful <- saveBinExp.exp(
		experiment, 
		dsn = dsn, 
		blockTabName = blockTabName,  
		breaksTabName = breaksTabName)
	
	if(saveSuccessful){
		saveBinExp.quantiles(
			experiment = experiment,
			binPCQuantiles = binPCQuantiles , 
			binQuantilesTabName = binQuantilesTabName,
			dsn=dsn)
	
	}
	

}
	

#----------------------------------------------------------------------------

extractPCBinQuantiles <- function(experiment, pcBinResults, quantiles = c(0.1,0.25,0.50, 0.75, 0.9)){
		
	numBins <- length(experiment$breaks)-1
	binColNames <- paste("Bin", seq(1,numBins), sep="")
	numQuantiles <- length(quantiles)

	obNums <- sort(unique(pcBinResults$obNum))
	obQls <- array(dim=c(numQuantiles*length(obNums),numBins+2))
	for (oN in obNums ){
		pcRes.oN <- pcBinResults[pcBinResults$obNum == oN, binColNames ]
		qls <- apply(pcRes.oN,2,quantile, quantiles)
		obQls[seq(1+(oN-1)*numQuantiles, oN * numQuantiles),] <- cbind(obNum=oN, quantile=quantiles, qls)
	}
	
	obQls <- cbind(ID=experiment$ID, obQls)
	colnames(obQls) <- c("ID","obNum", "quantiles", binColNames)	
	return(as.data.frame(obQls))
}

#----------------------------------------------------------------------------


saveBinExp.quantiles <- function(experiment, binPCQuantiles, binQuantilesTabName = "pcInRange_quant", dsn="modOffScores"){
		
	require(RODBC)
	chnl <- odbcConnect(dsn)
	sqlInfo <- checkIfIDExists(ID = experiment$ID, tableName = binQuantilesTabName, dsn=dsn)
	
	if(sqlInfo$IDExists) stop("ID already exists, either change ID number or experiment already run")
		

	if(sqlInfo$tableExists){
		# append
		sqlSave(channel= chnl, dat=binPCQuantiles, tablename=binQuantilesTabName, append=TRUE, rownames=FALSE)
	} else {
		#create table
		sqlSave(channel= chnl, dat=binPCQuantiles, tablename=binQuantilesTabName, append=FALSE, rownames=FALSE)
	}
	
	odbcClose(chnl)
}


#----------------------------------------------------------------------------

removeIDs <- function(IDsToKill, binQuantilesTabName=NULL, blockTabName="pcInRange_ExpDesign",  breaksTabName = "pcInRange_Breaks", dsn = "modOffScores" ){
	
	require(RODBC)
	chnl <- odbcConnect(dsn)
	N <- length(IDsToKill)
	
	qryString <- paste("DELETE FROM ", blockTabName, "WHERE ",  
		do.call(paste, as.list(paste(" ID = ",  IDsToKill[1:(N-1)], "OR"))),
		" ID = ", IDsToKill[N],";", sep=" ")
	
	sqlQuery(channel = chnl, query= qryString)
	
	qryString <- paste("DELETE FROM ", breaksTabName, "WHERE " , 
		do.call(paste, as.list(paste(" ID = ",  IDsToKill[1:(N-1)], "OR"))),
		" ID = ", IDsToKill[N],";", sep=" ")

	sqlQuery(channel = chnl, query= qryString)

	
	if(!is.null(binQuantilesTabName)){
		qryString <- paste("DELETE FROM ", binQuantilesTabName, "WHERE " , 
		do.call(paste, as.list(paste(" ID = ",  IDsToKill[1:(N-1)], "OR"))),
		" ID = ", IDsToKill[N],";", sep=" ")

		sqlQuery(channel = chnl, query= qryString)
		
	}
		
	odbcClose(chnl)
	
}


#----------------------------------------------------------------------------

getIDFromParams <- function(params, dsn="modOffScores", blockTabName="pcInRange_ExpDesign", breaksTabName= "pcInRange_Breaks", tolerancePC = 0.00001){
	# params is a list with the necessary params named as below
	require(RODBC)
	chnl <- odbcConnect(dsn)
	p<-params
	
	qryString <- paste(
		"select ID from ", blockTabName, 
		" where ",
		 getSqlTol("trthmean", tolerancePC, p$mean ) ," and ",
		 getSqlTol("trthvar", tolerancePC, p$var ) ," and ",
		 getSqlTol("trthweightL", tolerancePC, p$weightL ) ," and ",
		 getSqlTol("trthweightG", tolerancePC, p$weightG ) ," and ", 
 		 " xLen = ", p$xLen, " and ",
		 " seed = ", p$seed," and ",
		 " numObs = ", p$numObs," and ",					
		 " xMin = ", p$xMin, " and ",
		 " xMax = ", p$xMax, " and ",
		 " binQuantilesTabName = '", p$binQuantilesTabName,
		 "' ;",
		sep="" )
	
	IDs <- sqlQuery(channel=chnl, query=qryString)
	
	if(dim(IDs)[1]>1) stop("error:  ID is not uniquely defined by the parameters and experiment table - something has gone wrong")
	
	#IDmatchesBreaks <- NULL
	
	#for( ID in IDs[,"ID"]){
		
	#	IDbreaks <- sort(sqlQuery(channel=chnl, query = paste("select experimentbreaks from ", breaksTabName, " where ID = ", ID ," ;", sep="" ))[,"experimentbreaks"])
		
	#	IDmatchesBreaks <- c(IDmatchesBreaks , all(IDbreaks ==  params$breaks))
		
	#}
	
	#IDval <- IDs[,"ID"][IDmatchesBreaks]
	odbcClose(chnl)
	return(IDs[1,1])
}


#----------------------------------------------------------------------------
getBinGlideForIDandQuantile <- function(ID, quantile, binNum, binQuantilesTabName, dsn){
	require(RODBC)
	chnl <- odbcConnect(dsn)
	
	qGlides <- sqlQuery(
		channel=chnl,
		query = paste(
			"select obNum, bin", binNum, 
			" from ",binQuantilesTabName,
			" where ID = ", ID ,
			" and quantiles = ", quantile, " ;"	 , sep=""))
	
	qGlides <- qGlides[order(qGlides$obNum),]
	odbcClose(chnl)
	return(qGlides)
}


#----------------------------------------------------------------------------

getBinCrossingTime <- function(IDtruth, IDfcst,binNum,  q, p,binQuantilesTabName, dsn){
			
	# If I want to be q% confident that the observed bin proportion is outside of the 1-2(1-q) range:  see whether if falls either above the q% quantile or below the 1-q quantile.  
	
	# suppose we want a probability of p% that we'll be able to reject the forecast by time t with q% confidence - 
	
	# trth above fcst:  use the 1-p quantile against the q quantile
	# fcst above trth:  use the p quantile against the 1-q quantile
	
	ultBinPCTruth <- getUltBinProportion(
		ID = IDtruth, 
		binNum = binNum, 
		binQuantilesTabName = binQuantilesTabName, 
		dsn = dsn )
	
	ultBinPCfcst <- getUltBinProportion(
		ID = IDfcst, 
		binNum = binNum, 
		binQuantilesTabName = binQuantilesTabName, 
		dsn = dsn )
	
	if(ultBinPCfcst > ultBinPCTruth){
	# fcst above trth:  use the k quantile against the 1-x quantile
		qFcst <- 1 - q
		qTruth <- p
		diffScalar <- -1  
	} else {
	# trth above fcst:  use the 1-p quantile against the x quantile
		qFcst <- q
		qTruth <- 1-p
		diffScalar <- 1  
	}
	
	keyQl.exp <- getBinGlideForIDandQuantile(
		ID = IDfcst, 
		quantile = qFcst,
		binNum = binNum, 
		binQuantilesTabName = binQuantilesTabName, 
		dsn = dsn)[,2]   #this is "expected" or "forecast"
	
	keyQl.obs <- getBinGlideForIDandQuantile(
		ID = IDtruth, 
		quantile = qTruth,
		binNum = binNum, 
		binQuantilesTabName = binQuantilesTabName, 
		dsn = dsn)[,2]   #this is "observed" or "truth"
	

	diffQl <- (keyQl.obs - keyQl.exp)  #if t crossing < lnQ this will be positive for diffScalar =1 and negative otherwise
		
		lnQ <- length(diffQl)
		i <- 1
		
		repeat {
			if((all(diffScalar * diffQl[i:lnQ] > 0) | i==lnQ)) break
			i <- i + 1
			
		}	
		
		if (i ==  lnQ) i <- 10^(ceiling(log(lnQ, base=10))+1) -1
		
	return(i)	
}


#----------------------------------------------------------------------------
getUltBinProportion <- function(ID, binNum,binQuantilesTabName, dsn  ){
	
	binMedian <- getBinGlideForIDandQuantile(ID = ID, quantile = 0.5, binNum=binNum,  binQuantilesTabName=binQuantilesTabName, dsn=dsn )
	lastRow <-  dim(binMedian)[1]
	ultProportion <- binMedian[lastRow,2]
	
	return(ultProportion)
}


#----------------------------------------------------------------------------

getExperimentParams <- function(ID , blockTabName ,   
breaksTabName, dsn){
	
	require(RODBC)
	chnl <- odbcConnect(dsn=dsn)
	
	sqlQry <- paste("select * from ", blockTabName, " where ID = ", ID, " ;", sep=" ")
	bV <- sqlQuery(channel = chnl, query = sqlQry)

	sqlQry <- paste("select * from ", breaksTabName, " where ID = ", ID, " ;", sep=" ")
	bR <- sqlQuery(channel = chnl, query = sqlQry )
	
	xVals <- seq(bV$xMin, bV$xMax, length.out=bV$xLen)
	breaks <- sort(bR[,2])
	
	experiment <- list(
		truth = list(mean=bV$trthmean, var=bV$trthvar, weightL = bV$trthweightL, weightG = bV$trthweightG ),
		xVals = xVals,
		breaks = breaks,
		numObs = bV$numObs, 
		seed = bV$seed,
		ID = ID
		)
	
	return(experiment)
	
}


#----------------------------------------------------------------------------

plotPDFwithBins <- function(ID, blockTabName, breaksTabName, dsn, maxX=5, new=TRUE, labelBins=TRUE, ...){
	expr <- getExperimentParams(ID = ID, blockTabName = blockTabName, breaksTabName = breaksTabName, dsn = dsn)


	companyForecast <- getForecastPDF_weighted(
			xVals = expr$xVals, 
			mean=expr$truth$mean, 
			var=expr$truth$var, 
			weight.lnorm=expr$truth$weightL, weight.gamma=expr$truth$weightG)


	cFpdf <- approxfun(x=companyForecast$x, y=companyForecast$y)

	if(new){
	plot(companyForecast, type="l", xlim=c(0,maxX), xlab="", ylab="", ...)
	}


	for(i in 1:(length(expr$breaks)-1)){
	
		xRug <- seq(expr$breaks[i], min(expr$breaks[i+1], maxX), length.out=500)
		xB <- c(xRug, rev(xRug), xRug[1])
		yB <- c(rep(0, length(xRug)), cFpdf(rev(xRug)), 0 )
		polygon(x=xB, y=yB, col=i+1)
		if(labelBins){
		text(x=median(xRug), y=max(yB)/5, labels=paste("Bin ", i, sep=""))
		}
	}


}






