#-----------------------------------------------------------------------------------------------------------
#  READ ME

# This file produces all the graphics for the report triangleResults.pdf.   
# It draws in results from the tables tabled triangleInfDefcQuantiles*
# Variants covered by * relate to test performed -  see the file readMe_createDescription for details
# The script to produce the various graphics are shown below in the sections delimited with #----

# Graphics produced include:
#   exampleGlideQuantiles.pdf
#   forecastVsTruthExample.pdf
#	iDGlideExample.pdf
#	triangleResultsTimeCrossing.pdf
#	triangleResultsTimeCrossing512.pdf
#	triangleResultsTimeCrossing512.pdf
#	paste("squareResultsTimeCrossingmV512_*
#	paste("squareResultsTimeCrossingmV2048_s2_*
#	triTestSeeds_i14j9.pdf
#   <MULTIPLE SCORES GRAPHIC - TO BE NAMED>


#-----------------------------------------------------------------------------------------------------------
# SOURCE FILES FOR REQUIRED FUNCTIONS



source("modOffScores_sourcingCode.r")
source("triangle_sourcingCode.r")




scoreTypes <- c("ignorance", "properlinear","naivelinear", "powerrule","powerrule","powerrule","spherical","meansquarederror")

#"continuousrankedprobability"  excluded - too slow

scoreParams <- list(
NULL,
list(simpleIntegrate=TRUE),
NULL,
list(alpha=1.5, simpleIntegrate=TRUE),
list(alpha=2, simpleIntegrate=TRUE),
list(alpha=2.5, simpleIntegrate=TRUE),
list(simpleIntegrate=TRUE),
list(simpleIntegrate=TRUE)
)





#---------------------------------------------------------------------------------------
# PRODUCE EXPLANATORY GRAPHICS (you need to source in the required functions etc from above)
# 
#    exampleGlideQuantiles.pdf
#    forecastVsTruthExample.pdf
#	 iDGlideExample.pdf


N <- 2^9 # number of verifications
len <- 2^13
xVals <- seq(0.0001,15, length.out=len)
scoreType <- "ignorance"
scoreParam <- NULL
numGlides <- 2^10   # set to 10 for main runs
triangleSideDivisions <- 2^2
initialSeed <- 1236457
meanVal <- 1
varVal <- 0.65
resTabName <-  "triangleInfDefctQuantiles512"
probObsOsExp <- 0.75 
quantileExp <- 0.75

weights <- getTriangleArray(N=triangleSideDivisions)
# to avoid infinite scores with Pareto distribution
weights[weights[,1]==0  &  weights[,2]==0 & weights[,3]==1,] <- c(0.025, 0.025, 0.95)

truth=list(mean=1, var=0.65, weightL=weights[5,1], weightG=weights[5,2])
forecast=list(mean=1, var=0.65, weightL=weights[12,1], weightG=weights[12,2])
pdf("exampleGlideQuantilesFAB_below.pdf")
	plotGlidePathFAB(
		truth = truth, 
		forecast = forecast, 
		scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
		otherInfo = list(numVerifications = N, len = len,seedVal = initialSeed,
		rangeX0 = range(xVals)[1],rangeX1 = range(xVals)[2]),
		probObsOsExp = probObsOsExp ,
		quantileExp = quantileExp, 
		plotKeyQuantile = TRUE, 
		dsn = "modOffScores",
		resTabName = resTabName
	)

dev.off()

probObsOsExp <- 0.75
quantileExp <- 0.9

truth=list(mean=1, var=0.65, weightL=weights[14,1], weightG=weights[14,2])
forecast=list(mean=1, var=0.65, weightL=weights[8,1], weightG=weights[8,2])
pdf("exampleGlideQuantilesFAB_above.pdf")
	plotGlidePathFAB(
		truth = truth, 
		forecast = forecast, 
		scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
		otherInfo = list(numVerifications = N, len = len,seedVal = initialSeed,
		rangeX0 = range(xVals)[1],rangeX1 = range(xVals)[2]),
		probObsOsExp = probObsOsExp ,
		quantileExp = quantileExp, 
		plotKeyQuantile = TRUE, 
		dsn = "modOffScores",
		resTabName = resTabName
	)

dev.off()



compParms = c(meanVal, varVal, weights[5,1], weights[5,2])
trueParms = c(meanVal, varVal, weights[15,1], weights[15,2])

gR.sgl <- getGlideResults(
	compParms,
	trueParms,
	scoreType = "ignorance",
	scoreParam = NULL,
	xVals = xVals, 
	numGlides = 1,
	numVerifications = 2^13,
	start.seed = initialSeed+1, 
	len=len)


companyForecast <- getForecastPDF_weighted(
			xVals, 
			mean=compParms[1], 
			var=compParms[2], 
			weight.lnorm=compParms[3], weight.gamma=compParms[4])

trueDist <- getForecastPDF_weighted(
			xVals, 
			mean=trueParms[1], 
			var=trueParms[2], 
			weight.lnorm=trueParms[3], weight.gamma=trueParms[4])

require(Hmisc)
pdf("forecastVsTruthExampleFAB.pdf")
	plot(companyForecast, type="l", xlab="loss", ylab="density", col="blue")
	lines(trueDist, col="Red")
	legend(x=c(10,10), y=c(0.4,0.3), legend=c("forecast", "truth"), col=c("blue", "red"), 		lwd=1, cex=0.75, border=NULL, bty="n")
	subplot(plotTriangleWithTwoDots(compParms[3], compParms[4], trueParms[3], 					trueParms[4]), x=c(13,14.5), y=c(2.1,2.4))
dev.off()

pdf("iDGlideExampleFAB.pdf")
	plot(gR.sgl$fcstTime, gR.sgl$informationDeficit, type="l", 
		ylab="empirical information Deficit", xlab="time" , col="blue")
	eScoreFcst <- calcExpectedScore2("ignorance", forecastPDF = companyForecast, 
		underlyingPDF = companyForecast, scoreParams = NULL)
	eScoreTruth <-calcExpectedScore2("ignorance", forecastPDF = companyForecast, 
		underlyingPDF = trueDist, scoreParams = NULL) 
	abline(h= -(eScoreFcst - eScoreTruth), col="red" )
dev.off()




#------------------------------------------------------------------------------------
# TRIANGLE VALUES:  TIME OUT = 512 
#	triangleResultsTimeCrossing512.pdf


N <- 2^9 # number of verifications
len <- 2^13
xVals <- seq(0.0001,15, length.out=len)
scoreType <- "ignorance"
scoreParam <- NULL
numGlides <- 2^10   # set to 10 for main runs
triangleSideDivisions <- 2^2
initialSeed <- 1236457
meanVal <- 1
varVal <- 0.65
resTabName <-  "triangleInfDefctQuantiles512"
probObsOsExp <- 0.75 
quantileExp <- 0.75

resTab <- NULL
for (i in 1:15){ 
	
	#  This is now the Expected truth - i.e. i = the company forecast and this glide quantile set is what you would see if the company had correctly guessed the underlying distribution.   
	qExpectedGlides <- getQuantileData(
			truth=list(mean=1, var=0.65, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			forecast=list(mean=1, var=0.65, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = initialSeed,rangeX0  =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
	)
		
	
	
	for (j in seq(along = weights[,1])){
		
		print(paste(i, ":",j ))
		
		
		#  This is where the code is different from the truth as base case (and in fact the ONLY place it is different).  We are now looping over actual underlying distributions j - given the forecast i.  So the truth row below is now j (it varies) and the forecast row below is i (it is fixed whilst j varies).   Thus the roles of i and j are reversed below (compared to the glidePathQuantiles_analysis script see TruthAsBase folder).
		qObservedGlides <- getQuantileData(
			truth=list(mean=1, var=0.65, weightL=weights[j,1], 				
				weightG=weights[j,2]), 
			forecast=list(mean=1, var=0.65, weightL=weights[i,1],
				weightG=weights[i,2]), 
			scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = initialSeed, rangeX0 =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
		)
	
		

		
	# The 'truth' glide path  (qExpectedGlides) is actually the expected truth (distribution i);  whereas qObservedGlides gives the glide path you would see given the same forecast that underlies the expected truth (i) - but if the real truth is equal to distribution j.
	
	# You can calculate the expected truth glide path in advance
	# You can observe the actual empirical information deficit and note if it crosses various quantiles - crossing a quantile (k, say) says "there was only a k-chance of that occuring if our forecast is right. 
	#  I THINK THESE ARE THE SAME CONCEPT AS THE BANANA GRAPHS WE USE AT LLOYD'S !!!
	#  The point of this work is to say how quickly you'll know your forecast has a k-chance of being wrong.
	# the speend of certainty will depend on the true distribution j.  If this is very different from i you'll know quickly.
	# because j is random it can produce average values that look like i would produce - but this gets less and less likely as time passes and you add more and more values to the average.
	# the qFcstData glidepath shows you how likely the informationDeficit path is to be at a given point - at a particular time -if the true distirbution is j (but you are subtracting an expectation of i)
	
	# So you get a statement like "there is only a d(t)% chance of not being above the k(t)-quantile by time t if the true distribution is j and your forecast is i"  or 
	
	# "by time t, there is a d(t)% probability of being 1-k% confident that your choice is wrong if truth=j and forecast = i 
	
	
	# In summary the crossing point shows how quickly you can expect to be  
	
	
		if(dim(qObservedGlides)[1]>0 & dim(qExpectedGlides)[1]>0 ){
		tC <- getTimeCrossingFAB2(
			qObservedGlides = qObservedGlides, 
			qExpectedGlides = qExpectedGlides, 
			probObsOsExp = probObsOsExp ,
			quantileExp = quantileExp )		
		
		
		# so the true weight is not the truth - but the forecast
		# the fcst weight is not the fcst - but it is the true distribution
		# the tC variable should be interpreted as described above
		resRow <- data.frame(
			trueWeightL=weights[i,1], 				
			trueWeightG=weights[i,2],
			fcstWeightL=weights[j,1],
			fcstWeightG=weights[j,2],
			timeCrossing = tC
			)
		
		resTab <- rbind(resTab, resRow)
		} else {
			print(paste("combination not available: i=", i, "  j=", j, sep=""))
			
		}
		
	}   # next j - forecast

}  # next i - truth 



dim(resTab)
pdf(paste("triangleResultsFABTimeCrossing512", "_obs", probObsOsExp*100, "_exp", quantileExp*100,".pdf", sep=""))
par(mfrow=c(4,4))
for (i in 1:15){ 
	
	# Now the trueWeight is actually the forecast - the BLUE dot is the forecast
	# The numbers show the time (in score time intervals - i.e. the interval in time between the calculation of each score) at which there is a k% chance of being k% confident that your forecast is wrong if j <> i  
	
	trueWL <- weights[i,1]
	trueWG <- weights[i,2]
	triDat <- resTab[resTab$trueWeightL == trueWL & resTab$trueWeightG == trueWG, ]
	par(mar=c(1,1,1,1))
	plot(x=range(weights[,1])*1.1+c(-0.1,0), y = range(weights[,2])*1.1+c(-0.1,0), type="n", axes=FALSE, ann=FALSE, col="grey")
	
	if(i %in% c(1,2,3,4)) {
		mtext(side = 3, text=i, line=-0.5)
		if(i ==1) mtext(side=2, text=1, las=2)	
	}

	if( i == 5) mtext(side = 2, text=2, las=2)
	if( i == 9) mtext(side = 2, text=3, las=2)
	if( i == 13) mtext(side = 2, text=4, las=2)
	
	
	textWant <- !( (triDat$fcstWeightL == trueWL)  & (triDat$fcstWeightG == trueWG) ) 
	textCex <- (triDat$timeCrossing[textWant]==9999) * 0.4 +
		 (1-(triDat$timeCrossing[textWant]==9999))*0.7
	text(x=triDat$fcstWeightL[textWant], y = triDat$fcstWeightG[textWant], 
		labels=triDat$timeCrossing[textWant], cex=textCex)
	points(x=trueWL, y=trueWG, pch=19, col="blue")
}
dev.off()




pdf(paste("allGlidesTriFAB512_i", i, "j", j, "_obs", probObsOsExp*100, "_exp", quantileExp*100, ".pdf", sep=""))

for (i in 1:15){ 
 forecast=list(mean=1, var=0.65, weightL=weights[i,1], weightG=weights[i,2])
 for (j in 1:15){ 

    
	truth=list(mean=1, var=0.65, weightL=weights[j,1], weightG=weights[j,2])

	plotGlidePathFAB(
		truth = truth, 
		forecast = forecast, 
		scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
		otherInfo = list(numVerifications = N, len = len,seedVal = initialSeed,
		rangeX0 = range(xVals)[1],rangeX1 = range(xVals)[2]),
		probObsOsExp = probObsOsExp ,
		quantileExp = quantileExp, 
		plotKeyQuantile = TRUE, 
		dsn = "modOffScores",
		resTabName = resTabName
	)
	
	print(paste(i, j, sep=" : "))
 }
}

dev.off()





#---------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------
# SQUARE PLOT:  differing levels of mean and var  (TIME OUT = 2048)
#	paste("squareResultsTimeCrossingmV2048_s2_", keyQuantile*100, ".pdf"


# I'M NOT SURE THIS WAS SYMMETRIC?
# YES - THIS ONE WILL NEED TO BE RE-RUN WITH THE FORECAST AT THE CENTRE AND THE TRUTH VARYING...


N <- 2^11 # number of verifications
len <- 2^13
xVals <- seq(0.0001,15, length.out=len)
scoreType <- "ignorance"
scoreParam <- NULL
numGlides <- 2^10   # set to 10 for main runs
triangleSideDivisions <- 2^2
initialSeed <- 1236457
meanValCtr <- 1
varValCtr <- 0.65
resTabName <-  "triangleInfDefctQuantilesFABmV2048"

weights <- getTriangleArray(N=triangleSideDivisions)
# to avoid infinite scores with Pareto distribution
weights[weights[,1]==0  &  weights[,2]==0 & weights[,3]==1,] <- c(0.025, 0.025, 0.95)
probObsOsExp <- 0.75 
quantileExp <- 0.75 

i <- 8
j <- 8





resTab <- list()
k <- 1

# This is the expected truth
qExpectedGlides <- getQuantileData(
			truth=list(mean=meanValCtr, var=varValCtr, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			forecast=list(mean=meanValCtr, var=varValCtr, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = initialSeed,rangeX0  =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
	)

for (meanVal in c(1/1.2^3,1/1.2^2, 1/1.2,1,1.2,1.2^2) * meanValCtr){
	for (varVal in c(1/1.2^3,1/1.2^2, 1/1.2,1,1.2,1.2^2)* varValCtr){	
	
		print(paste(meanVal, varVal, k, sep=" : "))		
	
		# This is the quantiles with true truth - and the forecast as expected
		qObservedGlides <- getQuantileData(
			truth= list(mean=meanVal, var=varVal, weightL=weights[j,1],
				weightG=weights[j,2]),
			forecast=list(mean=meanValCtr, var=varValCtr, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = initialSeed, rangeX0 =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
		)
	
	
	

	
		
		if(dim(qObservedGlides)[1]>0 & dim(qExpectedGlides)[1]>0 ){
		tC <- getTimeCrossingFAB2(
			qObservedGlides = qObservedGlides, 
			qExpectedGlides = qExpectedGlides, 
			probObsOsExp = probObsOsExp ,
			quantileExp = quantileExp )	
		
		
		resRow <- data.frame(
			trueMean = meanValCtr,
			trueVar = varValCtr,
			fcstMean = meanVal,
			fcstVar = varVal,
			timeCrossing = tC
			)
		
		resTab[[k]] <- resRow
		} else {
			print(paste("combination not available: mean = ",meanVal, "  var = ", varVal, sep=""))
			
		}

	k <- k+1
	
	
	}   # next varVal

}  # next meanVal

resTab <- do.call(rbind, resTab)

pdf(paste("squareResultsTimeCrossingFABmV2048_s2_", i, "j", j, "_obs", probObsOsExp*100, "_exp", quantileExp*100,".pdf", sep=""))


	par(mar=c(3,3,1,1))
	plot(x=resTab$fcstMean, y = resTab$fcstVar, type="n", col="grey", axes=FALSE)
	axis(side=1, at=c(1/1.2^3,1/1.2^2, 1/1.2,1,1.2,1.2^2) * meanValCtr, cex.axis=0.5)
	axis(side=2, at=c(1/1.2^3,1/1.2^2, 1/1.2,1,1.2,1.2^2) * varValCtr, cex.axis=0.5)
	mtext("mean", side=1, line=2)
	mtext("var", side=2, line=2)
	textWant <- !( (resTab$fcstMean == meanValCtr)  & (resTab$fcstVar == varValCtr) ) 
	textCex <- (resTab$timeCrossing[textWant]==99999) * 0.4 + (1-(resTab
		$timeCrossing[textWant]==99999))*0.7
	text(x=resTab$fcstMean[textWant], y = resTab$fcstVar[textWant], labels=resTab
		$timeCrossing[textWant], cex=textCex)
	points(x=meanValCtr, y=varValCtr, pch=19, col="blue")
	

dev.off()

forecast=list(mean=1, var=0.65, weightL=weights[8,1], weightG=weights[8,2])

pdf(paste("allGlidesFAB_mV2048_i", i, "j", j, "_obs", probObsOsExp*100, "_exp", quantileExp*100, ".pdf", sep=""))

for (meanVal in c(1/1.2^3,1/1.2^2, 1/1.2,1,1.2,1.2^2) * meanValCtr){
 for (varVal in c(1/1.2^3,1/1.2^2, 1/1.2,1,1.2,1.2^2)* varValCtr){	
	truth=list(mean=meanVal, var=varVal, weightL=weights[8,1], weightG=weights[8,2])

	plotGlidePathFAB(
		truth = truth, 
		forecast = forecast, 
		scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
		otherInfo = list(numVerifications = N, len = len,seedVal = initialSeed,
		rangeX0 = range(xVals)[1],rangeX1 = range(xVals)[2]),
		probObsOsExp = probObsOsExp ,
		quantileExp = quantileExp, 
		plotKeyQuantile = TRUE, 
		dsn = "modOffScores",
		resTabName = resTabName
	)
	
	print(paste(meanVal, varVal, sep=" : "))
 }
}


dev.off()


pdf("illustrateDiffMV_FAB.pdf")
	pdfFcst <- getForecastPDF_weighted(xVals, mean=1, var=0.65, 0.5, 0.25)
	pdfTruth <- getForecastPDF_weighted(xVals, mean=1/1.2^2, var=0.65/1.2^2, 0.5, 0.25)
	plot(pdfTruth, col="red", type="l", xlab="loss", ylab="probability density")
	lines(pdfFcst, col="blue", lwd=1)
	legend(x=10, y=1.5,legend=c("truth", "forecast"), col=c("red", "blue"), lwd=1, cex=0.8)
dev.off()






#---------------------------------------------------------------------------------------
# TESTING LOTS OF DIFFERENT SEEDS
#	triTestSeeds_i14j9.pdf

# I'M NOT SURE THIS WAS SYMMETRIC?
# YES - THIS ONE WILL NEED TO BE RE-RUN WITH THE FORECAST AT THE CENTRE AND THE TRUTH VARYING...



N <- 2^9 # number of verifications
len <- 2^13
xVals <- seq(0.0001,15, length.out=len)
scoreType <- "ignorance"
scoreParam <- NULL
numGlides <- 2^10   # set to 10 for main runs
triangleSideDivisions <- 2^2
initialSeed <- 1236457
meanValCtr <- 1
varValCtr <- 0.65
resTabName <-  "triangleInfDefctQuantileFABSd512"

weights <- getTriangleArray(N=triangleSideDivisions)
# to avoid infinite scores with Pareto distribution
weights[weights[,1]==0  &  weights[,2]==0 & weights[,3]==1,] <- c(0.025, 0.025, 0.95)
probObsOsExp <- 0.75 
quantileExp <- 0.75 


chnl <- odbcConnect("modOffScores")
glSd <- sqlFetch(chnl, resTabName)
seeds <- unique(glSd$seedVal)


i <- 5
j <- 14

resTab <- list()
k <- 1


for (seedVal in seeds){


		#EXPECTED TRUTH
		qExpectedGlides <- getQuantileData(
			truth=list(mean=meanValCtr, var=varValCtr, weightL=weights[j,1], 				
				weightG=weights[j,2]), 
			forecast=list(mean=meanValCtr, var=varValCtr, weightL=weights[j,1], 				
				weightG=weights[j,2]), 
			scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = seedVal,rangeX0  =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
			)

	
		#roles of i and j interchanged
		qObservedGlides <- getQuantileData(
			truth=list(mean=meanValCtr, var=varValCtr, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			forecast=list(mean=meanValCtr, var=varValCtr, weightL=weights[j,1],
				weightG=weights[j,2]), 
			scoreInfo = list(scoreType="ignorance", scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = seedVal, rangeX0 =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
		)
	
		
		if(dim(qObservedGlides)[1]>0 & dim(qExpectedGlides)[1]>0 ){
		tC <- getTimeCrossingFAB2(
			qObservedGlides = qObservedGlides, 
			qExpectedGlides = qExpectedGlides, 
			probObsOsExp = probObsOsExp ,
			quantileExp = quantileExp )	
		
		
		resRow <- data.frame(
			seedVal = seedVal,
			timeCrossing = tC
			)
		
		resTab[[k]] <- resRow
		} else {
			print(paste("combination not available: mean = ",meanVal, "  var = ", varVal, sep=""))
			
		}

	k <- k+1		
	}  # next seedVal



resTab <- do.call(rbind, resTab)
oldPar <- par()


baseSeedNum <- resTab$seedVal == initialSeed
xCount <- seq(along = resTab$timeCrossing) 

pdf("triTestSeedsFAB_i9j14.pdf", width=2, height=8)
	set.seed(1235)
	plot(c(0.98,1.02), range(resTab$timeCrossing), type="n", ann=FALSE, axes=FALSE)
	abline(h=resTab$timeCrossing, col="grey")
	points( jitter(rep(1, length(xCount)-1), factor=1), resTab$timeCrossing[!baseSeedNum] ,pch=21)
	axis(side=2, at=resTab$timeCrossing, las=2, cex.axis=0.5)
	#mtext(side=2, text=resTab$timeCrossing[baseSeedNum],at=resTab$timeCrossing[baseSeedNum], col="blue", las=2)
	points((1:xCount)[baseSeedNum], resTab$timeCrossing[baseSeedNum], col="black", pch=19)
dev.off()




#-----------------------------------------------------------------------------------------------------------
# TESTING DIFFERENT SCORES

# I'M NOT SURE THIS WAS SYMMETRIC?
# YES - THIS ONE WILL NEED TO BE RE-RUN WITH THE FORECAST AT THE CENTRE AND THE TRUTH VARYING...



N <- 2^9#  2^9 # number of verifications
len <- 2^13
xVals <- seq(0.0001,15, length.out=len)
#scoreType <- "ignorance"
#scoreParam <- NULL
numGlides <- 2^10 #2^10   # set to 10 for main runs
triangleSideDivisions <- 2^2
initialSeed <- 1236457
meanVal <- 1
varVal <- 0.65
resTabName <-  "triangleInfDefctQuantilesFABScr512"
probObsOsExp <- 0.5 
quantileExp <- 0.75 



weights <- getTriangleArray(N=triangleSideDivisions)
# to avoid infinite scores with Pareto distribution
weights[weights[,1]==0  &  weights[,2]==0 & weights[,3]==1,] <- c(0.025, 0.025, 0.95)

chnl <- odbcConnect("modOffScores")
scr512 <- sqlFetch(chnl, resTabName)
scores <- as.character(unique(scr512$scoreType))


resTab <- list()
k <- 1

for (scrVal in scores){
 for (i in 15:15){ 
	
	
	qExpectedGlides <- getQuantileData(
			truth=list(mean=1, var=0.65, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			forecast=list(mean=1, var=0.65, weightL=weights[i,1], 				
				weightG=weights[i,2]), 
			scoreInfo = list(scoreType=scrVal, scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = initialSeed,rangeX0  =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
	)
		
	
	
	for (j in seq(along = weights[,1])){
		
		print(paste(i, ":",j ))
		
		qObservedGlides <- getQuantileData(
			truth=list(mean=1, var=0.65, weightL=weights[j,1], 				
				weightG=weights[j,2]), 
			forecast=list(mean=1, var=0.65, weightL=weights[i,1],
				weightG=weights[i,2]), 
			scoreInfo = list(scoreType=scrVal, scoreParam=NULL),
			otherInfo = list(numVerifications = N, len = len,
				seedVal = initialSeed, rangeX0 =range(xVals)[1],
				rangeX1 = range(xVals)[2]),
			dsn="modOffScores",
			resTabName=resTabName
		)
	
		

		
		
		if(dim(qObservedGlides)[1]>0 & dim(qExpectedGlides)[1]>0 ){
		tC <- getTimeCrossingFAB2(
			qObservedGlides = qObservedGlides, 
			qExpectedGlides = qExpectedGlides, 
			probObsOsExp = probObsOsExp ,
			quantileExp = quantileExp )	
		
		
		resRow <- data.frame(
			trueWeightL=weights[i,1], 				
			trueWeightG=weights[i,2],
			fcstWeightL=weights[j,1],
			fcstWeightG=weights[j,2],
			scoreType = scrVal,
			timeCrossing = tC
			)
		
		resTab[[k]] <- resRow
		} else {
			print(paste("combination not available: i=", i, "  j=", j, sep=""))
			
		}
		k <- k+1 		
	}   # next j - forecast

 }  # next i - truth

} # next score

resTab <- do.call(rbind, resTab)


dim(resTab)
pdf(paste("triangleResults_diffScoresFAB_i15", "_obs", probObsOsExp*100, "_exp", quantileExp*100, ".pdf", sep=""))

i <- 15
par(mfrow=c(2,2))
 for (scrVal in scores){

	trueWL <- weights[i,1]
	trueWG <- weights[i,2]
	triDat <- resTab[resTab$trueWeightL == trueWL & resTab$trueWeightG == trueWG & resTab$scoreType == scrVal, ]
	par(mar=c(1,1,1,1))
	plot(x=weights[,1], y = weights[,2], type="n", axes=FALSE, ann=FALSE, col="grey")
	textWant <- !( (triDat$fcstWeightL == trueWL)  & (triDat$fcstWeightG == trueWG) ) 
	textCex <- (triDat$timeCrossing[textWant]==9999) * 0.4 +
		 (1-(triDat$timeCrossing[textWant]==9999))*0.7
	text(x=triDat$fcstWeightL[textWant], y = triDat$fcstWeightG[textWant], 
		labels=triDat$timeCrossing[textWant], cex=textCex)
	points(x=trueWL, y=trueWG, pch=19, col="blue")
	title(main=scrVal)

 } # next score
dev.off()


