#--------------------------------------------------------------------------------
# assumes a sql database has been set up called scoreData (with the same dsn) to store the results in


#--------------------------------------------------------------------------------
# get required functions


source("plotWinningProportions.r")
source("requiredFunctions.r")
#--------------------------------------------------------------------------------
# Set output directory
root <- getwd()
setwd(  paste(root, "\outputs", sep=""))

#N <- 4  # once N=7 the true forecast wins 99.3%,  N=1 low wins,  N=2 low and med very close, N=3 true forecast wins 49%
M <- 10
#scoreType <- "ignorance"
eps <- 0.00000001

	# create 3 forecasts:  N(0, 1/sqrt(2))  N(0,1)    N(0,2)
	
	xVals <- list(
		xLow = c(
			seq(-6,-3-0.001, length.out=2^4), 
			seq(-3,-1-0.001, length.out=2^5),
			seq(-1,1-0.001,length.out=2^7), 
			seq(1,3-0.001, length.out=2^5), 
			seq(3,6, length.out=2^4)),
			
		xMed = c(
			seq(-6,-4-0.001, length.out=2^4), 
			seq(-4,-1.5-0.001, length.out=2^5),
			seq(-1.5,1.5-0.001,length.out=2^7), 
			seq(1.5,4-0.001, length.out=2^5), 
			seq(4,6, length.out=2^4)),	
	
		xHigh = c(
			seq(-6,-5-0.001, length.out=2^4), 
			seq(-5,-2-0.001, length.out=2^5),
			seq(-2,2-0.001,length.out=2^7), 
			seq(2,5-0.001, length.out=2^5), 
			seq(5,6, length.out=2^4))
	
		)
	
	
	fcsts <- list(
		fLow <- list(x=xVals$xLow,   y=dnorm(xVals$xLow, mean=0, sd= 1/sqrt(2)) ),
		fMed <- list(x=xVals$xMed,   y=dnorm(xVals$xMed, mean=0, sd= 1) ),
		fHigh <- list(x=xVals$xHigh,   y=dnorm(xVals$xHigh, mean=0, sd= sqrt(2)) )
	)


scoreNames <- c( "powerrule", "powerrule", "powerrule", "spherical","ignorance", "continuousrankedprobability", "meansquarederror","naivelinear")


scoreParams <- list(
	list(alpha=1.5, simpleIntegrate=TRUE), 
	list(alpha=2,simpleIntegrate=TRUE), 
	list(alpha=2.5,simpleIntegrate=TRUE), 
	list(simpleIntegrate=TRUE),
	NULL,
	list(simpleIntegrate=TRUE),
	list(simpleIntegrate=TRUE),
	NULL)



require(RODBC)
chnl <- odbcConnect("scoreData")

scoresWanted <- 1:8  


for(scrCounter in scoresWanted){
 for(N in 1:7){
  simResult <-  list()
  for (sim in 1:2^M){
		
	# sample 2^N points from N(0,1)
	obs <- rnorm(2^N, mean=0, sd=1)	

	# calc average score
	
	meanScore <- list()
	for(i in 1:3){
			scores <- NULL
			for(obVal in obs){

				scores <- c(scores, calcScore(scoreType = scoreNames[scrCounter], forecastPDF=fcsts[[i]], verification= obVal, scoreParams[[scrCounter]]))
	
			} # next obs
		meanScore[[i]] <- mean(scores)
		} # next obVal
	
		meanScore <- do.call(c, meanScore)

		minScr <- min(meanScore)
		# count 1 for the forecast that does best - 0 otherwise
		simResult[[sim]] <- (meanScore < minScr + eps)  &  (meanScore > minScr - eps)

		print(paste(scoreNames[scrCounter], " : ", N, " : ",   sim, sep=""))
		
	} # Next sim    repeat 2^M times - see the proportions
	simResult <- do.call(rbind, simResult)

	winningProportions <- apply(simResult, 2, mean)
	# save results to Database

	if(is.null(scoreParams[[scrCounter]]$alpha)) {
		 sP <- 99999
	} else {
		sP <- scoreParams[[scrCounter]]$alpha
	}


	rowToSave <- data.frame(scoreType=scoreNames[scrCounter], scoreParam = sP ,  M=M, N=N, fcstType=c("low", "med", "high"), winningProportions)


	sqlSave( channel=chnl, dat=rowToSave, tablename ="sparseGaussian" ,append=TRUE, rownames=FALSE)

 } # next N
 

} # next scoreType



odbcClose(chnl)

	
##############################################################################
# Analyse results


scoreNames <- c( "powerrule", "powerrule", "powerrule", "spherical","ignorance", "continuousrankedprobability", "meansquarederror","naivelinear")


scoreParams <- list(
	list(alpha=1.5, simpleIntegrate=TRUE), 
	list(alpha=2,simpleIntegrate=TRUE), 
	list(alpha=2.5,simpleIntegrate=TRUE), 
	list(simpleIntegrate=TRUE),
	NULL,
	list(simpleIntegrate=TRUE),
	list(simpleIntegrate=TRUE),
	NULL)




require(RODBC)
chnl <- odbcConnect("scoreData")
scrRes <- sqlFetch(chnl, sqtable="sparseGaussian")
odbcClose(chnl)

scoreType <- "ignorance"


for(i in 1:length(scoreNames)){

	sP <- scoreParams[[i]]$alpha
	if(is.null(sP)) sP <- 99999	

	sPnm <- ifelse(sP == 99999, "", paste("_",sP, sep=""))
	pdf(paste("winningProportions_",scoreNames[i], sPnm , "20150310.pdf", sep="" ))
		plotWinningProportions(scrRes=scrRes, scoreType =scoreNames[i], scoreParam=sP)
	dev.off()
}


resMed <- list()
for(i in 1:8){
	sP <- scoreParams[[i]]$alpha
	if(is.null(sP)) sP <- 99999	
	resMed[[i]] <- scrRes[scrRes$fcstType == "med"  & scrRes$scoreType == scoreNames[i] & 
		scrRes$scoreParam == sP,]
}

xx <- sort(as.numeric(unique(scrRes$N)))
yy <- seq(0,1, length.out=length(xx))
scrNameFull <- c( "powerrule1.5", "powerrule 2", "powerrule 2.5", "spherical","ignorance", "continuousrankedprobability", "meansquarederror","naivelinear")

pdf("winningProportions_N_0_1_allScores.pdf")
plot(x=xx, y=yy, type="n", xlab="log2(N)", sub="where N is number of observations", ylab=bquote(F[Perfect]))
legend(x=1, y=1, legend=scrNameFull, col=1:8, lwd=1, cex=0.9, border="white", bty="n")
for(i in 1:8){

	lines(resMed[[i]][,"winningProportions"], col=i)	
}


dev.off()









