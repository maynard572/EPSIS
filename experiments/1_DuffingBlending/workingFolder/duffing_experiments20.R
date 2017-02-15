###################################################################
# User settings

mainFolder <- "chooseAFolder"
resultSet <- "duffingExperiments"

scoreNames <- c( "powerrule", "powerrule", "powerrule", "spherical","ignorance", "properlinear", "continuousrankedprobability", "meansquarederror","naivelinear")


scoreParams <- list(list(alpha=1.5), list(alpha=2), list(alpha=2.5), rep(NULL,6))

# graphics for report
root <- "experiments20_"



#------------------------------------------------------------------
require(RODBC)
chnl <- odbcConnect(dsn="scoreData")
setwd( paste(mainFolder, "requiredFunctions", sep=""))
source("requiredFunctions.r")

setwd( paste(mainFolder, resultSet, sep=""))


#------------------------------------------------------------------

wantToGetScoresArray <- TRUE
if(wantToGetScoresArray){
	
# create initial sample values (will be given these in future)


sigmaU <- 0.1

len <- 2^11
sN <- scoreNames
sP <- scoreParams 
seed.Sm <- 1002
seed.v  <- 2002

sigmaMs <- c(seq(0.5, 0.8, by=0.1), seq(0.825,1.4, by=0.025), 1.5,1.6,1.7)/10 



appendVal <- FALSE  # if you are doing extra experiments you need to set this to true so it will add to the existing table….


scoreVals <- data.frame()


for (nStep in c(32)){ # either 8 or 32
	dufDat <- sqlFetch(channel=chnl, sqtable=paste("dufDat_", nStep,sep=""))

	
for (dufCol in c(5,17,20)){ #any number from 1 to 32
	
	S <- dufDat[,dufCol][1:2^13]
	S_V <- S[(1:2^12)]
	
	
for (mk in c(7)){
	M <- 2^mk   # the number of verifications
for (k2 in 0:9){
for (r in c(5,6,7)){
for (p in 0:0){



scoreVals.sgl <- getKernelSmoothedScores.dataframe(
	S_V = S_V,
	S_E = S[seq(2^12+1+2^r * p, 2^12 + 2^r +2^r * p)],
	sigmaU = sigmaU,
	M = M,
	len= len,
	sN = sN,
	sP = sP,
	seed.Sm = seed.Sm+k2,
	seed.v = seed.v,
	sigmaMs = sigmaMs)
	
scoreVals.sgl <- data.frame(sigmaU = sigmaU, nStep = nStep, initialisation = dufCol, seedIndex=k2, fcstSize=r, fcstPos=p, vPower = mk, scoreVals.sgl)

#scoreVals <- rbind(scoreVals, scoreVals.sgl)

sqlSave(channel=chnl, dat = scoreVals.sgl, tablename = paste(root, "scoreVals",sep=""), rownames=FALSE, append=appendVal)
appendVal <- TRUE # after the first time you append to the table

rm(scoreVals)
gc()


} # next p
} # next r
} # next k2
} # next mk
} # next dufCol
} # next nStep







getUGraph <- function(scoreType){
	t <- scoreVals.sgl
	t <- t[t$scoreType == scoreType,]
	t1 <- aggregate(t, by=list(t$sigmaM), FUN=mean)
	minScore <- getMin(x=t1$sigmaM, y=t1$scoreValue)
	plot(t1$sigmaM, t1$scoreValue, type='l')
	text(x=minScore$x, y=minScore$y, labels=paste("(",minScore$x,",", minScore$y,")"))
	title(main=scoreType)
}

# use append=TRUE if you want to add this to an existing table



} # end get scores array


if (!exists(x="scoreVals")) {
	scoreVals <- sqlFetch(channel = chnl, sqtable = paste(root, "scoreVals",sep=""))
	
}



# example for match
h <- data.frame(k=sample(c("a","b","c"), size=10, replace=TRUE), x=rnorm(10))
g<- data.frame(k=c("a","b","c"), v=c(5,4,3))
g[match(h$k, g$k), "v"]




for (dufCol in c(5,17,20)){

scoreValsS <- scoreVals[scoreVals$initialisation == dufCol,]

meanVals <- aggregate(
	scoreValsS$scoreValue, 
	by=list(
		seedIndex = scoreValsS$seedIndex, 
		fcstSize = scoreValsS$fcstSize,
		fcstPos = scoreValsS$fcstPos, 
		vPower = scoreValsS$vPower, 
		sigmaM = scoreValsS$sigmaM, 
		scoreType = scoreValsS$scoreType),
	FUN = mean999)

colnames(meanVals)[7] <- "scoreValue"

minVals <- getMinVals19(meanVals)	
minVals <- cbind(minVals, fcstID = paste(minVals$fcstSize, minVals$xmin, sep="_"))
forecastQuality <- sqlFetch(chnl, "forecastQuality")
forecastQuality <- cbind(forecastQuality, fcstID = paste(forecastQuality$fcstSize, forecastQuality$sigmaM, sep="_") )

minVals <- minVals[!duplicated(minVals),]


sN <- unique(as.character(minVals$scoreType))
scoreCombn <- combn(sN,2)
numComb <- dim(scoreCombn)[2]

wantEllipses <-FALSE
wantJitter <- FALSE
trueParm <- 0.1


for (fcstSize in 5:7){
 minValsS <- minVals[minVals$fcstSize==fcstSize,]
 pdf(paste(root,"scoreCompareDots_duffing_",dufCol,"fcstSize_",fcstSize,".pdf",sep=""))

for (inum in 1:numComb){

s1Name <- scoreCombn[1,inum]
s2Name <- scoreCombn[2, inum]

	#pch <- paste( 
	#	"(" , 
	#	minVals[minVals$scoreType== s1Name,"vPower"] ,
	#	",",
	#	minVals[minVals$scoreType== s1Name,"fcstSize"] ,
	#	")",
	#	 sep="")
	
	pch <- minVals[minVals$scoreType== s1Name,"seedIndex"] 
	
	if(wantEllipses){
		rd <- getRadii(radCol="vPower", s1Name=s1Name, s2Name=s2Name, trueParm=trueParm, FUN="median")
	} else {
		rd = list(rads1=NULL, rads2=NULL, radVals=NULL)
		
	}

	
	s1Val <- minValsS[minValsS$scoreType== s1Name,"xmin"]
	s2Val <- minValsS[minValsS$scoreType== s2Name,"xmin"]

	if(wantJitter){
		set.seed(12345) # to keep the jittering the same
		s1Val <- jitter(s1Val, factor=0.3)
		s2Val <- jitter(s2Val,factor=0.3)
	} 


	compareScores(
		s1 = s1Val,
		s2 = s2Val,
		s1Name = s1Name ,
		s2Name = s2Name ,
		trueParm = trueParm ,
		pch = pch ,
		forecastQuality = NULL,
		cex=0.5,
		rad1 = rd$rads1 ,
		rad2 = rd$rads2 ,
		radVals= rd$radVals,
		box=TRUE
		
		)	

	title(main=paste("Forecast size=2^",fcstSize, sep="" ))

} # next inum

dev.off()

} # next fcst size

} #next dufCol 


testFileName <- "texOutput"
for (dufCol in c(5,17,20)){
	scoreValsS <- scoreVals[scoreVals$initialisation == dufCol,]

	meanVals <- aggregate(
		scoreValsS$scoreValue, 
		by=list(
			seedIndex = scoreValsS$seedIndex, 
			fcstSize = scoreValsS$fcstSize,
			fcstPos = scoreValsS$fcstPos, 
			vPower = scoreValsS$vPower, 
			sigmaM = scoreValsS$sigmaM, 
			scoreType = scoreValsS$scoreType
			),
		FUN = mean999)

	colnames(meanVals)[7] <- "scoreValue"

	minVals <- getMinVals19(meanVals)	
	minVals <- cbind(minVals, fcstID = paste(minVals$fcstSize, minVals$xmin, sep="_"))
#	forecastQuality <- sqlFetch(chnl, "forecastQuality")
#	forecastQuality <- cbind(forecastQuality, fcstID = paste(forecastQuality$fcstSize, forecastQuality$sigmaM, sep="_") )

	minVals <- minVals[!duplicated(minVals),]

	for (fcstSize in 5:7){
		scoreResults <- minVals[minVals$fcstSize==fcstSize,]
		

		rowDescription <- paste("dufCol = ",dufCol, "  fcstSize = ", fcstSize, sep="")
		sN <- unique(as.character(scoreResults$scoreType))
		print(sN)
		
		getTestTex(
			rowDescription = rowDescription,
			sigmaU = 0.1,
			scoreResults = scoreResults,
			file=testFileName,
			sN = sN,
			sNPretty = sN)
	} # next fcst size

} #next dufCol 



testFileName <- "texOutputAllDuff"

	scoreValsS <- scoreVals

	meanVals <- aggregate(
		scoreValsS$scoreValue, 
		by=list(
			seedIndex = scoreValsS$seedIndex, 
			fcstSize = scoreValsS$fcstSize,
			fcstPos = scoreValsS$fcstPos, 
			vPower = scoreValsS$vPower, 
			sigmaM = scoreValsS$sigmaM, 
			scoreType = scoreValsS$scoreType,
			initialisation = scoreVals$initialisation),
		FUN = mean999)

	colnames(meanVals)[8] <- "scoreValue"

	minVals <- getMinVals19(meanVals)	
	minVals <- cbind(minVals, fcstID = paste(minVals$fcstSize, minVals$xmin, sep="_"))
#	forecastQuality <- sqlFetch(chnl, "forecastQuality")
#	forecastQuality <- cbind(forecastQuality, fcstID = paste(forecastQuality$fcstSize, forecastQuality$sigmaM, sep="_") )

	minVals <- minVals[!duplicated(minVals),]

	for (fcstSize in 5:7){
		scoreResults <- minVals[minVals$fcstSize==fcstSize,]
		

		rowDescription <- paste("  fcstSize = ", fcstSize, sep="")
		sN <- unique(as.character(scoreResults$scoreType))
		print(sN)
		
		getTestTex(
			rowDescription = rowDescription,
			sigmaU = 0.1,
			scoreResults = scoreResults,
			file=testFileName,
			sN = sN,
			sNPretty = sN)
	} # next fcst size








meanVals <- aggregate(
		scoreVals$scoreValue, 
		by=list(
			seedIndex = scoreVals$seedIndex, 
			fcstSize = scoreVals$fcstSize,
			fcstPos = scoreVals$fcstPos, 
			vPower = scoreVals$vPower, 
			sigmaM = scoreVals$sigmaM, 
			scoreType = scoreVals$scoreType,
			initialisation = scoreVals$initialisation
			),
		FUN = mean999)

	colnames(meanVals)[8] <- "scoreValue"

	minVals <- getMinVals19(meanVals)	
	minVals <- cbind(minVals, fcstID = paste(minVals$fcstSize, minVals$xmin, sep="_"))
#	forecastQuality <- sqlFetch(chnl, "forecastQuality")
#	forecastQuality <- cbind(forecastQuality, fcstID = paste(forecastQuality$fcstSize, forecastQuality$sigmaM, sep="_") )

	minVals <- minVals[!duplicated(minVals),]




jFac <- 3
fcstSizes <- unique(scoreVals$fcstSize)
for (fcstSize in fcstSizes){
	
xRange <- range(minVals$xmin)
pdf(paste("plotSigmaByScore20_fcstSizeSeed",fcstSize  ,".pdf", sep=""))
par(mfrow=c(3,3))
sNs <- as.character(unique(minVals$scoreType))
for (sN in sNs){

mV_sN <- minVals[minVals$scoreType == sN & minVals$fcstSize == fcstSize,]
inis <- sort(unique(mV_sN$initialisation))
numIni <- length(inis)

plot(x=xRange, y=c(0,numIni+1), type="n", xlab="Sigma", ylab="Duffing initialisation", axes=FALSE)
title(main=sN)
axis(1)
axis(2, at=seq(along=inis), labels= inis, las=1)
abline(v=trueParm <- 0.1, col="grey")

for (j in seq(along=inis)){
	xVals <-mV_sN[mV_sN$initialisation == inis[j],"xmin"] 
	seedIndex <- mV_sN[mV_sN$initialisation == inis[j],"seedIndex"] 
	abline(h=j, col="grey")
	#points(x=xVals, y=rep(j, length(xVals)),  col=(j/(numIni-1) - floor(j/(numIni-1)))*(numIni-1) +1, cex=0.5)
	
	text(x=jitter(xVals, factor=jFac), y=jitter(rep(j, length(xVals)), factor=jFac), labels=seedIndex,  col=1:length(xVals), cex=0.5)
	
}
}
dev.off()

}




	

