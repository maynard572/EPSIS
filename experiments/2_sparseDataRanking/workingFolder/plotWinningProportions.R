plotWinningProportions <- function(scrRes,  scoreType, scoreParam=99999){

	resHigh <- scrRes[scrRes$fcstType == "high"  & scrRes$scoreType == scoreType & 
		scrRes$scoreParam == scoreParam,]
	resMed <- scrRes[scrRes$fcstType == "med"  & scrRes$scoreType == scoreType & 
		scrRes$scoreParam == scoreParam,]
	resLow <- scrRes[scrRes$fcstType == "low"  & scrRes$scoreType == scoreType & 
		scrRes$scoreParam == scoreParam,]

	plot(resLow$N, resLow$winningProportions, ylim=c(0,1), type="n", pch=19, cex=0.5, axes=FALSE, ann=FALSE)
	axis(side=2, las=2)
	axis(side=1)
	
	sPnm <- ifelse(scoreParam == 99999, "", scoreParam)
	
	
	title(
		main = paste(scoreType, " ", sPnm, sep=""), 
		ylab=bquote(F[1]+F[2]+F[3]),  
		xlab="log2(N)",
		sub="where N is number of observations")

	polygon(
		x =  c(resLow$N, rev(resLow$N), resLow$N[1]),
		y =  c(rep(0, length(resLow$N)), rev(resLow$winningProportions), 0), col="red", border=NA)

	polygon(
		x =  c(resLow$N, rev(resLow$N), resLow$N[1]),
		y =  c(resLow$winningProportions, rev(resLow$winningProportions+resMed$winningProportions), 0), 
			col="blue", border=NA)

	polygon(
		x =  c(resLow$N, rev(resLow$N), resLow$N[1]),
		y =  c(resLow$winningProportions+resMed$winningProportions, 
			rev(resLow$winningProportions+resMed$winningProportions + resHigh$winningProportions), 0), 
				col="green", border=NA)


	cXLow <- max(min(sum(resLow$winningProportions), 3.5), 1.1)
	cXMed <- max(min(sum(resMed$winningProportions), 3.5), 1.1)
	cXHigh <- max(min(sum(resHigh$winningProportions), 3.5), 1.1)

	cYLow <- 0.3*resLow$winningProportions[1]
	cYMed <- resLow$winningProportions[1] + 0.5 * resMed$winningProportions[1]
	cYHigh <- resLow$winningProportions[1] + resMed$winningProportions[1] + 0.6 *  resHigh$winningProportions[1] 	


	text(x=cXLow, y=cYLow, labels="N(0, 1/sqrt(2))", pos=4)
	text(x=cXMed, y=cYMed, labels="N(0, 1)", col="white", pos=4)
	text(x=cXHigh, y=cYHigh, labels="N(0, sqrt(2)", pos=4)


}
