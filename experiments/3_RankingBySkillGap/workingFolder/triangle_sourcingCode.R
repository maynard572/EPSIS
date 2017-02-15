



# USER INPUTS

invPath <- "/Users/trevormaynard/Dropbox/Phd/modelOffice/experiments/reports/Investigations4"
runtimeGlobals <- list()
runtimeGlobals$dsn <-"ModelOffice5"    
runtimeGlobals$database <-"ModelOffice5"
runtimeGlobals$modelVersion <- "2.6"
expPath <- "/Users/trevormaynard/Dropbox/Phd/modelOffice/experiments/"
expBatch <- "Investigations4"

#########################################################################################
# READ IN REQUIRED FUNCTIONS AND SET UP WORKING DIRECTORY
require(RODBC)
getRFiles<-function(wd){
setwd(wd)
srcFiles<-dir()
for (file in srcFiles){
	source(file)
	}
}

addFwdSlash <- function(str){
	L <- nchar(str)
	if (!(substr(str,start=L,stop=L) == "/")) str <- paste(str, "/", sep="")
	return(str)
	}


# Read in score functions first
setwd("/Users/trevormaynard/Dropbox/Phd/robustIndex/SkillScores/requiredFunctions")
source("requiredFunctions.r")


# Now read in model office functions
invPath <- addFwdSlash(invPath)
expPath <- addFwdSlash(expPath)
R_WD <- paste(invPath, "R", sep="")
TexAndGraphics_WD <- paste(invPath, "TexAndGraphics", sep="")   # the directory into which all graphics are to be saved - and in which the main <report>.tex file is saved
getRFiles(paste("/Users/trevormaynard/Dropbox/Phd/modelOffice/modelVersions/model", runtimeGlobals$modelVersion, sep=""))
getRFiles("/Users/trevormaynard/Dropbox/Phd/modelOffice/resultsViewingFunctions/Model2")
getRFiles(R_WD)  


# Now read in forecastTest functions
setwd("/Users/trevormaynard/Dropbox/Phd/modelOffice/experiments/reports/inv4_forecastTests")
source("extractForecastFromRun_RequiredFunctions.r")


setwd("/Users/trevormaynard/Dropbox/Phd/modelOffice/Triangle")
source("triangle_RequiredFunctions.r")



