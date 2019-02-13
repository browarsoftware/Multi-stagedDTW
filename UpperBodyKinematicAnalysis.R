#Starting file for the analysis
#Change path here
path <- "e:\\Publikacje\\eecs 2018\\extended2\\UpperBodyKinematicAnalysis-master\\"

source(paste(path, "renderactor.R", sep =""))
source(paste(path, "rotatedata.R", sep =""))
source(paste(path, "calculatekinematic.R", sep =""))
source(paste(path, "dtwcomparator.R", sep =""))
                            
#source("e:\\Publikacje\\eecs 2018\\extended2\\UpperBodyKinematicAnalysis-master\\rotatedata.R")
#source("e:\\Publikacje\\eecs 2018\\extended2\\UpperBodyKinematicAnalysis-master\\calculatekinematic.R")
#source("e:\\Publikacje\\eecs 2018\\extended2\\UpperBodyKinematicAnalysis-master\\dtwcomparator.R")

################################################
#LOAD DATA
################################################

#refdata <- read.csv(paste(path, "data\\baerhug0.bvh.csv", sep =""))
#inputdata <- read.csv(paste(path, "data\\baerhug1.bvh.csv", sep =""))

#refdata <- read.csv(paste(path, "data\\dumbbellslifting0.bvh.csv", sep =""))
#inputdata <- read.csv(paste(path, "data\\dumbbellslifting1.bvh.csv", sep =""))


#refdata <- read.csv(paste(path, "data\\rightshouldertouch0.bvh.csv", sep =""))
#inputdata <- read.csv(paste(path, "data\\rightshouldertouch1.bvh.csv", sep =""))


#refdata <- read.csv(paste(path, "data\\risingrightarmfrontally0.bvh.csv", sep =""))
#inputdata <- read.csv(paste(path, "data\\risingrightarmfrontally1.bvh.csv", sep =""))


refdata <- read.csv(paste(path, "data\\righthandbetweenshoulders0.bvh.csv", sep =""))
inputdata <- read.csv(paste(path, "data\\righthandbetweenshoulders1.bvh.csv", sep =""))


#refdata <- read.csv(paste(path, "data\\shoulderelevation0.bvh.csv", sep =""))
#inputdata <- read.csv(paste(path, "data\\shoulderelevation1.bvh.csv", sep =""))



#left or right, depending of foot that remains on the ground, if you kick with right leg it should be left

leg <- "Left"
#leg <- "Right"

inputdataalignment <- rotatedata(inputdata, refdata, paste(leg, "Shoulder.Dx", sep = ""),paste(leg, "Shoulder.Dz", sep = ""),
                                 "RightShoulder.Dx","RightShoulder.Dz")
inputdataalignmentkinematic <- calculatekinematic(inputdataalignment, paste(leg, "Shoulder", sep = ""))
refdatakinematic <- calculatekinematic(refdata, paste(leg, "Shoulder", sep = ""))
inputdataalignmentkinematic <- aligninputandrefdata(inputdataalignmentkinematic, refdatakinematic, limbname = paste(leg, "Shoulder", sep = ""))
inputdataalignmentkinematicf <- generateFeatures(inputdataalignmentkinematic)
refdatakinematicf <- generateFeatures(refdatakinematic)



################################################
#BEGIN ANALYSIS
################################################
extremumtreshold <- 0.66
smoothSize <- 0.1
analyzerange <- ceiling(nrow(inputdataalignment) * smoothSize * 1) #sprawdzia, jaki powinien bya ten przedzia3!

footddf <- analyzedta(refdatakinematic, inputdataalignmentkinematic,
                      refdatakinematicf$RightHand, inputdataalignmentkinematicf$RightHand, FUN=euc.dist, smoothSize = smoothSize, 
                      c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), plottitle = "RightHand D", extremumtreshold = extremumtreshold, smoothSizeHelper = 0.1, 
                      whattodraw = "RightFoot", plotrgl = TRUE)

#plotsmoothingresults(footddf,"footddf", plotifnoextreams = TRUE, plotsmoothed = TRUE)

dftoanalyze <- tresholdresults(footddf, extremumtreshold)


plot(footddf$path1, footddf$path2, xlab = "Reference signal [sample id]", ylab = "Input signal [sample id]", main ="DTW alignment", type='l')


plot(footddf$data, xlab = "Time [10^-2 s]", ylab = "Distance [cm]", main ="DTW alignment function", col = 'black', type='l', ylim = c(min(footddf$data), max(footddf$data) * 1.5))
#lines(footddf$data, xlab = "Time [s^-100]", ylab = "Distance [cm]", main ="DTW alignment function")
lines(footddf$smoothdata, col = 'purple', lty = 2)
#points(footddf$smoothdata, col = 'blue')
for (a in 1:length(footddf$extremumid))
{
  points(footddf$extremumid[a], footddf$smoothdata[footddf$extremumid[a]], col = "cyan",pch = 1, cex = 2, lwd=3)
}
for (a in 1:length(dftoanalyze$extremumid))
{
  points(dftoanalyze$extremumid[a], dftoanalyze$smoothdata[dftoanalyze$extremumid[a]], col = "red",pch = 4, cex = 2, lwd=3)
}
legend(x= "topright", y=max(footddf$data) * 1.5, legend=c("Original", "Smoothed", "Maxima", "Maxima over treshold"), col=c("black", "purple", 'cyan','red'), lty=c(1,2,1), cex=0.8)


#this might render quite long, but looks nice :-)
rglplotanalyzedata(refdatakinematic, inputdataalignmentkinematic, refdatakinematicf$RightHand, inputdataalignmentkinematicf$RightHand, footddf$path1, footddf$path2,resultdata =   dftoanalyze, whattodraw = 'RightHand')

#right hand
summary(footddf$data)

#####################


plotsmoothingresults(dftoanalyze,"DTWaf of right hand trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE)
###################Shoulder
RightShoulderX <-analyzedta(refdatakinematicf, inputdataalignmentkinematicf, refdatakinematicf$ListRightShoulderX, inputdataalignmentkinematicf$ListRightShoulderX, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                   plottitle = "DTWaf of right shoulder X angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "RightShoulderX")
RightShoulderX <- tresholdresults(RightShoulderX, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightShoulderX, analyzerange),"DTWaf of right RightShoulder X trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

RightShoulderY <-analyzedta(refdatakinematicf, inputdataalignmentkinematicf, refdatakinematicf$ListRightShoulderY, inputdataalignmentkinematicf$ListRightShoulderY, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                            plottitle = "DTWaf of right shoulder Y angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "RightShoulderY")
RightShoulderY <- tresholdresults(RightShoulderY, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightShoulderY, analyzerange),"DTWaf of right RightShoulder Y trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

RightShoulderZ <-analyzedta(refdatakinematicf, inputdataalignmentkinematicf, refdatakinematicf$ListRightShoulderZ, inputdataalignmentkinematicf$ListRightShoulderZ, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                            plottitle = "DTWaf of right shoulder Z angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "RightShoulderZ")
RightShoulderZ <- tresholdresults(RightShoulderZ, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightShoulderZ, analyzerange),"DTWaf of right RightShoulder Z trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

summary(RightShoulderX$data)
summary(RightShoulderY$data)
summary(RightShoulderZ$data)

#########################Arm
RightArmX <-analyzedta(refdatakinematicf, inputdataalignmentkinematicf, refdatakinematicf$ListRightArmX, inputdataalignmentkinematicf$ListRightArmX, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                            plottitle = "DTWaf of right Arm X angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "RightArmX")
RightArmX <- tresholdresults(RightArmX, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightArmX, analyzerange),"DTWaf of right RightArm X trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

RightArmY <-analyzedta(refdatakinematicf, inputdataalignmentkinematicf, refdatakinematicf$ListRightArmY, inputdataalignmentkinematicf$ListRightArmY, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                            plottitle = "DTWaf of right Arm Y angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "RightArmY")
RightArmY <- tresholdresults(RightArmY, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightArmY, analyzerange),"DTWaf of right RightArm Y trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

RightArmZ <-analyzedta(refdatakinematicf, inputdataalignmentkinematicf, refdatakinematicf$ListRightArmZ, inputdataalignmentkinematicf$ListRightArmZ, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                            plottitle = "DTWaf of right Arm Z angle analysis",  path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, whattodraw = "RightArmZ")
RightArmZ <- tresholdresults(RightArmZ, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, RightArmZ, analyzerange),"DTWaf of right RightArm Z trajectory analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

#xxxx <- data.frame(unlist(inputdataalignmentkinematicf$ListRightArmZ), 
#                   unlist(inputdataalignmentkinematicf$ListRightArmY))
#yyyy <- data.frame(unlist(refdatakinematicf$ListRightArmZ), 
#                   unlist(refdatakinematicf$ListRightArmY))

summary(RightArmX$data)
summary(RightArmY$data)
summary(RightArmZ$data)

##################Elbow

elbowdf <- analyzedta(refdatakinematic, inputdataalignmentkinematic,refdatakinematicf$ListRightElbow, inputdataalignmentkinematicf$ListRightElbow, FUN=euc.dist1d, smoothSize = smoothSize, c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                      plottitle = "DTWaf of right knee angle analysis", path1 = footddf$path1, path2 = footddf$path2, extremumtreshold = extremumtreshold, smoothSizeHelper = 0.1, whattodraw = "RightElbow")
elbowdf <- tresholdresults(elbowdf, extremumtreshold)
plotsmoothingresults(alignextremum(dftoanalyze, elbowdf, analyzerange),"DTWaf of right elbow angle analysis", plotifnoextreams = TRUE, plotsmoothed = TRUE, ylab = "Angle [rad]")

summary(elbowdf$data)

######################################################################
