##############################################################################################################################
##R CODE FOR THE BRASS RELATIONAL GOMPERTZ MODEL OF FERTILITY
##
##EDDIE HUNSINGER, FEBRUARY 2011 (LAST UPDATED MARCH 2022)
##http://www.demog.berkeley.edu/~eddieh/
##
##EXAMPLE DATA IS LINKED, SO YOU SHOULD BE ABLE TO SIMPLY COPY ALL AND PASTE INTO R
##
##THERE IS NO WARRANTY FOR THIS CODE
##THIS CODE HAS NOT BEEN TESTED AT ALL-- PLEASE LET ME KNOW IF YOU FIND ANY PROBLEMS (edyhsgr@gmail.com)
##############################################################################################################################

##############################################################################################################################
#STEP 1: Read in and plot the fertility data (SampleFx is collected data, StandardFx is data to be used as the standard-- 
##so the goal is to adjust standard data (StandardFx) to meet the earliness and width of SampleFx)
##############################################################################################################################

fx<-read.table(file="https://github.com/AppliedDemogToolbox/Hunsinger_BrassRelationalGompertz/raw/master/Fx.csv",header=TRUE,sep=",")
agegrouplabels<-c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95-99","100+")

plot(fx$SampleFx,type="l",,ylim=c(0,.42),
main="Brass relational Gompertz model demonstration",
		xlab="age",
		ylab="fertility rate, proportional (sums to 1)",
		col="blue",lwd=5,axes=FALSE)
	axis(side=2,las=2)
	axis(side=1,at=1:length(fx$x),
		labels=agegrouplabels,cex.axis=.75,las=2)
points(fx$StandardFx,type="l",lty=2,col="red",lwd=5)
legend(x="topright",legend=c("sample","standard"),
       col=c("blue","red"),lty=c(1,2),lwd=c(5,5),bg="white")

##############################################################################################################################
#STEP 2: Enter the Brass Relational Gompertz function
##############################################################################################################################

#Note: minagegroupfert is the first age group that is non-zero, and maxagegroupfert is the last age group that is non-zero
Brass<-function(fx,fxBase,minagegroupfert,maxagegroupfert){
cumsumfx<-(cumsum(fx))
cumsumfxBase<-(cumsum(fxBase))
Yx<-log(-log(cumsumfx))
YxBase<-log(-log(cumsumfxBase))
Estimate<-lm(Yx[minagegroupfert:(maxagegroupfert-1)]~YxBase[minagegroupfert:(maxagegroupfert-1)])
Alpha<-Estimate$coefficients[1]
Beta<-Estimate$coefficients[2]
Brass<- data.frame(Beta=Beta,Alpha=Alpha)
return(Brass)}

##############################################################################################################################
#STEP 3: Using the Brass function, estimate the coefficients for adjusting the standard Fx data
##############################################################################################################################

BrassCoefficients<-Brass(fx$SampleFx,fx$StandardFx,3,10)
BrassCoefficients

##############################################################################################################################
#STEP 4: Adjust StandardFx with the estimated coefficients and plot the adjusted fx
##############################################################################################################################

cumsumfxStandard<-(cumsum(fx$StandardFx))
YxStandard<-log(-log(cumsumfxStandard))
YxStandardAdjusted<-NULL
for (i in 1:length(YxStandard)){YxStandardAdjusted[i]<-BrassCoefficients$Beta*YxStandard[i]+BrassCoefficients$Alpha}

cumsumfxAdjusted<-NULL
for (i in 1:length(YxStandard)){cumsumfxAdjusted[i]<-exp(-exp(YxStandardAdjusted[i]))}
fxAdjusted<-NULL
for (i in 2:length(YxStandard)){fxAdjusted[i]<-cumsumfxAdjusted[i]-cumsumfxAdjusted[i-1]}
fxAdjusted[1]<-0

points(fxAdjusted,type="l",lty=3,col="purple",lwd=5)
legend(x="topright",legend=c("sample","standard","adjusted standard"),
       col=c("blue","red","purple"),lty=c(1,2,3),lwd=c(5,5,5),bg="white")

##############################################################################################################################
#STEP 5: Plot some effects of change in Brass Alpha and Brass Beta
##############################################################################################################################

YxStandardHighAlpha<-NULL
for (i in 1:length(YxStandard)){YxStandardHighAlpha[i]<-(1)*YxStandard[i]+(.5)}
cumsumfxHighAlpha<-NULL
for (i in 1:length(YxStandard)){cumsumfxHighAlpha[i]<-exp(-exp(YxStandardHighAlpha[i]))}
fxHighAlpha<-NULL
for (i in 2:length(YxStandard)){fxHighAlpha[i]<-cumsumfxHighAlpha[i]-cumsumfxHighAlpha[i-1]}
fxHighAlpha[1]<-0

YxStandardLowAlpha<-NULL
for (i in 1:length(YxStandard)){YxStandardLowAlpha[i]<-(1)*YxStandard[i]+(-.5)}
cumsumfxLowAlpha<-NULL
for (i in 1:length(YxStandard)){cumsumfxLowAlpha[i]<-exp(-exp(YxStandardLowAlpha[i]))}
fxLowAlpha<-NULL
for (i in 2:length(YxStandard)){fxLowAlpha[i]<-cumsumfxLowAlpha[i]-cumsumfxLowAlpha[i-1]}
fxLowAlpha[1]<-0

YxStandardHighBeta<-NULL
for (i in 1:length(YxStandard)){YxStandardHighBeta[i]<-(1.5)*YxStandard[i]+(0)}
cumsumfxHighBeta<-NULL
for (i in 1:length(YxStandard)){cumsumfxHighBeta[i]<-exp(-exp(YxStandardHighBeta[i]))}
fxHighBeta<-NULL
for (i in 2:length(YxStandard)){fxHighBeta[i]<-cumsumfxHighBeta[i]-cumsumfxHighBeta[i-1]}
fxHighBeta[1]<-0

YxStandardLowBeta<-NULL
for (i in 1:length(YxStandard)){YxStandardLowBeta[i]<-(.5)*YxStandard[i]+(0)}
cumsumfxLowBeta<-NULL
for (i in 1:length(YxStandard)){cumsumfxLowBeta[i]<-exp(-exp(YxStandardLowBeta[i]))}
fxLowBeta<-NULL
for (i in 2:length(YxStandard)){fxLowBeta[i]<-cumsumfxLowBeta[i]-cumsumfxLowBeta[i-1]}
fxLowBeta[1]<-0

plot(fx$StandardFx,type="l",,ylim=c(0,.42),
main="Brass relational Gompertz model demonstration",
		xlab="age",
		ylab="fertility rate, proportional (sums to 1)",
		col="purple",lwd=5,axes=FALSE)
	axis(side=2,las=2)
	axis(side=1,at=1:length(fx$x),
		labels=agegrouplabels,cex.axis=.75,las=2)
points(fxHighAlpha,type="l",lty=2,col="red",lwd=5)
points(fxLowAlpha,type="l",lty=3,col="blue",lwd=5)
legend(x="topright",legend=c("standard","standard with 0.5 alpha","standard with -0.5 alpha"),
       col=c("purple","red","blue"),lty=c(1,2,3),lwd=c(5,5,5),bg="white")

plot(fx$StandardFx,type="l",ylim=c(0,.42),
main="Brass relational Gompertz model demonstration",
		xlab="age",
		ylab="fertility rate, proportional (sums to 1)",
		col="purple",lwd=5,axes=FALSE)
	axis(side=2,las=2)
	axis(side=1,at=1:length(fx$x),
		labels=agegrouplabels,cex.axis=.75,las=2)
points(fxHighBeta,type="l",lty=2,col="red",lwd=5)
points(fxLowBeta,type="l",lty=3,col="blue",lwd=5)
legend(x="topright",legend=c("standard","standard with 0.5 beta","standard with -0.5 beta"),
       col=c("purple","red","blue"),lty=c(1,2,3),lwd=c(5,5,5),bg="white")

#write.table(###, file="G:/###/###.csv", sep=",")

