##############################################################################################################################
##R CODE FOR THE BRASS RELATIONAL GOMPERTZ MODEL OF FERTILITY
##
##EDDIE HUNSINGER, FEBRUARY 2011 (LAST UPDATED DECEMBER 2018)
##http://www.demog.berkeley.edu/~eddieh/
##
##IF YOU WOULD LIKE TO USE, SHARE OR REPRODUCE THIS CODE, BE SURE TO CITE THE SOURCE
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

plot(fx$SampleFx,type="l",col="blue",lwd=5)
mtext(line=-14,text="Sample fx",font=2,cex=1,col="blue")
points(fx$StandardFx,type="l",col="red",lwd=5)
mtext(line=-10,text="Standard fx",font=2,cex=1,col="red")

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

points(fxAdjusted,type="l",col="purple",lwd=5)
mtext(line=-12,text="Standard fx Adjusted",font=2,cex=1,col="purple")

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

plot(fx$StandardFx,type="l",col="purple",lwd=5)
mtext(line=-12,text="Standard fx",font=2,cex=1,col="purple")
points(fxHighAlpha,type="l",col="red",lwd=5)
mtext(line=-10,text="High (0.5) Alpha fx",font=2,cex=1,col="red")
points(fxLowAlpha,type="l",col="blue",lwd=5)
mtext(line=-14,text="Low (-0.5) Alpha fx",font=2,cex=1,col="blue")

plot(fx$StandardFx,type="l",col="purple",lwd=5)
mtext(line=-12,text="Standard fx",font=2,cex=1,col="purple")
points(fxHighBeta,type="l",col="red",lwd=5)
mtext(line=-10,text="High (1.5) Beta fx",font=2,cex=1,col="red")
points(fxLowBeta,type="l",col="blue",lwd=5)
mtext(line=-14,text="Low (0.5) Beta fx",font=2,cex=1,col="blue")

#write.table(###, file="G:/###/###.csv", sep=",")

