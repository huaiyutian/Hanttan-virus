rm(list=ls())

#  select other packages
packages <- c("data.table", "tidyverse", "sf", "sp", "spdep","readxl",
              "dlnm", "tsModel", "hydroGOF","RColorBrewer", 
              "geofacet", "ggpubr", "ggthemes", "readxl", "dplyr", "mice","zoo", "splines")

# install.packages
# lapply(packages, install.packages, character.only = TRUE)

# load packages
lapply(packages, library, character.only = TRUE)

# library(hydroGOF)
# library(hydroTSM)
# 
# if (!require(devtools)) install.packages("devtools")
# library(devtools)
# install_github("hzambran/hydroTSM")
# 
# #hydroGOF
# package_url="https://cran.r-project.org/src/contrib/Archive/hydroGOF/hydroGOF_0.4-0.tar.gz"
# install.packages(package_url, repos = NULL, type = 'source')

source("D:/W/W/Hantaan/new/Hantaan fig file/data/fig3/clean_data.R")
##############################################################################

setwd("D:/W/W/Hantaan/new/Hantaan fig file/data/fig3/")
monthdata <- read_xlsx("Hu 2023.xlsx",sheet="month")

PR <- data.clean("monthdata")

PR <- data.process(PR)

PR <- data.control(PR)

set.seed(123)

cb.temp = crossbasis(PR$avg_temp, lag=24, argvar=list(fun="ns", knots = equalknots(PR$avg_temp, fun="ns", df=4)), arglag=list(fun="ns"))

cb.rainfall = crossbasis(PR$rainfall, lag=24, argvar=list(fun="ns", knots = equalknots(PR$rainfall, fun="ns", df=4)), arglag=list(fun="ns"))

cb.avg_patch_size = crossbasis(PR$avg_patch_size, lag=48, argvar=list(fun="ns", knots = equalknots(PR$avg_patch_size, fun="ns", df=4)), arglag=list(fun="ns"))

cb.RNdensity = crossbasis(PR$RNdensity, lag=24, argvar=list(fun="ns", knots = equalknots(PR$RNdensity, fun="ns", df=4)), arglag=list(fun="ns"))

cb.AAdensity = crossbasis(PR$AAdensity, lag=24, argvar=list(fun="ns", knots = equalknots(PR$AAdensity, fun="ns", df=4)), arglag=list(fun="ns"))

cb.RFdensity = crossbasis(PR$RFdensity, lag=24, argvar=list(fun="ns", knots = equalknots(PR$RFdensity, fun="ns", df=4)), arglag=list(fun="ns"))

model1 = glm(PR$AAdensity ~ cb.AAdensity + cb.RFdensity + cb.RNdensity 
             + cb.temp
             + cb.rainfall
             + cb.avg_patch_size
             + jan+feb+mar+apr+may+jun+jul+aug+sep+oct+nov+dec,
             family=gaussian(), PR)
summary(model1)


pred1.temp = crosspred(cb.temp, model1, cen=0, bylag=1)  
pred1.rainfall = crosspred(cb.rainfall, model1, cen=0, bylag=1) 
pred1.avg_patch_size = crosspred(cb.avg_patch_size, model1, cen=round(min(PR$avg_patch_size)), bylag=1)
pred1.AAdensity = crosspred(cb.AAdensity, model1, cen=0, bylag=1) 
pred1.RFdensity = crosspred(cb.RFdensity, model1, cen=0, bylag=1)
pred1.RNdensity = crosspred(cb.RNdensity, model1, cen=0, bylag=1)

dev.off()
par(mar = c(4, 4, 2, 1)) #set margins
par(mfrow=c(1,1))   

# png(file = "single_temp.png", width = 3000, height = 2500, res = 300)
plot(pred1.temp, "contour", xlab="MeanTemp", key.title=title("Effect"),cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title(xlab="avg_temp",ylab="Lag (months)",cex.main=2,cex.lab=1.5))


plot(pred1.rainfall, "contour", xlab="Rainfall", key.title=title("Effect"),cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title(xlab="Rainfall",ylab="Lag (months)",cex.main=2,cex.lab=1.5))


plot(pred1.avg_patch_size, "contour", xlab="avg_patch_size", key.title=title("Effect"),cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title(xlab="avg_patch_size",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.AAdensity, "contour", xlab="AAdensity", key.title=title("Effect"),cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title(xlab="AAdensity",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.RFdensity, "contour", xlab="RFdensity", key.title=title("Effect"),cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title(xlab="RFdensity",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.RNdensity, "contour", xlab="RNdensity", key.title=title("Effect"),cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title(xlab="RNdensity",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

plot(pred1.interaction, "contour", xlab="interaction", key.title=title("Effect"),cex.axis=2,
     plot.axes={axis(1,cex.axis=2)
       axis(2,cex.axis=2)},
     key.axes = axis(4,cex.axis=2),
     plot.title=title(xlab="AAdensity*Patch",ylab="Lag (months)",cex.main=2,cex.lab=1.5))

dev.off()
# par(mfrow = c(2, 2), oma = c(0.5, 0.5, 3, 0.5))
#

par(mar=c(3,4,2,1))
par(mfrow=c(2,3))
for (i in 1:6){
  plot(pred1.temp,"slices",lag=i,cumul=F,
       col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
       cex.axis=1.5,cex.lab=1,cex.main=1,lwd=1,
       ci.level=0.95,col=2,xlab="Avg_temp",
       ylab="Effect on AAdensity")
  mtext(paste("Lag",i,seq=''),side = 3,line =1,cex=0.5,font=2)
}
plot(pred1.temp,"slices",lag=6,cumul=F,
     col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=3,
     ci.level=0.95,col=2,xlab="Avg_temp",
     ylab="Effect on AAdensity")
mtext("Lag6",side = 3,line =1,cex=1.5,font=2)

par(mar=c(3,4,2,1))
par(mfrow=c(2,3))
for (i in 1:6){
  plot(pred1.rainfall,"slices",lag=i,cumul=F,
       col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
       cex.axis=1.5,cex.lab=1,cex.main=1,lwd=1,
       ci.level=0.95,col=2,xlab="Rainfall",
       ylab="Effect on AAdensity")
  mtext(paste("Lag",i,seq=''),side = 3,line =1,cex=0.5,font=2)
}
plot(pred1.rainfall,"slices",lag=1,cumul=F,
     col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=3,
     ci.level=0.95,col=2,xlab="Rainfall",
     ylab="Effect on AAdensity")
mtext("Lag1",side = 3,line =1,cex=1.5,font=2)

par(mar = c(3, 4, 2, 1)) 
par(mfrow=c(3,4))  
for (i in seq(1,12)){
  plot(pred1.avg_patch_size,"slices",lag=i,cumul=F,
       col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
       cex.axis=1,cex.lab=1,cex.main=1.5,lwd=1,
       ci.level=0.95,col=2,xlab="Avg_patch_size",
       ylab="Effect on AAdensity")
  mtext(paste("Lag",i,seq=''),side = 3,line =1,cex=0.5,font=2)
}

dev.off()
par(mar = c(4, 4, 3, 1))
par(mfrow=c(1,1))  
plot(pred1.avg_patch_size,"slices",lag=1,cumul=F,
     col="black", ci.arg=list(col = grey(0.85)),mgp=c(2.5,1,0),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=3,
     ci.level=0.95,col=2,xlab="Avg_patch_size",
     ylab="Effect on AAdensity")
mtext("Lag1",side = 3,line =1,cex=1.5,font=2)

par(mar = c(3, 4, 2, 1))
par(mfrow=c(2,3)) 
for (i in seq(1,6)){
  plot(pred1.AAdensity,"slices",lag=i,cumul=F,
       col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
       cex.axis=1,cex.lab=1,cex.main=1.5,lwd=1,
       ci.level=0.95,col=2,xlab="AA",
       ylab="Effect on AAdensity")
  mtext(paste("Lag",i,seq=''),side = 3,line =1,cex=0.5,font=2)
}
plot(pred1.AAdensity,"slices",lag=2,cumul=F,
     col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=3,
     ci.level=0.95,col=2,xlab="AA",
     ylab="Effect on CTdensity")
mtext("Lag6",side = 3,line =1,cex=1.5,font=2)

par(mar = c(3, 4, 2, 1))
par(mfrow=c(2,3)) 
for (i in seq(1,6)){
  plot(pred1.RFdensity,"slices",lag=i,cumul=F,
       col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
       cex.axis=1,cex.lab=1,cex.main=1.5,lwd=1,
       ci.level=0.95,col=2,xlab="AA",
       ylab="Effect on AAdensity")
  mtext(paste("Lag",i,seq=''),side = 3,line =1,cex=0.5,font=2)
}
plot(pred1.RFdensity,"slices",lag=3,cumul=F,
     col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=3,
     ci.level=0.95,col=2,xlab="CTdensity",
     ylab="Effect on AAdensity")
mtext("Lag3",side = 3,line =1,cex=1.5,font=2)


par(mar = c(3, 4, 2, 1)) 
par(mfrow=c(2,3)) 
for (i in seq(1,6)){
  plot(pred1.RNdensity,"slices",lag=i,cumul=F,
       col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
       cex.axis=1,cex.lab=1,cex.main=1.5,lwd=1,
       ci.level=0.95,col=2,xlab="AA",
       ylab="Effect on AAdensity")
  mtext(paste("Lag",i,seq=''),side = 3,line =1,cex=0.5,font=2)
}
plot(pred1.RNdensity,"slices",lag=1,cumul=F,
     col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=3,
     ci.level=0.95,col=2,xlab="RNdensity",
     ylab="Effect on AAdensity")
mtext("Lag1",side = 3,line =1,cex=1.5,font=2)

par(mar = c(3, 4, 2, 1))
par(mfrow=c(3,4)) 
for (i in seq(1,12)){
  plot(pred1.interaction,"slices",lag=6,cumul=F,
       col="black", ci.arg=list(col = grey(0.85)),mgp=c(3,1,0),
       cex.axis=1,cex.lab=1,cex.main=1.5,lwd=1,
       ci.level=0.95,col=2,xlab="Temp*Patch",
       ylab="Effect on CTdensity")
  mtext(paste("Lag",i,seq=''),side = 3,line =1,cex=0.5,font=2)
}

dev.off()



set.seed(123) 
sample_indices <- sample(1:nrow(PR), size = 0.2 * nrow(PR)) # 20% of the data is used for testing
test_data <- PR[sample_indices, ]
train_data <- PR[-sample_indices, ]

