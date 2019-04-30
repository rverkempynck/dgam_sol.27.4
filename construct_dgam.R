
rm(list = ls())

library(mapdata);
library(sp)
library(surveyIndex)

cmSize <- 1
years <- 2010:2018
# outFolder <- "."
genus <- "Solea"
bfamily <- "solea"
datafile <- "sole.RData"

## extracted dataset
load("final_dataset.RData")

## Area subsetting (column name "ICES_SUB" may depend on version of shape file)
# names(d[["HH"]])
# d <- subset(d, ICES_SUB %in% c("IVa","IVb","IVc"))
# d <- subset(d, SubArea %in% c("IVa","IVb","IVc"))
bubblePlot(d,scale=1/10)

## Do we need gear effects? (more than one?)
summary(d$Gear)

## impute missing depths
summary(d$Depth)
dmodel <- gam(log(Depth) ~ s(lon,lat,k=200),data=d[["HH"]])
selQ1c <- subset(d,is.na(Depth))
selQ1c$Depth <- 0; ## Guard against NA-error
d$Depth[is.na(d$Depth)] <- exp(predict(dmodel,newdata=selQ1c[[2]]))

## Check for enough age data in all years
xtabs(NoAtALK~Year+Age,data=d[["CA"]])
ages <- 0:7

####################
## Age-length key 
####################
mc.cores <- 2 ## Windows users should use mc.cores <- 1
## Declare settings for ALK model
mf  <-  NULL
ack <- TRUE;
useBICs <- TRUE;
varCofs <- FALSE;
maxKs <- 50;

d.before <- d

## If there are years where the youngest age group is not observed, the following is needed:
for(aa in ages){
  d  <- fixAgeGroup(d,age=aa,n=1,fun=ifelse(aa==min(ages),min,mean))
}

add.ALK<-function(d,ages){
  d[[1]] <- subset(d[[1]],Age>=min(ages))
  d <- addSpectrum(d,by=cmSize)
  ## split and fit by year
  d.ysplit  <-  split(d,d$Year)
  d.ALK <-  mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs,mc.cores=mc.cores)
  ## predict numbers-at-age
  d.Nage <- mclapply(d.ALK,predict,mc.cores=mc.cores)
  for(i in 1:length(d.ALK)) d.ysplit[[i]]$Nage=d.Nage[[i]];
  ## merge back into one object
  dd <- do.call("c",d.ysplit)
  dd    
}

d  <- add.ALK(d,ages)

## Make ctime : a numeric time variable 
d$ctime  <-  as.numeric(as.character(d$Year))

#######################
## Model
#######################

## Number of knots for spatial splines 
kvecP  <-  c(16,rep(100,length(ages)-1))
kvecZ  <-  kvecP/2

## Time-invariant spatial effect
modelsStatZ <- rep("Year+Gear+s(lon,lat,bs=c('tp'),k=kvecZ[a])+s(Depth,bs='ts',k=6)+offset(log(HaulDur))+s(TimeShotHour,k=5,bs='cc')",length(ages))
modelsStatP <- rep("Year+Gear+s(lon,lat,bs=c('tp'),k=kvecP[a])+s(Depth,bs='ts',k=6)+offset(log(HaulDur)) + s(TimeShotHour,k=5,bs='cc')",length(ages))

## Time-varying spatial effect
modelsNonStat <- rep("Year+Gear+te(ctime,lon,lat,d=c(1,2),bs=c('cs','tp'),k=c(4,40))+s(Depth,bs='ts',k=6)+offset(log(HaulDur))+s(TimeShotHour,k=5,bs='cc')",length(ages))

tknots  <- list(TimeShotHour=seq(0,24,length=5))

grid <- getGrid(d,nLon=30)
plot(grid[[1]],grid[[2]])

## Fit models
SI  <-  getSurveyIdx(d,ages=ages,myids=grid[[3]],cutOff=0.1,fam=rep("LogNormal",length(ages)),mc.cores=mc.cores,modelZ=modelsStatZ,modelP=modelsStatP,kvecP=kvecP,kvecZ=kvecZ,knotsZ=tknots,knotsP=tknots)

SI.ns  <-  getSurveyIdx(d,ages=ages,myids=grid[[3]],cutOff=0.1,fam=rep("LogNormal",length(ages)),mc.cores=mc.cores,modelZ=modelsNonStat,modelP=modelsNonStat,kvecP=kvecP,kvecZ=kvecZ,knotsZ=tknots,knotsP=tknots)

SI.ns.tw  <-  getSurveyIdx(d,ages=ages,myids=grid[[3]],cutOff=0.1,fam=rep("Tweedie",length(ages)),mc.cores=mc.cores,modelZ=modelsNonStat,modelP=modelsNonStat,kvecP=kvecP,kvecZ=kvecZ,knotsZ=tknots,knotsP=tknots)


final.index = SI.ns.tw

SI.alt = SI.ns$idx ## getSurveyIdxStratMean(d,ages)


## AIC / BIC
surveyIndex:::AIC.surveyIdx( SI )
surveyIndex:::AIC.surveyIdx( SI.ns )
surveyIndex:::AIC.surveyIdx( SI.ns.tw )

surveyIndex:::AIC.surveyIdx( SI, BIC=TRUE )
surveyIndex:::AIC.surveyIdx( SI.ns, BIC=TRUE )
surveyIndex:::AIC.surveyIdx( SI.ns.tw, BIC=TRUE )

mean(internalCons(SI$idx))
mean(internalCons(SI.ns$idx))
mean(internalCons(SI.ns.tw$idx))


## Plots

surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,1,1,1)),select=c("index"),plotByAge=FALSE)

## 2nd smooth term (depth)
surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,3,1,1)),select=c("2"),plotByAge=FALSE)

## 3rd smooth term (time of day)
surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,3,1,1)),select=c("3"),plotByAge=FALSE,ylim=c(-2 ,2),xlim=c(0,24))
dev.off()

## residuals
surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,3,1,1)),select=c("residuals"),plotByAge=FALSE)
surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,3,1,1)),select=c("fitVsRes"),plotByAge=FALSE)
surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,3,1,1)),select=c("resVsYear"),plotByAge=FALSE)
surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,3,1,1)),select=c("resVsShip"),plotByAge=FALSE)

## Map from stationary model
surveyIdxPlots(SI,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,1,1,1)),select=c("map"),plotByAge=FALSE,colors=rev(heat.colors(5)))


## maps from non-stationary model
pdf("maps.pdf",onefile=TRUE)
for(yy in years){
  surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,1,1,1)),select=c("map"),plotByAge=FALSE,colors=rev(heat.colors(5)),legend=FALSE,year=yy)
}
dev.off()

## Spatial residuals (positive model only)
surveyIdxPlots(final.index,d,alt.idx=SI.alt,myids=grid[[3]],par=list(mfrow=c(3,3),mar=c(4,3,1,1)),select=c("spatialResiduals"),plotByAge=FALSE,year=2014)

## Internal consistency plot using FLCore
##remotes::install_github("flr/FLCore")
library(FLCore)
tofli<-function(x){
  fli <- FLQuant(t(x),dimnames=list(age=colnames(x),year=rownames(x)))
  fli <- FLIndex(index=fli)
}
plot(tofli(final.index$idx),type="internal",use.rsq=FALSE)

## export to file
exportSI(final.index$idx,ages,years,toy=mean(d[["HH"]]$timeOfYear,na.rm=TRUE),file="index.txt",nam=paste(genus,bfamily," BTS + SNS - Delta-GAM; Last age is plus group, calculated",Sys.time()))

