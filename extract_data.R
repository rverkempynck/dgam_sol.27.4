
##install.packages("icesDatras")
##remotes::install_github("DTUAqua/DATRAS/DATRAS")
##remotes::install_github("casperwberg/surveyIndex/surveyIndex")

rm(list = ls())

library(DATRAS)
library(mapdata)
library(sp)

cmSize <- 1
years <- 2010:2018
genus <- "Solea"
bfamily <- "solea"
datafile <- "sole.RData"

if(!file.exists(datafile)){
    BTS <- getDatrasExchange("BTS",years=years,quarters=1:4,strict=FALSE)
    SNS <- getDatrasExchange("SNS",years=years,quarters=1:4,strict=FALSE)
    dAll <- c(BTS,SNS)
    dAll <-addSpatialData(dAll,"ICES_areas/ICES_Areas_20160601_cut_dense_3857.shp")
    d <- subset(dAll,Species==paste(genus,bfamily),Year %in% years,HaulVal=="V",StdSpecRecCode==1)
    dAll <- BTS <- SNS <- NULL
    save(d,file=datafile)
} else {
    load(datafile)
}

d <- subset(d, Country == "BEL" | (Gear == "BT8" & is.na(Rigging)) |
                Survey == "SNS")
xtabs( ~ Gear + Country, data = d[[2]])
xtabs( ~ Gear + Survey, data = d[[2]])

d <- addSpectrum(d, by = cmSize)
d <- addWeightByHaul(d)

bubblePlot(d,scale=1/10)

save(d, file = "final_dataset.RData")
