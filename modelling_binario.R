#### Species Distribution Modelling and Polygons ####

#rm(list=ls())
wd <- "~/Chamaecrista_conservation/" # main directory
wd_polygon <- "polygons/shape_files/" # where we save polygons
wd_current = "rasters/CHELSA_cur_V1_2B_r2_5m/2_5min/" # climatic rasters

library(sdm)
library(usdm)
library(raster)
library(dplyr)
library(maptools)
library(sf)
data(wrld_simpl)

# Loading cleaned points 
setwd(paste(wd))
points <- read.csv2("ser. Coriaceae clean.csv")[,c("scientificName", "decimalLongitude", "decimalLatitude")]
colnames(points) <- c("species","lon","lat")

##  1.  Loading climate layers   
setwd(paste(wd, "rasters/CHELSA_cur_V1_2B_r2_5m/2_5min", sep=""))
list.files(pattern="tif$") -> files
files[c(2,7,8,9,13,14,18,19)] -> files #selecting specific layers
sapply(files, raster, simplify = F) -> current
stack(current) -> current
current <- crop(current, extent(-60, -30, -25, -8.5)) #east of Brazil - campos rupestres region
cleaned_points <- points

##  2. SDM
geo <- cleaned_points[,c("species","lon","lat")]

# delete species with 3 points or less
geo <- geo[geo$species %in% names(which(table(geo$species) > 2)), ]
species <- names(which(table(geo$species)!=0))

# Modelling
Sys.time() -> start_time_total
for(i in 1:length(species)){
  Sys.time() -> start_time
  sp0 = species[i]
  sp1 <- geo[geo$species==sp0,]
  sp1$species <- 1
  coordinates(sp1) <- ~ lon + lat
  
  setwd(paste(wd, "distribution_models/1_current/", sep=""))
  d <- sdmData(species~., sp1, predictors = current, bg = list(n=1000)) # 1000 background points 
  m <- sdm(species~., d, methods =c("glm", "brt", "rf", "maxent"), # creating the models
           replication=c('boot'), n=5)


  Stat <- getEvaluation(m,stat=c('TSS', 'AUC'),opt=1, file=paste(sp0, "current.csv", sp="")) # evaluating the models
  write.csv(Stat, file=paste(sp0, "evaluation.csv", sp=""))
  
  
  x <- getVarImp(m, id=1, wtest="test.dep") # calculating relative importance of variaables in the models
  write.table(data.frame(), file=paste0(sp0, "_varimp.txt"))
  sink(paste0(sp0, "_varimp.txt"))
  print(x)
  sink()
  plot(x,'auc')
  
  en <- ensemble(m, current, paste(sp0, ".img", sep="") ,
                 setting=list(method='weighted',stat=c("AUC")) ) # ensemble ("fusion model")

  writeRaster(en, file=paste(sp0, "current", sep="_"), format="GTiff")
  
  setwd(paste(wd, "distribution_models/1_pdfs", sep=""))
  pdf(file=paste(sp0, "_roc", ".pdf", sep=""), width = 8, height = 5)
  roc(m, 1)
  title("glm model", adj=0, add=T)
  roc(m, 2)
  title("brt model", adj=0, add=T)
  roc(m, 3)
  title("rf model", adj=0, add=T)
  roc(m, 4)
  title("maxent model", adj=0, add=T)
  dev.off()
  
  pdf(file=paste(sp0, "_sdm", ".pdf", sep=""), width = 10, height = 5)
  plot(en, zlim=c(0,1))
  plot(sp1, pch=19, cex=0.1, col="tomato1", add=T)
  title(main=paste(sp0, "current distribution", sep=" "))
  dev.off()
  
  Sys.time() -> end_time
  print(c(species[i], "done!"))
  print(end_time-start_time)
}

Sys.time() -> end_time_total
print(end_time_total-start_time_total)

##############################################################

# Loading models just generated 
setwd(paste(wd, "distribution_models/1_current/", sep=""))
list.files(pattern="tif$") -> files
sapply(files, raster, simplify = F) -> current
sub("_current.tif", "", files) -> labs

setwd(paste(wd))

##############################################################

## 3. Transforming rasters into polygons 

# Establishing thresholds
thresholdC = 0.5

ranges_total <- data.frame()
Sys.time() -> start_time_total
for(u in 1:length(current)){
  Sys.time() -> start_time
  current[[u]]-> curr.bin
  
  curr.bin[current[[u]][] < thresholdC] <- 0
  curr.bin[current[[u]][] >= thresholdC] <- 1
  curr.bin[curr.bin[] == 0] <- NA
  
  if (is.na(which.min(curr.bin))==T) {
    print(c(labs[u], "sorry, no range...")) 
  }
  else {
    
    setwd(paste(wd, "polygons/plots", sep=""))
    sp1 <- geo[geo$species==labs[u],]
    coordinates(sp1) <- ~ lon + lat
    
    pdf(file=paste(labs[u],".pdf", sep="_"), width = 15, height=10)
    plot(curr.bin, col="green", legend=F)
    plot(wrld_simpl, add=T)
    plot(sp1, pch=19, cex=0.1, col="tomato1", add=T)
    title(main=paste(labs[u], "current distribution", sep=" "))
    
    dev.off() 
    
    # plot polygon
    
    setwd(paste(wd, "polygons/shape_files", sep=""))
    rasterToPolygons(curr.bin, dissolve=T) -> rangeC 
    
    if (is.null(rangeC)==T) {
      print(c(labs[u], "sorry, no range..."))
    }
    else {
      
      # shape
      writeSpatialShape(rangeC, paste(labs[u], "_current", sep=""))
      
      # total area in each time slice
      
      print(c(labs[u], "done!"))
      Sys.time() -> end_time
      print(end_time-start_time)
      
    }
  }
}  

### END
