### Biodiversity metric and Maps generation ###

#rm(list=ls())
wd <- "~/Chamaecrista_conservation/" # main directory
wd_polygon <- "polygons/shape_files/" # where we saved the polygons
wd_current <- "rasters/CHELSA_cur_V1_2B_r2_5m/2_5min/" # climatic rasters
wd_results <- "results_maps"

library(raster)
library(dplyr)
library(maptools)
library(sf)
library(rgdal)
library(picante)
data(wrld_simpl)

# Loading the tree
setwd(wd)
tree <- read.tree("final_tree.txt")

# Loading EDGE data frame
setwd(paste(wd))
phylotraits <- read.csv2("edge_values.csv")[1:19,1:10] 
phylotraits$Species <- sub(" ", "_", phylotraits$Species)

# Loading EDGE2 data frame
edge2_values <- read.csv("edge2_table.csv")
edge2_values$Species <- sub("_", " ", edge2_values$Species)

# Loading polygons
setwd(paste(wd, wd_polygon, sep=""))
list.files(pattern="_current.shp") -> files_cur
sapply(files_cur, readShapePoly, simplify = F) -> models
sub("_current.shp", "", files_cur) -> labs
labs <- sub(" ","_", labs)
names(models) <- labs

setwd(paste(wd, wd_current, sep=""))
mask <- raster("bio_1.tif")
mask <- crop(mask, extent(-60, -30, -25, -8.5))
mask[!is.na(mask)] <- 0 
plot(mask, axes = T, box = T, legend = T, main = "Neotropics")

## 1. Creating PD
rasters1 <- list() 
matrixPD <- matrix(0, ncol = 1, nrow = 1) # empty matrix for PD values
matrixPD <- data.frame(matrixPD)

# Creating the ocurrence matrix
for (i in 1:length(models)) {
  Sys.time() -> start_time
  r1 <- rasterize(models[[i]], mask, field=1)
  r1[is.na(r1)] <- 0
  resample(r1, mask) -> r1
  mask(r1, mask) -> rasters1[[i]]
  frame <- as.data.frame(values(rasters1[[i]]))
  tip <- sub(" ", "_", labs[i])
  colnames(frame) <- tip
  matrixPD <- cbind(matrixPD, frame)
  print(c(i, "in", length(models)))
  Sys.time() -> end_time
  print(end_time-start_time)
}
matrixPD[is.na(matrixPD)] <- 0

# Creating Species Richness map
names(rasters1) <- names(models)
rasters1 <- rasters1[which(names(rasters1) %in% tree$tip.label)]
res1 <- Reduce("+", rasters1, accumulate = TRUE)
l1<-length(res1)
speciesRichness <- res1[[l1]]
plot(speciesRichness, main = "Current SR") # Species richness map

# Creating the PD matrix
tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(matrixPD)))  
matrixPD <- matrixPD[,c(which(colnames(matrixPD) %in% tree$tip.label))]
pd_final <- pd(matrixPD, tree, include.root=TRUE)

pdMAP <- mask
values(pdMAP) <- pd_final[,1]

# Plotting
par(mfrow=c(1,2))
plot(speciesRichness, main = "SR")
plot(wrld_simpl, add = T)
plot(pdMAP, main = "PD")
plot(wrld_simpl, add = T)
dev.off()
setwd(paste(wd, wd_results, sep=""))
writeRaster(pdMAP, "PD.tif")
writeRaster(speciesRichness, "SR.tif")

# Regression and residuals for PD
linearRegression <- lm(pd_final[,1]~pd_final[,2])
abline(linearRegression, col="red") # adding trend line
summary(linearRegression)
res <- as.numeric(linearRegression$residuals)
resMap <- mask
values(resMap) <- res
plot(resMap, main="residual")
plot(wrld_simpl, add = T)
resMap[speciesRichness[]==0] <- NA
resMap[is.na(speciesRichness[])] <- NA
writeRaster(resMap, "Res_PD.tif")

## 2. Creating ED, ER and EDGE maps
setwd(paste(wd, wd_results, sep=""))
r1 <- list()
metric <- colnames(phylotraits[c(4, 6:10)])

for (r in 1:length(metric)) {
  # MAP
  for(i in 1:nrow(phylotraits)){
    if(phylotraits$Species[i] %in% names(rasters1)){
      x0 <- rasters1[[which(names(rasters1)==phylotraits$Species[i])]]
      col <- metric[r]
      values(x0) <- values(x0) * phylotraits[,col][i]
      r1[[i]] <- x0
    }
    else {
      print(paste(phylotraits$Species[i], "não existe!"))
    }
  }
  r1 <- r1[!sapply(r1, is.null)]
  final <- calc(stack(r1), sum)
  writeRaster(final, paste(col,".tif", sep = ""))
  
  # Residual map
  res2 <- lm(getValues(final) ~ getValues(speciesRichness), na.action=na.exclude)
  r <- mask
  r[] <- residuals.lm(res2)
  r[speciesRichness[]==0] <- NA
  r[is.na(speciesRichness[])] <- NA
  writeRaster(r, paste("residual_", col, ".tif", sep = ""))
  
}

# 3. Creating EDGE2 maps

r2 <- list() # cria lista vazia para ser "alimentada pelo loop"

# MAP
for(i in 1:nrow(edge2_values)){
  if(edge2_values$Species[i] %in% names(rasters1)){
    x0 <- rasters1[[which(names(rasters1)==edge2_values$Species[i])]]
    values(x0) <- values(x0) * edge2_values$EDGE[i]# mudando os valores do raster
    r1[[i]] <- x0
  }
  else {
    print(paste(edge2_values$Species[i], "não existe!"))
  }
}
r2
r2 <- r2[!sapply(r2, is.null)]

edge2 <- Reduce("+", r2, accumulate = TRUE)
l1<-length(edge2)
writeRaster(edge2[[l1]], "EDGE2.tif")

# Residual
res3 <- lm(getValues(edge2) ~ getValues(speciesRichness), na.action=na.exclude)
r <- mask
r[] <- residuals.lm(res3)
r[speciesRichness[]==0] <- NA
r[is.na(speciesRichness[])] <- NA
writeRaster(r, "residual_EDGE2.tif", overwrite=T) 

### END
