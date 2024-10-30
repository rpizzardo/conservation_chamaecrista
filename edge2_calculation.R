rm(list=ls())
## Calculating GE2 and EDGE2 for Chamaecrista species
# May 2023
library(ape)

wd <- "/Users/raquel.pizzardo/Documents/EDGE manual 6"
setwd(wd)

### GE2 ### See Gumbs 2023
# generate 1,000 GE2 values distributed across all Red List categories 
# and assign a distribution of GE2 values to each RL category to capture uncertainty

GE.2.calc <- function(pext){
  require(geiger)
  treesim <- sim.bdtree(n=10000)
  iucn <- sample(1:5, size=length(treesim$tip.label), replace=TRUE)
  data <- data.frame(species=treesim$tip.label, pext=pext[iucn])
  data <- data[order(data$pext),]
  data$rank <- seq_len(nrow(data))
  rank <- c(0, with(data, tapply(rank, pext, median)))
  pext <- c(0, pext)
  rank.sq <- rank^2; rank.cub <- rank^3; rank.qu <- rank^4; rank.quu <- rank^5
  model <- lm(pext ~ rank + rank.sq + rank.cub + rank.qu)
  data$rank.sq <- data$rank^2; data$rank.cub <- data$rank^3; data$rank.qu <- data$rank^4; data$rank.quu <- data$rank^5
  data$rank.pext <- predict(model, data)
  data$rank.pext[data$rank.pext <= 0] <- 0.0001
  data$rank.pext[data$rank.pext >= 1] <- 0.9999
  pext.LC <- data.frame(RL.cat = "LC", pext =data$rank.pext[data$pext == pext[2]])
  pext.NT <- data.frame(RL.cat = "NT", pext =data$rank.pext[data$pext == pext[3]])
  pext.VU <- data.frame(RL.cat = "VU", pext =data$rank.pext[data$pext == pext[4]])
  pext.EN <- data.frame(RL.cat = "EN", pext =data$rank.pext[data$pext == pext[5]])
  pext.CR <- data.frame(RL.cat = "CR", pext =data$rank.pext[data$pext == pext[6]])
  return(rbind(pext.CR,pext.EN, pext.VU, pext.NT, pext.LC))
}

EDGE.2.calc <- function(tree, pext){
  require(phylobase)
  require(data.table)
  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }
  tree_dat <- data.frame(Species = as.character(unique(tree$tip.label)),
                         TBL = tree$edge.length[sapply(c(1:length(tree$tip.label)),
                                                       function(x,y) which(y==x),y=tree$edge[,2])], 
                         pext = NA, ED = NA, EDGE = NA)
  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)
  tree <- as(tree, "phylo4")
  names(pext) <- c("species","pext")
  for(i in 1:length(tree_dat$Species)){
    tree_dat$pext[i] <- pext$pext[pext$species == tree_dat$Species[i]]
  }
  nodes <- descendants(tree, rootNode(tree), "all")
  for(i in 1:length(nodes)){
    tips <- descendants(tree, nodes[i], "tips")
    tips <- names(tips)
    tipscores <- which(pext$species %in% tips)
    tree@edge.length[which(tree@edge[,2] == nodes[i])] <- edgeLength(tree, nodes[i])*prod(as.numeric(pext$pext[tipscores]))
  }
  for(i in 1:length(tree_dat$Species)){
    tree_dat$EDGE[i] <- sum(tree@edge.length[which(tree@edge[,2] %in% ancestors(tree, 
                                                                                which(tipLabels(tree) == tree_dat$Species[i]), "ALL"))], na.rm=T)
    tree_dat$ED[i] <- tree_dat$EDGE[i] / as.numeric(tree_dat$pext[i])
  }
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)
  edge.res <- list(tree_dat,tree,ePD.dat)
  return(edge.res)
}


# calculating ge2
pext <- rev(c(0.97, 0.97/2, 0.97/4,0.97/8,0.97/16))
ge2 <- GE.2.calc(pext)


### EDGE2 ### See Gumbs 2023
# provide phylogenetic tree and dataframe with two columns: 
# the first comprising species names, the second comprising their associated GE2 scores (between 0 and 1)
# function returns three objects: 
# 1. dataframe with terminal branch length, GE2, ED2 and EDGE2 scores for each species
# 2. expected PD loss tree
# 3. PD and expected PD loss in MY for the clade

tree_c <- read.tree("final_tree.txt")
data <- read.csv("edge2_calc.csv")
data_to_edge <- data
results <- list()
for(i in 1:1000) {
  for(species_index in 1:nrow(data)){
    one_ge_value <- data[,2][species_index]
    if(one_ge_value=="EN") {
      one_sample <- round(sample(ge2$pext[which(ge2$RL.cat=="EN")], 1),3)
    }
    if(one_ge_value=="VU") {
      one_sample <- round(sample(ge2$pext[which(ge2$RL.cat=="VU")], 1),3)
    }
    if(one_ge_value=="LC") {
      one_sample <- round(sample(ge2$pext[which(ge2$RL.cat=="LC")], 1),3)
    }
    data_to_edge$iucn[species_index]<-one_sample
  }
  #colnames(data_to_edge)[2] <- "GE2"
  one_calc <- EDGE.2.calc(tree=tree_c, pext=data_to_edge)
  results[[i]] <- one_calc[[1]]
  cat(i, "\r")
}

all_results <- do.call(rbind, results)

#Edn