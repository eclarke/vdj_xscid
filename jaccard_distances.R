##
# Abundance-weighted Jaccard index
# 
# Erik Clarke [ecl@mail.med.upenn.edu]
##

library(data.table)
library(reshape2)
library(ape)
library(scatterplot3d)
library(pheatmap)
library(vegan)

CreateSpeciesDataframe <- function() {
    load('all_processed.Rdata')
  # the key is based on the V and J genes, as well as the sample (accn+repl)
  ap <- subset(ap, !is.na(VGeneName) & !is.na(JGeneName))
  ap$key <- paste(ap$accn, ap$replicate, ap$VGeneName, ap$JGeneName)
  genes <- data.table(ap)
  setkey(genes, key)
  # the genes data table aggregates the raw and normalized copy numbers by
  # species and sample
  genes <- genes[, list(raw.copy = sum(filt.copy), 
                        normalized.copy=sum(normalizedCopy), 
                        VGeneName=VGeneName, DGeneName=DGeneName, 
                        JGeneName=JGeneName, accn=accn, 
                        replicate=replicate), by=key]
  genes$species <- paste(genes$VGeneName, genes$JGeneName, sep='_')
  genes$sample <- paste(genes$accn, genes$replicate, sep='_')
  # we recast the genes data table in order to get the total number of
  # species per sample (using the raw.copy data because I don't know/trust
  # whatever Adaptive is doing with the normalized data; it may not be 
  # comparable across samples)
  species <- dcast(genes, species~sample, value.var='raw.copy',
                   fun.aggregate=mean)
  rownames(species) <- species$species
  species$species <- NULL
  species[is.na(species)] <- 0
  save(species, file="species_by_sample.RData")
  return(species)
}

JaccardDistance <- function(s1, s2, data, weighted) {
  # Computes the (optionally abundance-weighted) Jaccard distance between two 
  # samples.
  #
  # Args:
  #   s1: The first sample. Needs to be the name of a column in the dataframe.
  #   s2: The second sample. Also a column name.
  #   data: The data to use. The columns should have abundance data if using
  #         an abundance-weighted Jaccard distance.
  #   weighted: If TRUE, weigh the Jaccard distance by the abundance of each
  #             species.
  #
  # Returns:
  #   The Jaccard distance between sample 1 and sample 2.
  a <- rownames(data[data[, s1] > 0 & data[, s2] == 0, ])
  b <- rownames(data[data[, s1] == 0 & data[, s2] > 0, ])
  c <- rownames(data[data[, s1] > 0 & data[, s2] > 0, ])

  if (weighted) { fn <- sum } else { fn <- length }

  wa <- fn(data[a, s1])
  wb <- fn(data[b, s2])
  wc <- fn(data[c, c(s1, s2)])

  distance = (wa + wb) / (wa + wb + wc)

  return(distance)
}


PairwiseJaccardDistance <- function(sample.names, data, weighted, 
                                    na.diag=FALSE) {
  # Calculates the pairwise Jaccard distance between a list of samples.
  # 
  # Args:
  #   sample.names: A list of samples to compare against each other
  #   data: The dataframe to use. Same requirements as JaccardDistance fn.
  #   weighted: If TRUE, weigh the Jaccard distance by abundance
  #   na.diag: If TRUE, set the diagonal of the returned matrix to NA. This is
  #            useful for preventing the identity diagonal from skewing 
  #            heatmap color scales.
  #
  # Returns:
  #   An data frame containing sample information and an n by n distance
  #   matrix, where n is the length of the sample vector.
  accns <- substr(sample.names, 1, 8)
  status <- ifelse(grepl("GTSP", accns), "patient", "control")
  samples <- sample.names
  distance <- sapply(sample.names, 
                     function(x) {
                      sapply(sample.names, JaccardDistance, s2=x, 
                             data=data, weighted=weighted)
                      })

  diag(distance) <- ifelse(na.diag, NA, diag(distance))
  df <- list(accn=accns, status=status, sample=samples, 
             distance=distance)
  return(df)
}

DistancePCoA  <- function(df) {
  # Performs principal coordinate analysis on the distance matrix calculated
  # by the PairwiseJaccardDistance function.
  # 
  # Args:
  #   df: The data frame returned by PairwiseJaccardDistance
  #
  # Returns:
  #   The PCoA object with a dataframe at $df that includes both the axes and
  #   the sample information per row.
  
  pc <- pcoa(df$distance)
  df <- data.frame(df[1:3])
  pc$df <- cbind(pc$vectors, df)
  return(pc)
}

PCPlot3D <- function(pc, color.by='status', ...) {
  # Plots the first three principal coordinates of a distance matrix
  # calculated by the DistancePCoA function.
  #
  # Args:
  #   df: The data frame returned by DistancePCoA
  #   color.by: If 'status', colors according to whether a sample is control
  #             or patient. If 'accn', colors according to GTSP/CTRL accession
  #             number.
  #
  # Returns:
  #   Plots the resulting 3D scatterplot.
  if (color.by == 'status') {
    pc$df$pcolor <- ifelse(pc$df$status == "patient", "red", "blue")
  } else if (color.by == 'accn') {
    ## assigns a unique color to each accession number
    a <- unique(pc$df$accn)
    ra <- rainbow(length(a))
    names(ra) <- a
    pc$df$pcolor <- sapply(pc$df$accn, function(x) ra[[x]])
  }

  PC1 <- pc$df[, 'Axis.1']
  PC2 <- pc$df[, 'Axis.2']
  PC3 <- pc$df[, 'Axis.3']

  plot3d(PC1, PC2, PC3, col=pc$df$pcolor, type='h', ...)
  planes3d(c(0,0,1), col='grey', alpha=0.07)
  text3d(PC1, PC2, PC3, rownames(pc$df), pos=4, cex=0.5, col=pc$df$pcolor)
}

DistanceHeatmap <- function(df, clustering=F, na.diag=T, ...) {

  if (na.diag) diag(df$distance) <- NA
  if (clustering) {
    pheatmap(df$distance, ...)
  } else {
    pheatmap(df$distance, cluster_rows=F, cluster_cols=F, ...)
  }
}
