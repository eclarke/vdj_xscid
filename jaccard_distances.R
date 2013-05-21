##
# Abundance-weighted Jaccard index
# 
# Erik Clarke [ecl@mail.med.upenn.edu]
##

library(data.table)
library(reshape2)
library(ape)
library(scatterplot3d)

## Change this to true to re-create the data structure
if (F) {
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
    species[is.na(species)] <- 0
    save(species, file="species_by_sample.RData")
} else {
    load("species_by_sample.RData")
}

#comparisons <- combn(colnames(species)[-1], 2)

samples.all <- colnames(species)[-1]
accn.all <- sapply(samples.all, substr, 1, 8)
samples.patients <- samples.all[5:length(samples.all)]
accn.patients <- accn.all[5:length(accn.all)]

jaccardDistance <- function(s1, s2, data, weighted) {
    a <- data$species[data[, s1] > 0 & data[, s2] == 0]
    b <- data$species[data[, s1] == 0 & data[, s2] > 0]
    c <- data$species[data[, s1] > 0 & data[, s2] > 0]

    if (weighted) { fn <- sum } else { fn <- length }

    wa <- fn(data[a, s1])
    wb <- fn(data[b, s2])
    wc <- fn(data[c, c(s1, s2)])

    distance = (wa + wb) / (wa + wb + wc)

    return(distance)
}

##
# This function simply calls the jaccardDistance function over all
# pairwise combinations of the supplied species.names vector.
# If weighted = T, it uses the copy numbers to abundance-weight the
# Jaccard distance.
##
pairwiseOverSamples <- function(sample.names, weighted, na.diag=FALSE) {
    accns <- substr(sample.names, 1, 8)
    status <- ifelse(grepl("GTSP", accns), "patient", "control")
    samples <- sample.names
    results <- sapply(sample.names, 
                      function(x) {
                        sapply(sample.names, jaccardDistance, s2=x, 
                               data=species, weighted=weighted)
                        })

    diag(results) <- ifelse(na.diag, NA, diag(results))
    res.df <- list(accn=accns, status=status, sample=samples, 
                         results=results)
    return(res.df)
}

pw.awj.all <- pairwiseOverSamples(samples.all, weighted=T)
pw.uwj.all <- pairwiseOverSamples(samples.all, weighted=F)
pw.awj.patients <- pairwiseOverSamples(samples.patients, weighted=T)
pw.uwj.patients <- pairwiseOverSamples(samples.patients, weighted=F)

pc.scatter3d <- function(res.df, color='status', ...) {
    pc <- pcoa(res.df$results)
    df <- data.frame(res.df[1:3])
    pc$df <- cbind(pc$vectors, df)
    if (color == 'status'){
        pc$df$pcolor <- ifelse(pc$df$status == "patient", "red", "blue")
    } else if (color == 'accn') {
        a <- unique(df$accn)
        ra <- rainbow(length(a))
        names(ra) <- a
        print(ra)
        pc$df$pcolor <- sapply(df$accn, function(x) {cat(x); ra[[x]]})
    }

    PC1 <- pc$df[,'Axis.1']
    PC2 <- pc$df[,'Axis.2']
    PC3 <- pc$df[,'Axis.3']
    s3d <- scatterplot3d(PC1, PC2, PC3, color=pc$df$pcolor, pch=19, type='h',
                         lty.hplot=2, scale.y=0.75, ...)
    s3d.coords <- s3d$xyz.convert(PC1, PC2, PC3)
    text(s3d.coords$x, s3d.coords$y, labels=rownames(pc$df), pos=4, cex=0.5)
}



