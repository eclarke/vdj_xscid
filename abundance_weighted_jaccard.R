##
# Abundance-weighted Jaccard
# 
# Erik Clarke [ecl@mail.med.upenn.edu]
##

library(data.table)
library(reshape2)

## Change this to true to re-create the data structure
if (TRUE) {
    load('all_processed.Rdata')
    # the key is based on the V and J genes, as well as the sample (accn+repl)
    ap$key <- paste(ap$accn, ap$replicate, ap$VGeneName, ap$JGeneName)
    genes <- data.table(all.processed)
    setkey(genes, key)
    # the genes data table aggregates the raw and normalized copy numbers by
    # species and sample
    genes <- genes[, list(raw.copy = sum(filt.copy), 
                          normalized.copy=sum(normalizedCopy), 
                          VGeneName=VGeneName, DGeneName=DGeneName, 
                          JGeneName=JGeneName, accn=accn, 
                          replicate=replicate), by=key]
    # we recast the genes data table in order to get the total number of
    # species per sample (using the raw.copy data because I don't know/trust
    # whatever Adaptive is doing with the normalized data and may not be comp.
    # across samples)
    species.by.sample <- dcast(genes, species~sample, value.var='raw.copy',
                               fun.aggregate=mean)
    save(species.by.sample, file="species_by_sample.RData")
    
} else {
    load("species_by_sample.RData")
}