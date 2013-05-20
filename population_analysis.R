library(parallel)
if (F) {
  load('all_processed_data.RData')
  library(data.table)
  
  # isolate just the processed data from Adaptive
  ad.f <- subset(df.filt.all, !is.na(filt.copy))
  ad.f$key <- paste(ad.f$accn, ad.f$seq, sep="_")
  # collapse the replicate copy number data 
  cn <- data.table(ad.f)
  cn <- cn[, list(filt.copy=mean(filt.copy), normalizedCopy=mean(normalizedCopy),
                  normalizedFrequency=mean(normalizedFrequency), 
                  rawFrequency=mean(rawFrequency), accn=accn), by=key]
  # or, collapse the data according to V and J genes (D is too inconsistent)
  ad.f$key <- paste(ad.f$accn, ad.f$replicate, ad.f$VGeneName, ad.f$JGeneName)
  genes <- data.table(ad.f)
  setkey(genes, key)
  genes <- genes[, list(filt.copy = sum(filt.copy), 
                        normalized.copy=sum(normalizedCopy), 
                        VGeneName=VGeneName, DGeneName=DGeneName, 
                        JGeneName=JGeneName, accn=accn, 
                        replicate=replicate), by=key]

  accns <- unique(ad.f$accn)
} else if (T) {
  load("data_collapsed_by_VJgenes.RData")
}

##
# Performs the specified test between one accn/replicate and another
# using data from the specified data frame and column. 
# By default, the function returns just the resultant p-value; this can
# be modified by specifying full=T (although then the resulting matrix from
# `test.all` may be a bit difficult to parse).
# `...` specifies arguments to pass to the given test.
##
.test <- function(accn1, accn2, repl1, repl2, df, column = "normalizedCopy",
                  test = ks.test, full=F, t.value = "p.value", ...) {
    print(paste(accn1, accn2, repl1, repl2, sep=" "))
    x <- subset(df, accn == accn1 & replicate == repl1 & !is.na(column)
                )[, column]
    y <- subset(df, accn == accn2 & replicate == repl2 & !is.na(column)
                )[, column]
    t <- test(x, y, ...)
    if (full) return(t) else return(t[[t.value]])
}

##
# Calls the `.test` function on a pairwise comparison between the specified
# accession and every other. The `repl` argument specifies the replicate
# used to select the values representative of those accessions.
##
.by.accn <- function(accn, df, repl, column, test, ...) {
    sapply(accns, .test, df=df, accn2=accn, repl1=repl, repl2=repl, 
           column=column, test=test, ...)
}

##
# Calls the `.test` function on a pairwise comparison between the specified
# sample replicate and every other replicate of that sample.
##
.by.repl <- function(repl, df, repls, accn, column, test, ...) {
    sapply(repls, function(x) .test(accn1 = accn, accn2 = accn, repl1=x, 
                                    repl2 = repl, df=df, column=column, 
                                    test=test, ...))
}

.by.all <- function(x, all.samples, df, column, test, cl, t.value, full, ...) {
  cat("X", x, '\n')
  print(all.samples)
  parApply(cl = cl, X=all.samples, MARGIN=1, FUN=function(y) {
    cat("Y", y, '\n')
    .test(accn1 = x[[1]], accn2 = y[[1]], repl1 = x[[2]], repl2 = y[[2]],
          df = df, column = column, test = test, t.value = t.value, full, ...) })
}

##
# Performs a pairwise comparison using the specified test against either all
# samples or all replicates in a sample, using the specified data frame and 
# column containing the data.
# 
# test:     the statistical test to use
# df:       the data frame containing the data (usually ad.f, 
#               from when the script was loaded)
# by:       if "accn", compare samples using replicate 1. if anything else,
#               compare replicates from a given accn
# column:   the column of the given data frame that contains the data points
# accn:     the accn (sample) whose replicates are being compared.
#               required if `by` != "accn".
# full:     return the full results of the test
# t.value:  if not returning the full results, what named attribute of the 
#               test results to return
# 
##
test.all <- function(test, df, by = "accn", column = "normalizedCopy", 
                     accn=NULL, full=F, t.value="p.value", ...) {
    if (by == "accn") {
        results <- as.matrix(cbind(sapply(accns, .by.accn, df=df, repl="1", 
                             column=column, test=test, full, t.value, ...)))
    } else if (by == "repl") {
        stopifnot(accn != NULL)
        repls <- unique(ad.f$replicate[ad.f$accn == accn])
        results <- as.matrix(cbind(sapply(repls, .by.repl, repls=repls,
                             df=df, accn=accn, column=column, test=test, 
                             full, t.value, ...)))
    } else if (by == "all") {
      # not implemented, use parallel version
    }
    return(results)

}

##
# Exactly the same as test.all, but accepts a cluster to run the tests in
# parallel
##
par.test.all <- function(test, df, by = "all", column = "filt.copy", 
                     accn=NULL, full=F, t.value="p.value", ...) {
  cat("making cluster...")
  cl <- makeCluster(16)
  clusterEvalQ(cl, source("population_analysis.R"))
  if (by == "accn") {
    results <- as.matrix(cbind(parSapply(cl, accns, .by.accn, df=df, repl="1", 
                                      column=column, test=test, full, t.value, ...)))
  } else if (by == "repl") {
    stopifnot(accn != NULL)
    repls <- unique(ad.f$replicate[ad.f$accn == accn])
    results <- as.matrix(cbind(parSapply(cl, repls, .by.repl, repls=repls,
                                      df=df, accn=accn, column=column, test=test, 
                                      full, t.value, ...)))
  } else if (by == "all") {
    cat("running all by all..\n")
    all.s <- unique.data.frame(subset(ad.f, select=c('accn', 'replicate')))
    all.s <- all.s[order(all.s$accn, all.s$replicate),]
    rownames(all.s) <- paste(all.s$accn, all.s$replicate, sep="_")
    results <- apply(X=all.s, FUN=.by.all, MARGIN=1, cl=cl, all.samples=all.s, 
                     df=df, column=column, test=test, t.value=t.value,
                     full=full, ...)
  }
  stopCluster(cl)
  return(results)
}

