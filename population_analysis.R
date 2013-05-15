load('all_data.RData')
library(data.table)

# isolate just the processed data from Adaptive
ad.f <- subset(all.data, !is.na(filt.copy))
ad.f$key <- paste(ad.f$accn, ad.f$seq, sep="_")
# collapse the replicate copy number data 
cn <- data.table(ad.f)
cn <- cn[, list(filt.copy=mean(filt.copy), normalizedCopy=mean(normalizedCopy),
         normalizedFrequency=mean(normalizedFrequency), 
         rawFrequency=mean(rawFrequency), accn=accn), by=key]

accns <- unique(ad.f$accn)


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
    cat(accn1, accn2, repl1, repl2)
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
    } else {
        stopifnot(accn != NULL)
        repls <- unique(ad.f$replicate[ad.f$accn == accn])
        results <- as.matrix(cbind(sapply(repls, .by.repl, repls=repls,
                             df=df, accn=accn, column=column, test=test, 
                             full, t.value, ...)))
    }
    return(results)

}