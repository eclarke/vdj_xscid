##----------------------------------------------------------------------------
##  Data processing functions for unfiltered VDJ sequencing data
##  Erik Clarke (ecl@mail.med.upenn.edu)
##----------------------------------------------------------------------------
source("misc_functions.R")
path  <- 	"data/patients/processed/unfiltered"
files <-  list.files(path = path, pattern = glob2rx("*.adap.txt"))


## Function definitions
##----------------------------------------------------------------------------
createMetadata <- function(string) {
  r <- substr(string, 9, 9)
  # currently cannot handle replicates past "c". if more, make better soltn
  stopifnot(r %in% c("a", "b", "c"))
  replicate <- ifelse(r == "a", 1,
                      ifelse(r == "b", 2, 3))
  patient <- substr(string, 1, 8)
  path <- file.path(path, string)
  accn <- paste(substr(string, 1, 8), replicate, sep='.')
  return(data.frame(path=path, patient=patient, replicate=replicate, accn=accn))
}

readSeqData <- function(file) {
  dat.tmp <- read.delim(as.character(file), sep="", header=F, 
                        stringsAsFactors=F)
  meta    <- as.list(metadata[metadata$path == file, ])
  return(data.frame(seq=dat.tmp[,1], copy=dat.tmp[,2], meta))  
}


## Script
##----------------------------------------------------------------------------
## Data importing
cat("Collecting paths to data... ")
metadata  <- do.call(rbind, lapply(files, createMetadata))
cat("done.\n")
cat("Reading sequence data... ")
df <- do.call(rbind, lapply(metadata$path, readSeqData))
df$seq <- as.character(df$seq)
cat("done.\n")

## Identify singletons appearing across replicates
singletons <- df[df$copy == 1, c("seq", "patient")]
across.replicates <- which(duplicated(singletons))
repl.singletons <- singletons[across.replicates, ]

# Output results
n.singletons.replicated <- length(across.replicates)
n.singletons <- length(singletons$seq)
proportion.singletons <- n.singletons.replicated/n.singletons
cat("\n\t", "Total singletons across replicates:", n.singletons.replicated,
    "\n\t", "Proportion of total singletons:", n.singletons.replicated, "/", 
    n.singletons, "=", proportion.singletons, "\n\n")
print(summary(repl.singletons))


