require(tidyverse)

################################################################################
# Construction of randomized protein sequences
################################################################################
# Human aafrequencies were taken from:
# Bios 81(1):22-31 (2010), Tsuji et al.

aacids <- c("A","C","D","E","F","G","H","I","K","L",
            "M","N","P","Q","R","S","T","V","W","Y")
aafreq <- c(6.97, 2.31, 4.70, 6.84, 3.79, 6.60, 2.59, 4.42, 5.62, 10.02,
            2.21, 3.60, 6.11, 4.64, 5.68, 8.10, 5.31, 6.08, 1.32, 2.74)
names(aafreq) <- aacids
nonTrypticAAs <- aacids[!aacids %in% c("K","R")]
nonTaafreq <- aafreq[-which(names(aafreq) %in% c("K","R"))]

readFasta <- function(fastaFile) {
    fasta <- read_file(file=fastaFile)
    fasta.names <- unlist(str_extract_all(fasta, "(?<=>).+(?=\\n)"))
    fasta.names <- str_replace_all(fasta.names, "\\n", "")
    fasta <- unlist(str_split(fasta, ">.+(?=\\n)"))[-1]
    fasta <- str_replace_all(fasta, "\\n", "")
    names(fasta) <- fasta.names
    return(fasta)
}

makeDecoySeq <- function(proteinSeq, multiplier=1) {
    proteinLength <- nchar(proteinSeq)
    Mstart <- str_sub(proteinSeq, 1, 1) == "M"
    KRposition <- unlist(str_locate_all(proteinSeq, "[KR](?!P)"))
    KRposition <- map(1:multiplier, function(x) KRposition + proteinLength * (x - 1))
    KRposition <- do.call(c, KRposition)
    proteinLength <- proteinLength * multiplier
    if (Mstart) {
        decoySeq <- 
            paste0(c("M", sample(nonTrypticAAs, 
                                 proteinLength - 1, 
                                 replace=T, 
                                 prob=nonTaafreq)), collapse="")
    } else {
        decoySeq <- 
            paste0(sample(nonTrypticAAs, 
                          proteinLength, 
                          replace=T, 
                          prob=nonTaafreq), collapse="")
    }
    for (KR in KRposition) {
        str_sub(decoySeq, KR, KR) <- 
            sample(c("K","R"), 1, prob=aafreq[c("K","R")])
    }
    return(decoySeq)
}

makeDecoyDB <- function(fwdDB, multiplier=1) {
    #fwdDB is named character vector.  Each element is protein sequence.
    decDB <- map_chr(fwdDB, makeDecoySeq, multiplier)
    names(decDB) <- str_c(">", names(fwdDB))
    return(decDB)
}

wrap60 <- function(stringIn) {
    numberLines <- str_length(stringIn) %/% 60 + 1
    for (i in seq(numberLines)) {
        str_sub(stringIn, i * 60 + i, i * 60 + (i - 1)) <- "\n"
    }
    return(stringIn)
}

makeDecoyFile <- function(fwdDBfile, multiplier=1) {
    fwdDB <- readFasta(fwdDBfile)
    nameReg <- paste0("_dec", multiplier, "\\1")
    decDBfile <- str_replace(fwdDBfile, "(\\.[^\\.]+)$", nameReg)
    decDB <- makeDecoyDB(fwdDB, multiplier)
    decDB <- map(decDB, wrap60)
    flatDB <- paste0(map2_chr(names(decDB), decDB, paste, sep="\n"), collapse="\n")
    write_file(flatDB, decDBfile)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    fwdFile <- args[1]
    multiplier <- as.integer(args[2])
    makeDecoyFile(fwdFile, multiplier)
}

main()



