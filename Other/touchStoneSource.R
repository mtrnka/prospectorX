setClass(Class="PPsearchCompareXL",
         slots=c(
             parentDir="character",
             dataTable="data.frame",
             peaklistDir="character",
             fastaDB="character",
             modulFile="list",
             chainMap="list",
             caTracedPDB="data.frame",
             reductionState="factor",
             expSource="character")
)

setMethod("initialize", "PPsearchCompareXL",
          function(.Object, 
                   directory, 
                   dataFile, 
                   modFile, 
                   pdbDir=NA_character_, 
                   pdbFile=NA_character_, 
                   chainMapFile=NA_character_, 
                   redState=factor("spectralMatchPairs"), 
                   expSource=NA_character_){
              .Object@expSource <- expSource
              levels(.Object@reductionState) <- c("spectralMatchPairs",
                                                  "xlinkedResPairs",
                                                  "xlinkedPepPairs")
              .Object@reductionState = redState
              cat("*** Reading Prospector XL Search Table *** \n")
              .Object@parentDir <- directory
              setwd(directory)
              .Object@dataTable <- readProspectorXLOutput(dataFile)
              cat("*** Reading Modules *** \n")
              .Object@modulFile <- readModuleFile(modFile)
              .Object@dataTable <- populateModules(.Object@dataTable,
                                                   .Object@modulFile)
              cat("*** Calculating Decoys *** \n")
              .Object@dataTable <- calculateDecoys(.Object@dataTable)
              cat("*** Calculating Xlinked Pairs *** \n")
              .Object@dataTable <- calculatePairs(.Object@dataTable)
              .Object@dataTable <- assignXLinkClass(.Object@dataTable)
              .Object@dataTable <- calculatePercentMatched(.Object@dataTable)
              .Object@dataTable <- calculatePeptideLengths(.Object@dataTable)
              .Object@dataTable <- lengthFilter(.Object@dataTable, minLen = 4, maxLen = 25)
              .Object@dataTable <- scoreFilter(.Object@dataTable, minScore = 0)
              if (!is.na(pdbFile)) {
                  cat("*** Calculating Distances *** \n")
                  chainFile <- paste(pdbDir,"/",chainMapFile,sep="")
                  pdbFile <- paste(pdbDir,"/",pdbFile,sep="")
                  .Object@chainMap <- readChainMap(chainFile)
                  .Object@caTracedPDB <- parsePDB(pdbFile)
                  .Object@dataTable <- measureDistances(.Object@dataTable,
                                                        .Object@caTracedPDB,
                                                        .Object@chainMap)
              }
              #Add validity check on .Object
              #validObject(.Object)
              return(.Object)
          }
)

setMethod("show", "PPsearchCompareXL",
          function(object){
              cat("*** Class: PPsearchCompareXL *** \n")
              cat("*** Truncated Data table *** \n")
              print(head(object@dataTable, 5))
          }
)

setGeneric("getSearchTable", function(object) {
    standardGeneric("getSearchTable")
}
)

setMethod("getSearchTable", "PPsearchCompareXL",
          function(object){
              return(object@dataTable)
          }
)

setGeneric("bestResPair", function(object, ...) {
    standardGeneric("bestResPair")
}
)

setMethod("bestResPair", "PPsearchCompareXL",
          function(object, classifier=c("Score.Diff","Score",
                                        "Peptide.2","Peptide.1",
                                        "MSMS.Info")) {
              datTab <- object@dataTable
              datTabR <- object@dataTable
              xlinks <- unique(datTab$xlinkedResPair)
              for (param in classifier) {
                  datTabR <- reduceTo(datTabR, "xlinkedResPair", param)
              }
              cat(c(nrow(datTabR), length(xlinks),
                    sum(datTabR$num), nrow(datTab), "\n"), sep="\t")
              object@reductionState = factor("xlinkedResPairs")
              object@dataTable <- datTabR
              return(object)
          }
)

setGeneric("bestPepPair", function(object, ...) {
    standardGeneric("bestPepPair")
}
)

setMethod("bestPepPair", "PPsearchCompareXL",
          function(object, classifier=c("Score.Diff","Score",
                                        "Peptide.2","Peptide.1",
                                        "MSMS.Info")) {
              datTab <- object@dataTable
              datTabR <- object@dataTable
              xlinks <- unique(datTab$xlinkedPepPair)
              for (param in classifier) {
                  datTabR <- reduceTo(datTabR, "xlinkedPepPair", param)
              }
              cat(c(nrow(datTabR), length(xlinks),
                    sum(datTabR$num), nrow(datTab), "\n"), sep="\t")
              return(datTabR)
          }
)

reduceTo <- function(datTab, category, classifier) {
    if (class(datTab[[classifier]]) == "factor") {
        datTab[[classifier]] = as.character(datTab[[classifier]])
    }
    bestScores <- aggregate(datTab[[classifier]] ~ datTab[[category]], data=datTab, max)
    names(bestScores) <- c(category, classifier)
    outTab <- merge(bestScores, datTab, by=c(category, classifier))
    if (!("num" %in% colnames(outTab))) {
        xlinks <- table(datTab[[category]])
        outTab$num <- sapply(outTab[[category]], function(x) xlinks[x])
    }
    if ("Area.Pk.1" %in% colnames(outTab) & "Area.Pk.2" %in% colnames(outTab)) {
        xlinks <- aggregate(cbind(datTab$Area.Pk.1, datTab$Area.Pk.2) ~ 
                                datTab[[category]], data=datTab, sum)
        row.names(xlinks) <- xlinks[,1]
        outTab$Area.Pk.1 <- sapply(outTab[[category]], function(x) xlinks[x, 2])
        outTab$Area.Pk.2 <- sapply(outTab[[category]], function(x) xlinks[x, 3])
    }
    if ("Int.Pk.1" %in% colnames(outTab) & "Int.Pk.2" %in% colnames(outTab)) {
        xlinks <- aggregate(cbind(datTab$Int.Pk.1, datTab$Int.Pk.2) ~ 
                                datTab[[category]], data=datTab, sum)
        row.names(xlinks) <- xlinks[,1]
        outTab$Int.Pk.1 <- sapply(outTab[[category]], function(x) xlinks[x, 2])
        outTab$Int.Pk.2 <- sapply(outTab[[category]], function(x) xlinks[x, 3])
    }
    return(outTab)
}

readProspectorXLOutput <- function(inputFile){
    inFile <- file(inputFile,open="r")
    header <- readLines(inFile,n=1)
    dataTable <- readLines(inFile,n=-1)
    close(inFile)
    #Clean up header:
    header <- lineSplit(header)
    header <- gsub("[[:space:]]",".",header)
    header <- gsub("#","Num",header)
    acc_pos <- grep("Acc",header)
    header[acc_pos[1]] <- "Acc.1"
    header[acc_pos[2]] <- "Acc.2"
    spec_pos <- grep("Species",header)
    header[spec_pos[1]] <- "Species.1"
    header[spec_pos[2]] <- "Species.2"
    prot_pos <- grep("Protein",header)
    header[prot_pos[1]] <- "Protein.1"
    header[prot_pos[2]] <- "Protein.2"
    #Clean up data:
    dataTable <- rbind(sapply(as.list(dataTable),lineSplit))
    dataTable <- as.data.frame(t(dataTable),stringsAsFactors=F)
    names(dataTable) <- header
    #Add stuff to remove useless columns from data table.
    return(cleanTypes(dataTable))
}

readModuleFile <- function(modFile) {
    inFile <- file(modFile, open="r")
    header <- readLines(inFile, n=1)
    header <- lineSplit(header)
    dataTable <- character(5)
    while (T) {
        line <- readLines(inFile, n=1)
        if (length(line)==0) {break}
        line <- lineSplit(line)
        if (length(line)==0) {break}
        else if (length(line) < 5) {length(line)=5}
        dataTable <- rbind(dataTable,line)}
    close(inFile)
    row.names(dataTable) <- NULL
    dataTable <- as.data.frame(dataTable[-1,],stringsAsFactors=F)
    names(dataTable) <- header
    dataTable <- cleanTypes(dataTable)
    modules <- list()
    for (subunit in unique(dataTable$Subunit)) {
        conv <- dataTable[dataTable$Subunit==subunit,
                          c("Range_low","Range_high","Module")]
        modules[[subunit]] <- conv[order(conv$Range_low),]
    }
    return(modules)
}

readChainMap <- function(chainFile) {
    chainMap <- read.table(chainFile,header=F,sep="\t",stringsAsFactors=F)
    names(chainMap) <- c("Subunit","Chain")
    chains <- list()
    for (su in unique(chainMap$Subunit)) {
        chains[[su]]=unlist(chainMap[chainMap$Subunit==su,"Chain"])
    }
    return(chains)
}

parsePDB <- function(pdbAcc) {
    require(bio3d)
    if (grepl("\\.cif$",pdbAcc)) {
        struct <- read.cif(pdbAcc)
    } else if (grepl("\\.pdb$",pdbAcc)) {
        struct <- read.pdb(pdbAcc)
    } else {
        stop("*** structure file is neither pdb nor cif ***")
    }
    ca.ind <- atom.select(struct,"calpha")
    struct <- struct$atom[ca.ind$atom,]
    return(struct)
}

euclideanDistance <- function(parsedPDB, residue1, chain1, residue2, chain2) {
    if (chain1 == "Dec" | chain2 == "Dec") {
        return(sample(0:250, 1, replace=F))
    }
    coord1 <- parsedPDB[parsedPDB$chain==chain1 & parsedPDB$resno==residue1,
                        c("x","y","z")]
    if (nrow(coord1)==0) coord1 <- NA
    coord2 <- parsedPDB[parsedPDB$chain==chain2 & parsedPDB$resno==residue2,
                        c("x","y","z")]
    if (nrow(coord2)==0) coord2 <- NA
    distance <- sqrt(sum((coord1 - coord2)**2))
    return(round(distance,2))
}

multiEuclideanDistance <- function(parsedPDB, residue1, chain1, residue2, chain2) {
    minDist <- NA
    minChainPair <- NA
    for (chains1 in chain1) { for (chains2 in chain2) {
        dist <- euclideanDistance(parsedPDB, residue1, chains1, residue2, chains2)
        if (is.na(minDist) | dist < minDist) {
            minDist <- dist
            if (chains1 <= chains2) {
                minChainPair <- paste(chains1, chains2, sep="")
            } else {
                minChainPair <- paste(chains2, chains1, sep="")
            }
        }
    }}
    return(minDist)
    #    return(list(minDist, minChainPair))
}

measureDistances <- function(searchTable, parsedPDB, chainMap) {
    chainLookup <- function(proteinName) {
        if (grepl("r[0-9]\\_",proteinName)) {
            return ("Dec")
        } else if (is.null(chainMap[[proteinName]])) {
            return("")
        } else {chain <- chainMap[[proteinName]]
        return(chain)
        }
    }
    pdbList <- list(parsedPDB)
    chains1 <- vapply(searchTable$Protein.1, chainLookup, FUN.VALUE="")
    chains2 <- vapply(searchTable$Protein.2, chainLookup, FUN.VALUE="")
    distance <- mapply(
        multiEuclideanDistance,
        pdbList,
        searchTable$XLink.AA.1,
        chains1,
        searchTable$XLink.AA.2,
        chains2)
    searchTable$distance <- distance
    return(searchTable)
}

populateModules <- function(searchTable, moduleDefinitions) {
    modules <- list(moduleDefinitions)
    pep1.modul <- mapply(
        assignModule,
        searchTable$XLink.AA.1,
        searchTable$Acc.1,
        modules)
    pep2.modul <- mapply(
        assignModule,
        searchTable$XLink.AA.2,
        searchTable$Acc.2,
        modules)
    levs <- unlist(sapply(moduleDefinitions, function(x) c(x[3])))
    names(levs) <- NULL
    levs <- unique(levs)
    searchTable$Modul.1 <- factor(unlist(pep1.modul),levels=levs)
    searchTable$Modul.2 <- factor(unlist(pep2.modul),levels=levs)
    return(searchTable)
}

assignModule <- function(seqPosition,protein,modList) {
    if (protein %in% names(modList)) {
        n <- findInterval(as.integer(seqPosition),
                          modList[[as.character(protein)]]$Range_low)
        if (n != 0) {
            return(modList[[protein]]$Module[n])
        } else return(NA)
    }
    else return(NA)
}

calculateDecoys <- function(searchTable) {
    searchTable$Decoy <- "Target"
    decoyReg <- "(^r[[:digit:]]_|^dec|^DECOY)"
    searchTable[grepl(decoyReg,searchTable$Acc.1) |
                    grepl(decoyReg,searchTable$Acc.2),
                "Decoy"] <- "Decoy"
    searchTable[grepl(decoyReg,searchTable$Acc.1) &
                    grepl(decoyReg,searchTable$Acc.2),
                "Decoy"] <- "DoubleDecoy"
    searchTable$Decoy2 <- "Target"
    searchTable[grepl("Decoy",searchTable$Decoy), "Decoy2"] <- "Decoy"
    searchTable$Decoy2 <- factor(searchTable$Decoy2, levels=c("Decoy","Target"))
    searchTable$Decoy <- factor(searchTable$Decoy, levels=c("DoubleDecoy","Decoy","Target"))
    return(searchTable)
}

calculatePairs <- function(searchTable){
    searchTable$Res.1 <- paste(searchTable$XLink.AA.1, searchTable$Acc.1, sep=".")
    searchTable$Res.2 <- paste(searchTable$XLink.AA.2, searchTable$Acc.2, sep=".")
    searchTable$xlinkedResPair <- ifelse(searchTable$Res.1 <= searchTable$Res.2, 
                                         paste(searchTable$Res.1, searchTable$Res.2, sep="::"),
                                         paste(searchTable$Res.2, searchTable$Res.1, sep="::"))
    searchTable$xlinkedResPair <- as.factor(searchTable$xlinkedResPair)
    searchTable$xlinkedPepPair <- ifelse(searchTable$DB.Peptide.1 <= searchTable$DB.Peptide.2,
                                         paste(searchTable$DB.Peptide.1, searchTable$DB.Peptide.2, sep="::"),
                                         paste(searchTable$DB.Peptide.2, searchTable$DB.Peptide.1, sep="::"))
    searchTable$xlinkedPepPair <- as.factor(searchTable$xlinkedPepPair)
    return(searchTable)
}

assignXLinkClass <- function(searchTable){
    searchTable$xlinkClass <- "intraProtein"
    searchTable[searchTable$Protein.1 != searchTable$Protein.2, "xlinkClass"] <- "interProtein"
    return(searchTable)
}

calculatePercentMatched <- function(searchTable) {
    if ("Num.Pks" %in% names(searchTable) & "Num.Unmat" %in% names(searchTable)) {
    num.Matched <- searchTable$Num.Pks - searchTable$Num.Unmat
    searchTable$percMatch <- num.Matched / searchTable$Num.Pks
    }
    return(searchTable)
}

calculatePeptideLengths <- function(searchTable) {
    searchTable$Len.Pep.1 <- nchar(searchTable$DB.Peptide.1)
    searchTable$Len.Pep.2 <- nchar(searchTable$DB.Peptide.2)
    return(searchTable)
}

lengthFilter <- function(searchTable, minLen, maxLen) {
    searchTable <- searchTable[searchTable$Len.Pep.2 >= minLen & searchTable$Len.Pep.1 >= minLen,]
    searchTable <- searchTable[searchTable$Len.Pep.2 <= maxLen & searchTable$Len.Pep.1 <= maxLen,]
    return(searchTable)
}

scoreFilter <- function(searchTable, minScore=0) {
    searchTable <- searchTable[searchTable$Sc.1 >= minScore & searchTable$Sc.2 >= minScore,]
    return(searchTable)
}


lineSplit <- function(line) {
    l <- nchar(line)
    g <- gregexpr("\t",line)[[1]]
    gl <- g[length(g)]
    if (l == gl) {
        return(c(unlist(strsplit(line,"\t")),""))
    } else {
        return(unlist(strsplit(line,"\t")))
    }
}

cleanTypes <- function(dataTable) {
    return(as.data.frame(
        lapply(dataTable, function(x) {
            type.convert(x, as.is=TRUE)}),
        stringsAsFactors=FALSE))
}

calculateFDR <- function(datTab, threshold=0, scalingFactor=10) {
    datTab <- datTab[datTab$dvals >= threshold, ]
    fdrTable <- table(datTab$Decoy)
    if (is.na(fdrTable["DoubleDecoy"])) {fdrTable["DoubleDecoy"] <- 0}
    if (is.na(fdrTable["Decoy"])) {fdrTable["Decoy"] <- 0}
    ffTT <- fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2)
    ftTT <- (fdrTable[["Decoy"]] / scalingFactor) - 
        (2 * fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2))
    fdr <- (ffTT + ftTT) / fdrTable[["Target"]]
    #fdr <- (fdrTable[["Decoy"]] - fdrTable[["DoubleDecoy"]]) /
    #    (scalingFactor * fdrTable[["Target"]])
    #print(fdrTable)
    #return(sprintf("FDR: %0.3f", fdr))
    return(fdr)
}

violationRate <- function(datTab, threshold, maxDistance = 35) {
    datTab <- datTab[datTab$Decoy == "Target" & datTab$dvals >= threshold, ]
    violTab <- table(datTab$distance > maxDistance)
    viol <- violTab[2]
    good <- violTab[1]
    return(viol / (good + viol))
}


# calculateFDRbin <- function(datTab, breaks, scalingFactor=10) {
#     datTab$group <- as.numeric(cut(datTab$dvals, breaks))
#     result <- as.numeric(by(datTab, datTab$group, calculateFDR, threshold=-100))
#     return(abs(result))
# }


setGeneric("prepPairs", function(object, threshold, classifier) {
    standardGeneric("prepPairs")
}
)

setMethod("prepPairs", "PPsearchCompareXL",
          function(object, threshold, classifier){
              object@dataTable <- thresholdResults(object@dataTable, threshold, classifier)
              object@dataTable <- removeDecoys(object@dataTable)
              return(object)
          }
)

removeDecoys <- function(datTab) {
    datTabR <- datTab[datTab$Decoy=="Target",]
    return(datTabR)
}

thresholdResults <- function(datTab, threshold, classifier="Score.Diff") {
    datTabR <- datTab[datTab[[classifier]] >= threshold,]
    return(datTabR)
}

setGeneric("splitSources", function(object, sourceTab) {
    standardGeneric("splitSources")
}
)

setMethod("splitSources", "PPsearchCompareXL",
          #sourceTab is a list that maps Fractions onto expSources
          #st <- list("source1" = c("Fraction1, Fraction2, ...), "source2" = c(Fract))
          function(object, sourceTab){
              datTab <- object@dataTable 
              splitResult <- lapply(sourceTab, function (x) {
                  datTab2 <- datTab[datTab$Fraction %in% x, ]
                  obj <- object
                  obj@dataTable <- datTab2
                  obj }
              )
              for (result in names(splitResult)) {
                  splitResult[[result]]@expSource <- result
              }
              return(splitResult)
          }
)

setGeneric("renameProtein", function(object, oldName, newName) {
    standardGeneric("renameProtein")
}
)

setMethod("renameProtein", "PPsearchCompareXL",
          function(object, oldName, newName){
              object@dataTable <- 
                  .renameProtein(object@dataTable, oldName, newName)
              object@modulFile <- 
                  .renameModTable(object@modulFile, oldName, newName)
              return(object)
          }
)

.renameProtein <- function(datTab, oldName, newName) {
    datTab$Acc.1 <- gsub(oldName, newName, datTab$Acc.1)
    datTab$Acc.2 <- gsub(oldName, newName, datTab$Acc.2)
    datTab <- calculatePairs(datTab)
    return(datTab)
}

.renameModTable <- function(modTab, oldName, newName) {
    names(modTab) <- gsub(oldName, newName, names(modTab))
    return(modTab)
}

setGeneric("renumberProtein", function(object, protein, shift) {
    standardGeneric("renumberProtein")
}
)

setMethod("renumberProtein", "PPsearchCompareXL",
          function(object, protein, shift){
              dt <- object@dataTable
              dt[dt$Acc.1==protein,"XLink.AA.1"] <-
                  dt[dt$Acc.1==protein,"XLink.AA.1"] + shift
              dt[dt$Acc.2==protein,"XLink.AA.2"] <-
                  dt[dt$Acc.2==protein,"XLink.AA.2"] + shift
              dt <- populateModules(dt, object@modulFile)
              dt <- calculatePairs(dt)
              if (is.data.frame(object@caTracedPDB)) {
                  dt <- measureDistances(dt,object@caTracedPDB,object@chainMap)
              }
              object@dataTable <- dt
              return(object)
          }
)


setGeneric("compareTwo", function(data1, data2) {
    standardGeneric("compareTwo")
}
)

setMethod("compareTwo", "PPsearchCompareXL",
          function(data1, data2){
              #Check that reduction state of two dataTables is the same...
              #Maybe several different methods for doing the comparison:
              if ((data1@reductionState != "xlinkedResPairs") | 
                  (data2@reductionState != "xlinkedResPairs")) {
                  stop("This function only works on xlinkedResPair")
              }
              columnsToKeep <- c("xlinkedResPair","XLink.AA.1","Acc.1",
                                 "Modul.1","XLink.AA.2","Acc.2","Modul.2",
                                 "Score.Diff","Xlinker","num","distance")
              dt2 <- data2@dataTable[,columnsToKeep]
              dt1 <- data1@dataTable[,columnsToKeep]
              dt3 <- merge(dt1,dt2,by=columnsToKeep[1:7],all=T)
              dt3$Score.Diff <- apply(dt3[,c("Score.Diff.x","Score.Diff.y")],1, max, na.rm=T)
              dt3$distance <- apply(dt3[,c("distance.x","distance.y")],1, mean, na.rm=T)
#              dt3$num.x[is.na(dt3$num.x)] <- 1
#              dt3$num.y[is.na(dt3$num.y)] <- 1
#              dt3$num <- log(dt3$num.x / dt3$num.y, 2)
#              dt3$num <- scale(dt3$num)
              dt3$num.x[is.na(dt3$num.x)] <- 0
              dt3$num.y[is.na(dt3$num.y)] <- 0
              dt3$num <- dt3$num.x + dt3$num.y
              dt3 <- dt3[,-c(8:15)]
              data1@dataTable <- dt3
              return(data1)
          }
)



setGeneric("pairPlot", function(object, color="lightseagreen",
                                removeMods=NA_character_, displayEmpty=T) {
    standardGeneric("pairPlot")
}
)

setMethod("pairPlot", "PPsearchCompareXL",
          function(object, color, removeMods, displayEmpty){
              .pairPlot(object@dataTable, object@modulFile, color, 
                        removeMods, displayEmpty)
          }
)

.rejiggerMods <- function(modTab) {
    #modTab <- lapply(modTab, function(x) x[-nrow(x),])
    dummy <- function(x) {
        first <- x[1,]
        first[,2] <- first[,1]
        return(rbind(first, x))
    }
    modTab <- lapply(modTab, dummy)
    modTab <- do.call(rbind, modTab)
    modTab$Acc.1 <- gsub("\\..+$","", row.names(modTab))
    row.names(modTab) <- NULL
    names(modTab) <- c("Range_low","XLink.AA.1","Module","Acc.1")
    modTab <- modTab[,c(2,4)]
    modIndex <- expand.grid(1:nrow(modTab),1:nrow(modTab))
    modTab <- cbind(modTab[modIndex$Var1,], modTab[modIndex$Var2,])
    names(modTab) <- c("XLink.AA.1","Acc.1","XLink.AA.2","Acc.2")
    return(modTab)
}

.pairPlot <- function(datTab, modTab, color="lightseagreen", 
                      removeMods=NA_character_, displayEmpty=T) {
    #Function takes a dataTable which is reduced to Residue Pairs and plot
    #Dot for each unique crosslink with size equal to num of instances.  The
    #color can reflect the source of the crosslink to facilitate comparisons.
    #Basic skeleton of the function that will get built up later.  Ultimately
    #would like this to work on the XLSearchOutput object, which probably needs
    #flag to show if it is residue pair reduced etc.
    require(ggplot2)
    require(colorspace)
    datTab <- removeModule(datTab, removeMods)
    datTabR <- datTab
    tempXLinkAA1 <- datTab$XLink.AA.1
    tempAcc1 <- datTab$Acc.1
    datTabR$XLink.AA.1 <- datTab$XLink.AA.2
    datTabR$Acc.1 <- datTab$Acc.2
    datTabR$XLink.AA.2 <- tempXLinkAA1
    datTabR$Acc.2 <- tempAcc1
    datTab <- rbind(datTab, datTabR)
    if (!displayEmpty) {
        modTab <- modTab[names(modTab) %in% 
                             unique(c(datTab$Acc.1, datTab$Acc.2))]
    }
    modTab <- .rejiggerMods(modTab)
    gg <- ggplot(datTab, aes(XLink.AA.1, XLink.AA.2)) +
        #        geom_point(aes(size=num,alpha=Score.Diff,col=distance >35)) +
        geom_vline(aes(xintercept=XLink.AA.1), modTab, color="palevioletred2",
                   linetype="21", size=0.75, alpha=0.5) +
        geom_hline(aes(yintercept=XLink.AA.2), modTab, color="palevioletred2",
                   linetype="21", size=0.75, alpha=0.5) +
        facet_grid(Acc.2 ~ Acc.1, space="free", scales="free", 
                   as.table=F, switch="both") + 
        theme(panel.background = element_rect(fill=rgb(0.98,0.98,0.98)),
              panel.spacing = unit(0.5, units="line"),
              panel.border = element_rect(color = "grey", fill = NA, size = 1),
              panel.grid.major = element_line(color = "grey", size=.5),
              panel.grid.minor = element_line(color = "grey", size=.25),
              aspect.ratio=1) +#legend.position="none") +
        scale_x_continuous(breaks=seq(0,max(datTab$XLink.AA.1),by=50),
                           minor_breaks=seq(0,max(datTab$XLink.AA.1),by=25)) +
        scale_y_continuous(breaks=seq(0,max(datTab$XLink.AA.2),by=50),
                           minor_breaks=seq(0,max(datTab$XLink.AA.2),by=25)) +
        scale_size_area(limits=c(1,70)) +
        scale_alpha(limits=c(10,55), range=c(0.8,1))
    
    if (color=="byDistance") {gg <- gg + geom_point(aes(size=num,alpha=Score.Diff,col=distance <35))
    }   else if (color=="byQuant") {gg <- gg + geom_point(aes(size=numSC, alpha=Score.Diff,col=num)) +
        #scale_color_gradientn(breaks=c(-10,-3,-2,-1,0,1,2,3,10),colors=diverge_hsv(8))
        scale_color_gradient2(low=muted("blue"),high=muted("red"),mid="grey70",midpoint=0,limits=c(-1.5,1.5))
    }   else {gg <- gg + geom_point(aes(size=num,alpha=Score.Diff),col=color)
    }
    plot(gg)
}

removeModule <- function(datTab, modules) {
    datTab <- datTab[!(datTab$Modul.1 %in% modules | datTab$Modul.2 %in% modules),]
    return(datTab)
}

setGeneric("parseCrosslinker", function(object) {
    standardGeneric("parseCrosslinker")
}
)

setMethod("parseCrosslinker", "PPsearchCompareXL",
          function(object){
              datTab <- object@dataTable
              datTab <- .parseCrosslinker(datTab)
              object@dataTable <- datTab
              return(object)
          }
)

.parseCrosslinker <- function(datTab) {
    datTab$Xlinker <- "BS3"
    datTab[grepl("\\:2H\\(12\\)",datTab$Peptide.1) | 
               grepl("\\:2H\\(12\\)",datTab$Peptide.2),"Xlinker"] <- "BS3hvy"
    datTab$Xlinker <- as.factor(datTab$Xlinker)
    return(datTab)
}

setGeneric("buildClassifier", function(object) {
    standardGeneric("buildClassifier")
}
)

setMethod("buildClassifier", "PPsearchCompareXL",
          function(object){
              require(e1071)
              datTab <- object@dataTable
              ind <- sample(nrow(datTab),nrow(datTab)/2)
              train <- datTab[ind,]
              test <- datTab[-1 * ind,]
              params <- c("percMatch","Score.Diff","z","Rk.2","Rk.1","Len.Pep.1","Len.Pep.2")
              wghts <- numeric(0)
              wghts["Target"] <- table(test$Decoy2)[1] / sum(table(test$Decoy2),na.rm=T)
              wghts["Decoy"] <- table(test$Decoy2)[2] / sum(table(test$Decoy2),na.rm=T)
              fit <- svm(train$Decoy2 ~., 
                         subset(train, select=params),
                         kernel="linear",
                         cost=0.01,
                         tolerance=0.01,
                         class.weights=wghts
              )
              p <- predict(fit, subset(datTab,select=params),decision.values=T)
              datTab$dvals <- attr(p,"decision.values")
              print(paste("training weights:", round(wghts,2)))
              tab <- table(datTab$Decoy2, datTab$dvals > 0)
              if (tab[1] < tab[3]) {
                  datTab$dvals <- -1 * datTab$dvals
              }
              print(tab)
              print(paste("specificity:", round(tab[1]/(tab[1]+tab[3]),2)))
              object@dataTable <- datTab
              return(object)
          }
)

setGeneric("makeFilteredPeaklists", function(object, inputPeaklistDir, outputPeaklistDir) {
    standardGeneric("makeFilteredPeaklists")
}
)

setMethod("makeFilteredPeaklists", "PPsearchCompareXL",
          function(object, inputPeaklistDir, outputPeaklistDir){
              cwd <- getwd()
              setwd(inputPeaklistDir)
              datTab <- object@dataTable
              fileNames <- unique(datTab$Fraction)
              for (file in fileNames) {
                  subTab <- datTab[datTab$Fraction == file,"MSMS.Info"]
                  spectra <- sapply(subTab, extractSpecFromMGFpd, file)
                  outFileName <- paste(outputPeaklistDir, "/", file, sep="")
                  spectra <- spectra[order(subTab)]
                  cat(spectra, file=outFileName, sep="\n")
              }
              datTab <- formatMSViewerFile(datTab)
              outFileName <- paste(outputPeaklistDir, "/", "filteredPeakList.txt", sep="")
              write.table(datTab,outFileName,sep="\t",quote=F,row.names=F)
              setwd(cwd)
          }
)

formatMSViewerFile <- function(datTab) {
    require(stringr)
    datTab$percMatch <- round(datTab$percMatch * 100, 2)
    datTab$dvals <- round(datTab$dvals, 2)
    datTab$Spectrum <- 1
    datTab <- datTab[order(datTab$Fraction, datTab$MSMS.Info),]
    dupeDF <- datTab[,c("Fraction", "RT")]
    dupeDF$RT <- round(dupeDF$RT * 60, 0)
    dupes <- which(duplicated(dupeDF))
    datTab[dupes, "Spectrum"] <- 2
    datTab$Fraction <- gsub("(.*)\\.[^.]+$", "\\1", datTab$Fraction)
    names(datTab) <- gsub("Peptide\\.([[:digit:]])", "Peptide \\1", names(datTab))
    names(datTab) <- gsub("m.z", "m/z", names(datTab))
    annoyingColumns <- str_which(names(datTab), "(Int|Dec)[a-z]{2}\\.[[1-2]]")
    datTab <- datTab[, -annoyingColumns]
    datTab <- datTab[order(datTab$dvals, decreasing = T),]
    return(datTab)
}

extractSpecFromMGFpd <- function(scanNo, mgfFile) {
    require(readr)
    print(paste(scanNo, mgfFile, sep="\t"))
    startLine <- getLineNos(mgfFile, scanNo) - 5
    spectrum <- read_lines(mgfFile, skip= startLine -1, n_max=10000)
    endLine <- grep("END IONS", spectrum)[1]
    spectrum <- paste(c(spectrum[1:endLine],""), sep="", collapse="\n")
    return(spectrum)
}

getLineNos <- function(mgfFile, scanNo) {
    scanQuote <- shQuote(paste("SCANS=", scanNo, "(\\r|\\n)", sep=""))
    grepOut <- system2("egrep", args=c("-n", scanQuote, mgfFile), stdout=T)
    lineNo <- gsub(":SCANS=.+", "", grepOut)
    return(as.integer(lineNo))
}


