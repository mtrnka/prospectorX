print("touchStone module loaded")

readProspectorXLOutput <- function(inputFile){
    dataTable <- read_tsv(inputFile)
    header <- names(dataTable) %>%
        str_replace("_[[0-9]]$", "") %>%
        str_replace_all("[[:space:]]",".") %>%
        str_replace_all("#", "Num")
    acc_pos <- str_which(header, "Acc")
    header[acc_pos] <- c("Acc.1", "Acc.2")
    spec_pos <- str_which(header, "Species")
    header[spec_pos] <- c("Species.1", "Species.2")
    prot_pos <- str_which(header, "Protein.Name")
    header[prot_pos] <- c("Protein.1", "Protein.2")
    if ("Protein.MW" %in% header) {
        mw_pos <- str_which(header, "Protein.MW")
        header[mw_pos] <- c("MW.1", "MW.2")
    }
    names(dataTable) <- header
    if (!"Spectrum" %in% names(dataTable)) {
        dataTable$Spectrum <- 1
    }
    if (!"distance" %in% names(dataTable)) {
        dataTable$distance <- NA_real_
    }
    dataTable <- calculateDecoys(dataTable)
    dataTable <- calculatePairs(dataTable)
    dataTable <- assignXLinkClass(dataTable)
    dataTable <- calculatePercentMatched(dataTable)
    dataTable <- calculatePeptideLengths(dataTable)
    dataTable <- lengthFilter(dataTable, minLen = 3, maxLen = 35)
    dataTable <- scoreFilter(dataTable, minScore = 0)
    return(dataTable)
}

bestResPair <- function(datTab, classifier=c("Score.Diff", "Score", "MSMS.Info")) {
    datTabR <- datTab
    xlinks <- unique(datTab$xlinkedResPair)
    for (param in classifier) {
        if (param %in% names(datTabR)) {
            datTabR <- reduceTo(datTabR, "xlinkedResPair", param)
        }
    }
    cat(c(nrow(datTabR), length(xlinks),
          sum(datTabR$num), nrow(datTab), "\n"), sep="\t")
    return(datTabR)
}

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

# readModuleFile <- function(modFile) {
#     inFile <- file(modFile, open="r")
#     header <- readLines(inFile, n=1)
#     header <- lineSplit(header)
#     dataTable <- character(5)
#     while (T) {
#         line <- readLines(inFile, n=1)
#         if (length(line)==0) {break}
#         line <- lineSplit(line)
#         if (length(line)==0) {break}
#         else if (length(line) < 5) {length(line)=5}
#         dataTable <- rbind(dataTable,line)}
#     close(inFile)
#     row.names(dataTable) <- NULL
#     dataTable <- as.data.frame(dataTable[-1,],stringsAsFactors=F)
#     names(dataTable) <- header
#     dataTable <- cleanTypes(dataTable)
#     modules <- list()
#     for (subunit in unique(dataTable$Subunit)) {
#         conv <- dataTable[dataTable$Subunit==subunit,
#                           c("Range_low","Range_high","Module")]
#         modules[[subunit]] <- conv[order(conv$Range_low),]
#     }
#     return(modules)
# }

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
    print("***parsePDB***")
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
        #print(paste(residue1,chains1,residue2,chains2,sep="\t"))
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
#    searchTable <- searchTable[searchTable$Decoy == "Target",]
    print("***measureDistances***")
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
    print("*** measured distances on pdb file ***")
    return(searchTable)
}

# populateModules <- function(searchTable, moduleDefinitions) {
#     print("***populate Modules***")
#     modules <- list(moduleDefinitions)
#     pep1.modul <- mapply(
#         assignModule,
#         searchTable$XLink.AA.1,
#         searchTable$Acc.1,
#         modules)
#     pep2.modul <- mapply(
#         assignModule,
#         searchTable$XLink.AA.2,
#         searchTable$Acc.2,
#         modules)
#     levs <- unlist(sapply(moduleDefinitions, function(x) c(x[3])))
#     names(levs) <- NULL
#     levs <- unique(levs)
#     searchTable$Modul.1 <- factor(unlist(pep1.modul),levels=levs)
#     searchTable$Modul.2 <- factor(unlist(pep2.modul),levels=levs)
#     return(searchTable)
# }
# 
# assignModule <- function(seqPosition,protein,modList) {
#     if (protein %in% names(modList)) {
#         n <- findInterval(as.integer(seqPosition),
#                           modList[[as.character(protein)]]$Range_low)
#         if (n != 0) {
#             return(modList[[protein]]$Module[n])
#         } else return(NA)
#     }
#     else return(NA)
# }

assignModules <- function(searchTable, moduleFile) {
    modTab <- read_tsv(moduleFile)
    Modul.1 <- searchTable %>% 
        left_join(modTab, by=c("Acc.1"="Subunit")) %>%
        filter(is.na(Range_low) | (XLink.AA.1 >= Range_low & XLink.AA.1 <= Range_high)) %>%
        pull("Module")
    Modul.2 <- searchTable %>% 
        left_join(modTab, by=c("Acc.2"="Subunit")) %>%
        filter(is.na(Range_low) | (XLink.AA.2 >= Range_low & XLink.AA.2 <= Range_high)) %>%
        pull("Module")
    searchTable$Modul.1 <- as_factor(Modul.1)
    searchTable$Modul.2 <- as_factor(Modul.2)
    return(searchTable)
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
    searchTable <- searchTable %>% add_count(xlinkedResPair, name="numCSM")
    return(searchTable)
}

assignXLinkClass <- function(searchTable) {
    intraProteinLinks <- searchTable$Acc.1 == searchTable$Acc.2
    interProteinLinks <- searchTable$Acc.1 != searchTable$Acc.2
    s1 <- searchTable$Start.1
    s2 <- searchTable$Start.2
    e1 <- searchTable$End.1
    e2 <- searchTable$End.2
    condition1 <- (s1 >= s2) & (s1 <= e2)
    condition2 <- (e1 >= s2) & (e1 <= e2)
    interHomomerLinks <- condition1 | condition2
    searchTable$xlinkClass <- ifelse(intraProteinLinks,
                                     ifelse(interHomomerLinks,
                                     "interProtein, homomeric",
                                     "intraProtein"),
                                     "interProtein, heteromeric")
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

# lineSplit <- function(line) {
#     l <- nchar(line)
#     g <- gregexpr("\t",line)[[1]]
#     gl <- g[length(g)]
#     if (l == gl) {
#         return(c(unlist(strsplit(line,"\t")),""))
#     } else {
#         return(unlist(strsplit(line,"\t")))
#     }
# }

cleanTypes <- function(dataTable) {
    return(as.data.frame(
        lapply(dataTable, function(x) {
            type.convert(x, as.is=TRUE)}),
        stringsAsFactors=FALSE))
}

calculateFDR <- function(datTab, threshold=-100, classifier="dvals", scalingFactor=10) {
    datTab <- datTab[datTab[[classifier]] >= threshold, ]
    fdrTable <- table(datTab$Decoy)
    if (is.na(fdrTable["DoubleDecoy"])) {fdrTable["DoubleDecoy"] <- 0}
    if (is.na(fdrTable["Decoy"])) {fdrTable["Decoy"] <- 0}
    if (is.na(fdrTable["Target"])) {fdrTable["Target"] <- 0}
    ffTT <- fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2)
    ftTT <- (fdrTable[["Decoy"]] / scalingFactor) - 
        (2 * fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2))
    fdr <- (ffTT + ftTT) / fdrTable[["Target"]]
    return(fdr)
}

calculateHits <- function(datTab, threshold=-100, classifier="dvals") {
    datTab <- datTab[datTab[[classifier]] >= threshold, ]
    intra <- sum(str_count(datTab$xlinkClass, "intraProtein"))
    inter <- sum(str_count(datTab$xlinkClass, "interProtein"))
    return(c("thresh"=threshold, "intra"=intra, "inter"=inter))
}

calculateGT <- function(datTab, threshold=-100, classifier="dvals") {
    datTab <- datTab[datTab[[classifier]] >= threshold & 
                         datTab$Decoy=="Target", ]
    gtTable <- table(datTab$groundTruth)
    if (is.na(gtTable["FALSE"])) {gtTable["FALSE"] <- 0}
    if (is.na(gtTable["TRUE"])) {gtTable["TRUE"] <- 0}
    gt <- gtTable["FALSE"] / sum(gtTable)
    return(gt)
}

violationRate <- function(datTab, threshold, classifier="Score.Diff", maxDistance = 35) {
    datTab <- datTab[datTab[[classifier]] >= threshold,]
    dists <- datTab$distance
    dists <- dists[!is.na(dists)]
    nAbove <- sum(dists > maxDistance)
    nBelow <- sum(dists <=maxDistance)
    return(nAbove / (nAbove + nBelow))
}

generateErrorTable <- function(datTab, classifier="dvals",
                               errorFUN=calculateFDR, ...) {
    class.max <- ceiling(max(datTab[[classifier]]))
    class.min <- floor(min(datTab[[classifier]]))
    if (classifier=="dvals") {
        range.spacing = 0.1
    } else if (classifier=="Score.Diff") {
        range.spacing = 0.25
    } else {
        range.spacing = abs(class.max - class.min) / 100
    }
    class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
    error.rates <- unlist(lapply(class.range, function(threshold) {
        errorFUN(datTab, threshold=threshold, classifier=classifier, ...)
    }))
    error.rates[is.na(error.rates)] <- 0
    #plot(error.rates ~ class.range, type="l")
    
    num.hits <- map_dfr(class.range, function(threshold) {
        calculateHits(removeDecoys(datTab), threshold=threshold, classifier=classifier, ...)
    })
    num.hits$fdr <- error.rates
    num.hits$total <- num.hits$inter + num.hits$intra
    return(num.hits)
}

findThreshold <- function(datTab, targetER=0.05, minThreshold=-5,
                          classifier="dvals", errorFUN=calculateFDR, ...) {
    num.hits <- generateErrorTable(datTab, classifier, errorFUN, ...)
    error.rates <- num.hits$fdr
    class.range <- num.hits$thresh
    error.diff <- (error.rates - targetER)
    error.cross <- logical(length(class.range))
    for (i in 1:(length(class.range) - 1)) {
        error.cross[i] <- error.diff[i] > 0 & error.diff[i+1] <= 0
    }
    if (sum(error.cross)==0 & error.diff[1] < 0) {
        threshold=minThreshold
    } else if (sum(error.cross)==0 & error.diff[1] > 0) {
        threshold=class.max
        print("no acceptable threshold found")
    } else {
        firstCross <- which(error.cross)[1] + 1
        threshold <- class.range[firstCross]
    }
    if (threshold < minThreshold) {
        threshold = minThreshold
    }
    return(list(threshold, errorFUN(datTab, threshold, ...)))
}

removeDecoys <- function(datTab) {
    datTabR <- datTab[datTab$Decoy=="Target",]
    return(datTabR)
}

thresholdResults <- function(datTab, threshold, classifier="Score.Diff") {
    datTabR <- datTab[datTab[[classifier]] >= threshold,]
    return(datTabR)
}

renameProtein <- function(datTab, oldName, newName) {
    datTab$Acc.1 <- gsub(oldName, newName, datTab$Acc.1)
    datTab$Acc.2 <- gsub(oldName, newName, datTab$Acc.2)
    datTab <- calculatePairs(datTab)
    return(datTab)
}

renameModTable <- function(modTab, oldName, newName) {
    names(modTab) <- gsub(oldName, newName, names(modTab))
    return(modTab)
}

renumberProtein <-function(dt, protein, shift) {
    dt[dt$Acc.1==protein,"XLink.AA.1"] <-
        dt[dt$Acc.1==protein,"XLink.AA.1"] + shift
    dt[dt$Acc.2==protein,"XLink.AA.2"] <-
        dt[dt$Acc.2==protein,"XLink.AA.2"] + shift
    dt <- populateModules(dt, object@modulFile)
    dt <- calculatePairs(dt)
    if (is.data.frame(object@caTracedPDB)) {
        dt <- measureDistances(dt,object@caTracedPDB,object@chainMap)
    }
    return(dt)
}

# setGeneric("pairPlot", function(object, color="lightseagreen",
#                                 removeMods=NA_character_, displayEmpty=T) {
#     standardGeneric("pairPlot")
# }
# )
# 
# setMethod("pairPlot", "PPsearchCompareXL",
#           function(object, color, removeMods, displayEmpty){
#               .pairPlot(object@dataTable, object@modulFile, color, 
#                         removeMods, displayEmpty)
#           }
# )

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
#    datTab <- removeModule(datTab, removeMods)
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
                geom_point(aes(size=numCSM,alpha=Score.Diff,col=color)) +
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
    
    if (color=="byDistance") {
        gg <- gg + geom_point(aes(size=numCSM,alpha=Score.Diff,col=distance <35))
    }   else if (color=="byQuant") {
        gg <- gg + geom_point(aes(size=numCSM, alpha=Score.Diff,col=numCSM)) +
        #scale_color_gradientn(breaks=c(-10,-3,-2,-1,0,1,2,3,10),colors=diverge_hsv(8))
        scale_color_gradient2(low=muted("blue"), high=muted("red"), mid="grey70", midpoint=0, limits=c(-1.5,1.5))
    }   else {
        gg <- gg + geom_point(aes(size=numCSM, alpha=Score.Diff), col=color)
    }
    plot(gg)
}

removeModule <- function(datTab, modules) {
    datTab <- datTab[!(datTab$Modul.1 %in% modules | datTab$Modul.2 %in% modules),]
    return(datTab)
}

parseCrosslinker <- function(datTab) {
    datTab$Xlinker <- "BS3"
    datTab[grepl("\\:2H\\(12\\)",datTab$Peptide.1) | 
               grepl("\\:2H\\(12\\)",datTab$Peptide.2),"Xlinker"] <- "BS3hvy"
    datTab$Xlinker <- as.factor(datTab$Xlinker)
    return(datTab)
}

buildClassifier <- function(datTab) {
    datTab$massError <- abs(datTab$ppm)
    if (nrow(datTab) > 30000) {
        sampleNo <- 15000
    } else {
        sampleNo <- nrow(datTab)/2
    }
    ind <- sample(nrow(datTab),sampleNo)
    train <- datTab[ind,]
    test <- datTab[-1 * ind,]
    params <- c("Score.Diff","z","Score","numCSM","massError", "Rk.2","Rk.1")
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
    datTab$dvals <- as.numeric(attr(p,"decision.values"))
    print(paste("training weights:", round(wghts,2)))
    tab <- table(datTab$Decoy2, datTab$dvals > 0)
    if (tab[1] < tab[3]) {
        datTab$dvals <- -1 * datTab$dvals
        tab <- table(datTab$Decoy2, datTab$dvals > 0)
    }
    print(tab)
    print(paste("specificity:", round(tab[1]/(tab[1]+tab[3]),2)))
    return(datTab)
}

# setGeneric("makeFilteredPeaklists", function(object, inputPeaklistDir, outputPeaklistDir) {
#     standardGeneric("makeFilteredPeaklists")
# }
# )
# 
# setMethod("makeFilteredPeaklists", "PPsearchCompareXL",
#           function(object, inputPeaklistDir, outputPeaklistDir){
#               cwd <- getwd()
#               setwd(inputPeaklistDir)
#               datTab <- object@dataTable
#               fileNames <- unique(datTab$Fraction)
#               for (file in fileNames) {
#                   subTab <- datTab[datTab$Fraction == file,"MSMS.Info"]
#                   spectra <- sapply(subTab, extractSpecFromMGFpd, file)
#                   outFileName <- paste(outputPeaklistDir, "/", file, sep="")
#                   spectra <- spectra[order(subTab)]
#                   cat(spectra, file=outFileName, sep="\n")
#               }
#               datTab <- formatMSViewerFile(datTab)
#               outFileName <- paste(outputPeaklistDir, "/", "filteredPeakList.txt", sep="")
#               write.table(datTab,outFileName,sep="\t",quote=F,row.names=F)
#               setwd(cwd)
#           }
# )

formatMSViewerFile <- function(datTab) {
    require(stringr)
    if ("percMatch" %in% names(datTab)) {
        datTab$percMatch <- round(datTab$percMatch * 100, 2)
    }
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

generateMSViewerLink <- function(path, fraction, rt, z, peptide.1, peptide.2, spectrum) {
    if(!str_detect(fraction, "\\.[[a-z]]$")) {
        fraction <- paste0(fraction, ".mgf")
    }
    instrumentType <- switch(
        str_extract(fraction, "(?<=FTMSms2)[[a-zA-Z]]+"),
        "ESI-Q-high-res",
        ethcd = "ESI-EThcD-high-res",
        etd = "ESI-ETD-high-res",
        hcd = "ESI-Q-high-res" #Add other instrument types
    )
    linkType <- str_extract(peptide.1, "(?<=\\(\\+).+?(?=\\))")
    templateVals[6] <- file.path(path, fraction)
    templateVals[9] <- instrumentType
    templateVals[18] <- rt
    templateVals[20] <- spectrum
    templateVals[21] <- z
    templateVals[22] <- z
    templateVals[23] <- peptide.1
    templateVals[25] <- peptide.2
    templateVals[28] <- linkType
    templateVals <- url_encode(templateVals)
    zipped <- str_c(templateKeys[1:28], templateVals[1:28], sep="=", collapse="&")
    str_c('<a href=\"', zipped, '\" target=\"_blank\">Spectrum</a>')
}

generateCheckBoxes <- function(datTab) {
    datTab$selected <- str_c('<input type="checkbox" names="row', datTab$xlinkedResPair, '"value="', datTab$xlinkedResPair, '">',"")
    return(datTab)
}

formatXLTable <- function(datTab) {
    annoyingColumns <- str_which(names(datTab), "(Int|Dec)[a-z]{2}\\.[[1-2]]")
    datTab <- datTab[, -annoyingColumns]
    datTab <- datTab %>%
        select(-starts_with("Res"),
               -starts_with("Decoy"),
               -starts_with("Num\\."),
               -massError,
               -xlinkedPepPair)
    if (sum(!is.na(datTab$distance)) == 0) {
        datTab <- datTab %>%
            select(-distance)
    }
    datTab <- datTab[order(datTab$dvals, decreasing = T),]
    datTab <- datTab %>% select(any_of(c("selected", "link", "xlinkedResPair", 
                                "dvals", "distance", "m.z", "z", "ppm",
                                "DB.Peptide.1", "DB.Peptide.2", "Score", "Score.Diff",
                                "Sc.1", "Rk.1", "Sc.2", "Rk.2", "Acc.1",
                                "XLink.AA.1", "Protein.1", "Modul.1", "Species.1",
                                "Acc.2", "XLink.AA.2", "Protein.2", "Modul.2", "Species.2",
                                "xlinkClass", "Len.Pep.1", "Len.Pep.2",
                                "Peptide.1", "Peptide.2", "NumCSM",
                                "Fraction", "RT", "MSMS.Info")))
    if ("Protein.1" %in% names(datTab) & "Protein.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Protein.1 = factor(Protein.1), 
                                    Protein.2 = factor(Protein.2))
        fct_unify(list(datTab$Protein.1, datTab$Protein.2))
    }
    if ("Modul.1" %in% names(datTab) & "Modul.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Modul.1 = factor(Modul.1), 
                                    Modul.2 = factor(Modul.2))
        fct_unify(list(datTab$Modul.1, datTab$Modul.2))
    }
    if ("Species.1" %in% names(datTab) & "Species.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Species.1 = factor(Species.1), 
                                    Species.2 = factor(Species.2))
        fct_unify(list(datTab$Species.1, datTab$Species.2))
    }
    if ("Acc.1" %in% names(datTab) & "Acc.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Acc.1 = factor(Acc.1), 
                                    Acc.2 = factor(Acc.2))
        fct_unify(list(datTab$Acc.1, datTab$Acc.2))
    }
    if ("xlinkClass" %in% names(datTab)) {
        datTab <- datTab %>% mutate(xlinkClass = factor(xlinkClass))
    }
    if ("Fraction" %in% names(datTab)) {
        datTab$Fraction <- gsub("(.*)\\.[^.]+$", "\\1", datTab$Fraction)
        datTab <- datTab %>% mutate(xlinkClass = factor(Fraction))
    }
    if ("percMatch" %in% names(datTab)) {
        datTab$percMatch <- round(datTab$percMatch * 100, 2)
    }
    if ("dvals" %in% names(datTab)) {
        datTab$dvals <- round(datTab$dvals, 2)
        names(datTab) <- gsub("dvals", "SVM.score", names(datTab))
    }
    if ("Rk.1" %in% names(datTab) & "Rk.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Rk.1 = as.integer(Rk.1), 
                                    Rk.2 = as.integer(Rk.2))
    }
    if ("XLink.AA.1" %in% names(datTab) & "XLink.AA.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(XLink.AA.1 = as.integer(XLink.AA.1), 
                                    XLink.AA.2 = as.integer(XLink.AA.2))
    }
    if ("Len.Pep.1" %in% names(datTab) & "Len.Pep.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Len.Pep.1 = as.integer(Len.Pep.1), 
                                    Len.Pep.2 = as.integer(Len.Pep.2))
    }
    if ("z" %in% names(datTab)) {
        datTab <- datTab %>% mutate(z = as.integer(z))
    }
    if ("MSMS.Info" %in% names(datTab)) {
        datTab <- datTab %>% mutate(MSMS.Info = as.integer(MSMS.Info))
    }
    datTab <- datTab %>% mutate(xlinkedResPair = fct_drop(xlinkedResPair))
    return(datTab)
}

XVisOutput <- function(datTab) {
#    Protein1 <- str_c("sp", pull(datTab, Acc.1), pull(datTab, Species.1), sep = "|")
#    Protein2 <- str_c("sp", pull(datTab, Acc.2), pull(datTab, Species.2), sep = "|")
    Protein1 <- pull(datTab, Protein.1)
    Protein2 <- pull(datTab, Protein.2)
    AbsPos1 <- pull(datTab, XLink.AA.1)
    AbsPos2 <- pull(datTab, XLink.AA.2)
    IDScore <- pull(datTab, SVM.score)
    xVisXL <- data.frame("Protein1" = Protein1, "Protein2" = Protein2, 
                         "AbsPos1" = AbsPos1, "AbsPos2"= AbsPos2, "Id-Score" = IDScore)
}

getLineNos <- function(mgfFile, scanNo) {
    scanQuote <- shQuote(paste("SCANS=", scanNo, "(\\r|\\n)", sep=""))
    grepOut <- system2("egrep", args=c("-n", scanQuote, mgfFile), stdout=T)
    lineNo <- gsub(":SCANS=.+", "", grepOut)
    return(as.integer(lineNo))
}

calculatePrecMass <- function(pepMass, charge) {
    precMass <- pepMass * charge - 1.007276 * charge
    return(precMass)
}

getIndicesForFile <- function(file) {
    mgfFile <- readLines(file)
    scanNo <- na.omit(str_match(mgfFile, "TITLE=Scan\\s([[0-9]]+)"))
    scanNo <- as.integer(scanNo[,2])
    masterScanNo <- na.omit(str_match(mgfFile, "Master_Scan_Number=([[0-9]]+)"))
    masterScanNo <- as.integer(masterScanNo[,2])
    pepMass <- na.omit(str_match(mgfFile, "PEPMASS=([[0-9]]+\\.[[0-9]]+)"))
    pepMass <- as.numeric(pepMass[,2])
    charge <- na.omit(str_match(mgfFile, "CHARGE=([[0-9]])\\+"))
    charge <- as.integer(charge[,2])
    extract <- data.frame(scanNo=scanNo, masterScanNo=masterScanNo, 
                          pepmass=pepMass, charge=charge)
    return(extract)
}

createMasterScanFile <- function(peaklistDir, 
                                 ms3PAVApeaklist, 
                                 ms2PAVApeaklistCID,
                                 ms2PAVApeaklistETD=NA) {

    ms3masterscansCID <- getIndicesForFile(file.path(peaklistDir,
                                                     ms3PAVApeaklist))
    names(ms3masterscansCID) <- c("ms3ScanNo", "ms2cidScanNo", 
                                  "pepmass.ms3", "charge.ms3")
    ms2masterscansCID <- getIndicesForFile(file.path(peaklistDir,
                                                     ms2PAVApeaklistCID))
    if (!is.na(ms2PAVApeaklistETD)) {
        ms2masterscansETD <- getIndicesForFile(file.path(peaklistDir,
                                                         ms2PAVApeaklistETD))
        ms2masterscans <- full_join(ms2masterscansCID, ms2masterscansETD,
                                    by=c("masterScanNo", "pepmass","charge"))
        names(ms2masterscans) <- c("ms2cidScanNo", "ms1MasterScanNo",
                                   "pepmass.ms2", "charge.ms2", "ms2etdScanNo")
    } else {
        ms2masterscans <- ms2masterscansCID
        names(ms2masterscans) <- c("ms2cidScanNo", "ms1MasterScanNo",
                                   "pepmass.ms2", "charge.ms2")
    }
    ms2masterscans <- full_join(ms2masterscans, ms3masterscansCID, 
                                by="ms2cidScanNo")
    ms2masterscans$precMass.ms2 <- map2_dbl(ms2masterscans$pepmass.ms2, 
                                            ms2masterscans$charge.ms2,
                                            calculatePrecMass)
    if (!is.na(ms2PAVApeaklistETD)) {
        ms2masterscans <- ms2masterscans[,c(2,1,5:6,3:4,9,7:8)]
    } else {
        ms2masterscans <- ms2masterscans[,c(2,1,5,3:4,8,6:7)]
    }
    return(ms2masterscans)
}

readMS3results <- function(searchResultsDir, ms3SearchTable, masterScanFile) {
    ms3results <- read_tsv(ms3SearchTable, skip=2) %>% 
        select(c("DB Peptide", "Peptide", "Protein Mods",
               "Fraction", "RT", "Spectrum", "MSMS Info", "Score",
               "Expect", "Protein Name", "Acc #", "Species"))
    names(ms3results) <- c("DB.Peptide.ms3", "Peptide.ms3", "Protein.Mods.ms3",
                           "Fraction.ms3", "RT.ms3", "Spectrum.ms3", "ms3ScanNo", 
                           "Score.ms3", "Expect.ms3", "Protein.ms3", "Acc.ms3",
                           "Decoy.ms3")
    ms3results$Decoy.ms3[ms3results$Decoy.ms3 != "DECOY"] <- "Target"
    ms3results <- left_join(ms3results, masterScanFile, 
                            by = "ms3ScanNo")
    ms3results$precMass.ms3 <- map2_dbl(ms3results$pepmass.ms3,
                                        ms3results$charge.ms3, calculatePrecMass)
    return(ms3results)
}

processMS3xlinkResults <- function(ms3SearchTable,
                                   projectDir,
                                   pdbDir,
                                   chainMapFile,
                                   pdbFile,
                                   modFile,
                                   preProcessFunction=function(x) x) {
    
    ms3results <- ms3SearchTable
    n_ms3 <- ms3results %>% group_by(ms2cidScanNo) %>% nest()
    n_ms3 <- n_ms3 %>% mutate(scanGroupSize = map_int(data, ~nrow(.)))
    
    ms3crosslinksFlat <- map2_dfr(n_ms3$data, n_ms3$ms2cidScanNo, ms3ionFragFind)
    ms3crosslinksFlat <- preProcessFunction(ms3crosslinksFlat)
    ms3crosslinksFlat$Score.Diff <- ms3crosslinksFlat$Sc.2
    ms3crosslinksFlat <- calculatePairs(ms3crosslinksFlat)
    ms3crosslinksFlat <- assignXLinkClass(ms3crosslinksFlat)
    #ms3crosslinksFlat <- calculatePercentMatched(ms3crosslinksFlat)
    #ms3crosslinksFlat <- calculatePeptideLengths(ms3crosslinksFlat)
    ms3crosslinksFlat <- lengthFilter(ms3crosslinksFlat, minLen = 4, maxLen = 25)
    ms3crosslinksFlat <- scoreFilter(ms3crosslinksFlat, minScore = 0)
    
    # Modify this class or make new class to read in the MS3 results and all 
    # necessary peaklists and flatten.
    t1 <- new(Class="PPsearchCompareXL",
              directory=projectDir,
              dataFile="emptyTemplate.txt",
              modFile=modFile,
              pdbDir=pdbDir, 
              pdbFile=pdbFile, 
              chainMapFile=chainMapFile
    )
    t1@dataTable <- ms3crosslinksFlat
    t1@dataTable <- nces(t1@dataTable, t1@caTracedPDB, t1@chainMap)
    t1@dataTable <- populateModules(t1@dataTable, t1@modulFile)
    t1@dataTable <- calculateDecoys(t1@dataTable)
    return(t1)
}

getMassError <- function(targetMass, testMass) {
    # reports mass error in ppm
    massDiff <- 1e6 * (targetMass - testMass) / targetMass
    return(massDiff)
}

ms3ionFragFind <- function(scanGroupTibble, masterScanNo) {
    ms2Mass <- scanGroupTibble$precMass.ms2[1]
    scanGroupSize <- nrow(scanGroupTibble)
    if (scanGroupSize < 2) {
        return(data.frame())
    }
    pairedIons <- combn(scanGroupTibble$precMass.ms3, 2)
    pairedIndex <- combn(1:scanGroupSize, 2)
    summedPairs <- apply(pairedIons, 2, function(x) sum(x) + 18.01002)
    matchesMS2 <- unlist(lapply(summedPairs, getMassError, ms2Mass))
    CSMs <- which(abs(matchesMS2) <= 25) # +/-25 ppm mass error for ms2 ions
    if (length(CSMs) == 0) return(data.frame())
    placeHolder <- vector("list", length(CSMs))
    index = 1
    for (CSM in CSMs) {
        massError <- matchesMS2[CSM]
        peptidePair <- pairedIndex[,CSM]
        flatMS3pair <- flattenMS3pair(scanGroupTibble[peptidePair[1],],
                                      scanGroupTibble[peptidePair[2],],
                                      masterScanNo)
        flatMS3pair$ppm <- massError
        if (is.na(flatMS3pair$XLink.AA.1) | is.na(flatMS3pair$XLink.AA.2)) {
            flatMS3pair <- data.frame()
        }
        placeHolder[[index]] <- flatMS3pair
        index = index + 1
    }
    return(do.call(rbind, placeHolder))
}

flattenMS3pair <- function(ms3CSM1, ms3CSM2, masterScanNo=NA) {
    if (ms3CSM1$Expect.ms3 > ms3CSM2$Expect.ms3) {
        tempCSM <- ms3CSM1
        ms3CSM1 <- ms3CSM2
        ms3CSM2 <- tempCSM
    }
    flatCSM <- tibble(ms2cidScanNo=masterScanNo,
                      MSMS.Info.1=ms3CSM1$`ms3ScanNo`,
                      MSMS.Info.2=ms3CSM2$`ms3ScanNo`,
                      DB.Peptide.1=ms3CSM1$`DB.Peptide.ms3`,
                      DB.Peptide.2=ms3CSM2$`DB.Peptide.ms3`,
                      Len.Pep.1=nchar(ms3CSM1$`DB.Peptide.ms3`),
                      Len.Pep.2=nchar(ms3CSM2$`DB.Peptide.ms3`),
                      Peptide.1=ms3CSM1$Peptide.ms3,
                      Peptide.2=ms3CSM2$Peptide.ms3,
                      Score=ms3CSM1$Score.ms3 + ms3CSM2$Score.ms3,
                      Score.Diff=NA,
                      Fraction=ms3CSM1$Fraction.ms3,
                      RT.1=ms3CSM1$RT.ms3,
                      RT.2=ms3CSM2$RT.ms3,
                      Spectrum.1=ms3CSM1$Spectrum.ms3,
                      Spectrum.2=ms3CSM2$Spectrum.ms3,
                      ppm=NA,
                      Elemental.Composition=NA,
                      Sc.1=ms3CSM1$Score.ms3,
                      Sc.2=ms3CSM2$Score.ms3,
                      Expect.1=ms3CSM1$Expect.ms3,
                      Expect.2=ms3CSM2$Expect.ms3,
                      Protein.1=ms3CSM1$`Protein.ms3`,
                      Protein.2=ms3CSM2$`Protein.ms3`,
                      Acc.1=ms3CSM1$`Acc.ms3`,
                      Acc.2=ms3CSM2$`Acc.ms3`,
                      # This needs to change to support decoys soon.
                      Decoy.1=ms3CSM1$Decoy.ms3,
                      Decoy.2=ms3CSM2$Decoy.ms3,
                      # Parsing of the Xlinked residues is still simplistic
                      # and does not account for site localization ambiguity
                      XLink.AA.1=str_match(ms3CSM1$`Protein.Mods.ms3`, 
                                           "Xlink:DSSO_[[as]]_fragment@([[0-9]]+)")[2],
                      XLink.AA.2=str_match(ms3CSM2$`Protein.Mods.ms3`, 
                                           "Xlink:DSSO_[[as]]_fragment@([[0-9]]+)")[2]
    )
    return(flatCSM)
}

processMS2xlinkResultsBoosted <- function(ms2searchResults, ms3searchTable, dev=F) {
    #ms2search results are the S04 object and they should not be processed beyond CSMs ideally.
    #ms3searchTable is the tsv d by running readMS3results
    
    datTab <- ms2searchResults@dataTable
    ms2rescue <- left_join(datTab, ms3searchTable, by=c("MSMS.Info"="ms2etdScanNo"))
    ms2rescue$outcome <- NA
    ms2rescue$outcome[ms2rescue$DB.Peptide.1 == ms2rescue$`DB.Peptide.ms3`] <- "agree Pep1"
    ms2rescue$outcome[ms2rescue$DB.Peptide.2 == ms2rescue$`DB.Peptide.ms3`] <- "agree Pep2"
    ms2rescue$outcome[!is.na(ms2rescue$`DB.Peptide.ms3`) &
                          ms2rescue$`DB.Peptide.ms3`!= ms2rescue$DB.Peptide.1 &
                          ms2rescue$`DB.Peptide.ms3`!= ms2rescue$DB.Peptide.2] <- "disagree"
    print(table(ms2rescue$outcome))
    if (dev) {return(ms2rescue)}
    datTab <- left_join(datTab,
                        ms2rescue[which(ms2rescue$outcome=="agree Pep2"),
                                  c("MSMS.Info","DB.Peptide.ms3", "Peptide.ms3", 
                                    "Protein.Mods.ms3", "Fraction.ms3", "RT.ms3", 
                                    "Spectrum.ms3", "ms3ScanNo", "Score.ms3", "Expect.ms3", 
                                    "ms2cidScanNo", "pepmass.ms3", "charge.ms3", 
                                    "precMass.ms3")],
                        by="MSMS.Info")
    datTab[is.na(datTab$ms3ScanNo), "Expect.ms3"] <- 1
    datTab <- datTab %>% mutate(Score.Diff.Boosted = Score.Diff + -log(Expect.ms3))
    ms2searchResults@dataTable <- datTab
    return(ms2searchResults)
}

# 
# peptideGroups <- list(
#     Group1=c(
#         "SDKNR",
#         "KLINGIR",
#         "KFDNLTK",
#         "FIKPILEK",
#         "APLSASMIKR",
#         "NPIDFLEAKGYK",
#         "LPKYSLFELENGR",
#         "TEVQTGGFSKESILPK"),
#     Group2=c(
#         "VKYVTEGMR",
#         "FDNLTKAER",
#         "DFQFYKVR",
#         "YDENDKLIR",
#         "MIAKSEQEIGK",
#         "HKPENIVIEMAR",
#         "TILDFLKSDGFANR",
#         "KIECFDSVEISGVEDR",
#         "YVNFLYLASHYEKLK"),
#     Group3=c(
#         "LSKSR",
#         "DKPIR",
#         "KDLIIK",
#         "MKNYWR",
#         "KGILQTVK",
#         "NSDKLIAR",
#         "DDSIDNKVLTR"),
#     Group4=c(
#         "KLVDSTDK",
#         "IEKILTFR",
#         "KAIVDLLFK",
#         "VLSAYNKHR",
#         "IEEGIKELGSQILK",
#         "SSFEKNPIDFLEAK",
#         "SNFDLAEDAKLQLSK",
#         "HSLLYEYFTVYNELTKVK"),
#     Group5=c(
#         "KVTVK",
#         "EKIEK",
#         "VITLKSK",
#         "QLKEDYFK",
#         "QLLNAKLITQR",
#         "GGLSELDKAGFIK",
#         "MDGTEELLVKLNR"),
#     Group6=c(
#         "EVKVITLK",
#         "KPAFLSGEQK",
#         "ENQTTQKGQK",
#         "KTEVQTGGFSK",
#         "VVDELVKVMGR",
#         "LESEFVYGDYKVYDVR",
#         "MLASAGELQKGNELALPSK",
#         "NFMQLIHDDSLTFKEDIQK",
#         "VLPKHSLLYEYFTVYNELTK"),
#     Group7=c(
#         "KMIAK",
#         "ESILPKR",
#         "DLIIKLPK",
#         "FKVLGNTDR",
#         "SEQEIGKATAK",
#         "AIVDLLFKTNR",
#         "LKTYAHLFDDK",
#         "VNTEITKAPLSASMIK",
#         "YDEHHQDLTLLKALVR"),
#     Group8=c(
#         "KDWDPK",
#         "QQLPEKYK",
#         "KVLSMPQVNIVK",
#         "MTNFDKNLPNEK",
#         "QITKHVAQILDSR",
#         "KSEETITPWNFEEVVDK",
#         "KNGLFGNLIALSLGLTPNFK",
#         "SKLVSDFR"),
#     Group9=c(
#         "LKSVK",
#         "IIKDK",
#         "DWDPKK",
#         "LKGSPEDNEQK",
#         "VLSMPQVNIVKK",
#         "LENLIAQLPGEKK",
#         "LIYLALAHMIKFR",
#         "YPKLESEFVYGDYK"),
#     Group10=c(
#         "VPSKK",
#         "VTVKQLK",
#         "EDYFKK",
#         "VKYVTEGMR",
#         "GKSDNVPSEEVVK",
#         "LEESFLVEEDKK",
#         "QEDFYPFLKDNR"),
#     Group11=c(
#         "GQKNSR",
#         "AGFIKR",
#         "GYKEVK",
#         "VMKQLK",
#         "KDFQFYK",
#         "LVDSTDKADLR",
#         "SDNVPSEEVVKK",
#         "KNLIGALLFDSGETAEATR"),
#     Group12=c(
#         "HSIKK",
#         "DKQSGK",
#         "NLPNEKVLPK",
#         "QSGKTILDFLK",
#         "MNTKYDENDK",
#         "SVKELLGITIMER",
#         "TYAHLFDDKVMK",
#         "FNASLGTYHDLLKIIK")
# )

# Max number of crosslinks:
#peptideGroups %>% map_dbl(function(x) length(x)**2) %>% sum

compareFDRs <- function(datTab, classifier="dvals", scalingFactor=10, ...) {
    class.max <- ceiling(max(datTab[[classifier]]))
    class.min <- floor(min(datTab[[classifier]]))
    if (classifier=="dvals") {
        range.spacing = 0.1
    } else if (classifier=="Score.Diff") {
        range.spacing = 0.25
    } else {
        range.spacing = abs(class.max - class.min) / 100
    }
    class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
    FDRs <- unlist(lapply(class.range, function(threshold) {
        calculateFDR(datTab, threshold=threshold, classifier=classifier, scalingFactor=scalingFactor)
    }))
    FDRs[is.na(FDRs)] <- 0
    GTs <- unlist(lapply(class.range, function(threshold) {
        calculateGT(datTab, threshold=threshold, classifier=classifier)
    }))
    GTs[is.na(GTs)] <- 0
    plot(GTs ~ FDRs, type="l", lwd=2, ...)
    abline(a = 0, b=1, lt=2, col="red")
}


fdrPlots <- function(datTab, scalingFactor = 10, cutoff = 0, classifier="dvals") {
    datTab <- as.data.frame(datTab)
    # if (classifier=="dvals") {
    #     stepSize = 0.5
    # } else if (classifier=="Score.Diff") {
    #     stepSize = 1
    # } else {
    #     stepSize = 0.5
    # }
    minValue <- floor(min(datTab[[classifier]], na.rm = T))
    maxValue <- ceiling(max(datTab[[classifier]], na.rm = T))
    stepSize = mfloor((maxValue - minValue) / 50, 0.25)
    decoyHist <- hist(datTab[datTab$Decoy=="Decoy", classifier], 
                      breaks=seq(minValue, maxValue, by=stepSize), 
                      plot=F)
    TD <- decoyHist$counts
    doubleDecoyHist <- hist(datTab[datTab$Decoy=="DoubleDecoy", classifier], 
                            breaks=seq(minValue, maxValue, by=stepSize), 
                            plot=F)
    DD <- doubleDecoyHist$counts
    doubleDecoyHist$counts = DD / (scalingFactor ** 2)
    decoyHist$counts = (TD /scalingFactor - 2 * DD/ (scalingFactor**2))
    targetHist <- hist(datTab[datTab$Decoy=="Target", classifier],
                       breaks=seq(minValue, maxValue, by=stepSize),
                       plot=F)
    diffHist <- targetHist
    diffHist$counts <- targetHist$counts - decoyHist$counts - doubleDecoyHist$counts
    displayLim =c(minValue,maxValue)
    plot(targetHist, col="lightblue", xlim=displayLim, 
         xlab="SVM Score", main="FDR plot")
    plot(decoyHist, add=T, col="salmon")
    plot(doubleDecoyHist, add=T, col="goldenrod1")
    abline(v=cutoff, lwd=4, lt=2, col="red")
    legend("topright", c("Target", "Decoy", "DoubleDecoy"), 
           fill = c("lightblue", "salmon", "goldenrod1"),
           bty="n")
}    

assignGroups <- function(datTab, pgroups=peptideGroups) {
    datTab$group.A <- map_chr(datTab$DB.Peptide.1, function(x) {
        groupMembership <- map_lgl(pgroups, function(group) {
            x %in% group
        })
        if (sum(groupMembership)==0) {
            return(NA)
        } else if (sum(groupMembership) > 1) {
            return(paste(which(groupMembership), collapse="."))
        } else {
            return(which(groupMembership))
        }
    })
    datTab$group.B <- map_chr(datTab$DB.Peptide.2, function(x) {
        groupMembership <- map_lgl(pgroups, function(group) {
            x %in% group
        })
        if (sum(groupMembership)==0) {
            return(NA)
        } else if (sum(groupMembership) > 1) {
            return(paste(which(groupMembership), collapse="."))
        } else {
            return(which(groupMembership))
        }
    })
    datTab[which(datTab$group.A == "2.10" & datTab$group.B == "2"), "group.A"] <- "2"
    datTab[which(datTab$group.B == "2.10" & datTab$group.A == "2"), "group.B"] <- "2"
    datTab$group.A[datTab$group.A=="2.10"] <- "10"
    datTab$group.B[datTab$group.B=="2.10"] <- "10"
    datTab <- datTab %>% mutate(groundTruth=group.A==group.B)
    datTab$groundTruth[which(is.na(datTab$groundTruth))] <- FALSE
    datTab$crosslink <- ifelse(datTab$DB.Peptide.1 < datTab$DB.Peptide.2,
                               paste(datTab$DB.Peptide.1, datTab$DB.Peptide.2, sep=":"),
                               paste(datTab$DB.Peptide.2, datTab$DB.Peptide.1, sep=":")
    )                                     
    return(datTab)
}

removeWeirdDeadEnds <- function(datTab) {
    datTab <- datTab %>% filter(!str_detect(Peptide.1, "Xlink:DSS[[23]]"),
                                !str_detect(Peptide.2, "Xlink:DSS[[23]]"))
    return(datTab)
}

senspe <- function(threshold, datTab, classifier="dvals") {
    datTab <- datTab[datTab$Decoy=="Target",]
    P <- sum(datTab$groundTruth)
    N <- sum(!datTab$groundTruth)
    datTab <- datTab[datTab[[classifier]] >= threshold, ]
    TP <- sum(datTab$groundTruth)
    FP <- sum(!datTab$groundTruth)
    FN <- P-TP
    TN <- N-FP
    sensitivity <- TP / P
    specificity <- TN / N
    fdr <- FP / P
    return(data.frame(thresh=threshold, sens=sensitivity, spec=specificity, fdr=fdr))
}

calculateGroundTruthFDR <- function(datTab) {
    results <- datTab %>% filter(Decoy=="Target") %>%
        group_by(Fraction, groundTruth) %>% 
        count() %>% 
        pivot_wider(names_from=groundTruth, values_from=n)
    names(results)[which(names(results)=="TRUE")] <- "Correct"
    names(results)[which(names(results)=="FALSE")] <- "Incorrect"
    results[is.na(results)] <- 0
    results <- results %>% 
        mutate(FDR = round((100 * Incorrect)/(Incorrect + Correct),2))
    means <- tibble(Fraction="means:", 
                    Incorrect=mean(results$Incorrect, na.rm=T),
                    Correct=mean(results$Correct, na.rm=T), 
                    FDR=mean(results$FDR, na.rm=T))
    sds <- tibble(Fraction="sds:", 
                  Incorrect=sd(results$Incorrect, na.rm=T),
                  Correct=sd(results$Correct, na.rm=T), 
                  FDR=sd(results$FDR, na.rm=T))
    results <- bind_rows(results, means, sds)
    return(results)
}

reduceToXlinks <- function(datTab) {
    datTab <- datTab %>% 
        group_by(Fraction) %>% 
        nest() %>% mutate(data=map(data, bestResPairDT)) %>%
        unnest(cols=c(data))
    return(datTab)
}

reduceToXlinkHacky <- function(datTab, classifier=classifier) {
    quoClass <- enquo(classifier)
    datTab <- datTab %>%
        group_by(Fraction) %>%
        nest() %>% mutate(data=map(data, bestResPairHackDT, !! quoClass)) %>%
        unnest(cols=c(data))
    return(datTab)
}

bestResPairHackDT <- function(datTab, classifier=dvals){
    quoClass <- enquo(classifier)
    datTabR <- datTab %>% group_by(xlinkedResPair) %>%
        filter(!! quoClass==max(!! quoClass)) %>%
        ungroup()
    return(datTabR)
}

getRandomCrosslinks <- function(pdbFile, numCrosslinks) {
    pdbFile <- pdbFile[pdbFile$resid=="LYS",]
    xl1 <- sample(nrow(pdbFile), numCrosslinks, replace=T)
    xl2 <- sample(nrow(pdbFile), numCrosslinks, replace=T)
    randomDists <- mapply(multiEuclideanDistance, list(pdbFile), 
                          pdbFile[xl1, "resno"], 
                          pdbFile[xl1, "chain"],
                          pdbFile[xl2, "resno"],
                          pdbFile[xl2, "chain"]
    )
    return(randomDists)
}

massErrorPlot <- function(massErrors, lowThresh, highThresh, lowPlotRange, highPlotRange) {
    hist(massErrors, col="slategray3", 
         xlim=c(lowPlotRange, highPlotRange), 
         breaks=seq(lowPlotRange, highPlotRange, 0.5),
         xlab = "Mass Error (ppm)",
         main = "Precursor Mass Deviation")
    abline(v = c(lowThresh, highThresh), lwd = 4, lt = 2, col = "red")
}

distancePlot <- function(targetDists, randomDists, threshold) {
    targetDists.size <- length(!is.na(targetDists))
    randomDists.size <- length(!is.na(randomDists))
    maxDist <- ceiling(max(c(randomDists, targetDists), na.rm=T)) + 5
    dists <- hist(targetDists, breaks=seq(-5, maxDist, 5), plot=F)
    r.dists <- hist(randomDists, breaks=seq(-5, maxDist, 5), plot=F)
    r.dists$counts <- r.dists$counts * (targetDists.size / randomDists.size)
    col1 <- rgb(141, 90, 151, alpha = 150, maxColorValue = 255)
    col2 <- rgb(184, 235, 208, alpha = 150, maxColorValue = 255)
    ylim.value <- ceiling(max(c(dists$counts, r.dists$counts)))
    plot(r.dists, col=col1, ylim=c(0, ylim.value),
         xlab = expression(paste("C", alpha, "-C", alpha, " Distance (Å)")),
         main = "Crosslink Violations")
    plot(dists, col=col2, add=T)
    abline(v=threshold, lwd=4, lt=2, col="orangered")
    legend("topright", c("experimental", "random distribution"), fill = c(col2, col1),
           bty="n")
}

numHitsPlot <- function(num.hits, threshold) {
    num.hits <- num.hits %>%
        pivot_longer(cols=c(-fdr,-thresh),
                     names_to="crosslinkClass",
                     values_to="numHits")
    p <- num.hits %>% ggplot(aes(x=fdr, y=numHits)) +
        geom_line(aes(col=crosslinkClass)) +
        theme_bw() +
        theme(legend.position="top") +
        scale_color_brewer(type="qual", palette="Set1") +
        xlab("FDR") +
        ylab("Num Hits")
    p <- p + geom_vline(xintercept = threshold,
                        col="red", size=1.5)
    print(p)
}

#findThreshold(test, targetER=0.01)
#numHitsPlot(ErrorTable(test), 0.01)





