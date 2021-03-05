print("touchStone module loaded")
proton <- 1.007276

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
        return(NA)
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
        if (grepl("r[0-9]\\_", proteinName)) {
            return ("Dec")
        } else if (grepl("decoy", proteinName)) {
            return ("Dec")
        } else if (is.null(chainMap[[proteinName]])) {
            return("")
        } else {chain <- chainMap[[proteinName]]
        return(chain)
        }
    }
    pdbList <- list(parsedPDB)
    Acc.1 <- as.character(searchTable$Acc.1)
    Acc.2 <- as.character(searchTable$Acc.2)
    chains1 <- vapply(Acc.1, chainLookup, FUN.VALUE="")
    chains2 <- vapply(Acc.2, chainLookup, FUN.VALUE="")
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

assignModules <- function(searchTable, moduleFile) {
    modTab <- read_tsv(moduleFile) %>% group_by(Subunit) %>% nest()
    Modul.1 <- searchTable %>% 
        mutate(XLink.AA.1 = ifelse(XLink.AA.1 == 0, 1, XLink.AA.1)) %>%
        left_join(modTab, by=c("Acc.1"="Subunit")) %>% 
        mutate(data = map2(data, XLink.AA.1, function(df, x) {
            if (is.null(df)) {
                return(NULL)
            } else {
                df$inRange = map2_lgl(df$Range_low, df$Range_high, between, x=x)
                if (sum(df$inRange, na.rm=T) >= 1) {
                    df <- filter(df, inRange) %>% slice(1)
                } else if (sum(is.na(df$inRange)) >= 1) {
                    df <- filter(df, is.na(inRange)) %>% slice(1)
                } else {
                    df <- df %>% slice(1)
                }
                return(df)
            }
        })) %>% 
        unnest(data, keep_empty=T) %>%
        pull(Module) %>%
        replace_na("unknown")
    Modul.2 <- searchTable %>% 
        mutate(XLink.AA.1 = ifelse(XLink.AA.2 == 0, 1, XLink.AA.2)) %>%
        left_join(modTab, by=c("Acc.2"="Subunit")) %>% 
        mutate(data = map2(data, XLink.AA.2, function(df, x) {
            if (is.null(df)) {
                return(NULL)
            } else {
                df$inRange = map2_lgl(df$Range_low, df$Range_high, between, x=x)
                if (sum(df$inRange, na.rm=T) >= 1) {
                    df <- filter(df, inRange) %>% slice(1)
                } else if (sum(is.na(df$inRange)) >= 1) {
                    df <- filter(df, is.na(inRange)) %>% slice(1)
                } else {
                    df <- df %>% slice(1)
                }
                return(df)
            }
        })) %>% 
        unnest(data, keep_empty=T) %>%
        pull(Module) %>%
        replace_na("unknown")
    searchTable$xlinkedModulPair <- ifelse(Modul.1 <= Modul.2,
                                           paste(Modul.1, Modul.2, sep="::"),
                                           paste(Modul.2, Modul.1, sep="::"))
    searchTable$xlinkedModulPair <- as.factor(searchTable$xlinkedModulPair)
    mods <- unique(c(Modul.1, Modul.2)) %>% str_sort()
    searchTable$Modul.1 <- factor(Modul.1, levels = mods)
    searchTable$Modul.2 <- factor(Modul.2, levels = mods)
    searchTable <- searchTable %>% add_count(xlinkedModulPair, name="numMPSM")
    return(searchTable)
}

calculateDecoys <- function(searchTable) {
    searchTable$Decoy <- "Target"
    decoyReg <- "(^r[[:digit:]]_|^dec|^DECOY)"
    searchTable[grepl(decoyReg, searchTable$Acc.1) |
                    grepl(decoyReg, searchTable$Acc.2),
                "Decoy"] <- "Decoy"
    searchTable[grepl(decoyReg, searchTable$Acc.1) &
                    grepl(decoyReg, searchTable$Acc.2),
                "Decoy"] <- "DoubleDecoy"
    searchTable$Decoy2 <- "Target"
    searchTable[grepl("Decoy", searchTable$Decoy), "Decoy2"] <- "Decoy"
    searchTable$Decoy2 <- factor(searchTable$Decoy2, levels=c("Decoy","Target"))
    searchTable$Decoy <- factor(searchTable$Decoy, levels=c("DoubleDecoy","Decoy","Target"))
    return(searchTable)
}

calculatePairs <- function(searchTable){
    searchTable$xlinkedProtPair <- ifelse(searchTable$Decoy2 == "Decoy",
                                          ifelse(searchTable$Protein.1 <= searchTable$Protein.2,
                                          paste("decoy", searchTable$Protein.1, searchTable$Protein.2, sep="::"),
                                          paste("decoy", searchTable$Protein.2, searchTable$Protein.1, sep="::")),
                                          ifelse(searchTable$Acc.1 <= searchTable$Acc.2,
                                          paste(searchTable$Acc.1, searchTable$Acc.2, sep="::"),
                                          paste(searchTable$Acc.2, searchTable$Acc.1, sep="::")))
    accs <- unique(c(searchTable$Acc.1, searchTable$Acc.2)) %>% str_sort()
    searchTable <- searchTable %>% mutate(Acc.1 = factor(Acc.1, levels = accs),
                                          Acc.2 = factor(Acc.2, levels = accs))
    searchTable$Res.1 <- paste(searchTable$XLink.AA.1, searchTable$Acc.1, sep=".")
    searchTable$Res.2 <- paste(searchTable$XLink.AA.2, searchTable$Acc.2, sep=".")
    searchTable$xlinkedResPair <- ifelse(searchTable$Res.1 <= searchTable$Res.2, 
                                         paste(searchTable$Res.1, searchTable$Res.2, sep="::"),
                                         paste(searchTable$Res.2, searchTable$Res.1, sep="::"))
    searchTable$xlinkedResPair <- as.factor(searchTable$xlinkedResPair)
    searchTable$xlinkedProtPair <- as.factor(searchTable$xlinkedProtPair)
    searchTable <- searchTable %>% add_count(xlinkedResPair, name="numCSM")
    searchTable <- searchTable %>% add_count(xlinkedProtPair, name="numPPSM")
    return(searchTable)
}

assignXLinkClass <- function(searchTable) {
    intraProteinLinks <- searchTable$Protein.1 == searchTable$Protein.2
    interProteinLinks <- searchTable$Protein.1 != searchTable$Protein.2
    if (sum(c("Start.1", "Start.2", "End.1", "End.2") %in% names(searchTable)==4)) {
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
    } else {
        searchTable$xlinkClass <- ifelse(intraProteinLinks, "intraProtein", "interProtein")
    }
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

calculateDecoyFractions <- function(datTab, scalingFactor=10) {
    fdrTable <- table(datTab$Decoy)
    if (is.na(fdrTable["DoubleDecoy"])) {fdrTable["DoubleDecoy"] <- 0}
    if (is.na(fdrTable["Decoy"])) {fdrTable["Decoy"] <- 0}
    if (is.na(fdrTable["Target"])) {fdrTable["Target"] <- 0}
    ffTT <- fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2)
    ftTT <- (fdrTable[["Decoy"]] / scalingFactor) - 
        (2 * fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2))
    TT <- fdrTable[["Target"]]
    return(c("TT"=TT, "ftTT"=ftTT, "ffTT"=ffTT))
}

calculateFDR <- function(datTab, threshold=-100, classifier="dvals", scalingFactor=10) {
    datTab <- datTab[datTab[[classifier]] >= threshold, ]
    decoyFractions <- calculateDecoyFractions(datTab, scalingFactor)
    fdr <- (decoyFractions["ffTT"] + decoyFractions["ftTT"]) / decoyFractions["TT"]
    names(fdr) <- NULL
    return(fdr)
}

calculateSeparateFDRs <- function(datTab, thresholdIntra = -100, thresholdInter = -100, 
                                  classifier="dvals", scalingFactor=10) {
    intraHits <- datTab[str_detect(datTab$xlinkClass, "intraProtein"),]
    interHits <- datTab[str_detect(datTab$xlinkClass, "interProtein"),]
    intraHits <- intraHits[intraHits[[classifier]] >= thresholdIntra, ]
    interHits <- interHits[interHits[[classifier]] >= thresholdInter, ]
    intraFractions <- calculateDecoyFractions(intraHits, scalingFactor)
    interFractions <- calculateDecoyFractions(interHits, scalingFactor)
    decoyFractions <- intraFractions + interFractions
    fdr <- (decoyFractions["ffTT"] + decoyFractions["ftTT"]) / decoyFractions["TT"]
    names(fdr) <- NULL
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
    # if (classifier=="dvals") {
    #     range.spacing = 0.1
    # } else if (classifier=="Score.Diff") {
    #     range.spacing = 0.25
    # } else {
    range.spacing = abs(class.max - class.min) / 200
    # }
    class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
    error.rates <- unlist(lapply(class.range, function(threshold) {
        errorFUN(datTab, threshold=threshold, classifier=classifier, ...)
    }))
    error.rates[is.na(error.rates)] <- 0
    #plot(error.rates ~ class.range, type="l")
    
    num.hits <- map_dfr(class.range, function(threshold) {
        calculateHits(removeDecoys(datTab), threshold=threshold, classifier=classifier)
    })
    num.hits$fdr <- error.rates
    num.hits$total <- num.hits$inter + num.hits$intra
    num.hits <- num.hits %>% filter(total > 0)
    return(num.hits)
}

generateErrorTableSeparate <- function(datTab, classifier="dvals",
                                       errorFUN=calculateFDR, ...) {
    intraHits <- generateErrorTable(filter(datTab, str_detect(xlinkClass, "intraProtein")),
                                    classifier=classifier, errorFUN=errorFUN, ...)
    interHits <- generateErrorTable(filter(datTab, str_detect(xlinkClass, "interProtein")),
                                    classifier=classifier, errorFUN=errorFUN, ...)
    max.fdr <- mmax(max(intraHits$fdr, interHits$fdr), 0.05)
    fdr.spacing <- max.fdr / 50
    fdr.seq <- seq(0, max.fdr, fdr.spacing)
    
    fdr.seq %>% map_dfr(function(i) {
        intra <- intraHits %>% slice(which.min(abs(fdr - i))) %>%
            pull(intra)
        inter <- interHits %>% slice(which.min(abs(fdr - i))) %>%
            pull(inter)
        tibble(thresh = NA, intra, inter, fdr = i, total = intra + inter)
    })
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
        threshold <- minThreshold
    } else if (sum(error.cross)==0 & error.diff[1] > 0) {
        threshold <- max(class.range)
        print("no acceptable threshold found")
    } else {
        firstCross <- which(error.cross)[1] + 1
        threshold <- class.range[firstCross]
    }
    if (threshold < minThreshold) {
        threshold <- minThreshold
    }
    return(list("globalThresh"=threshold, "correspondingFDR" = errorFUN(datTab, threshold, classifier, ...)))
}

findSeparateThresholds <- function(datTab, targetER=0.05, minThreshold=-5,
                                   classifier="dvals", errorFUN=calculateFDR, ...) {
    interTab <- datTab %>% 
        filter(str_detect(xlinkClass, "inter"))
    interThresh <- findThreshold(interTab, targetER, minThreshold, classifier, errorFUN, ...)[[1]]
    intraTab <- datTab %>% 
        filter(str_detect(xlinkClass, "intra"))
    intraThresh <- findThreshold(intraTab, targetER, minThreshold, classifier, errorFUN, ...)[[1]]
    return(list("intraThresh"=intraThresh, "interThresh"=interThresh))
}

calculateSeparateHits <- function(datTab, thresholds = c("intraThresh"=-100, "interThresh"=-100), classifier="dvals") {
    datTab.intra <- datTab[datTab[[classifier]] >= thresholds["intraThresh"], ]
    intra <- sum(str_count(datTab.intra$xlinkClass, "intraProtein"))
    datTab.inter <- datTab[datTab[[classifier]] >= thresholds["interThresh"], ]
    inter <- sum(str_count(datTab.inter$xlinkClass, "interProtein"))
    return(tibble("thresh"=thresholds["interThresh"], "intra"=intra, "inter"=inter))
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

# renumberProtein <-function(dt, protein, shift) {
#     dt[dt$Acc.1==protein,"XLink.AA.1"] <-
#         dt[dt$Acc.1==protein,"XLink.AA.1"] + shift
#     dt[dt$Acc.2==protein,"XLink.AA.2"] <-
#         dt[dt$Acc.2==protein,"XLink.AA.2"] + shift
#     dt <- assignModules(dt, modulFile)
#     dt <- calculatePairs(dt)
#     if (is.data.frame(object@caTracedPDB)) {
#         dt <- measureDistances(dt,object@caTracedPDB,object@chainMap)
#     }
#     return(dt)
#}

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

buildClassifierExperimental <- function(datTab, params, scoreName) {
    datTab$massError <- abs(datTab$ppm - mean(datTab$ppm))
    datTab$charge <- as.factor(datTab$z)
    if (nrow(datTab) > 30000) {
        sampleNo <- 15000
    } else {
        sampleNo <- nrow(datTab)/2
    }
    ind <- sample(nrow(datTab),sampleNo)
    train <- datTab[ind,]
    test <- datTab[-1 * ind,]
    #params <- c("percMatched", "Score.Diff", "charge", "massError", "numCSM","numPPSM", "xlinkClass")
    #params <- c("Score.Diff", "percMatched", "massError", "charge", "numPPSM", "xlinkClass", "numCSM", "pep2.longest.run")
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
    datTab[[scoreName]] <- as.numeric(attr(p,"decision.values"))
    print(paste("training weights:", round(wghts,2)))
    tab <- table(datTab$Decoy2, datTab[[scoreName]] > 0)
    if (tab[1] < tab[3]) {
        datTab[[scoreName]] <- -1 * datTab[[scoreName]]
        tab <- table(datTab$Decoy2, datTab[[scoreName]] > 0)
    }
    print(tab)
    print(paste("specificity:", round(tab[1]/(tab[1]+tab[3]),2)))
    return(datTab)
}

buildClassifierMS3 <- function(datTab) {
    datTab$massError <- abs(datTab$ppm)
    if (nrow(datTab) > 30000) {
        sampleNo <- 15000
    } else {
        sampleNo <- nrow(datTab)/2
    }
    ind <- sample(nrow(datTab),sampleNo)
    train <- datTab[ind,]
    test <- datTab[-1 * ind,]
    params <- c("Sc.2","z","Sc.1","numCSM","massError")
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
    # if (tab[1] < tab[3]) {
    #     datTab$dvals <- -1 * datTab$dvals
    #     tab <- table(datTab$Decoy2, datTab$dvals > 0)
    # }
    print(tab)
    print(paste("specificity:", round(tab[1]/(tab[1]+tab[3]),2)))
    return(datTab)
}

generateMSViewerLink <- function(path, fraction, z, peptide.1, peptide.2, spectrum, 
                                 instrumentType = NA, outputType = "HTML") {
    if(!str_detect(fraction, "\\.[[a-z]]+$")) {
        fraction <- paste0(fraction, ".mgf")
    }
    if (is.na(instrumentType)) {
        instrumentType <- switch(
            #        str_extract(fraction, "(?<=FTMSms2)[[a-zA-Z]]+"),
            str_extract(fraction, "[[a-z]]+(?=\\.[[a-z]]+$)"),
            "ESI-Q-high-res",
            ethcd = "ESI-EThcD-high-res",
            etd = "ESI-ETD-high-res",
            hcd = "ESI-Q-high-res" #Add other instrument types
        )
    }
    linkType <- str_extract(peptide.1, "(?<=\\(\\+).+?(?=\\)([[A-Z]]|$|-))")
    templateVals["output_type"] <- outputType
    templateVals["data_filename"] <- file.path(path, fraction)
    templateVals["instrument_name"] <- instrumentType
    templateVals["scan_number"] <- spectrum
    templateVals["max_charge"] <- z
    templateVals["msms_precursor_charge"] <- z
    templateVals["sequence"] <- peptide.1
    templateVals["sequence2"] <- peptide.2
    templateVals["link_search_type"] <- linkType
    templateNames <- names(templateVals)
    templateVals <- url_encode(templateVals)
    names(templateVals) <- templateNames
    zipped <- str_c(names(templateVals), templateVals, sep="=", collapse="&")
    if (outputType == "HTML") {
        return(str_c('<a href=\"', zipped, '\" target=\"_blank\">Spectrum</a>'))
    } else {
        return(zipped)
    }
}

generateMSViewerLink.ms3 <- function(path, fraction, rt, z, peptide, spectrum) {
    if(!str_detect(fraction, "\\.[[a-z]]$")) {
        fraction <- paste0(fraction, ".mgf")
    }
    instrumentType <- switch(
        #        str_extract(fraction, "(?<=FTMSms2)[[a-zA-Z]]+"),
        str_extract(fraction, "[[a-z]]+(?=\\.[[a-z]]+$)"),
        "ESI-Q-high-res",
        ethcd = "ESI-EThcD-high-res",
        etd = "ESI-ETD-high-res",
        hcd = "ESI-Q-high-res",
        cid = "ESI-ION-TRAP-low-res"
    )
    templateVals.ms3[6] <- file.path(path, fraction)
    templateVals.ms3[9] <- instrumentType
    templateVals.ms3[18] <- rt
    templateVals.ms3[20] <- spectrum
    templateVals.ms3[21] <- z
    templateVals.ms3[22] <- z
    templateVals.ms3[23] <- peptide
    templateVals.ms3 <- url_encode(templateVals.ms3)
    zipped <- str_c(templateKeys.ms3[1:25], templateVals.ms3[1:25], sep="=", collapse="&")
    str_c('<a href=\"', zipped, '\" target=\"_blank\">', peptide, '</a>')
}

generateCheckBoxes <- function(datTab) {
    datTab$selected <- str_c('<input type="checkbox" names="row', datTab$xlinkedResPair, '"value="', datTab$xlinkedResPair, '">',"")
    return(datTab)
}

readMSProductInfo <- function(datTab) {
    require(rvest)
    require(progress)
#    msvFiles <- "/var/lib/prospector/seqdb/web/results/msviewer/A/V/AV11111111/Z20201215hcd"
    msvFiles <- "/var/lib/prospector/seqdb/web/results/msviewer/A/V/AV22222222/Z20201215ethcd"
    pb <- progress_bar$new(
        format = " reading spectral match [:bar] :percent eta: :eta",
        total = nrow(datTab)
        )
    ms.product.info <-
        pmap_chr(list(msvFiles, datTab$Fraction, datTab$RT, datTab$z, 
                      datTab$Peptide.1, datTab$Peptide.2, datTab$Spectrum, 
                      "Tab delimited text"), generateMSViewerLink) %>%
        map(function(msvLink) {
            spec.html <- read_html(msvLink)
            spec.node <- html_node(spec.html, xpath = '//*[@id="centerbody"]')
            spec.table <- read_tsv(html_text(spec.node))
            pb$tick()
            return(spec.table)
        })
    return(ms.product.info)
}

getPercentMatched <- function(ms.product.info) {
    intensities <- map_dfr(ms.product.info, function(x) {
        matchedIntensity <- sum((!is.na(x$`Peptide #`)) * x$Intensity, na.rm=T)
        totalIntensity <- sum(x$Intensity, na.rm=T)
        percentMatched <- matchedIntensity / totalIntensity
        sums <- x %>% group_by(`Peptide #`) %>%
            summarize(percIntensity = sum(Intensity, na.rm=T) / totalIntensity, .groups="drop") %>%
            pivot_wider(names_from = `Peptide #`, values_from = percIntensity, names_prefix = "percInt.pep")
        sums$percMatched <- percentMatched
        sums$percDiag.pep1 <- x %>% 
            filter(str_detect(`Ion Type`, "^MH[[\\*\\#]]"), `Peptide #` == "1") %>% 
            pull(Intensity) %>% 
            sum(na.rm=T) / totalIntensity
        sums$percDiag.pep2 <- x %>% 
            filter(str_detect(`Ion Type`, "^MH[[\\*\\#]]"), `Peptide #` == "2") %>% 
            pull(Intensity) %>% 
            sum(na.rm=T) / totalIntensity
        return(sums)
    })
    replaceVals <- rep(0.0, ncol(intensities))
    names(replaceVals) <- names(intensities)
    intensities <- replace_na(intensities, as.list(replaceVals))
    return(intensities)
}

longestBackboneIonSeries <- function(ms.product.info) {
    map_dfr(ms.product.info, function(spectrumTable) {
        spectrumTable <- spectrumTable %>% 
            mutate(`Ion Type` = str_replace(`Ion Type`, "^([[yb]]).+", "\\1")) %>%
            group_by(`Peptide #`, `Ion Type`) %>% 
            arrange(Index) %>% 
            nest()
        pep1.longest <- spectrumTable %>%
            filter(`Peptide #` == 1) %>%
            mutate(pep1.longest.run = map_dbl(data,
                                              function(x) {
                                                  ions <- unique(x$Index)
                                                  runLength <- max(rle(ions - lag(ions))$lengths, na.rm=T)
                                                 # runLength[is.infinite(runLength)] <- 0
                                              })
            ) %>%
            pull(pep1.longest.run) %>% max(na.rm=T)
        pep1.longest[is.infinite(pep1.longest)] <- 0
        pep2.longest <- spectrumTable %>%
            filter(`Peptide #` == 2) %>%
            mutate(pep2.longest.run = map_dbl(data,
                                              function(x) {
                                                  ions <- unique(x$Index)
                                                  runLength <- max(rle(ions - lag(ions))$lengths, na.rm=T)
                                                  #runLength[is.infinite(runLength)] <- 0
                                              })
            ) %>%
            pull(pep2.longest.run) %>% max(na.rm=T)
        pep2.longest[is.infinite(pep2.longest)] <- 0
        return(tibble(pep1.longest.run = pep1.longest, 
                      pep2.longest.run = pep2.longest))
    })
}

formatXLTable <- function(datTab) {
    annoyingColumns <- str_which(names(datTab), "(Int|Dec)[a-z]{2}\\.[[1-2]]")
    datTab <- datTab[, -annoyingColumns]
    datTab <- datTab %>%
        select(-starts_with("Res"),
               -starts_with("Decoy"),
               -starts_with("Num\\."),
               -any_of("massError"))
    if (sum(!is.na(datTab$distance)) == 0) {
        datTab <- datTab %>%
            select(-any_of("distance"))
    }
    datTab <- datTab[order(datTab$dvals, decreasing = T),]
    datTab <- datTab %>% select(any_of(c("selected", "link", "link.1", "link.2", "xlinkedResPair",
                                         "xlinkedProtPair", "xlinkedModulPair",
                                "dvals", "SVM.new", "distance", "m.z", "z", "ppm",
                                "DB.Peptide.1", "DB.Peptide.2", "Score", "Score.Diff",
                                "Sc.1", "Rk.1", "Sc.2", "Rk.2", "Acc.1",
                                "XLink.AA.1", "Protein.1", "Modul.1", "Species.1",
                                "Acc.2", "XLink.AA.2", "Protein.2", "Modul.2", "Species.2",
                                "xlinkClass", "Len.Pep.1", "Len.Pep.2",
                                "Peptide.1", "Peptide.2", "numCSM", "numPPSM", "numMPSM",
                                "Fraction", "RT", "MSMS.Info")))
    if ("Protein.1" %in% names(datTab) & "Protein.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Protein.1 = factor(Protein.1), 
                                    Protein.2 = factor(Protein.2))
    }
    if ("Modul.1" %in% names(datTab) & "Modul.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Modul.1 = factor(Modul.1), 
                                    Modul.2 = factor(Modul.2))
    }
    if ("Species.1" %in% names(datTab) & "Species.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Species.1 = factor(Species.1), 
                                    Species.2 = factor(Species.2))
    }
    if ("Acc.1" %in% names(datTab) & "Acc.2" %in% names(datTab)) {
        datTab <- datTab %>% mutate(Acc.1 = factor(Acc.1), 
                                    Acc.2 = factor(Acc.2))
    }
    if ("xlinkClass" %in% names(datTab)) {
        datTab <- datTab %>% mutate(xlinkClass = factor(xlinkClass))
    }
    if ("Fraction" %in% names(datTab)) {
        datTab$Fraction <- gsub("(.*)\\.[^.]+$", "\\1", datTab$Fraction)
        datTab <- datTab %>% mutate(Fraction = factor(Fraction))
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
    precMass <- pepMass * charge - proton * charge
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

createMasterScanFile <- function(ms3PAVApeaklist, 
                                 ms2PAVApeaklistCID,
                                 ms2PAVApeaklistETD=NA) {

    ms3masterscansCID <- getIndicesForFile(ms3PAVApeaklist)
    names(ms3masterscansCID) <- c("ms3ScanNo", "ms2cidScanNo", 
                                  "pepmass.ms3", "charge.ms3")
    ms2masterscansCID <- getIndicesForFile(ms2PAVApeaklistCID)
    if (!is.na(ms2PAVApeaklistETD)) {
        ms2masterscansETD <- getIndicesForFile(ms2PAVApeaklistETD)
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

readMS3results <- function(ms3SearchTable){
    ms3results <- read_tsv(ms3SearchTable, skip=2) %>% 
        select(c("DB Peptide", "Peptide", "Protein Mods",
               "Fraction", "RT", "Spectrum", "MSMS Info", "Score",
               "Expect", "Protein Name", "Acc #", "Species"))
    names(ms3results) <- c("DB.Peptide.ms3", "Peptide.ms3", "Protein.Mods.ms3",
                           "Fraction.ms3", "RT.ms3", "Spectrum.ms3", "ms3ScanNo", 
                           "Score.ms3", "Expect.ms3", "Protein.ms3", "Acc.ms3",
                           "Species.ms3")
    ms3results$Decoy.ms3 <- ifelse(str_detect(ms3results$Acc.ms3, "decoy"),
                                   "DECOY", "Target")
    ms3results <- ms3results %>% 
#        filter(str_detect(Protein.Mods.ms3, "Xlink:DSSO_[[as]]_fragment"))
        filter(str_detect(Protein.Mods.ms3, "XL:A-[[a-zA-Z\\(\\)]]+"))
    return(ms3results)
}

addMasterScanInfo <- function(ms3results, masterScanFile) {
   ms3results <- left_join(ms3results, masterScanFile,
                           by = "ms3ScanNo")
   ms3results$precMass.ms3 <- map2_dbl(ms3results$pepmass.ms3,
                                       ms3results$charge.ms3, calculatePrecMass)
   return(ms3results)
}

processMS3xlinkResults <- function(ms3results) {
    n_ms3 <- ms3results %>% group_by(ms2cidScanNo) %>% nest()
    n_ms3 <- n_ms3 %>% mutate(scanGroupSize = map_int(data, ~nrow(.)))
    ms3crosslinksFlat <- map2_dfr(n_ms3$data, n_ms3$ms2cidScanNo, ms3ionFragFind, tol=50)
    ms3crosslinksFlatMis <- map2_dfr(n_ms3$data, n_ms3$ms2cidScanNo, ms3ionFragFind, tol=50, mis=proton)
    ms3crosslinksFlatMis2 <- map2_dfr(n_ms3$data, n_ms3$ms2cidScanNo, ms3ionFragFind, tol=50, mis=2*proton)
    ms3crosslinksFlatMisM1 <- map2_dfr(n_ms3$data, n_ms3$ms2cidScanNo, ms3ionFragFind, tol=50, mis=-1*proton)
    ms3crosslinksFlat <- rbind(ms3crosslinksFlatMisM1, ms3crosslinksFlatMis2, ms3crosslinksFlatMis, ms3crosslinksFlat)
    ms3crosslinksFlat$Score.Diff <- ms3crosslinksFlat$Sc.2
    ms3crosslinksFlat <- assignXLinkClass(ms3crosslinksFlat)
    ms3crosslinksFlat <- lengthFilter(ms3crosslinksFlat, minLen = 3, maxLen = 35)
    ms3crosslinksFlat <- scoreFilter(ms3crosslinksFlat, minScore = 0)
    ms3crosslinksFlat <- calculateDecoys(ms3crosslinksFlat)
    ms3crosslinksFlat <- calculatePairs(ms3crosslinksFlat)
    if (!"distance" %in% names(ms3crosslinksFlat)) {
        ms3crosslinksFlat$distance <- NA_real_
    }
    return(ms3crosslinksFlat)
}

getMassError <- function(targetMass, testMass) {
    # reports mass error in ppm
    massDiff <- 1e6 * (targetMass - testMass) / targetMass
    return(massDiff)
}

ms3ionFragFind <- function(scanGroupTibble, masterScanNo, tol = 25, mis = 0) {
    H2O <- 18.01002
    ms2Mass <- scanGroupTibble$precMass.ms2[1] - mis
    scanGroupSize <- nrow(scanGroupTibble)
    if (scanGroupSize < 2) {
        return(data.frame())
    }
    scanGroupTibble$alkene <- str_count(scanGroupTibble$Peptide.ms3, fixed("XL:A-Alkene"))
    scanGroupTibble$thiol <- str_count(scanGroupTibble$Peptide.ms3, fixed("XL:A-Thiol(Unsaturated)"))
    scanGroupTibble$sulfenic <- str_count(scanGroupTibble$Peptide.ms3, fixed("XL:A-Sulfenic"))
    
#    scanGroupTibble$alkene <- str_detect(scanGroupTibble$Protein.Mods.ms3, "Alkene@[[0-9\\|]]+")
#    scanGroupTibble$thiol <- str_detect(scanGroupTibble$Protein.Mods.ms3, "Thiol\\(Unsaturated\\)@[[0-9\\|]]+")
#    scanGroupTibble$sulfenic <- str_detect(scanGroupTibble$Protein.Mods.ms3, "Sulfenic@[[0-9\\|]]+")
    summedPairs <- scanGroupTibble$precMass.ms3 %o% rep(1, scanGroupSize) +
        rep(1, scanGroupSize) %o% scanGroupTibble$precMass.ms3
    massErrors.AT <- getMassError(ms2Mass, summedPairs + H2O)
    matchesMS2 <- abs(massErrors.AT) <= tol
    allowed <- scanGroupTibble$alkene %o% scanGroupTibble$thiol
    ATpairs <- which(matchesMS2 * allowed == 1)
    massErrors.AS <- getMassError(ms2Mass, summedPairs)
    matchesMS2 <- abs(massErrors.AS) <= tol
    allowed <- scanGroupTibble$alkene %o% scanGroupTibble$sulfenic
    ASpairs <- which(matchesMS2 * allowed == 1)
    numMatches <- length(c(ATpairs, ASpairs))
    dimensions <- dim(allowed)
    if (numMatches==0) return(data.frame())
    placeHolder <- vector("list", numMatches)
    index = 1
    for (CSM in ATpairs) {
        massError <- massErrors.AT[CSM]
        peptidePair <- arrayInd(CSM, dimensions)
        flatMS3pair <- flattenMS3pair(scanGroupTibble[peptidePair[,1],],
                                      scanGroupTibble[peptidePair[,2],],
                                      masterScanNo,
                                      frag1="Alkene", frag2="Thiol\\(Unsaturated\\)")
        flatMS3pair$ppm <- massError
        missing <- is.na(flatMS3pair$XLink.AA.1) | is.na(flatMS3pair$XLink.AA.2)
        flatMS3pair <- flatMS3pair[!missing, ]
        if (nrow(flatMS3pair)==0) {
            flatMS3pair <- data.frame()
        }
        placeHolder[[index]] <- flatMS3pair
        index = index + 1
    }
    for (CSM in ASpairs) {
        massError <- massErrors.AS[CSM]
        peptidePair <- arrayInd(CSM, dimensions)
        flatMS3pair <- flattenMS3pair(scanGroupTibble[peptidePair[,1],],
                                      scanGroupTibble[peptidePair[,2],],
                                      masterScanNo,
                                      frag1="Alkene", frag2="Sulfenic")
        flatMS3pair$ppm <- massError
        missing <- is.na(flatMS3pair$XLink.AA.1) | is.na(flatMS3pair$XLink.AA.2)
        flatMS3pair <- flatMS3pair[!missing, ]
        if (nrow(flatMS3pair)==0) {
            flatMS3pair <- data.frame()
        }
        placeHolder[[index]] <- flatMS3pair
        index = index + 1
    }
    return(do.call(rbind, placeHolder))
}

flattenMS3pair <- function(ms3CSM1, ms3CSM2, masterScanNo=NA,
                           frag1="Alkene", frag2="Thiol\\(Unsaturated\\)") {
    if (ms3CSM1$Expect.ms3 > ms3CSM2$Expect.ms3) {
        tempCSM <- ms3CSM1
        ms3CSM1 <- ms3CSM2
        ms3CSM2 <- tempCSM
        tempFrag <- frag1
        frag1 <- frag2
        frag2 <- tempFrag
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
                      Decoy.1=ms3CSM1$Decoy.ms3,
                      Decoy.2=ms3CSM2$Decoy.ms3,
                      z=ms3CSM1$charge.ms2,
                      z.1=ms3CSM1$charge.ms3,
                      z.2=ms3CSM2$charge.ms3,
                      AA.1=str_match_all(ms3CSM1$`Protein.Mods.ms3`, 
                                          str_c("XL:A-", frag1, "@([[0-9|]]+)"))[[1]][,2],
                      AA.2=str_match_all(ms3CSM2$`Protein.Mods.ms3`, 
                                          str_c("XL:A-", frag2, "@([[0-9|]]+)"))[[1]][,2],
                      XLink.AA.1=str_c(AA.1, sep="|"),
                      XLink.AA.2=str_c(AA.2, sep="|")
                      )
    return(flatCSM)
}

processMS3xlinkResultsSingleFile <- function(scResult, ms3file, ms2file) {
    masterScanList <- createMasterScanFile(ms3file, ms2file)
    ms3 <- addMasterScanInfo(scResult, masterScanList)
    datTab <- processMS3xlinkResults(ms3)
    return(datTab)
}

processMS3xlinkResultsMultiFile <- function(scResults, ms3files, ms2files) {
    #scResults is search compare output for MS3 search of all files
    #ms3files and ms2files are character vectors with PAVA generated file names
    searchCompareMS3 <- readMS3results(scResults)
    searchCompareMS3 <- searchCompareMS3 %>% group_by(Fraction.ms3) %>% nest()
    baseFiles <- str_replace(searchCompareMS3$Fraction.ms3, "_[[A-Z]]+MSms[[0-9]][[a-z]]+", "")
    numMS3Files <- length(ms3files)
    numMS2Files <- length(ms2files)
    if ((numMS3Files != numMS2Files) | (numMS3Files != nrow(searchCompareMS3))) {
        print("incorrect input")
        return()
    }
    datTab <- map_dfr(baseFiles, function(bf) {
        print(bf)
        p1 <- searchCompareMS3 %>% filter(str_detect(Fraction.ms3, bf)) %>% unnest(cols=c(data))
        p2 <- ms3files[str_which(ms3files, bf)]
        p3 <- ms2files[str_which(ms2files, bf)]
        processMS3xlinkResultsSingleFile(p1, p2, p3)
    })
    return(datTab)
}


# processMS2xlinkResultsBoosted <- function(ms2searchResults, ms3searchTable, dev=F) {
#     #ms2search results are the S04 object and they should not be processed beyond CSMs ideally.
#     #ms3searchTable is the tsv d by running readMS3results
#     
#     datTab <- ms2searchResults@dataTable
#     ms2rescue <- left_join(datTab, ms3searchTable, by=c("MSMS.Info"="ms2etdScanNo"))
#     ms2rescue$outcome <- NA
#     ms2rescue$outcome[ms2rescue$DB.Peptide.1 == ms2rescue$`DB.Peptide.ms3`] <- "agree Pep1"
#     ms2rescue$outcome[ms2rescue$DB.Peptide.2 == ms2rescue$`DB.Peptide.ms3`] <- "agree Pep2"
#     ms2rescue$outcome[!is.na(ms2rescue$`DB.Peptide.ms3`) &
#                           ms2rescue$`DB.Peptide.ms3`!= ms2rescue$DB.Peptide.1 &
#                           ms2rescue$`DB.Peptide.ms3`!= ms2rescue$DB.Peptide.2] <- "disagree"
#     print(table(ms2rescue$outcome))
#     if (dev) {return(ms2rescue)}
#     datTab <- left_join(datTab,
#                         ms2rescue[which(ms2rescue$outcome=="agree Pep2"),
#                                   c("MSMS.Info","DB.Peptide.ms3", "Peptide.ms3", 
#                                     "Protein.Mods.ms3", "Fraction.ms3", "RT.ms3", 
#                                     "Spectrum.ms3", "ms3ScanNo", "Score.ms3", "Expect.ms3", 
#                                     "ms2cidScanNo", "pepmass.ms3", "charge.ms3", 
#                                     "precMass.ms3")],
#                         by="MSMS.Info")
#     datTab[is.na(datTab$ms3ScanNo), "Expect.ms3"] <- 1
#     datTab <- datTab %>% mutate(Score.Diff.Boosted = Score.Diff + -log(Expect.ms3))
#     ms2searchResults@dataTable <- datTab
#     return(ms2searchResults)
# }
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

# compareFDRs <- function(datTab, classifier="dvals", scalingFactor=10, ...) {
#     class.max <- ceiling(max(datTab[[classifier]]))
#     class.min <- floor(min(datTab[[classifier]]))
#     if (classifier=="dvals") {
#         range.spacing = 0.1
#     } else if (classifier=="Score.Diff") {
#         range.spacing = 0.25
#     } else {
#         range.spacing = abs(class.max - class.min) / 100
#     }
#     class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
#     FDRs <- unlist(lapply(class.range, function(threshold) {
#         calculateFDR(datTab, threshold=threshold, classifier=classifier, scalingFactor=scalingFactor)
#     }))
#     FDRs[is.na(FDRs)] <- 0
#     GTs <- unlist(lapply(class.range, function(threshold) {
#         calculateGT(datTab, threshold=threshold, classifier=classifier)
#     }))
#     GTs[is.na(GTs)] <- 0
#     plot(GTs ~ FDRs, type="l", lwd=2, ...)
#     abline(a = 0, b=1, lt=2, col="red")
# }

fdrPlots <- function(datTab, scalingFactor = 10, cutoff = 0, classifier="dvals",
                     minValue = floor(min(datTab[[classifier]], na.rm = T)),
                     maxValue = ceiling(max(datTab[[classifier]], na.rm = T)),
                     addLegend = T,
                     title = "FDR plot",
                     xlabel = "SVM Score") {
    datTab <- as.data.frame(datTab)
    stepSize = mmax((maxValue - minValue) / 50, 0.25)
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
         xlab=xlabel, main=NA)
    plot(decoyHist, add=T, col="salmon")
    plot(doubleDecoyHist, add=T, col="goldenrod1")
    abline(v=cutoff, lwd=2, lt=1, col="red")
    if (addLegend) {
        legend("topright", c("Target", "Decoy", "DoubleDecoy"), 
               fill = c("lightblue", "salmon", "goldenrod1"),
               bty="n")
    }
    title(title, adj=0)
}    

# assignGroups <- function(datTab, pgroups=peptideGroups) {
#     datTab$group.A <- map_chr(datTab$DB.Peptide.1, function(x) {
#         groupMembership <- map_lgl(pgroups, function(group) {
#             x %in% group
#         })
#         if (sum(groupMembership)==0) {
#             return(NA)
#         } else if (sum(groupMembership) > 1) {
#             return(paste(which(groupMembership), collapse="."))
#         } else {
#             return(which(groupMembership))
#         }
#     })
#     datTab$group.B <- map_chr(datTab$DB.Peptide.2, function(x) {
#         groupMembership <- map_lgl(pgroups, function(group) {
#             x %in% group
#         })
#         if (sum(groupMembership)==0) {
#             return(NA)
#         } else if (sum(groupMembership) > 1) {
#             return(paste(which(groupMembership), collapse="."))
#         } else {
#             return(which(groupMembership))
#         }
#     })
#     datTab[which(datTab$group.A == "2.10" & datTab$group.B == "2"), "group.A"] <- "2"
#     datTab[which(datTab$group.B == "2.10" & datTab$group.A == "2"), "group.B"] <- "2"
#     datTab$group.A[datTab$group.A=="2.10"] <- "10"
#     datTab$group.B[datTab$group.B=="2.10"] <- "10"
#     datTab <- datTab %>% mutate(groundTruth=group.A==group.B)
#     datTab$groundTruth[which(is.na(datTab$groundTruth))] <- FALSE
#     datTab$crosslink <- ifelse(datTab$DB.Peptide.1 < datTab$DB.Peptide.2,
#                                paste(datTab$DB.Peptide.1, datTab$DB.Peptide.2, sep=":"),
#                                paste(datTab$DB.Peptide.2, datTab$DB.Peptide.1, sep=":")
#     )                                     
#     return(datTab)
# }

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

bestResPair <- function(datTab, classifier=dvals){
    quoClass <- enquo(classifier)
    datTab <- datTab %>% group_by(xlinkedResPair) %>%
        filter(!! quoClass==max(!! quoClass)) %>%
        ungroup()
    return(datTab)
}

bestProtPair <- function(datTab, classifier=dvals){
    quoClass <- enquo(classifier)
    datTab <- datTab %>% group_by(xlinkedProtPair) %>%
        filter(!! quoClass==max(!! quoClass)) %>%
        ungroup()
    return(datTab)
}

bestModPair <- function(datTab, classifier=dvals){
    quoClass <- enquo(classifier)
    datTab <- datTab %>% group_by(xlinkedModulPair) %>%
        filter(!! quoClass==max(!! quoClass)) %>%
        ungroup()
    return(datTab)
}

bestResPairFraction <- function(datTab, classifier=classifier) {
    quoClass <- enquo(classifier)
    datTab <- datTab %>%
        group_by(Fraction) %>%
        nest() %>% mutate(data=map(data, bestResPair, !! quoClass)) %>%
        unnest(cols=c(data))
    return(datTab)
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
    if (abs(lowPlotRange) < 10 & abs(highPlotRange) < 10) {
        binw <- 0.5
    } else {
        binw <- 1
    }
    hist(massErrors, col="slategray3", 
         xlim=c(lowPlotRange, highPlotRange), 
         breaks=seq(lowPlotRange, highPlotRange, binw),
         xlab = "Mass Error (ppm)",
         main = NA)
    title("Precursor Mass Deviation", adj=0)
    abline(v = c(lowThresh, highThresh), lwd = 2, lt = 1, col = "red")
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
         main = NA)
    title("Crosslink Violations", adj=0)
    plot(dists, col=col2, add=T)
    abline(v=threshold, lwd=2, lt=1, col="red")
    legend("topright", c("experimental", "random distribution"), fill = c(col2, col1),
           bty="n")
}

numHitsPlot <- function(num.hits, threshold=0) {
    num.hits <- num.hits %>%
        filter(!is.infinite(fdr), total >= 0.025 * max(total, na.rm=T)) %>%
        pivot_longer(cols=c(-fdr,-thresh),
                     names_to="crosslinkClass",
                     values_to="numHits") %>%
        filter(!is.na(numHits))
    p <- num.hits %>% ggplot(aes(x=fdr, y=numHits)) +
        geom_path(aes(col=crosslinkClass), na.rm=T, size=1.25) +
        theme_bw() +
        theme(legend.position="right") +
        scale_color_brewer(type="qual", palette="Paired") +
        xlab("FDR") +
        ylab("Num Hits") +
        ggtitle("Num Hits vs FDR") +
        theme(plot.title = element_text(face = "bold", size = 14),
              )
    p <- p + geom_vline(xintercept = threshold,
                        col="red", size=1.25)
    print(p)
}

bestResPair <- function(datTab, classifier=dvals){
    quoClass <- enquo(classifier)
    datTab <- datTab %>% group_by(xlinkedResPair) %>%
        filter(!! quoClass==max(!! quoClass)) %>%
        ungroup()
    return(datTab)
}

summarizeProtData <- function(datTab) {
    datTab <- datTab %>% filter(Decoy=="Target")
    datTab <- bestResPair(datTab)
    datTab$Acc.1 <- fct_drop(datTab$Acc.1, only="decoy")
    datTab$Acc.2 <- fct_drop(datTab$Acc.2, only="decoy")
    matrixCounts <- datTab %>%
        group_by(Acc.1, Acc.2) %>%
        summarize(accCounts = sum(numCSM), .groups = "drop") %>%
        complete(Acc.1, Acc.2) %>% 
        unique() %>%
        pivot_wider(names_from = Acc.2, values_from = accCounts) %>%
        column_to_rownames("Acc.1") %>%
        as.matrix()
    matrixCounts[is.na(matrixCounts)] <- 0
    matrixCounts <- matrixCounts + t(matrixCounts)
    diag(matrixCounts) <- diag(matrixCounts) / 2
    accNames <- rownames(matrixCounts)
    matrixCounts <- clearAboveDiag(matrixCounts)
    matrixCounts <- as_tibble(matrixCounts)
    matrixCounts$Acc.1 <- accNames
    matrixCounts <- matrixCounts %>%
        pivot_longer(cols = -Acc.1, names_to = "Acc.2", values_to = "counts")
    return(matrixCounts)
}

summarizeModuleData <- function(datTab) {
    datTab <- datTab %>% filter(Decoy=="Target")
    datTab <- bestResPair(datTab)
    # datTab$Acc.1 <- fct_drop(datTab$Acc.1, only="decoy")
    # datTab$Acc.2 <- fct_drop(datTab$Acc.2, only="decoy")
    matrixCounts <- datTab %>%
        group_by(Modul.1, Modul.2) %>%
        summarize(modCounts = sum(numCSM), .groups = "drop") %>%
        complete(Modul.1, Modul.2) %>% 
        unique() %>%
        pivot_wider(names_from = Modul.2, values_from = modCounts) %>%
        column_to_rownames("Modul.1") %>%
        as.matrix()
    matrixCounts[is.na(matrixCounts)] <- 0
    matrixCounts <- matrixCounts + t(matrixCounts)
    diag(matrixCounts) <- diag(matrixCounts) / 2
    modNames <- rownames(matrixCounts)
    matrixCounts <- clearAboveDiag(matrixCounts)
    matrixCounts <- as_tibble(matrixCounts)
    matrixCounts$Modul.1 <- modNames
    matrixCounts <- matrixCounts %>%
        pivot_longer(cols = -Modul.1, names_to = "Modul.2", values_to = "counts")
    return(matrixCounts)
}

summaryPlot <- function(summarizedData) {
    accPlot <- summarizedData %>% 
        ggplot(aes(Acc.1, Acc.2, size=counts, col=counts)) + 
        geom_point(na.rm=T) +
        scale_size(range = c(-1, 12)) +
        scale_color_viridis_c(option="D") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle=90, hjust=1))
    print(accPlot)
}

clearAboveDiag <- function(sqMatrix) {
    n <- dim(sqMatrix)[1]
    m <- dim(sqMatrix)[2]
    if (n != m) return("only works for square matrices")
    a <- 1 + (1:n-1)*n
    b <- (1:n-1)*(n+1)
    abovePositions <- unlist(map2(a, b, seq)[2:n])
    sqMatrix[abovePositions] <- NA
    return(sqMatrix)
}

makeXiNetFile <- function(datTab) {
    datTab <- datTab %>% select(dvals, Acc.1, Acc.2, XLink.AA.1, XLink.AA.2)
    names(datTab) <- c("Score", "Protein1", "Protein2", "LinkPos1", "LinkPos2")
    return(datTab)
}