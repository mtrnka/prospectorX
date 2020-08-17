library(jsonlite)
# Prospector JSON contain -nan for coverage.  JSON format desired null instead.
system2("sed", args = c("'s/-nan/null/g'", "testCache.json", ">", "testCache2.json"))
json <- fromJSON("testCache2.json", simplifyVector = F)

samples <- json$samples
id <- samples$id
#proteins <- id$proteins
crosslinks <- id$crosslinks
xl <- crosslinks$`search$1`
xlp <- xl$crosslinked_proteins
xlp1 <- xlp[[1]]
xlp1pep <- xlp1$peptides

map_df(xlp1pep, as_data_frame)


test <- as.list(letters[1:5])
test[[3]] <- as.list(LETTERS[1:3])
as.character(purrr::flatten(test))
