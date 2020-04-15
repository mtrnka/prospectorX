fetchInstanceList <- function() {
   instances <- gce_list_instances()$items
   machinesType <- instances$machineType
   machineType <- str_extract(machinesType, "(?<=\\/)[[a-z]][[0-9]]+\\-.+?$")
   status <- instances$status
   instanceTable <- data.frame(names = instances$name,
                               types = machineType,
                               status = status,
                               stringsAsFactors = F)
   return(instanceTable)
}


