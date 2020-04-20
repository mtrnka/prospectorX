fetchInstanceList <- function() {
   instances <- gce_list_instances()$items
   machinesType <- instances$machineType
   machineType <- str_extract(machinesType, "(?<=\\/)[[a-z]][[0-9]]+\\-.+?$")
   status <- instances$status
   instanceTable <- data.frame(ID = seq_along(instances$name),
                               names = instances$name,
                               types = machineType,
                               status = status,
                               stringsAsFactors = F)
   return(instanceTable)
}

serverGCEbutton <- function(inputName, GCRfunc) {
   observeEvent(inputName, {
      instances <- fetchInstanceList()
      print(instances)
      selected <- as.integer(input$instanceNo)
      GCRfunc(instances$names[selected])
      gceConnection <- prospX()
      if (!is.null(gceConnection)) {
         instances$status[positionInList] <- "CONNECT-PROSPX"
         positionInList <- str_which(instances$names,
                                     gceConnection$name)
      }
   })
}   
