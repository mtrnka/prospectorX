library(googleComputeEngineR)
gce_get_project()
gce_list_instances()
prospX <- gce_vm("prospx-demo")
prospX <- gce_vm_stop("prospx-demo")

test <- "echo foo; echo bar"
launchSearch <- "/home/mtrnka/prospector.5.23.0/web/cgi-bin/mssearch.cgi -f /home/mtrnka/P181010_Phil/batchtags/prospX.xml"


prospX
View(gce_list_machinetype()$items)
View(gce_list_images(image_project="prospector-236322")$items)

prospXL <- gce_vm_create(name = "prospx",
                        predefined_type = "n1-standard-1",
                        image_family = "debian-9"
                        #        metadata = list(<key> = <value?)
)

out <- gce_ssh(prospX, runDemoScript, capture_text = T)
prospX <- gce_vm("prospx-demo")
gce_ssh(prospX, runDemoScript, capture_text = "gceRunOutput.txt")
runDemoScript <- "bash ./runDemo.sh"
