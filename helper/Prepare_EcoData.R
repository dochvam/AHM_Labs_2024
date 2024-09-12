
tryCatch(library(devtools), 
         error = function(e) {
           install.packages("devtools")
         })

tryCatch(library(EcoData), 
         error = function(e) {
           devtools::install_github(repo = "TheoreticalEcology/EcoData", 
                                    dependencies = T, build_vignettes = T)
         })

library(EcoData)


### Modify the data for Lab 1


