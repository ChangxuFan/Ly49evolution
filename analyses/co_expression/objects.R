source("~/R_packages/common/base.R")
source("~/R_packages/common/sc.R")
library(scFanc)
SAMPLES.vivier.2022 <- c(paste0("Gut_rep", 1:3), paste0("Liver_rep", 1:3), 
                         paste0("SG_rep", 1:2), paste0("Spleen_rep", 1:2))

KLRA.genes <- paste0("Klra", c(1, 3:10, "11-ps", "13-ps", "14-ps", 17))
