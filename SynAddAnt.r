library(readxl)
library(tidyverse)


synergy <- read_excel("~/life-history-traits-analysis/Tai/SynAddAnt.xlsx")
synergy


synergy$treatment <- factor(synergy$treatment, levels = c("MP", "PFOS", "PFOA", "MP_PFOA", "PFOS_PFOA", "MP_PFOS", "MP_PFOS_PFOA"))

synergy$Genotypes <- factor(synergy$Genotypes, levels = c("LRV-0-1", "LR2-36-01", "NULL-LRV-0-1", "NULL-LR2-36-01"))

synergy$trait <- factor(synergy$trait, levels = c("Age_maturity", "Size_maturity", "Fecundity", "Interval_brood" ))

