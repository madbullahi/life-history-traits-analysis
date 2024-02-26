library(readxl)
library(tidyverse)


synergy <- read_excel("~/life-history-traits-analysis/Tai/SynAddAnt.xlsx")
synergy


synergy$treatment <- factor(synergy$treatment, levels = c("MP", "PFOS", "PFOA", "MP_PFOA", "PFOS_PFOA", "MP_PFOS", "MP_PFOS_PFOA"))

synergy$Genotypes <- factor(synergy$Genotypes, levels = c("LRV-0-1", "LR2-36-01", "NULL-LRV-0-1", "NULL-LR2-36-01"))

synergy$trait <- factor(synergy$trait, levels = c("Age_maturity", "Size_maturity", "Fecundity", "Interval_brood" ))


# create the plot

SAA <-  ggplot(
  synergy,
  aes(
    x = treatment,
    y = SMD,
    fill = Genotypes
  )
) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(ymin = SMD - sd, ymax = SMD + sd),
    position = position_dodge(0.9),
    width = 0
  ) +
  scale_x_discrete(expand = c(0, 0.3)) +
  facet_wrap(~trait, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "top") +
  theme(panel.grid = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("black", "grey", "red", "blue")) 

SAA

# t-test
t-test <-  read_excel("~/life-history-traits-analysis/Tai/")
t.test(Age_maturity ~ Genotypes, data = synergy)

