
PHE_data <-  read_excel(file.choose(), sheet = "Sheet1") # read the data

PHE_data # print the data

# remove NA values

PHE_data <- PHE_data[complete.cases(PHE_data),]

PHE_data # print the data

# subet the variable Age_maturity, size_maturity, fecundity, and interval_brood

PHE_data_traits <- PHE_data[,c("Replicate","Genotype","Treatment","Lake_Phase","Age_maturity", "Size_maturity", "Fecundity", "Interval_brood")]

PHE_data_traits # print the data

# log transformation of the variables

PHE_data_traits <- log(PHE_data_traits[,5:8] +0.5)
PHE_data_traits # print the data


# manova

manova_results <- manova(cbind(Age_maturity, Size_maturity, Fecundity, Interval_brood) ~ Treatment * Lake_Phase, data = PHE_data)
summary(manova_results) # print the results

# interpret the results

print(summary.aov(manova_results))

# post hoc test [ perform the Tukey HSD test]

library(emmeans) #  load the emmeans package

emmeans_results <- emmeans(manova_results, c("Treatment", "Lake_Phase"))

Tukey_result <-  pairs(emmeans_results, adjust = "tukey") # print the results

Tukey_result
 
# save the results

write.csv(emmeans_results, "emmeans_results.csv")

# save the pairwise comparison results

write.csv(Tukey_result, "Tukey_result.csv")


# visualise the results

library(ggplot2) # load the ggplot2 package
library(car)

emmip(emmeans_results, ~ Treatment | Lake_Phase, CIs = TRUE) # plot the results

# save the plot

ggsave("emmeans_results.png")

