library(readxl)

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


# After running the manova on the dataset, the emmeans package allows you:
# 1. calculate the estimate marginnal means (EMMS) for each combination of the independent variables, while accounting for the other variables in the model.
# 2. perform pairwise comparisons between the EMMs to determine which groups differ significantly from each other.
# 3. visualise the results using interacctiion plots or other graphics.

# load packages
library(rstatix)

# calculate the mahalanobis distance

mahalanobis_distance <- mahalanobis(cbind(PHE_data$Age_maturity,PHE_data$Size_maturity,PHE_data$Fecundity,PHE_data$Interval_brood), 
                                    colMeans(cbind(PHE_data$Age_maturity, PHE_data$Size_maturity, PHE_data$Fecundity, PHE_data$Interval_brood)), 
                                    cov(cbind(PHE_data$Age_maturity, PHE_data$Size_maturity, PHE_data$Fecundity, PHE_data$Interval_brood)))

# identify the outliers

outliers <- PHE_data[mahalanobis_distance > qchisq(0.975, ncol(cbind(PHE_data$Age_maturity,
                                                                     PHE_data$Size_maturity,
                                                                     PHE_data$Fecundity,
                                                                     PHE_data$Interval_brood))),]
print(outliers)

# save the outliers

write.csv(outliers, "outliers.csv")


#  calculate the multivariate centroid

centroid <- aggregate(cbind(PHE_data$Age_maturity, PHE_data$Size_maturity, PHE_data$Fecundity, PHE_data$Interval_brood) ~ Treatment * Lake_Phase, data = PHE_data, FUN = mean)

# plot the 3D scatter plot with centroids
library(scatterplot3d)
scatterplot3d(PHE_data$Age_maturity, PHE_data$Size_maturity, PHE_data$Fecundity, PHE_data$Interval_brood, pch = 16, 
            main = "3D Scatter Plot with Centroids")
points(centroid$Age_maturity, centroid$Size_maturity, centroid$Fecundity, pch = 19, cex = 2)

# save the plot

ggsave("3D_Scatter_Plot_with_Centroids.png")
 



# Emmeans for only the Lake-Phase
 emmeans_LakePhase <- emmeans(manova_results, "Lake_Phase")

 Tukey_LakePhase <- pairs(emmeans_LakePhase, adjust = "tukey")

 # save the results
 
 write.csv(emmeans_LakePhase, "emmeans_LakePhase.csv")
write.csv(Tukey_LakePhase, "Tukey_LakePhase.csv") 
