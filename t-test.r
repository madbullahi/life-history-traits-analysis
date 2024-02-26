# create vectors for the null and observed values: Age at maturity (LRV-0-1) (note: this is an hypothetical data)

null_values <-  c(7,
                  8,
                  8,
                  9
)

observed_values <- c(6,
                     7,
                     7,
                     8
)


# Perform uparies t-test

t.test(null_values, observed_values)


# alternative approach

# calculate the mean and standard devaiaation of observed values and predicted values:

predicted_mean <- -1.4999

predicted_sd <- 0.5

observed_mean <- 2.1499995

observed_sd <- 1.4142136

# create the vectors for the observed and predicted values

predicted_values <- rnorm(4, mean = predicted_mean, sd = predicted_sd)

observed_values <- rnorm(4, mean = observed_mean, sd = observed_sd)

# perform un[paired t-test]

t_test_result <- t.test(predicted_values, observed_values)

t_test_result
