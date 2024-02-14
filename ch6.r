1 + 1
2 + 5 - 3
?sqrt()


require(stats)
require(graphics)
xx <- -9:9
plot(xx, sqrt(xx), col="red")
lines(spline(xx, sqrt(abs(xx)), n=101), col="pink")

# convert to upper case
toupper("Hello World")

# is 3 greater than 2?

3 > 2


my_result <- sqrt(abs(-9 ^ 3))
my_result


my_int <- 12L
my_int


??arrange()


update.packages()


# Update R 

UpdateR()

# for exploring and visualising data

library(tidyverse)

# for reading Excel files

library(readxl)

# for descriptive statistics

library(psych)

# for writing data to Excel

library(writexl)

is.vector(1:10)


my_number <- 8.2
sqrt(my_number)


is.vector(my_number) # check if my_number is a vector

length(my_number) # check the length of my_number


my_numbers <- c(5, 6, 7, 8, 9, 10)

is.vector(my_numbers) # check if my_number is a vector


length(my_numbers) # check the length of my_number


str(my_numbers) # check the structure of my_numbers


sqrt(my_numbers) # square root of my_numbers


roster_names <- c("John", "Paul", "George", "Ringo")

toupper(roster_names) # convert to upper case


my_vec <- c("A", 2, "C", 4, "E")

str(my_vec) # check the structure of my_vec

# Get the third element of roster_names

roster_names[3]

roster_names[1:3]

roster_names[2:length(roster_names)]

# Get the second and fifth element

roster_names[c(2, 4)]

my_numbers


sqrt(my_numbers[c(2:length(my_numbers))]) #

sqrt(my_numbers[2:3])
my_numbers

# from Excel Tables to R dataframe
# lets build and print a dataframe from scratch.

# create a dataframe

roster <- data.frame(
  Name = c("John", "Paul", "George", "Ringo"),
  Age = c(22, 23, 24, 25),
  Height = c(5.8, 5.9, 6.0, 5.7),
  Weight = c(160, 165, 170, 155)
)
 roster
 # a dataframe is a list of vectors of equal length(i.e a collection of vectors)
 
 data()

 head(iris) 

 # lets confirm that iris is really a dataframe
 
 is.data.frame(iris)


