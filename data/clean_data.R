setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")

library(mlogit)

# read in data
data <- read.csv("../../data/full_data.csv")

png("population.png", width = 600, height = 600)
rpar(pty = "s")
hist(data$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()

png("white.png", width = 600, height = 600)
rpar(pty = "s")
hist(data[data$Race == 1, ]$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()