setwd("~/Documents/GitHub/BO-BART/data")

library(mlogit)

# read in data
data <- read.csv("full_data.csv")

png("hist_pop.png", width = 450, height = 450)
Rpar(pty = "s")
hist(data$Total_person_income/1000, breaks = 30, xlab = "Total personal income (thousands of US$)", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()

png("white.png", width = 600, height = 600)
Rpar(pty = "s")
hist(data[data$Race == 1, ]$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()

png("edu.png", width = 700, height = 500)
par(mfrow = c(1,2), pty = "s")
hist(data[data$Education <= 16 & data$Education >= 1, ]$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
hist(data[data$Education <= 24 & data$Education >= 17, ]$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()