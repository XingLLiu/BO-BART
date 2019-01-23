<<<<<<< HEAD
#setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")
# read in data
<<<<<<< HEAD
#data <- read.csv("selected_data.csv")
=======
data <- read.csv("full_data.csv")
>>>>>>> 0b3958b002346881c0bcf1cd1c5179a9ef39ef6f

set.seed(1223)

cleanData <- data[-which(is.na(data$Mobility)), ]
cleanData <- cleanData[-which(is.na(cleanData$Employment)), ]
cleanData <- cleanData[-which(is.na(cleanData$Own_child)), ]
cleanData <- cleanData[-which(cleanData$Total_person_income <= 0), ]

write.csv(cleanData, file = "full_data.csv")
=======
setwd("~/Documents/GitHub/BO-BART/data")
>>>>>>> augustDev

library(mlogit)

# read in data
data <- read.csv("full_data.csv")

png("hist_pop.png", width = 450, height = 450)
Rpar(pty = "s")
hist(data$Total_person_income/1000, breaks = 30, xlab = "Total personal income (thousands of US$)", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()

<<<<<<< HEAD
<<<<<<< HEAD
=======
Rpar(pty = "s")
hist(data$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
>>>>>>> 0b3958b002346881c0bcf1cd1c5179a9ef39ef6f
=======
png("white.png", width = 600, height = 600)
Rpar(pty = "s")
hist(data[data$Race == 1, ]$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()

png("edu.png", width = 700, height = 500)
par(mfrow = c(1,2), pty = "s")
hist(data[data$Education <= 16 & data$Education >= 1, ]$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
hist(data[data$Education <= 24 & data$Education >= 17, ]$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
dev.off()
>>>>>>> augustDev
