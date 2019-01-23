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

ref <- sample(4076, 2038, replace = FALSE)

train <- cleanData[ref[1:1019], ]
candidate <- cleanData[ref[1020:2038], ]

write.csv(train, file = "train.csv")
write.csv(candidate, file = "candidate.csv")

<<<<<<< HEAD
=======
Rpar(pty = "s")
hist(data$Total_person_income, breaks = 30, xlab = "Income", ylab = "Frequency", main = NULL, xaxs = "i", yaxs = "i")
>>>>>>> 0b3958b002346881c0bcf1cd1c5179a9ef39ef6f
