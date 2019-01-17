#setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")
# read in data
#data <- read.csv("selected_data.csv")

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

