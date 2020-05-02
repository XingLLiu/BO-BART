populationData <- read.csv("../../data/full_data.csv")

cols <- ncol(populationData)
index <- populationData[, 1]
Total_person_income <- populationData[, cols]

populationData <- sapply(populationData[, 2:(cols-1)], as.factor)
populationData <- cbind(index, populationData, Total_person_income)
populationData <- data.frame(populationData)

train <- c()
ref <- c()
candidate <- populationData

set.seed(123)


for (i in 2:(cols-1)) {

  level <- levels(populationData[, i])
  nlevel <- length(level)

  for (j in 1:nlevel){

    ss <- candidate[which(candidate[, i] == level[j]), ]
    n <- sample(nrow(ss), 1)
    train <- rbind(train, ss[n, ])
    candidate <- candidate[-which(candidate[, 1] == ss[n, 1]), ]
    }
}

write.csv(train[, -1], "train2.csv")
write.csv(candidate[, -1], "candidate2.csv")


