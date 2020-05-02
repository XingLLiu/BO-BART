pop <- read.csv("../../data/full_data.csv")
pop <- sapply(pop[2:9], as.factor)
summary(pop)

MC <- read_csv("MC.csv")
BART <- read_csv("BART.csv")

MC <- MC[-c(1:44), -1]
BART <- BART[-c(1:44), -c(1,10)]

n <- nrow(MC)

MC_summary <- matrix(rep(0, 14), 2, 7, dimnames = list(c("Male", "Female"), 
                                                       c("Edu16", "Edu18", "Edu19", "Edu20", 
                                                         "Edu21", "Edu22", "Other")))
BART_summary <- matrix(rep(0, 14), 2, 7, dimnames = list(c("Male", "Female"), 
                                                         c("Edu16", "Edu18", "Edu19", "Edu20", 
                                                           "Edu21", "Edu22", "Other")))

sex <- c(1, 2)
edu <- c(16, 18, 19, 20, 21, 22)

for (i in 1:2){
  
  for (j in 1:6){
    
    MC_summary[i,j] <- sum(MC$Sex==sex[i] & MC$Education==edu[j])
    BART_summary[i,j] <- sum(BART$Sex==sex[i] & BART$Education==edu[j])
  }
  
  MC_summary[i, 7] <- sum(MC$Sex==sex[i]) - sum(MC_summary[i, 1:6])
  BART_summary[i, 7] <- sum(BART$Sex==sex[i]) - sum(BART_summary[i, 1:6])
}

MC_summary <- MC_summary/n*100
BART_summary <- BART_summary/n*100
