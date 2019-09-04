filename      <- commandArgs(trailingOnly = T)[1]
years_per_gen <- as.numeric(commandArgs(trailingOnly = T)[2])
mu            <- as.numeric(commandArgs(trailingOnly = T)[3])

rate          <- read.table(filename)
rate[,1]      <- rate[,1] * years_per_gen
rate          <- rate[(rate[,1] <= 2e6),]

cat(mean(abs(rate[,2] - mu))/mu < 0.1, "\n")


