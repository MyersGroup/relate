
argv <- commandArgs(trailingOnly = T)

filename <- paste("relate_",argv[1],".coal", sep = "")

epoche <- read.table(filename, skip = 1, nrow = 1)
coal <- read.table(filename, skip = 2)
rate <- read.table(paste("relate_",argv[1],"_avg.rate", sep = ""))

d <- dim(coal)[2]
coal[3:d] <- coal[3:d]/rate[,2]*1.25e-8
coal[which(is.na(coal))] <- 0

write("group1",filename)
write.table(epoche, filename, row.names = F, col.names = F, quote = F, append = T)
write.table(coal, filename, row.names = F, col.names = F, quote = F, append = T)


