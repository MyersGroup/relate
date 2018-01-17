
pop <- as.matrix(commandArgs(trailingOnly=TRUE))

label <- pop[1]
for(i in 2:length(pop)){
  label <- paste(label, pop[i], sep = "")
}

sample <- read.table("data/1000GP_Phase3_sub.sample",skip = 1)
filename <- paste("1000GP_Phase3_",label,".sample", sep = "")

index <- which(sample[,2] %in% pop)
write("ID POP GROUP SEX", filename)
write.table(sample[index,], filename, row.names = F, col.names = F, append = T, quote=FALSE)

