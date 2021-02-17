
#check if packages are installed
#list.of.packages <- c("ggplot2", "grid", "gridExtra")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(grid)
library(gridExtra)

argv <- commandArgs(trailingOnly = T)
filename <- argv[1]
years_per_gen <- as.numeric(argv[2])

#read in population size
groups   <- as.matrix(read.table(paste(filename, ".coal", sep = ""), nrow = 1))
t        <- years_per_gen*t(as.matrix(read.table(paste(filename, ".coal", sep = ""), skip = 1, nrow = 1)))
pop_size <- data.frame(time = numeric(0), pop_size = numeric(0), groups = numeric(0))
num_pops <- round(sqrt(dim(read.table(paste(filename, ".coal", sep = ""), skip = 2))[1]))
#num_pops <- (-1 + sqrt(1 + 8 * num_pops))/2

for(p1 in 1:num_pops){
  for(p2 in 1:p1){
    c        <- as.matrix(read.table(paste(filename, ".coal", sep = ""), skip = (p1-1) * num_pops + p2 + 1, nrow = 1))[-c(1:2)]
    str      <- rep(paste(groups[p1]," - ",groups[p2], sep = ""),length(c))
    pop_size <- rbind(pop_size, data.frame(time = t, pop_size = 0.5/c, groups = str))
  }
}
pop_size$time[which(pop_size$time > 1e7)] <- 1e7

#plot
p1 <- ggplot(pop_size) + 
  geom_step(aes(time, pop_size, color = groups, linetype = groups), lwd = 1.2) +
  scale_x_continuous(limits = c(1e3,1e7), trans="log10") + annotation_logticks(sides = "bl") +  
  scale_y_continuous(trans ="log10") +
  ylab("population size") +
  xlab("years ago") +

pdf(paste(filename,".pdf",sep = ""), width = 10, height = 8)
plot(p1)
dev.off()


