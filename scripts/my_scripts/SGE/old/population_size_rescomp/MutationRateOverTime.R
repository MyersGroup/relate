
library(reshape2)
library(ggplot2)
source("../../../arg_inference/Rfunctions/multiplot.R")
source("../../../arg_inference/Rfunctions/PublicationTheme.R")

args = commandArgs(trailingOnly=TRUE)
pop <- args[1]

coln <- as.matrix(read.table(paste("relate_",pop,".rate", sep=""),nrows = 1))
rate <- as.matrix(read.table(paste("relate_",pop,".rate", sep=""), skip = 1))
colnames(rate) <- c("time",coln)
rate <- as.data.frame(rate)

#rate[,1] <- (rate[,1] + c(rate[-1,1],exp(0.7*22)))/2

rate <- melt(rate, id.vars = "time")

downstream_base <- rep('X',dim(rate)[1])
upstream_base   <- rep('X',dim(rate)[1])
ancestral_base  <- rep('X',dim(rate)[1])
flanking_nucleotides <- rep('X', dim(rate)[1])
CpG                  <- rep('X', dim(rate)[1])
mutation        <- rep('X',dim(rate)[1])

foo <- mapply(function(i){
  
  split <- strsplit(as.character(rate$variable[i]),split = "")[[1]]
  upstream_base[i] <<- split[1]
  downstream_base[i]   <<- split[2]
  ancestral_base[i]  <<- split[3]
  if(ancestral_base[i] == 'C' && downstream_base[i] == 'G'){
    CpG[i] <<- "yes"    
  }else{ 
    CpG[i] <<- "no"  
  }
  flanking_nucleotides[i] <<- paste(split[1], split[2], collapse = "", sep=".")
  mutation[i]        <<- paste(split[3:5], collapse="")
  
}, 1:dim(rate)[1])

rate$downstream_base <- downstream_base
rate$upstream_base <- upstream_base
rate$ancestral_base <- ancestral_base
rate$mutation <- mutation
rate$CpG <- factor(CpG, levels=c("yes", "no"))
rate$flanking_nucleotides <- flanking_nucleotides

mean_rate <- read.table(paste("relate_",pop,"_avg.rate", sep=""))
mean_rate <- data.frame(time = mean_rate[,1], value = mean_rate[,2])

p1 <- ggplot(data=rate) + 
  #geom_hline(aes(yintercept = 1.25e-8), lwd = 0.9, linetype = 10) +
  #geom_line(aes(x = 28*time, y = value, color = flanking_nucleotides), lwd = 1.1) + 
  geom_line(aes(x = 28*time, y = value, color = upstream_base, linetype = downstream_base, group = flanking_nucleotides), lwd = 0.9) + 
  geom_line(data=mean_rate, aes(x = 28*time, y = value), col = "black", lwd = 0.9) + 
  scale_linetype_manual(values=c(1,2,3,10)) +
  facet_grid(mutation ~ .) +
  ylab("mutation rate") +
  xlab("time in years") +
  guides(color=guide_legend(title="upstream base "), linetype=guide_legend(title="downstream base ")) +
  theme_Publication(base_size = 25, legend.position = "bottom", panel.spacing = unit(2, "lines")) + theme(legend.title = element_text(size= rel(0.8))) +
  scale_y_continuous_sci_Publication(limits = c(1e-10,1e-6), trans="log10") + annotation_logticks(side = "l") +
  scale_x_continuous_sci_Publication(limits = c(1e2,1e7), breaks = c(10^{2:7}), trans="log10") + 
  scale_color_manual(values = c("red", "blue", "orange", "magenta"))

mean_rate$value[which(is.na(mean_rate$value) == T)] <- 1
mean_rate$value[which(mean_rate$value == 0)]        <- 1
rate$value   <- rate$value/mean_rate$value

p2 <- ggplot(data=rate) + 
  #geom_hline(aes(yintercept = 1.25e-8), lwd = 0.9, linetype = 10) +
  #geom_line(aes(x = 28*time, y = value, color = flanking_nucleotides), lwd = 1.1) + 
  geom_line(aes(x = 28*time, y = value, color = upstream_base, linetype = downstream_base, group = flanking_nucleotides), lwd = 0.9) + 
  scale_linetype_manual(values=c(1,2,3,10)) +
  facet_grid(mutation ~ .) +
  ylab("normalized mutation rate") +
  xlab("time in years") +
  guides(color=guide_legend(title="upstream base "), linetype=guide_legend(title="downstream base ")) +
  theme_Publication(base_size = 25, legend.position = "bottom", panel.spacing = unit(2, "lines")) + theme(legend.title = element_text(size= rel(0.8))) +
  scale_y_continuous_Publication(limits = c(1e-2,1e2), trans = "log10") + 
  scale_x_continuous_sci_Publication(limits = c(1e2,1e7), breaks = c(10^{2:7}), trans="log10") +
  #scale_color_brewer(type = 'qual', palette = 'GnBu')
  scale_color_manual(values = c("red", "blue", "orange", "magenta"))


#rate         <- rate[which(is.na(rate$value) == F),]
unique_time <- unique(rate$time)
rate_tmp <- subset(rate, rate$time < 1e6/28)
epoche_length <- unique_time
epoche_length <- epoche_length[-1] - epoche_length[-length(epoche_length)]
epoche_length <- data.frame(time = unique_time[-length(unique_time)], dur = epoche_length)
total_time <- sum(epoche_length$dur[which(epoche_length$time <= max(unique(rate_tmp$time)))])
for(t in unique(rate_tmp$time)){
  rate_tmp$value[which(rate_tmp$time == t)] <- rate_tmp$value[which(rate_tmp$time == t)] * epoche_length$dur[which(epoche_length$time == t)]/total_time
}
mean_by_type <- aggregate(rate_tmp$value, list(rate_tmp$flanking_nucleotides, rate_tmp$mutation), sum) 
#mean_by_type[,3] <- mean_by_type[,3]/min(unique_time[which(unique_time >= 1e6/28)]) 
#rate_now <- subset(rate, rate$time == 0)

foo <- mapply(function(l){
  index <- which(mean_by_type[,1] == rate$flanking_nucleotide[l])
  index <- intersect(index, which(mean_by_type[,2] == rate$mutation[l]))
  rate$value[l] <<- rate$value[l]/mean_by_type[index,3]
  #index <- which(rate_now$flanking_nucleotide == rate$flanking_nucleotide[l])
  #index <- intersect(index, which(rate_now$mutation == rate$mutation[l]))
  #rate$value[l] <<- rate$value[l]/rate_now$value[index]
},1:dim(rate)[1])

p6 <- ggplot(data=rate) + 
  #geom_hline(aes(yintercept = 1.25e-8), lwd = 0.9, linetype = 10) +
  #geom_line(aes(x = 28*time, y = value, color = flanking_nucleotides), lwd = 1.1) + 
  geom_line(aes(x = 28*time, y = value, color = upstream_base, linetype = downstream_base, group = flanking_nucleotides), lwd = 0.9) + 
  scale_linetype_manual(values=c(1,2,3,10)) +
  facet_grid(mutation ~ .) +
  ylab("normalized mutation rate") +
  xlab("time in years") +
  guides(color=guide_legend(title="upstream base "), linetype=guide_legend(title="downstream base ")) +
  theme_Publication(base_size = 25, legend.position = "bottom", panel.spacing = unit(2, "lines")) + theme(legend.title = element_text(size= rel(0.8))) +
  scale_y_continuous_Publication(limits = c(0,2)) + 
  scale_x_continuous_sci_Publication(limits = c(1e2,1e7), breaks = c(10^{2:7}), trans="log10") +
  #scale_color_brewer(type = 'qual', palette = 'GnBu')
  scale_color_manual(values = c("red", "blue", "orange", "magenta"))



file_plot <- paste("MutationRateThroughTime_",pop,".pdf",sep="")
pdf(file_plot,
    width     = 24,
    height    = 14,
    pointsize = 7,
    #family="jsMath-cmr10"
)

p1 <- p1 + ggtitle("(a)") + theme(plot.title = element_text(hjust = -0.1))
p6 <- p6 + ggtitle("(b)") + theme(plot.title = element_text(hjust = -0.1))

multiplot(p1,p6, cols=2)

dev.off()

