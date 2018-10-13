
library(reshape2)
library(ggplot2)
source("../../arg_inference/Rfunctions/multiplot.R")
source("../../arg_inference/Rfunctions/PublicationTheme.R")

GetRate <- function(pop){

  coln <- as.matrix(read.table(paste("result_",pop,"/relate_",pop,".rate", sep=""),nrows = 1))
  rate <- as.matrix(read.table(paste("result_",pop,"/relate_",pop,".rate", sep=""), skip = 1))
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
  rate$pop <- substr(pop,1,3)

  return(rate)

}

GetMeanRate <- function(pop){

  mean_rate <- read.table(paste("result_",pop,"/relate_",pop,"_avg.rate", sep=""))
  mean_rate <- data.frame(time = mean_rate[,1], value = mean_rate[,2]/3, pop = substr(pop,1,3))

  return(mean_rate)

}

PlotMutationRate2 <- function(rate, mean_rate){

  p1 <- ggplot(data=rate) + 
    geom_line(aes(x = 28*time, y = value, color = upstream_base, linetype = downstream_base, group = flanking_nucleotides), lwd = 0.9) + 
    geom_line(data=mean_rate, aes(x = 28*time, y = value), col = "black", lwd = 0.9) + 
    scale_linetype_manual(values=c(1,2,3,10)) +
    facet_grid(mutation ~ pop) +
    ylab("mutation rate") +
    xlab("years ago") +
    guides(color=guide_legend(title="upstream base "), linetype=guide_legend(title="downstream base ")) +
    theme_Publication(base_size = 45, legend.position = "bottom", panel.spacing = unit(6, "lines")) + theme(legend.title = element_text(size= rel(0.8))) +
    scale_y_continuous_sci_Publication(limits = c(1e-10,1e-6), trans="log10") + annotation_logticks(side = "b") +
    scale_x_continuous_sci_Publication(limits = c(1e3,1e6), breaks = c(10^{3:6}), trans="log10") + 
    scale_color_manual(values = c("red", "blue", "orange", "magenta"))

  #median_rate  <- aggregate(rate$value, list(rate$time), median)
  #colnames(median_rate) <- c("time", "value")
  mean_rate$value[which(is.na(mean_rate$value) == T)] <- 1
  mean_rate$value[which(mean_rate$value == 0)] <- 1
  mapply(function(x){

         #divide rate$value by entry in mean_rate that matches pop and time
         x1 <- rate$pop[x]
         x2 <- rate$time[x]
         rate$value[x] <<- rate$value[x]/mean_rate$value[intersect(which(mean_rate$pop == x1), which(mean_rate$time == x2))]

  },1:dim(rate)[1])


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
  mean_by_type <- aggregate(rate_tmp$value, list(rate_tmp$flanking_nucleotides, rate_tmp$mutation, rate_tmp$pop), sum) 
  #mean_by_type[,3] <- mean_by_type[,3]/min(unique_time[which(unique_time >= 1e6/28)]) 
  #rate_now <- subset(rate, rate$time == 0)
  #print(head(mean_by_type))

  foo <- mapply(function(l){
                index <- which(mean_by_type[,1] == rate$flanking_nucleotide[l])
                index <- intersect(index, which(mean_by_type[,2] == rate$mutation[l]))
                index <- intersect(index, which(mean_by_type[,3] == rate$pop[l]))
                rate$value[l] <<- rate$value[l]/mean_by_type[index,4]
                #index <- which(rate_now$flanking_nucleotide == rate$flanking_nucleotide[l])
                #index <- intersect(index, which(rate_now$mutation == rate$mutation[l]))
                #rate$value[l] <<- rate$value[l]/rate_now$value[index]
  },1:dim(rate)[1])

  p2 <- ggplot(data=rate) + 
    #geom_hline(aes(yintercept = 1.25e-8), lwd = 0.9, linetype = 10) +
    #geom_line(aes(x = 28*time, y = value, color = flanking_nucleotides), lwd = 1.1) + 
    geom_line(aes(x = 28*time, y = value, color = upstream_base, linetype = downstream_base, group = flanking_nucleotides), lwd = 0.9) + 
    scale_linetype_manual(values=c(1,2,3,10)) +
    facet_grid(mutation ~ pop) +
    ylab("normalized mutation rate") +
    xlab("years ago") +
    guides(color=guide_legend(title="upstream base "), linetype=guide_legend(title="downstream base ")) +
    theme_Publication(base_size = 45, legend.position = "bottom", panel.spacing = unit(6, "lines")) + theme(legend.title = element_text(size= rel(0.8))) +
    scale_y_continuous_Publication(limits = c(0.5,2)) + annotation_logticks(side = "b") + 
    scale_x_continuous_sci_Publication(limits = c(1e3,1e6), breaks = c(10^{3:6}), trans="log10") +
    #scale_color_brewer(type = 'qual', palette = 'GnBu')
    scale_color_manual(values = c("red", "blue", "orange", "magenta"))

  return(list(p1 = p1, p2 = p2))

}

rate <- function(pop){

  rate <- GetRate(paste(pop[1],pop[1],sep = ""))
  for(i in 2:length(pop)){
    rate <- rbind(rate, GetRate(paste(pop[i],pop[i],sep = "")))
  }

  return(rate)

}

mean_rate <- function(pop){

  mean_rate <- GetMeanRate(paste(pop[1],pop[1],sep = ""))
  for(i in 2:length(pop)){
    mean_rate <- rbind(mean_rate, GetMeanRate(paste(pop[i],pop[i],sep = "")))
  }

  return(mean_rate)

}

#pop <- c("GBR","FIN","CHB","JPT","YRI","LWK")
pop <- c("GBR","CHB","YRI")
rate <- rate(pop)
mean_rate <- mean_rate(pop)

p <- PlotMutationRate2(rate, mean_rate)

file_plot <- paste("MutationRateThroughTime.pdf",sep="")
pdf(file_plot,
    width     = 30,
    height    = 40,
    pointsize = 7,
    )

p1 <- p[[1]] + ggtitle("(a)") + theme(plot.title = element_text(hjust = -0.02))
p2 <- p[[2]] + ggtitle("(b)") + theme(plot.title = element_text(hjust = -0.02))

grid.arrange(ggplotGrob(p1),ggplotGrob(p2),layout_matrix = rbind(c(1),c(2)))

dev.off()
