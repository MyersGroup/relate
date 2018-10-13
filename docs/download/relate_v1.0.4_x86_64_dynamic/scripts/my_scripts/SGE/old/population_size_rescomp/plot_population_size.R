library(RColorBrewer)
library(ggplot2)
library(reshape2)
require(gridExtra)
require(grid)
source("../../../arg_inference/Rfunctions/PublicationTheme.R")
source("../../../arg_inference/Rfunctions/multiplot.R")

pop <- commandArgs(trailingOnly=TRUE)[1]

pop_labels <- as.matrix(read.table(paste("relate_",pop,"_bypop.coal", sep=""), nrows = 1))
epoche <- as.numeric(as.matrix(read.table(paste("relate_",pop,"_bypop.coal", sep=""), skip = 1, nrows = 1)))
coal <- as.matrix(read.table(paste("relate_",pop,"_bypop.coal", sep=""), skip = 2))

if(dim(coal)[1] == 1){
 
  pop_indices <- coal[c(1,2)] + 1
  coal <- as.matrix(1/coal[-c(1,2)])/2
  pop_size <- data.frame(time = 28* epoche, pop_size = coal, type = paste(pop_labels[pop_indices[1]],pop_labels[pop_indices[2]], sep = "-"))
  print(pop_size)

}else{
  
  pop_indices <- coal[,c(1,2)] + 1
  coal <- as.matrix(1/coal[,-c(1,2)])/2
  pop_size <- data.frame(time = 28* epoche, pop_size = coal[1,], type = paste(pop_labels[pop_indices[1,1]],pop_labels[pop_indices[1,2]], sep = "-"))
  for(i in 2:dim(coal)[1]){
    if(pop_indices[i,1] < pop_indices[i,2]){
      label1 <- pop_labels[pop_indices[i,1]]
      label2 <- pop_labels[pop_indices[i,2]]
    }else{
      label1 <- pop_labels[pop_indices[i,2]]
      label2 <- pop_labels[pop_indices[i,1]]
    }
    pop_size <- rbind(pop_size, data.frame(time = 28* epoche, pop_size = coal[i,], type = paste(label1,label2, sep = "-")))
  }

}
pop_size <- melt(pop_size, id.vars = c("time","type"))

p1 <- ggplot(data = pop_size) + geom_step(aes(time, value, color = type), lwd = 1.2) +
  theme_Publication(legend.position="bottom", base_size = 26) + 
  scale_x_continuous_sci_Publication(limits = c(1e2,1e7), trans = "log10") +
  scale_y_continuous_sci_Publication(limits = c(1e2,1e6), trans = "log10") + annotation_logticks() +
  ylab("population size") +
  xlab("time in years") + 
  scale_color_custom()


###################################################

coal=system("ls *byhap.coal", intern = T)
N=system("head -1 *.arg", intern=T)[2]
N=as.numeric(unlist(strsplit(unlist(N), "[^0-9]+")))[2]
samples <- paste("../1000GP_Phase3_",pop,".sample",sep="")

print(coal)
print(N)

epoches <- read.table(coal[1],skip=1,nrow=1)
coal <- read.table(coal[1],skip = 2)
samples <- as.matrix(read.table(samples[1], skip = 1)[,2])
N <- dim(samples)[1] * 2
samples <- as.vector(t(cbind(samples,samples)))
usamples <- sort(unique(samples))

get_matrix <- function(ep){

  M <- matrix(0,N,N)
  M[lower.tri(M, diag=F)] <- coal[,ep]
  M <- M + t(M)
  
  order <- order(samples)
  M <- M[order,order]
  char <- paste(round(28*epoches[ep]), "-", round(28*epoches[ep+1]))
  
  return(cbind(melt(M),char))
  
}

df <- get_matrix(6)

sorted_samples <- sort(samples)
n1 <- min(which(sorted_samples != sorted_samples[1]))
sample_segment <- data.frame( x1 = c(0,0), x2 = c(0,0), y1 = c(0,n1-1), y2 = c(n1-1,N), color = c("pop1","pop2"))
print(sample_segment)
print(head(df))

p2 <- ggplot(df) + geom_tile(aes(x = Var1, y = Var2, fill = value)) +
             geom_segment(data = sample_segment, aes(x = x1, y = y1, xend = x2, yend = y2, color = color), lwd = 5) +
             theme_Publication(legend.position = "bottom") + guides(fill = F) +
             theme(axis.line=element_blank(),axis.text.x=element_blank(),
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none",
                   panel.border = element_rect(colour = NA),
                   plot.margin=unit(c(0,0,0,0),"mm"),
                   plot.title = element_text(hjust = 0.05, vjust = 0.1)) +
             scale_x_continuous_Publication(limits = c(0,N)) +
             scale_y_continuous_Publication(limits = c(0,N)) + 
             scale_fill_gradientn( colors = c("lightblue","darkgreen"), trans = "log10", na.value = "blue") +
             #scale_fill_continuous(low = "orange", high = "yellow") + 
             scale_color_manual(values = c("#000000", "#F8766D")) + 
             ggtitle("2500 years ago")

#xmin <- 0.6
#xmax <- 0.98
#ymin <- 0.6
#ymax <- 0.98
#g2 <- ggplotGrob(p2)
#plot3 <- p1 + annotation_custom(grob = g2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

vp <- viewport(width = 0.3, height = 0.3, x = 0.75, y = 0.36)

#png("population_structure.png", width = 1000, height = 1000)
#plot(p1)
#plot(p2, vp = vp)
#dev.off()

pdf("population_structure.pdf", width = 10, height = 10)
plot(p1)
plot(p2, vp = vp)
dev.off()

