library(RColorBrewer)
library(ggplot2)
library(reshape2)
require(gridExtra)
require(grid)
#library(ggmap)
source("../../arg_inference/Rfunctions/PublicationTheme.R")
source("../../arg_inference/Rfunctions/multiplot.R")

get_plot <- function(pop, ep){

  dir=paste("result_",pop,"/", sep = "")

  pop_labels <- as.matrix(read.table(paste(dir,"relate_",pop,"_bypop.coal", sep=""), nrows = 1))
  epoche <- as.numeric(as.matrix(read.table(paste(dir,"relate_",pop,"_bypop.coal", sep=""), skip = 1, nrows = 1)))
  coal <- as.matrix(read.table(paste(dir,"relate_",pop,"_bypop.coal", sep=""), skip = 2))
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
  pop_size <- melt(pop_size, id.vars = c("time","type"))

  p1 <- ggplot(data = pop_size) + geom_step(aes(time, value, color = type), lwd = 1.2) +
    theme_Publication(base_size = 40, legend.position="bottom") + 
    scale_x_continuous_sci_Publication(limits = c(1e2,1e7), trans = "log10") +
    scale_y_continuous_sci_Publication(limits = c(1e2,1e6), trans = "log10") + annotation_logticks() +
    ylab("population size") +
    xlab("years ago") + 
    scale_color_custom()


  ###################################################

  coal=system(paste("ls ",dir,"*byhap.coal",sep=""), intern = T)
  samples <- paste("1000GP_Phase3_",pop,".sample", sep="")

  print(coal)
  print(samples)

  epoches <- read.table(coal[1],skip=1,nrow=1)
  coal <- read.table(coal[1],skip = 2)
  samples <- as.matrix(read.table(samples[1], skip = 1)[,2])
  samples <- as.vector(t(cbind(samples,samples)))
  N <- length(samples)
  print(N)
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

  df <- get_matrix(ep)

  sorted_samples <- sort(samples)
  n1 <- min(which(sorted_samples != sorted_samples[1]))
  sample_segment <- data.frame( x1 = c(-1,-1), x2 = c(-1,-1), y1 = c(0,n1-1), y2 = c(n1-1,N), color = c("pop1","pop2"))

  print(c(min(df$value[which(df$value > 0)], na.rm = T), max(df$value, na.rm = T)))
  p2 <- ggplot(df) + geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    geom_segment(data = sample_segment, aes(x = x1, y = y1, xend = x2, yend = y2, color = color), lwd = 5) +
    theme_Publication(base_size = 32, legend.position = "bottom") + guides(fill = F) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.border = element_rect(colour = NA),
          plot.margin=unit(c(0,0,0,0),"mm"),
          plot.title = element_text(hjust = 0.05, vjust = 0.1, size = rel(0.9))) +
    scale_x_continuous_Publication(limits = c(-1,N)) +
    scale_y_continuous_Publication(limits = c(0,N)) + 
    scale_fill_gradientn(colors = c("red", "yellow", "blue"), limit = c(min(df$value[which(df$value > 0)], na.rm = T), max(df$value, na.rm = T)), trans = "log10", na.value = "red") +
    #scale_fill_continuous(low = "orange", high = "yellow") + 
    scale_color_manual(values = c("#000000", "#F8766D")) 

  if(ep == 17)  p2 <- p2 + ggtitle("120,000 years ago")
  if(ep == 16)  p2 <- p2 + ggtitle("100,000 years ago")
  if(ep == 16)  p2 <- p2 + ggtitle("100,000 years ago")
  if(ep == 13)  p2 <- p2 + ggtitle("30,000 years ago")
  if(ep == 9)  p2 <- p2 + ggtitle("10,000 years ago")
  if(ep == 7)  p2 <- p2 + ggtitle("5000 years ago")
  if(ep == 5)  p2 <- p2 + ggtitle("2500 years ago")

  xmin <- 5.1 
  xmax <- 6.9
  ymin <- 2.2
  ymax <- 3.9
  g2 <- ggplotGrob(p2)
  plot3 <- p1 + annotation_custom(grob = g2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

  return(plot3)

}
#vp <- viewport(width = 0.3, height = 0.3, x = 0.75, y = 0.36)

p1 <- ggplotGrob(get_plot("CHBYRI",17) + ggtitle("(a)"))
p2 <- ggplotGrob(get_plot("GBRYRI",17) + ggtitle("(b)"))
p3 <- ggplotGrob(get_plot("CHBGBR",13) + ggtitle("(c)"))
p4 <- ggplotGrob(get_plot("CHBJPT",9) + ggtitle("(d)"))
p5 <- ggplotGrob(get_plot("GBRFIN",9) + ggtitle("(e)"))
p6 <- ggplotGrob(get_plot("LWKYRI",13) + ggtitle("(f)"))

#pdf("population_structure.pdf", width = 30, height = 20)
png("population_structure.png", width = 3000, height = 2000)

grid.arrange(rbind(cbind(p1, p2, p3), cbind(p4, p5, p6)))

dev.off()
