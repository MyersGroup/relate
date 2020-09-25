if(as.numeric(version$major) < 3 || (as.numeric(version$major) == 3 && as.numeric(version$minor) < 3.1)) {
    stop("Please update your R version to at least 3.3.1.")
}
while(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
while(!require(cowplot)) install.packages("cowplot", repos = "http://cran.us.r-project.org")

##################################################################################

#plotting parameters

#plot height and width
height              <- 30
width               <- 40
#ratio of tree to poplabels
ratio               <- c(6,2)

#tree linewidth
tree_lwd            <- 3
#mutation size
mut_size            <- 8
#size of | indicating population label
poplabels_shapesize <- 10 
#population label text size
poplabels_textsize  <- 50


##################################################################################

filetype <- function(path){
    f = file(path)
    ext = summary(f)$class
    close.connection(f)
    ext
}

TreeView <- function(filename_plot, years_per_gen, ...){

  plotcoords      <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
  plotcoords[3:4] <- plotcoords[3:4] * years_per_gen

  p <- ggplot()          + geom_segment(data = subset(plotcoords, seg_type != "m"), aes(x = x_begin, xend = x_end, y = y_begin, yend = y_end), ...) +
                            #geom_point(data = subset(plotcoords, seg_type == "m"), aes(x = x_begin, y = y_begin), color = "black", size = 2) +
                            theme(text = element_text(size=30),
                                  axis.line=element_blank(),axis.text.x=element_blank(),
                                  axis.ticks=element_blank(),
                                  axis.title.x=element_blank(),
                                  axis.title.y=element_blank(),legend.position="top",
                                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                  panel.grid.minor=element_blank(),plot.background=element_blank(), 
                                  strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"), 
                                  legend.key = element_blank(),
                                  legend.key.width= unit(3, "line"),
                                  legend.key.height= unit(1.5, "line"),
                                  legend.text=element_text(size=35),
                                  strip.text = element_text(face="bold"), plot.margin = margin(t = 0, r = 20, b = 60, l = 30, unit = "pt")) + 
                             scale_x_continuous(limits = c(0, max(plotcoords$x_begin)+1)) + scale_y_continuous(limits = c(0, max(plotcoords$y_end)))

  return(p)

}

AddMutations <- function(filename_plot, years_per_gen, ...){

  plotcoords      <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
  plotcoords[3:4] <- plotcoords[3:4] * years_per_gen

  mut_on_branches     <- read.table(paste(filename_plot,".plotcoords.mut", sep = ""), header = T)
  all_mut_on_branches <- table(mut_on_branches[,2])

  muts <- subset(plotcoords, seg_type == "m")
  mut_on_branches <- read.table(paste(filename_plot,".plotcoords.mut", sep = ""), header = T)
  for(i in 1:nrow(all_mut_on_branches)){
    index                   <- which(mut_on_branches$branchID == rownames(all_mut_on_branches)[i])
    mut_on_branches[index,"branchID"] <- paste(rownames(all_mut_on_branches)[i], 1:length(index))
  }
  all_mut_on_branches <- table(muts$branchID)
  for(i in 1:nrow(all_mut_on_branches)){
    index                   <- which(muts$branchID == rownames(all_mut_on_branches)[i])
    muts[index,"branchID"]  <- paste(rownames(all_mut_on_branches)[i], 1:length(index))
  }
  muts <- merge(muts, mut_on_branches, by = "branchID")
  ord           <- order(muts$pos, decreasing = F)
  muts          <- muts[ord,]
  muts          <- cbind(muts, id = 1:nrow(muts))
  muts$id       <- muts$id - muts$id[min(which(muts$pos >= snp))]

  #if(filetype(filename_mut) != "gzfile"){
    #mutation_file <- as.data.frame(fread(filename_mut, sep = ";"))
  #}else{
    #mutation_file <- as.data.frame(fread(paste("zcat",filename_mut), sep = ";"))
  #}
  mutation_file <- read.table(filename_mut, header = T, sep = ";")  
  colnames(mutation_file)[2] <- "pos"

  muts <- merge(muts, mutation_file[,c("pos","is_flipped")], by = "pos")
  muts$is_flipped[muts$is_flipped == 0] <- "unflipped"
  muts$is_flipped[muts$is_flipped == 1] <- "flipped"
  muts$is_flipped <- factor(muts$is_flipped, levels = c("unflipped", "flipped"))

  p <- geom_point(data = muts, aes(x = x_begin, y = y_begin, color = is_flipped), ...)

  return(p)
  #geom_text(data = muts, aes(x = x_begin + 6, y = y_begin, label = id), size = 8) +

}


PopLabels <- function(filename_plot, filename_poplabels, text_size = 100, ...){

  plotcoords <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
  poplabels  <- read.table(filename_poplabels, header = T)[,2:4]

  tips <- subset(plotcoords, seg_type == "t")
	if(all(is.na(poplabels[,3])) || any(poplabels[,3] != 1)){
		tips <- cbind(tips, population = poplabels[ceiling((tips$branchID+1)/2),1], region = poplabels[ceiling((tips$branchID+1)/2),2])
	}else{
		tips <- cbind(tips, population = poplabels[tips$branchID+1,1], region = poplabels[tips$branchID+1,2])
	}
  unique_region <- unique(poplabels[,2])
  p <- ggplot() + geom_point(data = tips, aes(x = x_begin, y = population, color = population), ...) +
                   theme(text = element_text(size=text_size),
                                  axis.line=element_blank(),axis.text.x=element_blank(),
                                  axis.ticks=element_blank(),
                                  axis.title.x=element_blank(),
                                  axis.title.y=element_blank(),legend.position="bottom",
                                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                  panel.grid.minor=element_blank(),plot.background=element_blank(), 
                                  strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                  strip.text = element_text(face="bold"), plot.margin = margin(t = 0, r = 20, b = 60, l = 60, unit = "pt")) +
                   guides(color = F, shape = F) + 
                   scale_x_continuous(limits = c(0, max(tips$x_begin)+1))

  return(p)

}

#################################################################################

argv <- commandArgs(trailingOnly = T)

PATH_TO_RELATE       <- argv[1] 
filename_haps        <- argv[2]
filename_sample      <- argv[3]
filename_poplabels   <- argv[4]
filename_anc         <- argv[5]
filename_mut         <- argv[6]
years_per_gen        <- as.numeric(argv[7])
snp                  <- argv[8]
filename_plot        <- argv[9]

##################################################################################

# run RelateTreeView to extract tmp files for plotting tree
system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode TreeView --anc ", filename_anc," --mut ", filename_mut," --snp_of_interest ", as.integer(snp)," -o ", filename_plot))
system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode MutationsOnBranches --anc ", filename_anc," --mut ", filename_mut," --haps ", filename_haps," --sample ", filename_sample," --snp_of_interest ", as.integer(snp)," -o ", filename_plot))

# plot tree
p1 <- TreeView(filename_plot, years_per_gen, lwd = tree_lwd) + 
      AddMutations(filename_plot, years_per_gen, size = mut_size) #+ scale_y_continuous(trans = "log10") 
      
# some modifications to theme
p1 <- p1 + theme(axis.text.y = element_text(size = rel(2.3)), 
                 legend.title = element_text(size = rel(1)), 
                 legend.text = element_text(size = rel(1)))  +  
           scale_color_manual(labels = c("unflipped", "flipped"), values = c("red", "blue"), drop = FALSE) +
           guides(color = guide_legend(nrow = 2, title = "")) 

# plot population labels
p2 <- PopLabels(filename_plot, filename_poplabels, text_size = poplabels_textsize, size = poplabels_shapesize, shape = "|")

##################################################################################

pdf(paste(filename_plot,".pdf", sep = ""), height = height, width = width)
plot_grid(p1,p2, rel_heights = ratio, labels = "", align = "v", ncol = 1)
dev.off()

# delete tmp files for plotting tree
system(paste("rm ",filename_plot, ".plotcoords", sep = ""))
system(paste("rm ",filename_plot, ".plotcoords.mut", sep = ""))
