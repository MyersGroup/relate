if(as.numeric(version$major) < 3 || (as.numeric(version$major) == 3 && as.numeric(version$minor) < 3.1)) {
    stop("Please update your R version to at least 3.3.1.")
}

while(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
while(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
while(!require(cowplot)) install.packages("cowplot", repos = "http://cran.us.r-project.org")

##################################################################################

#plotting parameters

#plot height and width
height              <- 30
width               <- 40
#ratio of tree to poplabels
ratio               <- c(6,1.5)

#tree linewidth
tree_lwd            <- 3
#mutation size
mut_size            <- 10
#size of | indicating population label
poplabels_shapesize <- 10 
#population label text size
poplabels_textsize  <- 30

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

  branches_below_mut <- read.table(paste(filename_plot,".plotcoords.mut", sep = ""), header = T)
  plotcoords$belowmut <- FALSE
  plotcoords$belowmut[plotcoords$branchID %in% branches_below_mut$branchID] <- TRUE

  p <- ggplot()          + geom_segment(data = subset(plotcoords, seg_type != "m"), aes(x = x_begin, xend = x_end, y = y_begin, yend = y_end, color = belowmut), ...) +
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
                             scale_x_continuous(limits = c(0, max(plotcoords$x_begin)+1))

  return(p)

}

AddMutations <- function(filename_plot, years_per_gen, ...){

	plotcoords      <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
	plotcoords[3:4] <- plotcoords[3:4] * years_per_gen

	mut_on_branches     <- read.table(paste(filename_plot,".plotcoords.mut", sep = ""), header = T)

  muts <- subset(plotcoords, seg_type == "v" | seg_type == "t")
	muts <- merge(mut_on_branches, muts, by = "branchID")

  muts %>% group_by(branchID) %>% mutate(y_begin = (1:length(y_begin)) * (max(y_end) - min(y_begin))/(1+length(y_begin)) + min(y_begin), y_end = y_begin ) -> muts

	p <- geom_point(data = muts, aes(x = x_begin, y = y_begin), color = "#6564db", ...)

	return(p)

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
                                  axis.line=element_blank(),
																	axis.text.x=element_blank(),
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
filename_dist        <- argv[7]
years_per_gen        <- as.numeric(argv[8])
snp                  <- argv[9]
filename_plot        <- argv[10]

##################################################################################
if(1){
# run RelateTreeView to extract tmp files for plotting tree
system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode TreeView --anc ", filename_anc," --mut ", filename_mut," --snp_of_interest ", as.integer(snp)," -o ", filename_plot))
if(filename_dist != "nodist"){
  system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode BranchesBelowMutation --anc ", filename_anc," --mut ", filename_mut," --dist ",filename_dist," --snp_of_interest ", as.integer(snp)," -o ", filename_plot))
}else{
  system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode BranchesBelowMutation --anc ", filename_anc," --mut ", filename_mut," --snp_of_interest ", as.integer(snp)," -o ", filename_plot))
}
}
# plot tree
p1 <- TreeView(filename_plot, years_per_gen, lwd = tree_lwd) 

if(1){
if(filename_dist != "nodist"){
  system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode MutationsOnBranches --anc ", filename_anc," --mut ", filename_mut," --haps ", filename_haps," --sample ", filename_sample," --dist ",filename_dist," --snp_of_interest ", as.integer(snp)," -o ", filename_plot))
}else{
  system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode MutationsOnBranches --anc ", filename_anc," --mut ", filename_mut," --haps ", filename_haps," --sample ", filename_sample," --snp_of_interest ", as.integer(snp)," -o ", filename_plot))
}
}

p1 <- p1 + AddMutations(filename_plot, years_per_gen, size = mut_size) #+ scale_y_continuous(trans = "log10") 
      
# some modifications to theme
p1 <- p1 + theme(axis.text.y = element_text(size = rel(2.3)), 
                 legend.title = element_text(size = rel(1)), 
                 legend.text = element_text(size = rel(1)))  +  
           scale_shape_manual(labels = c("unflipped", "flipped"), values = c(16, 17), drop = FALSE) +
           scale_color_manual(labels = c("ancestral", "derived"), values = c("black", "#f24236"), drop = FALSE) +
           guides(color = guide_legend(nrow = 2, title = ""), shape = guide_legend(nrow = 2, title = ""))
           
# plot population labels
p2 <- PopLabels(filename_plot, filename_poplabels, text_size = poplabels_textsize, size = poplabels_shapesize, shape = "|")

##################################################################################

pdf(paste(filename_plot,".pdf", sep = ""), height = height, width = width)
plot_grid(p1,p2, rel_heights = ratio, labels = "", align = "v", ncol = 1)
dev.off()

# delete tmp files for plotting tree
system(paste("rm ",filename_plot, ".plotcoords", sep = ""))
system(paste("rm ",filename_plot, ".plotcoords.mut", sep = ""))
