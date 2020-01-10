options(scipen=999)
# Assumes that input is of form
#  //(some text)
#  segsites: (num_segsites)
#  position: (positions)
#  10101110 etc

#Input filename
infile     <- commandArgs(trailingOnly = T)[1]
#output filename (without file extension)
outfile    <- commandArgs(trailingOnly = T)[2]
#if nsites argument is specified, multiply this to pos, otherwise multiply 1
if(length(commandArgs(trailingOnly = T)) == 3){
  nsites       <- as.numeric(commandArgs(trailingOnly = T)[3])
}else if(length(commandArgs(trailingOnly = T)) > 3 || length(commandArgs(trailingOnly = T)) < 2){

  cat("###################################################\n")
  cat("Usage: Rscript ms2haps.R infile.ms outfile nsite\n\n")
  cat("infile.ms: Input filename with file extension.\n")
  cat("outfile: Output filename withouth file extension.\n")
  cat("nsites: (Optional) Number of simulated sites. Default value is 1. This is multiplied to the positions.\n")
  cat("###################################################\n")

  quit(status=1)
}else{
  nsites <- 1
}

# read file as vector of strings
infile_lines <- readLines(infile)
startlines   <- which(startsWith(infile_lines, prefix = "//"))

if(length(startlines) == 0){
  cat("No input in ms format.\n")

  quit(status = 1)
}

# if only one realisation is in the file:
if(length(startlines) == 1){

  pos  <- unlist(strsplit(infile_lines[startlines+2], split = " "))
  pos  <- round(as.numeric(as.matrix(pos[-1])) * nsites)

  seq  <- infile_lines[(startlines+3):length(infile_lines)]
  seq  <- seq[seq!=""] #remove empty lines at the end
  # convert seq to a matrix
  seq  <- t(matrix(as.numeric(unlist(strsplit(seq, split = ""))), ncol = length(pos), byrow = T))

  # find sites with multiple mutations
  dupl <- which(duplicated(pos))
  if(length(dupl) > 0){
    dupl <- sort(unique(c(dupl, dupl-1)))
    pos  <- pos[-dupl]
    seq  <- seq[-dupl,]
  }
  if(length(pos) == 0){
    cat("No segsites\nBP have to be integers! (Use third argument)")
    quit(status = 1)
  }

  N <- dim(seq)[2]
  L <- dim(seq)[1]

  ##### write outfile.sample #####
  write.table(c("ID_1 ID_2 missing"), paste(outfile,".sample", sep = ""), row.names = F, col.names = F, quote = F)
  write.table(c("0 0 0"), paste(outfile,".sample", sep = ""), row.names = F, col.names = F, quote = F, append = T)
  sample <- as.data.frame(matrix(0,N/2,3))
  sample[,1] <- paste("UNR",1:(N/2), sep = "")
  sample[,2] <- paste("UNR",1:(N/2), sep = "")
  sample[,3] <- 0
  write.table(sample, paste(outfile,".sample", sep = ""), row.names = F, col.names = F, quote = F, append = T)

  ##### write outfile.haps #####
  haps <- cbind(rep(1,L), paste("SNP", 1:L, sep = ""), pos, "A", "T", seq)
  write.table(haps, paste(outfile,".haps", sep = ""), quote = F, row.names = F, col.names = F)

}

# if more than 1 realisations are in the file, output these as different "chromosomes"
if(length(startlines) > 1){

  startlines <- c(startlines, length(infile_lines)+1)

  for(i in 1:(length(startlines)-1)){

    pos  <- unlist(strsplit(infile_lines[startlines[i]+2], split = " "))
    pos  <- round(as.numeric(as.matrix(pos[-1])) * nsites)

    seq  <- infile_lines[(startlines[i]+3):(startlines[i+1]-1)]
    seq  <- seq[seq!=""] #remove empty lines at the end 
    # convert seq to a matrix
    seq  <- t(matrix(as.numeric(unlist(strsplit(seq, split = ""))), ncol = length(pos), byrow = T))

    # find sites with multiple mutations
    dupl <- which(duplicated(pos))
    if(length(dupl) > 0){
      dupl <- sort(unique(c(dupl, dupl-1)))
      pos  <- pos[-dupl]
      seq  <- seq[-dupl,]
    }
    if(length(pos) == 0){
      cat("No segsites\nBP have to be integers! (Use third argument)")
      quit(status = 1)
    }

    N <- dim(seq)[2]
    L <- dim(seq)[1]

    ##### write outfile.sample #####
    write.table(c("ID_1 ID_2 missing"), paste(outfile, "_chr", i, ".sample", sep = ""), row.names = F, col.names = F, quote = F)
    write.table(c("0 0 0"), paste(outfile, "_chr", i, ".sample", sep = ""), row.names = F, col.names = F, quote = F, append = T)
    sample <- as.data.frame(matrix(0,N/2,3))
    sample[,1] <- paste("UNR",1:(N/2), sep = "")
    sample[,2] <- paste("UNR",1:(N/2), sep = "")
    sample[,3] <- 0
    write.table(sample, paste(outfile, "_chr", i, ".sample", sep = ""), row.names = F, col.names = F, quote = F, append = T)

    ##### write outfile.haps #####
    haps <- cbind(rep(i,L), paste("SNP", 1:L, sep = ""), pos, "A", "T", seq)
    write.table(haps, paste(outfile, "_chr", i, ".haps", sep = ""), quote = F, row.names = F, col.names = F)

  }

}
