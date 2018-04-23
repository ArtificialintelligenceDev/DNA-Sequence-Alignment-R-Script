# Genomic Annotation R Script
# Joseph Santos 04-20-2018 [ORF Finding]

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(Biostrings)
library(seqinr) 

# Input Sequence.
AB003468 <- readDNAStringSet("AB003468.fasta")
AB003468 <- as.character(AB003468)

# Find start codons (START is ATG).
matchPattern("ATG",AB003468)

# Put search sequence into a gen var (increase reusability).
sequence <- AB003468

# Create start and stop codon variables.
start_codon <- "ATG"
stop_codons <- c("TGA","TAA","TAG")

# Create variables to store the results in.
start_pos <- c()
revstart_pos <- c()
stop_pos <- c()
revstop_pos <- c()

# Searching sequence for START and STOP codons.
matches <- matchPattern(start_codon,sequence)

# Finds all instances of start_codon in sequence and puts the results into matches.
start_pos <- c(start_pos, start(matches))

# Use reverseComplement to help us search the other 3 frames. 
revmatches <- matchPattern(reverseComplement(DNAString(start_codon)),sequence)
revstart_pos <- c(revstart_pos, start(revmatches))

# Sort the results.
start_pos <- sort(start_pos)
revstart_pos <- sort(revstart_pos,decreasing = TRUE)

# FOR loop will look for each possible value of STOP codons found in stop_codons.
for(codon in stop_codons){
  matches <- matchPattern(codon, sequence)
  stop_pos <- c(stop_pos, start(matches))
  revmatches <- matchPattern(reverseComplement(DNAString(codon)),sequence)
  revstop_pos <- c(revstop_pos, start(revmatches))
}

# sort STOP positions.
stop_pos <- sort(stop_pos)
revstop_pos <- sort(revstop_pos,decreasing = TRUE)

# K is threshold
k <- 150
lengths <- vector(mode = "numeric")

# Hold the location of the STOPS in each reading frame
stop_pointers <- c(0,0,0)
count <- 0

for(current_start in start_pos){
  frame <- (current_start%%3) +1
  stop_pointer <- stop_pointers[frame]
  if(stop_pointer <= length(stop_pos)&&(stop_pointer == 0 || stop_pos[stop_pointer] < current_start)){
    stop_pointer <- stop_pointer +1
  
  while((stop_pointer <= length(stop_pos)) && ((stop_pos[stop_pointer]<= current_start)
                                               || (((stop_pos[stop_pointer]%%3)+1)!=frame))){
    stop_pointer <- stop_pointer +1
  }
  stop_pointers[frame] <- stop_pointer
  
  
  if(stop_pointer <= length(stop_pos)){
    if((stop_pos[stop_pointer] + 2 - current_start +1) > k){
      count <- count +1
      print(count)
      print("Frame:")
      print(frame)
      print("Start:")
      print(current_start)
      print("Stop:")
      print(stop_pos[stop_pointer])
      print("Length:")
      lengths <- c(lengths, (stop_pos[stop_pointer] + 2 - current_start +1))
      print(stop_pos[stop_pointer] + 2 - current_start +1)
      print("Sequence:")
      print(subseq(sequence, current_start, stop_pos[stop_pointer]+2))
      }}
}
}

# Reverse 3 frames.

revstop_pointers <- c(0,0,0)

for(current_revstart in revstart_pos){
  current_revstart <- current_revstart + 2
  frame <-
    
    
    
    
    
    
    
    
  stop_pointer <- stop_pointers[frame]
  if(stop_pointer <= length(stop_pos)&&(stop_pointer == 0 || stop_pos[stop_pointer] < current_start)){
    stop_pointer <- stop_pointer +1
    
    while((stop_pointer <= length(stop_pos)) && ((stop_pos[stop_pointer]<= current_start)
                                                 || (((stop_pos[stop_pointer]%%3)+1)!=frame))){
      stop_pointer <- stop_pointer +1
    }
    stop_pointers[frame] <- stop_pointer
    
    
    if(stop_pointer <= length(stop_pos)){
      if((stop_pos[stop_pointer] + 2 - current_start +1) > k){
        count <- count +1
        print(count)
        print("Frame:")
        print(frame)
        print("Start:")
        print(current_start)
        print("Stop:")
        print(stop_pos[stop_pointer])
        print("Length:")
        lengths <- c(lengths, (stop_pos[stop_pointer] + 2 - current_start +1))
        print(stop_pos[stop_pointer] + 2 - current_start +1)
        print("Sequence:")
        print(subseq(sequence, current_start, stop_pos[stop_pointer]+2))
      }}
  }
}

