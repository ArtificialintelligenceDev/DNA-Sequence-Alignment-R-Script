# Sequence Alignment
# Joseph Santos

#Get latest Bioconductor packages.

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("msa")

# Libraries
install.packages("seqinr")
install.packages("phangorn")
library(Biostrings)
library(seqinr)
library(msa)
library(ape)
library(phangorn)

# Read prok.fasta into a new variable called prokaryotes.
prokaryotes <- read.fasta(file="prok.fasta",seqtype="DNA")

# Read the first two sequences in prokaryotes into seq1 and seq2 as character strings.

seq1<-as.character(prokaryotes[[1]])
seq1=paste(seq1,collapse = "")

seq2<-as.character(prokaryotes[[2]])
seq2=paste(seq2,collapse = "")

# Align seq1 and seq2 using the default settings of Biostrings.
pairalign <- pairwiseAlignment(pattern = seq2,subject = seq1)

# convert the alignment to FASTA.
pairalignString = BStringSet( c( toString( subject(pairalign) ), toString(pattern(pairalign))))
writeXStringSet(pairalignString,"aligned.txt",format = "FASTA")

# Load the data
coxgenes <- read.fasta(file = "cox1multi.fasta",seqtype = "AA")
cox1 <- as.character(coxgenes[[1]])
cox2 <- as.character(coxgenes[[2]])

# Incredibly simple dotplot
dotPlot(cox1,cox2,main="Human vs Mouse Cox1 Doptplot")

# Second dotplot (windowed)
dotPlot(cox1, cox2, wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

# Residues 1 -> 100 of the two sequences.
dotPlot(cox1[1:100], cox2[1:100], wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 first 100 AA Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

# Read files as character sets.
coxAA <- readAAStringSet("cox1multi.fasta")
prokDNA <- readDNAStringSet("prok.fasta")

#MSA routine (aligns the cox1 protein sequences using CLUSTALW).
coxAligned <- msa(coxAA)

# Align the DNA sequences.
prokAligned <- msa(prokDNA)

# Save the prokaryote sequence alignment as a stringset.
prokAlignStr=as(prokAligned,"DNAStringSet")

# Writes the stringset you created above as a FASTA file in your working directory.
writeXStringSet(prokAlignStr,file="prokAligned.fasta")

# Export the amino acid alignment contained within coxAligned as FASTA.
coxAlignStr = as(coxAligned,"AAStringSet")
writeXStringSet(coxAlignStr,file="coxAligned.fasta")

write.phylip(coxAligned,"coxAligned.phylip")

# Convert prokaryotic alignment to seqinr format.
prokAligned2 <-msaConvert(prokAligned,type = "seqinr::alignment")

# Generate a distance matrix.
prokdist <- dist.alignment(prokAligned2, "identity")

# neighbor-joining distance tree.
prokTree <- nj(prokdist)
plot(prokTree)

# Generate phylogenetic data (PhyDat) for parsimony tree.
prokAligned3 <- msaConvert(prokAligned, type="phangorn::phyDat")

# Pratchet generates parsimony trees.
ParsTree <- pratchet(prokAligned3)

# Calculate the likelihood using the pml function.
fit <- pml(prokTree, prokAligned3)

# Optimize this tree using optim.pml.(Jukes-Cantor model)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")

bootstrapped <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=FALSE, control = pml.control(trace=0))

# Visulaize bootstrapped variable.
plotBS(midpoint(fitJC$tree), bootstrapped, p = 50, type="p")
