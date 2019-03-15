cd tools
git clone git://github.com/gpertea/gclib
git clone git://github.com/gpertea/gffread
cd gffread
make release

#using gffread function to make the transcript fasta file
#gffread can also be used to generate a FASTA file with the DNA sequences for all transcripts in a GFF file. For this operation a fasta file with the genomic sequences have to be provided as well. For example, one might want to extract the sequence of all transfrags assembled from a Cufflinks assembly session. This can be accomplished with a command line like this:
#don't think I will need this

tools/gffread/gffread -w /mnt/lustre/users/k1632479/polyester/Mus_musculus.GRCm38_transcripts.fa -g /mnt/lustre/users/k1632479/grcm38/Ensembl/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa /mnt/lustre/users/k1632479/grcm38/Ensembl/Mus_musculus.GRCm38.95.gtf


#for polyester I am going to load 
#install polyester and module load R
module load general/R/3.5.0
R

#source("http://bioconductor.org/biocLite.R")
#biocLite("polyester")
library("polyester")
library("Biostrings")
library("stringr")


EnsemblfastaFile <- "/mnt/lustre/users/k1632479/polyester/Mus_musculus.GRCm38_transcripts.fa"
fastaFile <- readDNAStringSet(EnsemblfastaFile)
fastapath = "/mnt/lustre/users/k1632479/polyester/Mus_musculus.GRCm38_transcripts.fa"

#count the number of transcripts found in the gff fasta file
numtx = count_transcripts(fastapath)
tNames = word(names(readDNAStringSet(fastapath)))


set.seed(4)

fold_change_values = sample(c(0.5, 1, 2), size=2*numtx, prob=c(0.05, 0.9, 0.05), replace=TRUE)
fold_changes = matrix(fold_change_values, nrow=numtx, ncol=6)
#work out dimensions
colnames(fold_changes) <-c("X13799X1_1", "X13799X1_2", "X13799X1_3", "X13799X2_1", "X13799X2_2", "X13799X2_3" )
rownames(fold_changes) <- tNames


#read.gtf(NMD_transcript, filename="/mnt/lustre/users/k1632479/grcm38/Ensembl/Mus_musculus.GRC38.95_only_NMD.gtf")
#NMD_transcriptnames <- as.data.frame((NMD_transcript@ev$gtf$transcript_id))
#colnames(NMD_transcriptnames) <-c("t_name")
NMD_transcriptnames <- read.table("/mnt/lustre/users/k1632479/polyester/NMD_transcriptnames", quote="\"", comment.char="")

mean(fold_changes)
#[1] 1.025513

for (i in c(1:numtx)) {
  if (any(str_detect(NMD_transcriptnames$t_name, as.vector(rownames(fold_changes))[i])) ==TRUE ){
  fold_changes[i,1] <- 0
fold_changes[i,2] <- 0 
fold_changes[i,3] <- 0 
fold_changes[i,4]<- 100 
fold_changes[i,5] <- 100 
fold_changes[i,6] <- 100 }}

write.table(fold_changes, file = "/mnt/lustre/users/k1632479/polyester/fold_changes.txt")

simulate_experiment(fastaFile, numreps=c(3,3), meanmodel=TRUE,fold_changes=fold_changes, outdir='/mnt/lustre/users/k1632479/polyester/', transcriptid=tNames, seed=12, error_model='illumina6', bias='rnaf')

#distr normal
#'rnaf' represents positional bias that arises in protocols using RNA fragmentation
#error model
#no gcbias as human cell lines
#frag_GC_bias not using as adds
#meanmol- sets read per transcripts as a function of transcription length using linear model 

