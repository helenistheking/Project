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

fold_change_values = sample(c(1, 2, 3), size=2*numtx, prob=c(0.05, 0.9, 0.05), replace=TRUE)
fold_changes = matrix(fold_change_values, nrow=numtx, ncol=3)
#work out dimensions
rownames(fold_changes) <- tNames


#read.gtf(NMD_transcript, filename="/mnt/lustre/users/k1632479/grcm38/Ensembl/Mus_musculus.GRC38.95_only_NMD.gtf")
#NMD_transcriptnames <- as.data.frame((NMD_transcript@ev$gtf$transcript_id))
#colnames(NMD_transcriptnames) <-c("t_name")
NMD_transcriptnames <- read.table("/mnt/lustre/users/k1632479/polyester/NMD_transcriptnames", quote="\"", comment.char="")

mean(fold_changes)
#[1] 1.025513

for (i in c(1:numtx)) {
  if (any(str_detect(NMD_transcriptnames$t_name, as.vector(rownames(fold_changes))[i])) ==TRUE ){
  fold_changes[i,1] <- 1
fold_changes[i,2]<- 100 }}

mean(fold_changes)

write.table(fold_changes, "/mnt/lustre/users/k1632479/polyester/fold_changes.txt")

colnames(fold_changes) <- NULL
rownames(fold_changes) <- NULL

simulate_experiment(EnsemblfastaFile, numreps=c(3,3),fold_changes= fold_changes, readlen=100, reads_per_transcript=readspertx, outdir="/mnt/lustre/users/k1632479/polyester/", distr="empirical", error_model="illumina5", bias="rnaf")

##output of results 

#[1] 2.000585
#[1] 3.610394

#to run the simulate_experiment call when the fold changes has already been created






library("polyester")
library("Biostrings")
library("stringr")

#has to be a simplified two column matrix
fold_changes<- data.matrix(read.table("/mnt/lustre/users/k1632479/polyester/fold_changes.txt")[,1:2])
rownames(fold_changes) <- NULL

EnsemblfastaFile <- "/mnt/lustre/users/k1632479/polyester/Mus_musculus.GRCm38_transcripts.fa"
fasta_File <- readDNAStringSet(EnsemblfastaFile)

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM

#replace all 0 counts with 1
fastaFile_nozero <- replace(width(fasta_File), width(fasta_File) == 0, 1)

readspertx = round(1 * fastaFile_nozero / 100)

readspertx <- replace(readspertx, readspertx == 0, 1)


simulate_experiment(EnsemblfastaFile, numreps=c(3,3),fold_changes=fold_changes, reads_per_transcript=readspertx, paired=FALSE, outdir="/mnt/lustre/users/k1632479/polyester/simulatedread", distr="empirical", error_model="illumina5", bias="rnaf")


####polyester cannot handle the size of the fasta file ect.
##so doing one repeat of both and trying to run
#if this does not work then im going to try to seperate the log fold changes



#distr normal
#'rnaf' represents positional bias that arises in protocols using RNA fragmentation
#error model
#no gcbias as human cell lines
#frag_GC_bias not using as adds
#meanmol- sets read per transcripts as a function of transcription length using linear model 

#remove the sequencing model 


###after error 

> lengths <- sapply(fastaFile, FUN=length)
> min(lengths)
[1] 9
> w <- which(lengths > 500)
> fastaFile_new <- fastaFile[w]
> simplified_fold_changes_new <- simplified_fold_changes[w]
> simplified_fold_changes_new <- simplified_fold_changes[w,]
> readspertx_new <- readspertx[w]
> writeXStringSet(fastaFile_new, 'fastaFile_new.fa')
tNames_new = word(names(readDNAStringSet('fastaFile_new.fa')))
simulate_experiment('fastaFile_new.fa', numreps=c(3,3),fold_changes= simplified_fold_changes_new,readlen=100, reads_per_transcript=readspertx_new, outdir="/mnt/lustre/users/k1632479/polyester/", transcriptid=tNames_new, )


