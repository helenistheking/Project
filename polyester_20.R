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
fold_changes_20 = matrix(fold_change_values, nrow=numtx, ncol=3)
#work out dimensions
rownames(fold_changes_20) <- tNames


#read.gtf(NMD_transcript, filename="/mnt/lustre/users/k1632479/grcm38/Ensembl/Mus_musculus.GRC38.95_only_NMD.gtf")
#NMD_transcriptnames <- as.data.frame((NMD_transcript@ev$gtf$transcript_id))
#colnames(NMD_transcriptnames) <-c("t_name")
NMD_transcriptnames <- read.table("/mnt/lustre/users/k1632479/polyester/NMD_transcriptnames", quote="\"", comment.char="")

mean(fold_changes)
#[1] 1.025513

for (i in c(1:numtx)) {
  if (any(str_detect(NMD_transcriptnames$t_name, as.vector(rownames(fold_changes_20))[i])) ==TRUE ){
  fold_changes_20[i,1] <- 1
fold_changes_20[i,2]<- 20 }}

mean(fold_changes)

write.table(fold_changes, "/mnt/lustre/users/k1632479/polyester/fold_changes_20.txt")

colnames(fold_changes_20) <- NULL
rownames(fold_changes_20) <- NULL



#has to be a simplified two column matrix
fold_changes_20<- data.matrix(read.table("/mnt/lustre/users/k1632479/polyester/fold_changes_20.txt")[,1:2])
rownames(fold_changes_20) <- NULL

EnsemblfastaFile <- "/mnt/lustre/users/k1632479/polyester/Mus_musculus.GRCm38_transcripts.fa"
fasta_File <- readDNAStringSet(EnsemblfastaFile)

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM

#replace all 0 counts with 1
fastaFile_nozero <- replace(width(fasta_File), width(fasta_File) == 0, 1)

readspertx = round(20 * fastaFile_nozero / 100)



simulate_experiment(EnsemblfastaFile, numreps=c(1,1), paired=FALSE, fold_changes= fold_changes_20, reads_per_transcript=readspertx, outdir="/mnt/lustre/users/k1632479/polyester/simulated_reads_20", distr="empirical", error_model="illumina5", bias="rnaf")
