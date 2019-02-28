# need to seperate twice with three repeats

#edit the pheno_data and the names
PHENO_DATA.csv
sample,cell_type,treatment
X13799X1_1,ESC,DMSO
X13799X1_2,ESC,DMSO
X13799X1_3,ESC,DMSO
X13799X2_1,ESC,CHX
X13799X2_2,ESC,CHX
X13799X2_3,ESC,CHX


#create ballgown object 

module load general/R/3.5.0

#load ballgown in R
#set the working directory

setwd("/mnt/lustre/users/k1632479/ESC_NSC/STAR/stringtie/ballgown/ballgown_input_files/")

library(ballgown)
#library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#need to run ballgown to get the object for running polyester

pheno_data = read.table(file ="PHENO_DATA.csv", header = TRUE, sep = ",")


#read in expresion data created by string tie
sample_full_path <- paste("/mnt/lustre/users/k1632479/ESC_NSC/STAR/stringtie/ballgown/ballgown_input_files",pheno_data[,1], sep = '/')
sample_full_path <- paste(sample_full_path,"", sep = '/')


bg = ballgown(samples=as.vector(sample_full_path),pData=pheno_data)


#Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)
#Load gene names for lookup later in the tutorial
bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])
#Pull the gene_expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg_filt))
#View the first five rows of data (all columns) in one of the dataframes created
head(gene_expression)

save(bg, file="bg.Rda")

#move to Rproject on my desktop
#with this file preloaded Mus_musculus.GRC38.95_only_NMD.gt


# need to find the way to summarise the FPKM values before change and after
differential_expression_test <- stattest(bg, feature = "transcript", covariate = "cell_type", adjustvars = c("treatment"), getFC=TRUE, meas = "FPKM")
FPKM_qval_mean <- mean( $qval, na.rm = TRUE)



#creating a list of transcripts that are NMD related- and the all have got a LFC change of 0  and 2 respectively 


install.packages("refGenome")
library(refGenome)

#reading in the GTF containing the NMD transcripts
NMD_transcript <- ensemblGenome()
read.gtf(NMD_transcript, filename="Mus_musculus.GRC38.95_only_NMD.gtf")

NMD_transcriptnames <- as.data.frame(matrix(0 ,nrow=length(c(NMD_transcript@ev$gtf$transcript_id)), ncol=6),row.names=NMD_transcript@ev$gtf$transcript_id)
colnames(NMD_transcriptnames) <-c("X13799X1_1", "X13799X1_2", "X13799X1_3", "X13799X2_1", "X13799X2_2", "X13799X2_3" )
NMD_transcriptnames[,4:6] <- 2


#edit the FPKM values using the NMD_transcript_names loop
#see below NEED to edit


selected_transcript<- as.data.frame(cbind( transcript$X13799X1, transcript$X13799X1, transcript$X13799X1, transcript$X13799X2, transcript$X13799X2, transcript$X13799X2 ), row.names=transcript$transcript_id, stringAsFactors=FALSE) 


#using the mean LFC of DESeq2
#<!-- > mean(abs(na.omit(resLFC$log2FoldChange))) -->
#<!-- [1] 0.1966057 -->

#apply only to 2,3,5,6 columns random number variety 
multi.fun <- function(x) {  
y <- runif(1, 0.95, 1.05)
round(x*y, digits=0) }

selected_transcript <- data.frame(selected_transcript[,1], lapply(selected_transcript[,2:3], multi.fun), selected_transcript[,4], lapply(selected_transcript[,5:6], multi.fun))

#merge the NMD_transcripts

library("stringr")

for (i in c(2:length(row.names(selected_transcript)))){
if (any(str_detect(row.names(NMD_transcriptnames), row.names(selected_transcript[i,])) == TRUE))
{selected_transcript[i,] <- NMD_transcriptnames["row.names(selected_transcript[i,])",]}
}





head_selected_transcript <- head(selected_transcript)
head_NMD_transcriptnames <- head(NMD_transcriptnames)
head_selected_transcript[2,] <- head_NMD_transcriptnames["row.names(selected_transcript[2,])",]


#the desired output is selected_transcript


differential_expression_test <- stattest(bg, feature = "transcript", covariate = "cell_type", adjustvars = c("treatment"), getFC=TRUE, meas = "FPKM")
difference <- FPKM_qval_mean-mean( $qval, na.rm = TRUE)
#check change 
