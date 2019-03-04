Candidate_genes_ESC <- read_excel("~/Desktop/Saqi Lab/Candidate_genes_ESC.xls", range = "A1:C65536")
library(dplyr)
library(readxl)

#due to just being a one repeat, it is almost irrelevant

colnames(Candidate_genes_ESC)<- c("gene_id","ESC_DMSO","ESC_CHX")
Candidate_genes_ESC <- Candidate_genes_ESC %>% mutate( LogFC = log10(( ESC_CHX - ESC_DMSO ) /ESC_DMSO ))
Candidate_genes_ESC <- Candidate_genes_ESC %>% filter( LogFC >=0 )
Candidate_genes_ESC <- Candidate_genes_ESC %>% filter( Inf>LogFC )



Candidate_genes_NMD_list <- read.csv("~/Desktop/Saqi Lab/Candidate_genes_ESC.csv", header=FALSE)


Candidate_genes_NMD_list<- select(Candidate_genes_NMD_list, V1)
colnames(Candidate_genes_NMD_list) <- c("gene_id")
final_Candidate_genes_ESC <- merge.data.frame(Candidate_genes_NMD_list, Candidate_genes_ESC, by="gene_id" )
final_Candidate_genes_ESC <- final_Candidate_genes_ESC %>% arrange(desc(LogFC))
final_Candidate_genes_ESC <- unique(final_Candidate_genes_ESC)                                                                   

write.csv(final_Candidate_genes_ESC, file="candidate_genes_ESC_upreg_NMD.csv")
