#move to Rproject on my desktop
#with this file preloaded Mus_musculus.GRC38.95_only_NMD.gt

library(ballgown)
#library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)


# need to find the way to summarise the FPKM values before change and after
load("bg.Rda")
load("bg_original.Rda")
pData(bg)
differential_expression_test <- stattest(bg, feature = "transcript", covariate = "treatment", meas = "FPKM", getFC = TRUE)

#outputs computationally singular so cannot do differential expression analysis



#creating a list of transcripts that are NMD related- and the all have got a LFC change of 0  and 2 respectively 
#install.packages("refGenome")
# library(refGenome)

#reading in the GTF containing the NMD transcripts
# NMD_transcript <- ensemblGenome()
# read.gtf(NMD_transcript, filename="Mus_musculus.GRC38.95_only_NMD.gtf")
# 
# NMD_transcriptnames <- as.data.frame(matrix(0 ,nrow=length(c(NMD_transcript@ev$gtf$transcript_id)), ncol=7))
# colnames(NMD_transcriptnames) <-c("t_name", "X13799X1_1", "X13799X1_2", "X13799X1_3", "X13799X2_1", "X13799X2_2", "X13799X2_3" )
# NMD_transcriptnames$t_name <- NMD_transcript@ev$gtf$transcript_id
# NMD_transcriptnames[,5:7] <- 100


#edit the FPKM values using the NMD_transcript_names loop
#see below NEED to edit


# selected_transcript<- as.data.frame(cbind( transcript$X13799X1, transcript$X13799X1, transcript$X13799X1, transcript$X13799X2, transcript$X13799X2, transcript$X13799X2 ), row.names=transcript$transcript_id, stringAsFactors=FALSE) 


#using the mean LFC of DESeq2
#<!-- > mean(abs(na.omit(resLFC$log2FoldChange))) -->
#<!-- [1] 0.1966057 -->

#apply only to 2,3,5,6 columns random number variety 
multi.fun <- function(x) {  
  y <- runif(1, 0.90, 1.10)
  round(x*y, digits=6) }
#had to change run if range from 0.95 to 0.90 and 1.05 to 1.10 so that there is enough variety to run stattest

#looked at the code in github- polyester only uses FPKM values from the transcript section of  bg object
bg@expr[["trans"]][["FPKM.X13799X1_2"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X1_2"]], multi.fun)
bg@expr[["trans"]][["FPKM.X13799X1_3"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X1_3"]], multi.fun)
bg@expr[["trans"]][["FPKM.X13799X2_2"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X2_2"]], multi.fun)
bg@expr[["trans"]][["FPKM.X13799X2_3"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X2_3"]], multi.fun)

differential_expression_test <- stattest(bg1, feature = "transcript", covariate = "treatment", meas = "FPKM", getFC = TRUE)
FPKM_qval_mean <- mean(differential_expression_test$qval, na.rm = TRUE)
FPKM_fc_mean <- mean(differential_expression_test$fc, na.rm = TRUE)


library("stringr")

#log(100)=2 (base 10) so that why change FPKM values!

for (i in c(1:length(bg@expr[["trans"]][["FPKM.X13799X1_1"]]))) {
  if (any(str_detect(NMD_transcriptnames$t_name, bg@expr[["trans"]][["t_name"]][i])) == TRUE){
  bg@expr[["trans"]][["FPKM.X13799X1_1"]][i] <- 0 
bg@expr[["trans"]][["FPKM.X13799X1_2"]][i] <- 0 
bg@expr[["trans"]][["FPKM.X13799X1_3"]][i] <- 0 
bg@expr[["trans"]][["FPKM.X13799X2_1"]][i] <- 100 
bg@expr[["trans"]][["FPKM.X13799X2_2"]][i] <- 100 
bg@expr[["trans"]][["FPKM.X13799X2_3"]][i] <- 100 }
}

bg@expr[["trans"]][["FPKM.X13799X1_2"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X1_2"]], multi.fun)
bg@expr[["trans"]][["FPKM.X13799X1_3"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X1_3"]], multi.fun)
bg@expr[["trans"]][["FPKM.X13799X2_2"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X2_2"]], multi.fun)
bg@expr[["trans"]][["FPKM.X13799X2_3"]] <- sapply(bg@expr[["trans"]][["FPKM.X13799X2_3"]], multi.fun)

#the desired output is selected_transcript

differential_expression_test_2 <- stattest(bg, feature = "transcript", covariate = "treatment", meas = "FPKM", getFC = TRUE)
FPKM_qval_mean_2 <- mean(differential_expression_test_2$qval, na.rm = TRUE)
FPKM_fc_mean_2 <- mean(differential_expression_test_2$fc, na.rm = TRUE)

# FPKM_fc_mean
# 11143382226
# FPKM_fc_mean_2
# 2.967385e+18
# FPKM_qval_mean
# 0.9998665
# FPKM_qval_mean_2
# 0.9998592
