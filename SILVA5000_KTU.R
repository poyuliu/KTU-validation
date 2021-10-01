##### README #####
# Download testing dataset from https://doi.org/10.5281/zenodo.5535693 or https://github.com/poyuliu/KTU-validation
# 1. 5000-sequences SILVA 132 database subset: SILVA5000_DNA.fasta
# 2. 7-Levels taxonomic information of subset sequences: SILVA5000_tax.txt
# Save these two files to the local or working directory
#
# Download an R function from https://github.com/poyuliu/KTU-validation 
# 1. An R function for rearranging taxonomy format: taxrearrange.R
# Save these two files to the local or working directory
##### README #####

#install package from github:
library(devtools)
install_github("poyuliu/KTU")

#load KTU package
library(KTU)

# run KTU clustering
kluster <- klustering(repseq = "SILVA5000_DNA.fasta",write.fasta = TRUE,cores = 4) # load "SILVA5000_DNA.fasta" from your local or working directory
saveRDS(kluster,file="kluster.RDS") # output the clustering result to your local ot working directory 

# evaluate KTU clustering: (1) number of clusters, (2) ANI% of each KTU, and (3) cosine divergence of each KTU
evaluation <- KTU::KTUsim.eval(klusterRDS = "kluster.RDS",ASVfasta = "SILVA5000_DNA.fasta") # load "kluster.RDS" and "SILVA5000_DNA.fasta" from your local or working directory
saveRDS(evaluation,file="ktu_eval.RDS") # output the evaluation result to your local ot working directory 

# evaluate the reproducibility of the KTU algorithm by repeating KTU clustering 
kluster2 <- klustering(repseq = "~/ManuScript_KTU/data/KTU_revise2/SILVA_validation/SILVA5000_DNA.fasta",write.fasta = TRUE,cores = 4)
saveRDS(kluster2,file="kluster2.RDS")
kluster3 <- klustering(repseq = "~/ManuScript_KTU/data/KTU_revise2/SILVA_validation/SILVA5000_DNA.fasta",write.fasta = TRUE,cores = 4)
saveRDS(kluster3,file="kluster3.RDS")
kluster4 <- klustering(repseq = "~/ManuScript_KTU/data/KTU_revise2/SILVA_validation/SILVA5000_DNA.fasta",write.fasta = TRUE,cores = 4)
saveRDS(kluster4,file="kluster4.RDS")
kluster5 <- klustering(repseq = "~/ManuScript_KTU/data/KTU_revise2/SILVA_validation/SILVA5000_DNA.fasta",write.fasta = TRUE,cores = 4)
saveRDS(kluster5,file="kluster5.RDS")

# annotate the taxonomic information
source("taxrearrange.R") # load an R function for rearranging taxonomy format
silvatax <- read.delim("SILVA5000_tax.txt",header = F) # load "SILVA5000_tax.txt" from your local or working directory
length(unique(silvatax$V2)) # check how many unique taxonomic names (2673)
silvatax <- tax.rearrange(silvatax,file = F,sep = ";")  # rearranging semicolon-separated taxonomy format to 7 columns

# sort table by cluster numbers
silvatax <- as.data.frame(silvatax[match(names(kluster$clusters),rownames(silvatax)),])
silvatax$kluster <- kluster$clusters
silvatax <- silvatax[order(silvatax$kluster),]
names(silvatax)[1:7] <-c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")

# evaluate taxonomy consistency of each taxonomic level
consensus=c(5:10)/10 # 6 criteria (≧50%, ≧60%, ≧70%, ≧80%, ≧90%, and 100% consensus candidate nomenclatures within a KTU) for each taxonomy level 
consensus.7lv <- list()
for(z in 1:6){
  xxx <- matrix(data = NA,ncol = 7,nrow = max(silvatax$kluster))
  
  for(i in 1:max(silvatax$kluster)){
    xx <- silvatax[which(silvatax$kluster==i), ]
    xx[, -8] <- apply(xx[, -8], 2, as.character)
    xy <- apply(xx, 2, function(x) names(table(x))[which(prop.table(table(x)) >= consensus[z])])
    xxx[i, 1:7] <- c(xy[[1]][1], xy[[2]][1], xy[[3]][1],xy[[4]][1], xy[[5]][1], xy[[6]][1], xy[[7]][1])
    xxx[i, 7] <- ifelse(!is.na(xxx[i, 7]) & is.na(xxx[i,5]), NA, xxx[i, 7])
    for (j in 7:2) xxx[i, j] <- ifelse(!is.na(xxx[i, j]) & is.na(xxx[i, j - 1]), NA, xxx[i, j])
  }
  
  consensus.7lv[[z]] <- round((colSums(!is.na(xxx))/max(silvatax$kluster))*100,2) # summarize how many KTUs have consistent taxonomy under a consensus criterium
}

consensus.7lv <- as.data.frame(do.call(rbind,consensus.7lv))
rownames(consensus.7lv) <- consensus
colnames(consensus.7lv) <- c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")

# output results to local or working directory
write.csv(silvatax,"silvatax.csv")
write.csv(consensus.7lv,"consensus7lv.csv")

