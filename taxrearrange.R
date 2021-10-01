# Taxonomy table rearrange
tax.rearrange <- function(taxonomytsv,file=TRUE,sep="; "){
  if(file==TRUE){
    fileinput <- read.delim(taxonomytsv)
    taxa <- as.character(fileinput[,2])
  } else if(file==FALSE) taxa <- as.character(taxonomytsv[,2])
  taxa <- strsplit(taxa,sep)
  for(i in 1:length(taxa)){
    if(length(taxa[[i]]) < 7){
      if(7-length(taxa[[i]]) ==1){
        taxa[[i]][7] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==2){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==3){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==4){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
        taxa[[i]][4] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==5){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
        taxa[[i]][4] <- "Unassigned"
        taxa[[i]][3] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==6){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
        taxa[[i]][4] <- "Unassigned"
        taxa[[i]][3] <- "Unassigned"
        taxa[[i]][2] <- "Unassigned"
      } 
    }
  }
  
  taxa <- as.data.frame(do.call(rbind,taxa))
  if(file==TRUE){
    taxa <- cbind(Feature.ID=fileinput[,1],taxa)
  } else if(file==FALSE) taxa <- cbind(Feature.ID=taxonomytsv[,1],taxa)
  taxa <- apply(taxa, 2, as.character)
  rownames(taxa) <- taxa[,1]
  taxa <- taxa[,-1]
  return(taxa)
}