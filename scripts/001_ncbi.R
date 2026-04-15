################################################################
### 001 - LER AS SEQUENCIAS DO NCBI E SELECIONAR AS AMOSTRAS ###
################################################################

rm(list = ls())

if (!require("seqinr", quietly = TRUE))
  install.packages("seqinr")
if (!require("stringr", quietly = TRUE))
  install.packages("stringr")

library(seqinr)
library(stringr)

n <- 100  # Tamanho da amostra

delta_fasta <- read.fasta(file = "delta.fasta", seqtype = "DNA", as.string = TRUE, strip.desc = TRUE)

delta <- data.frame(
  accession = character(),
  variante = character(), 
  pais = character(),
  estado = character(),
  tamanho = integer(),
  stringsAsFactors = FALSE
)

for (i in 1:length(delta_fasta)) {
  delta <- rbind(delta, data.frame(
    accession = trimws(strsplit(getAnnot(delta_fasta[[i]]), "\\|")[[1]][1]),
    variante = strsplit(getAnnot(delta_fasta[[i]]), "\\|")[[1]][2],
    pais = strsplit(getAnnot(delta_fasta[[i]]), "\\|")[[1]][3],
    estado = strsplit(getAnnot(delta_fasta[[i]]), "\\|")[[1]][4],
    tamanho = strsplit(getAnnot(delta_fasta[[i]]), "\\|")[[1]][5]))
}

omicron_fasta <- read.fasta(file = "omicron.fasta", seqtype = "DNA", as.string = TRUE, strip.desc = TRUE)

omicron <- data.frame(
  accession = character(),
  variante = character(), 
  pais = character(),
  estado = character(),
  tamanho = integer(),
  stringsAsFactors = FALSE
)

for (i in 1:length(omicron_fasta)) {
  omicron <- rbind(omicron, data.frame(
    accession = trimws(strsplit(getAnnot(omicron_fasta[[i]]), "\\|")[[1]][1]),
    variante = strsplit(getAnnot(omicron_fasta[[i]]), "\\|")[[1]][2], 
    pais = strsplit(getAnnot(omicron_fasta[[i]]), "\\|")[[1]][3],
    estado = strsplit(getAnnot(omicron_fasta[[i]]), "\\|")[[1]][4],
    tamanho = strsplit(getAnnot(omicron_fasta[[i]]), "\\|")[[1]][5]))
}

#identificando quantas sequencias teremos na amostra com pais!='USA' e do 'USA' com estado informado
n_outros <- length(which(omicron$pais!='USA'))
n_usa <- length(which(omicron$estado!=''))

#selecionando as amostras
#delta
delta_n <- c(sample(which(delta$pais!='USA'), n_outros, replace = FALSE), sample(which(delta$estado!=''), n_usa, replace = FALSE), 
               sample(which(delta$pais=='USA' & delta$estado==''), n-n_outros-n_usa, replace = FALSE))

#omicron
omicron_n <- c(which(omicron$pais!='USA'), which(omicron$estado!=''), 
               sample(which(omicron$pais=='USA' & omicron$estado==''), n-n_outros-n_usa, replace = FALSE))

amostra <-rbind(delta[delta_n,], omicron[omicron_n,])
write.csv(amostra, "amostra.csv", row.names = FALSE)

for (i in 1:n) {
  write.fasta(sequences = toupper(delta_fasta[which(names(delta_fasta)==amostra[i,1])]), names = getAnnot(delta_fasta[which(names(delta_fasta)==amostra[i,1])]), file.out = paste0("seq/seq",i,".fasta"))
  write.fasta(sequences = toupper(omicron_fasta[which(names(omicron_fasta)==amostra[i+n,1])]), names = getAnnot(omicron_fasta[which(names(omicron_fasta)==amostra[i+n,1])]), file.out = paste0("seq/seq",i+n,".fasta"))
}