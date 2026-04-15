##########################################################################################
### 002 - LER AS POSICOES INICIAIS E FINAIS DO GENE SPIKE EM CADA SEQUENCIA DA AMOSTRA ### 
### E RETORNAR ARQUIVOS FASTA SOMENTE COM PARTE DA SEQUENCIA REFERENTE A ELE           ###
##########################################################################################

if (!require("seqinr", quietly = TRUE))
  install.packages("seqinr")
if (!require("stringr", quietly = TRUE))
  install.packages("stringr")
if (!require("readxl", quietly = TRUE))
  install.packages("readxl")

library(seqinr)
library(stringr)
library(readxl)

spike <- as.matrix(read_excel("amostra.xls"))[,c(1,6:8)]

for (i in 1:n) {
  write.fasta(sequences = str_sub(toupper(delta_fasta[which(delta$accession==as.character(spike[i,1]))]),
    spike[i,2],spike[i,3]), names = paste0(str_sub(getAnnot(delta_fasta[which(delta$accession==as.character(
    spike[i,1]))]), 1,str_length(getAnnot(delta_fasta[which(delta$accession==as.character(spike[i,1]))]))-5),spike[i,4]), 
    file.out = paste0("spike/spike",i,".fasta"))
  write.fasta(sequences = str_sub(toupper(omicron_fasta[which(omicron$accession==as.character(spike[i+n,1]))]),
    spike[i+n,2],spike[i+n,3]), names = paste0(str_sub(getAnnot(omicron_fasta[which(omicron$accession==as.character(
    spike[i+n,1]))]), 1,str_length(getAnnot(omicron_fasta[which(omicron$accession==as.character(spike[i+n,1]))]))-5),
    spike[i+n,4]), file.out = paste0("spike/spike",i+n,".fasta"))
}