##########################################################################################
### 003 - LER OS RESULTADOS DOS ALINHAMENTOS DAS SEQUENCIAS E APRTIR DAI SEGUIR COM AS ###
### ETAPAS DE ANALISE EXPLORATRIA DOS DADOS E INFERENCIA                               ### 
##########################################################################################

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
if (!require("dendextend", quietly = TRUE))
  install.packages("dendextend")

library(dplyr)
library(ggplot2)
library(pheatmap)
library(dendextend)

alinhamento <- read.csv("alinhamento.csv", sep=";")
alinhamento <- cbind(alinhamento,completa_pdist=alinhamento$completa_mismatches/(alinhamento$completa_matches+
                                                                                   alinhamento$completa_mismatches+alinhamento$completa_gap_openings+alinhamento$completa_gap_extentions))
alinhamento <- cbind(alinhamento,completa_jc69=-3/4*log(1-4/3*alinhamento$completa_pdist))
alinhamento <- cbind(alinhamento,completa_k2p=-1/2*log(1-2*((alinhamento$completa_mismatches-alinhamento$completa_gap_openings)/
                                                              (alinhamento$completa_matches+alinhamento$completa_mismatches+alinhamento$completa_gap_openings+alinhamento$completa_gap_extentions))-
                                                         alinhamento$completa_gap_openings/(alinhamento$completa_matches+alinhamento$completa_mismatches+alinhamento$completa_gap_openings+
                                                                                              alinhamento$completa_gap_extentions))-1/4*log(1-2*alinhamento$completa_gap_openings/(alinhamento$completa_matches+
                                                                                                                                                                                     alinhamento$completa_mismatches+alinhamento$completa_gap_openings+alinhamento$completa_gap_extentions)))
alinhamento <- cbind(alinhamento,spike_pdist=alinhamento$spike_mismatches/(alinhamento$spike_matches+
                                                                             alinhamento$spike_mismatches+alinhamento$spike_gap_openings+alinhamento$spike_gap_extentions))
alinhamento <- cbind(alinhamento,spike_jc69=-3/4*log(1-4/3*alinhamento$spike_pdist))
alinhamento <- cbind(alinhamento,spike_k2p=-1/2*log(1-2*((alinhamento$spike_mismatches-alinhamento$spike_gap_openings)/
                                                           (alinhamento$spike_matches+alinhamento$spike_mismatches+alinhamento$spike_gap_openings+alinhamento$spike_gap_extentions))-
                                                      alinhamento$spike_gap_openings/(alinhamento$spike_matches+alinhamento$spike_mismatches+alinhamento$spike_gap_openings+
                                                                                        alinhamento$spike_gap_extentions))-1/4*log(1-2*alinhamento$spike_gap_openings/(alinhamento$spike_matches+
                                                                                                                                                                         alinhamento$spike_mismatches+alinhamento$spike_gap_openings+alinhamento$spike_gap_extentions)))
alinhamento <- cbind(alinhamento,tempo2=abs(alinhamento$seq1_data_coleta-alinhamento$seq2_data_coleta))

# Histogramas por grupo
grupo1 <- alinhamento[which(alinhamento$grupo=="inter-variantes"),]
grupo2 <- alinhamento[which(alinhamento$grupo=="intra-delta"),]
grupo3 <- alinhamento[which(alinhamento$grupo=="intra-omicron"),]

# série completa
par(mfrow=c(3,3), oma = c(2, 0, 3, 0))

hist(grupo1$completa_pdist,probability=TRUE,col="lightblue",main="Inter-Variantes",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo1$completa_pdist),sd=sd(grupo1$completa_pdist))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo2$completa_pdist,probability=TRUE,col="lightgreen",main="Intra-Delta",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo2$completa_pdist),sd=sd(grupo2$completa_pdist))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo3$completa_pdist,probability=TRUE,col="lightyellow",main="Intra-Omicron",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo3$completa_pdist),sd=sd(grupo3$completa_pdist))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo1$completa_jc69,probability=TRUE,col="lightblue",main="Inter-Variantes",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo1$completa_jc69),sd=sd(grupo1$completa_jc69))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo2$completa_jc69,probability=TRUE,col="lightgreen",main="Intra-Delta",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo2$completa_jc69),sd=sd(grupo2$completa_jc69))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo3$completa_jc69,probability=TRUE,col="lightyellow",main="Intra-Omicron",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo3$completa_jc69),sd=sd(grupo3$completa_jc69))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo1$completa_k2p,probability=TRUE,col="lightblue",main="Inter-Variantes",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo1$completa_k2p),sd=sd(grupo1$completa_k2p))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo2$completa_k2p,probability=TRUE,col="lightgreen",main="Intra-Delta",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo2$completa_k2p),sd=sd(grupo2$completa_k2p))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo3$completa_k2p,probability=TRUE,col="lightyellow",main="Intra-Omicron",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo3$completa_k2p),sd=sd(grupo3$completa_k2p))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

mtext("Série Completa\n p-distance", side = 3, outer = TRUE, cex = 1, font = 2)
mtext("Jukes-Cantor\n\n\n\n\n\n\n\n\n\n\n\n\n\n Kimura 2-Parameter", side = 1, outer = TRUE, cex = 1, font = 2, line = -23)

# gene spike
par(mfrow=c(3,3), oma = c(2, 0, 3, 0))

hist(grupo1$spike_pdist,probability=TRUE,col="lightblue",main="Inter-Variantes",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo1$spike_pdist),sd=sd(grupo1$spike_pdist))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo2$spike_pdist,probability=TRUE,col="lightgreen",main="Intra-Delta",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo2$spike_pdist),sd=sd(grupo2$spike_pdist))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo3$spike_pdist,probability=TRUE,col="lightyellow",main="Intra-Omicron",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo3$spike_pdist),sd=sd(grupo3$spike_pdist))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo1$spike_jc69,probability=TRUE,col="lightblue",main="Inter-Variantes",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo1$spike_jc69),sd=sd(grupo1$spike_jc69))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo2$spike_jc69,probability=TRUE,col="lightgreen",main="Intra-Delta",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo2$spike_jc69),sd=sd(grupo2$spike_jc69))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo3$spike_jc69,probability=TRUE,col="lightyellow",main="Intra-Omicron",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo3$spike_jc69),sd=sd(grupo3$spike_jc69))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo1$spike_k2p,probability=TRUE,col="lightblue",main="Inter-Variantes",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo1$spike_k2p),sd=sd(grupo1$spike_k2p))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo2$spike_k2p,probability=TRUE,col="lightgreen",main="Intra-Delta",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo2$spike_k2p),sd=sd(grupo2$spike_k2p))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

hist(grupo3$spike_k2p,probability=TRUE,col="lightyellow",main="Intra-Omicron",xlab="",ylab="")
dados <- rnorm(1000,mean=mean(grupo3$spike_k2p),sd=sd(grupo3$spike_k2p))
x_valores <- seq(min(dados),max(dados),length=100)
y_valores <- dnorm(x_valores,mean=mean(dados),sd=sd(dados))
lines(x_valores,y_valores,col="red",lwd=2)

mtext("Gene Spike\n p-distance", side = 3, outer = TRUE, cex = 1, font = 2)
mtext("Jukes-Cantor\n\n\n\n\n\n\n\n\n\n\n\n\n\n Kimura 2-Parameter", side = 1, outer = TRUE, cex = 1, font = 2, line = -23)


# Q-Q Plot por grupo

# série completa
par(mfrow=c(3,3), oma = c(2, 0, 3, 0))

qqnorm(grupo1$completa_pdist,main="Inter-Variantes",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo1$completa_pdist,"pnorm",mean(grupo1$completa_pdist),sd(grupo1$completa_pdist))$p.value,6),")"),ylab="")
qqline(grupo1$completa_pdist,col="red",lwd=2)

qqnorm(grupo2$completa_pdist,main="Intra-Delta",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo2$completa_pdist,"pnorm",mean(grupo2$completa_pdist),sd(grupo2$completa_pdist))$p.value,6),")"),ylab="")
qqline(grupo2$completa_pdist,col="red",lwd=2)

qqnorm(grupo3$completa_pdist,main="Intra-Omicron",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo3$completa_pdist,"pnorm",mean(grupo3$completa_pdist),sd(grupo3$completa_pdist))$p.value,6),")"),ylab="")
qqline(grupo3$completa_pdist,col="red",lwd=2)

qqnorm(grupo1$completa_jc69,main="Inter-Variantes",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo1$completa_jc69,"pnorm",mean(grupo1$completa_jc69),sd(grupo1$completa_jc69))$p.value,6),")"),ylab="")
qqline(grupo1$completa_jc69,col="red",lwd=2)

qqnorm(grupo2$completa_jc69,main="Intra-Delta",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo2$completa_jc69,"pnorm",mean(grupo2$completa_jc69),sd(grupo2$completa_jc69))$p.value,6),")"),ylab="")
qqline(grupo2$completa_jc69,col="red",lwd=2)

qqnorm(grupo3$completa_jc69,main="Intra-Omicron",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo3$completa_jc69,"pnorm",mean(grupo3$completa_jc69),sd(grupo3$completa_jc69))$p.value,6),")"),ylab="")
qqline(grupo3$completa_jc69,col="red",lwd=2)

qqnorm(grupo1$completa_k2p,main="Inter-Variantes",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo1$completa_k2p,"pnorm",mean(grupo1$completa_k2p),sd(grupo1$completa_k2p))$p.value,6),")"),ylab="")
qqline(grupo1$completa_k2p,col="red",lwd=2)

qqnorm(grupo2$completa_k2p,main="Intra-Delta",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo2$completa_k2p,"pnorm",mean(grupo2$completa_k2p),sd(grupo2$completa_k2p))$p.value,6),")"),ylab="")
qqline(grupo2$completa_k2p,col="red",lwd=2)

qqnorm(grupo3$completa_k2p,main="Intra-Omicron",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo3$completa_k2p,"pnorm",mean(grupo3$completa_k2p),sd(grupo3$completa_k2p))$p.value,6),")"),ylab="")
qqline(grupo3$completa_k2p,col="red",lwd=2)

mtext("Série Completa\n p-distance", side = 3, outer = TRUE, cex = 1, font = 2)
mtext("Jukes-Cantor\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n Kimura 2-Parameter", side = 1, outer = TRUE, cex = 1, font = 2, line = -22)

# gene spike
par(mfrow=c(3,3), oma = c(2, 0, 3, 0))

qqnorm(grupo1$spike_pdist,main="Inter-Variantes",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo1$spike_pdist,"pnorm",mean(grupo1$spike_pdist),sd(grupo1$spike_pdist))$p.value,6),")"),ylab="")
qqline(grupo1$spike_pdist,col="red",lwd=2)

qqnorm(grupo2$spike_pdist,main="Intra-Delta",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo2$spike_pdist,"pnorm",mean(grupo2$spike_pdist),sd(grupo2$spike_pdist))$p.value,6),")"),ylab="")
qqline(grupo2$spike_pdist,col="red",lwd=2)

qqnorm(grupo3$spike_pdist,main="Intra-Omicron",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo3$spike_pdist,"pnorm",mean(grupo3$spike_pdist),sd(grupo3$spike_pdist))$p.value,6),")"),ylab="")
qqline(grupo3$spike_pdist,col="red",lwd=2)

qqnorm(grupo1$spike_jc69,main="Inter-Variantes",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo1$spike_jc69,"pnorm",mean(grupo1$spike_jc69),sd(grupo1$spike_jc69))$p.value,6),")"),ylab="")
qqline(grupo1$spike_jc69,col="red",lwd=2)

qqnorm(grupo2$spike_jc69,main="Intra-Delta",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo2$spike_jc69,"pnorm",mean(grupo2$spike_jc69),sd(grupo2$spike_jc69))$p.value,6),")"),ylab="")
qqline(grupo2$spike_jc69,col="red",lwd=2)

qqnorm(grupo3$spike_jc69,main="Intra-Omicron",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo3$spike_jc69,"pnorm",mean(grupo3$spike_jc69),sd(grupo3$spike_jc69))$p.value,6),")"),ylab="")
qqline(grupo3$spike_jc69,col="red",lwd=2)

qqnorm(grupo1$spike_k2p,main="Inter-Variantes",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo1$spike_k2p,"pnorm",mean(grupo1$spike_k2p),sd(grupo1$spike_k2p))$p.value,6),")"),ylab="")
qqline(grupo1$spike_k2p,col="red",lwd=2)

qqnorm(grupo2$spike_k2p,main="Intra-Delta",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo2$spike_k2p,"pnorm",mean(grupo2$spike_k2p),sd(grupo2$spike_k2p))$p.value,6),")"),ylab="")
qqline(grupo2$spike_k2p,col="red",lwd=2)

qqnorm(grupo3$spike_k2p,main="Intra-Omicron",xlab=paste0("teste KS\n (p-valor = ",round(ks.test(grupo3$spike_k2p,"pnorm",mean(grupo3$spike_k2p),sd(grupo3$spike_k2p))$p.value,6),")"),ylab="")
qqline(grupo3$spike_k2p,col="red",lwd=2)

mtext("Gene Spike\n p-distance", side = 3, outer = TRUE, cex = 1, font = 2)
mtext("Jukes-Cantor\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n Kimura 2-Parameter", side = 1, outer = TRUE, cex = 1, font = 2, line = -22)

#teste kruskal-Wallis
kruskal.test(completa_pdist ~ grupo, data = alinhamento)
kruskal.test(completa_jc69 ~ grupo, data = alinhamento)
kruskal.test(completa_k2p ~ grupo, data = alinhamento)
kruskal.test(spike_pdist ~ grupo, data = alinhamento)
kruskal.test(spike_jc69 ~ grupo, data = alinhamento)
kruskal.test(spike_k2p ~ grupo, data = alinhamento)

#box-plot
# série completa
par(mfrow=c(3,1), oma = c(2, 0, 3, 0))
boxplot(completa_pdist ~ grupo, data=alinhamento,main="p-distance",col = c("lightblue","lightgreen","lightyellow"),xlab="", ylab="")
boxplot(completa_jc69 ~ grupo, data=alinhamento,main="Jukes-Cantor",col = c("lightblue","lightgreen","lightyellow"),xlab="", ylab="")
boxplot(completa_k2p ~ grupo, data=alinhamento,main="Kimura 2-Parameter",col = c("lightblue","lightgreen","lightyellow"),xlab="", ylab="")
mtext("Série Completa", side = 3, outer = TRUE, cex = 1, font = 2)

# gene spike
par(mfrow=c(3,1), oma = c(2, 0, 3, 0))
boxplot(spike_pdist ~ grupo, data=alinhamento,main="p-distance",col = c("lightblue","lightgreen","lightyellow"),xlab="", ylab="")
boxplot(spike_jc69 ~ grupo, data=alinhamento,main="Jukes-Cantor",col = c("lightblue","lightgreen","lightyellow"),xlab="", ylab="")
boxplot(spike_k2p ~ grupo, data=alinhamento,main="Kimura 2-Parameter",col = c("lightblue","lightgreen","lightyellow"),xlab="", ylab="")
mtext("Gene Spike", side = 3, outer = TRUE, cex = 1, font = 2)

#heatmap
nomes<-c("OL757809.1","ON483559.1","OM961293.1","OL351442.1","ON909201.1","OM190651.1","OL958521.1",
         "OM692491.1","OM721708.1","OL966477.1","OK190628.1","MZ558098.1","OL958468.1","OM090161.1","ON507019.1",
         "OL989400.1","ON485812.1","OK190620.1","OM961262.1","OL989407.1","OL989392.1","OL351383.1","OL966486.1",
         "OK104624.1","OL958479.1","OK189643.1","OL533312.1","OK190314.1","OL929501.1","MZ996543.1","MZ979005.1",
         "OK396712.1","OL532798.1","OL416456.1","OL535485.1","MZ983858.1","OL776371.1","OL525066.1","OK600503.1",
         "MZ450170.1","OK040851.1","OL533077.1","OL527012.1","PP219952.1","OK654484.1","OP397223.1","OM244655.1",
         "OM246163.1","OM191603.1","OM524823.1","OK507690.1","OM689003.1","OM246247.1","OM153932.1","OK653873.1",
         "OL411226.1","PP142372.1","OM244445.1","PP142752.1","OM046899.1","OK385946.1","PP143202.1","OM524909.1",
         "OM245956.1","OM244619.1","PP299883.1","OL924833.1","OK652202.1","PP143434.1","PP220871.1","OK657000.1",
         "OP397072.1","OM244624.1","OM246028.1","OM244757.1","OM246049.1","OM244701.1","PP142907.1","OK383519.1",
         "OP397330.1","OL955122.1","OM244612.1","OM246166.1","OM246034.1","OK658366.1","OL541835.1","OM246133.1",
         "OM153467.1","OK158384.1","PP143706.1","OP397426.1","OK507706.1","OK652281.1","OK386497.1","OM246500.1",
         "OP397014.1","MZ726510.1","OM054471.1","OM244584.1","OM244636.1","PQ577963.1","PQ567389.1","PQ567392.1",
         "PQ567494.1","PQ281310.1","PQ281378.1","PQ281386.1","OR725981.1","OR725982.1","OR725983.1","OR717328.1",
         "OQ954650.1","OQ052586.1","OQ059020.1","OQ059024.1","OQ059027.1","OQ059028.1","OP183419.1","OP183425.1",
         "OP183428.1","OP183433.1","OP183435.1","OP183436.1","OP183440.1","OP183447.1","OP183454.1","OQ315917.1",
         "ON637175.1","ON384194.1","OM396818.1","OM332107.1","OM297839.1","OM263440.1","OM196622.1","OM098592.1",
         "OM080504.1","OM080566.1","OM080680.1","OM080705.1","OM080710.1","OM038373.1","OM038376.1","OL579980.1",
         "ON189634.1","ON313199.1","ON071365.1","OP879969.1","OP879987.1","OP419356.1","OP419299.1","ON230413.1",
         "ON313287.1","OP419274.1","ON313032.1","OP051818.1","ON071549.1","OM229419.1","ON305405.1","ON188712.1",
         "OP419315.1","OP419363.1","OP419272.1","OP419393.1","ON230809.1","ON189849.1","OP419385.1","OP419426.1",
         "OP419435.1","ON787303.1","ON118867.1","OP929937.1","OP983736.1","OM294238.1","ON313485.1","ON264105.1",
         "ON188719.1","OM390558.1","OP419449.1","ON189200.1","ON068774.1","ON966040.1","OP419335.1","ON305911.1",
         "ON003813.1","ON068785.1","ON188726.1","ON800468.1","ON313115.1","ON595182.1","OP419340.1","ON118588.1",
         "ON071109.1","ON188704.1","OP419295.1","OP419438.1","ON594963.1","OP880133.1","ON313238.1","ON208913.1",
         "ON226320.1")

heat<-data.frame(OL757809.1=numeric(200),ON483559.1=numeric(200),OM961293.1=numeric(200),
                 OL351442.1=numeric(200),ON909201.1=numeric(200),OM190651.1=numeric(200),OL958521.1=numeric(200),
                 OM692491.1=numeric(200),OM721708.1=numeric(200),OL966477.1=numeric(200),OK190628.1=numeric(200),
                 MZ558098.1=numeric(200),OL958468.1=numeric(200),OM090161.1=numeric(200),ON507019.1=numeric(200),
                 OL989400.1=numeric(200),ON485812.1=numeric(200),OK190620.1=numeric(200),OM961262.1=numeric(200),
                 OL989407.1=numeric(200),OL989392.1=numeric(200),OL351383.1=numeric(200),OL966486.1=numeric(200),
                 OK104624.1=numeric(200),OL958479.1=numeric(200),OK189643.1=numeric(200),OL533312.1=numeric(200),
                 OK190314.1=numeric(200),OL929501.1=numeric(200),MZ996543.1=numeric(200),MZ979005.1=numeric(200),
                 OK396712.1=numeric(200),OL532798.1=numeric(200),OL416456.1=numeric(200),OL535485.1=numeric(200),
                 MZ983858.1=numeric(200),OL776371.1=numeric(200),OL525066.1=numeric(200),OK600503.1=numeric(200),
                 MZ450170.1=numeric(200),OK040851.1=numeric(200),OL533077.1=numeric(200),OL527012.1=numeric(200),
                 PP219952.1=numeric(200),OK654484.1=numeric(200),OP397223.1=numeric(200),OM244655.1=numeric(200),
                 OM246163.1=numeric(200),OM191603.1=numeric(200),OM524823.1=numeric(200),OK507690.1=numeric(200),
                 OM689003.1=numeric(200),OM246247.1=numeric(200),OM153932.1=numeric(200),OK653873.1=numeric(200),
                 OL411226.1=numeric(200),PP142372.1=numeric(200),OM244445.1=numeric(200),PP142752.1=numeric(200),
                 OM046899.1=numeric(200),OK385946.1=numeric(200),PP143202.1=numeric(200),OM524909.1=numeric(200),
                 OM245956.1=numeric(200),OM244619.1=numeric(200),PP299883.1=numeric(200),OL924833.1=numeric(200),
                 OK652202.1=numeric(200),PP143434.1=numeric(200),PP220871.1=numeric(200),OK657000.1=numeric(200),
                 OP397072.1=numeric(200),OM244624.1=numeric(200),OM246028.1=numeric(200),OM244757.1=numeric(200),
                 OM246049.1=numeric(200),OM244701.1=numeric(200),PP142907.1=numeric(200),OK383519.1=numeric(200),
                 OP397330.1=numeric(200),OL955122.1=numeric(200),OM244612.1=numeric(200),OM246166.1=numeric(200),
                 OM246034.1=numeric(200),OK658366.1=numeric(200),OL541835.1=numeric(200),OM246133.1=numeric(200),
                 OM153467.1=numeric(200),OK158384.1=numeric(200),PP143706.1=numeric(200),OP397426.1=numeric(200),
                 OK507706.1=numeric(200),OK652281.1=numeric(200),OK386497.1=numeric(200),OM246500.1=numeric(200),
                 OP397014.1=numeric(200),MZ726510.1=numeric(200),OM054471.1=numeric(200),OM244584.1=numeric(200),
                 OM244636.1=numeric(200),PQ577963.1=numeric(200),PQ567389.1=numeric(200),PQ567392.1=numeric(200),
                 PQ567494.1=numeric(200),PQ281310.1=numeric(200),PQ281378.1=numeric(200),PQ281386.1=numeric(200),
                 OR725981.1=numeric(200),OR725982.1=numeric(200),OR725983.1=numeric(200),OR717328.1=numeric(200),
                 OQ954650.1=numeric(200),OQ052586.1=numeric(200),OQ059020.1=numeric(200),OQ059024.1=numeric(200),
                 OQ059027.1=numeric(200),OQ059028.1=numeric(200),OP183419.1=numeric(200),OP183425.1=numeric(200),
                 OP183428.1=numeric(200),OP183433.1=numeric(200),OP183435.1=numeric(200),OP183436.1=numeric(200),
                 OP183440.1=numeric(200),OP183447.1=numeric(200),OP183454.1=numeric(200),OQ315917.1=numeric(200),
                 ON637175.1=numeric(200),ON384194.1=numeric(200),OM396818.1=numeric(200),OM332107.1=numeric(200),
                 OM297839.1=numeric(200),OM263440.1=numeric(200),OM196622.1=numeric(200),OM098592.1=numeric(200),
                 OM080504.1=numeric(200),OM080566.1=numeric(200),OM080680.1=numeric(200),OM080705.1=numeric(200),
                 OM080710.1=numeric(200),OM038373.1=numeric(200),OM038376.1=numeric(200),OL579980.1=numeric(200),
                 ON189634.1=numeric(200),ON313199.1=numeric(200),ON071365.1=numeric(200),OP879969.1=numeric(200),
                 OP879987.1=numeric(200),OP419356.1=numeric(200),OP419299.1=numeric(200),ON230413.1=numeric(200),
                 ON313287.1=numeric(200),OP419274.1=numeric(200),ON313032.1=numeric(200),OP051818.1=numeric(200),
                 ON071549.1=numeric(200),OM229419.1=numeric(200),ON305405.1=numeric(200),ON188712.1=numeric(200),
                 OP419315.1=numeric(200),OP419363.1=numeric(200),OP419272.1=numeric(200),OP419393.1=numeric(200),
                 ON230809.1=numeric(200),ON189849.1=numeric(200),OP419385.1=numeric(200),OP419426.1=numeric(200),
                 OP419435.1=numeric(200),ON787303.1=numeric(200),ON118867.1=numeric(200),OP929937.1=numeric(200),
                 OP983736.1=numeric(200),OM294238.1=numeric(200),ON313485.1=numeric(200),ON264105.1=numeric(200),
                 ON188719.1=numeric(200),OM390558.1=numeric(200),OP419449.1=numeric(200),ON189200.1=numeric(200),
                 ON068774.1=numeric(200),ON966040.1=numeric(200),OP419335.1=numeric(200),ON305911.1=numeric(200),
                 ON003813.1=numeric(200),ON068785.1=numeric(200),ON188726.1=numeric(200),ON800468.1=numeric(200),
                 ON313115.1=numeric(200),ON595182.1=numeric(200),OP419340.1=numeric(200),ON118588.1=numeric(200),
                 ON071109.1=numeric(200),ON188704.1=numeric(200),OP419295.1=numeric(200),OP419438.1=numeric(200),
                 ON594963.1=numeric(200),OP880133.1=numeric(200),ON313238.1=numeric(200),ON208913.1=numeric(200),
                 ON226320.1=numeric(200),row.names=nomes)

for (i in 1:199) {
  for (j in (i+1):200) {
    heat[i,j] <- alinhamento$spike_k2p[which(alinhamento$seq1_accession==nomes[i] & alinhamento$seq2_accession==nomes[j])]   
  }
}

pheatmap(heat, cluster_rows = FALSE, cluster_cols = TRUE, main = "Heatmap Kimura 2-Parameter - Gene Spike")

#teste spearman
cor.test(grupo1$completa_pdist, grupo1$tempo2, method = "spearman")
cor.test(grupo1$completa_jc69, grupo1$tempo2, method = "spearman")
cor.test(grupo1$completa_k2p, grupo1$tempo2, method = "spearman")
cor.test(grupo1$spike_pdist, grupo1$tempo2, method = "spearman")
cor.test(grupo1$spike_jc69, grupo1$tempo2, method = "spearman")
cor.test(grupo1$spike_k2p, grupo1$tempo2, method = "spearman")

cor.test(grupo2$completa_pdist, grupo2$tempo2, method = "spearman")
cor.test(grupo2$completa_jc69, grupo2$tempo2, method = "spearman")
cor.test(grupo2$completa_k2p, grupo2$tempo2, method = "spearman")
cor.test(grupo2$spike_pdist, grupo2$tempo2, method = "spearman")
cor.test(grupo2$spike_jc69, grupo2$tempo2, method = "spearman")
cor.test(grupo2$spike_k2p, grupo2$tempo2, method = "spearman")

cor.test(grupo3$completa_pdist, grupo3$tempo2, method = "spearman")
cor.test(grupo3$completa_jc69, grupo3$tempo2, method = "spearman")
cor.test(grupo3$completa_k2p, grupo3$tempo2, method = "spearman")
cor.test(grupo3$spike_pdist, grupo3$tempo2, method = "spearman")
cor.test(grupo3$spike_jc69, grupo3$tempo2, method = "spearman")
cor.test(grupo3$spike_k2p, grupo3$tempo2, method = "spearman")

#data de coleta
boxplot(tempo2 ~ grupo, data=alinhamento,main="Tempo de coleta",col = c("lightblue","lightgreen","lightyellow"),xlab="", ylab="")

summary(alinhamento$tempo2[which(alinhamento$grupo=="inter-variantes")])
summary(alinhamento$tempo2[which(alinhamento$grupo=="intra-delta")])
summary(alinhamento$tempo2[which(alinhamento$grupo=="intra-omicron")])

#mismatches
summary(alinhamento$completa_mismatches[which(alinhamento$grupo=="inter-variantes")])
summary(alinhamento$completa_mismatches[which(alinhamento$grupo=="intra-delta")])
summary(alinhamento$completa_mismatches[which(alinhamento$grupo=="intra-omicron")])

summary(alinhamento$spike_mismatches[which(alinhamento$grupo=="inter-variantes")])
summary(alinhamento$spike_mismatches[which(alinhamento$grupo=="intra-delta")])
summary(alinhamento$spike_mismatches[which(alinhamento$grupo=="intra-omicron")])

#dendrograma
m = 200
matriz_distancia <- matrix(0, ncol=m, nrow=m)
vnome <- matrix("NA",ncol=m,nrow=1)
colnames(matriz_distancia) <- unique(c(unique(alinhamento$seq1_accession)[1:199],unique(alinhamento$seq2_accession)[1:199]))[c(1:(m/2),101:(100+m/2))]
rownames(matriz_distancia) <- unique(c(unique(alinhamento$seq1_accession)[1:199],unique(alinhamento$seq2_accession)[1:199]))[c(1:(m/2),101:(100+m/2))]
for (i in 1:m) {
  for (j in 1:m) {
    if (i>j)
      matriz_distancia[i,j] <- alinhamento$spike_k2p[which(alinhamento$seq1_accession==colnames(matriz_distancia)[j]&alinhamento$seq2_accession==rownames(matriz_distancia)[i])]
  }
}

for (i in 1:(m/2)) {
  vnome[i] <- paste0("D","")
  vnome[i+(m/2)] <- paste0("-","")
}

colnames(matriz_distancia) <- vnome
rownames(matriz_distancia) <- vnome
matriz_distancia <- dist(matriz_distancia)

hc <- hclust(matriz_distancia, method = "complete")
dend <- as.dendrogram(hc)

dend_final <- color_branches(dend, k = 2) %>%
  set("branches_lwd", 2) %>% # Ajustar a espessura das linhas
  set("labels_cex", 0.8) # Ajustar o tamanho dos rótulos

plot(dend_final,
     main = "Dendrograma Gene Spike",
     ylab = "Kimura 2-Parameter")