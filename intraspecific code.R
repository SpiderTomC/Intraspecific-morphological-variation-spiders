#initial loading of file
Multispid= read.csv("your file path")
intraspid=Multispid[,6:30]
intraspidlog=log(intraspid)
intraspiddist=dist(scale(intraspidlog))
intraspidpcoa=cmdscale(intraspiddist, k=10, eig=T)
plot(intraspidpcoa$points[,1],intraspidpcoa$points[,3],cex=0)
text(intraspidpcoa$points[,1],intraspidpcoa$points[,3],cex=0.6, labels= Multispid[,2])

#lda and probabilities
intralda<- lda(x=intraspidpcoa$points, grouping=Multispid[,2], CV=T)
intraspecies= c("Eriophora transmarina","Habronestes grahami","Latrodectus hasseltii","Nyssus albopunctatus", "Poltys noblei","Tasmanicosa godeffroyi","Triconephila plumipes","Venonia micarioides")
intraprob<-array()
for (i in 1:8){
  intraprob[i]<-mean(intralda$posterior[which(Multispid[,2]==intraspecies[i]),i],na.rm=T)
}

#size adjusted matrix + probabilities
adjintraspid= matrix(nrow = 145,ncol = 25)
for (i in 1:nrow(adjintraspid)) {
  rowratio<-mean(intraspid[,2],na.rm=T)/intraspid[i,2]
  for (j in 1:ncol(intraspid)){
    adjintraspid[i,j]<-intraspid[i,j]*rowratio
  }
}

adjintraspidlog=log(adjintraspid)
adjintraspiddist=dist(scale(adjintraspidlog))
adjintraspidpcoa=cmdscale(adjintraspiddist, k=10, eig=T)

adjintralda<- lda(x=adjintraspidpcoa$points, grouping=Multispid[,2], CV=T)
adjintraprob<-array()
for (i in 1:8){
  adjintraprob[i]<-mean(intralda$posterior[which(Multispid[,2]==intraspecies[i]),i],na.rm=T)
}


#wilcoxon tests for intraspecific loop

gr <- levels(Multispid[,2])
res <- matrix(NA, nrow= length(gr), ncol = length(gr), dimnames = list(gr,gr))

for (i in 1:ncol(res)){
  for (j in 1:nrow(res)){
    x<- wilcox.test(intraspidpcoa$points[,1][Multispid[,2] == gr[i]], intraspidpcoa$points[,1][Multispid[,2] == gr[j]])
    res[i,j] <- x$p.value
  }
}
res2 <- matrix(NA, nrow= length(gr), ncol = length(gr), dimnames = list(gr,gr))
for (i in 1:ncol(res2)){
  for (j in 1:nrow(res2)){
    x<- wilcox.test(intraspidpcoa$points[,2][Multispid[,2] == gr[i]], intraspidpcoa$points[,2][Multispid[,2] == gr[j]])
    res2[i,j] <- x$p.value
  }
}

res3 <- matrix(NA, nrow= length(gr), ncol = length(gr), dimnames = list(gr,gr))

for (i in 1:ncol(res3)){
  for (j in 1:nrow(res3)){
    x<- wilcox.test(intraspidpcoa$points[,3][Multispid[,2] == gr[i]], intraspidpcoa$points[,3][Multispid[,2] == gr[j]])
    res3[i,j] <- x$p.value
  }
}
#with adjusted
resintra <- matrix(NA, nrow= length(gr), ncol = length(gr), dimnames = list(gr,gr))

for (i in 1:ncol(resintra)){
  for (j in 1:nrow(resintra)){
    x<- wilcox.test(adjintraspidpcoa$points[,1][Multispid[,2] == gr[i]], adjintraspidpcoa$points[,1][Multispid[,2] == gr[j]])
    resintra[i,j] <- x$p.value
  }
}
resintra2 <- matrix(NA, nrow= length(gr), ncol = length(gr), dimnames = list(gr,gr))

for (i in 1:ncol(resintra2)){
  for (j in 1:nrow(resintra2)){
    x<- wilcox.test(adjintraspidpcoa$points[,2][Multispid[,2] == gr[i]], adjintraspidpcoa$points[,2][Multispid[,2] == gr[j]])
    resintra2[i,j] <- x$p.value
  }
}
resintra3 <- matrix(NA, nrow= length(gr), ncol = length(gr), dimnames = list(gr,gr))

for (i in 1:ncol(resintra3)){
  for (j in 1:nrow(resintra3)){
    x<- wilcox.test(adjintraspidpcoa$points[,3][Multispid[,2] == gr[i]], adjintraspidpcoa$points[,3][Multispid[,2] == gr[j]])
    resintra3[i,j] <- x$p.value
  }
}

#correlation tests (P. noblei, V. micar, E. trans, T. godeff, T. plum, H. grahami, N. albo, L. hassel against LatLong)
#axis 1
cor.test (intraspidpcoa$points[1:15,1],Multispid$Lat[1:15])
cor.test (intraspidpcoa$points[16:35,1],Multispid$Lat[16:35])
cor.test (intraspidpcoa$points[36:55,1],Multispid$Lat[36:55])
cor.test (intraspidpcoa$points[56:75,1],Multispid$Lat[56:75])
cor.test (intraspidpcoa$points[76:90,1],Multispid$Lat[76:90])
cor.test (intraspidpcoa$points[91:105,1],Multispid$Lat[91:105])
cor.test (intraspidpcoa$points[106:125,1],Multispid$Lat[106:125])
cor.test (intraspidpcoa$points[126:145,1],Multispid$Lat[126:145])
cor.test (intraspidpcoa$points[,1],Multispid$Lat)
#axis 2
cor.test (intraspidpcoa$points[1:15,2],Multispid$Lat[1:15])
cor.test (intraspidpcoa$points[16:35,2],Multispid$Lat[16:35])
cor.test (intraspidpcoa$points[36:55,2],Multispid$Lat[36:55])
cor.test (intraspidpcoa$points[56:75,2],Multispid$Lat[56:75])
cor.test (intraspidpcoa$points[76:90,2],Multispid$Lat[76:90])
cor.test (intraspidpcoa$points[91:105,2],Multispid$Lat[91:105])
cor.test (intraspidpcoa$points[106:125,2],Multispid$Lat[106:125])
cor.test (intraspidpcoa$points[126:145,2],Multispid$Lat[126:145])
cor.test (intraspidpcoa$points[,2],Multispid$Lat)
#axis 3
cor.test (intraspidpcoa$points[1:15,3],Multispid$Lat[1:15])
cor.test (intraspidpcoa$points[16:35,3],Multispid$Lat[16:35])
cor.test (intraspidpcoa$points[36:55,3],Multispid$Lat[36:55])
cor.test (intraspidpcoa$points[56:75,3],Multispid$Lat[56:75])
cor.test (intraspidpcoa$points[76:90,3],Multispid$Lat[76:90])
cor.test (intraspidpcoa$points[91:105,3],Multispid$Lat[91:105])
cor.test (intraspidpcoa$points[106:125,3],Multispid$Lat[106:125])
cor.test (intraspidpcoa$points[126:145,3],Multispid$Lat[126:145])
cor.test (intraspidpcoa$points[,3],Multispid$Lat)
#longtitude
#axis 1
cor.test (intraspidpcoa$points[1:15,1],Multispid$Long[1:15])
cor.test (intraspidpcoa$points[16:35,1],Multispid$Long[16:35])
cor.test (intraspidpcoa$points[36:55,1],Multispid$Long[36:55])
cor.test (intraspidpcoa$points[56:75,1],Multispid$Long[56:75])
cor.test (intraspidpcoa$points[76:90,1],Multispid$Long[76:90])
cor.test (intraspidpcoa$points[91:105,1],Multispid$Long[91:105])
cor.test (intraspidpcoa$points[106:125,1],Multispid$Long[106:125])
cor.test (intraspidpcoa$points[126:145,1],Multispid$Long[126:145])
cor.test (intraspidpcoa$points[,1],Multispid$Long)
#axis 2
cor.test (intraspidpcoa$points[1:15,2],Multispid$Long[1:15])
cor.test (intraspidpcoa$points[16:35,2],Multispid$Long[16:35])
cor.test (intraspidpcoa$points[36:55,2],Multispid$Long[36:55])
cor.test (intraspidpcoa$points[56:75,2],Multispid$Long[56:75])
cor.test (intraspidpcoa$points[76:90,2],Multispid$Long[76:90])
cor.test (intraspidpcoa$points[91:105,2],Multispid$Long[91:105])
cor.test (intraspidpcoa$points[106:125,2],Multispid$Long[106:125])
cor.test (intraspidpcoa$points[126:145,2],Multispid$Long[126:145])
cor.test (intraspidpcoa$points[,2],Multispid$Long)
#axis 3
cor.test (intraspidpcoa$points[1:15,3],Multispid$Long[1:15])
cor.test (intraspidpcoa$points[16:35,3],Multispid$Long[16:35])
cor.test (intraspidpcoa$points[36:55,3],Multispid$Long[36:55])
cor.test (intraspidpcoa$points[56:75,3],Multispid$Long[56:75])
cor.test (intraspidpcoa$points[76:90,3],Multispid$Long[76:90])
cor.test (intraspidpcoa$points[91:105,3],Multispid$Long[91:105])
cor.test (intraspidpcoa$points[106:125,3],Multispid$Long[106:125])
cor.test (intraspidpcoa$points[126:145,3],Multispid$Long[126:145])
cor.test (intraspidpcoa$points[,3],Multispid$Long)

#distance and variance explained for intraspecific
intracentroids<-array()
for (j in 1:length(intraspecies)){
  intrasample= sample(intraspidpcoa$points[which(Multispid[,2]==intraspecies[j]),1],nrow(intraspid[which(Multispid[,2]==intraspecies[j]),]))
  intracentroids[j]<- mean(intrasample)
}

intracentroids2<-array()
for (j in 1:length(intraspecies)){
  intrasample= sample(intraspidpcoa$points[which(Multispid[,2]==intraspecies[j]),2],nrow(intraspid[which(Multispid[,2]==intraspecies[j]),]))
  intracentroids2[j]<- mean(intrasample)
}
intracentroids3<-array()
for (j in 1:length(intraspecies)){
  intrasample= sample(intraspidpcoa$points[which(Multispid[,2]==intraspecies[j]),3],nrow(intraspid[which(Multispid[,2]==intraspecies[j]),]))
  intracentroids3[j]<- mean(intrasample)
}

intracentmat= matrix(cbind(intracentroids,intracentroids2, intracentroids3), nrow=8, ncol = 3)
intraclusters= kmeans(intraspidpcoa$points[,1:3], centers = intracentmat, nstart = 20)
#adjusted centroids
adjintracentroids<-array()
for (j in 1:length(intraspecies)){
  intrasample= sample(adjintraspidpcoa$points[which(Multispid[,2]==intraspecies[j]),1],nrow(intraspid[which(Multispid[,2]==intraspecies[j]),]))
  adjintracentroids[j]<- mean(intrasample)
}

adjintracentroids2<-array()
for (j in 1:length(intraspecies)){
  intrasample= sample(adjintraspidpcoa$points[which(Multispid[,2]==intraspecies[j]),2],nrow(intraspid[which(Multispid[,2]==intraspecies[j]),]))
  adjintracentroids2[j]<- mean(intrasample)
}
adjintracentroids3<-array()
for (j in 1:length(intraspecies)){
  intrasample= sample(adjintraspidpcoa$points[which(Multispid[,2]==intraspecies[j]),3],nrow(intraspid[which(Multispid[,2]==intraspecies[j]),]))
  adjintracentroids3[j]<- mean(intrasample)
}

adjintracentmat= matrix(cbind(adjintracentroids,adjintracentroids2, adjintracentroids3), nrow=8, ncol = 3)
adjtotintracentroids<- array()
for(j in 1:3) {
  totsample= sample(intraspidpcoa$points[,j],145)
  adjtotintracentroids[j]<- mean(totsample)
}
#distance calculations#
totintracentroids<- array()
for(j in 1:3) {
  totsample= sample(intraspidpcoa$points[,j],145)
  totintracentroids[j]<- mean(totsample)
}
allintradistance<-array()
for (i in 1:nrow(intraspid)){
  allintradistance[i]<- dist(rbind(intraspidpcoa$points[i,1:3],totintracentroids),method="euclidean")^2
}
sum(allintradistance)
#total variance explained per species
varexplained<-array()
indvspeciesdist<-array()
speciesdist<-array()
multispeciesdist<-array()
indvspeciessum<-array()
indvvarexplained<-array()
for (i in 1:nrow(intracentmat)){
  speciesrows<-which(Multispid[,2]==intraspecies[i])
  for (j in 1:length(speciesrows)){
    speciesdist[j]<-dist(rbind(intraspidpcoa$points[speciesrows[j],1:3],totintracentroids),method="euclidean")^2
    indvspeciesdist[j]<-dist(rbind(intraspidpcoa$points[speciesrows[j],1:3],intracentmat[i,]),method="euclidean")^2
  }
  indvspeciessum[i]<-sum(indvspeciesdist[1:length(speciesrows)])
  multispeciesdist[i]<- sum(speciesdist[1:length(speciesrows)])
  varexplained<-(1-sum(indvspeciessum)/sum(multispeciesdist))*100
  indvvarexplained[i]<- (1-sum(indvspeciesdist[1:length(speciesrows)])/sum(multispeciesdist[i]))*100
}

#adjusted matrix
adjvarexplained<-array()
adjindvspeciesdist<-array()
adjspeciesdist<-array()
adjmultispeciesdist<-array()
adjindvspeciessum<-array()
adjindvvarexplained<-array()
for (i in 1:nrow(adjintracentmat)){
  speciesrows<-which(Multispid[,2]==intraspecies[i])
  for (j in 1:length(speciesrows)){
    adjspeciesdist[j]<-dist(rbind(adjintraspidpcoa$points[speciesrows[j],1:3],adjtotintracentroids),method="euclidean")^2
    adjindvspeciesdist[j]<-dist(rbind(adjintraspidpcoa$points[speciesrows[j],1:3],adjintracentmat[i,]),method="euclidean")^2
  }
  adjindvspeciessum[i]<-sum(adjindvspeciesdist[1:length(speciesrows)])
  adjmultispeciesdist[i]<- sum(adjspeciesdist[1:length(speciesrows)])
  adjvarexplained<-(1-sum(adjindvspeciessum)/sum(adjmultispeciesdist))*100
  adjindvvarexplained[i]<- (1-sum(adjindvspeciesdist[1:length(speciesrows)])/sum(adjmultispeciesdist[i]))*100
}

#figures

png('Fig 01 Intramorphospace.png', width=700, height=700)
par(mar=c(4,4,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
plot(intraspidpcoa$points[,1],intraspidpcoa$points[,2],cex=0, xlab='PCoA axis 1',ylab='PCoA axis 2',xlim=c(-10,10),ylim=c(-3.5,3.5))
dataEllipse(intraspidpcoa$points[1:15,1],intraspidpcoa$points[1:15,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=1)
dataEllipse(intraspidpcoa$points[16:35,1],intraspidpcoa$points[16:35,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=16)
dataEllipse(intraspidpcoa$points[36:55,1],intraspidpcoa$points[36:55,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=16)
dataEllipse(intraspidpcoa$points[56:75,1],intraspidpcoa$points[56:75,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=1)
dataEllipse(intraspidpcoa$points[76:90,1],intraspidpcoa$points[76:90,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T,pch= 3)
dataEllipse(intraspidpcoa$points[91:105,1],intraspidpcoa$points[91:105,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#009e73",center.pch=NULL, add=T,pch=16)
dataEllipse(intraspidpcoa$points[106:125,1],intraspidpcoa$points[106:125,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#e79f00",center.pch=NULL, add=T,pch=16)
dataEllipse(intraspidpcoa$points[126:145,1],intraspidpcoa$points[126:145,2]*-1,levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#cc79a7",center.pch=NULL, add=T, pch=16)
legend(5.2,3.6,c("E. transmarina","H. grahami","L. hasseltii","P. noblei", "N. albopunctatus","T. plumipes","T. godeffroyi","V. micarioides"), col=c("#0072b2","#009e73","#cc79a7","#0072b2","#e79f00","#0072b2","#000000","#000000"),pch=c(16,16,16,1,16,3,1,16),cex=0.8,pt.cex=0.8,text.font=3)
dev.off()

#axis 1 vs 3
png('Fig 02 Intramorphospace.png', width=700, height=700)
par(mar=c(4,4,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
plot(intraspidpcoa$points[,1],intraspidpcoa$points[,3],cex=0, xlab='PCoA axis 1',ylab='PCoA axis 3',xlim=c(-10,9),ylim=c(-4,3))
dataEllipse(intraspidpcoa$points[1:15,1],intraspidpcoa$points[1:15,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=1)
dataEllipse(intraspidpcoa$points[16:35,1],intraspidpcoa$points[16:35,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=16)
dataEllipse(intraspidpcoa$points[36:55,1],intraspidpcoa$points[36:55,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=16)
dataEllipse(intraspidpcoa$points[56:75,1],intraspidpcoa$points[56:75,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=1)
dataEllipse(intraspidpcoa$points[76:90,1],intraspidpcoa$points[76:90,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T,pch= 3)
dataEllipse(intraspidpcoa$points[91:105,1],intraspidpcoa$points[91:105,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#009e73",center.pch=NULL, add=T,pch=16)
dataEllipse(intraspidpcoa$points[106:125,1],intraspidpcoa$points[106:125,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#e79f00",center.pch=NULL, add=T,pch=16)
dataEllipse(intraspidpcoa$points[126:145,1],intraspidpcoa$points[126:145,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#cc79a7",center.pch=NULL, add=T, pch=16)
legend(5.1,3.2,c("E. transmarina","H. grahami","L. hasseltii","P. noblei", "N. albopunctatus","T. plumipes","T. godeffroyi","V. micarioides"), col=c("#0072b2","#009e73","#cc79a7","#0072b2","#e79f00","#0072b2","#000000","#000000"),pch=c(16,16,16,1,16,3,1,16),cex=0.8,pt.cex=0.8)
dev.off()
#using size adjusted axis 1 vs. 2
png('Fig 03 Intramorphospace.png', width=700, height=700)
par(mar=c(4,4,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
plot(adjintraspidpcoa$points[,1],adjintraspidpcoa$points[,2],cex=0, xlab='PCoA axis 1',ylab='PCoA axis 2',xlim=c(-10,10),ylim=c(-6,6))
dataEllipse(adjintraspidpcoa$points[1:15,1],adjintraspidpcoa$points[1:15,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=1)
dataEllipse(adjintraspidpcoa$points[16:35,1],adjintraspidpcoa$points[16:35,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=16)
dataEllipse(adjintraspidpcoa$points[36:55,1],adjintraspidpcoa$points[36:55,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=16)
dataEllipse(adjintraspidpcoa$points[56:75,1],adjintraspidpcoa$points[56:75,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=1)
dataEllipse(adjintraspidpcoa$points[76:90,1],adjintraspidpcoa$points[76:90,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T,pch= 3)
dataEllipse(adjintraspidpcoa$points[91:105,1],adjintraspidpcoa$points[91:105,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#009e73",center.pch=NULL, add=T,pch=16)
dataEllipse(adjintraspidpcoa$points[106:125,1],adjintraspidpcoa$points[106:125,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#e79f00",center.pch=NULL, add=T,pch=16)
dataEllipse(adjintraspidpcoa$points[126:145,1],adjintraspidpcoa$points[126:145,2],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#cc79a7",center.pch=NULL, add=T, pch=16)
legend(5.5,6,c("E. transmarina","H. grahami","L. hasseltii","P. noblei", "N. albopunctatus","T. plumipes","T. godeffroyi","V. micarioides"), col=c("#0072b2","#009e73","#cc79a7","#0072b2","#e79f00","#0072b2","#000000","#000000"),pch=c(16,16,16,1,16,3,1,16),cex=0.8,pt.cex=0.8,text.font=3)
dev.off()
#size adjusted axis 1 vs. 3
png('Fig 04 Intramorphospace.png', width=700, height=700)
par(mar=c(4,4,2,2), mgp=c(2.5,0.75,0), cex=1.3, cex.lab=1.25, las=1)
plot(adjintraspidpcoa$points[,1],adjintraspidpcoa$points[,2],cex=0, xlab='PCoA axis 1',ylab='PCoA axis 3',xlim=c(-10,10),ylim=c(-5,5))
dataEllipse(adjintraspidpcoa$points[1:15,1],adjintraspidpcoa$points[1:15,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=1)
dataEllipse(adjintraspidpcoa$points[16:35,1],adjintraspidpcoa$points[16:35,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=16)
dataEllipse(adjintraspidpcoa$points[36:55,1],adjintraspidpcoa$points[36:55,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T, pch=16)
dataEllipse(adjintraspidpcoa$points[56:75,1],adjintraspidpcoa$points[56:75,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#000000",center.pch=NULL, add=T,pch=1)
dataEllipse(adjintraspidpcoa$points[76:90,1],adjintraspidpcoa$points[76:90,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#0072b2",center.pch=NULL, add=T,pch= 3)
dataEllipse(adjintraspidpcoa$points[91:105,1],adjintraspidpcoa$points[91:105,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#009e73",center.pch=NULL, add=T,pch=16)
dataEllipse(adjintraspidpcoa$points[106:125,1],adjintraspidpcoa$points[106:125,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#e79f00",center.pch=NULL, add=T,pch=16)
dataEllipse(adjintraspidpcoa$points[126:145,1],adjintraspidpcoa$points[126:145,3],levels=0.95,fill=T,fill.alpha=0.4,plot.points=T,col="#cc79a7",center.pch=NULL, add=T, pch=16)
legend(5,5,c("E. transmarina","H. grahami","L. hasseltii","P. noblei", "N. albopunctatus","T. plumipes","T. godeffroyi","V. micarioides"), col=c("#0072b2","#009e73","#cc79a7","#0072b2","#e79f00","#0072b2","#000000","#000000"),pch=c(16,16,16,1,16,3,1,16),cex=0.8,pt.cex=0.8,text.font=3)
dev.off()