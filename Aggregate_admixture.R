ADM=read.csv("Admixture.csv",row.names = 1)
head(ADM)

# Instruction to name bars:
barNaming <- function(vec) {
  retVec <- vec
  for (k in 2:length(vec)) {
    if (vec[k - 1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}


# Size of the plot - these are margins, also related to the end appearance of the plot
par(mar = c(1, 1, 1, 1))

#This is my selection of colors and display for K=8 
x=barplot(t(as.matrix(ADM[, 7:14])),
          col = rainbow(8),
          horiz = FALSE,
          las = 2,
          cex.axis = 0.1,
          cex.names = 0.65,
)
text(x,0.5,barNaming(ADM$Country),srt=45, col='white')

##consider phenotype


GROUPS=unique(ADM$Country);GROUPS

K=8
f=c();p=c();
for(g in GROUPS){
  print(g)
  Y=ADM[ADM$Country==g,7:14];head(Y); dim(Y)
  L<-min(10,dim(unique(Y))[1])
  print( L );print(g)
  N = length(Y[,1]);
  
  if(L<=2){
    dat_s=cbind(paste(g,0,sep='_'),Y)
    write.table(dat_s, "Medicago.csv", append=T ,col.names = FALSE, row.names = T, sep=',')
  }
  
  is_stop=FALSE;
  if(L>2){
    for(k in (L-1):2){
      print(k)
      if(!is_stop){
        
        kclus <- kmeans(Y,centers= k, iter.max=1000, nstart=10000)
        #  p=c(p,round(kclus$betweenss/kclus$totss*100))
        
        
        mind=10000;
        mink1=0; mink2=0;
        
        
        for(k1 in 1:(k-1)){
          for(k2 in (k1+1):k){
            KLD=0;
            #collapsed of k1 and k2
            Q=(kclus$centers[k1,]+kclus$centers[k2,])/2+0.0000001;
            P1= kclus$centers[k1,]+0.0000001;
            P2= kclus$centers[k2,]+0.0000001;
            KLD=sum(P1*log(P1/Q,2)+P2*log(P2/Q,2)+ Q*log(Q/P1,2)+Q*log(Q/P2,2))/2
            if(KLD<mind){
              mind=KLD;
              mink1=k1; mink2=k2;
            }#if
          }#over k2
        }#over k1
        
        f=c(f,mind)
        print(mind)
        if(mind>=0.04){ #stop collapsing
          is_stop=TRUE;
          print(paste("Number of clusters for ", g, k))
          dat_s=cbind(paste(g,kclus$cluster,sep='_'), Y)
          write.table(dat_s,  "Medicago.csv", append=T ,col.names = FALSE, row.names = T, sep=',')
          break;
        }#if
      } #if stop
    }#over k
    if(is_stop==FALSE){
      print(paste("Number of clusters for ", g, k-1))
      dat_s=cbind(paste(g,0,sep='_'), Y)
        write.table(dat_s,  "Medicago.csv",  append=T ,col.names = FALSE, row.names = T, sep=',')
    }#if stop is FALSE
  }#if L>2
}

ALL=unique(read.csv("Medicago.csv", header=F))
dim(ALL)
head(ALL)
names(ALL)=c("SAMPLE","CLUSTER",names(Y))
head(ALL)
dat_mean=aggregate(ALL[,3:(K+2)],by=list(ALL$CLUSTER),FUN='mean');head(dat_mean)
dat_sd=aggregate(ALL[,3:(K+2)],by=list(ALL$CLUSTER),FUN='sd')
dat_N=aggregate(ALL[,1],by=list(ALL$CLUSTER),FUN='length');head(dat_N)
dat_mean=merge(dat_mean,dat_N, by='Group.1')
dat_sd=merge(dat_sd,dat_N, by='Group.1')

write.csv(dat_mean,"Medicago.mean.csv")
write.csv(dat_sd,"Medicago.sd.csv")

d=dist(dat_mean[,2:(K+1)])
hc=hclust(d)

jpeg("medicago.ordered.clustered.jpg")

x=barplot(t(as.matrix(dat_mean[hc$order,2:(K+1)])), col=rainbow(K), border=NA,
         las=2)
text(x, 0.5,barNaming(dat_mean$Group.1[hc$order]), col='white',srt=45 )
dev.off()

  M=25
  ADM=read.csv('Admixture.csv',row.names = 1)
  TEST=subset(ADM, is.na(ADM$Lat))
  REF=subset(ADM, !is.na(ADM$Lat))
  Phen=read.csv('pheno.csv',row.names = 1); head(Phen)
  AMDPH=na.omit(merge(REF, Phen, by='row.names'))
  for(i in 1:dim(TEST)[1]){
    Y=TEST[i,7:14];Y
    DIST=c()
    for(j in 1:dim(AMDPH)[1]){
      Z=AMDPH[j,8:15];Z
      d=sum((Y-Z)^2);d
      DIST=c(DIST,d)
    }
    I=order(DIST)[1:M];I
    Ph=AMDPH$MaxSymptomScore[I];Ph
    Dist=DIST[I]+1.E-15;Dist
    w=min(Dist)/Dist;w
    predPh=sum(w*Ph)/sum(w);predPh
    predLat=sum(w*AMDPH$Lat[I])/sum(w);predLat
    predLon=sum(w*AMDPH$Long[I])/sum(w);predLon
  }