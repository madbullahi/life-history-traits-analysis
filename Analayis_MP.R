#####MANOVA##############

setwd("C:/Users/mxa1132-admin/OneDrive/Documents/life-history-traits-analysis")

library(readxl)
library(tidyverse)

MP <- read_excel(file.choose(), sheet = "MP")
MP
PFOS <- read_excel(file.choose(), sheet = "PFOS-Tia")
PFOS
PFOA <- read_excel(file.choose(), sheet = "PFOA-Tia")
PFOA
MP_PFOA <- read_excel(file.choose(), sheet = "MP_PFOA_Tia")
MP_PFOA
PFOS_PFOA <- read_excel(file.choose(), sheet = "PFOS_PFOA_Tia")
PFOS_PFOA


# removes NA values

MP <- MP[complete.cases(MP),]
PFOS <- PFOS[complete.cases(PFOS),]
PFOA <- PFOA[complete.cases(PFOA),]
MP_PFOA <- MP_PFOA[complete.cases(MP_PFOA),]
PFOS_PFOA <- PFOS_PFOA[complete.cases(PFOS_PFOA),]


# subset the variable Age_maturity, size_maturity, fecundity, and interval_brood

MP_traits <- MP[,c("Replicates","Genotypes","Treatment","Age_maturity", "Size_maturity", "Fecundity", "Interval_brood")]
PFOS_traits <- PFOS[,c("Replicates","Genotypes","Treatment","Age_maturity", "Size_maturity", "Fecundity", "Interval_brood")]
PFOA_traits <- PFOA[,c("Replicates","Genotypes","Treatment","Age_maturity", "Size_maturity", "Fecundity", "Invterval_brood")]
MP_PFOA_traits <- MP_PFOA[,c("Replicates","Genotypes","Treatment","Age_maturity", "Size_maturity", "Fecundity", "Interval_brood")]


# log transformation of the variables

MP_traits[4:7] <- log(MP_traits[,4:7] +0.5)
PFOS_traits[4:7] <- log(PFOS_traits[,4:7] +0.5)
PFOA_traits[4:7] <- log(PFOA_traits[,4:7] +0.5)
MP_PFOA_traits[4:7] <- log(MP_PFOA_traits[,4:7] +0.5)


#  manova

set.seed(123)
manova_results_MP <- manova(cbind(Age_maturity, Size_maturity, Fecundity, Interval_brood) ~ Treatment*Genotypes, data = MP_traits)
manova_results_PFOS <- manova(cbind(Age_maturity, Size_maturity, Fecundity, Interval_brood) ~ Treatment*Genotypes, data = PFOS_traits)
manova_results_PFOA <- manova(cbind(Age_maturity, Size_maturity, Fecundity, Invterval_brood) ~ Treatment*Genotypes, data = PFOA_traits)
manova_results_MP_PFOA <- manova(cbind(Age_maturity, Size_maturity, Fecundity, Interval_brood) ~ Treatment*Genotypes, data = MP_PFOA_traits)

# print the results

summary(manova_results_MP)
summary(manova_results_PFOS)
summary(manova_results_PFOA)
summary(manova_results_MP_PFOA)


## PTA analysis####

# convert the variables to factors

MP_traits$Treatment <- factor(MP_traits$Treatment, levels = c("control", "PET"))
MP_traits$Genotypes <- factor(MP_traits$Genotypes, levels = c("LRV-0-1", "LR2-36-01"))

#Levels

level <- as.factor(MP_traits[[3]])
taxa <- as.factor(MP_traits[[2]])
level
taxa

# create the matrix

Y <- as.matrix(MP_traits[,4:7])
Y
###############
taxaBylevel <- as.factor(paste(taxa, level)) # create the taxaBylevel variable)

print(taxaBylevel)
###########

Y <- prcomp(Y)$x # perform the PCA

#Basic data param,,eters (n=#spec, p=#ldmk, k=#dim)
n <- length(levels(taxa)) # create the n variable
n  #was n.taxa

p <- length(levels(level)) # create the p variable

p  #was n.level

k <- ncol(Y) # create the k variable
k #dimensions of phenotypic data

#####
# create the full model
#manova

lm.full <- lm(Y~ taxa * level, model = T, x = T, y = T, qr = T) # create the model

yhat.full <- predict(lm.full) # predict the model

# for resid randomization (used later)
lm.red <- lm(Y ~ taxa +level, model = T, x = T, y = T, qr = T) # create the model just for the covariate

yhat.red <- predict(lm.red) # predict the model

res.red <- resid(lm.red) # create the residuals

res.red
####################
lsmeans.obs <- NULL
for (i in 1:k) {
  mean.temp <- tapply(yhat.full[,i], taxaBylevel, mean)
  lsmeans.obs <- cbind(lsmeans.obs, mean.temp)  
}

lsmeans.obs

summary(manova(lm.full))


#############

arrayspecs <- function(A){
  n.m <-NULL # n.m stands for new means
  for(i in 1:n){
    temp <-as.matrix(A[((1+(i-1)*p):(i*p)),1:k])
    n.m <-cbind(n.m, temp)}
  trajectories <- array(n.m, dim = c(p, k, n))
}

### function calculating Trajectory attributes

# Pathlength distance
pathdist<-function(M) {as.matrix(dist(M))}
trajsize<-function(M){
  traj.pathdist<-array(0,dim=c(n,1))   		#loop across trajectories
  for (i in 1:n){
    temp<-pathdist(M[,,i])
    for (j in 1:(p-1)){
      traj.pathdist[i]<-traj.pathdist[i]+temp[j,j+1]
    }
  }
  traj.size.dist<-as.matrix(dist(traj.pathdist))		#TrajSizeDiff
}

#trajectory direction 
orient<-function(M) {(svd(var(M))$v[1:k,1])} 	#find orientation

trajorient<-function(M){
  traj.orient<-array(NA,dim=c(n,k))   		#loop across trajectories
  check.1<-array(NA,dim=c(n))
  for (i in 1:n){
    temp<-orient(M[,,i])
    traj.orient[i,]<-temp
    check.1[i]<-M[1,,i]%*%traj.orient[i,]  #check startingpoint location
    check.1[i]<-check.1[i]/abs(check.1[i])
    if(check.1[i]==-1) traj.orient[i,]<--1*traj.orient[i,]
  }
  options(warn=-1)				#b/c acos of 1 (e.g, on diagonal) yields warning
  traj.ang.diff<-(180/pi)*acos(traj.orient%*%t(traj.orient))
  #diag(traj.ang.diff)<-0
}

#trajectory shape
### GPA: following J.claude 2008; Morphometrics in R

trans <-function(A){scale(A, scale=F)} ##TRANSLATION


csize<-function(A)  ##CSIZE
{p<-dim(A)[1]
size<-sqrt(sum(apply(A,2,var))*(p-1))
list("centrois_size"=size, "scaled"= A/size)
}


mshape<-function(A){apply(A,c(1,2), mean)}  #meanshape

pPsup<- function(M1,M2){                    ## OPA rotation 1-->2
  k<-ncol(M1)
  Z1<-trans(csize(M1)[[2]])
  Z2<-trans(csize(M2)[[2]])
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$u; V<-sv$u;Delt<-sv$d
  sig<-sign(det(t(Z1)%*%Z2))
  Delt[k]<-sig*abs(Delt[k]); V[,k]<-sig*V[,k]
  Gam<-U%*%t(V)
  beta<-sum(Delt)
  list(Mp1=beta*Z1%*%Gam,Mp2=Z2,rotation=Gam,scale=beta,
       df=sqrt(1-beta^2))
}
pgpa<-function(A)
{p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]  
temp2<-temp1<-array(NA,dim=c(p,k,n)); Siz<-numeric(n)#; Qm2<-numeric(n)
for (i in 1:n)
{Acs<-csize(A[,,i])
Siz[i]<-Acs[[1]]
temp1[,,i]<-trans(Acs[[2]])}
Qm1<-dist(t(matrix(temp1,k*p,n)))
Q<-sum(Qm1); iter<-0
while (abs(Q)> 0.00001)
{for(i in 1:n){
  M<-mshape(temp1[,,-i])
  temp2[,,i]<-pPsup(temp1[,,i],M)[[1]]}
  Qm2<-dist(t(matrix(temp2,k*p,n)))
  Q<-sum(Qm1)-sum(Qm2)
  Qm1<-Qm2
  iter=iter+1
  temp1<-temp2}
list("rotated"=temp2,"it.number"=iter,"Q"=Q,"intereucl.dist"=Qm2,"mshape"=
       csize(mshape(temp2))[[2]],"cent.size"=Siz)
}
###############

## loop for GPA and shape distances
trajshape <- function(M){
  x<-pgpa(M)
  traj.shape.dist<-as.matrix(dist(x$intereucl.dist))
}
###############

## TrajectrorySummaryStat

sumstat<-function(M){
  M<-as.dist(M)
  x<-if(length(M)>1)(x=var(M)) else 0
  
}

###############

#########Main Loop########

traj.specs.obs<-arrayspecs(lsmeans.obs) #lsmeans.obs is the output of the main loop
trajsize.obs<-trajsize(traj.specs.obs)
trajdir.obs<-trajorient(traj.specs.obs)
diag(trajdir.obs)<-0 #b/c some NA/N values on diagonal)
trajshape.obs<-trajshape(traj.specs.obs)
sumstatsize.obs<-sumstat(trajsize.obs)
sumstatdir.obs<-sumstat(trajdir.obs)
sumstatshape.obs<-sumstat(trajshape.obs)

###############

#### PERMUTATION Procedure

permute<-1000
line<-nrow(Y)
PSize<-POrient<-PShape<-array(1,dim=c(n,n))
PSumSize<-PSumOrient<-PSumShape<-1
for(i in 1:permute){
  line.rand<-sample(line,replace=FALSE)
  res.temp<-cbind(line.rand,res.red)
  z<-(order(line.rand))
  res.temp2<-as.matrix(res.temp[z,])
  res.p<-res.temp2[,-1] # Rows of residual are now randomized
  y.rand<-yhat.red+res.p
  yhat.rand<-predict(lm(y.rand~taxa*level,model=T,x=T,y=T,qr=T))
  
  lsmeans.rand<-NULL
  for (j in 1:k){
    mean.temp<-tapply(yhat.rand[,j],taxaBylevel, mean)
    lsmeans.rand<-cbind(lsmeans.rand,mean.temp)
  }
  traj.specs.rand<-arrayspecs(lsmeans.rand)
  trajsize.rand<-trajsize(traj.specs.rand)
  trajdir.rand<-trajorient(traj.specs.rand)
  diag(trajdir.rand)<-0 #b/c some NA/N values on diagonal)
  trajshape.rand<-trajshape(traj.specs.rand)
  sumstatsize.rand<-sumstat(trajsize.rand)
  sumstatdir.rand<-sumstat(trajdir.rand)
  sumstatshape.rand<-sumstat(trajshape.rand)
  for(j in 1:n){
    for(jj in 1:n){
      PSize[j,jj]<-if(trajsize.rand[j,jj]>=trajsize.obs[j,jj])
        (PSize[j,jj]+1) else PSize[j,jj]
      POrient[j,jj]<-if(trajdir.rand[j,jj]>=trajdir.obs[j,jj])
        (POrient[j,jj]+1) else POrient[j,jj]
      PShape[j,jj]<-if(trajshape.rand[j,jj]>=trajshape.obs[j,jj])
        (PShape[j,jj]+1) else PShape[j,jj]
    }
  }
  PSumSize<-if(sumstatsize.rand>=sumstatsize.obs) (PSumSize+1) else PSumSize
  PSumOrient<-if(sumstatdir.rand>=sumstatdir.obs) (PSumOrient+1) else PSumOrient
  PSumShape<-if(sumstatshape.rand>=sumstatshape.obs) (PSumShape+1) else PSumShape
} #end of permutation loop
for (j in 1:n){
  for(jj in 1:n){
    PSize[j,jj]<-PSize[j,jj]/(permute+1)
    POrient[j,jj]<-POrient[j,jj]/(permute+1)
    PShape[j,jj]<-PShape[j,jj]/(permute+1)
  }
}
PSumSize<-PSumSize/(permute+1)
PSumOrient<-PSumOrient/(permute+1)
PSumShape<-PSumShape/(permute+1)

###############


trajsize.obs # trajectore magnitude, length of line
PSize
trajdir.obs # trajectory direction, distance between two dots
POrient 
#trajshape.obs # trajectory shape, distance between two lines
#PShape

###NOTE: summary statistics only valid for 2+ trajectories!!
sumstatsize.obs
PSumSize
sumstatdir.obs
PSumOrient
#sumstatshape.obs
#PSumShape


#############Plot Data and LS means
plot(Y[,1:2], type="n", xlab="PC1", ylab="PC2", asp=1)
YY=cbind(taxa,level,Y)
points(YY[,3:4],col=taxa, pch=c(1,19)[as.numeric(level)])
for(i in 1:n){
  for(j in 1:(p-1)){
    points(traj.specs.obs[(j:(j+1)),1,i],traj.specs.obs[(j:(j+1)),2,i],type="l",pch=21,col="black")
    
  }
}
for (i in 1:n){
  for(j in 2:(p-1)){
    points(traj.specs.obs[j,1,i],traj.specs.obs[j,2,i],pch=21,col="red", cex=1.5)
  }
}
for(i in 1:n){
  points(traj.specs.obs[1,1,i],traj.specs.obs[1,2,i],pch=21,col="blue", cex=1.5, bg="white")
}
for(i in 1:n){
  points(traj.specs.obs[p,1,i],traj.specs.obs[p,2,i],pch=21,col="black", cex=1.5, bg="black")
}