library(tsDyn)
library(data.table)
data("barry")
#Fruit and gases
x1<-as.data.table(read.csv(
  file = "C:/Users/bened/Desktop/Vrije/Lezioni/Dissertation/Dataset/Italy/Usati/FAOSTAT_veggiesdata_6-30-2021.csv"
))

x2<-as.data.table(read.csv(
  file = "C:/Users/bened/Desktop/Vrije/Lezioni/Dissertation/Dataset/Italy/Usati/methane.csv"
))

x3<-as.data.table(read.csv(
  file = "C:/Users/bened/Desktop/Vrije/Lezioni/Dissertation/Dataset/Italy/Usati/nitrous oxide.csv"
))


x4<-as.data.table(read.csv(
  file = "C:/Users/bened/Desktop/Vrije/Lezioni/Dissertation/Dataset/Italy/Usati/co2 world.csv"
))

x5<-as.data.table(read.csv(
  file = "C:/Users/bened/Desktop/Vrije/Lezioni/Dissertation/Dataset/Italy/Usati/FAOSTAT_metyeardata_10-17-2021.csv"
))

x6<-as.data.table(read.csv(
  file = "C:/Users/bened/Desktop/Vrije/Lezioni/Dissertation/Dataset/Italy/Usati/Rain_1961_2014.csv"
))



x1<-x1[29:55,'Value']
x2<-x2[2270:2296,'Total.including.LUCF..CH4.emissions..CAIT.']
x3<-x3[2270:2296,'Total.including.LUCF..N2O.emissions..CAIT.']
x4<-x4[10970:10996,'Annual.CO2.emissions']
x5<-x5[27:53,'Value']
x6<-x6[27:53,'Average.monthly.precipitation']


# R port of the gauss code by Luis Filipe Martins
y1 <- log(barry[,"dolcan"], base = exp(1))	
y2 <- log(barry[,"cpiUSA"], base = exp(1))	
y3 <- log(barry[,"cpiCAN"], base = exp(1))


y <- cbind(y1, y2, y3)
y <- cbind(x1, x2, x3,x4)


y <- as.data.table(y)

n <- nrow(y)
# Number of endogenous variables in the analysis
k <- ncol(y)
my.qr.eigen <- function(x){
  if(is.matrix(x) == FALSE){
    stop("\nPlease provide a quadratic matrix.\n")
  }
  x <- as.matrix(x)
  pQ <- diag(1, dim(x)[1]);
  # iterate 
  while(sum(abs(x[lower.tri(x, diag = FALSE)])) > 0.0000001){
    d <- qr(x);
    Q <- qr.Q(d);
    pQ <- pQ %*% Q;
    x <- qr.R(d) %*% Q; 
    
  }
  # Eigenvalues
  eig.values <- round(diag(x),5)
  # Eigenvectors
  eig.vec <- round(pQ,5)
  return(list(values = eig.values, vectors = eig.vec))
  
}
lag.data.table <- function(data, max.lag = 1){
  data[, shift(data, max.lag, give.names = TRUE)]
}

diffmatrix <- function(x,max.diff = 1,max.lag = 1){
  #Add if condition to make it possible to differentiate between matrix and vector                  
  if(is.vector(x) == TRUE ){
    myx <- embed(c(rep(NA,max.lag),diff(x,max.lag,max.diff)),max.diff)
    colnames(myx) <- paste("v1.d",max.diff, sep=".")
    return(myx)
  }
  
  else if(is.matrix(x) == TRUE){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
  }
  #if matrix the result should be 0, if only a vector it should be 1
  else if(as.integer(is.null(ncol(x))) == 0 ){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
    
  }
}

as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1))
# Removal of na.s
na.exclude(as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1)))

tvcoint <- function(y,p,m){
  # local ystar_1,ysub,y_1,dy,dylags,T,betau,resu,betav,resv,S00,S01,S10,S11,S00inv,S11inv,A,eval,evec,evs;
  y <- as.matrix(y)
  ystar_1 <- ycheb(y,varord = p,chebdim = m);		
  # dy is the same as the Z0 matrix in ca.jo
  # ylags would be the same as ZK
  # dylags is equivalent to Z1
  # checked the dimensions and content of all preliminary matrices against original gauss code
  dy <- na.exclude(as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1)))
  dylags <-  na.exclude(lag.data.table(dy, max.lag = p))
  dy <- dy[-(1:lags)]
  myT <- nrow(dylags);             
  dylags <- cbind(rep(1,myT), dylags) 
  # a= b/c
  # a = (c'c)^(-1)c'b
  # betau <- dy/dylags
  # betau <- (dylags'dylags)^-1 dylags dy
  betau <- as.matrix(solve(crossprod(as.matrix(dylags))) %*% crossprod(as.matrix(dylags), as.matrix(dy)))
  resu <- as.matrix(dy)-(as.matrix(dylags) %*% betau)     
  # betav <- ystar_1/dylags;     
  # betav <- (dylags'dylags)^-1 dylags ystar_1
  betav <-  as.matrix(solve(crossprod(as.matrix(dylags))) %*% crossprod(as.matrix(dylags), as.matrix(ystar_1)))
  resv <- ystar_1-(as.matrix(dylags) %*% as.matrix(betav));
  
  S00 <- (crossprod(resu))/ myT;        
  S01 <- (crossprod(resu,resv))/myT;        
  S10 <- t(S01);                
  S11 <- (crossprod(resv))/myT;       
  S00inv <- solve(S00);       
  S11inv <- solve(S11);       
  A <- S11inv %*% S10 %*% S00inv %*% S01; 
  # calculation of eigenvalues and -vectors with QR algorithm
  valeigen <- my.qr.eigen(A);
  evs <- cbind(valeigen$values,valeigen$vectors)[order(valeigen$values, decreasing = TRUE), ]
  evec = t(evs[,2:ncol(evs)])
  detS00 <- det(S00)
  output <- list("eigenvalues" = evs[, 1], "eigenvectors" = evec, "determinant" = detS00)
  return(output)
}

ycheb <- function(mydata,varord,chebdim){
  #  local i,yst,k,n,nn,ind;
  k <- ncol(mydata);
  nn <- nrow(mydata)
  yst <- matrix(nrow = nn-varord-1, ncol = (chebdim+1)*k)
  yst[, 1:k] <- as.matrix(mydata[(varord+1):(nn-1), ])
  if(chebdim == 0){
    return(yst);
  } else {
    n <- length(yst[, 1]);
    ind <- seq(from = varord+2, by = 1, length.out = n);
    i = 1;
    for(i in 1:chebdim){
      yst[, (i*k+1):((i+1)*k)] <- sqrt(2)*cos(i*pi*(ind-0.5)/n)*yst[,1:k]
    }
    return(yst);
  }
}

mmax <- round(n/10)	# maximum dimension of Chebishev Polynomials
p <- 1;			# VAR order
lrtvc <- matrix(nrow = mmax, ncol = k);	     # TVC Stat (row) m=1,...,mmax and (col) r=1,..,k
lrtvcpv <- matrix(nrow = mmax, ncol = k);	     # P-values using the Asymp Distr (chisquare(mrk))
betat <- matrix(nrow = n-p-1, ncol = k*mmax)  # b_t for r=1.  (row) i= observation; Cols 1 up to k= b1~...~bk (m=1); Cols k+1 up to 2k=b1~...~bk (m=2); ...
lnlikm <- matrix(nrow = mmax, ncol = k);	     # log likelihood for different m and r 
aic <- matrix(nrow = mmax, ncol = k);	     # AIC Akaike model selection criterion 
bic <- matrix(nrow = mmax, ncol = k);	     # BIC Schwarz msc 
hann <- matrix(nrow = mmax, ncol = k);	     # Hannan-Quinn msc 


return.tvcoint.0 <- tvcoint(y,p,0)
ll0=log(1 - return.tvcoint.0$eigenvalues, base = exp(1));

for(m in 1:mmax){
  result.tvcoint.m <- tvcoint(y,p,m); # OUT: Lambda; Eigenvectors q1...qr...q(m+1)k; det(S00)
  evect <- result.tvcoint.m$eigenvectors
  eval <- result.tvcoint.m$eigenvalues
  llm <- log(1-eval, base = exp(1));	
  ind <- seq(from = p+2, by = 1, length.out = n-p-1); 
  beta1sum <- matrix(nrow = m, ncol = length(ind));	
  beta2sum <- matrix(nrow = m, ncol = length(ind));	
  beta3sum <- matrix(nrow = m, ncol = length(ind));
  beta4sum <- matrix(nrow = m, ncol = length(ind));
  
  for(mm in 1:m){
    beta1sum[mm,] <- evect[k*mm+1,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));	
    beta2sum[mm,] <- evect[k*mm+2,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));	
    beta3sum[mm,] <- evect[k*mm+3,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));
    beta4sum[mm,] <- evect[k*mm+3,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));
  }
  
  betat[,(k*(m-1)+1):(k*m)] <- cbind(evect[(1:k), 1] + cbind(colSums(beta1sum), colSums(beta2sum), colSums(beta3sum), colSums(beta4sum)))		
  
  for(r in 1:k){
    lrtvc[m,r] <- (n-p-1)*sum(ll0[1:r]) - (n-p-1)*sum(llm[1:r]);	
    lrtvcpv[m,r] <- 1-pchisq(lrtvc[m,r],m*r*k)
    detm<-tvcoint(y,p,m)
    lnlikm[m,r] <- (log(r, base = exp(1))-k-k*log((2*pi), base = exp(1)))*(n-p-1)/2 - sum(llm[1:r])*(n-p-1)/2 - (log(detm$determinant, base = exp(1)))*(n-p-1)/2;				        
    npar <- (m+1)*k*r+r*k+k^2+(k+(p-1)*k^2);		
    aic[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*2/(n-p-1);					
    bic[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*(log(n-p-1, base= exp(1)))/(n-p-1);			
    hann[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*(log(log((n-p-1), base = exp(1)), base = exp(1)))*2/(n-p-1);	
  }
  
}

result.tvcoint.m <- tvcoint(y,p,m)

evect <- result.tvcoint.m$eigenvectors
eval <- result.tvcoint.m$eigenvalues
llm <- log(1-eval, base = exp(1));	
ind <- seq(from = p+2, by = 1, length.out = n-p-1); 
beta1sum <- matrix(nrow = m, ncol = length(ind));	
beta2sum <- matrix(nrow = m, ncol = length(ind));	
beta3sum <- matrix(nrow = m, ncol = length(ind));
beta4sum <- matrix(nrow = m, ncol = length(ind));
for(mm in 1:m){
  beta1sum[mm,] <- evect[k*mm+1,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));	
  beta2sum[mm,] <- evect[k*mm+2,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));	
  beta3sum[mm,] <- evect[k*mm+3,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));
  beta4sum[mm,] <- evect[k*mm+3,1]*sqrt(2)*cos(mm*pi*(t(ind)-0.5)/(n-p-1));
}
betat[,(k*(m-1)+1):(k*m)] <- cbind(evect[(1:k), 1] + cbind(colSums(beta1sum), colSums(beta2sum), colSums(beta3sum), colSums(beta4sum)))


#VECM AND JOHANSEN
#VECM AND JOHANSEN

library(tsDyn)
y_new<-y %>% 
  rename(
    Fruit = Value,
    CH4 = Total.including.LUCF..CH4.emissions..CAIT., C02=Annual.CO2.emissions, N2O=Total.including.LUCF..N2O.emissions..CAIT.
  )
cajo<-ca.jo(y_new,type="eigen",ecdet="none",K=2,spec="longrun",season=NULL,dumvar=NULL)
summary(cajo)
vecm<-VECM(y_new,2,r=3,include = "none",beta=NULL,estim = "ML")
summary(vecm)
print(vecm)

#STANDARDIZE DATA
library(dplyr)
y_norm <- y_new %>% mutate_at(c("Fruit", "N2O","CH4","C02"), ~(scale(.) %>% as.vector))
#NO CONST OR TREND
cajo<-ca.jo(y_norm,type="eigen",ecdet="none",K=2,spec="longrun",season=NULL,dumvar=NULL)
summary(cajo)
#TREND
cajo<-ca.jo(y_norm,type="eigen",ecdet="trend",K=2,spec="longrun",season=NULL,dumvar=NULL)
summary(cajo)
#CONST
cajo<-ca.jo(y_norm,type="eigen",ecdet="const",K=2,spec="longrun",season=NULL,dumvar=NULL)
summary(cajo)
#NO CONSTANT NOR TREND
vecm<-VECM(y_norm,1,r=1,include = "none",beta=NULL,estim = "ML")
summary(vecm)
print(vecm)
#CONST
vecm2<-VECM(y_norm,1,r=2,include = "both",beta=NULL,estim = "ML")
summary(vecm2)
print(vecm2)

#JOHANSEN USING THE RANK.TEST
rt<-rank.test(vecm2)
summary(rt)
print(rt)

rt<-rank.test(vecm2,type="eigen")
summary(rt)
print(rt)

#COINTEGRATION CONSIDERING ALL GASES AND CLIMATE VARIABLES
y <- cbind(x1, x2, x3,x4,x5,x6)
y <- as.data.table(y)

library(dplyr)
colnames(y)=c("Fruit","N2O","CH4","C02","Temp","Preci")
y_norm <- y %>% mutate_at(c("Fruit", "N2O","CH4","C02","Temp","Preci"), ~(scale(.) %>% as.vector))
#NO CONST OR TREND
cajo<-ca.jo(y_norm,type="eigen",ecdet="none",K=2,spec="longrun",season=NULL,dumvar=NULL)
summary(cajo)
#TREND
cajo<-ca.jo(y_norm,type="eigen",ecdet="trend",K=2,spec="longrun",season=NULL,dumvar=NULL)
summary(cajo)
#CONST
cajo<-ca.jo(y_norm,type="eigen",ecdet="const",K=2,spec="longrun",season=NULL,dumvar=NULL)
summary(cajo)
#NO CONSTANT NOR TREND
vecm<-VECM(y_norm,1,r=2,include = "none",beta=NULL,estim = "ML")
summary(vecm)
print(vecm)
#CONST
vecm2<-VECM(y_norm,1,r=1,include = "const",beta=NULL,estim = "ML")
summary(vecm2)
print(vecm2)

#JOHANSEN USING THE RANK.TEST
rt<-rank.test(vecm2)
summary(rt)
print(rt)

rt<-rank.test(vecm2,type="eigen")
summary(rt)
print(rt)
