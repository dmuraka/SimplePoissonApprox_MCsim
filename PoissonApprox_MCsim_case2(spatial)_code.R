
##########################################################################
########## Monte Carlo simulation experiments explained in
########## "Case 2: Model with spatial effects" section
########## in Murakami and Matsui (2021) Improved log-Gaussian approximation for
########## over-dispersed Poisson regression: application to spatial analysis of COVID-19
######################################################################

########## Function for generating data
# b   : True regression coefficients
# dpar: Dispersion parameter
# n   : Sample size

dgp_sp<-function(b, dpar, n){

  rqpois <- function(n, mu, dpar){
    rnbinom(n = n, mu = mu, size = mu/(dpar-1))
  }

  coords  <- cbind(rnorm(n),rnorm(n))
  meig    <- meigen(coords)
  sf      <- meig$sf
  ev      <- meig$ev
  r       <- rnorm(ncol(sf),sd=sqrt(ev))

  x	      <- cbind(1, rnorm(n), rnorm(n))
  xb      <- x %*% b
  sp_comp <- sd(xb)*scale(sf%*%r)
  y0      <- xb + sp_comp
  y       <- rqpois(n,mu=exp(y0),dpar=dpar)

  x	      <- as.data.frame(x)
  names(x)<- paste("var",1:length( b ),sep="")

  return(list(y=y, x=x, coords=coords, sp_comp=sp_comp))
}


#######################################################################################
#######################################################################################

#install.packages("mgcv")
#install.packages("spmoran")
library(mgcv);library(spmoran)
nsim  <- 1000       # Number of iteration
b     <- c(NA,2,0.5)# True regression coefficients (NA for intercept b[0])
dpar  <- 5          # True dispersion parameter
c     <- 0.5        # Parameter to avoid taking the logarithm of zero

Rmse1<-NULL         # RMSE for b[1]
Rmse2<-NULL         #          b[2]
Rmse3<-NULL         #          b[3]
Bias1<-NULL         # Bias for b[1]
Bias2<-NULL         #          b[2]
Bias3<-NULL         #          b[3]
Mean1_bse<-NULL     # Mean standard error for b[1]
Mean2_bse<-NULL     #          b[2]
Mean3_bse<-NULL     #          b[3]
Rmse_ranef<-NULL    # RMSE for the spatial random effects
for(const in c(-2,-1,0,1,2)){# True value for the intercept (b[1])
  for(n in c(50,200)){       # Sample size
    b[1]       <- const
    B1         <- matrix(0,nrow=nsim,ncol=4)
    B2         <- matrix(0,nrow=nsim,ncol=4)
    B3         <- matrix(0,nrow=nsim,ncol=4)
    Bse1       <- matrix(0,nrow=nsim,ncol=4)
    Bse2       <- matrix(0,nrow=nsim,ncol=4)
    Bse3       <- matrix(0,nrow=nsim,ncol=4)
    Rmse_ranef0<- matrix(0,nrow=nsim,ncol=4)
    for(i in 1:nsim){

      dum<-1
      ecount <-0
      while(dum==1){
        ######### Data generation
        dgp     <- dgp_sp( b = b, dpar = dpar, n = n)
        x       <- dgp$x
        y       <- dgp$y
        coords  <- dgp$coords
        sp_comp <- dgp$sp_comp

        ######### Poisson regression
        test1 <-try(mod1   <- gam(method="ML", family="poisson",y~x[,2]+x[,3]+s(coords[,1],coords[,2],m=2,bs="gp")))

        ######### Overdispersed Poisson regression
        test2 <-try(mod2   <- gam(method="ML", family="quasipoisson",y~x[,2]+x[,3]+s(coords[,1],coords[,2],m=2,bs="gp")))

        ######### Proposed approximation #####################
        zrat  <- sum(y==0)/n
        y3    <- log(y+0.5) - (1+0.5*zrat)/(y+0.5)
        test3 <- try(mod3   <- gam(method="ML", weights=y+0.5, y3~x[,2]+x[,3]+s(coords[,1],coords[,2],m=2,bs="gp")))

        ######### Taylor approximation
        y4    <- log(y+c)- c/(y+c)
        test4<-try(mod4     <- gam(method="ML", weights=y+0.5, y4~x[,2]+x[,3]+s(coords[,1],coords[,2],m=2,bs="gp")))

        test5<-(class(test1)[1]=="try-error")|(class(test2)[1]=="try-error")|(class(test3)[1]=="try-error")|(class(test4)[1]=="try-error")
        dum<-ifelse(test5,1,0)
        if(class(test5)[1]=="try-error"){
          ecount<-ecount + 1
          print(ecount)
        }
      }

      ######################## Summarize coefficient estimates
      B1[i,]<-c(coefficients(mod1)[1],coefficients(mod2)[1],coefficients(mod3)[1],coefficients(mod4)[1])
      B2[i,]<-c(coefficients(mod1)[2],coefficients(mod2)[2],coefficients(mod3)[2],coefficients(mod4)[2])
      B3[i,]<-c(coefficients(mod1)[3],coefficients(mod2)[3],coefficients(mod3)[3],coefficients(mod4)[3])

      ######################## Summarize coefficient standard errors
      Bse1[i,]<-c(summary(mod1)$se[1],summary(mod2)$se[1],summary(mod3)$se[1],summary(mod4)$se[1])
      Bse2[i,]<-c(summary(mod1)$se[2],summary(mod2)$se[2],summary(mod3)$se[2],summary(mod4)$se[2])
      Bse3[i,]<-c(summary(mod1)$se[3],summary(mod2)$se[3],summary(mod3)$se[3],summary(mod4)$se[3])

      ######################## Summarize spatial random effects estimates
      ranef <- cbind(predict(mod1,type="terms")[,3],predict(mod2,type="terms")[,3],predict(mod3,type="terms")[,3],predict(mod4,type="terms")[,3])
      Rmse_ranef0[i,]<-sqrt(colMeans((ranef-c(sp_comp))^2))
    }

    rmse1<-colSums((B1-b[1])^2/nsim)#### RMSE
    rmse2<-colSums((B2-b[2])^2/nsim)
    rmse3<-colSums((B3-b[3])^2/nsim)

    bias1<-colSums((B1-b[1])/nsim)  #### Bias
    bias2<-colSums((B2-b[2])/nsim)
    bias3<-colSums((B3-b[3])/nsim)

    mean1_bse<-colSums(Bse1/nsim)   #### Mean coefficient standard errors
    mean2_bse<-colSums(Bse2/nsim)
    mean3_bse<-colSums(Bse3/nsim)

    ################################### Summarize error statistics
    Rmse1<-rbind(Rmse1,c(c,n,const,rmse1))
    Rmse2<-rbind(Rmse2,c(c,n,const,rmse2))
    Rmse3<-rbind(Rmse3,c(c,n,const,rmse3))
    Bias1<-rbind(Bias1,c(c,n,const,bias1))
    Bias2<-rbind(Bias2,c(c,n,const,bias2))
    Bias3<-rbind(Bias3,c(c,n,const,bias3))
    Mean1_bse<-rbind(Mean1_bse,c(c,n,const,mean1_bse))
    Mean2_bse<-rbind(Mean2_bse,c(c,n,const,mean2_bse))
    Mean3_bse<-rbind(Mean3_bse,c(c,n,const,mean3_bse))

    Rmse_ranef<-rbind(Rmse_ranef,c(c,n,const,t(colMeans(Rmse_ranef0)) ))
    print(const)
  }
}

################################## Summarize results (RMSE)
rrr1<-as.data.frame(rbind(cbind(Rmse1[,1:3],b[2],b[3],Rmse1[,-(1:3)])))
rrr2<-as.data.frame(rbind(cbind(Rmse2[,1:3],b[2],b[3],Rmse2[,-(1:3)])))
rrr3<-as.data.frame(rbind(cbind(Rmse3[,1:3],b[2],b[3],Rmse3[,-(1:3)])))
names(rrr1)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")
names(rrr2)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")
names(rrr3)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")

rrr1      # RMSE of b[1]
rrr2      # RMSE of b[2]
rrr3      # RMSE of b[3]

################################## Summarize results (Bias)
bbb1<-as.data.frame(rbind(cbind(Bias1[,1:3],b[2],b[3],Bias1[,-(1:3)])))
bbb2<-as.data.frame(rbind(cbind(Bias2[,1:3],b[2],b[3],Bias2[,-(1:3)])))
bbb3<-as.data.frame(rbind(cbind(Bias3[,1:3],b[2],b[3],Bias3[,-(1:3)])))
names(bbb1)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")
names(bbb2)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")
names(bbb3)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")

bbb1      # Bias of b[1]
bbb2      # Bias of b[2]
bbb3      # Bias of b[3]

sss1_bse<-as.data.frame(rbind(cbind(Mean1_bse[,1:3],b[2],b[3],Mean1_bse[,-(1:3)])))
sss2_bse<-as.data.frame(rbind(cbind(Mean2_bse[,1:3],b[2],b[3],Mean2_bse[,-(1:3)])))
sss3_bse<-as.data.frame(rbind(cbind(Mean3_bse[,1:3],b[2],b[3],Mean3_bse[,-(1:3)])))
names(sss1_bse)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")
names(sss2_bse)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")
names(sss3_bse)<-c("c","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")

sss1      # Standard error of b[1]
sss2      # Standard error of b[2]
sss3      # Standard error of b[3]

rmse_rr <-as.data.frame(rbind(cbind(Rmse_ranef[,1:3],b[2],b[3],Rmse_ranef[,-(1:3)])))
names(rmse_rr)<-c("x","n","b0","b1","b2","poisson","od_poisson","proposed","taylor")

rmse_rr   # RMSE of the spatial random effects






