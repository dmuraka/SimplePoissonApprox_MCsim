
##########################################################################
########## Monte Carlo simulation experiments explained in
########## "Case 1: Basic over-dispersed Poisson regression model" section
########## in Murakami and Matsui (2021) Improved log-Gaussian approximation for
########## over-dispersed Poisson regression: application to spatial analysis of COVID-19
######################################################################

########## Function for generating data
# b   : True regression coefficients
# dpar: Dispersion parameter
# n   : Sample size

dgpFun   <-function( b, dpar, n ){
  
  rqpois <- function(n, mu, dpar) {
    rnbinom(n = n, mu = mu, size = mu/(dpar-1))
  }
  
  x1    <- rnorm(n)
  x2    <- rnorm(n)
  x     <- as.matrix(cbind(1,x1,x2))
  mu    <- exp(x%*%b)
  y     <- rqpois(n,mu=mu,dpar=dpar)
  
  return(list(x=x, y=y))
}


########################################################
########## Run Monte Carlo experiments #################
########################################################

library(MASS)
nsim  <- 1000       # Number of iteration
b     <- c(NA,2,0.5)# True regression coefficients (NA for intercept b[0])
dpar  <- 5          # True dispersion parameter
c     <- 0.5        # Parameter to avoid taking the logarithm of zero
# (not used in the proposed method)

Rmse1<-NULL         # RMSE for b[1]
Rmse2<-NULL         #          b[2]
Rmse3<-NULL         #          b[3]
Bias1<-NULL         # Bias for b[1]
Bias2<-NULL         #          b[2]
Bias3<-NULL         #          b[3]
Mean1_bse<-NULL     # Mean standard error for b[1]
Mean2_bse<-NULL     #          b[2]
Mean3_bse<-NULL     #          b[3]
for(n in c(50, 200)){############## Sample size
  for(const in c(-2,-1,0,1,2)){###### True value for the intercept (b[1])
    b[1]<- const
    B1  <-matrix(0,nrow=nsim,ncol=6)
    B2  <-matrix(0,nrow=nsim,ncol=6)
    B3  <-matrix(0,nrow=nsim,ncol=6)
    Bse1<-matrix(0,nrow=nsim,ncol=6)
    Bse2<-matrix(0,nrow=nsim,ncol=6)
    Bse3<-matrix(0,nrow=nsim,ncol=6)
    for(i in 1:nsim){
      
      test_error<-TRUE
      while(test_error){
        
        ###################### Data generation
        dgp     <- dgpFun( b = b, dpar = dpar, n = n)
        x       <- dgp$x
        y       <- dgp$y
        
        ###################### Negative binomial regression
        ###### (To avoid an error, NB model is estimated first)
        tres<-try(mod2  <- glm.nb(y~0+x) )
        test_error<- class(tres)[1]=="try-error"
        
      }
      
      ######################## Poisson regression
      mod1  <- glm(y~0+x, family = poisson(link="log"))
      
      ######################## Overdispersed Poisson regression
      mod3  <- glm(y~0+x, family = quasipoisson(link="log"))
      
      ######################## Proposed approximation ################
      zrat  <- sum(y==0)/n
      yy    <- log(y+0.5) - (1+0.5*zrat)/(y+0.5)
      mod4  <- lm(yy~0+x, weights=(y+0.5))
      xb4   <- predict(mod4)
      mu4   <- exp(xb4)
      sig4  <- sum( (y - mu4)^2 / mu4 )/(n-length(b))
      zz4   <- xb4 + (y - mu4)/mu4
      XXinv4<- solve(t(x)%*% (mu4*x) )
      beta4 <- XXinv4 %*% (t(x)%*% (mu4*zz4))
      bse4  <- sqrt(sig4)*sqrt(diag( XXinv4 ) )
      
      ######################## Taylor approximation
      yy   <-log(y+c) - c/(y+c)
      mod5  <-lm(yy~0+x, weights=(y+c))
      
      ######################## Basic log-normal approximation
      yy   <-log(y+c)
      mod6  <-lm(yy~0+x, weights=(y+c))
      
      ######################## Summarize coefficient estimates
      B1[i,]<-c(coefficients(mod1)[1],coefficients(mod2)[1],
                coefficients(mod3)[1],beta4[1], coefficients(mod5)[1], coefficients(mod6)[1])
      B2[i,]<-c(coefficients(mod1)[2],coefficients(mod2)[2],
                coefficients(mod3)[2],beta4[2], coefficients(mod5)[2], coefficients(mod6)[2])
      B3[i,]<-c(coefficients(mod1)[3],coefficients(mod2)[3],
                coefficients(mod3)[3],beta4[3], coefficients(mod5)[3], coefficients(mod6)[3])
      
      ######################## Summarize coefficient standard errors
      Bse1[i,]<-c(coefficients(summary(mod1))[1,2],coefficients(summary(mod2))[1,2],coefficients(summary(mod3))[1,2],bse4[1],
                  coefficients(summary(mod5))[1,2],coefficients(summary(mod6))[1,2])
      Bse2[i,]<-c(coefficients(summary(mod1))[2,2],coefficients(summary(mod2))[2,2],coefficients(summary(mod3))[2,2],bse4[2],
                  coefficients(summary(mod5))[2,2],coefficients(summary(mod6))[2,2])
      Bse3[i,]<-c(coefficients(summary(mod1))[3,2],coefficients(summary(mod2))[3,2],coefficients(summary(mod3))[3,2],bse4[3],
                  coefficients(summary(mod5))[3,2],coefficients(summary(mod6))[3,2])
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
    
    print(const)
  }
}

################################## Summarize results (RMSE)
rrr1<-as.data.frame(rbind(cbind(Rmse1[,1:3],b[2],b[3],Rmse1[,-(1:3)])))
rrr2<-as.data.frame(rbind(cbind(Rmse2[,1:3],b[2],b[3],Rmse2[,-(1:3)])))
rrr3<-as.data.frame(rbind(cbind(Rmse3[,1:3],b[2],b[3],Rmse3[,-(1:3)])))
names(rrr1)<-c("c","n","b0","b1","b2","pois","nb","qpois","ours","taylor","usual")
names(rrr2)<-c("c","n","b0","b1","b2","pois","nb","qpois","ours","taylor","usual")
names(rrr3)<-c("c","n","b0","b1","b2","pois","nb","qpois","ours","taylor","usual")

rrr1      # RMSE of b[1]
rrr2      # RMSE of b[2]
rrr3      # RMSE of b[3]

################################## Summarize results (Bias)
bbb1<-as.data.frame(rbind(cbind(Bias1[,1:3],b[2],b[3],Bias1[,-(1:3)])))
bbb2<-as.data.frame(rbind(cbind(Bias2[,1:3],b[2],b[3],Bias2[,-(1:3)])))
bbb3<-as.data.frame(rbind(cbind(Bias3[,1:3],b[2],b[3],Bias3[,-(1:3)])))
names(bbb1)<-c("c","n","b0","b1","b2","poisson","nb","od_poisson","proposed","taylor","loglinear")
names(bbb2)<-c("c","n","b0","b1","b2","poisson","nb","od_poisson","proposed","taylor","loglinear")
names(bbb3)<-c("c","n","b0","b1","b2","poisson","nb","od_poisson","proposed","taylor","loglinear")

bbb1      # Bias of b[1]
bbb2      # Bias of b[2]
bbb3      # Bias of b[3]

################################## Summarize results (Mean standard error)
sss1<-as.data.frame(rbind(cbind(Mean1_bse[,1:3],b[2],b[3],Mean1_bse[,-(1:3)])))
sss2<-as.data.frame(rbind(cbind(Mean2_bse[,1:3],b[2],b[3],Mean2_bse[,-(1:3)])))
sss3<-as.data.frame(rbind(cbind(Mean3_bse[,1:3],b[2],b[3],Mean3_bse[,-(1:3)])))
names(sss1)<-c("c","n","b0","b1","b2","poisson","nb","od_poisson","proposed","taylor","loglinear")
names(sss2)<-c("c","n","b0","b1","b2","poisson","nb","od_poisson","proposed","taylor","loglinear")
names(sss3)<-c("c","n","b0","b1","b2","poisson","nb","od_poisson","proposed","taylor","loglinear")

sss1      # Mean standard error of b[1]
sss2      # Mean standard error of b[2]
sss3      # Mean standard error of b[3]
