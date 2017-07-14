###########################################################################
#  Simple example for genomewide selection
#  Ridge regression and BLUP with data from the course notes of B. Hayes
###########################################################################

## Set working, input, and output directories
   # setwd("e:/SelectionTools-examples")
  st.input.dir  <- "data"
  st.output.dir <- "output"

  # st.script.dir <- "src/"
  # source(paste(st.script.dir,"/SelectionTools.R",sep=""))

## Start R package SelectionTools
  library("SelectionTools")

###########################################################################
# (1)  Ridge Regression
# (1.1) Calculations for estimating effects step by step
###########################################################################

st.read.marker.data      ("ex-gsbh-em",format="l") # data set "default"
st.read.performance.data ("ex-gsbh-ep")            # is used

x <- st.marker.data.statistics()

x$genotypes[1:10,1:12]

x$ind.list


( gs.build.Z()            )              # (1)
( gs.lambda.const ( lambda=0.001) )      # (2)
( gs.mme.coeff()  )                      # (3)
( gs.mme.rhs()    )                      # (4)
( gs.mme.solve()  )                      # (5)

( gs.esteff.rrwe.01 ( lambda = 0.001 ) ) #  (1) - (5)

# Check how good the model fits the data:
# Estimate the breeding values of the individuals of the estimation set

gs.plot.model.fit( estimation.set = "default")

###########################################################################
# (1.2) Predicting breeding values for a new data set with ridge regression
###########################################################################

st.read.marker.data ( filename = "ex-gsbh-em", # Located in 'st.ionput.dir'   
                      format   = "l",          # List format  
                      data.set = "e" )         # Name of the dataset
 
st.read.performance.data ( filename = "ex-gsbh-ep",
                           data.set = "e" )            

gs.esteff.rrwe.01 ( lambda = 0.001,             # Estimate the effects 
                    data.set="e"  )             # with ridge regression

st.read.marker.data  ( filename = "ex-gsbh-vm", # Marker data of a new
                       format   = "l",          # data set
                       data.set = "v"  )        # named "v"

c1 <- gs.estimate.gv ( estimation.set = "e",  # Data set with effects
                       prediction.set = "v" ) # Data set with marker data

c1

st.read.performance.data ( filename = "ex-gsbh-vp", # For validation
                           data.set = "v" )         # the actual phenotypes   

gs.plot.validation ( estimation.set = "e",      # Comapre model fit and   
                     validation.set = "v",      # prediction accuracy
                     title = "Rigde Regression" )

###########################################################################
# (1.3) REML 
###########################################################################

st.read.marker.data      ("ex-gsbh-em",format="l")
st.read.performance.data ("ex-gsbh-ep")

gs.build.Z()                      # (1) Design matrix
(gs.lambda.rmlc())                # (2) Determine lamda 
(gs.lambda.emstep(constvar= T))   #     One interation step
gs.mme.coeff()                    # (3)
gs.mme.rhs()                      # (4)
gs.mme.solve()                    # (5)

(gs.esteff.rmlc())                # (1)-(5) combined

gs.plot.model.fit( training.set = "default")

dev.new()
gs.plot.validation ( estimation.set="default",
                     validation.set="v",      
                     title = "REML" )         

# REML estimation with regress

st.read.marker.data      ("ex-gsbh-em",format="l")
st.read.performance.data ("ex-gsbh-ep")
gs.esteff.rmlc()
vc <- gs.lambda.emstep(constvar= T)[,"var"]
vc[1] # Residual variance
vc[2] # Genetic variance of one marker
vc[1] / vc[2] # Shrikage factor

library(regress)

y <- st.read.performance.data ("ex-gsbh-ep")[,-1]
Z <- as.matrix( gs.build.Z()[,-1],ncol=10) 
V <- Z %*% t(Z)
rg <- regress( y ~ 1 , ~ V)
summary(rg)                        # Variance components
(sh <- rg$sigma[2] /rg$sigma[1])   # Shrinkage factor

LA <- diag( ncol(Z) ) * sh         

MME.solve <- function (y, Z, LA) 
{
  X  <- matrix (1,ncol=1,nrow=nrow(Z))
  R1 <- cbind( t(X)%*%X , t(X)%*%Z )
  R2 <- cbind( t(Z)%*%X, (t(Z) %*% Z + LA))
  LHS <- rbind(R1, R2)
  RHS <- rbind ( t(X)%*%y, t(Z)%*%y )
  bhat <- ginv(LHS) %*% RHS
  return( bhat )
}            

beta.hat <- MME.solve (y, Z, LA )
beta.hat            


###########################################################################
# 1.3) Data from Hayes script
#      Bayes
###########################################################################

dta <- st.read.marker.data      ("ex-gsbh-em",format="l")
y   <- st.read.performance.data ("ex-gsbh-ep")

# Use gmp on Linux systems
# gs.set.precision.01(128)

( gs.esteff.bayes ())

dev.new()
gs.plot.model.fit( estimation.set = "default",      
                     title = "REML")

dev.new()
gs.plot.validation ( estimation.set="default",
                     validation.set="v",      
                     title = "Bayes" )         


###########################################################################
# Bayes estimation in R
# R version from B. Hayes
# for Cormaprison
###########################################################################

{# Prepare data 
  dta <- st.read.marker.data      ("ex-gsbh-em",format="l")
  y   <- st.read.performance.data ("ex-gsbh-ep")
  x   <- gs.build.Z()
}

x[1:5,1:5]
x <- as.matrix(x[,2:ncol(x)])
x[1:5,1:5]

y <- as.matrix(y[,2:ncol(y)])
              
{
  #initializatoin g, mu, ...
  nmkr = 10
  nid  = 325
  nit  = 1000
  #
  gs    = array(0,c(nit,nmkr))
  gVars = gs
  vars  = array(0,nit)
  mjus  = vars
  #
  g    = array(0.1, nmkr)
  mu   = 0.1
  gvar = array(0.1,nmkr)
  ones = array(1,nid)
  e    = array(0,nid)
  set.seed(1)
}

 for(k in 1:nit)
    {
      cat (sprintf("Run %i of %i\r",k,nit))
      # Sample varE from inv-chi-square posterior
      e   <-  y - x%*%g - mu
      s.e <-  sum(e*e)
      rnd <-  rchisq(1,nid)
      varE <- s.e / rnd 
      # Sample the mean from a normal posterior
      mu  =rnorm(1, (sum(y)-t(ones)%*%x%*%g)/nid, sqrt(varE/nid))
      #sample gvar from a normal posterior
      for(j in 1:nmkr){
        #Xu
        gvar[j]=(g[j]*g[j])/rchisq(1,0.998)
      }
      #sample g from a normal distribution
      z=array(0,nid)
      gold <- g
      for(j in 1:nmkr){
        gTmp=gold;
        gTmp[j]=0
        for(i in 1:nid){
          z[i]=x[i,j]
        }
        #calculate the mean of the distribution
        exp <- (t(z)%*%y-t(z)%*%x%*%gTmp-t(z)%*%ones*mu)/(t(z)%*%z+varE/gvar[j])
        sde <- sqrt(varE/(t(z)%*%z+varE/gvar[j]))
        g[j]=rnorm(1,exp,sde)
      }
      #store parameter values
      for (j in 1:nmkr){
        gs[k,j]=g[j]
        gVars[k,j]=gvar[j]
      }
      vars[k]=varE
      mjus[k]=mu
    }

ev <- 501:1000
eff <- rep(0,nmkr)
for (j in 1:nmkr) eff[j] <- mean (gs[ev,j])
data.frame(eff=round(eff,digits=5))
