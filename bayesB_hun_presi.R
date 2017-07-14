setwd("/home/friedhelm/Schreibtisch/praktikum_c/workshop/")

#Loading *so file
dyn.unload("bayesB_hung_nogsl.so")
dyn.load("bayesB_hung_nogsl.so")

#Defining function
baysB <-function(x, y, numit, propSegs, numMHCycles, chi_par) {
  
  gS <- matrix(0, nrow=numit, ncol=nmarkers)
  #muS <- matrix(0, c(numit))
  #vareS <- matrix(0, c(numit))
  gvareS <- matrix(0, nrow=numit, ncol=nmarkers)
  
  asdf<-.C("baysB", nrecords=as.integer(nrow(x)), nmarkers=as.integer(ncol(x)), numit=as.integer(numit), 
           numMHCycles = as.integer(numMHCycles), propSegs=as.numeric(propSegs), chi_par=as.numeric(chi_par),
           x=as.numeric(x), y= as.numeric(y), gStore=as.numeric(gS), gvarS= as.numeric(gvareS))
  store<-list(matrix(asdf$gStore,nrow=numit), matrix(asdf$gvarS,nrow=numit), asdf$muS, asdf$vareStore)
  names(store) <- c("gStore", "gvarStore", "muStore", "vareStore")
  
  return(store)
}



######read in file#####
numit <-1000
nmarkers<-3
propSegs <- 0.66
numMHCycles = 30
chi_par = 0.998;

x_c <- "/home/friedhelm/Schreibtisch/praktikum_c/workshop/xvec.inp"
y_c <- "/home/friedhelm/Schreibtisch/praktikum_c/workshop/yvec.inp"

x <-.C("read_in_x", as.character(x_c))


#####Parameters#####
numit <-1000
nmarkers<-3
propSegs <- 0.66
numMHCycles = 30
chi_par = 0.998;

#####Data#####
x <- matrix(scan("xvec.inp"), ncol=nmarkers, byrow=TRUE)
y <- matrix(scan("yvec.inp"), byrow=TRUE)

#####Execute BayesB####
starttime<-Sys.time()
a<-baysB(x,y,numit,propSegs,numMHCycles,chi_par)
endtime <-Sys.time()
duration_c <- endtime-starttime
duration_c

printout<-.C("printOut_test", as.integer(3))

plus<-.C("test1plus2", as.double(0))

numit_2 <- 30
b_wait <- matrix(ncol=3)

for (i in 1:numit_2){
  print(i)
  a<-baysB(x,y,numit,propSegs,numMHCycles,chi_par)
  g_calc_mean<-c(mean(a$gStore[100:numit,1]),mean(a$gStore[100:numit,2]),mean(a$gStore[100:numit,3]))
  b_wait <- rbind(b_wait, g_calc_mean)
  Sys.sleep(10)
}



#Mean of the marker effects
g_calc_mean<-c(mean(a$gStore[100:numit,1]),mean(a$gStore[100:numit,2]),mean(a$gStore[100:numit,3]))
g_calc_mean

#Data of new records
x_neu<- rbind(t(c(2,1,1)),t(c(1,1,1)),t(c(0,1,2)),t(c(1,2,1)),t(c(0,2,1)),
              t(c(1,1,1)),t(c(2,2,2)),t(c(2,2,1)),t(c(2,1,1)),t(c(0,1,2)))

TBV<-c(4,1,-4,1,-2,1,2,4,4,-4)s

GEBV <- x_neu%*%g_calc_mean
GEBV_sum <- matrix(nrow=nrow(x_neu),ncol=1)

for(i in 1:nrow(x_neu)){
  GEBV_sum[i]<-sum(GEBV[i,])
}

#####Plot TBV and result######
plot(TBV, cex = 1.25, pch=19, col = "darkgreen", ylim = c(-6,6))
title("BayesB")
points(GEBV_sum, pch=19, cex=1.25, col="red")
legend("right",c("TBV","GEBV"), col = c("darkgreen","red"), pch=19)
g_calc_mean
