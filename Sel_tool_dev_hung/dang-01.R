###########################################################################
# Example for genomic selection
###########################################################################

## Set working, input, and output directories

   st.input.dir  <- "input"
   st.output.dir <- "output"
   st.script.dir <- "src/"
   source(paste(st.script.dir,"/SelectionTools.R",sep=""))

   st.set.info.level(1)

###########################################################################
# Predict marker effects and genotypic values
###########################################################################

st.read.marker.data ( "ex-perez-t.mpo",   # Marker data of the training
                      format   = "m" ,   # set in matrix format
                      data.set = "t"  )  # Name of the data set

st.read.performance.data ( "ex-perez-yld.phe", # Yield
                           data.set ="t" )    

st.read.marker.data ( "ex-perez-p.mpo",  # Marker data of the prediction
                      format   = "m" ,  # set in matrix format
                      data.set = "p"  ) 

gs.esteff.rr ( method="rmlct", data.set="t" )               

gs.esteff.bayesB (hsq=0.8, data.set="t")


c1 <- gs.predict.genotypes ( training.set   = "t",  # Predict genotypic values 
                             prediction.set = "p" ) 
c1[1:20,]                                      

effekte <- gs.return.effects(data.set="t")

head(effekte)


###########################################################################
# Validate the predictions
###########################################################################


st.read.performance.data ( "ex-perez-yld.phe", # Yield of the lines of the 
                           data.set="p")      # prediction set

# pdf (file="ex-perez-01.pdf",width=8,height=4,pointsize=11)
gs.plot.validation ( training.set   ="t",
                     prediction.set ="p", title="blup")
# dev.off()

###########################################################################
# Prediction methods
###########################################################################

st.set.info.level (1)

gs.esteff.rr ( method="rmlct",  # REML, constant marker variances
               hsq=0.75,        # Start value for iterations
               data.set="t" )

gs.esteff.rr ( method="rrwet",  # Constant marker variances
               hsq=0.75,        # Preliminary estimate of the heritability
               data.set="t" )

gs.esteff.rr ( method="rmlat",  # Constant marker variances
               hsq=0.75,        # Preliminary estimate of the heritability
               data.set="t" )

gs.esteff.rr ( method="rmlc", hsq=0.75,      # Start value for iterations
               data.set="t" ) 
                                
gs.esteff.rr ( method="rrwe", hsq=0.9, data.set="t" )   
gs.esteff.rr ( method="rmlv", data.set="t" ) 
gs.esteff.rr ( method="rmla", data.set="t" ) 
gs.esteff.rr ( method="rrwa", hsq=0.9, data.set="t")  

library("regress")
gs.esteff.blup(data.set="t")

###########################################################################
# Compare accuracy of prediction methods
###########################################################################

st.read.marker.data      ( "ex-perez-t.mpo", format="m", data.set="t" )
st.read.performance.data ( "ex-perez-yld.phe", data.set="t" )    
st.read.marker.data      ( "ex-perez-p.mpo", format="m", data.set="p" )

st.read.performance.data ( "ex-perez-yld.phe", data.set="p")  
st.copy.marker.data ( "t1", source.data.set="t" )
st.copy.marker.data ( "t2", source.data.set="t" )

gs.esteff.rr ( method="rmlct", data.set="t1" )
gs.esteff.rr ( method="rmlat", data.set="t2" )               


pdf (file="ex-perez-01a.pdf",width=8,height=4,pointsize=11)
gs.plot.validation ( training.set   ="t1",
                     prediction.set ="p", title="RMLCT")
dev.off()

pdf (file="ex-perez-01b.pdf",width=8,height=4,pointsize=11)
gs.plot.validation ( training.set   ="t2",
                     prediction.set ="p", title="RMLAT")
dev.off()

###########################################################################
# Assess prediction accuracy with cross validation
###########################################################################

st.read.marker.data ( "ex-perez-c.mpo", format="m", data.set="c" )
st.read.performance.data ( "ex-perez-yld.phe", data.set ="c" )

r.rmlc <- gs.cross.validation ( data.set            = "c"         ,
                                estimation.method   = "rmlc"      ,
                                n.ts                = 306 - 306/5 ,
                                n.runs              = 100   )

r.rmla <- gs.cross.validation ( data.set             = "c"         ,
                                estimation.method   = "rmla"      ,
                                n.ts                = 306 - 306/5 ,
                                n.runs              = 100   )

mean(r.rmlc$cor)
mean(r.rmla$cor)

###########################################################################
# Comparison of effect sizes
###########################################################################

st.copy.marker.data ( "t1", "t" ); st.copy.marker.data ( "t2", "t" )
st.copy.marker.data ( "t3", "t" ); st.copy.marker.data ( "t4", "t" )
st.copy.marker.data ( "t5", "t" )

gs.esteff.rr ( method="rmlc", hsq=0.9, data.set="t1")  
gs.esteff.rr ( method="rrwe", hsq=0.9, data.set="t2")  
gs.esteff.rr ( method="rmla", hsq=0.9, data.set="t3")  
gs.esteff.rr ( method="rrwa", hsq=0.9, data.set="t4")  
gs.esteff.rr ( method="rmlv", hsq=0.9, maxiter=50000, data.set="t5")  


# pdf (file="ex-perez-02.pdf",width=8,height=8,pointsize=11)

par(mfrow=c(2,2))
   
gs.compare.effects (data.set.a="t1",  label.a ="RMLC",
                    data.set.b="t2",  label.b ="RRWE")

gs.compare.effects (data.set.a="t3",  label.a ="RMLA",
                    data.set.b="t4",  label.b ="RRWA")

gs.compare.effects (data.set.a="t1",  label.a ="RMLC",
                    data.set.b="t4",  label.b ="RMLA")

gs.compare.effects (data.set.a="t1",  label.a ="RMLC",
                    data.set.b="t5",  label.b ="RMLV")

# dev.off()

###########################################################################
# Comparison of computation time
# Notebook with Intel i7-3740QM with 32 GB RAM
###########################################################################

# > gs.esteff.rr ( method="rmlc", hsq=0.9, data.set="t1")  
# I (data set 't1'): Number of EM iterations: 48, precision: 0.000993
# I (data set 't1'): (Greatest change of a var. comp. in the last iteration:  0.099%)
# I: Elapsed time: 14.25s
 
# > gs.esteff.rr ( method="rrwe", hsq=0.9, data.set="t2")  
# I: Elapsed time: 0.15s
 
# > gs.esteff.rr ( method="rmla", hsq=0.9, data.set="t3")  
# I (data set 't3'): Number of EM iterations: 48, precision: 0.000993
# I (data set 't3'): (Greatest change of a var. comp. in the last iteration:  0.099%)
# I: Elapsed time: 14.23s
 
# > gs.esteff.rr ( method="rrwa", hsq=0.9, data.set="t4")  
# I: Elapsed time: 0.16s
 
# > gs.esteff.rr ( method="rmlv", hsq=0.9, maxiter=50000, data.set="t5")
# I (data set 't5'): Number of EM iterations: 48, precision: 0.000993
# I (data set 't5'): (Greatest change of a var. comp. in the last iteration:  0.099%)
# I (data set 't5'): Number of EM iterations: 3770, precision: 0.000996
# I (data set 't5'): (Greatest change of a var. comp. in the last iteration:  0.100%)
# I: Elapsed time: 18.58m
 
# > gs.esteff.rr ( method="rmlct", hsq=0.9, data.set="t1")
# I (data set 't1'): Number of EM iterations: 43, precision: 0.000962
# I (data set 't1'): (Greatest change of a var. comp. in the last iteration:  0.096%)
# I: Elapsed time: 0.15s
 
gs.esteff.rr ( method="rmlct", hsq=0.9, data.set="t1")  
gs.esteff.rr ( method="rrwet", hsq=0.9, data.set="t2")  
gs.esteff.rr ( method="rmlat", hsq=0.9, data.set="t3")  
library(regress); gs.esteff.blup ( data.set="t4")  
 
# > gs.esteff.rr ( method="rrwet", hsq=0.9, data.set="t2")
# I: Elapsed time: 0.01s

# > gs.esteff.rr ( method="rmlat", hsq=0.9, data.set="t3")
# I (data set 't3'): Number of EM iterations: 43, precision: 0.000962
# I (data set 't3'): (Greatest change of a var. comp. in the last iteration:  0.096%)
# I: Elapsed time: 0.14s
 
# > gs.esteff.blup ( data.set="t4")
# I: Elapsed time: 1.31s
 
###########################################################################
# Comparison of computation time
# Notebook with Intel i7-3740QM with 32 GB RAM
###########################################################################

# >   st.read.performance.data ( "large-1-gdc.p" )
# M (data set 'default'): No. of performance data: 931
# M (data set 'default'): No. of individuals: 931, no. of markers: 23988

# > gs.esteff.rr   (  method="rmlct" )
# I (data set 'default'): Number of EM iterations: 46, precision: 0.000948
# I (data set 'default'): (Greatest change of a var. comp. in the last iteration:  0.095%)
# I: Elapsed time: 6.49s

# > gs.esteff.rr   (  method="rmlat" )
# I (data set 'default'): Number of EM iterations: 46, precision: 0.000948
# I (data set 'default'): (Greatest change of a var. comp. in the last iteration:  0.095%)
# I: Elapsed time: 7.02s


