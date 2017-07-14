{
###########################################################################
# Example for genomic selection with sugar beet data
###########################################################################

library(SelectionTools)
  
   st.input.dir  <- "data"
   st.output.dir <- "output"
   st.script.dir <- "src/"
   source(paste(st.script.dir,"/SelectionTools.R",sep=""))
   
###########################################################################
# Load the data in a data set called "e"

{st.set.info.level(-2)
 st.read.marker.data ( "ex-gsr373.pop", format="l", data.set = "e" )
 st.restrict.marker.data ( NoAll.MAX = 2 ,  data.set  = "e")  
 st.restrict.marker.data ( MaMis.MAX = 0.2, data.set  = "e")  
 st.restrict.marker.data ( ExHet.MIN = 0.1, data.set  = "e")  
 st.restrict.marker.data ( InMis.MAX = 0.2, data.set  = "e")  
 st.read.performance.data ( "ex-gsr373sc.phe", data.set = "e"    )
 st.read.map ( "ex-gsr37.map", m.stretch=100, data.set="e" )
 st.read.marker.data ("ex-gsr377.pop", format="l", data.set="v" )
 st.restrict.marker.data ( NoAll.MAX = 2 ,  data.set  = "v")  
 st.restrict.marker.data ( MaMis.MAX = 0.2, data.set  = "v")  
 st.restrict.marker.data ( ExHet.MIN = 0.1, data.set  = "v")  
 st.restrict.marker.data ( InMis.MAX = 0.2, data.set  = "v")  
 st.read.performance.data ( "ex-gsr377sc.phe", data.set = "v" )
 st.read.map ( "ex-gsr37.map", m.stretch=100, data.set="v" ) 
 st.copy.marker.data ( target.data.set = "e1", source.data.set = "e" )
 st.copy.marker.data ( target.data.set = "e2", source.data.set = "e" )
 st.copy.marker.data ( target.data.set = "e3", source.data.set = "e" )
 st.copy.marker.data ( target.data.set = "e4", source.data.set = "e" )
 st.copy.marker.data ( target.data.set = "e5", source.data.set = "e" )
 st.set.info.level(0) }


###########################################################################
# Estimate effects with BLUP , constant variances

{ gs.esteff.rmlc ( data.set="e1" )
  gs.esteff.rrwe ( hsq=0.9, data.set="e2")
  gs.esteff.rmlv ( data.set="e3")
  gs.esteff.rmlr ( data.set="e4")
  gs.esteff.rrwr ( hsq=0.9, data.set="e5")   }


gs.compare.effects (data.set.a="e1", data.set.b="e5" )

gs.compare.effects (data.set.a="e4", data.set.b="e5" )
gs.compare.effects (data.set.a="e4", data.set.b="e5" )
      

gs.plot.model.fit ( estimation.set = "e2")

gs.compare.effects ( data.set.a="e4", data.set.b="e5" )

gs.compare.effects ( data.set.a="e1", data.set.b="e3" )
gs.compare.effects ( data.set.a="e1", data.set.b="e4" )
gs.compare.effects ( data.set.a="e3", data.set.b="e4" )
gs.compare.effects ( data.set.a="e3", data.set.b="e1" )


gs.plot.validation ( estimation.set="e1", validation.set="v")
gs.plot.validation ( estimation.set="e3", validation.set="v")
gs.plot.validation ( estimation.set="e4", validation.set="v")
gs.plot.validation ( estimation.set="e5", validation.set="v")



