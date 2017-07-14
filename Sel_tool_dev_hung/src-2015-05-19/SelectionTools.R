##################################################################
## Version information
##################################################################
st.script.date <- "2015/02/22"

##################################################################
## 0.)  Paths
##################################################################

d <- ls()
if ( ! ("st.script.dir" %in% d ) ) "st.script.dir" <<- ""
if ( ! ("st.data.dir"   %in% d ) ) "st.data.dir"   <<- ""
if ( ! ("st.input.dir"  %in% d ) ) "st.input.dir"  <<- ""
if ( ! ("st.output.dir" %in% d ) ) "st.output.dir" <<- ""

options(stringsAsFactors = FALSE) 

##################################################################
## 1.) Load compiled code
##################################################################

if (Sys.info()["sysname"] == "Linux") {
    if (Sys.info()["machine"] == "x86_64" ) {
      so.ext <- "-64.so"
    } else {
      so.ext <- ".so"
    }
  } else {
    so.ext <- ".dll"
  }

dyn.load (paste(st.script.dir,"SelectionTools",so.ext,sep=""))  

##################################################################
## 2.) R Code
##################################################################

## Base functions and genomic selection 
source(paste(st.script.dir,"genomic-selection.R",sep=""))

## Simulation
source (paste(st.script.dir,"plabsim.R",sep=""))     

## Optimize marker-assisted backcrossing
source(paste(st.script.dir,"optimize-mab.R",sep=""))

##################################################################
## 3.) Initializations
##################################################################

## Startup message
st.info(0,paste("SelectionTools",st.script.date))

## Random numbers for simulations
plabsim.init()

## Starting values
# st.set.num.threads(4)
set.info.level(0)

