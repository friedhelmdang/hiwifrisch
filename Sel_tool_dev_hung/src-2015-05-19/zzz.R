###
.onLoad <- function(lib, pkg) {

  library.dynam( "SelectionTools", pkg, lib )

  ## Startup message
  st.script.date <- "15.1.1"
  st.info(0,paste("SelectionTools",st.script.date))

  ## Random numbers for simulations
  plabsim.init()

  ## Starting values
  #  st.set.num.threads(4)
  set.info.level(0)

  ## Non defined paths
  if ( 0 == length(find("st.input.dir")) )  st.input.dir  <<- ""
  if ( 0 == length(find("st.output.dir")) ) st.output.dir <<- ""
  if ( 0 == length(find("st.data.dir")) )   st.data.dir   <<- ""
  
}
