
## Startup message

st.info <-  function( lev, msg )
  {
    c <- .C( "g_r_info",
             as.integer(lev),
             as.character(msg) );
  }

"st.set.num.threads" <- function( active )
{
  c  <-  .C( "gs_set_num_threads", as.integer(active) )
}

"st.set.openblas.threads" <- function( active )
{
  c  <-  .C( "gs_set_openblas_threads", as.integer(active) )
}


"st.get.num.threads" <- function()
{
  active <- 0;
  c  <-  .C( "gs_get_num_threads", t=as.integer(active) )
  return( c$t )
}

###################################################################################
# Genomic selection functions: gs.
###################################################################################

"gs.info" <-  function( lev, msg ){
    c <- .C("gs_r_info",
            as.integer(lev),
            as.character(msg)
            );
  }

"gs.set.num.threads" <- function( active )
{
  c  <-  .C( "gs_set_num_threads", as.integer(active) )
}

"gs.reset" <-  function()
{
    c <- .C("gs_reset");
  }

"gs.info" <-  function( lev, msg ){
    c <- .C("gs_r_info",
            as.integer(lev),
            as.character(msg)
            );
  }

"gs.set.info.level" <- function( level,data.set ="default" )
{
  c  <-  .C( "gs_set_info_level_GV",
            as.integer(level),
            as.character(data.set)   )
}

"gs.set.all.info.levels" <- function( level )
{
  c  <-  .C( "gs_set_all_info_levels",
            as.integer(level)   )
}


"gs.get.info.level" <- function ( data.set ="default" )
{
  level = 0
  c  <-  .C(  "gs_get_info_level_GV",
              level = as.integer(level),
              as.character(data.set)      )
  return (c$level)
}

"gs.test.mp" <- function(n = 0) {
  
  c  <-  .C(
    "gs_test_mp",
    as.integer(n)
  )
}


"st.set.matrix.ops" <- function(select.s="mt",select.p="st") {
  c <- .C( "st_set_matrix_ops",
           as.character(select.s),
           as.character(select.p)   )
  return(invisible(NULL))
}

"gs.start.timer" <- function (depth=1)
  {
    st.start.timer(depth)
  }


"gs.stop.timer" <- function ( depth=1, info.level=0 )
{
  st.stop.timer(depth,info.level)
 }


"gs.build.Z" <- function (
                              out.filename = "Z.matrix" ,
                              auxfiles     = TRUE,
                              data.set     = "default"
                             )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  X.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_build_Z_01_GV"          ,
          as.character(X.filename)   ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)      
          )
  X <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {

    h <- read.table(file=X.filename,header=F,nrows=1,stringsAsFactors=F)
    X <- read.table(file=X.filename,header=F,skip=1)
    unlink(X.filename)
    names(X) <- c("i",as.character(h))
  }
  gs.stop.timer (info.level=1)
  return( invisible(X) )
}


"gs.mme.coeff" <- function(
                              out.filename = "mme.coeff" ,
                              auxfiles     = TRUE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  A.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mme_coeff_01_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)      
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=F) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

"gs.build.V" <- function(
                              out.filename = "V.matrix" ,
                              auxfiles     = TRUE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  A.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_build_V_01_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)      
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=F) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}


"gs.mme.restcoeff" <- function( out.filename = "mme.coeff" ,
                                  auxfiles     = TRUE         )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  A.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mme_restcoeff_01_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval)
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=F) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

"gs.mme.rhs" <- function ( out.filename = "mme.rhs" ,
                              auxfiles     = TRUE      ,
                              data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  b.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mme_rhs_01_GV",
          as.character(b.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)        )
  b <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    b <- data.matrix( read.table(file=b.filename,header=F) )
    unlink(b.filename)
    colnames(b) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(b) )
}


"gs.mme.restrhs" <- function ( out.filename = "mme.rhs" ,
                                 auxfiles = TRUE              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  b.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mme_restrhs_01_GV",
          as.character(b.filename),
          as.integer(auxfiles),
          retval = as.integer(retval)  )
  b <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    b <- data.matrix( read.table(file=b.filename,header=F) )
    unlink(b.filename)
    colnames(b) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(b) )
}

"gs.mme.solve" <- function( out.filename = "mme.solution" ,
                               auxfiles     = TRUE           ,
                               data.set     = "default"      )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mme_solve_01_GV",
          as.character(x.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)        )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.mme.invcoeff" <- function( out.filename = "inv.mme" ,
                                  auxfiles = TRUE              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  o.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mme_invcoeff_01_GV",
          as.character(o.filename),
          as.integer(auxfiles),
          retval = as.integer(retval)  )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=o.filename,header=F) )
    unlink(o.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

## "gs.check.model.fit" <- function ( cor.filename = "breeding.values"  ,
##                                       auxfiles = TRUE                   )
## {
##   gs.start.timer()
##   o.filename <- paste(st.output.dir,op,cor.filename,sep="")
##   retval <- 0;
##   c <- .C ( "gs_check_model_fit_01_GV"  ,
##             as.character(o.filename)    ,
##             as.integer(auxfiles)        ,
##             retval = as.integer(retval) )
##   x <- NULL
##   if  ( (auxfiles) && (-2 != c$retval ) ) {
##     x <- data.matrix( read.table(file=o.filename,header=T) )
##     unlink(o.filename)
##   }
##   gs.stop.timer (info.level=1)
##   return( invisible(x) )
##}

## "gs.estimate.breeding.values" <- function (
##           mar.filename ,
##           cor.filename = "breeding.values",
##           auxfiles = TRUE                          )
## {
##   gs.start.timer()
##   if ( "" == st.input.dir )  ip <- "" else ip <- "/"
##   if ( "" == st.output.dir ) op <- "" else op <- "/"
##   m.filename <- paste(st.input.dir,ip,mar.filename,sep="")
##   o.filename <- paste(st.output.dir,op,cor.filename,sep="")
##   retval <- 0;
##   c <- .C ( "gs_estimate_breeding_values_01_GV",
##             as.character(m.filename),
##             as.character(o.filename),
##             as.integer(auxfiles),
##             retval = as.integer(retval))
##   x <- NULL
##   if  ( (auxfiles) && (-2 != c$retval ) ) {
##     x <- read.table(file=o.filename,header=F,skip=1,stringsAsFactors=F) 
##     unlink(o.filename)
##     if (2 == ncol(x)) names(x) <- c("i","yhat")
##     if (3 == ncol(x)) names(x) <- c("i","y","yhat")
##   }
##   gs.stop.timer (info.level=1)
##   return( invisible(x) )
## }


## "gs.validation" <- function ( mar.filename ,
##                                 per.filename ,
##                                 cor.filename = "predicted",
##                                 auxfiles = TRUE )
## {
##   gs.start.timer()
##   if ( "" == st.input.dir )  ip <- "" else ip <- "/"
##   if ( "" == st.output.dir ) op <- "" else op <- "/"
##   m.filename <- paste(st.input.dir,ip,mar.filename,sep="")
##   p.filename <- paste(st.input.dir,ip,per.filename,sep="")
##   o.filename <- paste(st.output.dir,op,cor.filename,sep="")
##   retval <- 0;
##   c <- .C ( "gs_validation_01_GV",
##             as.character(m.filename),
##             as.character(p.filename),
##             as.character(o.filename),
##             as.integer(auxfiles),
##             retval = as.integer(retval))
##   x <- NULL
##   if  ( (auxfiles) && (-2 != c$retval ) ) {
##     x <- data.matrix( read.table(file=o.filename,header=T) )
##     unlink(o.filename)
##   }
##   gs.stop.timer (info.level=1)
##   return( invisible(x) )
## }

"gs.single.marker.reg" <- function(  out.filename = "sm.reg" ,
                                        auxfiles     = TRUE      ,
                                        data.set     = "default"  )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_single_marker_reg_01_GV",
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)        )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.single.marker.aov" <- function(  alpha = 1 ,
                                     out.filename = "sm.reg" ,
                                     auxfiles     = TRUE      ,
                                     data.set     = "default"  )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_single_marker_aov_01_GV",
            as.double(alpha)           ,
            as.character(x.filename)    ,
            as.integer(auxfiles)        ,
            retval = as.integer(retval) ,
            as.character(data.set)        )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.lambda.const" <- function( lambda                  ,
                                  out.filename = "lambda" ,
                                  auxfiles = TRUE         ,
                                  data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
 c <- .C( "gs_lambda_const_01_GV",
           as.double(lambda),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval) ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.hsq" <- function( hsq,
                             out.filename = "lambda"  ,
                             auxfiles     = TRUE      ,
                             data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_const_02_GV",
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval) ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}



"gs.lambda.reg" <- function( hsq                      ,
                                out.filename = "lambda"  ,
                                auxfiles = TRUE          ,
                                data.set     = "default"  )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_reg_01_GV"        ,
           as.double(hsq)               ,
           as.character(x.filename)     ,
           as.integer(auxfiles)         ,
           retval = as.integer(retval)  ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.aov" <- function( hsq                      ,
                             alpha = 1                ,
                             out.filename = "lambda"  ,
                             auxfiles = TRUE          ,
                             data.set     = "default"  )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_aov_01_GV"        ,
           as.double(hsq)               ,
           as.double(alpha)             ,
           as.character(x.filename)     ,
           as.integer(auxfiles)         ,
           retval = as.integer(retval)  ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.emstep" <- function (constvar     = F        ,
                                   out.filename = "lambda" ,
                                   auxfiles     = TRUE,
                                   data.set     = "default"     )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_emstep_01_GV"    ,
           as.integer(constvar)        ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.lambda.rmlv" <- function( maxiter     = 1000       ,
                                 precision   = 0.001      ,
                                 out.filename = "lambda"  ,
                                 auxfiles = TRUE          ,
                                 data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_rmlv_01_GV"     ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.rmlc" <- function( maxiter     = 1000       ,
                                 precision   = 0.001      ,
                                 out.filename = "lambda"  ,
                                 auxfiles = TRUE          ,
                                 data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_rmlc_01_GV"     ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 

"gs.lambda.rmlr" <- function( maxiter     = 1000       ,
                                 precision   = 0.001      ,
                                 out.filename = "lambda"  ,
                                 auxfiles = TRUE          ,
                                 data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_rmlr_01_GV"     ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 

"gs.lambda.rmla" <- function( alpha        = 1        ,
                              maxiter      = 1000     ,
                              precision    = 0.001    ,
                              out.filename = "lambda" ,
                              auxfiles     = TRUE     ,
                              data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_rmla_01_GV"     ,
           as.double(alpha)           ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 



"gs.esteff.rrwr" <- function( hsq                       ,
                                 out.filename = "eff.rrwr" ,
                                 auxfiles = TRUE           ,
                                 data.set ="default"        )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rrwr_01_GV"      ,
           as.double(hsq)              ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

"gs.esteff.rrwa" <- function( hsq                       ,
                              alpha = 1                 ,
                              out.filename = "eff.rrwa" ,
                              auxfiles = TRUE           ,
                              data.set ="default"        )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rrwa_01_GV"      ,
           as.double(hsq)              ,
           as.double(alpha)            ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}



"gs.esteff.rrwe.01" <- function ( lambda                    ,
                                  out.filename = "eff.rrwe" ,
                                  auxfiles = TRUE           ,
                                  data.set ="default"         )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rrwe_01_GV",
           as.double(lambda),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

"gs.esteff.rrwe" <- function( hsq ,
                                out.filename = "eff.rrwe" ,
                                auxfiles = TRUE,
                                data.set ="default"              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rrwe_02_GV",
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}



"gs.set.allele.codes" <- function(    aa = 0,
                                      aA = 1,
                                      AA = 2,
                                      an = 0.5,
                                      nn = 1  ,
                                      nA = 1.5,
                                      data.set ="default" )
{
  c <- .C( "gs_set_allele_codes_01_GV",
           as.double(aa),
           as.double(aA),
           as.double(AA),
           as.double(an),
           as.double(nn),
           as.double(nA),
           as.character(data.set)   )
  return( invisible(NULL) )
}


"gs.build.esvs" <- function (es.size         ,
                                base.set       = "default" ,
                                estimation.set = "none"    ,
                                validation.set = "none"    ,
                                es_m = "cves_m" ,
                                es_p = "cves_p" ,
                                vs_m = "cvvs_m" ,
                                vs_p = "cvvs_p" ,
                                auxfiles = FALSE                )
{
  if ( "" == st.output.dir ) op <- "" else op <- "/"
    out.es_m    <- paste(st.output.dir,op,es_m,sep="") 
    out.es_p    <- paste(st.output.dir,op,es_p,sep="") 
    out.vs_m    <- paste(st.output.dir,op,vs_m,sep="") 
    out.vs_p    <- paste(st.output.dir,op,vs_p,sep="") 
    gs.start.timer()
    retval <- 0;
    c <- .C("gs_build_esvs_01_GV"           ,
            as.integer(es.size)          ,
            as.character(base.set)       ,            
            as.character(estimation.set) ,            
            as.character(validation.set) ,            
            as.character(out.es_m)       ,
            as.character(out.es_p)       ,
            as.character(out.vs_m)       ,
            as.character(out.vs_p)       ,
            as.integer  (auxfiles)       ,
            retval = as.integer(retval)   )
    gs.stop.timer (info.level=1)
  }

"gs.build.tsps" <- function (ts.size         ,
                             base.set       = "default" ,
                             training.set   = "none"    ,
                             prediction.set = "none"    ,
                             es_m = "cves_m" ,
                             es_p = "cves_p" ,
                             vs_m = "cvvs_m" ,
                             vs_p = "cvvs_p" ,
                             auxfiles = FALSE                )
{
  if ( "" == st.output.dir ) op <- "" else op <- "/"
    out.es_m    <- paste(st.output.dir,op,es_m,sep="") 
    out.es_p    <- paste(st.output.dir,op,es_p,sep="") 
    out.vs_m    <- paste(st.output.dir,op,vs_m,sep="") 
    out.vs_p    <- paste(st.output.dir,op,vs_p,sep="") 
    gs.start.timer()
    retval <- 0;
    c <- .C("gs_build_esvs_01_GV"           ,
            as.integer(ts.size)          ,
            as.character(base.set)       ,            
            as.character(training.set) ,            
            as.character(prediction.set) ,            
            as.character(out.es_m)       ,
            as.character(out.es_p)       ,
            as.character(out.vs_m)       ,
            as.character(out.vs_p)       ,
            as.integer  (auxfiles)       ,
            retval = as.integer(retval)   )
    gs.stop.timer (info.level=1)
  }


"gs.estimate.gv" <- function ( estimation.set ="default",
                               prediction.set ="default",
                               out.filename = "breeding.values",
                               auxfiles = TRUE                    )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  o.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C ( "gs_estimate_bv_01_GV",
            as.character(estimation.set)  ,
            as.character(prediction.set)  ,
            as.character(o.filename),
            as.integer(auxfiles)    ,
            retval = as.integer(retval)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- read.table(file=o.filename,header=F,skip=1,stringsAsFactors=F)
    unlink(o.filename)
    if (3 == ncol(x))  names(x) <- c("i","y","yhat")
    if (2 == ncol(x))  names(x) <- c("i","yhat")
  }
    gs.stop.timer (info.level=1)
  return( invisible(x) )
}
                                        


"gs.predict.genotypes" <- function ( training.set ="default",
                                     prediction.set ="default",
                                     out.filename = "breeding.values",
                                     auxfiles = TRUE                    )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  o.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C ( "gs_estimate_bv_01_GV",
            as.character(training.set)  ,
            as.character(prediction.set)  ,
            as.character(o.filename),
            as.integer(auxfiles)    ,
            retval = as.integer(retval)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- read.table(file=o.filename,header=F,skip=1,stringsAsFactors=F)
    unlink(o.filename)
    if (3 == ncol(x))  names(x) <- c("i","y","yhat")
    if (2 == ncol(x))  names(x) <- c("i","yhat")
  }
    gs.stop.timer (info.level=1)
  return( invisible(x) )
}




"gs.esteff.lsq" <- function ( out.filename = "eff.lsq" ,
                              auxfiles = TRUE           ,
                              data.set ="default"       )
{
  gs.start.timer()
  #
  Z <- gs.build.Z (data.set=data.set)
  r <-qr(as.matrix(Z[,2:ncol(Z)]))$rank
  n <- ncol(Z)-1
  if (r<n) { 
    m <- sprintf("Design matrix is not of full rank (%i < %i)",r,n)
    st.info(-2,m)
    gs.stop.timer (info.level=1)
    return (invisible(NULL))
  }
  #
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_lsq_01_GV",
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.testeff.lsq" <- function ( nperm = 0,
                               out.filename = "test.lsq" ,
                               auxfiles = TRUE           ,
                               data.set ="default"       )
{
  gs.start.timer()
  #
  Z <- gs.build.Z (data.set=data.set)
  r <-qr(as.matrix(Z[,2:ncol(Z)]))$rank
  n <- ncol(Z)-1
  if (r<n) { 
    m <- sprintf("Design matrix is not of full rank (%i < %i)",r,n)
    st.info(-2,m)
    gs.stop.timer (info.level=1)
    return (invisible(NULL))
  }
  #
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_testeff_lsq_01_GV",
           as.integer(nperm),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <-  read.table(file=x.filename,header=T)
#   x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.esteff.rmlc" <- function( maxiter     = 1000        ,
                              precision   = 0.001       ,
                              out.filename = "eff.rmlc" ,
                              auxfiles = TRUE           ,
                              data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rmlc_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.testeff.rmlc" <- function( nperm       = 1000        ,
                               maxiter     = 1000        ,
                               precision   = 0.001       ,
                               out.filename = "test.rmlc" ,
                               auxfiles = TRUE           ,
                               data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_testeff_rmlc_01_GV",
           as.integer(nperm),
           as.integer(maxiter),
           as.double(precision),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <-  read.table(file=x.filename,header=T)
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.esteff.rmlv" <- function(  maxiter     = 1000        ,
                                  precision   = 0.001       ,
                                  out.filename = "eff.rmlv" ,
                                  auxfiles = TRUE           ,
                                  data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rmlv_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

"gs.testeff.rmlv" <- function( nperm       = 1000        ,
                               maxiter     = 1000        ,
                               precision   = 0.001       ,
                               out.filename = "test.rmlv" ,
                               auxfiles = TRUE           ,
                               data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_testeff_rmlv_01_GV",
           as.integer(nperm),
           as.integer(maxiter),
           as.double(precision),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <-  read.table(file=x.filename,header=T)
#   x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.esteff.rmlr" <- function(  maxiter     = 1000        ,
                                  precision   = 0.001       ,
                                  out.filename = "eff.rmlr" ,
                                  auxfiles = TRUE           ,
                                  data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rmlr_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

"gs.esteff.rmla" <- function( alpha       = 1        ,
                              maxiter     = 1000        ,
                              precision   = 0.001       ,
                              out.filename = "eff.rmla" ,
                              auxfiles = TRUE           ,
                              data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rmla_01_GV",
           as.double(alpha),
           as.integer(maxiter),
           as.double(precision),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

"gs.esteff.bayes" <- function(  maxiter     = 1000        ,
                                   precision   = 0.001       ,
                                   gibbsburn   = 1000        ,
                                   gibbsiter   = 1000        ,
                                   out.filename = "eff.bayes" ,
                                   auxfiles = TRUE           ,
                                   data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_bayes_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.integer(gibbsburn),
           as.integer(gibbsiter),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.set.precision" <- function ( prec )
{
  c <- .C( "gs_set_precision_01_GV" ,
            as.integer(prec)        )
}


"gs.cross.validation" <- function ( estimation.method          ,
                                    n.ts                       ,
                                    n.runs                     ,
                                    hsq           = 0.9        ,
                                    alpha         = 1          , 
                                    maxiter       = 1000       ,
                                    precision     = 0.001      ,
                                    gibbsburn     = 1000       ,
                                    gibbsiter     = 1000       ,
                                    out.filename.r  = "crosval"  ,
                                    out.filename.u  = "effects"  ,
                                    auxfiles      = TRUE       ,
                                    data.set      = "default"  )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  r.filename <- paste(st.output.dir,op,out.filename.r,sep="")
  u.filename <- paste(st.output.dir,op,out.filename.u,sep="")
  retval <- 0;
  c <- .C( "gs_cross_validation_01_GV" ,
           as.character(estimation.method)        ,
           as.integer  (n.ts)  ,
           as.integer  (n.runs)        ,
           as.double   (hsq)           ,
           as.double   (alpha)         ,
           as.integer  (maxiter)       ,
           as.double   (precision)     ,
           as.integer  (gibbsburn),
           as.integer  (gibbsiter),
           as.character(r.filename)    ,
           as.character(u.filename)    ,
           as.integer  (auxfiles)      ,
           retval = as.integer(retval) ,
           as.character(data.set)       )
  outp <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <- data.matrix( read.table(file=r.filename,header=T) )
      unlink(r.filename)
      outp$cor <- x
      x <- data.matrix( read.table(file=u.filename,header=T) )
      unlink(u.filename)
      outp$eff <- x
  }
  
  gs.stop.timer (info.level=1)
  return( invisible(outp) )
}

"gs.restrict.marker.data.01" <- function ( outfile,
                                          traitfile = NULL,
                                          NoAll.MAX = 999,
                                          MaMis.MAX = 1,
                                          ExHet.MIN = 0,
                                          InMis.MAX = 1,
                                          data.set  = "default" )
{
    if ( "" == st.input.dir )  ip <- "" else ip <- "/"
    if ( "" == st.output.dir ) op <- "" else op <- "/"
    out.filename    <- paste(st.input.dir,op,outfile,sep="") 
    if ( !is.null(traitfile) )
      trait.filename  <- paste(st.input.dir,ip,traitfile,sep="") 
    else
      trait.filename = ""
    gs.start.timer()
    retval <- 0;
    c <- .C("gs_restrict_marker_data_01_GV"    ,
            as.character(out.filename)   ,
            as.character(trait.filename) ,
            as.integer(NoAll.MAX),
            as.double (MaMis.MAX),
            as.double (ExHet.MIN),
            as.double (InMis.MAX),
            retval = as.integer(retval),
            as.character(data.set)       )
    gs.stop.timer (info.level=1)
  }


###################################################################################
# Selection tools functions: st.
###################################################################################


"st.reset" <-  function(){
    c <- .C("gs_reset");
  }

"st.info" <-  function( lev, msg ){
    c <- .C(  "st_r_info",
              as.integer(lev),
              as.character(msg)    );
  }

"st.set.info.level" <- function( level )
{
  c  <-  .C( "st_set_info_level",
            as.integer(level)  )
}

"st.get.info.level" <- function ()
{
  level = 0
  c  <-  .C(  "st_get_info_level",
              level = as.integer(level)   )
  return (c$level)
}

"st.start.timer" <- function(depth=1)
  {
    if ( 0 > depth ) depth = 0 
    if ( depth > 5 ) depth = 4 else depth = depth-1
    
      c  <-  .C( "st_start_timer",
                 as.integer(depth),
                 as.double( proc.time()[3] ) )
  }

"st.stop.timer" <- function ( depth=1, info.level=0 )
{
    if ( 0 > depth ) depth = 0 
    if ( depth > 5 ) depth = 4 else depth = depth-1

    etime = 0
    c  <-  .C( "st_stop_timer",
               as.integer(depth),
               etime=as.double( etime ) )

    time <- proc.time()[3] - c$etime

    if (time > 3600){
      time <- time / 3600
      tunit <- "h"
    } else
    if (time > 60) {
      time <- time / 60
      tunit <- "m"
    } else
    {
      tunit <- "s"
    }
    st.info(info.level,paste("Elapsed time:",paste(round(time,digits=2),tunit,sep="")))
  }

#######

"st.read.marker.data" <- function ( filename,
                                    format,
                                    auxfiles  = FALSE,
                                    data.set  = "default" )
{
  st.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  dta.filename <- paste(st.input.dir,  ip, filename,sep="")
  ind.filename <- paste(st.output.dir, op ,filename,".ind",sep="")
  mar.filename <- paste(st.output.dir, op ,filename,".mar",sep="")
  gen.filename <- paste(st.output.dir, op ,filename,".gen",sep="")
  retval <- 0;

  if ("l" == format)
    {
      c <- .C("gs_read_marker_data_01_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }

    if ("m" == format)
    {
      c <- .C("gs_read_marker_data_02_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }

    if ("n" == format)
    {
      c <- .C("gs_read_marker_data_03_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }

    if ("t" == format)
    {
      c <- .C("gs_read_marker_data_04_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }
  
  outp <- NULL

  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      individual.list <- read.table(file=ind.filename,
                                    check.names = FALSE,
                                    header=TRUE,
                                    stringsAsFactors=F)
      outp$individual.list <- individual.list
      marker.list      <- read.table ( file=mar.filename,
                                       check.names = FALSE,
                                       header=TRUE,
                                       stringsAsFactors=F  )
      outp$marker.list <- marker.list
      genotypes        <- read.table ( file=gen.filename,
                                      check.names = FALSE,
                                      header=TRUE,
                                      stringsAsFactors=F)
      outp$genotypes   <- genotypes
      unlink(ind.filename)
      unlink(mar.filename)
      unlink(gen.filename)
    }
  
  st.stop.timer (info.level=1)
  return(invisible(outp))

}

"st.read.performance.data" <- function ( filename,
                                         out.filename = "y.vector",
                                         auxfiles     = TRUE,
                                         data.set     ="default"     )
{
  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  dta.filename <- paste ( st.input.dir,  ip, filename,sep="" )
  per.filename <- paste ( st.output.dir, op, out.filename,sep="" )
  unlink(per.filename)
  retval <- 0;
  c <- .C("gs_read_performance_data_01_GV",
          as.character(dta.filename),
          as.character(per.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  y <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      y <- read.table(file=per.filename,header=F,skip=1)
      unlink(per.filename)
      names(y) <- c("i","y")
    }
  gs.stop.timer (info.level=1)
  return( invisible(y) )
}

"st.return.performance.data" <- function ( out.filename = "y.vector",
                                           auxfiles     = TRUE,
                                           data.set     ="default"     )
{
  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  per.filename <- paste ( st.output.dir, op, out.filename,sep="" )
  unlink(per.filename)
  retval <- 0;
  c <- .C("gs_return_performance_data_01_GV",
          as.character(per.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  y <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      y <- read.table(file=per.filename,header=F,skip=1)  
      unlink(per.filename)
      names(y) <- c("i","y")
    }
  gs.stop.timer (info.level=1)
  return( invisible(y) )
}


"st.read.map" <- function ( filename,
                            skip = 0,
                            format     = "cpms", # or "mcp" 
                            m.stretch  = 1,
                            auxfiles   = TRUE,
                            data.set   = "default" )
{
  st.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  map.filename <- paste(st.input.dir,  ip, filename,sep="")
  out.filename <- paste(st.output.dir, op, filename,".srt",sep="")
  retval <- 0;
  c <- .C("gs_read_map_01_GV",
          as.character(map.filename),
          as.integer(skip),
          as.character(format),
          as.double (m.stretch),
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=F)
      unlink(out.filename)
      names(x) <- c("Chrom","Pos","Name","Class")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.get.map" <- function (  filename   = "map",
                            data.set   = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,filename,".srt",sep="")
  retval <- 0
  c <- .C("gs_get_map_01_GV",
          as.character(out.filename),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if    (-2 != c$retval ) 
    {
      x <-  read.table(file=out.filename,header=F,check.names=F)
      unlink(out.filename)
      names(x) <- c("Chrom","Pos","Name","Class")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.restrict.marker.data" <- function ( ind.list  = NULL,
                                        ind.file  = NULL,
                                        mar.list  = NULL,
                                        mar.file  = NULL,
                                        NoAll.MAX = 999,
                                        MaMis.MAX = 1,
                                        ExHet.MIN = 0,
                                        InMis.MAX = 1,
                                        data.set  = "default" )
{
  st.start.timer()

  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  
  if (  ! is.null (ind.file) && ! is.null (ind.list) ) {
      st.info(-2,"Provide either ind.file or ind.list")
      return (invisible(NULL));
    }

  if (  ! is.null (mar.file) && ! is.null (mar.list) ) {
      st.info(-2,"Provide either mar.file or mar.list")
      return (invisible(NULL));
    }
  
  ind.filename <- ""

  if ( ! is.null (ind.file) ) {
    ind.filename <- paste(st.input.dir,ip,ind.file,sep="")
  }

  if (  ! is.null( ind.list) ) {
      ind.file = paste("tmp-",round(runif(1,0,100000)),sep="") 
      ind.filename <- paste(st.input.dir,ip,ind.file,sep="")
      write.table (ind.list,file=ind.filename,quote=F,row.names=F,
                   col.names=F)
    }
  
  mar.filename <- ""
  
  if ( ! is.null( mar.file) ) {
    mar.filename <- paste(st.input.dir,ip,mar.file,sep="")
  }

  if (  ! is.null( mar.list) ) {
    mar.file = paste("tmp-",round(runif(1,0,100000)),sep="") 
    mar.filename <- paste(st.input.dir,ip,mar.file,sep="")
    write.table (mar.list,file=mar.filename,quote=F,row.names=F,
                 col.names=F)
  }
  
  retval <- 0;
  c <- .C("gs_restrict_marker_data_03_GV"    ,
          as.character(ind.filename),
          as.character(mar.filename),
          as.integer(NoAll.MAX),
          as.double (MaMis.MAX),
          as.double (ExHet.MIN),
          as.double (InMis.MAX),
          retval = as.integer(retval),
          as.character(data.set)       )
  
  if (  ! is.null( ind.list) ) unlink(ind.filename)
  if (  ! is.null( mar.list) ) unlink(mar.filename)

  st.stop.timer (info.level=1)

}

"st.copy.marker.data" <- function (  target.data.set,
                                     source.data.set )
{

  retval <- 0;
  c <- .C("gs_copy_marker_data_01_GV"   ,
          retval = as.integer(retval)   ,
          as.character(target.data.set) ,
          as.character(source.data.set)    )
}



"st.write.marker.data" <- function ( format     = "m",
                                     lfilename  = "tmp-l",
                                     mfilename  = "tmp-m",
                                     nfilename  = "tmp-n",
                                     ifilename  = "",
                                     f.ind      = 0,
                                     l.ind      = 0,
                                     add.inum   = ".1",
                                     auxfiles   = FALSE,
                                     data.set   = "default" )
{
  st.start.timer()
  #
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  # Data in list format #
  lpo.filename <- paste(st.output.dir,op,lfilename,".lpo",sep="")
  # Map for graphical genotypes #
  mmp.filename <- paste(st.output.dir,op,mfilename,".mmp",sep="")
  # Data matrix for graphical genotypes #
  mpo.filename <- paste(st.output.dir,op,mfilename,".mpo",sep="")
  # Genome paramteters for Plabsim #
  ngp.filename <- paste(st.output.dir,op,nfilename,".ngp",sep="")
  # Map for Plabsim #
  nmp.filename <- paste(st.output.dir,op,nfilename,".nmp",sep="")
  # Data matrix for Plabsim #
  npo.filename <- paste(st.output.dir,op,nfilename,".npo",sep="")
  # Individual list for sorting and selecting #
  ind.filename <- ""
  if ( "" != ifilename)
  ind.filename <- paste(st.output.dir,op,ifilename,sep="")
  
  retval <- 0;
  c <- .C("gs_write_marker_data_02_GV",
          as.character(format),
          as.character(lpo.filename),
          as.character(mmp.filename),
          as.character(mpo.filename),
          as.character(ngp.filename),
          as.character(nmp.filename),
          as.character(npo.filename),
          as.character(ind.filename),
          as.integer(f.ind),
          as.integer(l.ind),
          as.character(add.inum),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )

  outp <- NULL
  st.stop.timer (info.level=1)
  return( invisible(outp) )

}

"st.write.map" <- function (         mfilename  = "tmp-m",
                                     auxfiles   = TRUE,
                                     data.set   = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  mmp.filename <- paste(st.output.dir,op,mfilename,".mmp",sep="")

  retval <- 0;
  c <- .C( "gs_write_map_01_GV",
           as.character(mmp.filename),
           retval = as.integer(retval),
           as.character(data.set)          )

  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=mmp.filename,header=F)
      names(x) <- c("Chrom","Pos","Name","Class")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}



"st.recode.ref" <- function ( reference  = 0,
                              missing    = -1,
                              descending = FALSE,
                              data.set   = "default" )
{
  st.start.timer()
  retval <- 0;
  c <- .C("gs_recode_ref_01_GV",
          as.integer(reference),
          as.integer(missing),
          as.integer(descending),
          retval = as.integer(retval),
          as.character(data.set)
          )
  st.stop.timer (info.level=1)

}

"st.def.hblocks" <- function ( hap          = 1,
                               hap.unit     = 1,
                               hap.symbol   = "b",
                               out.filename = "blocks"  ,
                               auxfiles     = TRUE,
                               data.set     = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,out.filename,".hmap",sep="")
  retval <- 0;
  c <- .C("gs_def_hblocks_01_GV",
          as.integer(hap),
          as.integer(hap.unit),
          as.character(hap.symbol),
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=F)
      unlink(out.filename)
      names(x) <- c("Chrom","Pos","Name","Class","Markers")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.recode.hbc" <- function ( reference  = 0,
                             data.set   = "default" )
{
  st.start.timer()
  retval <- 0;
  c <- .C("gs_recode_hbc_01_GV",
          as.integer(reference),
          retval = as.integer(retval),
          as.character(data.set)
          )
  st.stop.timer (info.level=1)
}

"st.recode.hil" <- function ( data.set   = "default" )
{
  st.start.timer()
  retval <- 0;
  c <- .C("gs_recode_hil_01_GV",
          retval = as.integer(retval),
          as.character(data.set)
          )
  st.stop.timer (info.level=1)
}

"st.plot.ggt.src" <- function (
                        map   ,
                        pop   ,
                        f.lr  , 
                        f.tb  ,    
                        d.h   ,    
                        d.v   ,    
                        p.t   ,   
                        p.s   ,       
                        d.t   ,   
                        d.map ,   
                        c.nme ,   
                        i.nme ,   
                        color ,
                        z.min ,
                        z.max ,
                        cex.chrom,
                        cex.ind  ,
                        plt.png   ,
                        plt.pdf   ,
                        plt.name  ,
                        plt.width ,
                        plt.height,
                        plt.ptsize
                       )
{
  m <- map
  g <- pop
  #
  #
  close.screen(all=TRUE)
  
   if ( TRUE == plt.pdf )
     {
       pdf( file=plt.name, width=plt.width, height=plt.height, pointsize=plt.ptsize)
     }
   else if ( TRUE == plt.png )
     {
       png( file=plt.name, width=plt.width, height=plt.height, pointsize=plt.ptsize,
            res=300, units = "in")
     }
   else {
       # dev.new()
     }
  par (mar=c(0,0,0,0) )
  #
  # Number of individuals
  n.ind <- ncol(g)
  #
  clr <- color
  # Number of chromosomes
  n.chrom <- max(m[,1])
  # Chromosome names
  if ( ( 1 == length(c.nme) ) && ( "" == c.nme[1] ) ) {
    c.nme <- paste("Chr. ",1:n.chrom)
  }
  # Individual names
  if ( ( 1 == length(i.nme) ) && ( "" == i.nme[1] ) ) {
    i.nme <- names(g)
  }
  # Length of chromosomes
  {
    m.chrom <- 1:n.chrom
    l.chrom <- 1:n.chrom
    b.chrom <- vector("list",n.chrom)
    p.chrom <- vector("list",n.chrom)
    for (c in 1:n.chrom) {
      idx <- m[,1] == c
      m.chrom[c] <- sum(idx)
      m.pos <- m[idx,][,2]
      l.chrom[c] <- max(m.pos) - min(m.pos)
      m.brk <- 0:length(m.pos)
      for ( j in 2:length(m.pos) ) {
        m.brk[j] <- m.pos[j-1] + (m.pos[j] - m.pos[j-1]) / 2
      }
      m.brk[1+length(m.pos)] <- m.pos[length(m.pos)]
      b.chrom[[c]] <- m.brk
      p.chrom[[c]] <- m.pos
    }
  }
  # Width of displays
  ndc  <- 1 + n.chrom
  f.disp <- c(f.lr*mean(l.chrom),l.chrom)
  f.disp <- f.disp / sum (f.disp)
  s.disp <- f.disp
  for (kk in 2:ndc) s.disp[kk] <- s.disp[kk-1] +s.disp[kk]
  s.disp <- c(0,s.disp)
  #
  dcol <- matrix(0,ncol=2,nrow=ndc)
  for (kk in 1:ndc)
    {
      dcol[kk,1] = s.disp[kk];
      dcol[kk,2] = s.disp[kk+1];
    }

  # Height of displays
  ndr <- 2 + n.ind
  f.disp <- c ( f.tb, rep(1,n.ind) , f.tb)
  f.disp <- f.disp / sum (f.disp)
  s.disp <- f.disp
  for (kk in 2:ndr) s.disp[kk] <- s.disp[kk-1] +s.disp[kk]
  s.disp <- c(0,s.disp)
  #
  drow <- matrix(0,ncol=2,nrow=ndr)
  for (kk in 1:ndr)
    {
      drow[kk,1] = s.disp[kk];
      drow[kk,2] = s.disp[kk+1];
    }

  # Define the figures
  figs <- matrix (ncol=4,nrow=ndr*ndc )
  for (i in 1:ndr) for (c in 1:ndc)
    {
      idx  <- (i-1)*ndc + c
      figs[idx,] <- c(dcol[c,],drow[i,]  )
    }
  #
  sve      <- 1 - figs[,3]
  figs[,3] <- 1 - figs[,4]
  figs[,4] <- sve
  # 
  figs[,1] <- figs[,1] + d.h
  figs[,2] <- figs[,2] - d.h
  figs[,3] <- figs[,3] + d.v
  figs[,4] <- figs[,4] - d.v
  #
  figs <- round(figs,digits=5)
  #
  neu <- split.screen ( figs  )
  #
  # Plot chromosome.names and legend
  for (c in 1:n.chrom)
    {
      scr  <- (c+1)
      screen(scr)
      plot(0,0,t="n",axes=F, xlab="", ylab="")
      text (0,0,c.nme[c],cex=cex.chrom)
#
      scr <- (ndr-1) * ndc + (c+1)
      screen(scr)
      mi <- min(b.chrom[[c]])
      ma <- max(b.chrom[[c]])    
     if ( ( p.t == T) || ( p.s == T) )
       {
         plot(0,0, axes=F,t="n", xlab="", ylab="", xlim=c(mi,ma),ylim=c(-d.t,d.t) )
         lines( c(mi,ma),c(0,0))
       }
      if ( p.t == T)
        {
          for (mm in 1:length(p.chrom[[c]]) )
            lines (c(p.chrom[[c]][mm],p.chrom[[c]][mm]),c(0,d.t/2))
        }
      if ( p.s == T)
        {
          t <- seq(0,ma,d.map)
          text( t,rep(-d.t/2,length(t)),t  )
        }
      
     }
  # Plot individual labels
  for (i in 1:n.ind)
    {
      # Activate screen
      scr  <- (i)*ndc + (1)
      screen(scr)
      plot(0,0,t="n",axes=F, xlab="", ylab="")
      text (0,0,i.nme[i],cex=cex.ind)
    }
  # Plot all chromosomes
  for (i in 1:n.ind) for (c in 1:n.chrom)
    {
      # Activate screen
      scr  <- (i)*ndc + (c+1)
      screen(scr)
     #    
      p.m <- matrix( 0 ,
                   nrow=m.chrom[c] , ncol=2 , 
                   byrow=F                    )
      #
      idx <- m[,1] == c
      m.nme <- as.character ( m[idx,][,3] )
      idx2 <- row.names(g) %in% m.nme
      p.g <- as.character(g[idx2,i])
      #
      splt <- strsplit(p.g,"/")
      #
      for (p in 1:m.chrom[c]) {
        p.m[p,1] <- as.numeric(splt[[p]][1]) 
        p.m[p,2] <- as.numeric(splt[[p]][2]) 
      }
      image(
            x = b.chrom[[c]]    ,
            y = c(0,1,2)        ,
            z = p.m             ,
            axes = 0            ,
            zlim=c(z.min,z.max), 
            col = clr        
            )
      box()
      #
    }
  close.screen(all=TRUE)
 if (( TRUE == plt.pdf ) | ( TRUE == plt.png )) dev.off() 
}


"st.plot.ggt" <- function (
                        data.set   = "default",
                        ifilename  = "",
                        i.list     = "",
                        f.ind      = 0,
                        l.ind      = 0,
                        f.lr  = 2,      # Width of left description
                        f.tb  = 2,      # Hight of top/bottom description
                        d.h   = 0.001,  # Distance between chromosomes
                        d.v   = 0.001,  # Distance between individuals
                        p.t   = F    ,  # Plot marker tickmarks
                        p.s   = F    ,  # Plot cm scale
                        d.t   = 50   ,  # Lengh of marker tickmarks
                        d.map = 100  ,  # Distance between marker labels
                        c.nme = ""   ,  # Chromosome names
                        i.nme = ""   ,  # Individual names
                        cex.chrom = 1,  # Scaling factor chromosome names
                        cex.ind   = 1,  # Scaling factor individual names
                        color = c ("yellow","blue","red","green",
                              "wheat", "skyblue","tomato","palegreen",
                              "yellow4","darkblue","darkblue","darkgreen" ),
                        z.min        = 0,
                        z.max        = 12,
                        plt.pdf      = F,
                        plt.png      = F,
                        plt.fname    = "ggt",
                        plt.width    = 10,
                        plt.height   = 10,
                        plt.ptsize   = 10
                       )
{
  st.start.timer()
#
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"

  mfilename = paste("tmp-",round(runif(1,0,100000)),sep="") 
#
   mmp.filename  = st.outdir(paste(mfilename,".mmp",sep=""))
   mpo.filename  = st.outdir(paste(mfilename,".mpo",sep=""))

  if (T == plt.pdf )
    out.filename = paste(st.output.dir,op,plt.fname,".pdf",sep="")
  else if (T == plt.png )
    out.filename = paste(st.output.dir,op,plt.fname,".png",sep="")
#
  if ( "" != i.list[1] )  {
    ifilename = "tmp-inds"
    ind.filename <- paste(st.output.dir,op,ifilename,sep="")
    write.table ( i.list,
                  file=ind.filename,
                  quote=F,
                  row.names = F, col.names = F)
  } 

  if ( ( "" == i.list[1] ) && ( "" != ifilename )) {
    i.list <- read.table( paste(st.input.dir,ip,ifilename,sep=""),
                          stringsAsFactors=F)[,1]
    ind.filename <- paste(st.output.dir,op,ifilename,sep="")
    write.table ( i.list,
                  file=ind.filename,
                  quote=F,
                  row.names = F, col.names = F)
  } 
  
  st.write.marker.data ( format     = "m",
                         mfilename  = mfilename,
                         ifilename  = ifilename,
                         f.ind      = f.ind,
                         l.ind      = l.ind,
                         data.set   = data.set  )
                                        #
  if ( !file.exists(mmp.filename))
    {
      st.info(-2,"No linkage map loaded");
      return ;
    }

  retval <- 0;

  m <- read.table(mmp.filename,stringsAsFactors=F, check.names = FALSE,)
  p <- read.table(mpo.filename, check.names = FALSE)

  unlink(mmp.filename)
  unlink(mpo.filename)
  
  st.plot.ggt.src  (
                        map = m , 
                        pop = p , 
                        f.lr,  
                        f.tb,  
                        d.h,   
                        d.v,   
                        p.t,   
                        p.s,   
                        d.t,   
                        d.map, 
                        c.nme, 
                        i.nme, 
                        color,
                        z.min,
                        z.max,
                        cex.chrom,
                        cex.ind  ,
                        plt.png,
                        plt.pdf,
                        plt.name=out.filename,
                        plt.width,
                        plt.height,
                        plt.ptsize 
               )
     if (( TRUE == plt.pdf ) | ( TRUE == plt.png )) st.info(1,paste("GGT file: ",out.filename,sep=""))
  st.stop.timer (info.level=1)
}


"st.marker.data.statistics" <- function ( data.set ="default",
                                          filename ="marker.stats",
                                          ind = T,
                                          mar = T,
                                          gen = T,
                                          auxfiles = TRUE )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  ind.filename <- paste(st.output.dir,op,filename,".ind",sep="")
  mar.filename <- paste(st.output.dir,op,filename,".mar",sep="")
  gen.filename <- paste(st.output.dir,op,filename,".gen",sep="")
  retval <- 0;
  c <- .C("gs_marker_stats_01_GV",
          as.character(ind.filename),
          as.character(mar.filename),
          as.character(gen.filename),
          as.integer(ind),
          as.integer(mar),
          as.integer(gen),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )

  outp <- NULL

  if  ( (auxfiles) && (-2 != c$retval ) )
    {
       if (ind) {
        individual.list <- read.table ( file=ind.filename,
                                        header=TRUE,
                                        check.names=F,                                        
                                        stringsAsFactors=F)
        outp$individual.list <- individual.list
      }
      if (mar) {
        marker.list      <- read.table ( file=mar.filename,
                                         check.names=F,                                        
                                         header=TRUE,
                                         stringsAsFactors=F  )
      outp$marker.list <- marker.list
      }
      if (gen) {
        genotypes        <- read.table ( file=gen.filename,
                                         check.names=F,                                        
                                         header=TRUE,
                                         stringsAsFactors=F)
        outp$genotypes   <- genotypes
      }
      unlink(ind.filename)
      unlink(mar.filename)
      unlink(gen.filename)
     }
  
  st.stop.timer (info.level=1)
  return(invisible(outp))

}


st.plot.gene.diversity <- function (
   data.set ="default",
   plt.pdf      = F,
   plt.fname    = "gene.div",
   plt.width    = 10,
   plt.height   = 6,
   plt.ptsize   = 14
)
{
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  nfilename = paste("tmp-",round(runif(1,0,100000)),sep="") 
#
   smp.file = st.outdir(paste(nfilename,".nmp",sep=""))
   spo.file = st.outdir(paste(nfilename,".npo",sep=""))
   sgp.file = st.outdir(paste(nfilename,".ngp",sep=""))
   out.filename = paste(st.output.dir,op,plt.fname,".pdf",sep="")
#
#
  st.write.marker.data ( format     = "n",
                         nfilename  = nfilename,
                         data.set   = data.set  )
#
  if ( !file.exists(sgp.file) )
    {
      st.info(-2,"No linkage map loaded");
      return ;
    }
#
   reset.all()
   cl <- scan(sgp.file,quiet=T)
   set.genome.par( no.chrom=length(cl), no.hom=2, chrom.len=cl )
   linkage.map.load(smp.file)
   population.matrix.load(spo.file)
#
  unlink(smp.file)
  unlink(spo.file)
  unlink(sgp.file)
#
   l.names <- as.character(list.populations()$PopName)
   p.name <- "pop"
   remove.population(p.name)
   for (i in 1:length(l.names))
      append.population(p.name,l.names[i])
#
    mp <- get.map()
    gene.div <- NULL
    for (i in 1:nrow(mp)){
      af <- evaluate.allele(p.name,as.character(mp[i,]$name))
      gd <- 1 - sum( ( af$count/sum(af$count) )^2 )
      gene.div <- c(gene.div,gd)
    }
#
  mp$gene.div <- gene.div
#
  c.pos <- cl
  for ( i in 2:length(c.pos) ) c.pos[i] <- c.pos[i] + c.pos[i-1]
#
  mp$cpos <- mp$pos 
  for ( i in 1:nrow(mp) )
     if (mp$chrom[i] > 1)
        mp$cpos[i] <- mp$cpos[i] + c.pos[ mp$chrom[i]-1]
#
  if (T == plt.pdf)
    {
      pdf(
          file  = out.filename,
          width = plt.width, 
          height=plt.height, 
          pointsize=plt.ptsize)
    }
  else
   {
      # dev.new()
    }
  # Plot
  par(mar=c(4,4,1,1))
  plot(
     mp$cpos,
     mp$gene.div,
     type="l",
     lwd=2,
     col="blue",
     ylim=c(0,1),
     xaxt="n",
     xlab="Chromosome",
     ylab="Gene diversity"
  )
#
  c.pos <- c(0,c.pos)
  axis( 1,
       at=c.pos,
       labels=c(1:length(cl)," "),
       tck=-.01)
#
    # Chromosome borders
    for (i in 1:(length(c.pos))){
      x <- c.pos[i]
      y0 <- 0
      y1  <- 1
      lines(c(x,x),c(y0,y1))
    }
#
  # Gene diversity per chromosome
  gdiv.chr <- tapply(mp$gene.div,mp$chrom,mean)
   for (i in 1: (length(c.pos)-1) ){
    text((c.pos[i]+c.pos[i+1])/2,0.95, sprintf("%4.2f",gdiv.chr[i]))
  }
   if (T == plt.pdf) dev.off()

  ret.val <- data.frame( chrom = mp$chrom,
                         pos   = mp$pos,
                         name  = mp$name,
                         cpos  = mp$cpos,
                         gene.div=mp$gene.div);
  
  return(invisible(ret.val))
  
 }

"gs.return.effects" <- function( out.filename ="effects"    ,
                                 data.set     ="default" )
{
  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_return_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
  if  (-2 == c$retval ) return(invisible(NULL))
  x <- read.table(file=x.filename,header=T,stringsAsFactors=FALSE)
  unlink(x.filename)
  return (x)
}

"gs.set.effects" <- function( eff,
                              filename ="seteffects" ,
                              data.set ="default"  )
{
  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,filename,sep="")
  write.table(eff,file=x.filename,quote=FALSE,row.names=FALSE)
  retval <- 0;
  c <- .C( "gs_set_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
}

"gs.return.pvals" <- function( out.filename ="effects"    ,
                               data.set     ="default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_return_pvals_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
  if  (-2 == c$retval ) return(invisible(NULL))
  x <- read.table(file=x.filename,header=T,stringsAsFactors=FALSE)
  unlink(x.filename)
  return(x)
}

"gs.check.pvals" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C( "gs_check_pvals_01_GV"      ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
  return( as.logical(c$retval) )
}

gs.plot.effects <- function (data.set = "default",
                             absolute = TRUE,
                             alpha    = 0.05,
                             p.adjust.method = "none",
                             pvals = NULL,
                             xlim  = NULL,
                             ylim  = NULL,
                             plt   = TRUE,
                             col.sp = "blue", pch.sp = 19 ,
                             col.sn = "blue", pch.sn = 19  , 
                             col.ns = "blue", pch.ns = 1 , 
                             ... ) 
{
  gs.start.timer()
  {
    if ( "" == st.input.dir )  ip <- "" else ip <- "/"
    if ( "" == st.output.dir ) op <- "" else op <- "/"
    tmp.filename = paste("tmp-", round(runif(1, 0, 1e+05)), sep = "")
    auxfiles = TRUE
    x.filename <- paste(st.output.dir, op, tmp.filename, sep = "")
    retval <- 0
    pvals.av <- gs.check.pvals(data.set)
    if ( TRUE == pvals.av )
      {
        c <- .C("gs_return_pvals_01_GV", as.character(x.filename), 
                retval = as.integer(retval), as.character(data.set))
      }
    else
      {
        c <- .C("gs_return_effects_01_GV", as.character(x.filename), 
                retval = as.integer(retval), as.character(data.set))
      }
    if (-2 == c$retval)  return(invisible(NULL))
    x <- read.table(file = x.filename, header = T, 
                                check.names = F)
    if (TRUE == pvals.av)
      { x <- data.matrix(x)
        effect <- x[2:nrow(x),1]
        pvalue <- x[2:nrow(x),2]
      } else {
        effect <- x$effect[2:nrow(x)]
      }
    if (!is.null(pvals) )
      { 
        pvalue <- pvals
        pvals.av <- TRUE
      }
    unlink(x.filename)
    {
      nfilename = paste("tmp-", round(runif(1,0,1e+05)),sep = "")
      smp.file = st.outdir(paste(nfilename, ".nmp", sep = ""))
      spo.file = st.outdir(paste(nfilename, ".npo", sep = ""))
      sgp.file = st.outdir(paste(nfilename, ".ngp", sep = ""))
      st.write.marker.data(format = "n", nfilename = nfilename, 
                           data.set = data.set)
      if (!file.exists(sgp.file)) {
        st.info(-2, "No linkage map loaded")
        return (invisible(NULL))
      }
      reset.all()
      cl <- scan(sgp.file, quiet = T)
      set.genome.par(no.chrom = length(cl), no.hom = 2, chrom.len = cl)
      linkage.map.load(smp.file)
      unlink(smp.file)
      unlink(spo.file)
      unlink(sgp.file)
      mp <- get.map()
    }
    if (TRUE==absolute) mp$effect <- abs(effect) else  mp$effect <- effect
    if (TRUE == pvals.av) {
      mp$pvalue   <- pvalue
      mp$pvalue.a <- p.adjust(mp$pvalue,p.adjust.method) 
    }
    c.pos <- cl
    for (i in 2:length(c.pos)) c.pos[i] <- c.pos[i] + c.pos[i-1]
    for (i in 1:nrow(mp)) if (mp$chrom[i] > 1) 
      mp$pos[i] <- mp$pos[i] + c.pos[mp$chrom[i] - 1]
  }
  if (TRUE == plt)
      {
        par(mar = c(4, 4, 1, 1))
        if ( is.null(ylim) )
          ylim = c( min( c(min(mp$effect),0) ), max(mp$effect) )
        if ( is.null(xlim) )
          xlim = c( min( c(min(mp$pos),0) ), max(mp$pos) )
        plot(mp$pos, mp$effect, type = "n", 
             xaxt = "n", xlab = "Chromosome", ylab = "Effect size",
             ylim=ylim, xlim=xlim,...)
        c.pos <- c(0, c.pos)
        axis(1, at = c.pos, labels = c(1:length(cl), " "), tck = -0.01)
        grid()
        for (i in 1:(length(c.pos))) {
          x <- c.pos[i]
          lines( c(x,x), ylim )
        }
        if (F == absolute) lines ( xlim,c(0,0))
        tol = 1e-8
        if ( TRUE == pvals.av )
          {
            mp1   <- subset ( mp, pvalue.a <  alpha )
            mp1ap <- subset ( mp1, effect >   tol )
            mp1an <- subset ( mp1, effect <  -tol )
            mp1b  <- subset ( mp1, abs(effect) <= tol )
            mp2   <- subset ( mp, pvalue.a >= alpha )
            points(mp1b$pos,  mp1b$effect, type = "p", lwd = 2, col=col.ns, pch=pch.ns)
            points(mp2$pos,   mp2$effect,  type = "p", lwd = 2, col=col.ns, pch=pch.ns)
            points(mp1ap$pos, mp1ap$effect, type = "p", lwd = 2, col=col.sp, pch=pch.sp)
            points(mp1an$pos, mp1an$effect, type = "p", lwd = 2, col=col.sn, pch=pch.sn)
          }
        else
          {
            points(mp$pos, mp$effect, type = "p", lwd = 2, col = "blue")
          }
      }
    ret.val <- data.frame(chrom = mp$chrom, pos = mp$pos, name = mp$name, 
        cpos = mp$pos, effect = mp$effect)
    if ( TRUE == pvals.av ) ret.val$pvalue <- mp$pvalue
    gs.stop.timer(info.level = 1)
    return(invisible(ret.val))
}

gs.compare.effects <- function ( data.set.a,
                                 data.set.b,
                                 label.a      = data.set.a,
                                 label.b      = data.set.b   )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  tmp.filename = paste("tmp-",round(runif(1,0,100000)),sep="") 
  auxfiles     = TRUE
  x.filename <- paste(st.output.dir,op,tmp.filename,sep="")

  retval <- 0;
  c <- .C( "gs_return_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set.a)      )

  if (-2 == c$retval) return (invisible(NULL))
  
  x <- read.table(file=x.filename,header=T) 
  effects.a <- x$effect[2:nrow(x)]

  retval <- 0;
  c <- .C( "gs_return_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set.b)      )

  if (-2 == c$retval) return (invisible(NULL))
  
  x <- read.table(file=x.filename,header=T) 
  effects.b <- x$effect[2:nrow(x)]
  
  unlink(x.filename)

  max <- max ( c( abs(effects.a),abs(effects.b) ) ) 
  min <- -max;
  
  par(mar=c(4,4,1,1))
  plot( effects.b ~effects.a ,
        type="p",
        lwd=2,
        col="blue",
        ylim = c(min,max),
        xlim = c(min,max),
        xlab=label.a,
        ylab=label.b
  )
  lines (c(min,max),c(min,max))
#
  gs.stop.timer (info.level=1)
  return (invisible(NULL))
  
 }

gs.plot.model.fit <- function ( training.set ="default",
                                title        = ""         )
{
  c1 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = training.set )
  if ( is.null(c1) ) return (invisible(NULL) )
  #
  c1 <- c1[,2:3]
  #
  par(mar=c(4,4,1,1) )
  #
  l <- c(min(c1),max(c1))
  p1 <- l[1] + 0.65*(l[2]-l[1] ) ; p2 <- l[1] + 0.10*(l[2]-l[1] )
  p3 <- l[1] + 0.05*(l[2]-l[1] ) ; p4 <- l[1] + 0.90*(l[2]-l[1] )
  #
  plot(c1,ylim=l,xlim=l,col="blue")
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(c1)[2],
                   cor(c1,method="spearman")[2]),pos=4)
  text(p3,p4,sprintf("%s\nModel fit",title),pos=4)
  lines(l,l)
  lines ( c(mean(c1[,1]),mean(c1[,1])), c(min(c1[,2]),max(c1[,2])) )
  lines ( c(min(c1[,1]),max(c1[,1])), c(mean(c1[,2]),mean(c1[,2])) )
}

gs.plot.validation <- function ( training.set ,
                                 prediction.set ,
                                 title=""          )
{
  #
  c1 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = training.set )
  c2 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = prediction.set )
  #
  if ( 3 == ncol(c1) ) c1 <- c1[,2:3] else return (invisible(NULL) )
  if ( 3 == ncol(c2) ) c2 <- c2[,2:3] else
   {
     st.info (-2,"No phenotypic data in the prediction set")
     return (invisible(NULL) )
   }
  #
  if ( is.null(c1) || is.null(c1) ) return (invisible(NULL) )
  #
  par( mfcol=c(1,2) , mar=c(4,4,1,1) )
  #
  l <- c( min(c1,c2),max(c1,c2) )
  #
  p1 <- l[1] + 0.65*(l[2]-l[1] ) ; p2 <- l[1] + 0.10*(l[2]-l[1] )
  p3 <- l[1] + 0.05*(l[2]-l[1] ) ; p4 <- l[1] + 0.90*(l[2]-l[1] )
  #
  plot(c1,ylim=l,xlim=l,col="blue")
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(c1)[2],
                     cor(c1,method="spearman")[2]),pos=4)
  text(p3,p4,sprintf("%s\nModel fit",title),pos=4)
  lines(l,l)
  lines ( c(mean(c1[,1]),mean(c1[,1])), c(min(c1[,2]),max(c1[,2])) )
  lines ( c(min(c1[,1]),max(c1[,1])), c(mean(c1[,2]),mean(c1[,2])) )
  #
  plot(c2,ylim=l,xlim=l,col="blue")
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(c2)[2],
                     cor(c2,method="spearman")[2]),pos=4)
  text(p3,p4,sprintf("%s\nValidation",title),pos=4)
  lines(l,l)
  lines ( c(mean(c2[,1]),mean(c2[,1])), c(min(c2[,2]),max(c2[,2])) )
  lines ( c(min(c2[,1]),max(c2[,1])), c(mean(c2[,2]),mean(c2[,2])) )
}

st.plot.corr <- function ( x , y,
                          title        = "" ,
                          pch          = 16 ,
                          cex          = 1  ,
                          xlab         = "x",
                          ylab         = "y"  )
{
  par ( mar=c(4,4,1,1) )
  #
  lx <- c(min(x),max(x))
  ly <- c(min(y),max(y))
  p1 <- lx[1] + 0.65*(lx[2]-lx[1] ) ; p2 <- ly[1] + 0.10*(ly[2]-ly[1] )
  p3 <- lx[1] + 0.05*(lx[2]-lx[1] ) ; p4 <- ly[1] + 0.90*(ly[2]-ly[1] )
  #
  plot(x,y,xlab=xlab,ylab=ylab,type="n")
  points(x,y,pch=pch,cex=cex);
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(x,y),
                   cor(x,y,method="spearman")),pos=4)
  
  text(p3,p4,sprintf("%s",title),pos=4)
  lines ( c(lx[1],lx[2]) , c(ly[1],ly[2]))
  lines ( c(mean(x),mean(x)), c(min(y),max(y)) )
  lines ( c(min(x),max(x)), c(mean(y),mean(y)) )
}

st.plot.corr.l <- function ( x , y,
                          title        = "" ,
                          pch          = 16 ,
                          cex          = 1  ,
                          xlab         = "x",
                          ylab         = "y",
                          ylim,
                          xlim )
{
  par ( mar=c(4,4,1,1) )
  #
  lx <- c(min(x),max(x))
  ly <- c(min(y),max(y))
  p1 <- lx[1] + 0.65*(lx[2]-lx[1] ) ; p2 <- ly[1] + 0.10*(ly[2]-ly[1] )
  p3 <- lx[1] + 0.05*(lx[2]-lx[1] ) ; p4 <- ly[1] + 0.90*(ly[2]-ly[1] )
  #
  plot(x,y,xlab=xlab,ylab=ylab,type="n",xlim=xlim,ylim=ylim)
  points(x,y,pch=pch,cex=cex);
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(x,y),
                   cor(x,y,method="spearman")),pos=4)
  
  text(p3,p4,sprintf("%s",title),pos=4)
  lines ( c(lx[1],lx[2]) , c(ly[1],ly[2]))
  lines ( c(mean(x),mean(x)), c(min(y),max(y)) )
  lines ( c(min(x),max(x)), c(mean(y),mean(y)) )
}


st.genetic.distances.02 <- function ( measure = "mrd", # "mrd" "rd" "euc"
                                      format  = "l",   # "l" "m"
                                      filename  = "genetic.distances",
                                      auxfiles  = FALSE,
                                      data.set  = "default" )
{
   nfilename = paste("tmp-",round(runif(1,0,100000)),sep="") 
#
   smp.file = st.outdir(paste(nfilename,".nmp",sep=""))
   spo.file = st.outdir(paste(nfilename,".npo",sep=""))
   sgp.file = st.outdir(paste(nfilename,".ngp",sep=""))
   if ( "" == st.output.dir ) op <- "" else op <- "/"
   out.filename = paste(st.output.dir,op,filename,".gdi",sep="")
#
  st.write.marker.data ( format     = "n",
                         nfilename  = nfilename,
                         data.set   = data.set  )
  dta  <- read.table(spo.file)
#
  unlink(smp.file)
  unlink(spo.file)
  unlink(sgp.file)
#   
  distance  <- gd.genetic.distance (dta,measure=measure) 
#
  if  (auxfiles)
    {
      write.table ( distance,
                    file=out.filename,
                    quote = F,
                    row.names = F,
                    col.names = F        )
    }
  
  st.stop.timer (info.level=1)
   if ("m" == format)
     {
       p <- gd.data.parameters(dta)
       aa <- distance[,3]
       bb <- p$no.pop
       cc <- p$pop.list
       dd <- structure(aa, Size = bb, Labels = cc, Diag = FALSE,
                       Upper = T, method = "euclidean", class = "dist")
       return ( invisible(dd) )
     }
   else
     {
       return ( invisible(distance) )
     }
}

st.genetic.distances <- function ( measure = "mrd", # "mrd" "rd" "euc"
                                   format  = "l",   # "l" "m"
                                   filename  = "genetic.distances",
                                   auxfiles  = FALSE,
                                   data.set  = "default" )
{
   if ( "" == st.output.dir ) op <- "" else op <- "/"
   out.filename = paste(st.output.dir,op,filename,".gdi",sep="")
#   
   gs.cross.eval.gd ( dist=measure , data.set=data.set )
   c <- gs.cross.info ( data.set=data.set )
   distance <- data.frame (OTU1=c$P1Name,OTU2=c$P2Name,Measure=c$gd)

#
  if  (auxfiles)
    {
      write.table ( distance,
                    file=out.filename,
                    quote = F,
                    row.names = F,
                    col.names = F        )
    }
  
  st.stop.timer (info.level=1)
   if ("m" == format)
     {
       c <- st.marker.data.statistics ( ind=T,
                                        mar=F,
                                        gen=F,
                                        auxfiles=T,
                                        data.set=data.set)
       n <- c$individual.list$Name
       aa <- distance[,3]
       bb <- length(n)
       cc <- n
       dd <- structure(aa, Size = bb, Labels = cc, Diag = FALSE,
                       Upper = T, method = "euclidean", class = "dist")
       return ( invisible(dd) )
     }
   else
     {
       return ( invisible(distance) )
     }
}


"st.indir" <- function ( filename )
{
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  paste(st.input.dir,ip,filename,sep="")
}

"st.outdir" <- function ( filename )
{
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  paste(st.output.dir,op,filename,sep="")
}

"st.datadir" <- function ( filename )
{
   if ( "" == st.data.dir ) dp <- "" else dp <- "/"
   paste(st.data.dir,dp,filename,sep="")
}

"st.scriptdir" <- function ( filename )
{
   if ( "" == st.script.dir ) sp <- "" else sp <- "/"
   paste(st.script.dir,sp,filename,sep="")
}

"st.id" <- function ( filename )
{
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  paste(st.input.dir,ip,filename,sep="")
}

"st.od" <- function ( filename )
{
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  paste(st.output.dir,op,filename,sep="")
}

"st.dd" <- function ( filename )
{
   if ( "" == st.data.dir ) dp <- "" else dp <- "/"
   paste(st.data.dir,dp,filename,sep="")
}

"st.sd" <- function ( filename )
{
   if ( "" == st.script.dir ) sp <- "" else sp <- "/"
   paste(st.script.dir,sp,filename,sep="")
}


"gs.cross.eval.gd" <- function ( dist     = "rd",
                                 data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_gd_01_GV"       ,
          as.character(dist)          ,   
          retval = as.integer(retval) ,
          as.character(data.set)       )
  gs.stop.timer (info.level=1)
}


"gs.cross.info.gd" <- function( out.filename ="cross.info.gd",
                                data.set     ="default"        )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_cross_info_gd_01_GV"    ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)        )
  x <- NULL
  if  (-2 != c$retval ) {
    x <-  read.table(file=x.filename,header=T,stringsAsFactors=F)
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.cross.eval.mu" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_mu_01_GV"    ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.mi" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_mi_01_GV"    ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.ma" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_ma_01_GV"    ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}


"gs.cross.eval.va" <- function ( data.set     = "default"  ,
                                 pop.type     = "unlinked" ,
                                 t            = 0          ,
                                 map.function = "Haldane"    )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_va_01_GV"    ,
          as.character(pop.type)      ,
          as.integer  (t)             ,
          as.character(map.function)  ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.es" <- function ( data.set = "default",
                                 alpha    = 0.1      ,
                                 N        = 0        ,
                                 G        = 0          )
{
  gs.start.timer()

  if ( (0==N) || (0==G) ) {
    i <- dnorm(qnorm(1-alpha,0,1),0,1)/alpha
  } else {
    alpha   <- N/G
    i.alpha <- dnorm(qnorm(1-alpha,0,1),0,1)/alpha
    i       <- i.alpha - ((G-N)/(2*N*(G+1)*i.alpha))
  }

  retval <- 0;
  c <- .C("gs_cross_eval_es_01_GV"       ,
          as.double (i)                  ,
          retval = as.integer(retval)    ,
          as.character(data.set)          )
  gs.stop.timer (info.level=1)
}


"gs.cross.info" <- function( bestn  = 0,
                             sortby = "index",
                             out.filename ="cross.info" ,
                             data.set     ="default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_cross_info_01_GV"       ,
           as.integer  (bestn)         ,
           as.character(sortby)        ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  x <- NULL
  if  (-2 != c$retval ) {
    x <-  read.table(file=x.filename,header=T,stringsAsFactors=F)
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.chrom.stats" <- function ( filename = "chrom.stats",
                               auxfiles = TRUE,
                               data.set ="default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,filename,".chr",sep="")
  retval <- 0;
  c <- .C("gs_chrom_stats_01_GV",
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=F)
      unlink(out.filename)
      names(x) <- c("Chrom","NLoci","Length")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.calc.rf" <- function ( map.function = "Haldane",
                           auxfiles     = TRUE     ,
                           data.set     = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,map.function,".rf",sep="")
  retval <- 0;
  c <- .C("gs_calc_rf_01_GV"         ,
          as.character(map.function) ,
          as.character(out.filename) ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)       )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=F)
      unlink(out.filename)
      names(x) <- c("Chrom","Locus1","Locus2","RecFreq")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.calc.q" <- function ( pop.type = "DH",
                          t        = 0   , 
                          auxfiles = TRUE,
                          data.set ="default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,pop.type,".q",sep="")
  retval <- 0;
  c <- .C("gs_calc_q_01_GV",
          as.character(pop.type),
          as.integer(t),
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=F)
      unlink(out.filename)
      names(x) <- c("Chrom","Locus1","Locus2","q")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.get.simpop" <- function ( pop.name ,
                              data.set ="default" )
  {
   bg <- plabsim$allele.missing.indicator
   retval <- 0
   c <- .C("st_get_simpop_GV",
            as.character        (pop.name),
            as.integer          (bg),
            retval = as.integer (retval)  ,
            as.character        (data.set)   )
}

"st.set.simpop" <- function ( pop.name ,
                              data.set ="default" )
{
  st.write.marker.data ( format    = "n"      ,
                         nfilename = pop.name ,
                         data.set  = data.set   )
  smp.file <- st.outdir(paste(pop.name,".nmp",sep=""))
  spo.file <- st.outdir(paste(pop.name,".npo",sep="")) 
  sgp.file <- st.outdir(paste(pop.name,".ngp",sep="")) 
  #
  reset.all()
  cl <- scan(sgp.file,quiet =T)
  set.genome.par( no.chrom=length(cl), no.hom=2, chrom.len=cl )
  linkage.map.load(smp.file)
  population.matrix.load(spo.file,backcross=TRUE)
  #
  pp <- as.character(list.populations()$PopName)
  population.rename ( pp[length(pp)], pop.name )
  for (i in (length(pp)-1):1 ) 
       population.concat ( pop.name,pp[i])
  #
  unlink ( smp.file )
  unlink ( spo.file )
  unlink ( sgp.file )
}


"il.eval.lines" <- function ( allele    = 1,
                              chrom     = 0,
                              lower     = 0,
                              upper     = 0,
                              exclude   = 0,
                              data.set  = "default"  )

{
 
  #
  mar.stats <- st.marker.data.statistics(data.set=data.set,ind=T,mar=F,gen=F)
  #
  NoLines <- 0
  c  <- .C( "il_get_population_size_GV",
             NoLines = as.integer (NoLines) , 
             as.character         (data.set)   )
  #
  retval <- 0
  abs  <- numeric(c$NoLines)
  tot  <- numeric(c$NoLines)
  nseg <- numeric(c$NoLines)
  lseg <- numeric(c$NoLines)

  if ( !( ( length(chrom) == length(lower) ) &&
          ( length(lower) == length(upper) )    )  )
      {
        st.info(-2,"Incorrect target region borders")
        return ( invisible(NULL) )
      }

  if ( ( length(chrom) < c$NoLines ) ||
       ( length(lower) < c$NoLines ) ||
       ( length(upper) < c$NoLines )    )
      {
        chrom <- rep (chrom[1],c$NoLines)
        lower <- rep (lower[1],c$NoLines)
        upper <- rep (upper[1],c$NoLines)
      }
  
  c <- .C ( "il_eval_lines_GV",
           as.integer (allele),
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer(exclude),
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  
  if ( -2 == c$retval ) { return ( invisible(NULL) )}

  r <- data.frame ( name=mar.stats$individual.list$Name, nseg=c$nseg, lseg=c$lseg,
                   abs=c$abs,tot=c$tot,rel=(c$abs/c$tot) )
    
    return ( r )

}


il.ideal.library <- function ( n.c =  5  ,            # number of chromosomes 
                               l.c = 100 ,            # length of chromosomes
                               s.l =  20 ,            # segment length   
                               data.set = "default" )
{
  m <- NULL
  for (c in 1:n.c) 
  m <- rbind ( m,data.frame( chrom = rep(c,l.c),
                             pos   = 1:l.c,
                             name  = paste ("C",c,"MA",1:l.c,sep=""),
                             class = rep(paste ("C",c,sep=""),l.c) ))
  m.name = (paste("tmp-map",round(runif(1,0,100000)),sep="") )
  p.name = (paste("tmp-pop",round(runif(1,0,100000)),sep="") )
  write.table (m,file=st.outdir(m.name),quote=F,row.names=F,col.names=F,sep =" ")
  x <- matrix("AA",nrow=n.c*l.c,ncol=n.c*l.c/s.l)
  colnames(x) <- paste ( "IL", 1:(n.c*l.c/s.l),  sep="" )
  rownames(x) <- m$name
  i <- 1
  for (p in 1:(n.c*l.c)) {
    x[p,i] <- "CC"
    if (p/s.l==round(p/s.l)) i <- i+1
  }
  write.table (x,file=st.outdir(p.name),quote=F, sep =" ")
  sve <- st.input.dir
  st.input.dir <<- st.output.dir
  st.read.marker.data (p.name, format="m",auxfiles = FALSE,data.set=data.set)
  st.read.map ( m.name,data.set=data.set)
  st.input.dir <<- sve
  unlink (st.outdir(m.name))
  unlink (st.outdir(p.name))
}

il.overlapping.library <- function ( n.c =  5  ,            # number of chromosomes 
                                     l.c = 100 ,            # length of chromosomes
                                     s.l =  20 ,            # segment length   
                                     data.set = "default" )
{
  m <- NULL
  for (c in 1:n.c ) 
  m <- rbind ( m,data.frame( chrom = rep(c,l.c),
                             pos   = 1:l.c,
                             name  = paste ("C",c,"MA",1:l.c,sep=""),
                             class = rep(paste ("C",c,sep=""),l.c) ))
  m.name = (paste("tmp-map",round(runif(1,0,100000)),sep="") )
  p.name = (paste("tmp-pop",round(runif(1,0,100000)),sep="") )
  write.table (m,file=st.outdir(m.name),quote=F,row.names=F,col.names=F,sep =" ")
  x <- matrix("AA",nrow=n.c*l.c,ncol=n.c*(1+l.c/s.l))
  colnames(x) <- paste ( "IL", 1:(n.c*(1+l.c/s.l)),  sep="" )
  rownames(x) <- m$name
  i <- 1
  for (p in 1:(n.c*l.c)) {
    x[p,i] <- "CC"
    if (p/s.l==round(p/s.l)) i <- i+1
    if (p/l.c==round(p/l.c)) i <- i+1
  }
  i <- 2
  for (p in 1:(n.c*l.c)) {
    x[p,i] <- "CC"
    if (p/s.l==round(p/s.l)) i <- i+1
    if (p/l.c==round(p/l.c)) i <- i+1
  }

  write.table (x,file=st.outdir(p.name),quote=F, sep =" ")
  sve <- st.input.dir
  st.input.dir <<- st.output.dir
  st.read.marker.data (p.name, format="m",auxfiles = FALSE,data.set=data.set)
  st.read.map ( m.name,data.set=data.set)
  st.input.dir <<- sve
  unlink (st.outdir(m.name))
  unlink (st.outdir(p.name))
}


"il.eval.library" <- function ( rp.allele = 1,
                                dp.allele = 2,
                                chrom     = 0,
                                lower     = 0,
                                upper     = 0,
                                data.set  = "default"  )
{
  retval <- 0
  cov    <- 0
  dep    <- 0
  seglib <- 0
  c <- .C ( "il_eval_library_GV",
            cov    = as.numeric (cov   ),
            dep    = as.numeric (dep   ),
            seglib = as.numeric (seglib),
                     as.integer (dp.allele),
            retval = as.integer (retval), 
            as.character        (data.set) )
  if ( -2 == c$retval ) { return ( invisible(NULL) )}
  r <- c ( cov=c$cov , dep=c$dep , seglib=c$seglib )
  #
  NoLines <- 0
  c  <- .C( "il_get_population_size_GV",
           NoLines = as.integer (NoLines) , 
           as.character         (data.set)   )
  abs  <- numeric(c$NoLines)
  tot  <- numeric(c$NoLines)
  nseg <- numeric(c$NoLines)
  lseg <- numeric(c$NoLines)
  retval <- 1
  #
  c0 <- .C ( "il_eval_lines_GV",
           as.integer (dp.allele),            # Donor allele
           as.integer ( rep(0,c$NoLines) ) ,  # genome wide
           as.numeric ( rep(0,c$NoLines) ) ,
           as.numeric ( rep(0,c$NoLines) ) ,
           as.integer ( 0 ) , 
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c0$retval ) { return ( invisible(NULL) )}
  
  r <- c(r,
         dpg  = mean (c0$abs/c0$tot),
         ndps = mean (c0$nseg)      ,
         ldps = mean (c0$lseg)          )
  if ( 0 == chrom[1] ) {  return ( r ) }

  # Analysis if target regions are defined
  if ( ( length(chrom) != c$NoLines ) ||
       ( length(lower) != c$NoLines ) ||
       ( length(upper) != c$NoLines )    )
    {
      st.info(-2,"Incorrect length of target region borders")
      return ( invisible(NULL) )
    }
  #
  c1 <- .C ( "il_eval_lines_GV",
           as.integer (dp.allele),    # Donor allele
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer ( 0 ),          # in the target regions
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c1$retval ) { return ( invisible(NULL) )}
  #
  c2 <- .C ( "il_eval_lines_GV",
           as.integer (rp.allele),    # Recipient allele
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer ( 1 ),          # outside the target regions
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c2$retval ) { return ( invisible(NULL) )}
  #
  c3 <- .C ( "il_eval_lines_GV",
           as.integer (dp.allele),    # Donor allele
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer ( 1 ),          # outside the target regions
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c3$retval ) { return ( invisible(NULL) )}
  #
  r <- c(r,
         tr.dpg  = mean (c1$abs/c1$tot),
         nt.rpg  = mean (c2$abs/c2$tot),
         nt.ndps = mean (c3$nseg)      ,
         nt.ldps = mean (c3$lseg)          )
  #
  return ( r )

}

  st.select.phen <- function (p, n=0, t=0,decreasing=TRUE,nsmall=1) {
    p  <- p[ order( p$y, decreasing=decreasing), ]
    if (t>0)  {
      if (decreasing)  p <- p[p$y>t,] else p <- p[p$y<t,]
    }
    if ((n>0)&&(n<nrow(p))) {
      p  <- p[1:n,]
    }
    p$descr <- paste(p$i,format(p$y,nsmall=nsmall))
    return(p)
  }

st.simple.ggt.plot <- function ( data.set="default", ... ) {
    f <- 1; i <- l <- 73;
    nme <-st.marker.data.statistics(data.set=data.set)$individual.list$Name
    max <- length(st.marker.data.statistics(data.set=data.set)$individual.list$Name)
    m <- 1
    repeat {
      f.name <- sprintf ("%s_part%03i",data.set,m)
      if (l>max) {fk <- (max-f+3) / (i+2) } else {fk<- 1}
      st.plot.ggt ( ...,
                    f.ind=f,
                    l.ind=l,
                    data.set=data.set ,
                    plt.ptsize = 8,
                    plt.height = 12*fk,
                    plt.width = 8,
                    d.t = 50,
                    d.map = 100,
                    plt.pdf=T,
                    plt.fname= f.name)
      f <- f+i
      l <- l+i
      m <- m+1  
      if (f>max) break;
    }
  }

st.mark.alleles <- function ( markers = NULL,      # Provide markers directly
                              data.set="default", 
                              effect.set=NULL,     # Significant in another data set
                              target.set,          
                              hap.lst = NULL,      # Second data set uses haplotypes  
                              alpha=0.05,           
                              p.adjust.method="none",
                              replace = "3"
                             )
{
  gs.start.timer()
 # Determine marker if not provided 
  if ( is.null(markers) ) 
    {
    if ( is.null(effect.set) ) effect.set <- data.set
    #
    tmp.filename = paste("tmp-", round(runif(1, 0, 1e+05)), sep = "")
    auxfiles = TRUE
    if ( "" == st.output.dir ) op <- "" else op <- "/"
    x.filename <- paste(st.output.dir, op, tmp.filename, sep = "")
    retval <- 0
    pvals.av <- gs.check.pvals(effect.set)
    if ( FALSE == pvals.av ){
      st.info(-2,"No p-values estimated")
      return(invisible(NULL))
    }
    c <- .C("gs_return_pvals_01_GV", as.character(x.filename), 
            retval = as.integer(retval), as.character(effect.set))
    x <- read.table(file = x.filename, header = TRUE, 
        check.names=FALSE, stringsAsFactors=FALSE)
    unlink(x.filename)
    x <- x[-1,]
    tol = 1e-8
    x <- subset (x, abs(effect) > tol) 
    x$pvalue <- p.adjust(x$pvalue,p.adjust.method)
    markers <- substring(rownames(subset(x,pvalue<alpha)),1,7)
    #
    mark.trans <- NULL
    if ( !is.null(hap.lst) ) {
      for (i in 1:length(markers)) {
        m = as.character(hap.lst$Markers[ markers[i] == hap.lst$Name ])
        mark.trans <- c(mark.trans,unlist(strsplit(m,";")))
      }
      markers <- mark.trans
    }
  }
# Replace marker alleles 
    n <- round(runif(1,0,100000))
    m.file <- paste("tmp-",n,sep="") 
    #
    st.write.marker.data ( format = "m", mfilename=m.file,data.set=data.set )
    fpop <-  paste(st.od(m.file),".mpo",sep="")
    la <- read.table (fpop, stringsAsFactors=F )
    #
    hom.rep = paste(replace,"/",replace,sep="")
    het.rep = paste("1","/",replace,sep="")
    #
    idx <- row.names(la) %in% markers
    for (i in 1:length(idx)) 
      if (TRUE == idx[i]){
        i2 <-  la[i,] == "2/2";        la[i,i2] <- hom.rep
        i2 <-  la[i,] == "1/2";        la[i,i2] <- het.rep
        }
    #     
    write.table(la,file=fpop,quote=F)
# Load data set
    ifs <- st.input.dir
    st.input.dir <<- st.output.dir
    print(paste(st.input.dir,m.file,".mpo",sep=""))
    st.read.marker.data ( paste(m.file,".mpo",sep=""),
                          format="m",data.set=target.set,
                          auxfiles=FALSE)
    st.read.map ( paste(m.file,".mmp",sep=""),
                  data.set=target.set,auxfiles=FALSE)
# Remove temp files
    unlink(paste(st.id(m.file),".mpo",sep=""))
    unlink(paste(st.id(m.file),".mmp",sep=""))
# Set back global input directory
    st.input.dir <<- ifs
}


gs.write.pseff <- function (e,beta0=NA,file="ps.effects") {
  if ( is.na(beta0) )
    write ( e[1,4],file ) else write ( beta0,file )
  ps.effects <- rbind ( cbind (e[2:nrow(e),c(1,2)],
                               effect=e[2:nrow(e),4] / 2 ) ,
                        cbind ( marker= e[2:nrow(e),1],
                                allele= e[2:nrow(e),3],
                                effect=-e[2:nrow(e),4]/2 ) )
  write.table (ps.effects,file,quote=F,
               row.names=FALSE, col.names=FALSE , append=TRUE )
}

"gs.neg.effall" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C( "gs_neg_effall_01_GV"      ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
}

"gs.pos.effall" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C( "gs_pos_effall_01_GV"      ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
}

########################################################################################
# Check regress for estimation of the shrinkage factor BEGIN
########################################################################################

                                        
#  Code from the CRAN package package by
#  Davil Clifford, Peter McCullagh, HJ Auinger
#  Version 1.3-10 of 2013-04-18
#

gs.esteff.blup <- function (data.set="default")
  {
    st.start.timer(depth=2)
    i <- st.get.info.level()
    st.set.info.level(0)
    gs.build.Z ( auxfiles=F, data.set=data.set )
    V <- gs.build.V(data.set=data.set)
    y <- st.return.performance.data(data.set=data.set)$y
    muell <- capture.output ( rg <- regress( y ~ 1 , ~ V) )
    sh <- rg$sigma[2] /rg$sigma[1]
    gs.lambda.const ( lambda=sh, auxfiles=F,data.set=data.set )    
    gs.mme.coeff  ( auxfiles=F, data.set=data.set)
    gs.mme.rhs    ( auxfiles=F, data.set=data.set) 
    gs.mme.solve  ( auxfiles=F, data.set=data.set)
    st.set.info.level(i)
    gs.stop.timer (depth=2,info.level=1)
  }


BLUP <- function(model,RE=NULL) {
    ## model - output from regress

    ## RE - vector of names of random effects, should match names of
    ## model$sigma

    ## Returns the conditional mean, variance and variance covariance
    ## matrices of random effects given the data. Default is all
    ## random effects ignoring those associated with the identity
    ## matrix.

    if(length(RE)==0) RE <- setdiff(model$Vnames,"In")
    if(any(is.na(match(RE,model$Vnames)))) stop(paste("RE should be a subset of",model$Vnames))
    if(class(model)!="regress") stop("model should be of class regress")

    ## conditional expected values - all of them
    Wy <- model$W %*% (model$model[[1]] - model$fitted)
    Us <- list()
    for(ii in 1:length(model$Z)) {
        Us[[ii]] <- as.vector(model$sigma[ii] * t(model$Z[[ii]]) %*% Wy)
        names(Us[[ii]]) <- colnames(model$Z[[ii]])
    }
    names(Us) <- names(model$sigma)

    ## Conditional Covariance matrices - ignore In for now
    ks <- unlist(lapply(model$Z,ncol))
    COV <- diag(rep(model$sigma,ks))
    TT <- NULL
    for(ii in (1:length(model$Z))) TT <- cbind(TT,model$sigma[ii] * model$Z[[ii]])
    COV <- COV - t(TT) %*% model$W %*% TT

    indices <- list()
    start <- c(1,1+cumsum(ks))[1:length(ks)]
    end <- cumsum(ks)
    for(ii in 1:length(start)) indices[[ii]] <- start[ii]:end[ii]

    names(start) <- names(model$sigma)
    names(end) <- names(model$sigma)

    means <- unlist(Us)
    rownames(COV) <- colnames(COV) <- names(means)
    vars <- diag(COV)

    allInd <- NULL
    for(ii in 1:length(RE)) {
        ind <- match(RE[ii],model$Vnames)
        ##output[[ii]] <- list("Mean"=means[start[ind]:end[ind]],
        ##                     "Variance"=vars[start[ind]:end[ind]],
        ##                     "Covariance"=COV[start[ind]:end[ind],start[ind]:end[ind]])
        allInd <- c(allInd,start[ind]:end[ind])
    }
    ##names(output) <- RE
    output <- list("Mean"=means[allInd],
                   "Variance"=vars[allInd],
                   "Covariance"=COV[allInd,allInd])
    return(output)
}

SWsolve <- function(S,K,D,Dinv=NULL,b) {
    ## solve(a,b) where a has the form SKS' + D using the Sherman Morrison Woodbury identities

    if(is.matrix(K) & is.matrix(D) & !is.null(Dinv)) {
        ## Case 1 - all are matrices
        tSDi <- crossprod(S,Dinv)
        Kinv <- solve(K)
        ret <- solve(Kinv + tSDi %*% S, tSDi)
        ret <- Dinv - crossprod(tSDi,ret)
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }

    if(is.numeric(K) & !is.null(Dinv)) {
        tSDi <- crossprod(S,Dinv)
        ret <- solve(1/K * diag(ncol(S)) + tSDi %*% S, tSDi)
        ret <- Dinv - crossprod(tSDi,ret)
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }

    if(is.numeric(D) & is.matrix(K)) {
        ret <- 1/D * diag(nrow(S)) - 1/D^2 * S %*% solve(solve(K) + 1/D * crossprod(S),t(S))
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }

    if(is.numeric(K) & is.numeric(D)) {
        ret <- 1/D * diag(nrow(S)) - 1/D^2 * S %*% solve(1/K * diag(ncol(S)) + 1/D * crossprod(S),t(S))
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }
}


SWsolve2 <- function(Zlist,clist,b) {
    ## Invert a matrix of the form sum [ clist[i] tcrossprod(Zlist[[ii]])) ] using Sherman Morrison Woodbury Identities
    if(length(Zlist)!=(length(clist)-1)) stop()
    k <- length(Zlist)
    D <- clist[1] * tcrossprod(Zlist[[1]])
    diag(D) <- diag(D) + clist[k+1]
    Dinv <- SWsolve(Zlist[[1]],clist[1],clist[k+1])
    if(k==1) {
        if(!missing(b)) Dinv <- Dinv %*% b
        return(Dinv)
    }
    for(ii in 2:k) {
        Dinv <- SWsolve(Zlist[[ii]],clist[ii],D,Dinv)
        D <- D + clist[ii]*tcrossprod(Zlist[[ii]])
    }
    if(!missing(b)) Dinv <- Dinv %*% b
    return(Dinv)
}

## generalised inverse

ginv <- function (X, tol = sqrt(.Machine$double.eps))
{
  ## taken from library MASS
  if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
    stop("X must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X))
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive))
    array(0, dim(X)[2:1])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
                            t(Xsvd$u[, Positive, drop = FALSE]))
}

summary.regress <- function(object, ...) object

print.regress <- function(x, digits=3, fixed.effects=T, ...)
  {
      cat("Likelihood kernel: K = ")
      if(length(x$kernel) == 1){
          cat(max(sign(x$kernel), 0))
      } else if(!is.null(x$Kcolnames)) cat(x$Kcolnames, sep="+")

      cat("\n\nMaximized log likelihood with kernel K is ",round(x$llik,digits),"\n",sep=" ")
      indent.lin <- max(nchar(dimnames(x$beta)[[1]]))
      indent.var <- max(nchar(x$Vnames))
      indent <- max(indent.lin,indent.var)

      extra.space <- ""
      space.var <- extra.space
      for(i in 0:(indent-indent.var)) space.var <- paste(space.var," ",sep="")
      space.lin <- extra.space
      for(i in 0:(indent-indent.lin)) space.lin <- paste(space.lin," ",sep="")

      coefficients <- cbind(x$beta,x$beta.se)
      dimnames(coefficients)[[2]] <- c("Estimate","Std. Error")
      coefficients <- round(coefficients,digits)
      if(fixed.effects) {
          cat("\nLinear Coefficients:\n")
          row.names(coefficients) <- paste(space.lin,dimnames(x$beta)[[1]],sep="")
          print(coefficients)
          cat("\n")
      } else {
          cat("\nLinear Coefficients: not shown\n\n")
      }

      ## New version of regress automatically converts to the linear
      ## scale - as if pos was a vector of zeroes

      var.coefficients <- cbind(x$sigma,sqrt(diag(as.matrix(x$sigma.cov))))
      row.names(var.coefficients) <- paste(space.var,x$Vnames,sep="")
      dimnames(var.coefficients)[[2]] <- c("Estimate","Std. Error")
      var.coefficients <- round(var.coefficients,digits)
      cat("Variance Coefficients:\n")
      print(var.coefficients)
      cat("\n")
  }

regress <- function(formula, Vformula, identity=TRUE, kernel=NULL,
                    start=NULL, taper=NULL, pos, verbose=0, gamVals=NULL,
                    maxcyc=50, tol=1e-4, data,
                    fraction=NULL,print.level=NULL){

  ## Vformula can just be something like ~ V0 + V1
  ## or leave it out or Vformula=NULL
  ## assume its in the form ~ V1 + V2 + ... + Vn or missing or Vformula=NULL
  ## for random effects and random interactions for factors A and B include
  ## ~ A + B + I(A:B)

  if(!is.null(print.level)) {
    cat("\nWarning: print.level has been replaced by verbose and has been deprecated.\nIt will be removed in the next version of regress\n\n")
    verbose <- print.level
  }
  if(!is.null(print.level)) {
      cat("\nWarning: fraction has been deprecated and replaced by taper - a vector of values from 0 to 1 giving a unique fraction for each step of the Newton Raphson algorithm.\n")
  }

  if(verbose>9) cat("Extracting objects from call\n")
  if(missing(data)) data <- environment(formula)
  mf <- model.frame(formula,data=data,na.action=na.pass)
  mf <- eval(mf,parent.frame())
  y <- model.response(mf)

  model <- list()
  ##model$formula <- formula
  ##model$Vformula <- Vformula
  model <- c(model,mf)

  if(missing(Vformula)) Vformula <- NULL

  ## Find missing values in fixed part of the model :: Change Aui Mar 1 2012
  isNA <-  apply(is.na(mf), 1, any)

  if(!is.null(Vformula))
  {
      V <- model.frame(Vformula,data=data,na.action=na.pass)
      V <- eval(V, parent.frame())
      # find missings in random part of the model Aui Mar 1 2012
      mfr <- is.na(V)
      if(ncol(mfr) == 1){
        isNA <- isNA | mfr
      } else {
        isNA <- isNA | apply(mfr[,!apply(mfr,2,all)], 1, any) # use only columns of the matrix with some nonmissing values
      }
      rm(mfr)
      Vcoef.names <- names(V)
      V <- as.list(V)
      k <- length(V)
  } else {
      V <- NULL
      k <- 0
      Vcoef.names=NULL
  }

  if(ncol(mf)==1) mf <- cbind(mf,1)
  X <- model.matrix(formula, mf[!isNA,]) # Aui Mar 1 2012 account for missings in the random part

  y <- y[!isNA]
  n <- length(y)
  Xcolnames <- dimnames(X)[[2]]
  if(is.null(Xcolnames)) {
      Xcolnames <- paste("X.column",c(1:dim(as.matrix(X))[2]),sep="")
  }

  X <- matrix(X, n, length(X)/n)
  qr <- qr(X)
  rankQ <- n-qr$rank
  if(qr$rank) {
      X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
      Xcolnames <- Xcolnames[qr$pivot[1:qr$rank]]
  } else {
      cat("\nERROR: X has rank 0\n\n")
  }

  if(verbose>9) cat("Setting up kernel\n")

  if(missing(kernel)){
      K <- X
      colnames(K) <- Xcolnames
      reml <- TRUE
      kernel <- NULL
  } else {
      if(length(kernel)==1 && kernel>0){
          K <- matrix(rep(1, n), n, 1)
          colnames(K) <- c("1")
      }
      if(length(kernel)==1 && kernel<=0) {
          K <- Kcolnames <- NULL
          KX <- X
          rankQK <- n
      }
      if(length(kernel) > 1) {
          ##K is a matrix I hope :: Change Aui Mar 1 2012
          if(is.matrix(kernel)) {
              K <- kernel[!isNA,]
          } else {
              K <- model.frame(kernel, data=data, na.action=na.pass)
              K <- eval(K, parent.frame())
              if(ncol(K) == 1){
                dimNamesK <- dimnames(K)
                K <- K[!isNA, ]
                dimNamesK[[1]] <- dimNamesK[[1]][!isNA]
                K <- data.frame(V1 = K)
                dimnames(K) <- dimNamesK
              } else {
                K <- K[!isNA, ]
              }
              K <- model.matrix(kernel, K)
          }
      }
      reml <- FALSE
  }

  if(!is.null(K)){
      Kcolnames <- colnames(K)
      qr <- qr(K)
      rankQK <- n - qr$rank
      if(qr$rank == 0) K <- NULL else {
          K <- matrix(K[, qr$pivot[1:qr$rank]],n,qr$rank)
          Kcolnames <- Kcolnames[qr$pivot[1:qr$rank]]
          KX <- cbind(K, X) # Spanning K + X: Oct 12 2011
          qr <- qr(KX)
          KX <- matrix(KX[, qr$pivot[1:qr$rank]],n,qr$rank) # basis of K+X
      }
  }

  if(missing(maxcyc)) maxcyc <- 50
  if(missing(tol)) tol <- 1e-4
  delta <- 1

  if(verbose>9) cat("Removing parts of random effects corresponding to missing values\n")
  ## remove missing values
  for(i in 1:k)
  {
      if(is.matrix(V[[i]]))
      {
          V[[i]] <- V[[i]][!isNA, !isNA]
      }
      if(is.factor(V[[i]]))
      {
          V[[i]] <- V[[i]][!isNA]
      }
  }

  In <- diag(rep(1,n),n,n)

  if(identity) {
      ##if(k) for(i in k:1) V[[i+1]] <- V[[i]]
      ##V[[1]] <- In
      V[[k+1]] <- as.factor(1:n)
      names(V)[k+1] <- "In"
      k <- k+1

      ##Vcoef.names <- c("Id",Vcoef.names)
      Vcoef.names <- c(Vcoef.names,"In")
      Vformula <- as.character(Vformula)
      Vformula[1] <- "~"
      Vformula[2] <- paste(Vformula[2],"+In")
      Vformula <- as.formula(Vformula)
  }

  model <- c(model,V)
  model$formula <- formula
  model$Vformula <- Vformula

  ## specify which parameters are positive and which are negative
  ## pos = c(1,1,0) means first two parameters are positive, third is either
  if(!missing(pos)) pos <- as.logical(pos)
  if(missing(pos)) pos <- rep(FALSE,k)
  pos <- c(pos,rep(FALSE,k))
  pos <- pos[1:k]

  ## Sherman Morrison Woodbury identities for matrix inverses can be brought to bear here
  if(verbose>9) cat("Checking if we can apply the Sherman Morrison Woodbury identites for matrix inversion\n")
  if (all(sapply(V, is.factor)) & k>2 ) {  # Contribution by Hans Jurgen Auinger
      SWsolveINDICATOR <- TRUE
  } else SWsolveINDICATOR <- FALSE
  Z <- list()
  for (i in 1:length(V)) {
      if (is.factor(V[[i]])) {
          Vi <- model.matrix(~V[[i]] - 1)
          colnames(Vi) <- levels(V[[i]])
          Z[[i]] <- Vi
          V[[i]] <- tcrossprod(Vi)
      } else{
          Z[[i]] <- V[[i]]
      }
  }
  names(Z) <- names(V)

  ## So V is always a list of variance coavriance matrices, Z contains
  ## the model matrices of factors when we need to invoke the Sherman
  ## Woodbury identities

  ## Expected Fisher Information
  A <- matrix(rep(0, k^2), k, k)
  entries <- expand.grid(1:k,1:k)

  x <- rep(0,k)
  sigma <- c(1,rep(0, k-1))
  stats <- rep(0, 0)

  ## START ALGORITHM

  if(missing(taper)){
      taper <- rep(0.9, maxcyc)
      if(missing(start) && k>1) taper[1:2] <- c(0.5, 0.7)
  } else {
      taper <- pmin(abs(taper), 1)
      if((l <- length(taper)) < maxcyc) taper <- c(taper, rep(taper[l], maxcyc-l))
  }

  if(!is.null(start)) {
      ## pad start with zeros if required
      start <- c(start, rep(1,k))
      start <- start[1:k]
  }

  if(k>2 && is.null(start)) start <- rep(var(y,na.rm=TRUE),k)
  if(k==1 && is.null(start)) start <- var(y,na.rm=TRUE)

  if(is.null(start) && k==2) {
      if(missing(gamVals)) {
          gamVals <- seq(0.01,0.02,length=3)^2
          gamVals <- sort(c(gamVals,seq(0.1,0.9,length=3),1-gamVals))
          gamVals <- 0.5
      }
      if(length(gamVals)>1) {
          if(verbose>=1) cat("Evaluating the llik at gamma = \n")
          if(verbose>=1) cat(gamVals)
          if(verbose>=1) cat("\n")
          reg.obj <- reml(gamVals,y,X,V[[1]],V[[2]],verbose=verbose)
          llik <- reg.obj$llik
          llik <- as.double(llik)
          if(verbose>=2) cat(llik,"\n")
          gam <- gamVals[llik==max(llik)]
          gam <- gam[1]
          if(verbose>=2) cat("MLE is near",gam,"and llik =",max(llik),"there\n")
      }
      if(length(gamVals)==1) {
          ## go straight to the Newton Raphson at gamVals
          gam <- gamVals[1]
          reg.obj <- list(rms=var(y))
      }
      start <- c(1-gam,gam)*reg.obj$rms
      ## it tends to take huge steps when starting at gam=0.9999
      if(gam==0.9999) {
          taper[1] <- taper[1]/100
          maxcyc <- maxcyc*10
      }
      if(verbose>=1) cat(c("start algorithm at",round(start,4),"\n"))
  }

  if(is.null(start) & k>2) {
      ## Never gets here by default - but this could be implemented,
      ## though it does add on a few extra iterations at the
      ## start.... not necessary in basic examples

      LLvals <- NULL
      ## equal weights
      V2 <- V[[2]]
      for(ii in 3:k) V2 <- V2 + V[[ii]]
      LLvals <- c(LLvals,reml(0.5,y,X,V[[1]],V2)$llik)
      ## Most at one end
      V2 <- V[[1]] + V2 ## total of all Vs
      for(ii in 1:k) {
          V2 <- V2 - V[[ii]]
          LLvals <- c(LLvals,reml(0.75,y,X,V2,V[[ii]])$llik)
      }
      best <- which.max(LLvals)
      if(verbose) {
          cat("Checking starting points\n")
          cat("llik values of", LLvals, "\n")
      }
      if(best==1) {
          start <- rep(var(y,na.rm=TRUE),k)
      } else {
          start <- rep(0.25,k)
          start[best] <- 0.75
      }
  }

  sigma <- coef <- start
  ## reparameterise so everything will get into the correct spot after exp
  coef[pos] <- log(sigma[pos])
  coef[!pos] <- sigma[!pos]

  ## Set the memory requirements beforehand
  T <- vector("list", length=k)
  for(ii in 1:k) T[[ii]] <- matrix(NA,n,n)

  for(cycle in 1:maxcyc){

      ## Limit how far out we go on the logarithmic scale
      ind <- which(pos)
      if(length(ind)) {
          coef[ind] <- pmin(coef[ind],20)
          coef[ind] <- pmax(coef[ind],-20) ## so on regular scale everything is between exp(-20) and exp(20)
          sigma[ind] <- exp(coef[ind])
      }

      if(verbose>=1) {
          cat(cycle, "sigma =",sigma)
          ##cat(sigma)
      }

      ##Sigma <- matrix(0,dim(V[[1]])[1],dim(V[[1]])[2])
      ## if(verbose>9) cat("Sherman Morrison Woodbury",SWsolveINDICATOR,"\n")

      if(!SWsolveINDICATOR) {
          Sigma <- 0
          ## can we get rid of this loop?
          for(i in 1:k) Sigma <- Sigma + V[[i]]*sigma[i]

          W <- solve(Sigma,In)
      } else {
          W <- SWsolve2(Z[1:(k-1)],sigma)
      }

      if(is.null(K)) WQK <- W else {
          WK <- W %*% K
          WQK <- W - WK %*% solve(t(K)%*%WK, t(WK))
      }
      if(reml) WQX <- WQK else {
          WX <- W %*% KX        # including the kernel (Oct 12 2011)
          WQX <- W - WX %*% solve(t(KX)%*%WX, t(WX))
      }

      rss <- as.numeric(t(y) %*% WQX %*% y)

      ##if(verbose>9) cat("Sigma[1:5]",Sigma[1:5],"\n")
      ##if(verbose>9) cat("RSS",rss,"WQX[1:5]",WQX[1:5],"\n")
      sigma <- sigma * rss/rankQK
      coef[!pos] <- sigma[!pos]
      coef[pos] <- log(sigma[pos])
      WQK <- WQK * rankQK/rss
      WQX <- WQX * rankQK/rss
      rss <- rankQK ## looks bad but the rss is absorbed into WQK so the rss term comes out of eig below

      eig <- sort(eigen(WQK,symmetric=TRUE,only.values=TRUE)$values, decreasing=TRUE)[1:rankQK]
      if(any(eig < 0)){
          cat("error: Sigma is not pos def on contrasts: range(eig)=", range(eig), "\n")
          WQK <- WQK + (tol - min(eig))*diag(nobs)
          eig <- eig + tol - min(eig)
      }
      ldet <- sum(log(eig))
      llik <- ldet/2 - rss/2
      if(cycle == 1) llik0 <- llik
      delta.llik <- llik - llik0
      llik0 <- llik

      if(verbose && reml) cat(" resid llik =", llik,"\n")
      if(verbose && !reml) cat(" llik =", llik, "\n")

      if(verbose) cat(cycle, "adjusted sigma =",sigma)
      if(cycle>1) {
          if(verbose && reml) cat(" delta.llik =", delta.llik, "\n")
          if(verbose && !reml) cat(" delta.llik =", delta.llik, "\n")
      } else cat("\n")

      ## now the fun starts, derivative and expected fisher info
      ## the 0.5 multiple is ignored, it is in both and they cancel

      ##T <- list(NULL)
      x <- NULL

      ## derivatives are now D[[i]] = var.components[i]*V[[i]]
      var.components <- rep(1,k)
      ind <- which(pos)
      if(length(ind)) var.components[ind] <- sigma[ind]

      ## Slow part - order k n-squared
      if(!SWsolveINDICATOR) {
          if(identity) {
              T[[k]] <- WQK
              if(k>1) {
                  for(ii in (k-1):1) T[[ii]] <- WQK %*% V[[ii]]
              }
          } else {
              for(ii in 1:k) T[[ii]] <- WQK %*% V[[ii]]
          }
      } else {
          if(identity) {
              T[[k]] <- WQK
              if(k>1) {
                  for(ii in (k-1):1) T[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
              }
          } else {
              for(ii in 1:k) T[[ii]] <- tcrossprod(WQK %*% Z[[ii]],Z[[ii]])
          }

          ##if(k>=6) {
          ##    ## One line to do all - may be memory inefficient though
          ##    T <- lapply(Z,function(x) tcrossprod(WQ %*% x, x))
          ##} else {
          ##    if(k>=1) T[[1]] <- tcrossprod(WQ %*% Z[[1]], Z[[1]])
          ##    if(k>=2) T[[2]] <- tcrossprod(WQ %*% Z[[2]], Z[[2]])
          ##    if(k>=3) T[[3]] <- tcrossprod(WQ %*% Z[[3]], Z[[3]])
          ##    if(k>=4) T[[4]] <- tcrossprod(WQ %*% Z[[4]], Z[[4]])
          ##    if(k>=5) T[[5]] <- tcrossprod(WQ %*% Z[[5]], Z[[5]])
          ##}
      }


      x <- sapply(T,function(x) as.numeric(t(y) %*% x %*% WQX %*% y - sum(diag(x))))
      x <- x * var.components

      ## See nested for loops commented out below - evaluating the Expected Fisher Information, A
      ff <- function(x) sum(T[[x[1]]] * t(T[[x[2]]])) * var.components[x[1]] * var.components[x[2]]
      aa <- apply(entries,1,ff)
      A[as.matrix(entries)] <- aa

      ##for(i in 1:k)
      ##{
      ##if(identity && i==1 && pos[i]) {
      ##  T[[i]] <- WQ
      ##} else

      ##if(SWsolveINDICATOR==FALSE) {
      ##    T[[i]] <- WQ %*% V[[i]]
      ##} else {
      ##    T[[i]] <- WQ %*% tcrossprod(V[[i]])
      ##}
      ##x[i] <- as.numeric(t(y) %*% T[[i]] %*% WQ %*% y - sum(diag(T[[i]])))
      ##x[i] <- x[i]*var.components[i]

      ##for(j in c(1:i))
      ##  {
      ## expected fisher information
      ## ##A[j,i] <- Tr(T[[j]] %*% T[[i]])
      ## ##the Ts are not symmetric, hence the transpose below

      ##A[j,i] <- sum(T[[j]] * t(T[[i]]))
      ##A[j,i] <- A[j,i]*var.components[i]*var.components[j]
      ##A[i,j] <- A[j,i]
      ##}
      ##}

      stats <- c(stats, llik, sigma[1:k], x[1:k])
      if(verbose==-1) {
          ##cat(c(rllik1, rllik2, sigma[1:k], x[1:k]),"\n")
      }

      A.svd <- ginv(A)
      x <- A.svd %*% x

      if(qr(A)$rank < k){
          if(cycle==1) {
              if(verbose) {
                  cat("Warning: Non identifiable dispersion model\n")
                  ##print(round(A,6))
                  cat(sigma)
                  cat("\n")
              }
          }
      }

      ## end of newton-raphson step
      ## x is  -l'(sigma)/l''(sigma)
      ## hence we add instead of subtract
      ## taper controls the proportion of each step to take
      ## for some reason 0.5 works very well

      ##if(all(pos==1)) x <- sign(x) * pmin(abs(x),5) ## limit maximum shift we can take in one step to 5 units on log scale
      coef <- coef + taper[cycle] * x
      sigma[!pos] <- coef[!pos]
      sigma[pos] <- exp(coef[pos])

      ## check the change in llik is small
      if(cycle > 1 & abs(delta.llik) < tol*10) break
      if(max(abs(x)) < tol) break
  }

  ## Recompute Sigma at adjusted sigman values
  if(!SWsolveINDICATOR) {
      Sigma <- 0
      ## can we get rid of this loop?
      for(i in 1:k) Sigma <- Sigma + V[[i]]*sigma[i]

      W <- solve(Sigma,In)
  } else {
      W <- SWsolve2(Z[1:(k-1)],sigma)
  }


  if(cycle==maxcyc)
  {
      ## issue a warning
      if(verbose) cat("WARNING:  maximum number of cycles reached before convergence\n")
  }

  stats <- as.numeric(stats)
  stats <- matrix(stats, cycle, 2*k+1, byrow=TRUE)
  colnames(stats) <- c("llik",paste("s^2_", Vcoef.names, sep=""), paste("der_", Vcoef.names, sep=""))

  WX <- W %*% X
  XtWX <- crossprod(X,WX)
  cov <- XtWX
  cov <- solve(cov, cbind(t(WX),diag(1,dim(XtWX)[1])))
  beta.cov <- matrix(cov[,(dim(t(WX))[2]+1):dim(cov)[2]],dim(X)[2],dim(X)[2])
  cov <- matrix(cov[,1:dim(t(WX))[2]],dim(X)[2],dim(X)[1])

  beta <- cov %*% y
  beta <- matrix(beta,length(beta),1)
  row.names(beta) <- Xcolnames
  beta.se <- sqrt(abs(diag(beta.cov)))
  pos.cov <- (diag(beta.cov)<0)
  beta.se[pos.cov] <- NA
  beta.se <- matrix(beta.se,length(beta.se),1)
  row.names(beta.se) <- Xcolnames
  rms <- rss/rankQ

  fitted.values <- X%*%beta
  Q <- In -  X %*% cov

  predicted <- NULL
  predictedVariance <- NULL
  predictedVariance2 <- NULL
  if(identity) {
      gam <- sigma[k]  ## coefficient of identity, last variance term
      if(SWsolveINDICATOR) {
          ## Sigma is undefined
          Sigma <- 0
          for(i in 1:k)
          {
              Sigma <- Sigma + V[[i]] * sigma[i]
          }
      }
      predicted <- fitted.values + (Sigma - gam*In) %*% W%*%(y - fitted.values)
      ## predictedVariance <- Sigma - (Sigma - gam*In) %*% W %*% (Sigma - gam*In)
      ## really the last term should be transposed but they are square symmetric matrices here so....
      ## If you multiply it out and take the diagonal only it simplifies to....
      predictedVariance <- 2*gam - gam^2*diag(W)
      ## Variance of a new observation conditional on data (assuming
      ## known beta)
      predictedVariance2 <- diag(gam^2 * WX %*% beta.cov %*% t(WX))
      ## Additional variance of a new observation conditional on data taking into
      ## account variation in beta
  }

  ## scale dictated by pos
  sigma.cov <- (A.svd[1:k, 1:k] * 2)

  ## Last Step:  June 17th 2005
  ## Convert the estimates for the variance parameters, their standard
  ## errors etc to the usual scale

  FI <- A/2

  ## convert FI using pos
  FI.c <- matrix(0,dim(FI)[1],dim(FI)[2])
  FI.c <- FI / tcrossprod((sigma-1)*pos+1)
  ##for(i in 1:dim(FI)[1])
  ##  for(j in 1:dim(FI)[2])
  ##    FI.c[i,j] <- FI[i,j]/(((sigma[i]-1)*pos[i]+1)*((sigma[j]-1)*pos[j]+1))

  sigma.cov <- ginv(FI.c)
  names(sigma) <- Vcoef.names
  rownames(sigma.cov) <- colnames(sigma.cov) <- Vcoef.names

  result <- list(trace=stats, llik=llik, cycle=cycle, rdf=rankQ,
                 beta=beta, beta.cov=beta.cov, beta.se=beta.se,
                 sigma=sigma[1:k], sigma.cov=sigma.cov[1:k,1:k], W=W,
                 Q=Q, fitted=fitted.values, predicted=predicted,
                 predictedVariance=predictedVariance,
                 predictedVariance2=predictedVariance2,pos=pos,
                 Vnames=Vcoef.names, formula=formula,
                 Vformula=Vformula, Kcolnames=Kcolnames,
                 model=model,Z=Z, X=X, Sigma=Sigma)
  class(result) <- "regress"
  return(result)
}

## when two matrices are passed to regress this is also called
## to evaluate the REML at certain values of gamma and find a
## good place to start the regress algorithm

reml <- function(lambda, y, X, V0, V1,verbose=0){

  if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)

      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
  V0 <- V0[isNA==F,isNA==F]
  V1 <- V1[isNA==F,isNA==F]
  X <- X[isNA==F,]
  X <- as.matrix(X)

  qr <- qr(X)
  ##print(qr$rank)
  n <- dim(as.matrix(y))[1]
  In <- diag(1,n)

  X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  llik <- rep(0, length(lambda))
  if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]

  n <- dim(X)[1]
  if(missing(V0)) V0 <- diag(rep(1, n), n, n)
  rank <- n - qr$rank
  ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")
  for(i in 1:length(lambda))
    {
      if(verbose>=2) cat(lambda[i],"\n")

      Sigma <- (1-lambda[i])*V0 + lambda[i] * V1
      ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
      ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
      ##W <- chol2inv(cholesky)
      ##WX <- W %*% X
      WX <- solve(Sigma,cbind(X,In))
      W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
      WX <- WX[,1:dim(X)[2]]
      XtWX <- t(X)%*%WX
      WQ <- W - WX%*%solve(XtWX,t(WX))
      rss <- t(y) %*% WQ %*% y
      logdetrss <- sum(log(eigen(rss)$values[1:q]))
      eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
      ldet <- sum(log(eVals))
      llik[i] <- Re(ldet*q/2 - rank*logdetrss/2)
    }
  imax <- sort.list(-llik)[1]
  lambdamax <- lambda[imax]
  curv <- 0
  if(imax > 1 && imax < length(lambda)){
    delta <- (lambda[imax+1] - lambda[imax-1])/2
    slope <-  (llik[imax+1] - llik[imax-1])/2
    curv <- llik[imax-1] -2*llik[imax] + llik[imax+1]
    lambdamax <- lambdamax - slope/curv * delta
    curv <- -curv/delta^2
  }
  lamMax <- lambdamax
  Sigma <- (1-lamMax)*V0 + lamMax * V1
  ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
  ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
  ##W <- chol2inv(cholesky)
  ##WX <- W %*% X
  WX <- solve(Sigma,cbind(X,In))
  W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
  WX <- WX[,1:dim(X)[2]]
  XtWX <- t(X)%*%WX
  FItWX <- solve(XtWX,t(WX))
  WQ <- W - WX%*%FItWX
  rss <- t(y) %*% WQ %*% y
  beta <- FItWX %*% y

  list(llik=as.numeric(llik),rms=rss/rank, beta=beta, gamma=lambda, gamMax=lambdamax,W=W)
}

remlOptimize <- function(y, X, V0, V1,verbose=0,...){

  if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)

      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
  V0 <- V0[isNA==F,isNA==F]
  V1 <- V1[isNA==F,isNA==F]
  X <- X[isNA==F,]
  X <- as.matrix(X)

  qr <- qr(X)
  ##print(qr$rank)
  n <- dim(as.matrix(y))[1]
  In <- diag(1,n)

  X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]

  n <- dim(X)[1]
  if(missing(V0)) V0 <- diag(rep(1, n), n, n)
  rank <- n - qr$rank
  ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")

  f <- function(lambda,verbose=verbose) {
    if(verbose>=2) cat(lambda,"\n")
    Sigma <- (1-lambda)*V0 + lambda * V1
    ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
    ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
    ##W <- chol2inv(cholesky)
    ##WX <- W %*% X
    WX <- solve(Sigma,cbind(X,In))
    W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
    WX <- WX[,1:dim(X)[2]]
    XtWX <- t(X)%*%WX
    WQ <- W - WX%*%solve(XtWX,t(WX))
    rss <- t(y) %*% WQ %*% y
    logdetrss <- sum(log(eigen(rss)$values[1:q]))
    eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
    ldet <- sum(log(eVals))
    llik <- Re(ldet*q/2 - rank*logdetrss/2)
    llik
  }

  res <- optimize(f,interval=c(0,1),maximum=TRUE,verbose=verbose,...)
  lamMax <- res$maximum
  llikMax <- res$objective

  Sigma <- (1-lamMax)*V0 + lamMax * V1
  ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
  ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
  ##W <- chol2inv(cholesky)
  ##WX <- W %*% X
  WX <- solve(Sigma,cbind(X,In))
  W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
  WX <- WX[,1:dim(X)[2]]
  XtWX <- t(X)%*%WX
  FItWX <- solve(XtWX,t(WX))
  WQ <- W - WX%*%FItWX
  rss <- t(y) %*% WQ %*% y
  beta <- FItWX %*% y

  list(llik=as.numeric(llikMax),rms=rss/rank, beta=beta, gamMax=lamMax,W=W)
}

########################################################################################
# Check regress for estimation of the shrinkage factor END
########################################################################################
