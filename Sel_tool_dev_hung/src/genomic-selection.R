
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
                              auxfiles     = FALSE,
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
                                  auxfiles     = FALSE       )
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
                              auxfiles     = FALSE      ,
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
                                 auxfiles = FALSE              )
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
                               auxfiles     = FALSE           ,
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
                              hsq         = 0.9        ,
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
           as.double(hsq)             ,
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
                              hsq         = 0.9        ,
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
           as.double(hsq)             ,
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
                              hsq          = 0.9      ,
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
           as.double(hsq)             ,
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
                                 auxfiles = FALSE          ,
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
                              auxfiles = FALSE          ,
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
                                  auxfiles = FALSE          ,
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
                                auxfiles = FALSE,
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
                              hsq         = 0.9         ,
                              out.filename = "eff.rmlc" ,
                              auxfiles = FALSE          ,
                              data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rmlc_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
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

"gs.testeff.rmlc" <- function( nperm       = 1000         ,
                               maxiter     = 1000         ,
                               precision   = 0.001        ,
                               hsq         = 0.9          ,
                               out.filename = "test.rmlc" ,
                               auxfiles = TRUE            ,
                               data.set ="default"         )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_testeff_rmlc_01_GV",
           as.integer(nperm),
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
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
                               hsq         = 0.9         ,
                               out.filename = "eff.rmlv" ,
                               auxfiles = FALSE          ,
                               data.set ="default"       )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rmlv_01_GV",
           as.integer(maxiter),
           as.double(precision),
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

"gs.testeff.rmlv" <- function( nperm       = 1000        ,
                               maxiter     = 1000        ,
                               precision   = 0.001       ,
                               hsq         = 0.9         ,
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
           as.double(hsq),
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


"gs.esteff.rmlr" <- function( maxiter     = 1000        ,
                              precision   = 0.001       ,
                              hsq         = 0.9         ,       
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
           as.double(hsq),
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
                              hsq         = 0.9         ,
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
           as.double(hsq),
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


"gs.cross.validation.0" <- function ( estimation.method          ,
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

"st.def.hblocks" <- function ( hap          =  5,
                               hap.unit     = 99,
                               ld.threshold = 0.8,
                               ld.criterion = "none",
                               tolerance    = 0,
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
          as.double(ld.threshold),
          as.character(ld.criterion),
          as.integer(tolerance),
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

"st.recode.hil" <- function ( out.filename = "blocks"  ,
                              auxfiles     = TRUE,
                              data.set   = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,out.filename,".hall",sep="")
  retval <- 0;
  c <- .C("gs_recode_hil_01_GV",
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
      names(x) <- c("Block","AlleleNr","AlleleDef")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
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

gs.esteff.blup <- function (data.set="default")
  {
    if ( !("package:regress" %in% search() ) )
      {
        st.info(-2,"Package regress is required")
        return()
      }
    st.start.timer(depth=2)
    i <- st.get.info.level()
    st.set.info.level(0)
    gs.build.Z ( auxfiles=F, data.set=data.set )
    V <- gs.build.V(data.set=data.set)
    y <- st.return.performance.data(data.set=data.set)$y
    muell <- capture.output ( rg <- regress( y ~ 1 , ~ V) )
    sh <- rg$sigma[2] /rg$sigma[1]
    gs.lambda.const ( lambda=sh, auxfiles=F,data.set=data.set )    
    gs.mmet.coeff  ( auxfiles=F, data.set=data.set)
    gs.mmet.rhs    ( auxfiles=F, data.set=data.set) 
    gs.mmet.solve  ( auxfiles=F, data.set=data.set)
    st.set.info.level(i)
    gs.stop.timer (depth=2,info.level=1)
  }

################################################################################
# BEGIN: Ridge regression after transformation
################################################################################

"gs.mmet.coeff" <- function(
                              out.filename = "mme.coeff" ,
                              auxfiles     = FALSE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  A.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mmet_coeff_01_GV",
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

"gs.mmet.coeff.3" <- function(
                              out.filename = "mme.coeff" ,
                              auxfiles     = FALSE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  A.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mmet_coeff_03_GV",
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


"gs.mmet.rhs" <- function ( out.filename = "mme.rhs" ,
                              auxfiles     = FALSE      ,
                              data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  b.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mmet_rhs_01_GV",
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

"gs.mmet.solve" <- function( out.filename = "mme.solution" ,
                               auxfiles     = FALSE        ,
                               data.set     = "default"      )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C("gs_mmet_solve_01_GV",
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

"gs.lambda.rmlct" <- function( maxiter     = 1000       ,
                              precision    = 0.001      ,
                              hsq          = 0.9       ,
                              out.filename = "lambda"  ,
                              auxfiles = FALSE          ,
                              data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_rmlct_01_GV"    ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.double(hsq)             ,
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

"gs.lambda.rmlat" <- function( alpha        = 1        ,
                              maxiter      = 1000     ,
                              precision    = 0.001    ,
                              hsq          = 0.9       ,
                              out.filename = "lambda" ,
                              auxfiles     = TRUE     ,
                              data.set     = "default" )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_lambda_rmlat_01_GV"     ,
           as.double(alpha)            ,
           as.integer(maxiter)         ,
           as.double(precision)        ,
           as.double(hsq)              ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)       )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 


"gs.esteff.rr" <- function( method       = "rmlct",  
                            maxiter      = 1000        ,
                            precision    = 0.001       ,
                            hsq          = 0.80        ,
                            alpha        = 1           ,
                            out.filename = "eff.rr"    ,
                            auxfiles     = FALSE       ,
                            data.set     = "default"     )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_rr_01_GV"        ,
           as.character(method)        ,
           as.integer(maxiter)         ,
           as.double(precision)        ,
           as.double(hsq)              ,
           as.double(alpha)            ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=T) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
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
  c <- .C( "gs_cross_validation_02_GV" ,
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

################################################################################
# END: Ridge regression after transformation
################################################################################

##############################################################################
# BEGIN: LD based haplotypes                                           2015-04
##############################################################################
                                        
"st.calc.ld" <- function ( ld.measure   = "r2" ,
                           auxfiles     = TRUE ,
                           data.set     = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,"ld",sep="")
  retval <- 0;
  c <- .C("gs_calc_ld_01_GV"         ,
          as.character(ld.measure) ,
          as.character(out.filename) ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)       )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=F)
      unlink(out.filename)
      names(x) <- c("Chrom","Locus1","Locus2","Name1","Name2","LD")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

# For problems with large LD files
"st.calc.ld.2" <- function ( ld.measure   = "r2" ,
                             auxfiles     = TRUE ,
                             keep.ldfile  = FALSE,
                             data.set     = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,"ld",sep="")
  retval <- 0;
  c <- .C("gs_calc_ld_01_GV"         ,
          as.character(ld.measure) ,
          as.character(out.filename) ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)       )
  x <- NULL
  if (TRUE==keep.ldfile)  {return(invisible(NULL))}
  if  ( (auxfiles) && (-2 != c$retval ) )
      {
          y <- scan(file=out.filename,
                    what=list(integer(),integer(),integer(),character(),character(),double()),
                    quiet = TRUE)
          x <- data.frame(y)
          unlink(out.filename)
          names(x) <- c("Chrom","Locus1","Locus2","Name1","Name2","LD")
      }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}


st.conv.genable <- function(data.set)
{
  ge <- st.marker.data.statistics(data.set=data.set)$genotypes 
  mp <- st.get.map(data.set=data.set)
  ilf <- data.frame (name=mp$Name,
                     chr=mp$Chrom,
                     pos=sprintf("%15.0f",round(mp$Pos*10000)),
                     strand="+",
                     ge[,2:ncol(ge)])
  for (i in 5:ncol(ilf)) ilf[,i] <- gsub( "/","" ,               
                                    gsub( "1","A",               
                                    gsub( "2","C",               
                                    gsub( "3","G",               
                                    gsub( "4","T",               
                                    gsub("-1","-", ilf[,i])))))) 
  write.table(ilf,"tmp.illu",quote=F,row.names=F,col.names=T)
  write.table(data.frame(id=names(ilf)[5:ncol(ilf)],trt=0,sex=0),
              file="tmp.phe",quote=F,row.names=F,col.names=T)
  convert.snp.illumina(inf="tmp.illu", out="tmp.raw", strand="file")
  gAdta <- load.gwaa.data ( phe="tmp.phe", gen="tmp.raw", force=TRUE)
  unlink("tmp.phe")
  unlink("tmp.illu")
  unlink("tmp.raw")
  return(gAdta)
}

st.LDplot.ld <- function (ld,chrom)
  {
    l <- ld[ld$Chrom==chrom,]
    n <- nrow(l)
    d <- (1+sqrt(1+8*n))/2
    LD <- matrix(nrow=d,ncol=d)
    for (i in 1:n) LD[ l[i,2],l[i,3] ] <- l[i,6]
    return(LD)
  }

st.LDplot.map <- function(chrom, data.set)
  {
    mp <- st.get.map(data.set=data.set)
    m <- mp[mp$Chrom==chrom,]$Pos
    return(m)
  }

"st.set.hblocks" <- function ( haplotype.list ,
                               hap.symbol   = "s",
                               out.filename = "blocks"  ,
                               auxfiles     = TRUE,
                               data.set     = "default" )
{
  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  out.filename <- paste(st.output.dir,op,out.filename,".hmap",sep="")
  retval <- 0;
  #
    {
      s <- strsplit(haplotype.list$Markers,";")
      for (i in 1:length(s)) s[[i]] <-  rep(i,length(s[[i]]))
      hl <- unlist(s)
    } 
  #
  c <- .C("gs_set_hblocks_01_GV",
          as.integer(length(hl)),
          as.integer(hl),
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


##############################################################################
# END:   LD based haplotypes                                           2015-04
##############################################################################

##############################################################################
# Begin:   Genetid distances for factoials                             2015-11
##############################################################################


st.genetic.distances.fct <- function ( measure = "mrd", # "mrd" "rd" "euc"
                                       format  = "l",   # "l" "m"
                                       split   = 1,     # Last individual of the first set
                                       filename  = "genetic.distances",
                                       auxfiles  = FALSE,
                                       data.set  = "default" )
{
   if ( "" == st.output.dir ) op <- "" else op <- "/"
   out.filename = paste(st.output.dir,op,filename,".gdi",sep="")
#   
   gs.cross.eval.gd.fct ( dist=measure , split=split, data.set=data.set )
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

"gs.cross.eval.gd.fct" <- function ( dist     = "rd",
                                     split    = 1   ,
                                     data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_gd_fct_01_GV"    ,
          as.character(dist)          ,
          as.integer(split)           ,
          retval = as.integer(retval) ,
          as.character(data.set)       )
  gs.stop.timer (info.level=1)
}

##############################################################################
# End:  Genetic distances for factoials                                2015-11
##############################################################################

##############################################################################
# BAYES                                                              2017-07
##############################################################################

"gs.esteff.bayesB" <- function( numit ,
                                out.filename = "eff.bayesB" ,
                                auxfiles = FALSE,
                                data.set ="default"              )
{
  gs.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  x.filename <- paste(st.output.dir,op,out.filename,sep="")
  retval <- 0;
  c <- .C( "gs_esteff_bayesB_01_GV",
           as.double(numit),
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

##############################################################################
# ENDE: BAYES                                                         2017-07
##############################################################################
