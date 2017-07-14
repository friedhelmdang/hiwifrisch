/////////////////////////////////////////////////////////////////////////
// SeletionTools main source file genomic-selection.c                   //
// (c) Matthias Frisch                                                   //
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#define MAXTHREADS 256

#ifdef MT
   #include <omp.h> 
   #include <pthread.h>
#else
  #define omp_get_thread_num() 0
  #define omp_set_num_threads(thr) {}
  #define omp_get_num_threads() 0
  #define omp_set_dynamic(thr) {}
#endif

#ifdef use_openblas
  #include "cblas.h"
  #include "lapacke.h"
  #define use_blaslapack 1
#endif

#ifdef use_mkl
  #include "mkl.h"
  #define use_blaslapack 1
#endif

#ifdef MultPrec
  #include <gmp.h>
#endif

/* Maximum string size for general use */
#define gs_MAXSTR 4096

/* Maximum string size for names of markers and indviduals */
#define gs_MAXNME 128

/* Maximum number of alleles per marker */
#define gs_MAXALL 4096

/* Maximum line size in input files */
#define gs_MAXLNE 65536


/* Common errors */
#define ERR_F {gs_info(gv,-2,"File open error");*retval=-2;return;}
#define ERR_M {gs_info(gv,-2,"Memory allocation error");*retval=-2;return;}

/* Global variables and definitions  */
typedef unsigned long int gs_UNSIGNED ;
# define gs_UNSIGNED_MAX ULONG_MAX

// typedef __float128 gs_FLOAT ;
// typedef __float80 gs_FLOAT ;
typedef double gs_FLOAT ;

/* Type linkage map*/
typedef struct
{
    int    chrom;
    double pos;
    char   name [gs_MAXSTR];
    char   class[gs_MAXSTR];
} gs_mappoint_type;

/* Pair of two unsigned values */
typedef struct
{
  gs_UNSIGNED a;
  gs_UNSIGNED n;
} gs_usp;

/* Parameters of a cross of two individuals*/
typedef struct
{
  gs_UNSIGNED p1;
  gs_UNSIGNED p2;
  gs_FLOAT    gd;  // Genetic distance of the parents
  gs_FLOAT    mu;  // Expected value of the cross
  gs_FLOAT    mi;  // Minimum value of the cross
  gs_FLOAT    ma;  // Maximum value of the cross
  gs_FLOAT    va;  // Variance of the cross
  gs_FLOAT    es;  // Expectation of the selected fraction
} gs_cross_type;
 
/* Chromosomes and recombination frequencies */
typedef struct
{
  gs_UNSIGNED        NLocC;       // Number of loci on the chromosome
  gs_mappoint_type** first_loc;   // Pointer to the first locus 
  gs_UNSIGNED        first_loc_i; // Index of the first locus 
  gs_FLOAT*          rf;          // Matrix of recombination frequencies
  gs_FLOAT*          q;           // Matrix of conditional probabilites
} gs_chrom_type;


/* Data set */

typedef struct {
  
  char name[gs_MAXSTR];
  char msg [gs_MAXSTR];
  int info_level;

  /* Variables initialized by gs_read_markerdata */

  int MaxIL;        /* Maximum length of an individual name */
  int MaxML;        /* Maximum length of an marker name */

  gs_UNSIGNED  
    NoInd ,   /* Number of individuals                   */
    NoMar ;   /* Number of markers                       */

  char (*IndName) [gs_MAXNME]; /* Individual names*/ 
  char (*MarName) [gs_MAXNME]; /* Marker names*/ 

  int       * All1  ; /* Genotypic data */
  int       * All2  ; /* Genotypic data */
  gs_FLOAT  * genv  ; /* Genotypic value */

  /* Variables initiaized by gs_marker_stats*/

  gs_UNSIGNED (*NoAll)  ;             /* Number of allele per Marker */
  int (*All) [gs_MAXALL] ;          /* Allele names */
  gs_UNSIGNED (*AllCnt) [gs_MAXALL] ; /* Allele count */

  double (*Pic) ; /* Pic value */
  double (*Mis) ; /* Frequency of Missing values per marker */

  double (*IndMis) ; /* Frequency of missing markers per ind */

  /* Variables for gs_build_Z and gs_build_X */

  double* Z ;   /* Design matrix */
  double* y ;   /* Phenotype     */

  gs_UNSIGNED NoEff ;          /* Number of effects  */
  char (*EffNme) [gs_MAXNME] ; /* Marker names */ 
  int *EffAll ;                /* Allele number of the effect */
  int *CmpAll ;                /* Allele number complementary allele */

  double aa; double aA; double AA; /* Coding alleles */
  double an; double nn; double nA; /* Coding alleles */

  /* Variables for gs_build_A */
  /* to solve Ax=b */

  gs_UNSIGNED dimu ;
  gs_UNSIGNED dimA ;
  double*  A ;   
  double*  u ;   
  double*  b ;
  double*  As;      // Save the value of A and b when using 
  double*  bs ;     // in place matrix operations

  /* Variables for the transformed model */

  double*  GI ;  // G-star inverted: inv(ZZ')  
  double*  a ;   // Solution vector for the transformed model   

  /* p-value for test of u */
  double*  p;

  /* Shrinkage factor */
  double* ladi ;

  /* Variables for gs_validation */
  double* est ;

  /* Variables for gs_single_marker_reg_01 and _aov_01*/
  double* SMsqm ;
  double* SMsqe ;
  double* SMF   ;
  double* SMp   ;
  double* SMvar ;

  /* Variables for variance component estimation */
  double* sigma;
  double  sigma_e;

  /* Linkage map */  
  gs_UNSIGNED         NLoc;   // Number of loci in the map
  gs_mappoint_type  * map_S;  // Array for storage
  gs_mappoint_type ** map;    // Pointers for access
  gs_UNSIGNED       * mdl;    // Index of the data line for a mappoint

  /* Storage of recoded of haploblock alleles */
  int * Ref1  ; 
  int * Ref2  ; 

  /* Linkage map of haplotypes */
  gs_UNSIGNED         NHap  ;  // Number of haplotype blocks
  gs_mappoint_type  * Hap_S ;  // Array for storage
  gs_UNSIGNED  * h_u;          // Lower border of block 
  gs_UNSIGNED  * h_o;          // Upper border
  gs_UNSIGNED  * h_n;          // No. of loci in the block

  /* Parameters of a cross */
  gs_UNSIGNED      NCrs;   // Number of crosses
  gs_cross_type  * crs_S;  // Array for storage
  gs_cross_type  **crs;    // Array of pointers for access

  /* Chromosomes and recombination frequencies */
  gs_UNSIGNED NChrom;      // Number of chromosomes
  gs_chrom_type * chrom ;  // Array of chromosome parameters

  /* Temprary storage for iterations without free */
  double *t0,*t1,*t2,*t3,*t4,*t5,*t6,*t7,*t8,*t9;
  gs_UNSIGNED *u0;   
} gs_varset_type;


/* Standard variable set is accessed with pointer */
gs_varset_type GV =  {
  .name       = "default",    
  .msg        = "",
  .info_level = 0 ,
  .MaxIL      = 0, 
  .MaxML      = 0, 
  .NoInd      = 0,  
  .NoMar      = 0,  
  .IndName    = NULL,
  .MarName    = NULL,
  .All1       = NULL,
  .All2       = NULL,
  .genv       = NULL,
  .NoAll      = NULL,
  .All        = NULL,
  .AllCnt     = NULL,
  .Pic        = NULL,
  .Mis        = NULL,
  .IndMis     = NULL,
  .Z          = NULL, 
  .y          = NULL, 
  .NoEff      = 0,   
  .EffNme     = NULL,
  .EffAll     = NULL,
  .CmpAll     = NULL,
  .aa         = 0, 
  .aA         = 1,
  .AA         = 2,
  .an         = 0.5,
  .nn         = 1  ,
  .nA         = 1.5,
  .dimu       = 0,
  .dimA       = 0,
  .A          = NULL,   
  .u          = NULL,   
  .b          = NULL,
  .As         = NULL,
  .bs         = NULL,
  .GI         = NULL,
  .a          = NULL,
  .p          = NULL,
  .ladi       = NULL,
  .est        = NULL,
  .SMsqm      = NULL,
  .SMsqe      = NULL,
  .SMF        = NULL,
  .SMp        = NULL,
  .SMvar      = NULL,
  .sigma      = NULL,
  .sigma_e    = 0,
  .NLoc       = 0,
  .map_S      = NULL,
  .map        = NULL,
  .mdl        = NULL,
  .Ref1       = NULL,
  .Ref2       = NULL,
  .NHap       = 0,
  .Hap_S      = NULL,
  .h_u        = NULL,
  .h_o        = NULL,
  .h_n        = NULL,
  .NCrs       = 0,
  .crs_S      = NULL,
  .crs        = NULL,
  .NChrom     = 0,
  .chrom      = NULL,
  .t0=NULL, .t1=NULL, .t2=NULL, .t3=NULL, .t4=NULL, 
  .t5=NULL, .t6=NULL, .t7=NULL, .t8=NULL, .t9=NULL,
  .u0=NULL
} ;

gs_varset_type *gv = &GV; 
gs_varset_type *sh = &GV; 

gs_varset_type ES; gs_varset_type *es = &ES; int esi = 0;
gs_varset_type VS; gs_varset_type *vs = &VS; int vsi = 0;

gs_varset_type **gs_data = NULL;
gs_UNSIGNED      gs_ndta = 0;

/* Is the random number generator started? */
int rgs = 0;

/* Global constants */
const char * gs_tok_del = "; \r\n\t\""  ;
char * nix  = ""; char ** nixp  = &nix;
int    nein = 0 ; int  *  neinp = & nein;
double null = 0 ; double *  nullp = & null;

int Srv=0; int* Sretval = &Srv;
char*empty=""; char ** Sout_filename = &empty;
int auf = 0; int*Sauxfiles = &auf;
int mz = -2; int*imz = &mz;

///////////////////////////////////////////////////////////////////////////
// Distribution functions, Numerical Recipies (Start)                    //
///////////////////////////////////////////////////////////////////////////

void gs_info ( gs_varset_type *, int , char* );

void nrerror(char error_text[])
{
  gs_info(gv,-2,error_text);
  exit(1);
}

#define MAXIT 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double betacf(double a, double b, double x)
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) 
	  nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN

double gammln_(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double betai(double a, double b, double x)
{
	double bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln_(a+b)-gammln_(a)-gammln_(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double st_fprob (double f,  gs_UNSIGNED df1, gs_UNSIGNED df2)
/* Prob ( F(df1,df2) > f ) */
{
	double prob = betai(0.5*df2,0.5*df1,df2/(df2+df1*f));
	if (prob > 1.0) prob=2.0-prob;
	return (prob);
}

void st_fprob_GV ( double* f, int* df1, int* df2, double* prob)
{
  *prob = st_fprob( *f,	(gs_UNSIGNED)(*df1), (gs_UNSIGNED)(*df2) ); 
}

double st_tprob (double t,  gs_UNSIGNED df)
/* Prob ( |t(df)| > t ) */
{
  double prob = 0.5*betai(0.5*df,0.5,df/(df+t*t)) ;
  return (prob);
}

void st_tprob_GV ( double* t, int* df, double* prob)
{
  *prob = st_tprob( *t,	(gs_UNSIGNED)(*df) ); 
}

///////////////////////////////////////////////////////////////////////////
// Distribution functions, Numerical Recipies (End)                      //
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Build in linear Algebra Routines, Numerical Recipies (Start)          //
///////////////////////////////////////////////////////////////////////////

#define NRANSI
#define TINY 1.0e-20
#define NR_END 1
#define FREE_ARG char*

double *vector(gs_UNSIGNED nl, gs_UNSIGNED nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}
void free_vector(double *v, gs_UNSIGNED nl)
/* free a double vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

double **matrix(gs_UNSIGNED nrl, gs_UNSIGNED nrh, 
		gs_UNSIGNED ncl, gs_UNSIGNED nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  gs_UNSIGNED i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

void free_matrix(double **m, gs_UNSIGNED nrl, gs_UNSIGNED ncl)
/* free a double matrix allocated by dmatrix() */
{
	free( (FREE_ARG) (m[nrl]+ncl-NR_END) );
	free( (FREE_ARG) (m+nrl-NR_END)      );
}

gs_UNSIGNED *ivector(gs_UNSIGNED nl, gs_UNSIGNED nh)
/* allocate an gs_UNSIGNED vector with subscript range v[nl..nh] */
{
  gs_UNSIGNED *v;

  v=(gs_UNSIGNED*)malloc( (size_t) ( (nh-nl+1+NR_END) *
			              sizeof(gs_UNSIGNED)) );
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}

void free_ivector(gs_UNSIGNED *v, gs_UNSIGNED nl)
/* free an gs_UNSIGNED vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
double **convert_matrix(double *a, gs_UNSIGNED nrl, gs_UNSIGNED nrh, 
			gs_UNSIGNED ncl, gs_UNSIGNED nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	gs_UNSIGNED i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

void free_convert_matrix(double **b, gs_UNSIGNED nrl)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void ludcmp_st(double **a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *d)
{
  gs_UNSIGNED i,imax=0,j,k;
  double big,dum,sum,temp;
  double *vv;

  vv=vector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,1);
}

void gs_ludcmp_st(double *a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *d)
{
  double **aa = convert_matrix( a, 1,n, 1,n);

  ludcmp_st ( aa,  n, indx-1, d );
 
  free_convert_matrix (aa,1);
 
}

void lubksb_st(double **a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *b)
{
	gs_UNSIGNED i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void gs_lubksb_st(double *a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *b )
{
  double **aa = convert_matrix( a, 1,n, 1,n);

  lubksb_st ( aa, n, indx-1, b-1 );   

  free_convert_matrix (aa,1);

}


void ludcmp_mt(double **a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *d)
{
  gs_UNSIGNED i,imax=0,j,k;
  double big,dum,sum,temp;
  double *vv;

  vv=vector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0)  nrerror("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    #pragma omp parallel for private(i,k,sum,dum)  
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        #pragma omp critical
	{
	  big=dum;
	  imax=i;
	}
      }
    }
    if (j != imax) {
      #pragma omp parallel for private(k,dum) 
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      #pragma omp parallel for private(i)
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,1);
}


void lubksb_mt(double **a, gs_UNSIGNED n, gs_UNSIGNED *indx, double b[])
{
	gs_UNSIGNED i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii) {
		  for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		}
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


void gs_luinv_st (double *a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *inv)
{

  gs_UNSIGNED i,j;

  double *col = malloc ( n * sizeof(double) );
  
  for ( j=0 ; j < n ; j++) {
     for ( i=0 ; i < n ; i++ ) col[i] = 0.0 ;
     col[j] = 1.0 ;
     gs_lubksb_st ( a, n, indx, col );
     for ( i=0 ; i < n ; i++ ) inv[ i*n +j ] = col[i] ;
  }

  free(col);

}

void choldc_st ( double **a, gs_UNSIGNED  n, double p[] )
{
  
  gs_UNSIGNED i,j,k;
  double sum;
  
  for (i=1;i<=n;i++) {
    for (j=i;j<=n;j++) {
      for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
      if (i == j) {
	if (sum <= 0.0)
	  nrerror("choldc failed");
	p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }
}

void gs_choldc_st (double *a, double *b, gs_UNSIGNED n )
{
  gs_UNSIGNED i,j;

  double **aa = convert_matrix( a, 1,n, 1,n);
  double **bb = convert_matrix( b, 1,n, 1,n);
 
  choldc_st ( aa, n, b );
 
  // Write diagonal elements from b to a
  for ( i=0; i<n; i++ ) aa[i][i] = b[i]; 

  // Fill upper triangle of a with 0's
  for ( i=0; i<(n-1); i++ ) for (j=i+1; j<n; j++)  aa[i][j] = 0 ;   

  // Write  the transpose to b
  for ( i=0; i<n; i++ ) for (j=0; j<n; j++)  bb[i][j] = aa[j][i];

  free_convert_matrix (aa,1);
  free_convert_matrix (bb,1);
}

#undef TINY
#undef NRANSI
#undef NREND
#undef FREE_ARG


void (*ludcmp)() = ludcmp_st;
void (*lubksb)() = lubksb_st;

void gs_solve_Axb_stmt ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
			 double *A, double  *x, double  *b , 
			 int*   retval    )
{
    gs_UNSIGNED i;

    for ( i=0; i<nra ; i++) x [i] = b[i];

    double **a, *xx, d;
    gs_UNSIGNED *indx;

    indx=ivector(1,nra);
    a = convert_matrix( &A[0], 1,nra, 1,nca);
    xx = &x[0] - 1;

    ludcmp(a,nra,indx,&d);
    lubksb(a,nra,indx,xx);

    free_convert_matrix(a,1);
    free_ivector(indx,1);

    *retval = 0;
}

void gs_solve_Axb_st ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		       double *A, double  *x, double  *b , 
		       int* retval )
{
  ludcmp = ludcmp_st;
  gs_solve_Axb_stmt (nra,nca,A,x,b,retval);
}

void gs_solve_Axb_mt ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		       double *A, double  *x, double  *b, 
		       int* retval )
{
  ludcmp = ludcmp_st;
  gs_solve_Axb_stmt (nra,nca,A,x,b,retval);
}

///////////////////////////////////////////////////////////////////////////
// Linear Algebra Routines, Numerical Recipies (End)                     //
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Linear algebra without libraries                                      //
///////////////////////////////////////////////////////////////////////////

/* matrix matrix add */
void gs_mm_add ( gs_UNSIGNED nr, gs_UNSIGNED nc,
		 double *A, double* B, double *C )
{
  gs_UNSIGNED i, j;

  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
	C [ i*nc + j ] = A [ i*nc + j ] + B [ i*nc + j ] ;
}

/* scalar matrix add */
void gs_sm_add ( gs_UNSIGNED nr, gs_UNSIGNED nc,
		 double *A, double* B, double *C )
{
  gs_UNSIGNED i, j;

  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
	C [ i*nc + j ] = *A  + B [ i*nc + j ] ;
}

/* matrix matrix multiply */
void gs_mm_mul_st ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		    gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		    double *A, double  *B, double  *C )	
{
  gs_UNSIGNED i, j, k;

  if ( nca != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
  
  for (i = 0; i < nra; i++) 
    for (j = 0; j < ncb; j++)
      {
	C [ i*ncb + j ] = 0 ;
	for (k = 0; k < nca; k++) 
	  {
	    C [ i*ncb + j ] += A [ i*nca + k ] * B [ k*ncb + j ];
	  }
      }
}

void gs_mtm_mul_st ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
	             gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		     double *A, double  *B, double  *C )	
{
  gs_UNSIGNED i, j, k;

  if ( nca != ncb )
    {
      sprintf( gv->msg, 
	       "Incorrect dimensions for matrix multiplication: %lu %lu",
	       nca, ncb);
      gs_info(gv,-2,gv->msg);
      return;
    }
  
  for (i = 0; i < nra; i++) 
    for (j = 0; j < nrb; j++)
      {
	C [ i*nrb + j ] = 0 ;
	for (k = 0; k < nca; k++) 
	  {
	    C [ i*nrb + j ] += A [ i*nca + k ] * B [ j*ncb + k ];
	  }
      }
}

void gs_mm_mul_mt ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		    gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		    double *A, double  *B, double  *C )	
{

  if ( nca != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
 
  #pragma omp parallel
  {
    gs_UNSIGNED i, j, k;
    long double add;
 
    #pragma omp for
    for (i = 0; i < nra; i++) 
      for (j = 0; j < ncb; j++)
	{
	  add = 0 ;
	  for (k = 0; k < nca; k++) 
	    {
	      add += A [ i*nca + k ] * B [ k*ncb + j ];
	    }
	  C [ i*ncb + j ] = add;
	}
  }
}

/* matrix matrix multiply */
void gs_mm_mul_gF ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		    gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		    gs_FLOAT *A,  gs_FLOAT  *B,  gs_FLOAT  *C )	
{
  gs_UNSIGNED i, j, k;

  if ( nca != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
  
  for (i = 0; i < nra; i++) 
    for (j = 0; j < ncb; j++)
      {
	C [ i*ncb + j ] = 0 ;
	for (k = 0; k < nca; k++) 
	  {
	    C [ i*ncb + j ] += A [ i*nca + k ] * B [ k*ncb + j ];
	  }
      }
}

#ifdef MultPrec

/* matrix matrix multiply */
void gs_mm_mul_mpf ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		     mpf_t *A,  mpf_t  *B,  mpf_t  *C )	
{
  gs_UNSIGNED i, j, k;
  mpf_t t; mpf_init(t);

  if ( nca != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
  
  for (i = 0; i < nra; i++) 
    for (j = 0; j < ncb; j++)
      {
	mpf_set_d ( C[ i*ncb + j ] , 0 ) ;
	for (k = 0; k < nca; k++) 
	  {
	    mpf_mul ( t , A [ i*nca + k ] , B [ k*ncb + j ] );
	    mpf_add ( C[ i*ncb + j ] , C [ i*ncb + j ] ,t);
	  }
      }

  mpf_clear(t);
}

void gs_tmm_mul_mpf ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		      gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		      mpf_t *A,  mpf_t  *B,  mpf_t  *C )	
{
  gs_UNSIGNED i, j, k;
  mpf_t t; mpf_init(t);

  if ( nra != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
  
  for (i = 0; i < nca; i++) 
    for (j = 0; j < ncb; j++)
      {
	mpf_set_d ( C [ i*ncb + j ] , 0) ;
	for (k = 0; k < nra; k++) 
	  {
	    mpf_mul ( t , A [ k*nca + i ] , B [ k*ncb + j ] );
	    mpf_add ( C[ i*ncb + j ] , C [ i*ncb + j ] ,t);
	  }
       }

  mpf_clear(t);
}

#endif



/* transpose(matrix) matrix multiply */
void gs_tmm_mul_02 ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		     double *A, double  *B, double  *C )	
{
  gs_UNSIGNED i, j, k;

  if ( nra != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
  
  for (i = 0; i < nca; i++) 
    for (j = 0; j < ncb; j++)
      {
	C [ i*ncb + j ] = 0 ;
	for (k = 0; k < nra; k++) 
	  {
	    C [ i*ncb + j ] += A [ k*nca + i ] * B [ k*ncb + j ];
	  }
      }
}

/* transpose(matrix) matrix multiply */
void gs_tmm_mul_st ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		     double *A, double  *B, double  *C )	
{
  gs_UNSIGNED i, j, k;
  long double add;

  if ( nra != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
  
  for (i = 0; i < nca; i++) 
    for (j = 0; j < ncb; j++)
      {
	add =0;
	for (k = 0; k < nra; k++) 
	  {
	    add += A [ k*nca + i ] * B [ k*ncb + j ];
	  }
	C [ i*ncb + j ] = add ;

       }
}

/* transpose(matrix) matrix multiply */
void gs_tmm_mul_mt ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		     double *A, double  *B, double  *C )	
{

  if ( nra != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }

  #pragma omp parallel
  {
    gs_UNSIGNED i, j, k;
    long double add;
    
    #pragma omp for
    for (i = 0; i < nca; i++) 
      for (j = 0; j < ncb; j++)
	{
	  add =0;
	  for (k = 0; k < nra; k++) 
	    {
	      add += A [ k*nca + i ] * B [ k*ncb + j ];
	    }
	  C [ i*ncb + j ] = add ;
	}
  }
}

void gs_tmm_mul_gF ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		     gs_FLOAT *A, gs_FLOAT  *B, gs_FLOAT  *C )	
{
  gs_UNSIGNED i, j, k;
  gs_FLOAT add;

  if ( nra != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
  
  for (i = 0; i < nca; i++) 
    for (j = 0; j < ncb; j++)
      {
	add =0;
	for (k = 0; k < nra; k++) 
	  {
	    add += A [ k*nca + i ] * B [ k*ncb + j ];
	  }
	C [ i*ncb + j ] = add ;

       }
}

///////////////////////////////////////////////////////////////////////////
// // Linear algebra without libraries                                      //
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Linear Algebra
///////////////////////////////////////////////////////////////////////////


#ifdef use_gsl

 #include "gsl/gsl_blas.h"
 #include "gsl/gsl_linalg.h"

void gs_tmm_mul_gsl ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		      gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		      double *A, double  *B, double  *C )	
{
  if ( nra != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
 
   gsl_matrix_view AA = gsl_matrix_view_array(A, nra, nca);
   gsl_matrix_view BB = gsl_matrix_view_array(B, nrb, ncb);
   gsl_matrix_view CC = gsl_matrix_view_array(C, nca, ncb);

   gsl_blas_dgemm (CblasTrans, CblasNoTrans,
		   1.0, &AA.matrix, &BB.matrix,
  		   0.0, &CC.matrix);

}

void gs_mm_mul_gsl ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
		     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		     double *A, double  *B, double  *C )	
{
  if ( nca != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }
 
   gsl_matrix_view AA = gsl_matrix_view_array(A, nra, nca);
   gsl_matrix_view BB = gsl_matrix_view_array(B, nrb, ncb);
   gsl_matrix_view CC = gsl_matrix_view_array(C, nra, ncb);

   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		   1.0, &AA.matrix, &BB.matrix,
  		   0.0, &CC.matrix);

}

void gs_mtm_mul_gsl ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
	              gs_UNSIGNED nrb, gs_UNSIGNED ncb,
		      double *A, double  *B, double  *C )	
{
  if ( nca != ncb )
    {
      sprintf( gv->msg, 
	       "Incorrect dimensions for matrix multiplication: %lu %lu",
	       nca, ncb);
      gs_info(gv,-2,gv->msg);
      return;
    }
 
   gsl_matrix_view AA = gsl_matrix_view_array(A, nra, nca);
   gsl_matrix_view BB = gsl_matrix_view_array(B, nrb, ncb);
   gsl_matrix_view CC = gsl_matrix_view_array(C, nra, nrb);

   gsl_blas_dgemm (CblasNoTrans, CblasTrans,
		   1.0, &AA.matrix, &BB.matrix,
  		   0.0, &CC.matrix);

}


void gs_solve_Axb_gsl ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
			double *A, double  *x, double  *b, 
			int* retval  )
{
    int s;

    gsl_matrix_view AA = gsl_matrix_view_array (A, nra, nca);
    gsl_vector_view xx = gsl_vector_view_array (x, nra);
    gsl_vector_view bb = gsl_vector_view_array (b, nra);
     
    gsl_permutation * p = gsl_permutation_alloc (nra);
     
    gsl_linalg_LU_decomp (&AA.matrix, p, &s);
    gsl_linalg_LU_solve  (&AA.matrix, p, &bb.vector, &xx.vector);

    gsl_permutation_free(p);

    *retval = 0;
}


void gs_ludcmp_gsl(double *a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *d)
{
  int s;

  gsl_matrix_view AA = gsl_matrix_view_array (a, n, n);
  gsl_permutation * p = gsl_permutation_alloc (n);
     
  gsl_linalg_LU_decomp (&AA.matrix, p, &s);
 
  for (gs_UNSIGNED i = 0; i < n; i ++) 
    indx[i] = 1+(gs_UNSIGNED)gsl_permutation_get(p,(size_t)i) ;

  gsl_permutation_free(p);

  *d = (int)s;
}


void gs_lubksb_gsl(double *a, gs_UNSIGNED n, gs_UNSIGNED *indx, double *b )
{
    double *x = malloc ( n * sizeof(double) );
    for (gs_UNSIGNED i=0 ; i < n ; i++ ) x[i] = b[i] ;

    gsl_permutation * p = gsl_permutation_alloc (n);
    for (size_t i = 0; i < n; i ++) 
      p->data[i] = (size_t)indx[i] - 1;

    gsl_matrix_view AA = gsl_matrix_view_array (a, n, n);
    gsl_vector_view xx = gsl_vector_view_array (b, n);    // is returned!
    gsl_vector_view bb = gsl_vector_view_array (x, n);
          
    gsl_linalg_LU_solve  (&AA.matrix, p, &bb.vector, &xx.vector);

    gsl_permutation_free(p);
}

void gs_luinv_gsl ( double *a, gs_UNSIGNED n, gs_UNSIGNED *indx,double *i )
{

  gsl_permutation * p = gsl_permutation_alloc (n);
  for (size_t j = 0; j < n; j ++) 
    p->data[j] = (size_t)indx[j] - 1;

  gsl_matrix_view A = gsl_matrix_view_array (a, n, n);
  gsl_matrix_view I = gsl_matrix_view_array (i, n, n);
     
  gsl_linalg_LU_invert ( &A.matrix, p, &I.matrix);

  gsl_permutation_free(p);

}



void gs_svdcmp_gsl(double *a, double *v, double*s, gs_UNSIGNED m, gs_UNSIGNED n)

{                 // A = USV' ,  m >= n
                  // A: m x n     V: n x m   s: m 

  gsl_matrix_view AA = gsl_matrix_view_array (a, m, n);
  gsl_matrix_view VV = gsl_matrix_view_array (v, m, n);
  gsl_vector_view ss = gsl_vector_view_array (s, m);

  gsl_vector *w = gsl_vector_alloc (m);
       
  gsl_linalg_SV_decomp (&AA.matrix, &VV.matrix , &ss.vector, w);

  gsl_vector_free(w);
 
}

void gs_svdsol_gsl(double *u, double *v, double*s, gs_UNSIGNED m, gs_UNSIGNED n,
		   double *b, double *x)
{ 

  gsl_matrix_view UU = gsl_matrix_view_array (u, m, n);
  gsl_matrix_view VV = gsl_matrix_view_array (v, m, n);
  gsl_vector_view ss = gsl_vector_view_array (s, m);
  gsl_vector_view bb = gsl_vector_view_array (b, m);
  gsl_vector_view xx = gsl_vector_view_array (x, m);
       
  gsl_linalg_SV_solve ( &UU.matrix, &VV.matrix,  &ss.vector,  
			&bb.vector, &xx.vector);
}

#endif


#ifdef use_blaslapack

void gs_tmm_mul_blaslapack ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
			     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
			     double *A, double  *B, double  *C )	
{
  if ( nra != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }


  cblas_dgemm (CblasRowMajor, // op(A) x op(B)
	       CblasTrans,    // op(A): m x k 
	       CblasNoTrans,  // op(B): k x n,  C: m x n
	       (int)nca,  // m: Number of rows in op(A) (and C)
	       (int)ncb,  // n: Number of columns in op(B) (and C)
	       (int)nra,  // k: Number of columns in op(A), rows in op(B)
	       1.0, 
	       A,
	       (int)nca,  // lda: Row stride of A
	       B, 
	       (int)ncb,  // ldb: Row stride of B
	       0.0, 
	       C, 
	       (int)ncb); // ldc: Row stride of C
}

void gs_mtm_mul_blaslapack ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
			     gs_UNSIGNED nrb, gs_UNSIGNED ncb,
			     double *A, double  *B, double  *C )	
{
  if ( nca != ncb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }


  cblas_dgemm (CblasRowMajor, // op(A) x op(B)
	       CblasNoTrans,  // op(A): m x k 
	       CblasTrans,  // op(B): k x n,  C: m x n
	       (int)nra,  // m: Number of rows in op(A) (and C)
	       (int)nrb,  // n: Number of columns in op(B) (and C)
	       (int)nca,  // k: Number of columns in op(A), rows in op(B)
	       1.0, 
	       A,
	       (int)nca,  // lda: Row stride of A
	       B, 
	       (int)ncb,  // ldb: Row stride of B
	       0.0, 
	       C, 
	       (int)nrb); // ldc: Row stride of C
}

void gs_mm_mul_blaslapack ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
			    gs_UNSIGNED nrb, gs_UNSIGNED ncb,
			    double *A, double  *B, double  *C )	
{
  if ( nca != nrb )
    {
      gs_info(gv,-2,"Incorrect dimensions for matrix multiplication");
      return;
    }

  cblas_dgemm (CblasRowMajor, // op(A) x op(B)
	       CblasNoTrans,  // op(A): m x k 
	       CblasNoTrans,  // op(B): k x n,  C: m x n
	       (int)nra,  // m: Number of rows in op(A) (and C)
	       (int)ncb,  // n: Number of columns in op(B) (and C)
	       (int)nrb,  // k: Number of columns in op(A), rows in op(B)
	       1.0, 
	       A,
	       (int)nca,  // lda: Row stride of A
	       B, 
	       (int)ncb,  // ldb: Row stride of B
	       0.0, 
	       C, 
	       (int)ncb); // ldc: Row stride of C
}


void gs_solve_Axb_blaslapack ( gs_UNSIGNED nra, gs_UNSIGNED nca, 
			       double *A, double  *x, double  *b, 
			       int* retval  )
{
                           // Solve A x = B for x
                           // A: LDA x N
                           // b: LDB x 1
    lapack_int 
      n    = (lapack_int)nra, // 
      nrhs = 1,               // Number of right hand sides 
      lda  = (lapack_int)nra, // Row stride of A
      ldb  = 1  ,             // Row stride of B
      info;
    
    gs_UNSIGNED dummy; dummy = nca ; // just to avoid compiler warning

    lapack_int ipiv[(lapack_int)nra];

    for (gs_UNSIGNED i=0; i<nra ; i++) x[i] = b[i];

    info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv,
			  x, ldb );

    if( info > 0 ) {
      gs_info(gv,-1,"Singular coefficient matrix in the MME (dgesv)");
      *retval = -2;
      return;
    } 

    *retval = 0;
}



void gs_ludcmp_blaslapack(double *a, gs_UNSIGNED nn, gs_UNSIGNED *indx, double *d)
{
 
    lapack_int 
      n = (lapack_int)nn  ,  
      m = (lapack_int)nn  ,  
      lda = (lapack_int)nn, 
      info;

    lapack_int ipiv[(lapack_int)n];
    
    info = LAPACKE_dgetrf(  LAPACK_ROW_MAJOR, m, n, a, lda, ipiv) ;

    if( info > 0 ) { gs_info(gv,-1,"Singular coefficient matrix in the MME (dgetrf)");}
 
    for (gs_UNSIGNED i=0; i<nn; i++) indx[i] = ipiv[i];

    *d = (double) info;

}

void gs_lubksb_blaslapack(double *a, gs_UNSIGNED nn, gs_UNSIGNED *indx, double *b )
{

    lapack_int 
      n    = (lapack_int)nn,   // 
      nrhs = 1,                // Number of right hand sides 
      lda  = (lapack_int)nn,   // Row stride of A
      ldb  = 1  ,              // Row stride of B
      info;

    lapack_int ipiv[(lapack_int)n];

    for (lapack_int i=0; i<n; i++)  ipiv[i] = (lapack_int)indx[i] ;

    info = LAPACKE_dgetrs( LAPACK_ROW_MAJOR, 'N', n, nrhs, a, lda, ipiv, b, ldb);
}

void gs_luinv_blaslapack (double *a, gs_UNSIGNED nn, gs_UNSIGNED *indx, double *i)
{
 
    lapack_int 
      n = (lapack_int)nn  ,  
      lda = (lapack_int)nn, 
      info;

    lapack_int ipiv[(lapack_int)n];

    for (lapack_int i=0; i<n; i++)  ipiv[i] = (lapack_int)indx[i] ;
  
    info = LAPACKE_dgetri(  LAPACK_ROW_MAJOR,  n, a, lda, ipiv) ;

    memcpy ( i , a , (size_t) n*n * sizeof(double) );
}


#endif


///////////////////////////////////////////////////////////////////////////
// Lapack wrappers for linear algebra
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Chosing linear algebra routines                                       //
///////////////////////////////////////////////////////////////////////////

int  st_matrix_ops_no       = 4 ;
char *st_matrix_names[]    = {"st","mt","gsl","blaslapack"};
int  st_matrix_availabe[5] = {1,1,
  #ifdef use_gsl 
   1,
  #else
   0,
  #endif 
  #ifdef use_blaslapack
   1
  #else
   0
  #endif 

};

#define arg1 gs_UNSIGNED,gs_UNSIGNED,gs_UNSIGNED, gs_UNSIGNED,double*,double*,double*
#define arg2 gs_UNSIGNED,gs_UNSIGNED,double*,double*,double*,int*
#define arg3 double*,gs_UNSIGNED,gs_UNSIGNED*,double*

// Array of pointers to avaliable linear algebra functions
void (*gs_tmm_mul_funcs[])(arg1)  = {
  gs_tmm_mul_st,
  gs_tmm_mul_mt,
  #ifdef use_gsl 
    gs_tmm_mul_gsl,
  #else
    NULL,
  #endif 
  #ifdef use_blaslapack 
    gs_tmm_mul_blaslapack 
  #else
    NULL
  #endif 
} ;


void (*gs_mm_mul_funcs[])(arg1)    = {
  gs_mm_mul_st,
  gs_mm_mul_mt,
  #ifdef use_gsl 
    gs_mm_mul_gsl,
  #else
    NULL,
  #endif 
  #ifdef use_blaslapack 
    gs_mm_mul_blaslapack
  #else
   NULL
  #endif 
} ;

void (*gs_mtm_mul_funcs[])(arg1)    = {
  gs_mtm_mul_st,
  gs_mtm_mul_st,
  #ifdef use_gsl 
    gs_mtm_mul_st,
  #else
    NULL,
  #endif 
  #ifdef use_blaslapack 
    gs_mtm_mul_blaslapack
  #else
   NULL
  #endif 
} ;


void (*gs_solve_Axb_funcs[])(arg2)  = {
  gs_solve_Axb_st,
  gs_solve_Axb_mt,
  #ifdef use_gsl 
    gs_solve_Axb_gsl,
  #else
    NULL,
  #endif 
  #ifdef use_blaslapack 
    gs_solve_Axb_blaslapack,
  #else
    NULL
  #endif 
} ;

void (*gs_ludcmp_funcs[])(arg3)  = {
  gs_ludcmp_st,
  gs_ludcmp_st,
  #ifdef use_gsl 
    gs_ludcmp_gsl,
  #else
    NULL,
  #endif 
  #ifdef use_blaslapack 
    gs_ludcmp_blaslapack
  #else
    NULL
  #endif 
} ;

void (*gs_lubksb_funcs[])(arg3)  = {
  gs_lubksb_st,
  gs_lubksb_st,
  #ifdef use_gsl 
    gs_lubksb_st,  // gsl version is slower
  #else
    NULL,
  #endif 
  #ifdef use_blaslapack 
    gs_lubksb_blaslapack     // openblas version is slower
  #else
    NULL
  #endif 
} ;

void (*gs_luinv_funcs[])(arg3)  = {
  gs_luinv_st,
  gs_luinv_st,
  #ifdef use_gsl 
    gs_luinv_st,   // gsl version is slower
  #else
    NULL,
  #endif 
  #ifdef use_blaslapack 
    gs_luinv_blaslapack
  #else
    NULL
  #endif 
} ;


#undef arg1
#undef arg2
#undef arg2


// Pointer to the linear algebar routine in use
// Default settings

#ifdef use_blaslapack

void (*gs_tmm_mul)()   = gs_tmm_mul_blaslapack;  
void (*gs_mm_mul)()    = gs_mm_mul_blaslapack;   
void (*gs_mtm_mul)()   = gs_mtm_mul_blaslapack;   
void (*gs_solve_Axb)() = gs_solve_Axb_blaslapack;
void (*gs_ludcmp)()    = gs_ludcmp_blaslapack;
void (*gs_lubksb)()    = gs_lubksb_st;
void (*gs_luinv)()     = gs_luinv_blaslapack;

#endif

#ifdef use_gsl

void (*gs_tmm_mul)()   = gs_tmm_mul_gsl;  
void (*gs_mm_mul)()    = gs_mm_mul_gsl;  
void (*gs_mtm_mul)()   = gs_mtm_mul_gsl;   
void (*gs_solve_Axb)() = gs_solve_Axb_gsl;
void (*gs_ludcmp)()    = gs_ludcmp_gsl;
void (*gs_lubksb)()    = gs_lubksb_st;
void (*gs_luinv)()     = gs_luinv_st;

#endif

#if ( (! defined use_gsl ) && ( ! defined use_blaslapack ) )

void (*gs_tmm_mul)()   = gs_tmm_mul_mt;  
void (*gs_mm_mul)()    = gs_mm_mul_mt;   
void (*gs_mtm_mul)()   = gs_mtm_mul_st;   
void (*gs_solve_Axb)() = gs_solve_Axb_mt;
void (*gs_ludcmp)()    = gs_ludcmp_st;
void (*gs_lubksb)()    = gs_lubksb_st;
void (*gs_luinv)()     = gs_luinv_st;

#endif


// Pointer to the linear algebra routine for in parallel regions
// Currently not used.
void (*gs_tmm_mul_p)()   = gs_tmm_mul_st;  
void (*gs_mm_mul_p)()    = gs_mm_mul_st;   
void (*gs_mtm_mul_p)()   = gs_mtm_mul_st;   
void (*gs_solve_Axb_p)() = gs_solve_Axb_st;
void (*gs_ludcmp_p)()    = gs_ludcmp_st;
void (*gs_lubksb_p)()    = gs_lubksb_st;
void (*gs_luinv_p)()     = gs_luinv_st;


void st_set_matrix_ops (char** select_s, char** select_p)
{

  int found_s = 0 , found_p= 0 ;
  for (int i=0; i<st_matrix_ops_no; i++)
    {
      if ( ( 0 == strcmp ( *select_s, st_matrix_names[i] ) ) &&
	   ( 1 == st_matrix_availabe[i]))
	{
	  gs_tmm_mul   = gs_tmm_mul_funcs[i];
	  gs_mm_mul    = gs_mm_mul_funcs[i];
	  gs_mtm_mul   = gs_mtm_mul_funcs[i];
	  gs_solve_Axb = gs_solve_Axb_funcs[i];
	  gs_ludcmp    = gs_ludcmp_funcs[i];
	  gs_lubksb    = gs_lubksb_funcs[i];
	  gs_luinv     = gs_luinv_funcs[i];
	  found_s = 1;
	}
      if ( ( 0 == strcmp ( *select_p, st_matrix_names[i] ) ) &&
	   ( 1 == st_matrix_availabe[i]))
	{
	  gs_tmm_mul_p   = gs_tmm_mul_funcs[i];
	  gs_mm_mul_p    = gs_mm_mul_funcs[i];
	  gs_mtm_mul_p   = gs_mtm_mul_funcs[i];
	  gs_solve_Axb_p = gs_solve_Axb_funcs[i];
	  gs_ludcmp_p    = gs_ludcmp_funcs[i];
	  gs_lubksb_p    = gs_lubksb_funcs[i];
	  gs_luinv_p     = gs_luinv_funcs[i];
	  found_p = 1;
	}
    }
  if ( 0==found_s){
    sprintf( gv->msg, "Matrix operations '%s' not available", *select_s);
    gs_info(gv,-2,gv->msg);
  }  
  else {
    sprintf( gv->msg, "Matrix operations '%s'", *select_s);
    gs_info(gv,0,gv->msg);
  }

  /*
  if ( 0==found_p) {
    sprintf( gv->msg, "Matrix operations %s is not available", *select_p);
    gs_info(gv,0,gv->msg);
  }  
  else {
    sprintf( gv->msg, "Matrix operations for parallel operation %s", *select_p);
    gs_info(gv,0,gv->msg);
  }
  */

}

///////////////////////////////////////////////////////////////////////////
// Chosing linear algebra routines                                       //
///////////////////////////////////////////////////////////////////////////


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void gs_set_num_threads (int * activate)
{

  if ( *activate > MAXTHREADS)  *activate = MAXTHREADS ;
  if ( *activate < 1 )          *activate = 1 ;
  omp_set_num_threads (*activate);
}

void gs_get_num_threads (int * activate)
{
  #pragma omp parallel
  { 
    *activate = omp_get_num_threads ();
  }
  sprintf( gv->msg, "Number of threads %i", *activate);
  gs_info(gv,0,gv->msg);
}

void gs_reset( gs_varset_type *gv )
{
  gs_UNSIGNED  i;

  gv->msg[0]     = '\0';
  gv->MaxIL      = 0; 
  gv->MaxML      = 0; 
  gv->NoInd      = 0;  
  gv->NoMar      = 0;  
  gv->IndName    = realloc (gv->IndName,0); gv->IndName = NULL;
  gv->MarName    = realloc (gv->MarName,0); gv->MarName = NULL;
  gv->All1       = realloc (gv->All1   ,0); gv->All1    = NULL;
  gv->All2       = realloc (gv->All2   ,0); gv->All2    = NULL;
  gv->genv       = realloc (gv->genv   ,0); gv->genv    = NULL;
  gv->NoAll      = realloc (gv->NoAll  ,0); gv->NoAll   = NULL;
  gv->All        = realloc (gv->All    ,0); gv->All     = NULL;
  gv->AllCnt     = realloc (gv->AllCnt ,0); gv->AllCnt  = NULL;
  gv->Pic        = realloc (gv->Pic    ,0); gv->Pic     = NULL;
  gv->Mis        = realloc (gv->Mis    ,0); gv->Mis     = NULL;
  gv->IndMis     = realloc (gv->IndMis ,0); gv->IndMis  = NULL;
  gv->Z          = realloc (gv->Z      ,0); gv->Z       = NULL; 
  gv->y          = realloc (gv->y      ,0); gv->y       = NULL; 
  gv->NoEff      = 0;  
  gv->EffNme     = realloc (gv->EffNme ,0); gv->EffNme  = NULL;
  gv->EffAll     = realloc (gv->EffAll ,0); gv->EffAll  = NULL;
  gv->CmpAll     = realloc (gv->CmpAll ,0); gv->CmpAll  = NULL;
  gv->aa         = 0; 
  gv->aA         = 1;
  gv->AA         = 2;
  gv->an         = 0.5; 
  gv->nn         = 1;
  gv->nA         = 1.5;
  gv->dimu       = 0; 
  gv->dimA       = 0; 
  gv->A          = realloc (gv->A      ,0); gv->A       = NULL;   
  gv->u          = realloc (gv->u      ,0); gv->u       = NULL;   
  gv->b          = realloc (gv->b      ,0); gv->b       = NULL;
  gv->As         = realloc (gv->As     ,0); gv->As      = NULL;
  gv->bs         = realloc (gv->bs     ,0); gv->bs      = NULL;
  gv->GI         = realloc (gv->GI     ,0); gv->GI      = NULL;   
  gv->a          = realloc (gv->a      ,0); gv->a       = NULL;   
  gv->p          = realloc (gv->p      ,0); gv->p       = NULL;
  gv->ladi       = realloc (gv->ladi   ,0); gv->ladi    = NULL;
  gv->est        = realloc (gv->est    ,0); gv->est     = NULL;
  gv->SMsqm      = realloc (gv->SMsqm  ,0); gv->SMsqm   = NULL;
  gv->SMsqe      = realloc (gv->SMsqe  ,0); gv->SMsqe   = NULL;
  gv->SMF        = realloc (gv->SMF    ,0); gv->SMF     = NULL;
  gv->SMp        = realloc (gv->SMp    ,0); gv->SMp     = NULL;
  gv->SMvar      = realloc (gv->SMvar  ,0); gv->SMvar   = NULL;
  gv->sigma      = realloc (gv->sigma,0)  ; gv->sigma   = NULL;
  gv->sigma_e    = 0;
  gv->NLoc       = 0;
  gv->map_S      = realloc (gv->map_S,0)  ; gv->map_S   = NULL;
  gv->map        = realloc (gv->map,  0)  ; gv->map     = NULL;
  gv->mdl        = realloc (gv->mdl,  0)  ; gv->mdl     = NULL;
  gv->Ref1       = realloc (gv->Ref1   ,0); gv->Ref1    = NULL;
  gv->Ref2       = realloc (gv->Ref2   ,0); gv->Ref2    = NULL;
  gv->NHap       = 0;
  gv->Hap_S      = realloc (gv->Hap_S,0);  gv->Hap_S    = NULL;
  gv->h_u        = realloc (gv->h_u,0)  ;  gv->h_u      = NULL;
  gv->h_o        = realloc (gv->h_o,0)  ;  gv->h_o      = NULL;
  gv->h_n        = realloc (gv->h_n,0)  ;  gv->h_n      = NULL;
  gv->NCrs       = 0;
  gv->crs_S      = realloc (gv->crs_S,0)  ; gv->crs_S   = NULL;
  gv->crs        = realloc (gv->crs,  0)  ; gv->crs     = NULL;
  for (i = 0; i < gv->NChrom; i ++ ) { 
    gv->chrom[i].rf = realloc (gv->chrom[i].rf,0);
    gv->chrom[i].q  = realloc (gv->chrom[i].q ,0);  
  }
  gv->NChrom     = 0;
  gv->chrom      = realloc (gv->chrom, 0)  ; gv->chrom  = NULL;
  gv->t0 = realloc(gv->t0,0);gv->t0=NULL;
  gv->t1 = realloc(gv->t1,0);gv->t1=NULL;
  gv->t2 = realloc(gv->t2,0);gv->t2=NULL;
  gv->t3 = realloc(gv->t3,0);gv->t3=NULL;
  gv->t4 = realloc(gv->t4,0);gv->t4=NULL;
  gv->t5 = realloc(gv->t5,0);gv->t5=NULL;
  gv->t6 = realloc(gv->t6,0);gv->t6=NULL;
  gv->t7 = realloc(gv->t7,0);gv->t7=NULL;
  gv->t8 = realloc(gv->t8,0);gv->t8=NULL;
  gv->t9 = realloc(gv->t9,0);gv->t9=NULL;
  gv->u0 = realloc(gv->u0,0);gv->u0=NULL;
}

void gs_init( gs_varset_type *gv )
{
  gv->name[0]    = '\0';
  gv->msg[0]     = '\0';
  gv->info_level = 0;
  gv->MaxIL      = 0; 
  gv->MaxML      = 0; 
  gv->NoInd      = 0;  
  gv->NoMar      = 0;  
  gv->IndName    = NULL;
  gv->MarName    = NULL;
  gv->All1       = NULL;
  gv->All2       = NULL;
  gv->genv       = NULL;
  gv->NoAll      = NULL;
  gv->All        = NULL;
  gv->AllCnt     = NULL;
  gv->Pic        = NULL;
  gv->Mis        = NULL;
  gv->IndMis     = NULL;
  gv->Z          = NULL; 
  gv->y          = NULL; 
  gv->NoEff      = 0;   			                       
  gv->EffNme     = NULL;
  gv->EffAll     = NULL;
  gv->CmpAll     = NULL;
  gv->aa         = 0; 
  gv->aA         = 1;
  gv->AA         = 2;
  gv->an         = 0.5; 
  gv->nn         = 1;
  gv->nA         = 1.5;
  gv->dimu       = 0;                   
  gv->dimA       = 0; 
  gv->A          = NULL;   
  gv->u          = NULL;   
  gv->b          = NULL;
  gv->As         = NULL;   
  gv->bs         = NULL;
  gv->GI         = NULL;   
  gv->a          = NULL;   
  gv->p          = NULL;   
  gv->ladi       = NULL;
  gv->est        = NULL;
  gv->SMsqm      = NULL;
  gv->SMsqe      = NULL;
  gv->SMF        = NULL;
  gv->SMp        = NULL;
  gv->SMvar      = NULL;
  gv->sigma      = NULL;
  gv->sigma_e    = 0;
  gv->NLoc       = 0;
  gv->map_S      = NULL;
  gv->map        = NULL;
  gv->mdl        = NULL;
  gv->Ref1       = NULL;
  gv->Ref2       = NULL;
  gv->NHap       = 0;
  gv->Hap_S      = NULL;
  gv->h_u        = NULL;
  gv->h_o        = NULL;
  gv->h_n        = NULL;
  gv->NCrs       = 0;
  gv->crs_S      = NULL;
  gv->crs        = NULL;
  gv->NChrom     = 0;
  gv->chrom      = NULL;
  gv->t0=NULL;  gv->t1=NULL;  gv->t2=NULL;  gv->t3=NULL;  gv->t4=NULL; 
  gv->t5=NULL;  gv->t6=NULL;  gv->t7=NULL;  gv->t8=NULL;  gv->t9=NULL;
  gv->u0=NULL;
}

/*--------------*/

void gs_info ( gs_varset_type *gv,
	       int level, 
	       char* message )
{
    if ( gv->info_level >= level )
    {
        if ( 1 <= level )        Rprintf( "I (data set '%s'): %s\n",gv->name, message );
        else if (  0 == level )  Rprintf( "M (data set '%s'): %s\n",gv->name, message );
        else if ( -1 == level )  Rprintf( "W (data set '%s'): %s\n",gv->name, message );
        else if ( -2 == level )  Rprintf( "E (data set '%s'): %s\n",gv->name, message );
    }
}

void gs_info_r ( gs_varset_type *gv,
		 int level, 
		 char* message )
{
    if ( gv->info_level >= level )
    {
        if ( 1 <= level )        Rprintf( "I (data set '%s'): %s\r",gv->name, message );
        else if (  0 == level )  Rprintf( "M (data set '%s'): %s\r",gv->name, message );
        else if ( -1 == level )  Rprintf( "W (data set '%s'): %s\r",gv->name, message );
        else if ( -2 == level )  Rprintf( "E (data set '%s'): %s\r",gv->name, message );
    }
}


void gs_r_info(int* level, char** message)
{
  gs_info(gv, *level, *message );
}

void gs_set_info_level(gs_varset_type *gv , int *level )
{
  if ( *level < -2 )  *level = -2;
  if ( *level > 2 )   *level =  2;
  gv->info_level = *level;
  sprintf( gv->msg          ,
	   "Info level %i"  ,
	   gv->info_level     );
  gs_info(gv,2,gv->msg);
}

void gs_get_info_level ( gs_varset_type *gv , int *level )
{
  *level = gv->info_level;
}

gs_varset_type* gs_fdta (char* set_name)
{
  gs_UNSIGNED i;
  gs_varset_type* sp = NULL;

  /* Default data set */
  if ( 0 == strcmp ( set_name , "default" ) ) 
    sp = gv;

  if (NULL == sp)
    {
      /* User defined data set is available */
      for ( i = 0; i < gs_ndta; i ++ )
	if ( 0 == strcmp ( set_name , gs_data[i]->name ) )
	  sp = gs_data[i];
    }

  if (NULL == sp)
    {
      /* New dataset */  
      gs_ndta++;
      gs_data = realloc(gs_data,gs_ndta*sizeof(gs_varset_type*));
      gs_data[gs_ndta-1] = malloc(sizeof(gs_varset_type));
      sp = gs_data[gs_ndta-1];
      
      gs_init(sp);
      strcpy (sp->name,set_name);
      sp->info_level = gv->info_level;
      sprintf(sp->msg,"New data set generated '%s'",sp->name);
      gs_info(sp,0,sp->msg);
    }

  sprintf(sp->msg,"Working on data set '%s'",sp->name);
  gs_info(sp,2,sp->msg);
  
  return sp;

} 


void gs_get_info_level_GV ( int *level,
			    char** set_name )
{
  gs_get_info_level( gs_fdta(*set_name) , level );
}

void gs_set_info_level_GV ( int *level,
			    char** set_name )
{
  gs_set_info_level( gs_fdta(*set_name), level );
}

void gs_set_all_info_levels ( int *level )
{
  gs_set_info_level( gv, level );
  for ( gs_UNSIGNED i = 0 ; i < gs_ndta; i++)
    gs_set_info_level(gs_data[i], level );
}


/* qsort pointers to string comparison function*/
int gs_string_cmp(const void *a, const void *b) 
{ 
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    return strcmp(*ia, *ib);
} 

/* qsort strings  comparison function*/
int gs_strcmp(const void *a, const void *b) 
{ 
    const char *ia = (const char *)a;
    const char *ib = (const char *)b;
    return strcmp(ia, ib);
} 

/* qsort int comparison function */ 
int gs_intcmp(const void *a, const void *b) 
{ 
    const int *ia = (const int *)a; 
    const int *ib = (const int *)b;
    return *ia  - *ib; 
} 

/* qsort short comparison function */ 
int gs_shortcmp(const void *a, const void *b) 
{ 
    const short *ia = (const short *)a; 
    const short *ib = (const short *)b;
    return *ia  - *ib; 
} 

/* qsort double comparison function */ 
int gs_dbl_cmp(const void *a, const void *b) 
{ 
    const double *ia = (const double *)a; 
    const double *ib = (const double *)b;
    return (*ia > *ib) - (*ia < *ib); 
}

void gs_marker_stats_01 ( gs_varset_type* , char** , char** , char** ,
			  int* , int* , int* , int* , int* );

long long int getline (char** , size_t* , FILE* ) ;

void gs_remove_newline(char*);

void gs_determine_alleles (int*, int*, char*);

void gs_read_marker_data_01 ( gs_varset_type *gv ,
			      char** dta_filename,
			      char** ind_filename,
			      char** mar_filename,
			      char** gen_filename,
			      int*   auxfiles,
			      int*   retval        )
{
  gs_UNSIGNED  
    IndSize = 0, /* Lenght of individual name in input file */
    MarSize = 0, /* Length of marker name in input file     */
    AllSize = 0, /* Lenght of allele name in input file     */
    NoLines = 0; /* Number of lines in the input file       */

  FILE *fp;

  char*  line = NULL;   char** linep = &line;
  size_t n    = 0;      size_t* np   = &n;

  char *p ;

  gs_reset(gv);

  /* Determine size of the data set */
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;

  getline ( linep , np , fp ) ;

  while (!feof(fp))
    {
      if ( 0 < getline ( linep , np , fp ) )
	{
	  p = strtok (line,gs_tok_del);
          if ( NULL == p ) {
	    sprintf( gv->msg,"Error in line %lu of the data file.",2 + NoLines);
	    gs_info(gv,-2,gv->msg); *retval=-2;return;
	  }
	  if ( strlen(p) > IndSize) IndSize = strlen(p);

	  p = strtok ('\0',gs_tok_del);
          if ( NULL == p ) {
	    sprintf( gv->msg,"Error in line %lu of the data file.",2 + NoLines);
	    gs_info(gv,-2,gv->msg); *retval=-2;return;
	  }
	  if ( strlen(p) > MarSize) MarSize = strlen(p);

	  p = strtok ('\0',gs_tok_del);
	  
	  if ( (NULL !=p) && (strlen(p) > AllSize)) AllSize = strlen(p);

	  strcpy(line,"\0");
	  NoLines++; 
	}
    }
  fclose(fp);

  if ( (int)IndSize > gv->MaxIL) gv->MaxIL = IndSize;
  gv->MaxML = MarSize;

  /* Messages on data size*/
  sprintf ( gv->msg, "%lu %s", NoLines, "Datalines" );
  gs_info ( gv, 1, gv->msg );
  sprintf(gv->msg,
	  "Description length: Ind %lu, Mar %lu, All %lu", 
	  IndSize, MarSize, AllSize);
  gs_info(gv,1,gv->msg);
  
  /* Read in data */

  gs_UNSIGNED i,j,idx;

  static char (*ind) [1+IndSize] = NULL;
  static char (*mar) [1+MarSize] = NULL;
  static char (*all) [1+AllSize] = NULL;

  ind = realloc ( ind, NoLines*sizeof(char[1+IndSize]) );
  mar = realloc ( mar, NoLines*sizeof(char[1+MarSize]) );
  all = realloc ( all, NoLines*sizeof(char[1+AllSize]) );
  if ( (NULL == ind) || (NULL == mar) || (NULL == all) ) ERR_M;

  /* Read in data lines */
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline ( linep , np , fp ) ;
  for (j = 0; j < NoLines; j++)
    {
      getline ( linep , np , fp ) ;
      gs_remove_newline ( line ) ;
      p = strtok (line,gs_tok_del);
      strcpy(ind[j],p);
      p = strtok ('\0',gs_tok_del);
      strcpy(mar[j],p);
      p = strtok ('\0',gs_tok_del);
      if (NULL==p) strcpy(all[j],"--");
      else strcpy(all[j],p);
    }
  fclose(fp);
  free(line);

  /* Determine individual and marker names */

  /* char *a [NoLines];   */

  char *(*a);
  a = malloc ( NoLines*sizeof(char*) ) ;
  if (NULL == a) ERR_M;

  for (j = 0 ; j < NoLines; j++)  a[j] = ind[j];
  qsort( a, NoLines, sizeof(char *), gs_string_cmp );

  gv->NoInd  = 1;
  gv->IndName = realloc ( gv->IndName, 1*sizeof(char[gs_MAXNME]) );
  strcpy ( gv->IndName[0], a[0] );

  for (j = 1 ; j < NoLines; j++) 
    {
      if ( 0 != strcmp ( a[j] , a[j-1] ) )
	{
	  gv->NoInd++;
	  gv->IndName = realloc(gv->IndName,
				gv->NoInd * sizeof(char[gs_MAXNME]));
	  strcpy( gv->IndName[gv->NoInd-1] , a[j] );
	}
    }

  for (j = 0 ; j < NoLines; j++)  a[j] = mar[j];
  qsort( a, NoLines, sizeof(char *), gs_string_cmp );

  gv->NoMar  = 1;
  gv->MarName = realloc ( gv->MarName, 1*sizeof(char[gs_MAXNME]) );
  strcpy ( gv->MarName[0], a[0] );

  for (j = 1; j < NoLines; j++)
    {
      if ( 0 != strcmp ( a[j] , a[j-1] ) )
	{
	  gv->NoMar++;
	  gv->MarName = realloc(gv->MarName,
				gv->NoMar * sizeof(char[gs_MAXNME]));
	  strcpy( gv->MarName[gv->NoMar-1] , a[j] );
      }
    }

  free(a);

  /* Allocate memory for the genotypic data */
  gv->All1 =  realloc(gv->All1, gv->NoMar * gv->NoInd * sizeof(int));
  gv->All2 =  realloc(gv->All2, gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == gv->All1)  || (NULL == gv->All2) ) ERR_M;

  for (i=0; i < gv->NoMar; i++)
    for (j=0; j < gv->NoInd; j++)
      {
	idx = i*gv->NoInd + j;
	gv->All1[idx] = -1;
	gv->All2[idx] = -1;
      }

  /* Read data from list to matrix */

  #pragma omp parallel 
  {
  int a=-1 ,b=-1;  int* za=&a;  int* zb=&b;
  gs_UNSIGNED i,j,idx;

  #pragma omp for 
  for (long long int k = 0; k < (long long int)NoLines; k++){
    for (i = 0; i < gv->NoMar; i++){
      if ( 0 == strcmp (mar[k], gv->MarName[i]) )
	{
	  for (j=0; j <gv->NoInd; j++){
	    if ( 0 == strcmp (ind[k], gv->IndName[j])) 
	      {
		gs_determine_alleles (za,zb,all[k]);
		idx = i*gv->NoInd + j;
		if ( -1 == gv->All1[idx] ) 
		  gv->All1[idx] = a;
		gv->All2[idx] = b;
	      }
	  }
	} /* if mar[k] == MarName[i]  */
    } /* for i = 0 .. gv->NoMar    */
  } /* for k = 0 .. NoLines  */

  } /* parallel */

  free(mar); mar = NULL;
  free(ind); ind = NULL;
  free(all); all = NULL;

  int Srv=0; int* Sretval = &Srv;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       auxfiles,auxfiles,auxfiles,
		       auxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;

}


void gs_read_marker_data_01_GV ( char** dta_filename,
			      char** ind_filename,
			      char** mar_filename,
			      char** gen_filename,
			      int*   auxfiles,
			      int*   retval,
			      char** set_name)
{
  gs_read_marker_data_01 ( gs_fdta(*set_name),
			dta_filename,
			ind_filename,
			mar_filename,
			gen_filename,
			auxfiles,
			retval         ) ;

}


/* Design matrix of type 1: */
/* biallelic markers */
/* one column per marker */
/* design matrix consists of counts of the greater allele */
/* greater allele is the 'effect allele' */
/* missing equals heterozyous */

/* effect allele:        b */
/* complementray allele: a */
/* everything else:      c */

/* c is missing or alleles that are present in the prediction set but not */
/* the estimation set */

/* Values of the design matrix */

/* a:    0        */
/* c:   1/2       */
/* b:    1        */

/* aa   ab    bb  */
/* 0     1     2  */

/* ac   cc    bc   */
/* 1/2   1    3/2  */

/* Storage of the complementray allele is required to distinguish  */
/* between  */
/* aa ac cc and */
/* ab bc */

void gs_build_Z_01 ( gs_varset_type *gv ,
		     char** out_filename,
		     int*   auxfiles,
		     int*   retval)
{

  gs_UNSIGNED i, m, idx_all, idx_Z;
  FILE *outfp;

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    *retval = -2; return ;
  }

  /* Allocate memory for the design matrix */
  gv->Z =  realloc( gv->Z, gv->NoInd * gv->NoMar * sizeof(double) );

  /* Check wheter EffAlls are already defined */
  int EffAllDefined = 1;
  if ( NULL ==  gv->EffAll ) EffAllDefined = 0; 

  /* Allocate memory for the effect descriptions */
  gv->NoEff =  1 + gv->NoMar;
  gv->EffNme = realloc (gv->EffNme, gv->NoEff * sizeof(char[gs_MAXNME]) );
  gv->EffAll = realloc (gv->EffAll, gv->NoEff * sizeof(int) );
  gv->CmpAll = realloc (gv->CmpAll, gv->NoEff * sizeof(int) );

  if ( (NULL==gv->Z) || (NULL==gv->EffNme) || (NULL==gv->EffAll) 
       || (NULL==gv->CmpAll)) ERR_M;

  /* Names of effects */
  strcpy(gv->EffNme[0],"mu");
  gv->EffAll[0] = 0;
  gv->CmpAll[0] = 0;
  for ( m = 1; m < gv->NoEff; m++)
    strcpy(gv->EffNme[m],gv->MarName[m-1]);		       
  
  /* Initialize design matrix */

  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar; m++) {
      idx_Z   = i*gv->NoMar + m ;
      gv->Z[idx_Z] =  gv->aA;
    }

  int a, b, lar_all, sml_all,
    mono_warning ;

  mono_warning = 0;

  for ( m = 0; m < gv->NoMar; m++)
    { 
      if ( NULL == gv->u ) {
	// Use the larger allele as EffAll
	lar_all = 0;
	for ( i = 0; i < gv->NoInd ; i++){
	  idx_all = m*gv->NoInd + i ;
	  if (gv->All1[idx_all] > lar_all) lar_all = gv->All1[idx_all];
	  if (gv->All2[idx_all] > lar_all) lar_all = gv->All2[idx_all];
	}
	/* Determine the complementary allele for each marker */
	sml_all = lar_all ;
	for ( i = 0; i < gv->NoInd ; i++){
	  idx_all = m*gv->NoInd + i ;
	  if ( (gv->All1[idx_all] < lar_all) && (gv->All1[idx_all] > 0) )
	    sml_all = gv->All1[idx_all];
	  if ( (gv->All2[idx_all] < lar_all) && (gv->All2[idx_all] > 0) )
	    sml_all = gv->All2[idx_all];
	}
	/* Check whether marker is biallelic */
	if ( (lar_all==sml_all) || (lar_all<2)  ) mono_warning = 1;
        /* Determine genotype */
	for ( i = 0; i < gv->NoInd ; i++){
	  idx_all = m*gv->NoInd + i ;
	  if ( ( (gv->All1[idx_all] != lar_all) && 
		 (gv->All1[idx_all] != sml_all) && 
		 (gv->All1[idx_all] != -1 ) )          ||
	       ( (gv->All2[idx_all] != lar_all) && 
		 (gv->All2[idx_all] != sml_all) && 
		 (gv->All2[idx_all] != -1 ) )              )
	    {
	      sprintf( gv->msg, "More than two alleles: %s",gv->MarName[m] );
	      gs_info(gv,-2,gv->msg);
	      *retval=-2;return;
	    }
	}
	/* Save allele for which effects will be determined */
	gv->EffAll[m+1] = lar_all;
	gv->CmpAll[m+1] = sml_all;
      } else {
	// Use EffAll and CmpAll from previous run 
	lar_all = gv->EffAll[m+1] ;
	sml_all = gv->CmpAll[m+1] ;
      }

     /* Construct design matrix */
      for ( i = 0; i < gv->NoInd ; i++)
	{
	  idx_all = m*gv->NoInd + i ;

	  a = gv->All1[idx_all] ;
          b = gv->All2[idx_all] ;

	  idx_Z   = i*gv->NoMar + m ; 

	  if ( (a==lar_all) && (b==lar_all) ) { 
	    gv->Z[idx_Z] = gv->AA; 
	  }
	  else	{
	    if ( (a==sml_all) && (b==sml_all) )  {
	      gv->Z[idx_Z] = gv->aa;
	    }
	    else {
	      gv->Z[idx_Z] = gv->aA ;
	    }
	  }
	}
    }
  
if ( 1 == mono_warning  )
	{
	  sprintf( gv->msg, "Monomorphic markers in the design matrix" );
	  gs_info(gv,-1,gv->msg);
	}

   

  /* Write design matrix to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( m = 0; m < gv->NoMar; m++)
	fprintf(outfp,"%s.%i ",gv->EffNme[m+1],gv->EffAll[m+1]);
      fprintf(outfp,"\n");

      for ( i = 0; i < gv->NoInd ; i++)
	{
	  fprintf(outfp,"%s ",gv->IndName[i]);
	  for ( m = 0; m < gv->NoMar; m++)
	    {
	      idx_Z = i*gv->NoMar + m ;
	      fprintf(outfp,"%f ",gv->Z[idx_Z]);
	    }
	  fprintf(outfp,"\n");
	}

      fclose(outfp);  
  }

  *retval = 0; return;

}

void gs_build_Z_01_GV ( 
		       char** out_filename        ,
		       int*   auxfiles            ,
		       int*   retval              ,
		       char** set_name                               
		      )
{
  gs_build_Z_01 ( gs_fdta(*set_name) ,
		  out_filename,
		  auxfiles,
		  retval         );
}

void gs_neg_effall_01 (  gs_varset_type *gv ,
			 int*    retval         )
{
  gs_UNSIGNED m;
  int c;

  if ( NULL == gv->u ) {
      gs_info(gv,-2, "No effects estimated");
      *retval = -2; return ;
  }

  for ( m = 1; m < gv->NoEff; m++) {
    if ( gv->u[m] > 0) {
      c = gv->EffAll[m];
      gv->EffAll[m] = gv->CmpAll[m];
      gv->CmpAll[m] = c;
    }
    gv->u[m] = 0;
  }

  *retval = 0; return;
}

void gs_neg_effall_01_GV (  int*    retval,
			    char** set_name         )
{
  gs_neg_effall_01 ( gs_fdta(*set_name),
		     retval         ) ;
}

void gs_pos_effall_01 (  gs_varset_type *gv ,
			 int*    retval         )
{
  gs_UNSIGNED m;
  int c;

  if ( NULL == gv->u ) {
      gs_info(gv,-2, "No effects estimated");
      *retval = -2; return ;
  }

  for ( m = 1; m < gv->NoEff; m++) {
    if ( gv->u[m] < 0) {
      c = gv->EffAll[m];
      gv->EffAll[m] = gv->CmpAll[m];
      gv->CmpAll[m] = c;
    }
    gv->u[m] = 0;
  }

  *retval = 0; return;
}

void gs_pos_effall_01_GV (  int*    retval,
			    char** set_name         )
{
  gs_pos_effall_01 ( gs_fdta(*set_name),
		     retval         ) ;
}


long long int getline(char **, size_t *, FILE *);

void gs_remove_newline(char*);

void gs_read_performance_data_01 ( gs_varset_type *gv ,
				   char** dta_filename,
				   char** out_filename,
				   int*   auxfiles,
				   int*   retval         )
{
  gs_UNSIGNED  
    IndSize = 0, /* Lenght of individual name in input file  */
    PerSize = 0, /* Length of perfromance data in input file */
    NoLines = 0; /* Number of lines in the input file       */

  FILE *fp;
  FILE *outfp;

  char*  line = NULL;   char** linep = &line;
  size_t n    = 0;      size_t* np   = &n;

  char *p;
  long long int c;


  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  /* Determine size of the data set */

  if ( (fp = fopen(*dta_filename,"r")) == NULL ) ERR_F;
  getline (linep,np,fp);
  while (!feof(fp))
    {
      c = getline (linep,np,fp);
      if ( c > 0 ) 
	{
	  p = strtok (line,gs_tok_del);
	  if ( strlen(p) > IndSize) IndSize = strlen(p);
	  p = strtok ('\0',gs_tok_del);
	  if ( strlen(p) > PerSize) PerSize = strlen(p);
	  strcpy(line,"\0");
	  NoLines++;
	}
    }
  fclose(fp);

  /* Messages on data size*/

  sprintf(gv->msg,"%lu %s",NoLines, "Datalines");
  gs_info(gv,1,gv->msg);
  sprintf(gv->msg,
	  "Description length: Ind %lu, Performance %lu", 
	  IndSize, PerSize);
  gs_info(gv,1,gv->msg);
  
  /* Read in data */

  gs_UNSIGNED i,j,k,m,i_n,idx;

  static char (*ind) [1+IndSize] = NULL;
  static char (*per) [1+PerSize] = NULL;
  
  ind = realloc ( ind, NoLines*sizeof(char[1+IndSize]) );
  per = realloc ( per, NoLines*sizeof(char[1+PerSize]) );
  if ( ( NULL== per ) || ( NULL== per )) ERR_M;

  /* Read in data lines */
  if ( (fp = fopen(*dta_filename,"r")) == NULL ) ERR_F;
  getline (linep,np,fp);
  for (j = 0; j < NoLines; j++)
    {
      getline (linep,np,fp);
      p = strtok (line,gs_tok_del);
      strcpy(ind[j],p);
      p = strtok ('\0',gs_tok_del);
      strcpy(per[j],p);
    }
  fclose(fp);
  free(line);

  // Allocate memory
  gv->y = realloc( gv->y, gv->NoInd*sizeof(double) );
  if ( NULL== gv->y ) ERR_M;

  // Copy performance data
  for (i = 0; i < gv->NoInd; i++) 
    {
      gv->y[i] = -99999;
      for (k = 0; k < NoLines; k++)
	{
	  if ( 0 == strcmp ( gv->IndName[i] , ind[k] ) )
	    gv->y[i] = strtod(per[k],NULL);
	}	
    }

  // Count missing performane data
  gs_UNSIGNED miss_phen = 0;
  for (i = 0; i < gv->NoInd; i++)
    if (-99999 == gv->y[i] ) miss_phen++;
 
  sprintf ( gv->msg, "No. of performance data: %lu", gv->NoInd-miss_phen);
  gs_info(gv,0,gv->msg);
  
  if ( 0 < miss_phen)
    {
      // Reduce marker data if performance data are missing
      {
	/* Temporary pointers to the map */
	gs_UNSIGNED         NLoc  = gv->NLoc; 
	gs_mappoint_type  * map_S = gv->map_S; gv->map_S = NULL; 
	gs_mappoint_type ** map   = gv->map;   gv->map   = NULL;
	gs_UNSIGNED * mdl         = gv->mdl;   gv->mdl   = NULL;

	// New storage
	gs_UNSIGNED NoInd = gv->NoInd - miss_phen; 
	gs_UNSIGNED NoMar = gv->NoMar; 

	char (*IndName) [gs_MAXNME] = malloc(NoInd*sizeof(char[gs_MAXNME])); 
	char (*MarName) [gs_MAXNME] = gv->MarName; gv->MarName = NULL; 

	int * All1  =  malloc( NoMar * NoInd * sizeof(int) );
	int * All2  =  malloc( NoMar * NoInd * sizeof(int) );

	double *y = malloc( NoInd * sizeof(double) );

	if ( NULL==All1  || NULL==All2 || NULL==IndName || NULL==y ) ERR_M;

	// Copy marker names and phennotypes to new storage
	gs_UNSIGNED CountInd = 0;
	for ( i = 0; i < gv->NoInd; i++  )
	  if  (-99999 != gv->y[i])  
	    {
	      strcpy ( IndName[CountInd] , gv->IndName[i]);
	      y[CountInd] = gv->y[i];
	      CountInd ++;
	    }

	// Copy marker data matrix
	for ( m=0; m < gv->NoMar; m++ )
	  {
	    CountInd = 0;
	    for ( i = 0; i < gv->NoInd; i++  )
	      if  (-99999 != gv->y[i])
		{
		  i_n = m * NoInd     + CountInd;
		  idx = m * gv->NoInd + i;
		  All1[i_n] = gv->All1[idx]; 
		  All2[i_n] = gv->All2[idx];
		  CountInd ++;
		}
	  }
  
	// Replace old with new data  
	gs_reset(gv);
	gv->NoInd   = NoInd;        
	gv->NoMar   = NoMar; 
	gv->IndName = IndName;    
	gv->MarName = MarName; 
	gv->All1    = All1;
	gv->All2    = All2;
	gv->y       = y;

	/* Write back the map */
	gv->NLoc  = NLoc ;  gv->map_S = map_S;
	gv->map   = map  ;  gv->mdl   = mdl  ;

	/* Determine marker statistics */
	int Srv=0; int* Sretval = &Srv;
	char*empty="";
	char ** ind_filename = &empty; char ** mar_filename = &empty;
	char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
	gs_marker_stats_01 ( gv ,
			     ind_filename,
			     mar_filename,
			     gen_filename,
			     Sauxfiles,Sauxfiles,Sauxfiles,
			     Sauxfiles,
			     Sretval        );
	if (-2 == *Sretval) { *retval = -2; return; }
      }
    }

  /* Write observation vector to output file */
  if (1 == *auxfiles) {

    if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

    fprintf(outfp,"y\n");
    for ( i = 0; i < gv->NoInd ; i++)
      fprintf(outfp,"%s %f\n",gv->IndName[i], gv->y[i]);
    
    fclose(outfp);  
  }

  *retval = 0; return ;

}

void gs_read_performance_data_01_GV (  char** dta_filename,
				       char** out_filename,
				       int*   auxfiles,
				       int*   retval,
				       char** set_name   )
{
  gs_read_performance_data_01 ( gs_fdta(*set_name),
				dta_filename,
				out_filename,
				auxfiles,
				retval ) ;
}

void gs_return_performance_data_01 ( gs_varset_type *gv ,
				     char** out_filename,
				     int*   auxfiles,
				     int*   retval         )
{
  gs_UNSIGNED i;

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return;
  }

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  /* Write observation vector to output file */
  if (1 == *auxfiles) {

    FILE *outfp;

    if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

    fprintf(outfp,"y\n");
    for ( i = 0; i < gv->NoInd ; i++)
      fprintf(outfp,"%s %f\n",gv->IndName[i], gv->y[i]);
    
    fclose(outfp);  
  }

  *retval = 0; return ;

}

void gs_return_performance_data_01_GV ( char** out_filename,
					int*   auxfiles,
					int*   retval,
					char** set_name   )
{
  gs_return_performance_data_01 ( gs_fdta(*set_name),
				  out_filename,
				  auxfiles,
				  retval ) ;
}



void gs_lambda_const_01 ( gs_varset_type *gv ,
			  double* lambda, 
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  gs_UNSIGNED i,m;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if ( NULL== gv->ladi )  ERR_M;

  for ( i = 0 ; i < gv->NoMar ; i++)  gv->ladi[i] = * lambda;

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_const_01_GV ( double* lambda      , 
			     char**  out_filename,
			     int*    auxfiles    , 
			     int*    retval      ,
			     char**  set_name     )

{
  gs_lambda_const_01 ( gs_fdta(*set_name),
		       lambda, 
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}


void gs_lambda_const_02 ( gs_varset_type *gv ,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  gs_UNSIGNED i,m;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available, run 'build.Z' first");
    *retval = -2; return ;
  }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if ( NULL== gv->ladi )  ERR_M;

  double lambda =  ( 1 / *hsq -1 ) * (double)gv->NoMar; 

  for ( i = 0 ; i < gv->NoMar ; i++)  gv->ladi[i] = lambda;

  /* Write output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_const_02_GV ( double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval         ,
			     char** set_name         )
{
  gs_lambda_const_02 ( gs_fdta(*set_name),
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}

void gs_single_marker_aov_01 ( gs_varset_type *gv  ,
			       double* alpha       ,
			       char** out_filename ,
			       int*   auxfiles     , 
			       int*   retval         )
{
  /* single marker aov with a random model to estimate the
     variance component of one marker.
     0 1 2 design matrix
  */

  gs_UNSIGNED i,m,idx_Z;
  FILE *outfp;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2, "No Z matrix available, run 'build.Z' first");
    *retval = -2; return;
  }

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return;
  }

  /* Check the design matrix Z for 0,1,2 */
  int ext = 0; double zz;
  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar; m++) {
      idx_Z   = i*gv->NoMar + m ;
      zz = gv->Z[idx_Z] ;
      if ( (0!=zz)&&(1!=zz)&&(2!=zz) ) {ext = 1;break;} ;
    }
  
  if ( 1 == ext ) {
    gs_info(gv,-2,
      "Incorrect design matrix for single marker anova");
    *retval = -2; return;
  }

  /* Allocate memory for the results */
  gv->SMp   =  realloc( gv->SMp, gv->NoMar   * sizeof(double) );
  gv->SMF   =  realloc( gv->SMF, gv->NoMar   * sizeof(double) );
  gv->SMsqm =  realloc( gv->SMsqm, gv->NoMar * sizeof(double) );
  gv->SMsqe =  realloc( gv->SMsqe, gv->NoMar * sizeof(double) );
  gv->SMvar =  realloc( gv->SMvar, gv->NoMar * sizeof(double) );
  if ( ( NULL ==  gv->SMp)   || ( NULL ==  gv->SMF ) ||
       ( NULL ==  gv->SMsqm) || ( NULL ==  gv->SMsqe ) ||
       ( NULL ==  gv->SMvar )                             ) ERR_M;

  /* Initialize vector for p-values */
  for ( m = 0; m < gv->NoMar ; m++) gv->SMp [m] = 0;

  double k = (double) (gv->NoInd);
  double sy=0, syy=0,
    SQT=0, SQM=0, SQE=0, mm, kk;
  
  double sy0, sy1, sy2, 
         ny0, ny1, ny2;

  /* Calculate single marker anova */

  for ( i = 0; i < gv->NoInd; i++ )
    {
      sy  += gv->y[i] ;
      syy += gv->y[i] * gv->y[i] ;
    }
  SQT = syy - sy*sy/k;

  for ( m = 0; m < gv->NoMar ; m++) 
    {
      sy0=0; sy1=0; sy2=0;
      ny0=0; ny1=0; ny2=0;

      for ( i = 0; i < gv->NoInd; i++ )
	{
	  idx_Z   = i*gv->NoMar + m ;
	  if (0==gv->Z[idx_Z]) { ny0++; sy0 += gv->y[i]; } else
	  if (2==gv->Z[idx_Z]) { ny2++; sy2 += gv->y[i]; } else
	  if (1==gv->Z[idx_Z]) { ny1++; sy1 += gv->y[i]; } 
	}

      mm  =    ( 0==ny0 ? 0 : 1)
             + ( 0==ny1 ? 0 : 1)
	     + ( 0==ny2 ? 0 : 1 );

      SQE = syy - ( 0==ny0 ? 0 : sy0*sy0/ny0)
	        - ( 0==ny1 ? 0 : sy1*sy1/ny1)
    	        - ( 0==ny2 ? 0 : sy2*sy2/ny2) ;

      SQM = SQT - SQE ; 

      gv->SMsqm [m] = SQM ;
      gv->SMsqe [m] = SQE ;

      if ( (2>mm) || ((k-mm)<1) ) {
	  gv->SMvar[m] = FLT_MIN;
	  gv->SMF[m] = 0 ;
	  gv->SMp[m] = 1 ;
      } 
      else {
	  kk = ( 2/(mm-1) ) * ( k - (ny0*ny0 + ny1*ny1 + ny2*ny2)/k );
	  gv->SMvar [m] = (SQM/(mm-1) - SQE/(k-mm)) / kk;
	  if (gv->SMvar[m] <= 0 ) gv->SMvar[m] = FLT_MIN;
	  gv->SMF[m] = ( SQM/(mm-1) ) / ( SQT/(k-mm));
	  gv->SMp[m] = pf( gv->SMF[m] ,(int)(mm-1),(int)(k-mm),0,0);
	  if ( gv->SMp[m] > *alpha )  gv->SMvar[m] = FLT_MIN;
      }
      
    }

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"SQM SQE var F p\n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f %f %f %f\n",
		gv->MarName[m],
		gv->SMsqm [m],
		gv->SMsqe [m],
		gv->SMvar [m],
		gv->SMF[m],
		gv->SMp[m]
		);

      fclose(outfp);  
  }

  *retval = 0; return;

}

void gs_single_marker_aov_01_GV ( double* alpha       ,
				  char** out_filename,
				  int*   auxfiles    , 
				  int*   retval      ,
				  char** set_name     )
{
  gs_single_marker_aov_01 ( gs_fdta(*set_name) ,
			    alpha       ,
			    out_filename,
			    auxfiles, 
			    retval         );
}


void gs_lambda_aov_01 (  gs_varset_type *gv ,
			 double* hsq,
                	 double* alpha,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;

  double Svar = 0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    {
      gv->ladi[m] = ( ( 1 / *hsq) - 1 ) * Svar / gv->SMvar[m];
    }

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_aov_01_GV ( double* hsq           ,
			   double* alpha         ,
			   char**  out_filename  ,
			   int*    auxfiles      , 
			   int*    retval        ,
			   char** set_name        )
{
  gs_lambda_aov_01 ( gs_fdta(*set_name),
		     hsq,
		     alpha,
		     out_filename,
		     auxfiles, 
		     retval         ) ;
}


void gs_single_marker_reg_01 ( gs_varset_type *gv  ,
			       char** out_filename ,
			       int*   auxfiles     , 
			       int*   retval         )
{
  gs_UNSIGNED i,m,idx_Z;
  FILE *outfp;


  if ( NULL == gv->Z ) {
    gs_info(gv,-2, "No Z matrix available, run 'build.Z' first");
    *retval = -2; return;
  }

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return;
  }

  /* Allocate memory for the results */
  gv->SMp   =  realloc( gv->SMp, gv->NoMar   * sizeof(double) );
  gv->SMF   =  realloc( gv->SMF, gv->NoMar   * sizeof(double) );
  gv->SMsqm =  realloc( gv->SMsqm, gv->NoMar * sizeof(double) );
  gv->SMsqe =  realloc( gv->SMsqe, gv->NoMar * sizeof(double) );
  if ( ( NULL ==  gv->SMp)   || ( NULL ==  gv->SMF ) ||
       ( NULL ==  gv->SMsqm) || ( NULL ==  gv->SMsqe ) ) ERR_M;

  /* Initialize vector for p-values */
  for ( m = 0; m < gv->NoMar ; m++) gv->SMp [m] = 0;

  double k = (double) (gv->NoInd);
  double sy=0, syy=0, sx=0, sxx=0, sxy=0,
    b,
    SQT=0, SQM=0, SQE=0, F, p;
  
  /* Calculate single marker regression */

  for ( i = 0; i < gv->NoInd; i++ )
    {
      sy  += gv->y[i] ;
      syy += gv->y[i] * gv->y[i] ;
    }
  SQT = syy - sy*sy/k;

  for ( m = 0; m < gv->NoMar ; m++) 
    {
      sx=0; sxx=0; sxy=0;
      for ( i = 0; i < gv->NoInd; i++ )
	{
	  idx_Z   = i*gv->NoMar + m ;
	  sx  += gv->Z[idx_Z];
	  sxx += gv->Z[idx_Z] * gv->Z[idx_Z] ;
	  sxy += gv->Z[idx_Z] * gv->y[i]     ;
	}

      b = ( sxy - sx*sy/k ) / 
          ( sxx - sx*sx/k) ;
    
      SQM = b * (sxy - sx*sy/k);
      SQE = SQT - SQM;
      F   = SQM / (SQE/(k-2)) ;
      p   = pf(F,1,(int)gv->NoInd,0,0);

      gv->SMsqm [m] = SQM ;
      gv->SMsqe [m] = SQE ;
      gv->SMF   [m] = F   ;
      gv->SMp   [m] = p   ; 

    }

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"SQM SQE F p\n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f %f %f\n",
		gv->MarName[m],
		gv->SMsqm [m],
		gv->SMsqe [m],
		gv->SMF[m],
		gv->SMp[m]
		);

      fclose(outfp);  
  }

  *retval = 0; return;

}

void gs_single_marker_reg_01_GV ( char** out_filename,
				  int*   auxfiles    , 
				  int*   retval      ,
				  char** set_name     )
{
  gs_single_marker_reg_01 ( gs_fdta(*set_name) ,
			    out_filename,
			    auxfiles, 
			    retval         );
}


void gs_lambda_reg_01 (  gs_varset_type *gv ,
			 double* hsq,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_reg_01 ( gv, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;

  double SSQM = 0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    SSQM += gv->SMsqm[m] ;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] = ( ( 1 / *hsq) - 1 ) * SSQM / gv->SMsqm[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_reg_01_GV ( double* hsq           ,
			   char**  out_filename  ,
			   int*    auxfiles      , 
			   int*    retval        ,
			   char** set_name        )
{
  gs_lambda_reg_01 ( gs_fdta(*set_name),
		     hsq,
		     out_filename,
		     auxfiles, 
		     retval         ) ;
}


void gs_mme_coeff_01 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->y ) {
      gs_info(gv,-2,
        "No performance data available, run 'read.performance.data' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  /* Allocate memory for the A matrix */

  gv->dimA = 1 + gv->NoMar;
  gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
  if ( NULL == gv->A ) ERR_M;

  for ( i = 0; i < gv->dimA ; i++) 
    for ( j = 0; j < gv->dimA ; j++) 
	gv->A [ i*gv->dimA + j ] =  0;

  /* Initialize the 1 and I matrix */

  double *one_vec; 

  one_vec = malloc ( gv->NoInd * sizeof(double) );
  if ( NULL == one_vec ) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++) one_vec[i] = 1;

  /* Init part 1*/
  
  gv->A [0] = (double) gv->NoInd;

  /* Init part 2 and 3*/

  double *tmp1;
  tmp1 = malloc ( gv->NoMar * sizeof(double) );
  if ( NULL == tmp1 ) ERR_M;

  gs_tmm_mul(gv->NoInd, 1,
  	     gv->NoInd, gv->NoMar,
  	     one_vec, gv->Z, tmp1);

  for ( j = 0; j < gv->NoMar; j++)  
    {
      gv->A [     0*gv->dimA + (1+j) ] = tmp1[j];
      gv->A [ (j+1)*gv->dimA + 0     ] = tmp1[j];
    }

  /* Init part 4*/

  /* implementation of tmp2 as nr-matrix would save space */
  double *tmp2;
  tmp2 = malloc ( gv->NoMar * gv->NoMar *sizeof(double) );
  if ( NULL == tmp2 ) ERR_M;

  gs_tmm_mul(gv->NoInd, gv->NoMar,
  	     gv->NoInd, gv->NoMar,
  	     gv->Z, gv->Z, tmp2);
  
  /* Build A */
  for ( i = 0; i < gv->NoMar; i++)
    for ( j = 0; j < gv->NoMar; j++)
      {
	gv->A [ (1+i)*gv->dimA + (1+j) ] = tmp2[ i*gv->NoMar + j];
      }

  /* Store A without applied lambda */

  if (1) 
    {
      gv->As =  realloc( gv->As, gv->dimA * gv->dimA *sizeof(double) );
      if ( NULL == gv->As ) ERR_M;
      for ( i = 0; i < gv->dimA; i++)
	for ( j = 0; j < gv->dimA; j++)
	  {
	    gv->As [ i*gv->dimA + j ] = gv->A [ i*gv->dimA + j ] ;
	  }
    }

  /* Apply lambda */
  for ( i = 1; i < gv->dimA; i++)
	gv->A [ i*gv->dimA + i ]  += gv->ladi[i-1] ;

  free(one_vec);
  free(tmp1);
  free(tmp2);

  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mme_coeff_01_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_mme_coeff_01 (gs_fdta(*set_name),
		   out_filename,
		   auxfiles, 
		   retval        ) ;
}


void gs_build_V_01 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  double *tmp2;
  tmp2 = malloc ( gv->NoInd * gv->NoInd *sizeof(double) );
  if ( NULL == tmp2 ) ERR_M;

  gs_mtm_mul(gv->NoInd, gv->NoMar,
  	     gv->NoInd, gv->NoMar,
  	     gv->Z, gv->Z, tmp2);
  
  /* Write V matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->NoInd; i++) {
	  for ( j = 0; j < gv->NoInd; j++) {
	      idx_A   = i*gv->NoInd + j ;
	      fprintf(outfp,"%f ",tmp2 [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  free(tmp2);

  *retval = 0; return ;
}

void gs_build_V_01_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_build_V_01 (gs_fdta(*set_name),
		   out_filename,
		   auxfiles, 
		   retval        ) ;
}



void gs_mme_restcoeff_01 ( gs_varset_type* gv ,
			   char** out_filename,
			   int*   auxfiles, 
			   int*   retval        )
{

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  if ( NULL == gv->As ) 
    {
      /* No stored value availabe, Calculate */
      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }
  else
    {
      /* Restore A */
      gv->dimA = 1 + gv->NoMar;
      gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
      if ( NULL == gv->A ) ERR_M;

      memcpy ( gv->A , gv->As, (size_t)gv->dimA * 
                               (size_t)gv->dimA * sizeof(double));
      
      /* Apply lambda */
      for ( i = 1; i < gv->dimA; i++)
	gv->A [ i*gv->dimA + i ]  += gv->ladi[i-1] ;
    }

  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mme_restcoeff_01_GV ( char** out_filename,
			 int*   auxfiles, 
			 int*   retval        )
{
  gs_mme_restcoeff_01 (gv,
		       out_filename,
		       auxfiles, 
		       retval        ) ;
}




void gs_mme_rhs_01 ( gs_varset_type *gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  gs_UNSIGNED i,j, dim_b;
  FILE *outfp;

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available, run 'build.Z' first");
    *retval = -2; return ;
  }

  /* Allocate memory for the b and initialize */
  dim_b = 1 + gv->NoMar;
  gv->b =  realloc( gv->b, dim_b * sizeof(double) );
  if ( NULL == gv->b ) ERR_M;
  for ( i = 0; i < dim_b ; i++) gv->b [i] =  0;

  double *one_vec;
  one_vec = malloc ( gv->NoInd * sizeof(double) );
  if ( NULL == one_vec ) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++) one_vec[i] = 1;

  /* Init part 1*/
  
  double *tmp ;
  tmp = malloc ( 1 * sizeof(double) );
  if ( NULL == tmp ) ERR_M;

  gs_tmm_mul(gv->NoInd, 1,
  	     gv->NoInd, 1,
  	     one_vec, gv->y, tmp);

  gv->b [0] = tmp[0];

  /* Init part 2*/

  double *tmp2 ;
  tmp2 = malloc ( gv->NoMar *sizeof(double) );
  if ( NULL == tmp2 ) ERR_M;

  gs_tmm_mul(gv->NoInd, gv->NoMar,
  	     gv->NoInd, 1,
  	     gv->Z, gv->y, tmp2);

  for ( j = 0; j < gv->NoMar; j++)  
    gv->b [ 1+j ] = tmp2[j];

  /* Store b */
  if (1) 
    {
      gv->bs =  realloc( gv->bs, dim_b * sizeof(double) );
      if ( NULL == gv->bs ) ERR_M;
      for ( j = 0; j < dim_b; j++)  
	gv->bs[j] = gv->b[j];
    }

  free(one_vec);
  free(tmp);
  free(tmp2);

  /* Write b vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( i = 0; i < dim_b ; i++)
	fprintf(outfp,"%f \n",gv->b [i] );
      
      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mme_rhs_01_GV ( char** out_filename,
			int*   auxfiles, 
			int*   retval,
			char** set_name    )
{
  gs_mme_rhs_01 ( gs_fdta(*set_name) ,
		  out_filename,
		  auxfiles, 
		  retval        ) ;
}

void gs_mme_restrhs_01 ( gs_varset_type *gv ,
			 char** out_filename,
			 int*   auxfiles, 
			 int*   retval        )
{
  gs_UNSIGNED i, dim_b ;
  FILE *outfp;

  int Srv; int* Sretval = &Srv;

  dim_b = 1 + gv->NoMar;

  if ( NULL == gv->bs ) 
    {
      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }
  else
    {
      gv->b =  realloc( gv->b, dim_b * sizeof(double) );
      if ( NULL == gv->b ) ERR_M;

      memcpy ( gv->b, gv->bs, (size_t) dim_b );
    }

  /* Write b vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( i = 0; i < dim_b ; i++)
	fprintf(outfp,"%f \n",gv->b [i] );
      
      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mme_restrhs_01_GV ( char** out_filename,
			    int*   auxfiles, 
			    int*   retval        )
{
  gs_mme_restrhs_01 ( gv ,
		     out_filename,
		     auxfiles, 
		     retval        ) ;
}


void gs_mme_solve_01 ( gs_varset_type *gv ,
		      char** out_filename,
		      int*   auxfiles, 
		      int*   retval        )
{
  gs_UNSIGNED i;
  FILE *outfp;

  int Srv=0; int* Sretval = &Srv;

  if ( NULL == gv->A ) {
      gs_info(gv,-2, "No coefficient matrix available");
      *retval = -2; return ;
  }

  if ( NULL == gv->b )  {
      gs_info(gv,-2, "No right hand side available");
      *retval = -2; return ;
  }

  /* Allocate memory for the solution vector */
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  gs_solve_Axb ( gv->dimA, gv->dimA, gv->A, gv->u, gv->b, Sretval );

  free(gv->A); 
  gv->A =NULL;

  if (*Sretval==-2) {*retval = -2; return ;}

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
   
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect\n");
      for ( i = 0; i < gv->dimu ; i++)
	fprintf(outfp,"%s.%i %f \n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u [i]);

      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mme_solve_01_GV ( char** out_filename ,
			  int*   auxfiles     , 
			  int*   retval       ,
			  char** set_name      )
{
  gs_mme_solve_01 ( gs_fdta(*set_name) ,
		    out_filename,
		    auxfiles, 
	            retval        ) ;
}

void gs_mme_invcoeff_01 ( gs_varset_type *gv ,
			  char** out_filename,
			  int*   auxfiles, 
			  int*   retval        )
{ 
  // Function is only used by the R inteface
  // No use of linear algebra libraries

  if ( NULL == gv->A ) {
      gs_info(gv,-2, "No A matrix available, run 'build.A' first");
      *retval = -2; return ;
  }

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  {
    double **a, d;
    gs_UNSIGNED *indx, i, j;

    indx=ivector(1,gv->dimA);
    a = convert_matrix( &(gv->A[0]), 1, gv->dimA, 1, gv->dimA);

    double *col = vector(1,gv->dimA);
    double **y  = matrix(1,gv->dimA,1,gv->dimA);

    ludcmp(a,gv->dimA,indx,&d);

    for( j=1 ; j <= gv->dimA ; j++) {
      for ( i=1 ; i <= gv->dimA ; i++ ) col[i] = 0.0 ;
      col[j] = 1.0 ;
      lubksb(a,gv->dimA,indx,col);
      for ( i=1 ; i <= gv->dimA ; i++ ) y[i][j] = col[i];
    }
    
    for( j=1 ; j <= gv->dimA ; j++) 
      for ( i=1 ; i <= gv->dimA ; i++ ) 
	{
	  idx_A = (i-1)*gv->dimA + (j-1) ;
	  gv->A [idx_A] = y[i][j] ;
	}

    free_convert_matrix (a,1);
    free_ivector (indx,1);
    free_vector  (col,1);
    free_matrix  (y,1,1);
  }

  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mme_invcoeff_01_GV ( char** out_filename,
			     int*   auxfiles, 
			     int*   retval        )
{
  gs_mme_invcoeff_01 ( gv ,
		       out_filename,
		       auxfiles, 
		       retval        ) ;
}


void gs_restrict_marker_data_03 ( gs_varset_type *gv ,
				  char**  ind_fname,
				  char**  mar_fname,
				  int*    NoAllMAX,
				  double* MaMisMAX,
				  double* ExHetMIN,
				  double* InMisMAX, 
				  int*   retval        )
{
  gs_UNSIGNED i,j,m,idx,i_n;

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  /* (0) Prepare output lists */
  FILE *fp;
  char line[gs_MAXSTR];

  /* Default: Output all individuals*/
  gs_UNSIGNED out_i[gv->NoInd] ;  
  for ( i = 0; i < gv->NoInd; i++) out_i[i] = 1;

  /* Individuals file available */
  if (0 != strcmp(*ind_fname,"")) {
    for ( i = 0; i < gv->NoInd; i++) out_i[i] = 0;
    if ( NULL == (fp = fopen(*ind_fname,"r")) ) ERR_F;
    while (!feof(fp)) 
      if ( 0 < fscanf(fp,"%s",line) ) {
	if ( strlen(line) > gs_MAXSTR) {
	  gs_info ( gv,-2,"Decription too long");
	  *retval = -2; return ;
	}
	for ( i=0; i < gv->NoInd; i++ ) 
	  if ( 0 == strcmp ( gv->IndName[i] , line ) ) {
	    out_i[i] = 1;
	    break;
	  }
      }
    fclose(fp);
  }

  /* Default: Output all markers */
  gs_UNSIGNED out_m[gv->NoMar] ;  
  for ( j = 0; j < gv->NoMar; j++) out_m[j] = 1;

  /* Individuals file available */
  if (0 != strcmp(*mar_fname,"")) {
    for ( j = 0; j < gv->NoMar; j++) out_m[j] = 0;
    if ( NULL == (fp = fopen(*mar_fname,"r")) ) ERR_F;
    while (!feof(fp)) 
      if ( 0 < fscanf(fp,"%s",line) ) {
	if ( strlen(line) > gs_MAXSTR) {
	  gs_info ( gv,-2,"Decription too long");
	  *retval = -2; return ;
	}
	for ( j=0; j < gv->NoMar; j++ ) 
	  if ( 0 == strcmp ( gv->MarName[j] , line ) ) {
	    out_m[j] = 1;
	    break;
	  }
      }
    fclose(fp);
  }


  /* (1) Determine number of markers and individuals  */

  gs_UNSIGNED CountInd = 0, CountMar = 0;

  for ( i=0; i < gv->NoInd; i++  )
    if ( ( *InMisMAX  >= gv->IndMis[i] ) && 
         (1 == out_i[i])                    )
	CountInd ++;

  for ( m=0; m < gv->NoMar; m++ )
    if ( ( ( (gs_UNSIGNED)(*NoAllMAX) ) >= 
	   (-1 == gv->All[m][0] ? gv->NoAll[m]-1 : gv->NoAll[m] )
	  ) &&
	 ( *MaMisMAX >= gv->Mis[m] ) &&
	 ( *ExHetMIN <= gv->Pic[m] ) && 
         ( 1 == out_m[m] )                                         )
      CountMar ++;

  int no_marker_reduction = (int)(CountMar == gv->NoMar);

  /* (2) Allocate memory for the new data set */

  // Marker matrix

  gs_UNSIGNED NoInd = CountInd; 
  gs_UNSIGNED NoMar = CountMar; 

  char (*IndName) [gs_MAXNME] = malloc(NoInd*sizeof(char[gs_MAXNME])); 
  char (*MarName) [gs_MAXNME] = malloc(NoMar*sizeof(char[gs_MAXNME])); 

  int * All1  =  malloc( NoMar * NoInd * sizeof(int) );
  int * All2  =  malloc( NoMar * NoInd * sizeof(int) );

  if ((NULL==All1)||(NULL==All2)||(NULL==IndName)||(NULL==MarName)) 
    ERR_M;

  // Linkage map 

  gs_UNSIGNED         t_NLoc  = 0  ;  
  gs_mappoint_type  * t_map_S = NULL;
  gs_mappoint_type ** t_map   = NULL;
  gs_UNSIGNED       * t_mdl   = NULL;

  if (NULL != gv->map_S ) // only if a map is loaded
    {
      t_NLoc  = NoMar  ;  
      t_map_S = malloc ( t_NLoc*sizeof(gs_mappoint_type) );
      t_map   = malloc ( t_NLoc*sizeof(gs_mappoint_type*))  ;   
      t_mdl   = malloc ( t_NLoc*sizeof(gs_UNSIGNED))   ; 
      if ((NULL==t_map_S)||(NULL== t_map)||(NULL==t_mdl))  ERR_M;
    }

  // Phenotype data

  double *y = NULL;

   if (NULL != gv->y ) // only if a phenotypes are loaded
    {
      y = malloc( NoInd * sizeof(double) );
      if (NULL==y)  ERR_M;
    }


  /* (3) Copy matching data to new storage */

  // Individuals names
  CountInd = 0;
  for (i=0; i<gv->NoInd; i++  )
    if ( ( *InMisMAX  >= gv->IndMis[i] ) && 
         (1 == out_i[i])                    )
      {
	strcpy ( IndName[CountInd] , gv->IndName[i]);
	CountInd ++;
      }

  // Marker names

  CountMar = 0;
  for ( m=0; m < gv->NoMar; m++ )
    if ( ( ( (gs_UNSIGNED)(*NoAllMAX) ) >= 
	   (-1 == gv->All[m][0] ? gv->NoAll[m]-1 : gv->NoAll[m] )
	 ) &&
	 ( *MaMisMAX >= gv->Mis[m] ) &&
	 ( *ExHetMIN <= gv->Pic[m] ) && 
         ( 1 == out_m[m] )                                        
        )
      {
	strcpy ( MarName[CountMar] , gv-> MarName[m] );
	CountMar ++;
      }

  // Marker data

  CountMar = 0;
  for ( m=0; m < gv->NoMar; m++ )
    if (
        ( ( (gs_UNSIGNED)(*NoAllMAX) ) >= 
	  (-1 == gv->All[m][0] ? gv->NoAll[m]-1 : gv->NoAll[m] )
	) &&
	( *MaMisMAX >= gv->Mis[m] ) &&
	( *ExHetMIN <= gv->Pic[m] ) && 
         ( 1 == out_m[m] )           
       )
      {
	CountInd = 0;
	for (i=0; i<gv->NoInd; i++  )
	  if ( ( *InMisMAX  >= gv->IndMis[i] ) && 
	       (1 == out_i[i])                    )
	    {
	      i_n = CountMar * NoInd     + CountInd;
	      idx =        m * gv->NoInd + i;
	      All1[i_n] = gv->All1[idx]; 
	      All2[i_n] = gv->All2[idx];
	      CountInd ++;
	    }
	CountMar ++;
      }

  // Marker map

  if (NULL != gv->map_S )
    {
      CountMar = 0;
      for ( m=0; m < gv->NoMar; m++ )
	if ( ( ( (gs_UNSIGNED)(*NoAllMAX) ) >= 
	       (-1 == gv->All[m][0] ? gv->NoAll[m]-1 : gv->NoAll[m] )
	       ) &&
	     ( *MaMisMAX >= gv->Mis[m] ) &&
	     ( *ExHetMIN <= gv->Pic[m] ) && 
	     ( 1 == out_m[m] )                                        
	     )
	  {
	    memcpy ( t_map_S + CountMar,  gv->map_S + m, 
		     (size_t)sizeof(gs_mappoint_type));
	    t_map   [CountMar] =  &t_map_S[CountMar];
	    t_mdl   [CountMar] =  CountMar;
	    CountMar ++;
	  }
      t_NLoc = CountMar;
    }

  // Phenotypes

  if ( NULL != gv->y )
    {
      CountInd = 0;
      for (i=0; i<gv->NoInd; i++  )
	if ( ( *InMisMAX  >= gv->IndMis[i] ) && 
	     (1 == out_i[i])                    )
	  {
	    y[CountInd] = gv->y[i];
	    CountInd ++;
	  }
    }

  /* If only indivdiuals are skipped, then copy marker effects*/
  
  gs_UNSIGNED  t_NoEff                = gv->NoEff  ;         
  char      (* t_EffNme) [gs_MAXNME]  = gv->EffNme ;
  int        * t_EffAll               = gv->EffAll ;  
  double     * t_u      	      = gv->u      ;    

  if ( no_marker_reduction )
    { // prevent the reset of gv from freeing the saved ptrs
      gs_info(gv,1,"Effects retained");

      gv->EffNme = NULL; 
      gv->EffAll = NULL; 
      gv->u      = NULL; 
    }
  else
    {
      gs_info(gv,1,"Effects removed");
    }

  /* (3) Replace old with new data matrix */

  gs_reset(gv);

  /* Marker data */
  gv->NoInd = NoInd; 
  gv->NoMar = NoMar; 

  gv->IndName = IndName;  
  gv->MarName = MarName; 

  gv->All1 = All1;
  gv->All2 = All2;

  /* Map  */
  gv->NLoc  = t_NLoc  ;
  gv->map_S = t_map_S ;
  gv->map   = t_map   ;
  gv->mdl   = t_mdl   ;
  
  /* Phenotypes */
  gv->y = y; 

  /* Effects */
  if ( no_marker_reduction )
    {
      gv->NoEff  = t_NoEff  ;         
      gv->EffNme = t_EffNme ;
      gv->EffAll = t_EffAll ;  
      gv->u      = t_u      ;    
    }

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv;
  char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return;
}

void gs_restrict_marker_data_03_GV ( char**  ind_fname,
				     char**  mar_fname,
				     int*    NoAllMAX,
				     double* MaMisMAX,
				     double* ExHetMIN,
				     double* InMisMAX, 
				     int*    retval,
				     char**  set_name       )
{
  gs_restrict_marker_data_03 ( gs_fdta(*set_name) ,
			       ind_fname,
			       mar_fname,
			       NoAllMAX,
			       MaMisMAX,
			       ExHetMIN,
			       InMisMAX, 
			       retval        );
}


void gs_copy_marker_data_01 ( gs_varset_type *new ,
			      gs_varset_type *gv  ,
			      int*   retval         )
{
  gs_UNSIGNED i,m;

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  gs_reset(new);

  /* Marker matrix */

  new->NoInd = gv->NoInd; 
  new->NoMar = gv->NoMar; 

  new->IndName = realloc( new->IndName, 
			  new->NoInd * sizeof(char[gs_MAXNME]) ); 
  new->MarName = realloc( new->MarName, 
			  new->NoMar * sizeof(char[gs_MAXNME]) ); 
  new->All1   =  realloc( new->All1,  
		          new->NoMar *  new->NoInd * sizeof(int) );
  new->All2   =  realloc( new->All2,  
			  new->NoMar *  new->NoInd * sizeof(int) );

  if ( (NULL == new->All1)  || (NULL == new->All2) || 
       ( NULL == new->IndName )  || ( NULL == new->MarName ) ) ERR_M;

  for ( i=0; i < gv->NoInd; i++  )
	strcpy ( new->IndName[i] , gv->IndName[i]);

  for ( m=0; m < gv->NoMar; m++ )
	strcpy ( new->MarName[m] , gv->MarName[m] );

  int * p1 ; int * p2 ;

  for ( m=0; m < gv->NoMar; m++ )
    {
	  p1 = &(gv ->All1[m* gv->NoInd]); 
	  p2 = &(new->All1[m*new->NoInd]); 
	  memcpy( p2 , p1, (size_t)new->NoInd * sizeof(int) );
	  p1 = &(gv ->All2[m* gv->NoInd]); 
	  p2 = &(new->All2[m*new->NoInd]); 
	  memcpy( p2 , p1, (size_t)new->NoInd * sizeof(int) );
    }
  
  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv;
  char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( new ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  if (NULL != gv->map)  // Map
    {
      new->NLoc  = gv->NLoc;
      new->map_S = realloc ( new->map_S, 
			     new->NLoc*sizeof(gs_mappoint_type) );
      new->map   = realloc ( new->map,   
			     new->NLoc*sizeof(gs_mappoint_type*) );
      new->mdl   = realloc ( new->mdl,   
			     new->NLoc*sizeof(gs_UNSIGNED) );
      if ( NULL==new->map || NULL==new->map_S || NULL==new->mdl ) ERR_M; 

      memcpy( new->map_S ,
	      gv->map_S  , 
	      (size_t)new->NLoc * sizeof(gs_mappoint_type ) );

      memcpy( new->map ,
	      gv->map  , 
	      (size_t)new->NLoc * sizeof(gs_mappoint_type*) );

      memcpy( new->mdl ,
	      gv->mdl  , 
	      (size_t)new->NLoc * sizeof(gs_UNSIGNED) );
    }

  if (NULL != gv->y)  // Phenotypes
    {
      new->y = realloc( new->y, new->NoInd*sizeof(double) );
      if ( NULL== new->y ) ERR_M;
      memcpy( new->y ,
	      gv->y  , 
	      (size_t)new->NoInd * sizeof(double) );
    }

  *retval = 0; return;
}

void gs_copy_marker_data_01_GV ( int*    retval     ,
				 char**  set_name_t ,
				 char**  set_name_s   )
{
  gs_copy_marker_data_01 ( gs_fdta(*set_name_t) ,
			   gs_fdta(*set_name_s) , 
			   retval        );
}


void gs_estimate_bv_01 ( gs_varset_type *effect_set ,
			 gs_varset_type *marker_set ,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval          )
{
   if ( NULL == effect_set->u ) {
    gs_info(effect_set,-2,"No estimated effects available");
    *retval = -2; return;
   }
 
   if ( ( 0 == marker_set->NoInd) || ( 0 == marker_set->NoMar) ) {
     gs_info(marker_set,-2,"No marker data available");
     *retval = -2; return ;
   }
   gs_UNSIGNED i, m, idx_all, idx_Z;

   // The genetic value
   marker_set->est = realloc ( marker_set->est, 
			       marker_set->NoInd*sizeof(double) );
   if (NULL == marker_set->est) ERR_M;


   int by_matrix_multiplication = 1;

   // Check whether the same marker set is used 
   if ( marker_set->NoMar != effect_set->NoMar)
     {
       by_matrix_multiplication = 0;                // Not the same no. of markers
     }
   else {
     for ( m = 0; m < marker_set->NoMar; m++) {
       if ( 0 != strcmp(marker_set->MarName[m], effect_set->MarName[m]) ) {
	   by_matrix_multiplication = 0;          // Not the same markers
	   break;
	 }
     }
   }

   if ( marker_set != effect_set) by_matrix_multiplication = 0;

   if ( 1 == by_matrix_multiplication )
     {
       gs_info(effect_set,1,"Marker sets are identical");
       
       // Make design matrix if necessary
       if ( NULL == marker_set->Z) {

	 gs_info(effect_set,1,"Constructing design matrix");

	 marker_set->Z =  realloc( marker_set->Z, marker_set->NoInd * 
				   marker_set->NoMar * sizeof(double) );
	 int a, b, lar_all, sml_all;
	 for ( m = 0; m < marker_set->NoMar; m++)
	   {
	     // Determine the larges allele for the marker 
	     lar_all = effect_set->EffAll[m+1];
	     sml_all = effect_set->CmpAll[m+1] ;
	     
	     // Construct design matrix
	     for ( i = 0; i < marker_set->NoInd ; i++)
	       {
		 idx_all = m*marker_set->NoInd + i ;
		 
		 a = marker_set->All1[idx_all] ;
		 b = marker_set->All2[idx_all] ;
		 
		 idx_Z   = i*marker_set->NoMar + m ;

		 if ( (a==lar_all) && (b==lar_all) ) { 
		   marker_set->Z[idx_Z] = marker_set->AA; 
		 }
		 else	{
		   if ( (a==sml_all) && (b==sml_all) )  {
		     marker_set->Z[idx_Z] = marker_set->aa;
		   }
		   else {
		     if ( ((a==lar_all) && (b==sml_all)) ||
			  ((a==sml_all) && (b==lar_all))   ) {
		       marker_set->Z[idx_Z] = marker_set->aA ;
		     }
		     else {
		       if ( (a==lar_all) || (b==lar_all) ) {
			 marker_set->Z[idx_Z] = marker_set->nA ;
		       }
		       else {
			 if ( (a==sml_all) || (b==sml_all) ) {
			   marker_set->Z[idx_Z] = marker_set->an ;
			 }
			 else {
			   marker_set->Z[idx_Z] = marker_set->nn ;
			 }
		       }
		     }
		   }
		 }
	       }
	   }  // Build design matrix
       } // two different data sets?

       gs_mm_mul( marker_set->NoInd, marker_set->NoMar,
		  marker_set->NoMar,1,
		  marker_set->Z,
		  &(effect_set->u[1]),
		  marker_set->est);		 

       for ( i = 0; i < marker_set->NoInd; i++) 
	 marker_set->est[i] += effect_set->u[0];

     }
   else  // Not by matrix multiplication
     {
       gs_info(effect_set,1,"Different marker sets");

        #pragma omp parallel
        {
	  gs_UNSIGNED i,j,m,idx_all;
	  int a, b, lar_all, sml_all;
	  double count;
	  
          #pragma omp for
	  for ( i = 0; i < marker_set->NoInd; i++) 
	    {
	      marker_set->est[i] = effect_set->u[0];
	      for (j = 1; j < effect_set->dimu; j++) 
		{
		  int found = 0;
		  for ( m = 0; m < marker_set->NoMar; m++) 
		    {
		      if ( 0 == strcmp (marker_set->MarName[m],
					effect_set->EffNme[j]) )
			{
			  lar_all = effect_set->EffAll[j];
			  sml_all = effect_set->CmpAll[j] ;

			  idx_all = m*marker_set->NoInd + i ;
			  a = marker_set->All1[idx_all] ;
			  b = marker_set->All2[idx_all] ;

			  if ( (a==lar_all) && (b==lar_all) ) { 
			    count = effect_set->AA; 
			  }
			  else	{
			    if ( (a==sml_all) && (b==sml_all) )  {
			      count = effect_set->aa;
			    }
			    else {
			      if ( ((a==lar_all) && (b==sml_all)) ||
				   ((a==sml_all) && (b==lar_all))   ) {
				count = effect_set->aA ;
			      }
			      else {
				if ( (a==lar_all) || (b==lar_all) ) {
				  count = effect_set->nA ;
				}
				else {
				  if ( (a==sml_all) || (b==sml_all) ) {
				    count = effect_set->an ;
				  }
				  else {
				    count = effect_set->nn ;
				  }
				}
			      }
			    }
			  }

			  marker_set->est[i] += count * effect_set->u[j];
			  found = 1;
			  break;
			}
		    }
		  if (0 == found) {
		    count = effect_set->nn;
		    marker_set->est[i] += count * effect_set->u[j];
		  }
		}
	    }
	}  // parallel
     } // else by matrix multiplication
     
 /* Write estimated values to output file */
 if (1 == *auxfiles)
   {
     gs_UNSIGNED i;
     FILE *fp;
     if ( (fp = fopen(*out_filename,"w")) == NULL ) ERR_F;
     if (NULL == marker_set->y )
       {
	 fprintf(fp,"yhat \n");
	 for ( i = 0; i < marker_set->NoInd ; i++)
	   fprintf(fp,"%s %f\n",
		   marker_set->IndName[i],
		   marker_set->est[i]);
       }
     else
       {
	 fprintf(fp,"y yhat \n");
	 for ( i = 0; i < marker_set->NoInd ; i++)
	   fprintf(fp,"%s %f %f\n",
		   marker_set->IndName[i],
		   marker_set->y[i],
		   marker_set->est[i]);
       }

     fclose(fp);  
   }

 *retval = 0; return;
}

void gs_estimate_bv_01_GV ( char** effect_set ,
			    char** marker_set ,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval          )
{
  gs_estimate_bv_01 ( gs_fdta(*effect_set),
		      gs_fdta(*marker_set),
		      out_filename,
		      auxfiles, 
		      retval          ) ;

}

/* Estimate breeding values for the loaded marker set */
void gs_check_model_fit_01_GV ( char**  out_filename ,
				int*    auxfiles     , 
				int*    retval       )
{
  gs_estimate_bv_01 ( gv,
		      gv,
		      out_filename,
		      auxfiles, 
		      retval              ) ;

}

/* Estimate breeding values for an new marker set */
void gs_estimate_breeding_values_01_GV ( char**  in_marker,
					 char**  out_filename,
					 int*    auxfiles, 
					 int*    retval    )
{
  int Srv; int* Sretval = &Srv;
 
  gs_varset_type VS; gs_varset_type *valset = &VS;
  gs_init (valset);

  gs_read_marker_data_01 (valset,in_marker,nixp,nixp,nixp,neinp,Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }
                                 
  gs_estimate_bv_01 (gv,valset,out_filename,auxfiles,Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_reset(valset);

  *retval = 0; return;

}

void gs_validation_01_GV ( char**  in_marker,
			   char**  in_performance,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval     )
{

  int Srv; int* Sretval = &Srv;
 
  gs_varset_type VS; gs_varset_type *valset = &VS;
  gs_init (valset);

  gs_read_marker_data_01 (valset,in_marker,nixp,nixp,nixp,neinp,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }

  gs_read_performance_data_01 (valset,in_performance,nixp,neinp,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }

  gs_estimate_bv_01 (gv,valset,nixp,neinp,Sretval);
  if (-2 == *Sretval) {gs_reset(valset); *retval = -2; return; }


 /* Write estimated values to output file */
  if (1 == *auxfiles)
    {
      FILE *fp;
      if ( NULL == (fp = fopen(*out_filename,"w"))  ) 
	{gs_reset(valset);ERR_F;}
      fprintf(fp,"y yhat \n");
      for (gs_UNSIGNED i = 0; i < valset->NoInd ; i++)
	fprintf(fp,"%s %f %f\n",
		valset->IndName[i],
		valset->y[i],
		valset->est[i]);
      fclose(fp);  
    }

  gs_reset (valset);
  
  *retval = 0; return;
}

void gs_build_esvs_01 ( gs_varset_type *gv ,
			gs_varset_type *es ,
			gs_varset_type *vs ,
			int*    n_estimation ,
			char**  out_es_m,
			char**  out_es_p,
			char**  out_vs_m,
			char**  out_vs_p,
			int*    auxfiles, 
			int*    retval          )
{
  
  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
     gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
     *retval = -2; return ;
   }

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
     *retval = -2; return ;
   }

  gs_reset(es);
  gs_reset(vs);

   if (!rgs) { srand( (unsigned) time(NULL) ); rgs=1; }

   gs_UNSIGNED CVSiEst  = (gs_UNSIGNED)*n_estimation;

   gs_UNSIGNED i,j,k,m, idx_g, idx_e, idx_v, idx;

   /* Estimation set */
   gs_UNSIGNED *CVes = malloc ( gv->NoInd * sizeof(gs_UNSIGNED) );
   if (NULL == CVes) ERR_M;
 
   /* Random permutation*/

   int *p = malloc( gv->NoInd * sizeof(int) );
   if (NULL == p) {free(CVes);ERR_M;}

   /* Random permutaion of the indiviudals */
   /* Knuth Stanford Graph Base */
   for (int i = 0; i < (int)gv->NoInd; ++i) {
     int j = rand() % (i + 1);
     p[i] = p[j];
     p[j] = i;
   }
   /* Estimation set 1, validation set 0 */
   /* Take the first random individuals as estimation set*/
   for ( i = 0; i < gv->NoInd; i++ ) CVes[i]      = 0;
   for ( i = 0; i < CVSiEst; i++ )   CVes[ p[i] ] = 1;
   free(p);

   /* Copy resampled data */

   /* Individuals */
   es->NoInd = CVSiEst;
   vs->NoInd = gv->NoInd - es->NoInd;

   es->IndName = realloc(es->IndName, es->NoInd * sizeof(char[gs_MAXNME]));
   vs->IndName = realloc(vs->IndName, vs->NoInd * sizeof(char[gs_MAXNME]));
   if ( (NULL==es->IndName) || (NULL==vs->IndName) ) {free(CVes);ERR_M;}

   for ( i=0,j=0,k=0; i < gv->NoInd; i++)
     {
       if ( 1 == CVes[i] )
	 {
	   strcpy( es->IndName[j] , gv->IndName[i] );
	   j++;
	 }
       else
	 {
	   strcpy( vs->IndName[k] , gv->IndName[i] );
	   k++;
	 }
     }

   /* Markers*/
   es->NoMar = gv->NoMar;
   vs->NoMar = gv->NoMar;

   es->MarName = realloc ( es->MarName, es->NoMar * sizeof(char[gs_MAXNME]) );
   vs->MarName = realloc ( vs->MarName, vs->NoMar * sizeof(char[gs_MAXNME]) );
   if ( (NULL==es->MarName) || (NULL==vs->MarName) ) {free(CVes);ERR_M;}
   
   for ( m = 0; m < gv->NoMar; m ++)
     {
	  strcpy( es->MarName[m] , gv->MarName[m] );
	  strcpy( vs->MarName[m] , gv->MarName[m] );
     }

  /* Genotypic data */
  es->All1 =  realloc(es->All1, es->NoMar * es->NoInd * sizeof(int));
  es->All2 =  realloc(es->All2, es->NoMar * es->NoInd * sizeof(int));
  if ( (NULL == es->All1)  || (NULL == es->All2) ) {free(CVes);ERR_M;}

  vs->All1 =  realloc(vs->All1, vs->NoMar * vs->NoInd * sizeof(int));
  vs->All2 =  realloc(vs->All2, vs->NoMar * vs->NoInd * sizeof(int));
  if ( (NULL == vs->All1)  || (NULL == vs->All2) ) {free(CVes);ERR_M;}

  for ( m = 0; m < gv->NoMar; m++ ) 
    {
      for ( i=0,j=0,k=0; i < gv->NoInd; i++)
	{
	   idx_g = m*gv->NoInd + i;
	   idx_e = m*es->NoInd + j;
	   idx_v = m*vs->NoInd + k;

	   if ( 1 == CVes[i] )
	     {
	       es->All1[idx_e] = gv->All1[idx_g];
	       es->All2[idx_e] = gv->All2[idx_g];
	       j++;
	     }
	   else
	     {
	       vs->All1[idx_v] = gv->All1[idx_g];
	       vs->All2[idx_v] = gv->All2[idx_g];
	       k++;
	     }
 	 }
     }

   /* Note: further data as generated by read_marker_data may be required*/

   /* Performance data */
   es->y = realloc( es->y, es->NoInd*sizeof(double) );
   vs->y = realloc( vs->y, vs->NoInd*sizeof(double) );
   if ( ( NULL== es->y ) || ( NULL== vs->y ) ) {free(CVes);ERR_M;}

   for ( i=0,j=0,k=0; i < gv->NoInd; i++)
     {
       if ( 1 == CVes[i] )
	 { 
	   es->y[j] = gv->y[i];
	   j++;
	 }
       else
	 {
	   vs->y[k] = gv->y[i];
	   k++;
	 }
     }

   free(CVes);

  /* Write output files */
  if (1 == *auxfiles){

    FILE * fp1 = fopen(*out_es_m,"w") ;
    FILE * fp2 = fopen(*out_es_p,"w") ;
    FILE * fp3 = fopen(*out_vs_m,"w") ;
    FILE * fp4 = fopen(*out_vs_p,"w") ;
    if ( (NULL==fp1) || (NULL==fp2) || (NULL==fp3) || (NULL==fp4)) ERR_F;

      fprintf(fp1, "I;M;A\n");
      for ( i=0; i < es->NoInd; i++  )
	for ( m=0; m < es->NoMar; m++ ) 
	  {
	    idx = m*es->NoInd + i;
	    if  (-1 != es->All1[idx]) 
	      {
		fprintf(fp1, "%s;%s;%i\n",
			es->IndName[i], es->MarName[m], es->All1[idx] );
		if (es->All1[idx] != es->All2[idx])
		  fprintf(fp1, "%s;%s;%i\n",
			  es->IndName[i], es->MarName[m], es->All2[idx] );
	      }
	  }

      fprintf(fp2, "I;P\n");
      for ( i=0; i < es->NoInd; i++  )
	fprintf(fp2, "%s;%f\n", es->IndName[i], es->y[i] );

      fprintf(fp3, "I;M;A\n");
      for ( i=0; i < vs->NoInd; i++  )
	for ( m=0; m < vs->NoMar; m++ ) 
	  {
	    idx = m*vs->NoInd + i;
	    if  (-1 != vs->All1[idx]) 
	      {
		fprintf(fp3, "%s;%s;%i\n",
			vs->IndName[i], vs->MarName[m], vs->All1[idx] );
		if (vs->All1[idx] != vs->All2[idx])
		  fprintf(fp3, "%s;%s;%i\n",
			  vs->IndName[i], vs->MarName[m], vs->All2[idx] );
	      }
	  }

      fprintf(fp4, "I;P\n");
      for ( i=0; i < vs->NoInd; i++  )
	fprintf(fp4, "%s;%f\n", vs->IndName[i], vs->y[i] );

      fclose(fp1);  
      fclose(fp2);  
      fclose(fp3);  
      fclose(fp4);  

  }
  
  *retval = -0; return ;

}

void gs_build_esvs_01_GV ( int*    n_estimation ,
			   char**  base_set,
			   char**  estimation_set,
			   char**  validation_set,
			   char**  out_es_m,
			   char**  out_es_p,
			   char**  out_vs_m,
			   char**  out_vs_p,
			   int*    auxfiles, 
			   int*    retval          )
{

  /* Achtung instabil bei gleichen Ausgabefilenamen */

  /* First use of global validation and estimation set */
  if (0==esi) {gs_init(es);esi=1;}
  if (0==vsi) {gs_init(vs);vsi=1;}

  gs_varset_type* ees = es;
  gs_varset_type* vvs = vs;

  if ( 0 != strcmp ( "none" , *estimation_set ) )
    ees = gs_fdta(*estimation_set);

  if ( 0 != strcmp ( "none" , *validation_set ) )
    vvs = gs_fdta(*validation_set);

  gs_build_esvs_01 ( gs_fdta(*base_set) ,
		  ees ,
		  vvs ,
		  n_estimation ,
		  out_es_m,
		  out_es_p,
		  out_vs_m,
		  out_vs_p,
		  auxfiles, 
		  retval          );
}

void gs_esteff_rrwr_01 (  gs_varset_type *gv ,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_reg_01 ( gv, hsq, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwr_01_GV ( double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_esteff_rrwr_01 ( gs_fdta(*set_name),
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_esteff_rrwa_01 (  gs_varset_type *gv ,
			  double* hsq,
			  double* alpha,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_aov_01 ( gv, hsq, alpha, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwa_01_GV ( double* hsq,
                            double* alpha,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_esteff_rrwa_01 ( gs_fdta(*set_name),
		      hsq,
		      alpha,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_esteff_rrwe_01 (  gs_varset_type *gv ,
			  double* lambda, 
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_const_01 ( gv, lambda,nixp, neinp, retval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwe_01_GV ( double* lambda   , 
			    char**  out_filename,
			    int*    auxfiles , 
			    int*    retval   ,
			    char** set_name   )
{
  gs_esteff_rrwe_01 ( gs_fdta(*set_name),
		      lambda            ,
		      out_filename      ,
		      auxfiles          , 
		      retval             ) ;
}

void gs_esteff_rrwe_02 (  gs_varset_type *gv ,
		       double* hsq,
		       char**  out_filename,
		       int*    auxfiles, 
		       int*    retval         )
{
  int Srv; int* Sretval = &Srv;
	   
  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_const_02 ( gv, hsq, nixp, neinp, retval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rrwe_02_GV ( double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_esteff_rrwe_02 ( gs_fdta(*set_name),
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_return_effects_01 (  gs_varset_type *gv ,
			     char**  out_filename,
			     int*    retval         )
{
  gs_UNSIGNED i;
  FILE *outfp;
   
  if ( NULL == gv->u ) {
      gs_info(gv,-2, "No effects estimated");
      *retval = -2; return ;
  }

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;
  
  fprintf(outfp,"marker allele cmpallele effect\n");
  for ( i = 0; i < gv->dimu ; i++)
    fprintf(outfp,"%s %i %i %20.16f \n",
	    gv->EffNme[i],
	    gv->EffAll[i],
	    gv->CmpAll[i],
	    gv->u [i]);
  
  fclose(outfp);  

  *retval = 0; 
  return ;
}

void gs_return_effects_01_GV ( char**  out_filename,
			       int*    retval,
			       char** set_name         )
{
  gs_return_effects_01 ( gs_fdta(*set_name),
			 out_filename,
			 retval         ) ;
}


void gs_set_effects_01 (  gs_varset_type *gv ,
		          char**  filename,
		          int*    retval         )
{

  gs_UNSIGNED i,j,k;

  /* Check wheter marker data is loaded */
  if ( 0 == gv->NoMar ) {
      gs_info(gv,-2, "No marker data available");
      *retval = -2; return ;
  }

  /* Remove old design matrix and estimated genotypes */
  gv->Z   = realloc (gv->Z,0);   gv->Z   = NULL; 
  gv->est = realloc (gv->est,0); gv->est = NULL;

  /* Allocate memory for the effect descriptions */
  gv->NoEff =  1 + gv->NoMar;
  gv->EffNme = realloc (gv->EffNme, gv->NoEff * sizeof(char[gs_MAXNME]) );
  gv->EffAll = realloc (gv->EffAll, gv->NoEff * sizeof(int) );
  gv->CmpAll = realloc (gv->CmpAll, gv->NoEff * sizeof(int) );
  if ( (NULL==gv->EffNme) || (NULL==gv->EffAll) || (NULL==gv->CmpAll)) ERR_M;

  /* Names of effects */
  strcpy(gv->EffNme[0],"mu");
  for ( i = 1; i < gv->NoEff; i++)
    strcpy(gv->EffNme[i],gv->MarName[i-1]);		       

  /* Allocate memory for the solution vector */
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  /* Initialize effects with zeros */
  for ( i = 0; i < gv->dimu ; i++) {
    gv->EffAll[i] = 0;
    gv->u [i]     = 0;

  }

  {
    gs_UNSIGNED NoLines = 0; /* Number of lines in the input file       */

    FILE *fp;

    char*  line = NULL;   char** linep = &line;
    size_t n    = 0;      size_t* np   = &n;
    
    char *p ;

    if ( NULL == (fp = fopen(*filename,"r")) ) ERR_F;

    /* Headline */
    getline ( linep , np , fp ) ;

    /* Count lines */
    while (!feof(fp))
      {
	if ( 0 < getline ( linep , np , fp ) )
	  { 
            // Read in a data line

	    p = strtok (line,gs_tok_del);
	    if ( NULL == p ) 
	      {gs_info(gv,-2,"Error in map file");*retval=-2;return;}
	    if ( strlen(p) > gs_MAXNME)
	      {gs_info(gv,-2,"Marker name to long");*retval=-2;return;}
	    
	    p = strtok ('\0',gs_tok_del);
	    if  ( NULL == p )  
	      {gs_info(gv,-2,"Error in map file");*retval=-2;return;}
	    
	    p = strtok ('\0',gs_tok_del); 
	    if ( NULL == p ) {gs_info(gv,-2,"Error in map file");*retval=-2;return;}

	    p = strtok ('\0',gs_tok_del); 
	    if ( NULL == p ) {gs_info(gv,-2,"Error in map file");*retval=-2;return;}

	    strcpy(line,"\0");
	    NoLines++; 
	  }
      }
    fclose(fp);

    /* Allocate memory */
    char   (*EffNme) [1+gs_MAXNME] = NULL; /* Effect name */
    int     *EffAll                = NULL;  /* Effect allele */     
    int     *CmpAll                = NULL;  /* Complementary allele */     
    double  *Eff                   = NULL;  /* Effect size */ 
  
    EffNme = realloc ( EffNme, NoLines*sizeof(char[1+gs_MAXNME]) );
    EffAll = realloc ( EffAll, NoLines*sizeof(int) );
    CmpAll = realloc ( CmpAll, NoLines*sizeof(int) );
    Eff    = realloc ( Eff ,   NoLines*sizeof(double) );
    if ( (NULL == EffNme) || (NULL == EffAll) || (NULL == CmpAll)
	 || (NULL == Eff) ) ERR_M;

    /* Read in data lines */
    if ( NULL == (fp = fopen(*filename,"r")) ) ERR_F;
    getline ( linep , np , fp ) ;
    for (j = 0; j < NoLines; j++)
      {
	getline ( linep , np , fp ) ;
	gs_remove_newline ( line ) ;

	p = strtok (line,gs_tok_del);
        strcpy(EffNme[j],p);

	p = strtok ('\0',gs_tok_del);
	EffAll[j] = (int)atoi (p);

	p = strtok ('\0',gs_tok_del);
	CmpAll[j] = (int)atoi (p);

	p = strtok ('\0',gs_tok_del);
	Eff[j] = (double)atof(p);

      }
    fclose(fp);
    free(line);

    for (j = 0; j < NoLines; j++)
      {
	if ( 0 == strcmp(EffNme[j],gv->EffNme[j]) ) 
	  { /* same order of effect names */
	    gv->EffAll[j] = EffAll[j];
	    gv->CmpAll[j] = CmpAll[j];
	    gv->u[j]      = Eff[j] ;
	    
	  }
	else
	  { /* different order of effect names */
	    for (k = 0; k < (1 + gv->NoMar) ; k++)
	      {
		if ( 0 == strcmp(EffNme[j],gv->EffNme[k]) ) {
		  gv->EffAll[k] = EffAll[j];
		  gv->CmpAll[k] = CmpAll[j];
		  gv->u[k]      = Eff[j] ;
		  break;
		}
	      }
	  }
      }

    free(EffNme);
    free(EffAll);
    free(CmpAll);
    free(Eff   );

  }

  *retval = 0; 
  return ;
}

void gs_set_effects_01_GV ( char**  filename,
		            int*    retval,
		            char** set_name         )
{
  gs_set_effects_01 ( gs_fdta(*set_name),
		      filename,
		      retval         ) ;
}




void gs_return_pvals_01 (  gs_varset_type *gv ,
			   char**  out_filename,
			   int*    retval         )
{
  gs_UNSIGNED i;
  FILE *outfp;
   
  if ( NULL == gv->p ) {
      gs_info(gv,-2, "No p-values estimated");
      *retval = -2; return ;
  }

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;
  
  fprintf(outfp," effect pvalue \n");
  for ( i = 0; i < gv->dimu ; i++)
    fprintf(outfp,"%s.%i %f %f \n",
	    gv->EffNme[i],
	    gv->EffAll[i],
	    gv->u [i],
	    gv->p [i]);
  
  fclose(outfp);  

  *retval = 0; 
  return ;
}

void gs_return_pvals_01_GV ( char**  out_filename,
			     int*    retval,
			     char** set_name         )
{
  gs_return_pvals_01 ( gs_fdta(*set_name),
			 out_filename,
			 retval         ) ;
}

void gs_check_pvals_01 (  gs_varset_type *gv ,
			  int*    retval         )
{
  if ( NULL == gv->p ) *retval = 0; 
  else *retval = 1;
  return;
}

void gs_check_pvals_01_GV (  int*    retval,
			     char** set_name         )
{
  gs_check_pvals_01 ( gs_fdta(*set_name),
		      retval         ) ;
}


void gs_set_allele_codes_01 (  gs_varset_type *gv ,
			       double* aa         ,
			       double* aA         ,
			       double* AA         ,
			       double* an         ,
			       double* nn         ,
			       double* nA           )
{
  gv->aa         = * aa ; 
  gv->aA         = * aA ;
  gv->AA         = * AA ;
  gv->an         = * an ; 
  gv->nn         = * nn ;
  gv->nA         = * nA ;

  /* Remove old design matrix and estimated genotypes */
  gv->Z   = realloc (gv->Z,0);   gv->Z   = NULL; 
  gv->est = realloc (gv->est,0); gv->est = NULL;

}

void gs_set_allele_codes_01_GV (  double* aa         ,
				  double* aA         ,
				  double* AA         ,
				  double* an         ,
				  double* nn         ,
				  double* nA         ,
				  char** set_name  )
{

  gs_set_allele_codes_01 ( gs_fdta(*set_name) ,
			   aa                 ,
			   aA                 ,
			   AA                 ,
			   an                 ,
			   nn                 ,
			   nA                  ) ;
}




void gs_lambda_emstep_01 ( gs_varset_type *gv     ,
			   int*    constvar       ,
			   char**  out_filename   ,
			   int*    auxfiles       , 
			   int*    retval         )
{
  gs_UNSIGNED i,j,m;

  if (  NULL== gv->ladi ) {
    gs_info(gv,-2,"No lambda defined");
    *retval = -2; return ;
  }

  int Srv; int* Sretval = &Srv;

  gs_mme_restcoeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_restrhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* --- */

  /* Allocate memory for the solution vector */
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  /* Initialize the solution vector */
  for ( i = 0; i < gv->dimu ; i++) gv->u [i] =  gv->b[i];

  /* Trace of the inverse coefficient matrix */
  gv->t0 = realloc(gv->t0, gv->dimA * sizeof(double) );

  { // LU decomposition 
    double d;

    gv->u0  = realloc ( gv->u0, gv->dimA * sizeof(gs_UNSIGNED) );

    /* LU decomposition */
     gs_ludcmp ( gv->A, gv->dimA, gv->u0, &d);

    /* Solutions for u */
     gs_lubksb ( gv->A, gv->dimA, gv->u0, gv->u );

    /* Inverse */
     gv->t4 = realloc (  gv->t4, gv->dimA*gv->dimA * sizeof(double) );

    gs_luinv ( gv->A, gv->dimA, gv->u0 ,gv->t4 );

    /* Diagonal elements of the inverse coefficient matrix */
     for ( j=0 ; j < gv->dimA ; j++ ) 
       gv->t0[j] = gv->t4[ gv->dimA *j + j];
   
  }


  /* Allocate memory for the variace components */
  gv->sigma =  realloc( gv->sigma, gv->dimA * sizeof(double) );
  if ( NULL == gv->sigma ) ERR_M;

  gv->t1 = realloc (gv->t1, sizeof(double));
  gv->t2 = realloc (gv->t2, ( (gv->NoInd > gv->NoMar) ?
		      gv->NoInd : gv->NoMar    ) 
		    * sizeof(double)            );
  gv->t3 = realloc (gv->t3, ( (gv->NoInd > gv->NoMar) ?
		           gv->NoInd : gv->NoMar    ) 
			 * sizeof(double)            );

  for ( i=0 ; i<gv->NoInd ; i++) gv->t2[i] = 1;

  gs_tmm_mul (gv->NoInd,1 ,         // y'y
	      gv->NoInd,1 ,	     
	      gv->y       ,
	      gv->y       ,
	      gv->sigma ) ;

  gs_tmm_mul (gv->NoInd,1 ,         // 1'y
	      gv->NoInd,1 ,	     
	      gv->t2          ,
	      gv->y       ,
	      gv->t1          ) ;

  gs_tmm_mul (gv->NoInd, gv->NoMar, // Z'y	  	     
	      gv->NoInd, 1 ,
	      gv->Z        ,
	      gv->y        ,
	      gv->t2           ) ;

  gs_tmm_mul (gv->NoMar, 1 ,	     
	      gv->NoMar, 1 ,
	      1+(gv->u)    ,
	      gv->t2           ,
	      gv->t3           ) ;
  
  if (0 == *constvar) /* Variable variance of the markers     */
    { /* Residual is held constant, and needs to be estimated */
      /* in a previous run with constant variances            */
       gv->sigma[0] = gv->sigma_e ;
       for ( i=1; i < gv->dimA; i++)
	{
	  gv->sigma[i] = 
	     gv->u[i] * gv->u[i] + gv->sigma[0] * gv->t0[i]  ;
	}
    }
  else /* Constant variance across all markers              */
    {  /* Residual error is estimated in each iteration     */
      gv->sigma[0] -= gv->u[0] *gv->t1[0];
      gv->sigma[0] -= gv->t3[0];                     
      gv->sigma[0]  = gv->sigma[0] / (double)gv->NoInd;
      gv->sigma_e   = gv->sigma[0];

      gs_tmm_mul (gv->NoMar,1, gv->NoMar,1, 1+(gv->u), 1+(gv->u), gv->t2 ) ;

      double tr = 0; 
      for ( i=1; i < gv->dimA; i++ ) tr +=gv->t0 [ i ];

      for (i=1; i <gv->dimA; i++)
	gv->sigma[i] = (gv->t2[0] + gv->sigma[0] * tr ) / (double) gv->NoMar;

    }


  for (i=1; i <gv->dimA; i++)
    gv->ladi[i-1] = gv->sigma[0] / gv->sigma[i];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda var \n");
      fprintf(outfp,"err %f %f \n", 0.0, gv->sigma[0] );
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f\n",
		gv->MarName[m],
		gv->ladi[m] ,
		gv->sigma[1+m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_emstep_01_GV ( int*    constvar      ,
			      char**  out_filename  ,
			      int*    auxfiles      , 
			      int*    retval        ,
			      char** set_name       )

{
  gs_lambda_emstep_01 (  gs_fdta(*set_name),
			 constvar,
			 out_filename,
			 auxfiles, 
			 retval         ) ;
}

void gs_lambda_rml_01 (  gs_varset_type* gv ,
			 int*    maxiter,
			 double* precision,
			 int*    constvar,
			 int*    init_lambda,
			 double* hsq,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  gs_UNSIGNED m;
  int i;
  double qu, maxqu=0;
  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  double *vc_save = malloc (1+gv->NoMar*sizeof(double) );
  if (NULL==vc_save) ERR_M;

  for ( m = 0 ; m < 1+gv->NoMar ; m++) vc_save[m] = 0;

  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;


  if ( 1 == *init_lambda )
    {
      // Initialisation with 1 
      // for ( m = 0 ; m < gv->NoMar ; m++)  
      // gv->ladi[m] = 1 ;
       
      /* Initialization with constant lambda */
      double lambda =  ( 1 / *hsq -1 ) * (double)gv->NoMar; 
      for ( m = 0 ; m < gv->NoMar ; m++)  
	gv->ladi[m] = lambda;
    }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ; 
  if (-2 == *Sretval) { *retval = -2; return; }  

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }
      
  for (i=0; i < *maxiter; i ++) 
    {
      if (i>0)  for ( m = 0 ; m < 1+gv->NoMar ; m++) 
		  vc_save[m] = gv->sigma[m];

      gs_lambda_emstep_01 ( gv,constvar,nixp,neinp,Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      maxqu = 0;
      for ( m = 0 ; m < 1+gv->NoMar ; m++) {
	qu = fabs ( (vc_save[m] - gv->sigma[m]) / gv->sigma[m] );
	if (qu>maxqu) maxqu=qu;
      }
      if (maxqu < *precision) break; 
  }
  
  sprintf(gv->msg, "Number of EM iterations: %i, precision: %f",i,maxqu);
  gs_info(gv,1,gv->msg);
  sprintf(gv->msg, "(Greatest change of a var. comp. in the last iteration: %6.3f%%)",
	  maxqu*100);
  gs_info(gv,1,gv->msg);
  
  free(vc_save);

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmlc_01 (  gs_varset_type *gv ,
			  int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int _cc = 1; int* constvar    = &_cc;
  int _il = 1; int* init_lambda = &_il;

  gs_lambda_rml_01  ( gv,
		      maxiter,
		      precision,
		      constvar,
		      init_lambda,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_rmlv_01 (  gs_varset_type *gv ,
			  int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{

  int Srv =0; int* Sretval = &Srv;

  int _cc = 1; int* constvar = &_cc;
  int _il = 1; int* init_lambda = &_il;

  /* First run with constant variances to estimate residual */
  /* variancs                                              */
  _cc = 1 ;
  _il = 1;
  gs_lambda_rml_01  ( gv,
		      maxiter,
		      precision,
		      constvar,
		      init_lambda,
		      hsq,
		      nixp,
		      neinp, 
		      Sretval         ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* Second run with variable variances to estimate marker */
  /* variances                                             */
  _cc = 0; 
  _il = 1;
  gs_lambda_rml_01  ( gv,
		      maxiter,
		      precision,
		      constvar,
		      init_lambda,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_lambda_rmlv_01_GV ( int*    maxiter,
			    double* precision,
			    double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval         ,
			    char** set_name         )
{
  gs_lambda_rmlv_01 ( gs_fdta(*set_name),
		      maxiter,
		      precision,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_rmlc_01_GV ( int*    maxiter,
			    double* precision,
			    double* hsq,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval         ,
			    char** set_name         )
{
  gs_lambda_rmlc_01 ( gs_fdta(*set_name),
		      maxiter,
		      precision,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_lambda_rmlr_01 (  gs_varset_type *gv     ,
			  int*    maxiter        ,
			  double* precision      ,
			  double* hsq            ,
			  char**  out_filename   ,
			  int*    auxfiles       , 
			  int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_reg_01 ( gv, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlc_01 ( gv,
		      maxiter,precision,hsq,
		      out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  long double SSQM = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    SSQM += gv->SMsqm[m] ;

  SSQM = SSQM / (long double) gv->NoMar;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  ( SSQM / gv->SMsqm[m]  ) * gv->ladi[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmlr_01_GV ( int*    maxiter       ,
			    double* precision     ,
			    double* hsq            ,
			    char**  out_filename  ,
			    int*    auxfiles      , 
			    int*    retval        ,
			    char** set_name        )
{
  gs_lambda_rmlr_01 ( gs_fdta(*set_name),
		      maxiter,
		      precision,
		      hsq,
		      out_filename      ,
		      auxfiles          , 
		      retval            ) ;
}

void gs_lambda_rmla_01 (  gs_varset_type *gv     ,
			  double* alpha          ,
			  int*    maxiter        ,
			  double* precision      ,
			  double* hsq            ,
			  char**  out_filename   ,
			  int*    auxfiles       , 
			  int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlc_01 ( gv,
		      maxiter,precision,hsq,
		      out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  long double Svar = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  Svar = Svar / (long double) gv->NoMar;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  ( Svar / gv->SMvar[m]  ) * gv->ladi[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmla_01_GV ( double* alpha,
                            int*    maxiter       ,
			    double* precision     ,
			    double* hsq           ,
			    char**  out_filename  ,
			    int*    auxfiles      , 
			    int*    retval        ,
			    char** set_name        )
{
  gs_lambda_rmla_01 ( gs_fdta(*set_name),
		      alpha,
		      maxiter,
		      precision,
		      hsq,
		      out_filename      ,
		      auxfiles          , 
		      retval            ) ;
}


void gs_esteff_rmlc_01 ( gs_varset_type*, int*, double*,double*,
			 char**, int*, int*);



void gs_bayes_iterate_01 ( gs_varset_type* gv ,
			   int*    gibbsburn,
			   int*    gibbsiter,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval         )
{
  gs_UNSIGNED i,j,k,r,m,idx_Z;

  gs_FLOAT t[4];

  /* Temporary storage */
  gs_FLOAT *t1 = malloc ( gv->NoInd *  sizeof(gs_FLOAT) ) ;
  if (NULL==t1) ERR_M;

  gs_FLOAT *t2 = malloc ( gv->NoMar *  sizeof(gs_FLOAT) ) ;
  if (NULL==t2) ERR_M;

  gs_FLOAT *Zu = malloc ( gv->NoInd *  sizeof(gs_FLOAT) ) ;
  if (NULL==Zu) ERR_M;

  gs_FLOAT *z = malloc  ( gv->NoInd *  sizeof(gs_FLOAT) ) ;
  if (NULL==z) ERR_M;

  gs_FLOAT *gtmp = malloc ( gv->dimA *  sizeof(gs_FLOAT) ) ;
  if (NULL==gtmp) ERR_M;

  gs_FLOAT *gold = malloc ( gv->dimA *  sizeof(gs_FLOAT) ) ;
  if (NULL==gold) ERR_M;

  /* Vector of residuals */
  gs_FLOAT *e = malloc ( gv->NoInd *sizeof(gs_FLOAT) );
  if (NULL==e) ERR_M;

  /* Vector of ones */
  gs_FLOAT *o = malloc ( gv->NoInd *sizeof(gs_FLOAT) );
  if (NULL==o) ERR_M;
  for (i = 0; i < gv->NoInd; i++ ) o[i] = 1;

  /* Vector for storing the effects */
  gs_FLOAT *effects = malloc ( gv->dimA  *sizeof(gs_FLOAT) );
  if (NULL==effects) ERR_M;
  for (j = 0; j < gv->dimA; j++ ) effects[j] = 0;

  /* Z u and sigma as matrix in gs_FLOAT */
  gs_FLOAT *Z = malloc(gv->NoInd * gv->NoMar * sizeof(gs_FLOAT) );
  if (NULL==Z) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar; m++) {
      idx_Z   = i*gv->NoMar + m ;
      Z[idx_Z] = gv->Z[idx_Z];
    }
  gs_FLOAT *u =  malloc( gv->dimu * sizeof(gs_FLOAT) );
  if ( NULL == u ) ERR_M;
  for ( i = 0; i < gv->dimu ; i++) u [i] = gv->u [i];

  gs_FLOAT *sigma =  malloc( gv->dimA * sizeof(gs_FLOAT) );
  if ( NULL == sigma ) ERR_M;
  for ( i = 0; i < gv->dimA ; i++)  sigma[i] =  gv->sigma[i] ;

  gs_FLOAT *y =  malloc( gv->NoInd * sizeof(gs_FLOAT) );
  if ( NULL == y ) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++)  y[i] =  gv->y[i] ;

  
  // for (j = 0; j < gv->dimA; j++ )  u[j] = 0.0 , sigma[j] = 0.000001 ;  
  

  GetRNGstate();

  for (r=0; r < (gs_UNSIGNED) ( *gibbsburn + *gibbsiter ); r++) 
  {

    if ( round ( (double) r/100) == (double) r/100 ){
      sprintf ( gv->msg,
		"%6i/%i  se2: %5.2f,  b0: %5.2f",
		(int) r,
		(int) ( *gibbsburn + *gibbsiter ),
		(double)sigma[0],(double)u[0]  ) ;

      gs_info(gv,1,gv->msg);
    }

  /* Error variance */

  gs_mm_mul_gF ( gv->NoInd, gv->NoMar,            
		 gv->NoMar, 1        ,
		 Z,
		 &(u[1]),
		 Zu                    ); 

 for (i = 0; i < gv->NoInd; i++ )      
      e[i] = y[i]  -  Zu[i] - u[0];

  gs_tmm_mul_gF ( gv->NoInd, 1,        
		  gv->NoInd, 1, 
		  e,
		  e,
		  t1            );

  //gs_FLOAT a = t1[0];
  //gs_FLOAT b = rchisq( (double) gv->NoInd ) ;
  //sigma[0]  = a / b ;                            

 

  /* Mean */

  gs_mm_mul_gF ( gv->NoInd, gv->NoMar,            
		 gv->NoMar, 1        ,
		 Z,
		 &(u[1]),
		 Zu                    ); 

  gs_tmm_mul_gF (gv->NoInd, 1,
		 gv->NoInd, 1,
		 o,
		 Zu,
		 t1);

  gs_FLOAT sum_y = 0;
  for (i = 0; i < gv->NoInd; i++ ) sum_y += y[i];

  gs_FLOAT expmean =  ( sum_y - t1[0] ) / (gs_FLOAT) gv->NoInd;
  gs_FLOAT sdmean  =  sqrt( sigma[0]    / (gs_FLOAT) gv->NoInd ) ;
  u[0] = rnorm (expmean,sdmean);        
 
  /* Marker variances */

  for (j = 1; j < gv->dimA; j++ ) 
    {
      sigma[j] =  ( u[j] * u[j] ) / rchisq ( 1+4.012 ); // <----
    }
  
  
  /* Effects */

  for (j=1; j < gv->dimA; j++) gold[j] = u[j];

  for (j=1; j < gv->dimA; j++)
    {
      for (k=1; k < gv->dimA; k++) gtmp[k] = gold[k];
      gtmp[j] = 0 ; 

      for (i=0; i < gv->NoInd; i ++) {
	z[i] = Z[ (i*gv->NoMar) + (j-1) ];
      }

      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,1,
		     z,
		     y,
		     &(t[0]) );
      
      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,gv->NoMar,
		     z,
		     Z,
		     t2 );
      
      gs_mm_mul_gF( 1,gv->NoMar,
		    gv->NoMar,1,
		    t2,
		    &(gtmp[1]),
		    &(t[1]) );

      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,1,
		     z,
		     o,
		     &(t[2]) );

      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,1,
		     z,
		     z,
		     &(t[3]) );

      gs_FLOAT expeff = ( t[0] - t[1] - t[2]*u[0] ) /
	              ( t[3] + sigma[0] / sigma[j]);

      gs_FLOAT sdeff  = sqrt ( sigma[0] /
			   ( t[3] + sigma[0] / sigma[j]));

      u[j] = rnorm(expeff,sdeff);                      // <----

    }


  /* Store the effects */
  if ( r >=  (gs_UNSIGNED) *gibbsburn )
    for (j = 0; j < gv->dimA; j++ ) 
      effects[j]  += u[j];

  }

  /* Write back average effects */
  for (j = 0; j < gv->dimA; j++ ) {
    gv->u[j] = effects[j] / (gs_FLOAT)(*gibbsiter);
  }

  /* Free temporary variables */
  free ( t1   );
  free ( t2   );
  free ( Zu   );
  free ( z    );
  free ( gtmp );
  free ( gold );
  free ( e    );
  free ( o    );
  free ( effects );
  free ( Z    );
  free ( u    );
  free ( sigma);
  free ( y    );

  PutRNGstate();

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect\n");
      for ( i = 0; i < gv->dimu ; i++)
	fprintf(outfp,"%s.%i %f \n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u [i]);

      fclose(outfp);  
  }

  *retval = 0; return ;

}

#ifdef MultPrec

void gs_bayes_iterate_02 ( gs_varset_type* gv ,
			   int*    gibbsburn,
			   int*    gibbsiter,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval         )
{
  gs_UNSIGNED i,j,k,r,m,idx_Z;

  gs_info(gv,1,"Multiple precision version of bayes iteration");

  /* Temporary storage */

  mpf_t tt[4];
  for ( i=0; i<4; i++ ) mpf_init ( tt[i] );

  mpf_t *tt1 = malloc ( gv->NoInd *  sizeof(mpf_t) ) ;
  for ( i=0; i<gv->NoInd; i++ ) mpf_init ( tt1[i] );

  mpf_t *tt2 = malloc ( gv->NoMar *  sizeof(mpf_t) ) ;
  for ( i=0; i<gv->NoMar; i++ ) mpf_init ( tt2[i] );

  mpf_t *ZZu = malloc ( gv->NoInd *  sizeof(mpf_t) ) ;
  for ( i=0; i<gv->NoInd; i++ ) mpf_init ( ZZu[i] );

  mpf_t *zz  = malloc ( gv->NoInd *  sizeof(mpf_t) ) ;
  for ( i=0; i<gv->NoInd; i++ ) mpf_init ( zz[i] );

  mpf_t *ggtmp = malloc ( gv->dimA *  sizeof(mpf_t) ) ;
  for ( i=0; i<gv->dimA; i++ ) mpf_init ( ggtmp[i] );

  mpf_t *ggold = malloc ( gv->dimA *  sizeof(mpf_t) ) ;
  for ( i=0; i<gv->dimA; i++ ) mpf_init ( ggold[i] );

  mpf_t *ee = malloc ( gv->dimA *  sizeof(mpf_t) ) ;
  for ( i=0; i<gv->dimA; i++ ) mpf_init ( ee[i] );

  /* Vector of residuals */
  mpf_t *eo = malloc ( gv->NoInd *sizeof(mpf_t) );
  for ( i=0; i<gv->NoInd; i++ ) mpf_init ( eo[i] );

  /* Vector of ones */
  mpf_t *oo = malloc ( gv->NoInd *sizeof( mpf_t ) );
  for ( i=0; i<gv->NoInd; i++ ) mpf_init_set_d ( oo[i], 1 );

   /* Vector for storing the effects */
  mpf_t *eeffects = malloc ( gv->dimA  *sizeof( mpf_t ) );
  for (j = 0; j < gv->dimA; j++ ) mpf_init_set_d ( eeffects[j] , 1 );

  /* Z u and sigma as matrix in gs_FLOAT */
  mpf_t *ZZ = malloc(gv->NoInd * gv->NoMar * sizeof( mpf_t ) );
  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar; m++) {
      idx_Z   = i*gv->NoMar + m ;
      mpf_init_set_d ( ZZ[idx_Z] , gv->Z[idx_Z]);
    }

  mpf_t *uu =  malloc( gv->dimu * sizeof ( mpf_t ) );
  for ( i = 0; i < gv->dimu ; i++) mpf_init_set_d ( uu[i] , gv->u[i]);

  mpf_t *ssigma =  malloc( gv->dimA * sizeof ( mpf_t ) );
  for ( i = 0; i < gv->dimA ; i++) mpf_init_set_d ( ssigma[i] , gv->sigma[i]) ;

  mpf_t *yy =  malloc( gv->NoInd * sizeof ( mpf_t ) );
  for ( i = 0; i < gv->NoInd ; i++) mpf_init_set_d ( yy[i] ,  gv->y[i]) ;


  /* Temporary storage */

  gs_FLOAT t[4];

  gs_FLOAT *t1 = malloc ( gv->NoInd *  sizeof(gs_FLOAT) ) ;
  if (NULL==t1) ERR_M;

  gs_FLOAT *t2 = malloc ( gv->NoMar *  sizeof(gs_FLOAT) ) ;
  if (NULL==t2) ERR_M;

  gs_FLOAT *Zu = malloc ( gv->NoInd *  sizeof(gs_FLOAT) ) ;
  if (NULL==Zu) ERR_M;

  gs_FLOAT *z = malloc  ( gv->NoInd *  sizeof(gs_FLOAT) ) ;
  if (NULL==z) ERR_M;

  gs_FLOAT *gtmp = malloc ( gv->dimA *  sizeof(gs_FLOAT) ) ;
  if (NULL==gtmp) ERR_M;

  gs_FLOAT *gold = malloc ( gv->dimA *  sizeof(gs_FLOAT) ) ;
  if (NULL==gold) ERR_M;

  /* Vector of residuals */
  gs_FLOAT *e = malloc ( gv->NoInd *sizeof(gs_FLOAT) );
  if (NULL==e) ERR_M;

  /* Vector of ones */
  gs_FLOAT *o = malloc ( gv->NoInd *sizeof(gs_FLOAT) );
  if (NULL==o) ERR_M;
  for (i = 0; i < gv->NoInd; i++ ) o[i] = 1;

  /* Vector for storing the effects */
  gs_FLOAT *effects = malloc ( gv->dimA  *sizeof(gs_FLOAT) );
  if (NULL==effects) ERR_M;
  for (j = 0; j < gv->dimA; j++ ) effects[j] = 0;

  /* Z u and sigma as matrix in gs_FLOAT */
  gs_FLOAT *Z = malloc(gv->NoInd * gv->NoMar * sizeof(gs_FLOAT) );
  if (NULL==Z) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar; m++) {
      idx_Z   = i*gv->NoMar + m ;
      Z[idx_Z] = gv->Z[idx_Z];
    }
  gs_FLOAT *u =  malloc( gv->dimu * sizeof(gs_FLOAT) );
  if ( NULL == u ) ERR_M;
  for ( i = 0; i < gv->dimu ; i++) u [i] = gv->u [i];

  gs_FLOAT *sigma =  malloc( gv->dimA * sizeof(gs_FLOAT) );
  if ( NULL == sigma ) ERR_M;
  for ( i = 0; i < gv->dimA ; i++)  sigma[i] =  gv->sigma[i] ;

  gs_FLOAT *y =  malloc( gv->NoInd * sizeof(gs_FLOAT) );
  if ( NULL == y ) ERR_M;
  for ( i = 0; i < gv->NoInd ; i++)  y[i] =  gv->y[i] ;

  
  // for (j = 0; j < gv->dimA; j++ )  u[j] = 0.0 , sigma[j] = 0.000001 ;  
  

  GetRNGstate();

  for (r=0; r < (gs_UNSIGNED) ( *gibbsburn + *gibbsiter ); r++) 
  {

    if ( round ( (double) r/100) == (double) r/100 ){
      sprintf ( gv->msg,
		"%6i/%i  se2: %5.2f,  b0: %5.2f",
		(int) r,
		(int) ( *gibbsburn + *gibbsiter ),
		(double) mpf_get_d(ssigma[0]),
		(double) mpf_get_d(uu[0])            ) ;

      gs_info(gv,1,gv->msg);
    }

  /* Error variance */

                gs_mm_mul_gF ( gv->NoInd, gv->NoMar,            
			       gv->NoMar, 1        ,
			       Z,
			       &(u[1]),
			       Zu                    ); 

  gs_mm_mul_mpf ( gv->NoInd, gv->NoMar,            
		  gv->NoMar, 1        ,
		  ZZ,
		  &(uu[1]),
		  ZZu                    ); 


                 for (i = 0; i < gv->NoInd; i++ )      
		   e[i] = y[i]  -  Zu[i] - u[0];


 for (i = 0; i < gv->NoInd; i++ ) 
   {
     mpf_sub ( tt[1] , ZZu[i] , uu[0]    );

     // <------------------------------
     // <------------------------------
     // <------------------------------
  
     mpf_sub ( ee[i]  , yy[i] ,  tt[1] );

     Rprintf("%i %f %f %f\n",
	     i,
	     mpf_get_d (ee[i]),
	     mpf_get_d (yy[i]),
	     mpf_get_d (tt[1])
	     );

     // <------------------------------
     // <------------------------------
     // <------------------------------


   }     

                   gs_tmm_mul_gF ( gv->NoInd, 1,        
				   gv->NoInd, 1, 
				   e,
				   e,
				   t1            );
/*
   gs_tmm_mul_mpf ( gv->NoInd, 1,        
		    gv->NoInd, 1, 
		    ee,
		    ee,
		    &(tt[0])          );
*/

                                //gs_FLOAT a = t1[0];
                                //gs_FLOAT b = rchisq( (double) gv->NoInd ) ;
                                //sigma[0]  = a / b ;                            

   // mpf_set_d ( tt[1], rchisq( (double) gv->NoInd ) ) ;
   // mpf_div    ( ssigma[0], tt[0], tt[1] );


  /* Mean */

  gs_mm_mul_gF ( gv->NoInd, gv->NoMar,            
		 gv->NoMar, 1        ,
		 Z,
		 &(u[1]),
		 Zu                    ); 

  gs_tmm_mul_gF (gv->NoInd, 1,
		 gv->NoInd, 1,
		 o,
		 Zu,
		 t1);

  gs_FLOAT sum_y = 0;
  for (i = 0; i < gv->NoInd; i++ ) sum_y += y[i];

  gs_FLOAT expmean =  ( sum_y - t1[0] ) / (gs_FLOAT) gv->NoInd;
  gs_FLOAT sdmean  =  sqrt( sigma[0]    / (gs_FLOAT) gv->NoInd ) ;
  u[0] = rnorm (expmean,sdmean);        
 
  /* Marker variances */

  for (j = 1; j < gv->dimA; j++ ) 
    {
      sigma[j] =  ( u[j] * u[j] ) / rchisq ( 4.012 ); // <----
    }
  
  
  /* Effects */

  for (j=1; j < gv->dimA; j++) gold[j] = u[j];

  for (j=1; j < gv->dimA; j++)
    {
      for (k=1; k < gv->dimA; k++) gtmp[k] = gold[k];
      gtmp[j] = 0 ; 

      for (i=0; i < gv->NoInd; i ++) {
	z[i] = Z[ (i*gv->NoMar) + (j-1) ];
      }

      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,1,
		     z,
		     y,
		     &(t[0]) );
      
      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,gv->NoMar,
		     z,
		     Z,
		     t2 );
      
      gs_mm_mul_gF( 1,gv->NoMar,
		    gv->NoMar,1,
		    t2,
		    &(gtmp[1]),
		    &(t[1]) );

      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,1,
		     z,
		     o,
		     &(t[2]) );

      gs_tmm_mul_gF( gv->NoInd,1,
		     gv->NoInd,1,
		     z,
		     z,
		     &(t[3]) );

      gs_FLOAT expeff = ( t[0] - t[1] - t[2]*u[0] ) /
	              ( t[3] + sigma[0] / sigma[j]);

      gs_FLOAT sdeff  = sqrt ( sigma[0] /
			   ( t[3] + sigma[0] / sigma[j]));

      u[j] = rnorm(expeff,sdeff);                      // <----

    }


  /* Store the effects */
  if ( r >=  (gs_UNSIGNED) *gibbsburn )
    for (j = 0; j < gv->dimA; j++ ) 
      effects[j]  += u[j];

  }

  /* Write back average effects */
  for (j = 0; j < gv->dimA; j++ ) {
    gv->u[j] = effects[j] / (gs_FLOAT)(*gibbsiter);
  }

  /* Free temporary variables */
  free ( t1   );
  free ( t2   );
  free ( Zu   );
  free ( z    );
  free ( gtmp );
  free ( gold );
  free ( e    );
  free ( o    );
  free ( effects );
  free ( Z    );
  free ( u    );
  free ( sigma);
  free ( y    );

  PutRNGstate();

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect\n");
      for ( i = 0; i < gv->dimu ; i++)
	fprintf(outfp,"%s.%i %f \n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u [i]);

      fclose(outfp);  
  }

  *retval = 0; return ;

}


#endif

void(*gs_bayes_iterate_00)() = gs_bayes_iterate_01;



void gs_esteff_lsq_01 (  gs_varset_type *gv ,
			 char**  out_filename,
			 int*    auxfiles, 
			 int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  if ( (1+gv->NoMar) > gv->NoInd ) {
    gs_info(gv,-2,"p > n");  
    *retval = -2; return; }

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  double lambda = 0;
  double * lambdaptr = &lambda;

  gs_lambda_const_01(gv, lambdaptr, nixp, neinp, Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_lsq_01_GV ( char**  out_filename,
			   int*    auxfiles, 
			   int*    retval,
			   char** set_name         )
{
  gs_esteff_lsq_01 ( gs_fdta(*set_name),
		     out_filename,
		     auxfiles, 
		     retval         ) ;
}


void gs_lm_ftest_01 ( gs_varset_type *gv ,
		      char** out_filename,
		      int*   auxfiles, 
		      int*   retval        )
/* F-test for full rank model */
{
  gs_UNSIGNED i,j;
  FILE *outfp;

  /* Allocate memory for the p-values */
  double *Q_val  =  malloc( gv->dimu * sizeof(double) );
  double *F_val  =  malloc( gv->dimu * sizeof(double) );
  gv->p = realloc (gv->p, gv->dimu * sizeof(double) );
  if ( NULL==Q_val ||  NULL==F_val  ||  NULL==gv->p ) ERR_M;

  /* Allocate memory for the contrast vector */
  double *k  =  malloc( gv->dimu * sizeof(double) );
  if ( NULL == k ) ERR_M;

  /* A is X'X */
  gs_mme_restcoeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* (X'X)⁻1  is the inverse of A */
  double * ixpx = malloc( gv->dimu * gv->dimu * sizeof(double) );
  if ( NULL == ixpx ) ERR_M;

  /* Invert A */
  {
    double **a, **y, *col, d;
    gs_UNSIGNED *indx;

    indx=ivector(1,gv->dimA);
    col = vector(1,gv->dimA);
    a = convert_matrix( &(gv->A[0]), 1, gv->dimA, 1, gv->dimA);
    y = convert_matrix( &(ixpx[0]),  1, gv->dimA, 1, gv->dimA);

    ludcmp(a,gv->dimA,indx,&d);

    for ( j=1; j <= gv->dimA; j++ ){
      for ( i=1; i <= gv->dimA; i++) col[i] = 0.0;
      col[j] = 1.0 ;
      lubksb(a,gv->dimA,indx,col);
      for (i=1; i <= gv->dimA; i++) y[i][j] = col[i];
    }

    free_ivector(indx,1);
    free_vector(col,1);
    free_convert_matrix(a,1);
    free_convert_matrix(y,1);
  }

  free(gv->A);
  gv->A =NULL;

  /* b ( from Au=b) is X'y */
  gs_mme_restrhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_UNSIGNED  dfe;
  double kb, kxxk, yy, bxy, sse;
  double * kixtx = malloc ( gv->dimA *sizeof(double) );
  
  gs_tmm_mul (gv->NoInd,1,gv->NoInd,1,gv->y,gv->y,&yy);
  gs_tmm_mul (gv->dimu,1,gv->dimu,1,gv->u,gv->b,&bxy);
  sse = yy - bxy;
  dfe = gv->NoInd - (gv->NoMar+1);

  sprintf(gv->msg, "SSE = %9f, DFE = %lu",sse,dfe);
  gs_info(gv,1,gv->msg);

  /* Calculate p-values */
  for ( i = 0; i < gv->dimu ; i++) 
    {
      for ( j = 0; j < gv->dimu ; j++)  k[j]=0;
      k[i]=1;
      
      gs_tmm_mul ( gv->dimA,1, gv->dimA,1,        k, gv->u,  &kb );
      gs_tmm_mul ( gv->dimA,1, gv->dimA,gv->dimA, k, ixpx,   kixtx);
      gs_tmm_mul ( gv->dimA,1, gv->dimA,1,        kixtx, k,  &kxxk);
      Q_val[i] = kb * (1/kxxk) * kb ;
      F_val[i] = Q_val[i] / (sse/(double)dfe);
      gv->p[i] = st_fprob(F_val[i],1,dfe);
    }

  /* Write solution vector and p-values to output file */
  if (1 == *auxfiles) {

    char sig[5];
   
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect ssh f pvalue sig\n");
      for ( i = 0; i < gv->dimu ; i++)
	{
	  if (gv->p[i]<0.001) strcpy(sig,"***"); 
	  else if (gv->p[i]<0.01) strcpy(sig,"** "); 
	  else if (gv->p[i]<0.05) strcpy(sig,"*  ");
	  else if (gv->p[i]<0.10) strcpy(sig,".  ");
	  else  strcpy(sig,"");

	fprintf(outfp,"%s.%i %f %f %f %f \"%s\"\n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u[i],
		Q_val[i],
		F_val[i],
		gv->p[i],
		sig     );
	}

      fclose(outfp);  
  }

  free(Q_val); 
  free(F_val); 
  free(k); 
  free(ixpx); 
  free(kixtx);

  *retval = 0; return ;
}

void gs_esteff_rmlc_01 (  gs_varset_type *gv ,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlc_01 ( gv,
		      maxiter,precision,hsq,
		      nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmlc_01_GV (  int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmlc_01 ( gs_fdta(*set_name),
                      maxiter,   
		      precision, 
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_esteff_rmlv_01 (  gs_varset_type *gv ,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlv_01 ( gv,
		      maxiter,precision,hsq,
		      nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmlv_01_GV (  int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmlv_01 ( gs_fdta(*set_name),
                      maxiter,   
		      precision,
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_permtest_01 ( gs_varset_type *gv ,
		      char *  scheme       ,
		      int*    nperm,
		      int*    maxiter,
		      double* precision,
		      double* hsq,
		      char** out_filename,
		      int*   auxfiles, 
		      int*   retval        )
{

  if ( RAND_MAX < (int) gv->NoInd )  {
    gs_info(gv,-2,"More individuals than RAND_MAX");
    *retval = -2; return ;
  }

  // Start random number generator, if required
  if (!rgs) { srand( (unsigned) time(NULL) ); rgs=1; }

  // Estimated effects in the permutation runs  
  double * eff_st = malloc( (*nperm)* gv->dimu * sizeof(double));
  if (!eff_st) ERR_M;

  // Omit the intercept
  for (gs_UNSIGNED i =0; i < gv->dimu; i++)
    for (int r = 0; r < *nperm; r++)
      eff_st[(i*(*nperm))+r] = 0;

  int ntr = 0;
  #pragma omp parallel
  {
    ntr = omp_get_num_threads();
  }
  if (0==ntr) ntr=1;


#ifdef use_openblas
  int eins = 1;
  openblas_set_num_threads(eins); 
  goto_set_num_threads(eins); 
  omp_set_num_threads (ntr);
#endif

#ifdef use_mkl
  int eins = 1;
  mkl_set_num_threads(eins); 
  omp_set_num_threads (ntr);
#endif


  // Serialize because data set management uses global variables
  for (int t=0; t<ntr; t++)
    {
      char esname[12];
      sprintf(esname,"dta_t%03i",t);
      gs_set_info_level(sh,imz); 
      gs_varset_type *es = gs_fdta(esname);
      gs_copy_marker_data_01 (es,gv,Sretval);
      if (-2 == *Sretval) { *retval = -2;  }
    }

  // --- once for each thread 
  
  #pragma omp parallel 
  {
    gs_UNSIGNED i,j,idx1,idx2;
    int Srv; int* Sretval = &Srv;
 
    // Pointer to the data set
    int t = omp_get_thread_num();
    char esname[12];
    sprintf(esname,"dta_t%03i",t);
    gs_varset_type *es = gs_fdta(esname);

    // Random permutation
    int *p = malloc( gv->NoInd * sizeof(int) );
    if (!p) { *retval = -2;  }

    // Effect estimates for one marker
    double *e = malloc ( *nperm * sizeof(double) );
    if (!e) { *retval = -2;  }

    gs_UNSIGNED no_unsucc_rep = 0;

    // --- once for each permutation
    #pragma omp for schedule(dynamic)
    for (int r = 0; r < *nperm; r++) 
      {

	int success_rep = 1;
	do {

	for (i =1; i < gv->dimu; i++)
	  {

	    // Random permutaion of the indiviudals 
	    // Knuth Stanford Graph Base 
	    for (int i = 0; i < (int)gv->NoInd; ++i) {
	      int j = rand() % (i + 1);
	      p[i] = p[j];
	      p[j] = i;
	    }
	    
	    // Randomize the alleles of the i the marker
	    for (j = 0; j < gv->NoInd; ++j) {

	      idx1 = (i-1)*gv->NoInd + j;
	      idx2 = (i-1)*gv->NoInd + p[j];

	      es->All1[idx1] = gv->All1[idx2];
	      es->All2[idx1] = gv->All2[idx2];

	    }

	    if ( 0 == strcmp (scheme,"lsq") )
	      { 
		gs_esteff_lsq_01 ( es, nixp, neinp, Sretval ) ;
		if (-2 == *Sretval) { *retval = -2;  }
	      }
	    
	    if ( 0 == strcmp (scheme,"rmlc") )
	      { 
		gs_esteff_rmlc_01(es,maxiter,precision,hsq,nixp,neinp,Sretval);
		if (-2 == *Sretval) { *retval = -2;  }
	      }

	    if ( 0 == strcmp (scheme,"rmlv") )
	      {
		gs_esteff_rmlv_01(es,maxiter,precision,hsq,nixp,neinp,Sretval);
		if (-2 == *Sretval) { *retval = -2;  }
	      }

	    eff_st[(i*(*nperm))+r] = es->u[i];
	  }

	if (-2 == *Sretval) {
	  success_rep = 0;
	  no_unsucc_rep++;
	}

	} while ( ( 0 == success_rep ) && ( no_unsucc_rep < 1000) );

      }
 
    free(e);
    free(p);

  } // parallel

  // Sort the effects
  #pragma omp parallel for 
  for (gs_UNSIGNED i=0; i < gv->dimu; i++) 
    qsort( &eff_st[(i*(*nperm))] ,(*nperm),sizeof(double),gs_dbl_cmp);

  // Allocate memory for the p-values
  double *p_seq  =  malloc( gv->dimu * sizeof(double) );
  double *p_leq  =  malloc( gv->dimu * sizeof(double) );
  gv->p = realloc (gv->p, gv->dimu * sizeof(double) );
  if ( !p_seq ||  !p_leq  ||  !gv->p ) ERR_M;

  // Determine quantiles
  #pragma omp parallel 
  for (gs_UNSIGNED i = 0; i < gv->dimu; i++) 
    {
      gs_UNSIGNED r;
      for (r = 0; r < (gs_UNSIGNED)*nperm; r++)
	if (gv->u[i] < eff_st[(i*(*nperm))+r]) break;
      p_seq[i] = (double) r;

      for (r = 0; r < (gs_UNSIGNED)*nperm; r++)
	if (gv->u[i] > eff_st[(i*(*nperm))+(*nperm-1-r)]) break;
      p_leq[i] = (double) r;

      if ( p_seq[i] < p_leq[i] ) 
	gv->p[i] =  p_seq[i]/(double)(*nperm);
      else 
	gv->p[i] =  p_leq[i]/(double)(*nperm);
    }

#ifdef use_openblas
  openblas_set_num_threads(ntr); 
  goto_set_num_threads(ntr); 
#endif

#ifdef use_mkl
  mkl_set_num_threads(ntr); 
#endif


  /* Write solution vector and p-values to output file */
  if (1 == *auxfiles) {

    FILE *outfp;
    char sig[5];

    if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

    fprintf(outfp,"effect  p.seq p.leq pvalue sig\n");
    for (gs_UNSIGNED i = 0; i < gv->dimu ; i++)
      {
	if (gv->p[i]<0.001) strcpy(sig,"***"); 
	else if (gv->p[i]<0.01) strcpy(sig,"** "); 
        else if (gv->p[i]<0.05) strcpy(sig,"*  ");
	else if (gv->p[i]<0.10) strcpy(sig,".  ");
	else  strcpy(sig,"");

	fprintf(outfp,"%s.%i %f %f %f %f \"%s\"\n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u[i],
		p_seq[i],
		p_leq[i],
		gv->p[i],
		sig     );
      }
      fclose(outfp);  
  }

  free(p_seq);
  free(p_leq);
  free(eff_st);

  // gs_mt_matrix_ops();
  // gs_mt_lu_ops();

  *retval = 0; return ;

}


void gs_testeff_lsq_01 (  gs_varset_type *gv ,
			  int*    nperm,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_esteff_lsq_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  if ( 0 == *nperm)
    {
      gs_lm_ftest_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }
  else
    {
      char scheme[]="lsq";
      int mi=0; int* maxiter = &mi;
      double pr=0; double* precision = &pr;
      double _hsq = 0.9; double* hsq=&_hsq;
      gs_permtest_01 ( gv,scheme,nperm,maxiter,precision,hsq,
		       out_filename,auxfiles,retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  *retval = 0; return ;
}

void gs_testeff_lsq_01_GV ( int*    nperm,
			    char**  out_filename,
			    int*    auxfiles, 
			    int*    retval,
			    char** set_name         )
{
  gs_testeff_lsq_01 ( gs_fdta(*set_name),
		      nperm,
	              out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_testeff_rmlc_01 (  gs_varset_type *gv ,
			   int*    nperm,
			   int*    maxiter,
			   double* precision,
			   double* hsq,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_esteff_rmlc_01 ( gv, maxiter, precision,hsq, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  char scheme[]="rmlc";
  gs_permtest_01 (gv,scheme,nperm,maxiter,precision,hsq,out_filename,
		  auxfiles,retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}


void gs_testeff_rmlc_01_GV ( int*    nperm, 
			     int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_testeff_rmlc_01 ( gs_fdta(*set_name),
		       nperm,
		       maxiter,   
		       precision, 
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}

void gs_testeff_rmlv_01 (  gs_varset_type *gv ,
			   int*    nperm,
			   int*    maxiter,
			   double* precision,
			   double* hsq,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_esteff_rmlv_01 ( gv, maxiter, precision, hsq, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  char scheme[]="rmlv";
  gs_permtest_01 (gv,scheme,nperm,maxiter,precision,hsq,
		  out_filename,auxfiles,retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}


void gs_testeff_rmlv_01_GV ( int*    nperm, 
			     int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_testeff_rmlv_01 ( gs_fdta(*set_name),
		       nperm,
		       maxiter,   
		       precision, 
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}


void gs_esteff_rmlr_01 (  gs_varset_type *gv ,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlr_01 ( gv,
		      maxiter,precision,hsq,
		      nixp,neinp,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmlr_01_GV (  int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmlr_01 ( gs_fdta(*set_name),
                      maxiter,   
		      precision, 
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}

void gs_esteff_rmla_01 (  gs_varset_type *gv ,
			  double* alpha,
                          int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  gs_build_Z_01 ( gv , nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmla_01 ( gv,
		      alpha,
		      maxiter,precision,hsq,
		      nixp,neinp,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_rmla_01_GV (  double* alpha,
			     int*    maxiter,
			     double* precision,
			     double* hsq,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_rmla_01 ( gs_fdta(*set_name),
		      alpha,
                      maxiter,   
		      precision, 
		      hsq,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}



void gs_esteff_bayes_01 (  gs_varset_type *gv ,
			   int*    maxiter,
			   double* precision,
			   int*    gibbsburn,
			   int*    gibbsiter,
			   char**  out_filename,
			   int*    auxfiles, 
			   int*    retval         )
{
  int Srv; int* Sretval = &Srv;

  double _hsq = 0.9; double* hsq=&_hsq;

  gs_esteff_rmlv_01 ( gv,
		      maxiter,precision,hsq,
		      nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  (*gs_bayes_iterate_00) ( gv,
			   gibbsburn, gibbsiter,
			   out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;
}

void gs_esteff_bayes_01_GV (  int*    maxiter,
			     double* precision,
			     int*    gibbsburn,
			     int*    gibbsiter,
			     char**  out_filename,
			     int*    auxfiles, 
			     int*    retval,
			     char** set_name         )
{
  gs_esteff_bayes_01 ( gs_fdta(*set_name),
                      maxiter,   
		      precision, 
		      gibbsburn,
		      gibbsiter,
		      out_filename,
		      auxfiles, 
		      retval         ) ;
}


void gs_cross_validation_01 ( gs_varset_type *gv   ,
			      char ** scheme       ,
			      int*    n_estimation ,
			      int*    n_runs       , 
			      double* hsq          , 
			      double* alpha        ,
			      int*    maxiter      ,
			      double* precision    ,
			      int*    gibbsburn    ,
			      int*    gibbsiter    ,
			      char**  out_filename_r , 
			      char**  out_filename_u , 
			      int*    auxfiles     , 
			      int*    retval       )
 { 
   if ( RAND_MAX < (int) gv->NoInd )  {
       gs_info(gv,-2,"More individuals than RAND_MAX");
       *retval = -2; return ;
     }

   if ( *n_estimation > (int) gv->NoInd )  {
       gs_info(gv,-2,"Training is set larger than base population");
       *retval = -2; return ;
     }


   if (  !(    ( 0 == strcmp (*scheme,"rrwe") )
	     ||( 0 == strcmp (*scheme,"rrwr") )
	     ||( 0 == strcmp (*scheme,"rrwa") )
	     ||( 0 == strcmp (*scheme,"rmlc") )
	     ||( 0 == strcmp (*scheme,"rmlv") )
             ||( 0 == strcmp (*scheme,"rmlr") ) 
             ||( 0 == strcmp (*scheme,"rmla") ) 
             ||( 0 == strcmp (*scheme,"bayes") ) 
      )   )
     {
       gs_info(gv,-2,"Unknown estimation scheme");
       *retval = -2; 
       return;
     }

   int ntr = 0;
   #pragma omp parallel
   {
     ntr = omp_get_num_threads();
   }
   if (0==ntr) ntr=1;


#ifdef use_openblas
  int eins = 1;
  openblas_set_num_threads(eins); 
  goto_set_num_threads(eins); 
  omp_set_num_threads (ntr);
#endif

#ifdef use_mkl
  int eins = 1;
  mkl_set_num_threads(eins); 
  omp_set_num_threads (ntr);
#endif

   if (!rgs) { srand( (unsigned) time(NULL) ); rgs=1; }

   int CVNoRuns = *n_runs;  

   /* Results */
   double* CVcor = malloc ( CVNoRuns * sizeof(double) );
   if (NULL == CVcor) ERR_M;
   for (gs_UNSIGNED c = 0; c < (gs_UNSIGNED) CVNoRuns ; c++) CVcor[c] = 0;

   gs_FLOAT* CVest = malloc ( CVNoRuns * (1+gv->NoMar) *sizeof(gs_FLOAT) );


     /* --- once for each thread */
     #pragma omp parallel
     {
       int Srv; int* Sretval = &Srv;
 
       gs_varset_type ES; gs_varset_type *es = &ES;
       gs_init (es);

       gs_varset_type VS; gs_varset_type *vs = &VS;
       gs_init (vs);
       
       gs_UNSIGNED i;
       double cor;       
       double sy=0, syy=0, sx=0, sxx=0, sxy=0, k;

       #pragma omp for
       for (int r=0 ; r < CVNoRuns; r ++)  
	 /* --- once for each CVRun */
	 {

	   cor = 0;

	   /* Build estimation and validation set */
	   gs_build_esvs_01(gv,es,vs,n_estimation, 
			 nixp,nixp,nixp,nixp,neinp, Sretval);
	   if (-2 == *Sretval) { *retval = -2;  }
	   
	   /* Estimate effects in estimation set */
	   if ( 0 == strcmp (*scheme,"rrwe") )
	     {
	        gs_esteff_rrwe_02 ( es, hsq, nixp,neinp, Sretval ) ;
		if (-2 == *Sretval) { *retval = -2;  }
	     }
	   if ( 0 == strcmp (*scheme,"rrwr") )
	     {
	        gs_esteff_rrwr_01 ( es, hsq, nixp,neinp, Sretval ) ;
		if (-2 == *Sretval) { *retval = -2;  }
	     }
	   if ( 0 == strcmp (*scheme,"rrwa") )
	     {
                gs_esteff_rrwa_01 ( es, hsq, alpha, nixp,neinp, Sretval ) ;
		if (-2 == *Sretval) { *retval = -2;  }
	     }
	   if ( 0 == strcmp (*scheme,"rmlc") )
	     {
	       gs_esteff_rmlc_01(es,maxiter,precision,hsq,nixp,neinp,Sretval);
		if (-2 == *Sretval) { *retval = -2;  }
	     }
	   if ( 0 == strcmp (*scheme,"rmlv") )
	     {
	        gs_esteff_rmlv_01(es,maxiter,precision,hsq,nixp,neinp,Sretval);
		if (-2 == *Sretval) { *retval = -2;  }
	     }
	   if ( 0 == strcmp (*scheme,"rmlr") )
	     {
	        gs_esteff_rmlr_01(es,maxiter,precision,hsq,nixp,neinp,Sretval);
		if (-2 == *Sretval) { *retval = -2;  }
	     }
	   if ( 0 == strcmp (*scheme,"rmla") )
	     {
                gs_esteff_rmla_01(es,alpha,maxiter,precision,hsq,nixp,neinp,Sretval);
		if (-2 == *Sretval) { *retval = -2;  }
	     }
	   if ( 0 == strcmp (*scheme,"bayes") )
	     {
	        gs_esteff_bayes_01(es,maxiter,precision,gibbsburn,gibbsiter,
				   nixp,neinp,Sretval);
		if (-2 == *Sretval) { *retval = -2;  }
	     }

	   /* Save effect estimates*/
	   for (i =0; i < es->dimu; i++) CVest[(r*es->dimu)+i] = es->u [i];

	   /* Estimate breeding values in validation set */
	   gs_estimate_bv_01 (es,vs,nixp,neinp,Sretval);
	   if (-2 == *Sretval) { *retval = -2;}

	   /* Correlation estimated - observed */

	   sy=0, syy=0, sx=0, sxx=0, sxy=0;
	   k = (double) vs->NoInd;
	   
	   for ( i = 0; i < vs->NoInd; i++ )  {
	       sx  += vs->est[i] ;
	       sxx += vs->est[i] * vs->est[i] ;
	       sy  += vs->y[i] ;
	       syy += vs->y[i] * vs->y[i] ;
	       sxy += vs->est[i] * vs->y[i] ;
	     }
	   cor =    (      k*sxy - sx*sy            ) /
	       sqrt ( (k*sxx-sx*sx) * (k*syy-sy*sy) ) ;

	   CVcor[r] = cor;

	 } // for r
       
       gs_reset(es);
       gs_reset(vs);

     } // parallel

#if defined(use_openblas) && defined(MT)
  openblas_set_num_threads(ntr); 
  goto_set_num_threads(ntr); 
#endif

#if defined(use_mkl) && defined(MT)
  mkl_set_num_threads(ntr); 
#endif

     gs_UNSIGNED i,r;
     gs_FLOAT sum;

     /* Allocate memory for the solution vector */
     int Srv; int* Sretval = &Srv;
     gs_esteff_rrwe_02 ( gv, hsq, nixp,neinp, Sretval ) ;
     if (-2 == *Sretval) { *retval = -2;  }

     /* Average effects estimated in CV */
     for ( i=0 ; i < gv->dimu; i++)
       {
	 sum = 0;
	 for ( r=0; r < (gs_UNSIGNED)CVNoRuns ; r++)
	   sum += CVest[(r*gv->dimu)+i] ;
	 gv->u[i] = sum / (gs_FLOAT)CVNoRuns;
       }


   /* Write results to output file */
   if (1 == *auxfiles)
     {
       FILE *fp;
       if ( (fp = fopen(*out_filename_r,"w")) == NULL ) ERR_F;
       fprintf(fp,"cor \n");
       for (int c = 0; c < CVNoRuns ; c++)
       	 fprintf(fp,"%i %f\n",c, CVcor[c] );
       fclose(fp);  

       if ( (fp = fopen(*out_filename_u,"w")) == NULL ) ERR_F;

       fprintf(fp,"effect\n");
       for ( i = 0; i < gv->dimu ; i++)
	 fprintf(fp,"%s.%i %f \n",
		 gv->EffNme[i],
		 gv->EffAll[i],
		 gv->u [i]);
       fclose(fp);  
     }
   
   free(CVcor);
   free(CVest);

   *retval = 0; return;

 } 


void gs_cross_validation_01_GV ( char ** scheme       ,
				 int*    n_estimation ,
				 int*    n_runs       , 
				 double* hsq          , 
				 double* alpha        ,
				 int*    maxiter      ,
				 double* precision    ,
				 int*    gibbsburn,
				 int*    gibbsiter,
				 char**  out_filename_r , 
				 char**  out_filename_u , 
				 int*    auxfiles     , 
				 int*    retval       ,
				 char** set_name      )
 {
   gs_cross_validation_01 ( gs_fdta(*set_name),
			    scheme       ,
			    n_estimation ,
			    n_runs       , 
			    hsq          , 
			    alpha        ,
			    maxiter      ,
			    precision    ,
			    gibbsburn    ,
			    gibbsiter    ,
			    out_filename_r , 
			    out_filename_u , 
			    auxfiles     , 
			    retval           ) ;
 }

#ifdef MultPrec

void gs_set_precision_01 (const int prec)
{
  if (0==prec)  
    {
      gs_info(gv,0,"Standard precision");

      gs_bayes_iterate_00  = gs_bayes_iterate_01;
    }
  else
    {
      mpf_set_default_prec(prec);
      sprintf(gv->msg, "Precision set to %i bits", prec);
      gs_info(gv,0,gv->msg);

      gs_bayes_iterate_00  = gs_bayes_iterate_02;
    }
}
#endif

#ifndef MultPrec

void  gs_set_precision_01 (const int prec)
{
  int dummy = prec; dummy++; //avoid compile warning
  gs_info(gv,-1,"Multiple precision not supported");
}

#endif

void gs_set_precision_01_GV (int* prec)
{
  int p = *prec;
  gs_set_precision_01 (p) ;
}

/* Nov 2011 */
/* Functions for plotting graphical genotypes */

/* qsort pointers to mappoints comparison function*/
int gs_mappoint_cmp(const void *a, const void *b) 
{ 
  int retval = 0;

  const gs_mappoint_type **ia = (const gs_mappoint_type **)a;
  const gs_mappoint_type **ib = (const gs_mappoint_type **)b;

  if ( (*ia)->chrom > (*ib)->chrom ) {
    retval = 1;
  } 
  else {
    if ( (*ia)->chrom < (*ib)->chrom ) {
      retval = -1 ;
    } 
    else {
      if ( (*ia)->pos > (*ib)->pos ) {
	retval = 1;
      } 
      else {
	if ( (*ia)->pos < (*ib)->pos ) {
	  retval = -1 ;
	}
      }
    }
  }

  return retval;
} 


void gs_reset_all_but_marker_matrix (gs_varset_type *gv)
{
  gs_UNSIGNED NoInd = gv->NoInd; 
  gs_UNSIGNED NoMar = gv->NoMar; 

  char (*IndName) [gs_MAXNME] = gv->IndName; gv->IndName = NULL;  
  char (*MarName) [gs_MAXNME] = gv->MarName; gv->MarName = NULL; 

  int * All1  = gv->All1; gv->All1 = NULL;
  int * All2  = gv->All2; gv->All2 = NULL;

  gs_reset(gv);

  gv->NoInd = NoInd; 
  gv->NoMar = NoMar; 

  gv->IndName = IndName;  
  gv->MarName = MarName; 

  gv->All1 = All1;
  gv->All2 = All2;

}

void gs_reset_chrom_stats()
{
  gs_UNSIGNED c;

  for (c = 0; c < gv->NChrom; c ++ ) 
    {
      gv->chrom[c].rf = realloc (gv->chrom[c].rf, 0); 
      gv->chrom[c].q  = realloc (gv->chrom[c].q,  0);  
    }

  gv->chrom = realloc( gv->chrom, 0 ); gv->chrom = NULL; 

  gv->NChrom = 0;
}

void gs_chrom_stats_01 (gs_varset_type *gv ,
			char** out_filename,
			int* auxfiles,
			int* retval          )
{
  gs_UNSIGNED c, i;

  if ( 0 == gv->NLoc ) { 
      sprintf (gv->msg, "No linkage map was loaded");
      gs_info(gv,-2,gv->msg);
      *retval = -2;   return;
    }

  gs_reset_chrom_stats();

  /* Count chromosomes and memory allocation*/

  gv->NChrom = 1;
  for ( i = 1 ; i < gv->NLoc; i++ )
    if  (gv->map[i-1]->chrom != gv->map[i]->chrom) gv->NChrom++;

  gv->chrom = realloc( gv->chrom, gv->NChrom * sizeof (gs_chrom_type) );

  for (c = 0; c < gv->NChrom; c ++ ) 
    {
      gv->chrom[c].rf = NULL;
      gv->chrom[c].q  = NULL;
    }

  sprintf ( gv->msg, "No. of chromosomes: %lu", gv->NChrom );
  gs_info ( gv, 1, gv->msg );

  /* First locus of each chromosome and number of loci  */

  gv->chrom[0].NLocC = 1;
  gv->chrom[0].first_loc   = &(gv->map[0]);
  gv->chrom[0].first_loc_i = 0;
  for ( c=0, i=1 ; i < gv->NLoc; i++ ) 
    {
      if ( gv->map[i-1]->chrom == gv->map[i]->chrom )
	{
	  gv->chrom[c].NLocC ++ ;
	}
      else
	{ 
	  c ++ ;
	  gv->chrom[c].first_loc   = &(gv->map[i]);
	  gv->chrom[c].first_loc_i = i;
	  gv->chrom[c].NLocC = 1;
	}
    }


  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      gs_mappoint_type ** p;
 
      for ( c=0; c < gv->NChrom; c++ )
	{
	  p =  gv->chrom[c].first_loc ;
	  sprintf( gv->msg,"%lu %lu %f \n",
		   1+c,
		   gv->chrom[c].NLocC,
		   ( p[gv->chrom[c].NLocC-1] ) -> pos
		   );
	  fprintf(outfp,"%s",gv->msg);
	}
      
      fclose(outfp);
    }

  *retval=0; return;
}

void gs_chrom_stats_01_GV ( char**  out_filename,
			    int*    auxfiles,
			    int*    retval,
			    char**  set_name)
{
  gs_chrom_stats_01 ( gs_fdta(*set_name),
		      out_filename,
		      auxfiles,
		      retval         ) ;
}

gs_FLOAT gs_rf_haldane( gs_FLOAT d)
{
  return  ( (1 - exp(-2*d) ) /2  );
}

gs_FLOAT gs_rf_kosambi( gs_FLOAT d)
{
  return ( exp(4*d) - 1 ) / ( exp(4*d) + 1 ) / 2;
}

void gs_calc_rf_01 ( gs_varset_type *gv  ,
		     char** map_function ,
		     char** out_filename ,
		     int* auxfiles       ,
		     int* retval           )
{
  gs_UNSIGNED c, i, j, idx;
  gs_FLOAT d;
  gs_mappoint_type ** p;

  if ( 0 == gv->NChrom ) { 
      sprintf (gv->msg, "No chromosome statistics available");
      gs_info(gv,-2,gv->msg);
      *retval = -2;   return;
    }

  gs_FLOAT (*gs_rf) (gs_FLOAT d);

  if ( 0 == strcmp(*map_function,"Haldane"))
    {
      gs_rf = gs_rf_haldane;
    }

  else if ( 0 == strcmp(*map_function,"Kosambi"))
    {
      gs_rf = gs_rf_kosambi;
    }

  else 
    {
      gs_info(gv,-2,"Unknown map function");
      *retval = -2; return;
    }

  /* Allocate memory for recombination frequencies*/
  for (c = 0; c < gv->NChrom; c ++ ) 
    {
      gv->chrom[c].rf = realloc (gv->chrom[c].rf, 
	gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT) );
      gv->chrom[c].q  = realloc (gv->chrom[c].q , 
	gv->chrom[c].NLocC * gv->chrom[c].NLocC * sizeof(gs_FLOAT) ); 
    }

  for ( c=0; c < gv->NChrom; c++ )
    {
       p =  gv->chrom[c].first_loc ;
       for  ( i = 0 ; i < gv->chrom[c].NLocC ; i ++)
	 for  ( j = 0 ; j < gv->chrom[c].NLocC ; j ++)
	   {
	     d = fabs ( (p[i])->pos - (p[j])->pos );

	     idx = i * gv->chrom[c].NLocC + j;
	     (gv->chrom[c].rf) [idx] = gs_rf (d/100) ;
	   }
    }

  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( c=0; c < gv->NChrom; c++ )
	{
	  for  ( i = 0 ; i < ( gv->chrom[c].NLocC - 1 ) ; i ++)
	    for  ( j = i+1 ; j < gv->chrom[c].NLocC ; j ++)
	      {
		idx = i * gv->chrom[c].NLocC + j;

		sprintf( gv->msg,"%lu %lu %lu %f \n",
			 1+c,
			 1+i,
			 1+j,
			 (double) (gv->chrom[c].rf) [idx]     );
		fprintf(outfp,"%s",gv->msg);

	      }
	  
	}

      fclose(outfp);
    }

  *retval=0; return;

 }


void gs_calc_rf_01_GV (	char** map_function ,
			char** out_filename ,
			int*   auxfiles     ,
			int*   retval       ,
			char** set_name       )
{
  gs_calc_rf_01 ( gs_fdta(*set_name),
		  map_function      ,
		  out_filename      ,
		  auxfiles          ,
		  retval              ) ;
}

gs_FLOAT gs_q_dh (int* t, gs_FLOAT r)
{
  return 0.5 + ( 1 - 2*r) / 2 * pow ( (double)(1-r), (double)*t);
}

gs_FLOAT gs_q_ssd (int* t, gs_FLOAT r)
{
  return 0.5 + ( 1 - 2*r) / ( 2 - 4*r ) * pow ( (double)(1-r), 1+(double)(*t) );
}

void gs_calc_q_01 ( gs_varset_type *gv   ,
		    char**  pop_type     ,
		    int*    t            ,
		    char** out_filename  ,
		    int* auxfiles        ,
		    int* retval            )
{
  gs_UNSIGNED c, i, j, idx;

  if  ( ( 0 == gv->NChrom ) || ( NULL == gv->chrom[0].rf ) ) { 
      sprintf (gv->msg, "No recombination frequencies available");
      gs_info(gv,-2,gv->msg);
      *retval = -2;   return;
    }

  gs_FLOAT (*gs_q) (int* t, gs_FLOAT r);

  if ( 0 == strcmp(*pop_type,"DH"))
    {
      gs_q = gs_q_dh;
    }

  else if ( 0 == strcmp(*pop_type,"SSD"))
    {
      gs_q = gs_q_ssd;
    }

  else 
    {
      gs_info(gv,-2,"Unknown population type");
      *retval = -2; return;
    }

  
  for ( c=0; c < gv->NChrom; c++ )
    {
      for  ( i = 0 ; i < gv->chrom[c].NLocC ; i ++)
	for  ( j = 0 ; j < gv->chrom[c].NLocC ; j ++)
	  {
	    idx = i * gv->chrom[c].NLocC + j;
	    (gv->chrom[c].q)[idx] =  gs_q (t,(gv->chrom[c].rf)[idx]);
	  }
    }
  
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( c=0; c < gv->NChrom; c++ )
	{
	  for  ( i = 0 ; i < ( gv->chrom[c].NLocC - 1 ) ; i ++)
	    for  ( j = i+1 ; j < gv->chrom[c].NLocC ; j ++)
	      {
		idx = i * gv->chrom[c].NLocC + j;

		sprintf( gv->msg,"%lu %lu %lu %f \n",
			 1+c,
			 1+i,
			 1+j,
			 (double) (gv->chrom[c].q) [idx]     );
		fprintf(outfp,"%s",gv->msg);

	      }
	  
	}

      fclose(outfp);
    }

  *retval=0; return;

 }


void gs_calc_q_01_GV ( char**  pop_type     ,
		       int*    t            ,
		       char**  out_filename ,
		       int*    auxfiles     ,
		       int*    retval       ,
		       char**  set_name       )
{
  gs_calc_q_01 ( gs_fdta(*set_name),
		 pop_type,
		 t,
		 out_filename,
		 auxfiles,
		 retval         ) ;
}


void gs_replace_comma (char *p)
{
  char* q = strchr(p,',');
  if (NULL != q) *q = '.';
}

void gs_read_map_01 ( gs_varset_type *gv ,
		      char**  map_filename,
		      int*    skip,
		      char**  fmt,
		      double* map_stretch,
		      char**  out_filename,
		      int*    auxfiles,
		      int*    retval        )
{

  gs_UNSIGNED  
    MarSize = 0, /* Length of marker name in input file     */
    ClsSize = 0, /* Lenght of class name in input file      */
    NoLines = 0, /* Number of lines in the input file       */
    lne     = 0, m;

  FILE *fp;

  char*  line = NULL;   char** linep = &line;
  size_t n    = 0;      size_t* np   = &n;

  char *p ;

  if ( !( (0 == strcmp(*fmt,"cpms")) ||
  	  (0 == strcmp(*fmt,"mcp"))    ) )
    {
      gs_info(gv,-2,"Implemented map formats: cpms, mcp");
      *retval = -2; return ;
    }

  gs_reset_chrom_stats();

  /* (1) Read in the map */

  /* Determine number of lines in the map */
  if ( NULL == (fp = fopen(*map_filename,"r")) ) ERR_F;

  for (m=0; m < (gs_UNSIGNED)(*skip); m++ )
    getline ( linep , np , fp ) ;


  if (0 == strcmp(*fmt,"cpms"))
    {
      while (!feof(fp))
	{
	  if ( 0 < getline ( linep , np , fp ) )
	    {
	      p = strtok (line,gs_tok_del);
	      if (( NULL == p )   || ( 0 == (int)atoi(p)) )
		{gs_info(gv,-2,"Error in map file");*retval=-2;return;}

	      p = strtok ('\0',gs_tok_del);
	      if  ( NULL == p )  
		{gs_info(gv,-2,"Error in map file");*retval=-2;return;}

	      p = strtok ('\0',gs_tok_del); 
	      if ( NULL == p ) {gs_info(gv,-2,"Error in map file");*retval=-2;return;}
	      if ( strlen(p) > MarSize) MarSize = strlen(p); 

	      p = strtok ('\0',gs_tok_del);
	      if ( NULL == p ) {gs_info(gv,-2,"Error in map file");*retval=-2;return;}
	      if ( strlen(p) > ClsSize) ClsSize = strlen(p);
	      NoLines++;
	    }
	}
    }

  if (0 == strcmp(*fmt,"mcp"))
    {
      while (!feof(fp))
	{
	  if ( 0 < getline ( linep , np , fp ) )
	    {
	      p = strtok (line,gs_tok_del);
	      if ( NULL == p ) {gs_info(gv,-2,"Error in map file");*retval=-2;return;}
	      if ( strlen(p) > MarSize) MarSize = strlen(p); 

	      p = strtok ('\0',gs_tok_del);
	      if (( NULL == p ) || ( 0 == (int)atoi(p)) )
		{gs_info(gv,-2,"Error in map file");*retval=-2;return;}


	      p = strtok ('\0',gs_tok_del); 
	      if ( NULL == p ) 
		{gs_info(gv,-2,"Error in map file");*retval=-2;return;}

	      NoLines++;
	    }
	}
      ClsSize = 1;
    }


  fclose(fp);

  /* Messages on data size*/
  sprintf(gv->msg,
  	  "Description length: Mar %lu, Class %lu",
  	  MarSize, ClsSize);
  gs_info(gv,1,gv->msg);
  sprintf ( gv->msg, "No. of markers in map file: %lu", NoLines );
  gs_info ( gv, 1, gv->msg );
  if ( ( MarSize > gs_MAXSTR) || ( ClsSize > gs_MAXSTR) )
    {
      gs_info(gv,-2,"Decription too long");
      *retval = -2; return ;
    }
  
  /* Allocate memory for the map */

  gv->NLoc  = NoLines;
  gv->map_S = realloc (gv->map_S, gv->NLoc*sizeof(gs_mappoint_type) );
  gv->map   = realloc (gv->map,   gv->NLoc*sizeof(gs_mappoint_type*) );
  gv->mdl   = realloc (gv->mdl,   gv->NLoc*sizeof(gs_UNSIGNED) );
  if ( (NULL == gv->map) || (NULL == gv->map_S) || (NULL == gv->mdl) ) ERR_M; 

  /* Read in the map */
  if ( NULL == (fp = fopen(*map_filename,"r")) ) ERR_F;

  for (m=0; m < (gs_UNSIGNED)(*skip); m++ )
    getline ( linep , np , fp ) ;

  lne = 0;

  if (0 == strcmp(*fmt,"cpms"))
    {
      while (!feof(fp))
	{
	  if (  0 < getline ( linep , np , fp ) )
	    {
	      p = strtok (line,gs_tok_del);
	      gv->map_S[lne].chrom = (int)atoi (p);
	      
	      p = strtok ('\0',gs_tok_del);
	      gs_replace_comma(p);
	      gv->map_S[lne].pos   = *map_stretch * (double)atof(p);
	      
	      p = strtok ('\0',gs_tok_del);
	      strcpy(gv->map_S[lne].name,p);
	      
	      p = strtok ('\0',gs_tok_del);
	      strcpy(gv->map_S[lne].class,p);
	      
	      lne++;
	    }
	}
    }

  if (0 == strcmp(*fmt,"mcp"))
    {
      while (!feof(fp))
	{
	  if (  0 < getline ( linep , np , fp ) )
	    {
	      p = strtok (line,gs_tok_del);
	      strcpy(gv->map_S[lne].name,p);

	      p = strtok ('\0',gs_tok_del);
	      gv->map_S[lne].chrom = (int)atoi (p);
	      
	      p = strtok ('\0',gs_tok_del);
	      gs_replace_comma(p);
	      gv->map_S[lne].pos   = *map_stretch * (double)atof(p);

	      strcpy(gv->map_S[lne].class,"m");
	      
	      lne++;
	    }
	}
    }


  fclose(fp);
  free(line);

  /* Array of pointers to mappoints */
  for (lne=0 ; lne < NoLines ; lne++ )
      gv->map[lne] = &(gv->map_S[lne]);

  /* Sort and correct identical marker positions */
  int jitter = 1;

  while ( 1 == jitter )
    {
      jitter = 0;

      qsort ( gv->map, 
	      gv->NLoc, 
	      sizeof (gs_mappoint_type *), 
	      gs_mappoint_cmp);

      for (lne=1 ; lne < NoLines ; lne++ )
      	{
      	  if ( ( (*gv->map[lne-1]).pos   == (*gv->map[lne]).pos) &&
      	       ( (*gv->map[lne-1]).chrom == (*gv->map[lne]).chrom)   )
      	    {
      	      (*gv->map[lne]).pos += 0.00123456;
      	      jitter = 1;
      	    }

    /*   for (lne=1 ; lne < NoLines ; lne++ ) */
    /* 	{ */
    /* 	  if  ( 0 == gs_mappoint_cmp( &(gv->map[lne-1]) , &(gv->map[lne])) ) */
    /* 	    {  */
    /* 	      Rprintf("%i %f \n",  */
    /* 		      (*gv->map[lne]).chrom,  */
    /* 		      (*gv->map[lne]).pos); */
    /* 	      (*gv->map[lne]).pos += 0.00123456; */
    /* 	      jitter = 1; */
    /* 	    } */
    /* 	} */

	}
    }

  /* Find corresponding data lines in the genotype matrix */
  gs_UNSIGNED NoMarMatching = 0;
  for (lne=0; lne <  gv->NLoc; lne++) {
    gv->mdl[lne] = gs_UNSIGNED_MAX;
    for ( m=0; m < gv->NoMar; m++ ) {
      if ( 0 == strcmp ( gv->MarName[m] , (gv->map[lne])->name ) ) {
	gv->mdl[lne] = m;
	NoMarMatching++;
      }
    }
  }

  /* Check for duplicated marker names */
  gs_UNSIGNED l1,l2,duplicates;
  duplicates=0;
  for (l1=0; l1 < gv->NLoc-1; l1++) 
    {
      for (l2=l1+1; l2 < gv->NLoc; l2++)
	{
	  if ( 0 == strcmp ( (gv->map[l1])->name , (gv->map[l2])->name ) ) 
	    {
	      gv->mdl[l1] = gs_UNSIGNED_MAX;
	      NoMarMatching--;
	      sprintf( gv->msg, "Duplicate marker name %s", (gv->map[l1])->name );
	      gs_info(gv,-1,gv->msg);
	      duplicates=1;
	    }
	}
    }

  if (1==duplicates){
    sprintf ( gv->msg, "Remove duplicated markers from map" );
    gs_info(gv,-2,gv->msg);
    *retval = -2;
    return ;
  }


  /* (2) Generate a new data matrix */

  /* Temporary storage marker data */
  gs_UNSIGNED NoInd = gv -> NoInd;    // Note: No changes here!
  gs_UNSIGNED NoMar = NoMarMatching; 

  char (*MarName) [gs_MAXNME] = malloc(NoMar*sizeof(char[gs_MAXNME])); 

  int * All1  =  malloc( NoMar * NoInd * sizeof(int) );
  int * All2  =  malloc( NoMar * NoInd * sizeof(int) );
  
  // Phenotypic data
  double *y = gv->y; gv->y = NULL;

  if ( (NULL == All1)  || (NULL == All2) || ( NULL == MarName )) ERR_M;

  /* Temporary strorage map */
  gs_UNSIGNED NLoc          = NoMarMatching;
  gs_mappoint_type  * map_S = malloc ( NLoc*sizeof(gs_mappoint_type) );
  gs_mappoint_type ** map   = malloc ( NLoc*sizeof(gs_mappoint_type*));
  gs_UNSIGNED * mdl         = malloc ( NLoc * sizeof(gs_UNSIGNED) );
  if ( (NULL == map) || (NULL == mdl) ) ERR_M; 

  /* Write data to temporary storage */
  gs_UNSIGNED m_n = 0;                  // Row in the new data matrix
  int * p1 ;                            // Address of a row in old matrix
  int * p2 ;                            // Address of a row in new matrix
  for (lne=0; lne <  gv->NLoc; lne++)   // Rows of the complete map 
    {
      m = gv->mdl[lne] ;                // Row of the original data matrix
      if (gs_UNSIGNED_MAX != m)         // Marker data is available
	{
	  p1 = &(gv->All1[m*NoInd]);    // First allele
	  p2 = &(All1[m_n*NoInd]);
	  memcpy( p2 , p1, (size_t)NoInd * sizeof(int) );
	  p1 = &(gv->All2[m*NoInd]);    // Second allele
	  p2 = &(All2[m_n*NoInd]);
	  memcpy( p2 , p1, (size_t)NoInd * sizeof(int) );
	  strcpy ( MarName[m_n], gv->MarName[m] ); // Marker name

	  memcpy ( &(map_S[m_n]) , &(*(gv->map[lne])), sizeof(gs_mappoint_type));
	  map[m_n] = &(map_S[m_n]); // Pointer to locus in map storage
	  mdl[m_n] = m_n;           // Index of map data line

	  m_n++;
	}
    }
  
  /* Write from temporary to final storage */
  gv->NoMar = NoMar;
  free ( gv->MarName );  gv->MarName = MarName;
  free ( gv->All1 );     gv->All1 = All1;
  free ( gv->All2 );     gv->All2 = All2;

  gv->NLoc = NLoc;
  free ( gv->map );  gv->map   = map;
  free ( gv->map_S); gv->map_S = map_S;
  free ( gv->mdl );  gv->mdl   = mdl;
  
  gv->y = y;

  sprintf ( gv->msg, "No. of matching markers: %lu", gv->NLoc );
  gs_info ( gv, 1, gv->msg );

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv; char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  /* Write map file */
  FILE *outfp;
  if (1 == *auxfiles){

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

  for (lne=0; lne <  gv->NLoc; lne++)
    if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
      {
	sprintf(gv->msg,"%i %15.12f %s %s\n",
		(gv->map[lne])->chrom, 
		(gv->map[lne])->pos  ,
		(gv->map[lne])->name ,
		(gv->map[lne])->class
		);
	fprintf(outfp,"%s",gv->msg);
      }

  fclose(outfp);

  }
 
  *retval = 0; return ;

}

void gs_read_map_01_GV ( char**  map_filename,
			 int*    skip,
			 char**  fmt,
			 double* map_stretch,
			 char**  out_filename,
			 int*    auxfiles,
			 int*    retval,
			 char**  set_name)
{

  gs_read_map_01 ( gs_fdta(*set_name),
		   map_filename,
		   skip,
		   fmt,
		   map_stretch,
		   out_filename,
		   auxfiles,
		   retval         ) ;
}

void gs_get_map_01 ( gs_varset_type *gv ,
		     char**  out_filename,
		     int*    retval        )
{
  gs_UNSIGNED  lne = 0;
  FILE *outfp;

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;
  
  for (lne=0; lne <  gv->NLoc; lne++)
    if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
      {
	sprintf(gv->msg,"%i %15.12f %s %s\n",
		(gv->map[lne])->chrom, 
		(gv->map[lne])->pos  ,
		(gv->map[lne])->name ,
		(gv->map[lne])->class
		);
	fprintf(outfp,"%s",gv->msg);
      }

  fclose(outfp);
  *retval = 0; return ;
}

void gs_get_map_01_GV (	 char**  out_filename,
			 int*    retval,
			 char**  set_name)
{

  gs_get_map_01 ( gs_fdta(*set_name),
		  out_filename,
		  retval         ) ;
}



int minimum (int a, int b) { return (a < b) ? a : b; }

void gs_str_rem_us (char* strg) 
{
  gs_UNSIGNED i, ii = strlen(strg);
  for (i = 0; i< ii; i++) 
    {
      if (
	  ('.' == strg[i]) ||
	  ('-' == strg[i]) ||
	  ('/' == strg[i]) ||
	  ('_' == strg[i]) 
	  )
	strg[i] = 'x';
    }


}

void gs_recode_ref_01 ( gs_varset_type* gv , 
			int* reference,
			int* missing,
			int* descending,
			int* retval        )
/* Recoding to a reference indidivudal on individual marker basis */
{
  gs_UNSIGNED  i, j, idx, idx_R;
  int a,b,c,d;

  if ( *reference < 1 ) {
      gs_info(gv,-2, "Invalid reference individual");
      *retval = -2; return ;
  }

  int * Ref1 =  malloc(gv->NoMar * gv->NoInd * sizeof(int));
  int * Ref2 =  malloc(gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == Ref1)  || (NULL == Ref2) ) ERR_M;
  
  gs_UNSIGNED R = (gs_UNSIGNED) *reference - 1;
      
  // Recode alleles in local memory
  for (i=0; i < gv->NoMar; i++)
    {
      idx_R = i*gv->NoInd + R;
      a = gv->All1[idx_R];
      b = gv->All2[idx_R];
      
      for (j=0; j <gv->NoInd; j++)
	{
	  idx = i*gv->NoInd + j;
	  c = gv->All1[idx];
	  d = gv->All2[idx];

	  if ( ( -1 == a ) || ( -1 == b ) || 
	       ( -1 == c ) || ( -1 == d )    ) 
	    {
	      Ref1[idx] = *missing;
	      Ref2[idx] = *missing;
	    } 
	  else
	    {
	      if (0==*descending) {
		if (( (c==a)&&(d==b)) ||
		    ( (c==b)&&(d==a))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 1;
		  } 
		else if (
			 ( (c==a)||(c==b)) ||
			 ( (d==a)||(d==b))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 2;
		  } 
		else
		  {
		    Ref1[idx] = 2;
		    Ref2[idx] = 2;
		  }
	      } else {
		if (( (c==a)||(c==b)) &&
		    ( (d==a)||(d==b))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 1;
		  } 
		else if (
			 ( (c==a)||(c==b)) ||
			 ( (d==a)||(d==b))    ) 
		  {
		    Ref1[idx] = 1;
		    Ref2[idx] = 2;
		  } 
		else
		  {
		    Ref1[idx] = 2;
		    Ref2[idx] = 2;
		  }
	      }
	    }
	}
    }
  
  // Write back the recoded alleles
  free(gv->All1); gv->All1 = Ref1;
  free(gv->All2); gv->All2 = Ref2;

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv;
  char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return;

}

void gs_recode_ref_01_GV ( int* reference,
			   int* missing,
			   int* descending,
			   int* retval,
			   char**  set_name)
{
  gs_recode_ref_01 ( gs_fdta(*set_name),
		     reference,
		     missing,
		     descending,
		     retval            ) ;

}

void gs_def_hblocks_01 ( gs_varset_type* gv , 
			 int* hap,
			 int* hap_unit,
			 char** hap_symbol,
			 char** out_filename,
			 int*   auxfiles,
			 int* retval        )
/* Determining haplotype borders */
{

  if ( 0 == gv->NLoc ) 
    {
      sprintf (gv->msg, "No linkage map was loaded");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  gs_UNSIGNED i,l,u,o,n, fib, lib;

  gv->h_u   = realloc ( gv->h_u , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  gv->h_o   = realloc ( gv->h_o , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  gv->h_n   = realloc ( gv->h_n , gv->NLoc*sizeof(gs_UNSIGNED) )  ;
  if ( (NULL == gv->h_u)||(NULL == gv->h_o)||(gv->h_n == 0) ) ERR_M;

  gv->NHap  = 0;                    
  o = 0; // upper block border
  n = 0; // number of (observed) markers
  u = 0;
  fib = gs_UNSIGNED_MAX; // first in block
  lib = gs_UNSIGNED_MAX; // last in block

  if (0 == *hap_unit) // disjunct segments
    {

      gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_UNSIGNED i,m,cc,l1,l2,idx1,idx2;
      for ( cc=0 ; cc < gv->NChrom ; cc ++ ) {
	fib = gv->chrom[cc].first_loc_i;

	for ( m=1 ; m < gv->chrom[cc].NLocC ; m++) {
	  n++;
	  for ( i=0; i < gv->NoInd; i++) {
	    l2 = gv->chrom[cc].first_loc_i + m ;
	    l1 = l2-1;
	    
	    idx1 = l1*gv->NoInd + i;
	    idx2 = l2*gv->NoInd + i;
	    
	    if (  (gv->All1[idx1] != gv->All1[idx2]) ||
		  (gv->All2[idx1] != gv->All2[idx2])    )
	      {
		lib = fib + n -1;
		gv->h_u [gv->NHap] = fib;
		gv->h_o [gv->NHap] = lib;
		gv->h_n [gv->NHap] = n;
		gv->NHap++;
		n = 0; 
		fib = lib + 1;
		break;
	      }
	  }
	}

	lib = fib + n;
	gv->h_u [gv->NHap] = fib;
	gv->h_o [gv->NHap] = lib;
	gv->h_n [gv->NHap] = n;
	gv->NHap++;
	n = 0; 
      }
    }
  else
    {
      if (1 == *hap_unit) // Marker count
	{
	  do {
	    // find blocks
	    do {
	      if ( gv->mdl[o] < gs_UNSIGNED_MAX ) {
		if  ( gs_UNSIGNED_MAX == fib )  fib = o;
		lib = o;
		n++;
	      }
	      o++; 
	    }  while  (  ( o < gv->NLoc) &&                         
			 !( (int)n > (*hap-1) ) &&
			 !( gv->map[o-1]->chrom !=  gv->map[o]->chrom ) );
	    
	    if (n>0) {
	      gv->h_u [gv->NHap] = fib;
	      gv->h_o [gv->NHap] = lib;
	      gv->h_n [gv->NHap] = n;
	      gv->NHap++;
	      
	      n = 0; 
	      fib = gs_UNSIGNED_MAX; 
	      lib = gs_UNSIGNED_MAX;
	    }
	    
	  } while  ( o < (gv->NLoc-1) ) ;
	} 
      else // Segment length
	{
	  do {
	    // find blocks
	    do {
	      if ( gv->mdl[o] < gs_UNSIGNED_MAX ) {
		if  ( gs_UNSIGNED_MAX == fib )  fib = o;
		lib = o;
		n++;
	      }
	      o++; 
	    }  while  (  ( o < gv->NLoc) &&                         
			 ! ( (double)*hap < 
			     ( gv->map[lib]->pos -  gv->map[fib]->pos ) ) &&
			 !( gv->map[o-1]->chrom !=  gv->map[o]->chrom ) );
	    
	    if (n>0) {
	      gv->h_u [gv->NHap] = fib;
	      gv->h_o [gv->NHap] = lib;
	      gv->h_n [gv->NHap] = n;
	      gv->NHap++;
	      
	      n = 0; 
	      fib = gs_UNSIGNED_MAX; 
	      lib = gs_UNSIGNED_MAX;
	    }
	  } while  ( o < (gv->NLoc-1) ) ;
	}
    }

  // Determine haploblock map
  gv->Hap_S = realloc (gv->Hap_S, gv->NHap*sizeof(gs_mappoint_type) );
  if (NULL == gv->Hap_S)  ERR_M;
  for (l = 0 ; l < gv->NHap; l++)
    {
      gv->Hap_S[l].chrom =   gv->map[ gv->h_u[l] ]->chrom ;
      gv->Hap_S[l].pos   = ( gv->map[ gv->h_o[l] ]->pos + 
			     gv->map[ gv->h_u[l] ]->pos   ) / 2;
      sprintf ( gv->Hap_S[l].name,"%1s%06lu",&(*hap_symbol[0]),l);
      sprintf ( gv->Hap_S[l].class,"b");
    }

  /* Write individual names to .hap file */
  FILE *outfp;
  char line [gs_MAXSTR];

  if (1 == *auxfiles){

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

	for ( l=0; l < gv->NHap; l++)
	  {
	    strcpy(line,"");
	    for ( i = gv->h_u[l] ; i <= gv->h_o[l] ; i++ ) 
	      {
		if ( gv->mdl[i] < gs_UNSIGNED_MAX )
		  {
		    strcat( line, gv->MarName[ gv->mdl[i] ] );
		    strcat( line, ";" );
		  }
	      }
	    sprintf( gv->msg,"%i %f %s %s %s\n",
		     gv->Hap_S[l].chrom, 
		     gv->Hap_S[l].pos,
		     gv->Hap_S[l].name ,
		     gv->Hap_S[l].class,
		     line
		     );
	    fprintf(outfp,"%s",gv->msg);
	  }
	fclose(outfp);
  }

}

void gs_def_hblocks_01_GV ( int*   hap,
			    int*   hap_unit,
			    char** hap_symbol,
			    char** out_filename,
			    int*   auxfiles,    
			    int*   retval,
			    char** set_name)
{
  gs_def_hblocks_01 ( gs_fdta(*set_name),
		      hap,
		      hap_unit,
		      hap_symbol,
		      out_filename,
		      auxfiles,
		      retval            ) ;
}



void gs_swap_block_marker (gs_varset_type* gv )
{

  gs_UNSIGNED i;

  // Marker (=haploblock) data
  gs_UNSIGNED NoInd = gv->NoInd; 
  gs_UNSIGNED NoMar = gv->NHap; 

  char (*IndName) [gs_MAXNME] = gv->IndName; gv->IndName = NULL;
  char (*MarName) [gs_MAXNME] = malloc(NoMar*sizeof(char[gs_MAXNME])); 
  for ( i = 0 ; i < gv->NHap; i++ ) 
    strcpy (MarName[i], gv->Hap_S[i].name );

  int * All1  =  gv->Ref1; gv->Ref1 = NULL;
  int * All2  =  gv->Ref2; gv->Ref2 = NULL;

  // Map data
  gs_UNSIGNED         NLoc  = gv->NHap; 
  gs_mappoint_type  * map_S = gv->Hap_S; gv->Hap_S=NULL; 
  gs_mappoint_type ** map   = malloc ( NLoc*sizeof(gs_mappoint_type*) );
  gs_UNSIGNED * mdl         = malloc ( NLoc * sizeof(gs_UNSIGNED) );

  for (i=0 ; i < NLoc ; i++ ) 
    {
      map[i] = &(map_S[i]);
      mdl[i] = i;
    }

  // Phenotypes
  double *y = gv->y ;
  gv->y = NULL ;       // prevent from being freed

  gs_reset(gv);

  gv->NoInd = NoInd; 
  gv->NoMar = NoMar; 

  gv->IndName = IndName;  
  gv->MarName = MarName; 

  gv->All1 = All1;
  gv->All2 = All2;

  gv->NLoc  = NLoc ;
  gv->map_S = map_S;
  gv->map   = map;
  gv->mdl   = mdl;

  gv->y     = y;

}

void gs_recode_hbc_01 ( gs_varset_type* gv , 
			int* reference,
			int* retval        )
/* Recoding to a haplotype blocks for BC individuals */
{

  if ( 0 == gv->NHap ) 
    {
      sprintf (gv->msg, "No haplotype blocks were defined");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  gs_UNSIGNED  l, i, j, idx, idx_R;

  // Allocate memory

  gv->Ref1 = realloc(gv->Ref1,gv->NoMar * gv->NoInd * sizeof(int));
  gv->Ref2 = realloc(gv->Ref2,gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == gv->Ref1)  || (NULL == gv->Ref2) ) ERR_M;

  // Determine hayplotypes
  gs_UNSIGNED R = (gs_UNSIGNED) *reference - 1;
  gs_UNSIGNED idx_H, n0,n1,n2;
  int a,b,c,d;

  for ( l = 0 ; l < gv->NHap ; l++ ) {

    for (j = 0; j < gv->NoInd ; j ++)
      {
	n0 = 0; n1 = 0; n2 = 0; 
	idx_H = l*gv->NoInd + j;
	
	for ( i = gv->h_u[l] ; i <= gv->h_o[l] ; i++ ) 
	  if ( gv->mdl[i] < gs_UNSIGNED_MAX )
	    {
	      
	      idx_R = i*gv->NoInd + R;
	      a = gv->All1[idx_R];
	      b = gv->All2[idx_R];
	      
	      idx = i*gv->NoInd + j;
	      c = gv->All1[idx];
	      d = gv->All2[idx];
	      
	      if ( !(                             // not missing
		     ( -1 == a ) || ( -1 == b ) || 
		     ( -1 == c ) || ( -1 == d )
						  ) ) 
		{
		  if (                          // two common alleles
		      ( (c==a)||(c==b)) &&
		      ( (d==a)||(d==b))    ) 
		    {
		      n2++ ;
		    } 
		  else if (                     // one common allele
			   ( (c==a)||(c==b)) ||
			   ( (d==a)||(d==b))    ) 
		    {
		      n1++ ;
		    } 
		  else
		    {
		      n0++;
		    }
		} 
 	    }
	
	gv->Ref1[idx_H] =  1; // Can originate from the RP
	gv->Ref2[idx_H] =  1; // 

	if (n0>0)                // At least in one marker in the block
	  {                      // both alleles can not 
	    gv->Ref1[idx_H] = 2; // originate from 
	    gv->Ref2[idx_H] = 2; // the RP
	  }  
	else if  (n1>0)  
	  {
	    gv->Ref1[idx_H] = 1;
	    gv->Ref2[idx_H] = 2;
	  }
 
      }
  }

  /* Use the haplotype data as markers */
  gs_swap_block_marker (gv);

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv; char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }
  
  *retval = 0; return;
}


void gs_recode_hbc_01_GV ( int*   reference,
			   int*   retval,
			   char** set_name)
{
  gs_recode_hbc_01 ( gs_fdta(*set_name),
		     reference,
		     retval            ) ;

}

/* qsort haplotype alleles compariosn*/
int gs_hall_cmp(const void *a, const void *b) 
{ 

  const gs_usp *ia = (const gs_usp *)a;
  const gs_usp *ib = (const gs_usp *)b;

  if ( (ia)->n < (ib)->n ) return (1);
  if ( (ia)->n > (ib)->n ) return (-1);
  return (0);
} 

void gs_recode_hil_01 ( gs_varset_type* gv , 
			int* retval        )
/* Recoding to a haplotype blocks for inbred lines */
{

  if ( 0 == gv->NHap ) 
    {
      sprintf (gv->msg, "No haplotype blocks were defined");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  gs_UNSIGNED  i, j, idx, idx_T;

  gv->Ref1 = realloc(gv->Ref1,gv->NoMar * gv->NoInd * sizeof(int));
  gv->Ref2 = realloc(gv->Ref2,gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == gv->Ref1)  || (NULL == gv->Ref2) ) ERR_M;

  // Transpose genotypic data
  int * AT1 = malloc ( gv->NoMar * gv->NoInd * sizeof(int));
  int * AT2 = malloc ( gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == AT1)  || (NULL == AT2) ) ERR_M;

  for (i=0; i < gv->NoMar; i++)
    for (j=0; j < gv->NoInd; j++) 
      {
	idx   = i*gv->NoInd + j;
	idx_T = j*gv->NoMar + i;
	AT1[idx_T] = gv->All1[idx] ;
	AT2[idx_T] = gv->All2[idx] ;
      }

  { // Parallel for each Haplotyp-Block
    gs_UNSIGNED i, j, l, k, idx_H, idx_T;
    size_t n;
    int* P1;
    int  found;

    // Pointers to Haplotypes
    int* Hall[2*gv->NoInd];
    gs_UNSIGNED NoHall;

    // Number of occurence of a haplotypes
    gs_UNSIGNED OcHall[2*gv->NoInd];

    for ( l = 0 ; l < gv->NHap ; l++ ) 
      {
	for (k = 0; k < 2*gv->NoInd; k++) OcHall[k] = 0; 
	NoHall = 0;
	// First haplo-allele is the first homologoue of the first
	// individual
	i = gv->h_u [l];            
	j = 0;
	idx_T = j*gv->NoMar + i;
	Hall [NoHall]   = &(AT1[idx_T]);
	OcHall [NoHall] ++;
	NoHall ++;
	// Size of a haplotype block
	n = (size_t) sizeof(int) * (size_t) (gv->h_n [l]);

	for (j = 0; j < gv->NoInd ; j ++)
	  {
	    // Index of the haploblock in the output matrix
	    idx_H = l*gv->NoInd + j;
	    // Adress of the beginning of the l-th block of the j-th
	    // Individual in the transposed matrix
	    idx_T = j*gv->NoMar + i;
	    P1    = &(AT1[idx_T]);

	    // Search whether the haplo-allel was already observed
	    found = 0;
	    for (k = 0; k < NoHall; k ++) 
	      {
		if ( 0 == memcmp( Hall [k] , P1 , n) )
		  {
		    gv->Ref1[idx_H] =  k; 
		    OcHall [k] ++;
		    found = 1;
		    break;
		  }
	      }
	    
	    // Haplo-Allele Not found
	    if (0 == found)
	      {
		Hall [NoHall]  = P1;
		gv->Ref1[idx_H] =  NoHall; 
		OcHall [NoHall] ++;
		NoHall++;
	      }
	    
	    // The same for the second homologue
	    P1    = &(AT2[idx_T]);
	    found = 0;
	    for (k = 0; k < NoHall; k ++) 
	      {
		if ( 0 == memcmp( Hall [k] , P1 , n) )
		  {
		    gv->Ref2[idx_H] =  k; 
		    OcHall [k] ++;
		    found = 1;
		    break;
		  }
	      }
	    if (0 == found)
	      {
		Hall [NoHall]  = P1;
		gv->Ref2[idx_H] =  NoHall; 
		OcHall [NoHall] ++;
		NoHall++;
	      }
	    
	  } // All individuals

	// Recoding the haplotypes acording to frequency
	gs_usp ar[NoHall];
	for (k = 0; k < NoHall; k ++) {
	  ar[k].a = k;
	  ar[k].n = OcHall [k];
	}
	qsort( ar, NoHall, sizeof(gs_usp), gs_hall_cmp );
	
	for (j = 0; j < gv->NoInd ; j ++)
	  {
	    idx_H = l*gv->NoInd + j;

	    for (k = 0; k < NoHall; k ++) 
	      if ( (gs_UNSIGNED)gv->Ref1[idx_H] == ar[k].a ) {
		gv->Ref1[idx_H] = 1+k; 
		break;
	      }
	    for (k = 0; k < NoHall; k ++) 
	      if ( (gs_UNSIGNED)gv->Ref2[idx_H] == ar[k].a ) {
		gv->Ref2[idx_H] = 1+k; 
		break;
	      }

	    }

      } // All haplotype blocks
  
  } // Parallel for each haplotyp-block

  free(AT1);
  free(AT2);

  /* Use the haplotype data as markers */
  gs_swap_block_marker (gv);

  /* Determine marker statistics */
  int Srv=0; int* Sretval = &Srv; char*empty="";
  char ** ind_filename = &empty; char ** mar_filename = &empty;
  char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       Sauxfiles,Sauxfiles,Sauxfiles,
		       Sauxfiles,
		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }
  
  *retval = 0; return;

}

void gs_recode_hil_01_GV ( int*   retval,
			   char** set_name)
{
  gs_recode_hil_01 ( gs_fdta(*set_name),
		     retval            ) ;
}


void gs_write_marker_data_02 ( gs_varset_type *gv ,
			       char** type ,
			       char** lpo_filename,
			       char** mmp_filename,
			       char** mpo_filename, 
			       char** ngp_filename, 
			       char** nmp_filename, 
			       char** npo_filename, 
			       char** ind_filename,
			       int*   f_ind,
			       int*   l_ind,
			       char** add_inum,
			       int*   auxfiles,
			       int*   retval )
{

  if ( NULL == gv->All1 ) 
    {
      sprintf (gv->msg, "No marker data available");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  gs_UNSIGNED lne, last_lne=0, i, j, jj, k, idx;
  char line [gs_MAXSTR];
  int m,n,a,b,c;

  /* List of indices of the individuals to be writtten out */
  /* Either: f_ind to l_ind                                */
  /* OR: listed in the file with ind_filename              */
  gs_UNSIGNED Npli = 0;             // No output indiv
  gs_UNSIGNED pli_idx[gv->NoInd] ;  // Indexes of output indiv

  /* Default: Output all individuals*/
  Npli = gv->NoInd;
  for ( jj = 0; jj < Npli; jj++) pli_idx[jj] = jj;

  /* Definition of output individuals by f_ind to l_ind */
  /* Start and stop values for output individuals */
  if ( (0 != *l_ind) || (0 != *f_ind) )
    {
      int f = *f_ind;
      int l = *l_ind;

      if ( (int)gv->NoInd < f ) f = (int)gv->NoInd ;
      if ( (int)gv->NoInd < l ) l = (int)gv->NoInd ;

      if ( 0 < f )  f--; else if ( 0 > f) f = 0;

      if (( 1 > l ) || ( (int)gv->NoInd < l) ) 
	l = (int)gv->NoInd;

      jj =0;
      for ( j=f; (int)j < l; j++ ) 
	{
	  pli_idx[jj] = j;
	  jj++;
	}
      Npli = jj;
    }

  /* Definition of output indivdiuals by names */
  /* List of output individuals */
  if (0 != strcmp(*ind_filename,""))
    {
      FILE *fp;

      if ( NULL == (fp = fopen(*ind_filename,"r")) ) ERR_F;
      j = 0;                           // Number of names read in
      while (!feof(fp))
	{
	  if ( 0 < fscanf(fp,"%s",line) )
	    {
	      if ( strlen(line) > gs_MAXSTR) {
		gs_info ( gv,-2,"Decription too long");
		*retval = -2; return ;
	      }

	      for ( k=0; k < gv->NoInd; k++ ) 
		{
		  if ( 0 == strcmp ( gv->IndName[k] , line ) ) 
		    {
		      pli_idx[j] = k;
		      j++;
		      break;
		    }
		}
	    }
	}
      fclose(fp);
      Npli = j;
    }

  FILE *outfp;
  
  int* A1 = gv->All1; 
  int* A2 = gv->All2; 

  if (1 == *auxfiles) {}

  if (NULL == gv->map)
    {
      // If no map is available then use all markers  
      gv->NLoc  = gv->NoMar;
      gv->mdl   = realloc (gv->mdl, gv->NLoc*sizeof(gs_UNSIGNED) );
      if ( (NULL == gv->mdl) ) ERR_M; 
      for ( i=0; i < gv->NoMar; i++ ) gv->mdl[i] = i;
    }

  if  ( ( 'l' == *type[0] ) || ( 'a' == *type[0] ) ) // List format output
    {
      gs_UNSIGNED m;

      if ( 'a' == *type[0] )   // Append to file
	{
	  if ( NULL == (outfp = fopen(*lpo_filename,"a")) ) ERR_M;
	}
      else                     // New file
	{
	  if ( NULL == (outfp = fopen(*lpo_filename,"w")) ) ERR_M;
	  fprintf(outfp, "I;M;A\n");
	}

      for (i=0; i<gv->NoInd; i++  )
	for ( m=0; m < gv->NoMar; m++ )
	  {
	      idx = m*gv->NoInd + i;
	      if  (-1 != gv->All1[idx])
		{
		  fprintf(outfp, "%s;%s;%i\n",
			  gv->IndName[i],
			  gv->MarName[m], 
			  gv->All1[idx] );
		  if (gv->All1[idx] != gv->All2[idx])
		    fprintf(outfp, "%s;%s;%i\n",
			    gv->IndName[i],
			    gv->MarName[m], 
			    gv->All2[idx] );
		}
	    }
      fclose(outfp);
    }

  if  ( 'm' == *type[0] )    // Matrix fomat output for ggts
    {               
      // Map data
      if (NULL != gv->map)   
	{ 
	  if ( (outfp = fopen(*mmp_filename,"w")) == NULL ) ERR_F;
	  for (lne=0; lne <  gv->NLoc; lne++){
	    if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
	      {
		sprintf( gv->msg,"%i %15.12f %s %s\n",
			 (gv->map[lne])->chrom, 
			 (gv->map[lne])->pos  ,
			 (gv->map[lne])->name ,
			 (gv->map[lne])->class
			 );
		fprintf(outfp,"%s",gv->msg);
	      }
	  }
	  fclose(outfp);
	}

      // Marker data
      if ( (outfp = fopen(*mpo_filename,"w")) == NULL ) ERR_F;
      for ( jj = 0 ; jj < Npli; jj++ )
	{
	  j = pli_idx[jj];
	  fprintf(outfp,"%s ",gv->IndName[j]);
	}
      fprintf(outfp,"\n" );
      for (lne=0; lne <  gv->NLoc; lne++){
	if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
	  {
	    i = gv->mdl[lne];
	    fprintf(outfp,"%s ",gv->MarName[i]);
	    for ( jj = 0 ; jj < Npli; jj++ )
	      {
		j = pli_idx[jj];
		idx = i*gv->NoInd + j;
		fprintf( outfp, "%i/%i ", A1[idx], A2[idx] ) ;
	      }
	    fprintf(outfp,"\n");
	  }
      }
      fclose(outfp);
      
    } // type = m

  if  ( 'n' == *type[0] ) // NTSys format for Plabsim
    {

      // Map data
      if (NULL != gv->map)
	{ // Linkage map in M
	  if ( (outfp = fopen(*nmp_filename,"w")) == NULL ) ERR_F;
	  for (lne=0; lne <  gv->NLoc; lne++){
	    if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
	      {
		sprintf( gv->msg,"%i %15.12f %s %s\n",
			 (gv->map[lne])->chrom, 
			 (gv->map[lne])->pos/100  ,
			 (gv->map[lne])->name ,
			 (gv->map[lne])->class
			 );
		fprintf(outfp,"%s",gv->msg);
	      }
	  }
	  fclose(outfp);

	  // Chromosome lengths 
	  if ( (outfp = fopen(*ngp_filename,"w")) == NULL ) ERR_F;
	  last_lne = 0;
	  for (lne=0; lne <  gv->NLoc; lne++){
	    if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
	      {
		if ( (lne>0) && 
		     ( (gv->map[lne])->chrom > (gv->map[last_lne])->chrom)) 
		  {
		    fprintf(outfp,"%9.4f\n",
			    0.01 + ( gv->map[last_lne])->pos/100);
		  }
		last_lne = lne; // last map point with marker data
	      }
	  }
	  // Last chromosome
	  fprintf(outfp,"%9.4f\n",  0.01+(gv->map[last_lne])->pos/100);
	  fclose(outfp);
	} // Map data 

      // Plabsim population type 3 format (NTSys)
      if ( (outfp = fopen(*npo_filename,"w")) == NULL ) ERR_F;
      for ( jj = 0 ; jj < Npli; jj++ )
	{
	  j = pli_idx[jj];
	  strcpy (line,gv->IndName[j]);
	  gs_str_rem_us(line);
	  fprintf(outfp,"%s%s ",line,*add_inum);
	}
      fprintf(outfp,"\n" );

      for (lne=0; lne <  gv->NLoc; lne++) {

	if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
	  {
	    i = gv->mdl[lne];
	    m = gv->NoAll[i];  // No. of alleles 
	    
	    // Write out alleles
	    for ( n=0; n < m; n++ )
	      if (-1 != gv->All[i][n])  // only non-missing
		{
		  fprintf ( outfp,"%s.%i ", gv->MarName[i],gv->All[i][n] );
		  for ( jj = 0 ; jj < Npli; jj++ )
		    {
		      j = pli_idx[jj];
		      idx = i*gv->NoInd + j;
		      a = A1[idx]; 
		      b = A2[idx];
		      c = ((gv->All[i][n]==a)||(gv->All[i][n]==b))?1:0;
		      fprintf( outfp, "%i ", c);
		    }
		  fprintf(outfp,"\n");
		}
	  }
      }
      fclose(outfp);
    }

  *retval = 0;   return ;
}

void gs_write_marker_data_02_GV ( 
				 char** type ,
				 char** lpo_filename,
				 char** mmp_filename,
				 char** mpo_filename, 
				 char** ngp_filename, 
				 char** nmp_filename, 
				 char** npo_filename, 
				 char** ind_filename,
				 int*   f_ind,
				 int*   l_ind,
				 char** add_inum,
				 int*   auxfiles,
				 int*   retval,         
				 char** set_name)
{

  gs_write_marker_data_02 ( gs_fdta(*set_name),
			    type ,          
			    lpo_filename,   
			    mmp_filename,   
			    mpo_filename,   
			    ngp_filename,   
			    nmp_filename,   
			    npo_filename,   
			    ind_filename,   
			    f_ind,          
			    l_ind,          
			    add_inum,       
			    auxfiles,       
			    retval         ) ;
}

void gs_write_map_01 ( gs_varset_type *gv  ,
		       char** mmp_filename ,
		       int*   retval        )
{

  if (NULL == gv->map)   
    {
      sprintf (gv->msg, "No linkage map was loaded");
      gs_info(gv,-2,gv->msg);
      *retval = -2; 
      return;
    }

  gs_UNSIGNED lne;
  
  FILE* outfp;

  if ( (outfp = fopen(*mmp_filename,"w")) == NULL ) ERR_F;

  for (lne=0; lne <  gv->NLoc; lne++)
    {
      if ( gv->mdl[lne] < gs_UNSIGNED_MAX )
	{
	  sprintf( gv->msg,"%i %15.12f %s %s\n",
		   (gv->map[lne])->chrom, 
		   (gv->map[lne])->pos  ,
		   (gv->map[lne])->name ,
		   (gv->map[lne])->class
		   );
	  fprintf(outfp,"%s",gv->msg);
	}
    }
  fclose(outfp);
  
  *retval = 0;   return ;
}

void gs_write_map_01_GV ( char** mmp_filename,
			  int*   retval,         
			  char** set_name      )
{

  gs_write_map_01 ( gs_fdta(*set_name),
		    mmp_filename,   
		    retval             ) ;
}


void gs_determine_alleles (int* a, int* b, char* p)
{
  const char  sep = '/'  ;
  char str [gs_MAXSTR];
  char *q1;
  char *q2;

  *a = -1;
  *b = -1;

  if ( 0 == strlen(p) ) return;

  strcpy(str,p);
  q1 = str;
 
  if ( 2 == strlen(str) )
    {
      q2 = &str[1];
    }
  else
    {
      q2  = strchr(str,sep);
      if (NULL != q2)
	{
	  q2[0] = '\0';
	  q2++;
	}
      else 
	{
	  q2=q1;
	}
    }
  
  *a = (int) atoi(q1);
  if ( 0 == (*a) )
    {
           if ( (q1[0] == 'A') || (q1[0] == 'a') )  *a = 1;
      else if ( (q1[0] == 'C') || (q1[0] == 'c') )  *a = 2;
      else if ( (q1[0] == 'G') || (q1[0] == 'g') )  *a = 3;
      else if ( (q1[0] == 'T') || (q1[0] == 't') )  *a = 4;
      else if ( (q1[0] == 'D') || (q1[0] == 'd') )  *a = 5;
      else  *a = -1;
    }

  *b = (int) atoi(q2);
  if ( 0 == (*b) )
    {
           if ( (q2[0] == 'A') || (q2[0] == 'a') )  *b = 1;
      else if ( (q2[0] == 'C') || (q2[0] == 'c') )  *b = 2;
      else if ( (q2[0] == 'G') || (q2[0] == 'g') )  *b = 3;
      else if ( (q2[0] == 'T') || (q2[0] == 't') )  *b = 4;
      else if ( (q2[0] == 'D') || (q2[0] == 'd') )  *b = 5;
      else  *b = -1;
    }  
}

void gs_remove_newline(char*p)
{
  char * pp = strchr(p,'\n');
  if ( NULL != pp ) *pp = '\0'; 
  pp = strchr(p,'\r');
  if ( NULL != pp ) *pp = '\0'; 
}

long long int getline(char **lineptr, size_t *n, FILE *stream) {

  char *p ;
  size_t offset;
  char c;
  int chunk = 32000;

  if ( ( NULL == *lineptr ) || ( 0 ==*n ) )
    {
      *lineptr = realloc ( *lineptr, chunk * sizeof(char) );
      if ( NULL == *lineptr ) return -1;
      *n = chunk;
    }

  c = fgetc(stream);
  p = *lineptr;
  while(c != EOF) 
    {
      offset = p - *lineptr;
      if ( offset > (*n - 1)) 
	{
	  *n = *n + chunk;
	  *lineptr = realloc(*lineptr, *n * sizeof(char) );
	  if (*lineptr == NULL) return -1;
	  p = *lineptr + offset;
	}
      *p++ = c;
      if (c == '\n')  break;
      c = fgetc(stream);
    }

  if ( '\r' == (*p-1) ) 
    {
      --p;
      *p = '\n';
    }  

  *p++ = '\0';

  return p - *lineptr - 1;
}

void gs_read_marker_data_02 ( gs_varset_type *gv ,
			      char** dta_filename,
			      char** ind_filename,
			      char** mar_filename,
			      char** gen_filename,
			      int*   auxfiles,
			      int*   retval        )
{

  FILE *fp;
  
  char*  line = NULL;   char** linep = &line;
  size_t n    = 0;      size_t* np   = &n;

  char *p , *q ;
  long long int c;

  gs_reset(gv);

  /* Check for skipfirst */

  gs_UNSIGNED skipfirst=0 , f=0 , s=0 ;

  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;

  getline ( linep , np , fp ) ;
  gs_remove_newline ( line ) ;
  p = line ;
  q = strtok ( line , gs_tok_del ) ;
  while (q != NULL) {
    f++ ; p = q + 1 ; q = strtok ( '\0' , gs_tok_del ) ;
  }  ;

  getline ( linep , np , fp ) ;
  gs_remove_newline ( line ) ;
  p = line ;
  q = strtok ( line , gs_tok_del ) ;
  while (q != NULL) {
    s++ ; p = q + 1 ; q = strtok ( '\0' , gs_tok_del ) ;
  }  ;

  fclose(fp);


  if (s == f) skipfirst = 1;

  /* Determine Individuals */
  gv->NoInd = 0;
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline (linep,np,fp);
  gs_remove_newline(line);
  p = strtok(line,gs_tok_del);
  if ( 1 == skipfirst) p = strtok('\0',gs_tok_del);

  while (p != NULL) 
    {
      if ( strlen(p) > (gs_MAXNME-1) )
	{
	  gs_info(gv,-2,"Individual name too long");
	  *retval = -2; return ;
	}

      /* New individual*/
      gv->NoInd++;
      gv->IndName = realloc(gv->IndName, gv->NoInd * sizeof(char[gs_MAXNME]));
      strcpy( gv->IndName[gv->NoInd-1] , p );
      /* Next token */
      p = strtok ( '\0' , gs_tok_del );
    }

  /* Determine markers */
  gv->NoMar = 0;
  while (!feof(fp))
    {
      c = getline (linep,np,fp);
      if ( c > 0 ) 
  	{
	  gs_remove_newline(line);
  	  p = strtok (line,gs_tok_del);
	  if ( NULL == p )
	    {
	      gs_info(gv,-2,"Empty line in input file");
	      *retval = -2; return ;
	    }
 	  if ( strlen(p) > (gs_MAXNME-1) )
	    {
	      gs_info(gv,-2,"Marker name too long");
	      *retval = -2; return ;
	    }
	  /* New marker*/
  	  gv->NoMar++;
  	  gv->MarName = realloc ( gv->MarName, gv->NoMar * sizeof(char[gs_MAXNME]));
  	  strcpy( gv->MarName[gv->NoMar-1] , p );
  	}
    }

  fclose(fp);

  /* Allocate memory for the genotypic data */
  gs_UNSIGNED i, j, idx;
  int a=-1 ,b=-1;  int* za=&a;  int* zb=&b;

  gv->All1 =  realloc(gv->All1, gv->NoMar * gv->NoInd * sizeof(int));
  gv->All2 =  realloc(gv->All2, gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == gv->All1)  || (NULL == gv->All2) ) ERR_M;

    /* Read in marker data */ 
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline (linep,np,fp);

  for (i=0; i < gv->NoMar; i++)
    {
      getline (linep,np,fp);             // fgets (line,32000,fp); 
      gs_remove_newline(line);
      p = strtok(line,gs_tok_del);

      for (j=0; j < gv->NoInd; j++)
	{
	  p = strtok ( '\0' , gs_tok_del );
	  gs_determine_alleles (za,zb,p);
	  idx = i*gv->NoInd + j;
	  gv->All1[idx] = a;
	  gv->All2[idx] = b;
	}
    }
  
  free (line);
  fclose(fp);

  *retval = 0;

  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       auxfiles,auxfiles,auxfiles,
		       auxfiles,
		       retval        );

   return ;

}

void gs_read_marker_data_02_GV ( char** dta_filename,
				 char** ind_filename,
				 char** mar_filename,
				 char** gen_filename,
				 int*   auxfiles,
				 int*   retval,
				 char** set_name)
{
  gs_read_marker_data_02 ( gs_fdta(*set_name),
			   dta_filename,
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   auxfiles,
			   retval         ) ;

}

void gs_read_marker_data_04 ( gs_varset_type *gv ,
			      char** dta_filename,
			      char** ind_filename,
			      char** mar_filename,
			      char** gen_filename,
			      int*   auxfiles,
			      int*   retval        )
{

  FILE *fp;
  char sep = ';' ;
  
  char*  line = NULL;   char** linep = &line;
  size_t n    = 0;      size_t* np   = &n;

  char *p, *q;
  long long int c;

  gs_reset(gv);

  /* Check for skipfirst */

  gs_UNSIGNED skipfirst=0 , f=0 , s=0 ;

  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;

  getline ( linep, np, fp ) ;
  gs_remove_newline ( line ) ;
  q = line-1 ;
  do {
    f++ ;  p = q + 1 ; q = strchr ( p , sep ) ;
  } while (q != NULL) ;

  getline ( linep, np, fp ) ;
  gs_remove_newline ( line ) ;
  q = line-1 ;
  do {
    s++ ;  p = q + 1 ;  q = strchr ( p , sep ) ;
  } while (q != NULL) ;

  fclose(fp);

  if (s == f) skipfirst = 1;

  /* Determine markers */

  gv->NoMar= 0;
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline (linep,np,fp);
  gs_remove_newline(line);
  p = strtok(line,gs_tok_del);
  if ( 1 == skipfirst) p = strtok('\0',gs_tok_del);
 
  while (p != NULL) 
    {
      if ( strlen(p) > (gs_MAXNME-1) ) {
	  gs_info(gv,-2,"Marker name too long");
	  *retval = -2; return ;
      }

      /* New marker*/
      gv->NoMar++;
      gv->MarName = realloc(gv->MarName,gv->NoMar * sizeof(char[gs_MAXNME]));
      strcpy( gv->MarName[gv->NoMar-1] , p );

      /* Next token */
      p = strtok ( '\0' , gs_tok_del );
    }

  /* Determine individuals */
  gv->NoInd = 0;
  while (!feof(fp))
    {
      c = getline (linep,np,fp);
      if ( c > 0 ) {
	  gs_remove_newline(line);
  	  p = strtok (line,gs_tok_del);
	  if ( NULL == p )
	    {
	      gs_info(gv,-2,"Empty line in input file");
	      *retval = -2; return ;
	    }
	  if ( strlen(p) > (gs_MAXNME-1) )
	    {
	      gs_info(gv,-2,"Individual name too long");
	      *retval = -2; return ;
	    }
	  /* New individual*/
	  gv->NoInd++;
	  gv->IndName = realloc ( gv->IndName, gv->NoInd * sizeof(char[gs_MAXNME]));
	  strcpy( gv->IndName[gv->NoInd-1] , p );
      }
    }

  fclose(fp);

  /* Allocate memory for the genotypic data */
  gs_UNSIGNED i, j, idx;
  int a=-1 ,b=-1;  int* za=&a;  int* zb=&b;

  gv->All1 =  realloc(gv->All1, gv->NoMar * gv->NoInd * sizeof(int));
  gv->All2 =  realloc(gv->All2, gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == gv->All1)  || (NULL == gv->All2) ) ERR_M;

  /* Read in marker data */ 
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline (linep,np,fp);

  for (j=0; j < gv->NoInd; j++)
    {
      getline (linep,np,fp);   
      gs_remove_newline(line);
      p = strchr(line,sep); p++;

      for (i=0; i < gv->NoMar; i++)
	{
	  q = strchr ( p , sep );
	  if (NULL != q) *q = '\0';
	  // Rprintf("%lu %lu ",i,j); 

	  gs_determine_alleles (za,zb,p);

	  idx = i*gv->NoInd + j;
	  gv->All1[idx] = a;
	  gv->All2[idx] = b;
	  
	  p = q+1;
	}
    }
  
  free (line);
  fclose(fp);

  *retval = 0;

  gs_marker_stats_01 ( gv ,
		       ind_filename,
		       mar_filename,
		       gen_filename,
		       auxfiles,auxfiles,auxfiles,
		       auxfiles,
		       retval        );

   return ;

}

void gs_read_marker_data_04_GV ( char** dta_filename,
				 char** ind_filename,
				 char** mar_filename,
				 char** gen_filename,
				 int*   auxfiles,
				 int*   retval,
				 char** set_name)
{
  gs_read_marker_data_04 ( gs_fdta(*set_name),
			   dta_filename,
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   auxfiles,
			   retval         ) ;

}


void gs_marker_stats_01 ( gs_varset_type *gv ,
			  char** ind_filename,
			  char** mar_filename,
			  char** gen_filename,
			  int*   ind,    
			  int*   mar,    
			  int*   gen,    
			  int*   auxfiles,
			  int*   retval        )
{

  /* Messages on number of markers and individuals */
  sprintf(gv->msg,
  	  "No. of individuals: %lu, no. of markers: %lu",
  	  gv->NoInd, gv->NoMar);
  gs_info(gv,0,gv->msg);

    gs_UNSIGNED i,j,k,m,idx,idx_all;
    int already_defined;

  /* Determine alleles  */

  gv->NoAll = realloc ( gv->NoAll, gv->NoMar* sizeof(gs_UNSIGNED) );
  gv->All   = realloc ( gv->All,   gv->NoMar* sizeof(int[gs_MAXALL]));
  if ( (NULL == gv->NoAll)  || (NULL == gv->All) ) ERR_M;

  for ( m = 0; m < gv->NoMar; m++)
    {
      gv->NoAll[m] = 0;

      for ( i = 0; i < gv->NoInd ; i++)
  	{
  	  idx_all = m*gv->NoInd + i ;

  	  already_defined = 0;
  	  for (k = 0; k < gv->NoAll[m]; k++)
  	    if ( gv->All1[idx_all] == gv->All[m][k] ) already_defined = 1;
  	  if (!already_defined){
  	    gv->All [m] [ gv->NoAll[m] ] = gv->All1[idx_all] ;
  	    gv->NoAll[m]++;
  	    if (gv->NoAll[m] == gs_MAXALL) {
  	      gs_info(gv,-2,"To many alleles");
  	      *retval = -2; return ;
  	    }
  	  }
	  
  	  already_defined = 0;
  	  for (k = 0; k < gv->NoAll[m]; k++)
  	    if ( gv->All2[idx_all] == gv->All[m][k] ) already_defined = 1;
  	  if (!already_defined){
  	    gv->All [m] [ gv->NoAll[m] ] = gv->All2[idx_all] ;
  	    gv->NoAll[m]++;
  	    if (gv->NoAll[m] == gs_MAXALL) {
  	      gs_info(gv,-2,"To many alleles per marker");
  	      *retval = -2; return ;
  	    }

  	  }
  	}
      qsort( gv->All[m], gv->NoAll[m], sizeof(int), gs_intcmp);
    }

  /* Determine allele frequencies  */
  gv->AllCnt = realloc ( gv->AllCnt, gv->NoMar*sizeof(gs_UNSIGNED[gs_MAXALL]));
  if (NULL == gv->AllCnt)  ERR_M;

  for ( m = 0; m < gv->NoMar; m++)
    for (k = 0; k < gv->NoAll[m]; k++)
      {
  	gv->AllCnt [m][k] = 0;
  	for ( i = 0; i < gv->NoInd ; i++) {
  	  idx_all = m*gv->NoInd + i ;
  	  if ( gv->All1[idx_all] == gv->All[m][k]) gv->AllCnt[m][k]++;
  	  if ( gv->All2[idx_all] == gv->All[m][k]) gv->AllCnt[m][k]++;
  	}
      }

  /* Marker statistics: missing and pic */

  gv->Pic = realloc ( gv->Pic, gv->NoMar* sizeof(double) );
  gv->Mis = realloc ( gv->Mis, gv->NoMar* sizeof(double) );
  if ( (NULL == gv->Pic)  || (NULL == gv->Mis) ) ERR_M;

  double observed;
  for ( m = 0; m < gv->NoMar; m++)
    {

     /* percentage missing per marker*/
      if ( -1 == gv->All[m][0])
  	gv->Mis[m] = ( (double) gv->AllCnt[m][0] ) / 2 / (double) gv->NoInd;
      else
  	gv->Mis[m] = 0;

      /* expected heterozygosity */
      gv->Pic[m] = 1;
      observed  = (double)( 2 * gv->NoInd ) ;
      if ( -1 == gv->All[m][0])
  	observed  -= gv->AllCnt[m][0];
      else
  	gv->Pic[m] -= pow ( (double)gv->AllCnt[m][0]/observed , 2);

      if ( observed > 0) {
  	for (k = 1; k < gv->NoAll[m]; k++)
  	  gv->Pic[m] -= pow ( (double)gv->AllCnt[m][k]/observed , 2);
      }
      else {
  	gv->Pic[m] = 0;
      }
    }

  /* Individual statistics: missing */
  gv->IndMis = realloc ( gv->IndMis, gv->NoInd*sizeof(double) );
  if (NULL == gv->IndMis)  ERR_M;

  gs_UNSIGNED count_miss ;
  for (i=0; i<gv->NoInd; i++  )
    {
      count_miss = 0;
      for (m=0; m < gv->NoMar; m++)
  	{
  	  idx = m*gv->NoInd + i;
  	  if (-1 == gv->All1[idx]) count_miss++;
  	}
      gv->IndMis[i] = (double)count_miss / (double)gv->NoMar;
    }

  FILE *outfp;

  /* Write individual names to .ind file */
  if ( (1 == *auxfiles) && (1 == *ind) ){

    if ( (outfp = fopen(*ind_filename,"w")) == NULL ) ERR_F;

    fprintf(outfp,"%s %s\n", "Name", "InMis");
    for (i=0; i < gv->NoInd; i++){
      sprintf(gv->msg,"%s %f\n",
  	      gv->IndName[i],
  	      gv->IndMis[i]);
      fprintf(outfp,"%s",gv->msg);
    }
    fclose(outfp);

  }

  /* Write marker names to .mar file */
  if ( (1 == *auxfiles) && (1 == *mar) ){

    if ( (outfp = fopen(*mar_filename,"w")) == NULL ) ERR_F;

    gs_UNSIGNED k, l, no_obs_all = 0;
    int obs_all[gs_MAXALL];
    int already_defined;

    for ( m=0; m < gv->NoMar; m++)
      for (k = 0; k < gv->NoAll[m]; k++)
  	{
  	  already_defined = 0;
  	  for (l = 0; l < no_obs_all; l++)
  	    if ( gv->All[m][k] == obs_all[l] ) already_defined = 1;
  	  if (!already_defined){
  	    obs_all[no_obs_all] = gv->All[m][k];
  	    no_obs_all++;
  	  }
  	}
    qsort( obs_all, no_obs_all, sizeof(int), gs_intcmp);
    
    fprintf(outfp,"Name NoAll MaMis ExHet " );
    for ( l=0; l<no_obs_all; l++) {
      if (-1 == obs_all[l]) fprintf(outfp,"AM ");
      else fprintf(outfp," A%i", obs_all[l]);
    }
    fprintf(outfp,"\n");
    
    for (m=0; m < gv->NoMar; m++)
      {
  	sprintf(gv->msg,"%s %lu %5.3f %5.3f ",
  		gv->MarName[m],
  		(-1 == gv->All[m][0] ? gv->NoAll[m]-1 : gv->NoAll[m] ),
  		gv->Mis[m],
  		gv->Pic[m]
  		);
  	fprintf(outfp,"%s",gv->msg);
  	for ( l=0; l<no_obs_all; l++) {
  	  gs_UNSIGNED ac = 0;
  	  for (k = 0; k < gv->NoAll[m]; k++)
  	    if ( obs_all[l] ==gv->All[m][k])  ac = gv->AllCnt[m][k];
  	  fprintf(outfp," %lu", ac );
  	}
  	fprintf(outfp,"\n");
      }
    fclose(outfp);
  }
 
  /* Write genotypes to .gen file */
  if ((1 == *auxfiles) && (1 == *gen) ) {
 
    if ( (outfp = fopen(*gen_filename,"w")) == NULL ) ERR_F;

    fprintf(outfp,"%s ","Mar/Ind");
    for (j=0; j < gv->NoInd; j++)
      fprintf(outfp,"%s ",gv->IndName[j]);
    fprintf(outfp,"\n" );
    for (i=0; i < gv->NoMar; i++)
      {
  	fprintf(outfp,"%s ",gv->MarName[i]);
  	for (j=0; j <gv->NoInd; j++){
  	  idx = i*gv->NoInd + j;
  	  fprintf(outfp,
  		  "%i/%i ",
  		  gv->All1[idx],
  		  gv->All2[idx]);
  	}
  	fprintf(outfp,"\n");
      }
	
    fclose(outfp);
    
  } 

  *retval = 0; return ;

}

void gs_marker_stats_01_GV ( 
			     char** ind_filename,
			     char** mar_filename,
			     char** gen_filename,
			     int*   ind,
			     int*   mar,
			     int*   gen,
			     int*   auxfiles,
			     int*   retval,
			     char** set_name)
{
  gs_marker_stats_01 (     gs_fdta(*set_name),
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   ind,
			   mar,
			   gen,
			   auxfiles,
			   retval         ) ;

}

void gs_read_marker_data_03 ( gs_varset_type *gv ,
			      char** dta_filename,
			      char** ind_filename,
			      char** mar_filename,
			      char** gen_filename,
			      int*   auxfiles,
			      int*   retval        )
{

  FILE *fp;

  char*  line = NULL;   char** linep = &line;
  size_t n    = 0;      size_t* np   = &n;

  char *p;
  long long int c;

  gs_reset(gv);

  /* Determine Individuals */
  gv->NoInd = 0;
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline (linep,np,fp);   
  gs_remove_newline(line);
  p = strtok(line,gs_tok_del);

  while (p != NULL) 
    {
      if ( strlen(p) > (gs_MAXNME-1) )
	{
	  gs_info(gv,-2,"Individual name too long");
	  *retval = -2; return ;
	}

      /* New individual*/
      gv->NoInd++;
      gv->IndName = realloc(gv->IndName,
			gv->NoInd * sizeof(char[gs_MAXNME]));
      strcpy( gv->IndName[gv->NoInd-1] , p );
      /* Next token */
      p = strtok ( '\0' , gs_tok_del );
    }

  /* Determine size of the data set */
  gs_UNSIGNED NoLines = 0;
  gs_UNSIGNED MarSize = 0;
  while (!feof(fp))
    {
      c = getline (linep,np,fp);
      if ( c > 0 ) 
	{
	  p = strtok (line,gs_tok_del);
	  if ( strlen(p) > MarSize) MarSize = strlen(p);
	  NoLines++;
	}
    }
  fclose(fp);

  gv->MaxML = MarSize;

  /* Messages on data size*/
  sprintf ( gv->msg, "%lu %s", NoLines, "Datalines" );
  gs_info ( gv, 1, gv->msg );
  sprintf (gv->msg, "Description length: Mar %lu",  MarSize);
  gs_info (gv,1,gv->msg);

  /* Read in marker names */

  gs_UNSIGNED i,j,m,idx;
  
  static char (*mar) [1+MarSize] = NULL;
  mar = realloc ( mar, NoLines*sizeof(char[1+MarSize]) );
  if (NULL == mar) ERR_M;

  /* Read in marker names */
 
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline (linep,np,fp);
  for (j = 0; j < NoLines; j++)
    {
      fgets(line,gs_MAXLNE-1,fp);
      p = strtok (line,".");
      strcpy(mar[j],p);
    }

  fclose(fp);

  /* Marker names */

  gv->NoMar  = 1;
  gv->MarName = realloc ( gv->MarName, 1*sizeof(char[gs_MAXNME]) );
  strcpy ( gv->MarName[0], mar[0] );

  for (j = 1; j < NoLines; j++)
    {
      if ( 0 != strcmp ( mar[j] , mar[j-1] ) )
  	{
  	  gv->NoMar++;
  	  gv->MarName = realloc(gv->MarName,
  				gv->NoMar * sizeof(char[gs_MAXNME]));
  	  strcpy( gv->MarName[gv->NoMar-1] , mar[j] );
      }
    }

  /* Allocate memory for the genotypic data */
  
  gv->All1 =  realloc(gv->All1, gv->NoMar * gv->NoInd * sizeof(int));
  gv->All2 =  realloc(gv->All2, gv->NoMar * gv->NoInd * sizeof(int));
  if ( (NULL == gv->All1)  || (NULL == gv->All2) ) ERR_M;

  for (m=0; m < gv->NoMar; m++)
    for (i=0; i < gv->NoInd; i++)
  	{
  	  idx = m*gv->NoInd + i;
    	  gv->All1[idx] = -1;
  	  gv->All2[idx] = -1;
  	}

  /* Read in marker data */ 
  int conv,al;
  if ( NULL == (fp = fopen(*dta_filename,"r")) ) ERR_F;
  getline (linep,np,fp);

  for (j = 0; j < NoLines; j++)
    {
      getline (linep,np,fp); 
      p  = strtok (line,".");
      al = atoi (strtok (NULL," ")); 
      for ( m=0; m<gv->NoMar; m++) 
	if (0 == strcmp(gv->MarName[m],p)) break;
      for (i=0; i < gv->NoInd; i++)
	{
	  conv = atoi ( strtok (NULL," "));
	  if (1==conv) 
	    {
	      idx = m*gv->NoInd + i;
	      if ((-1 == gv->All1[idx])&&(-1 == gv->All2[idx]))
		{
		  gv->All1[idx] = (int)al;
		  gv->All2[idx] = (int)al;
		}
	      else
		{
		  gv->All2[idx] = (int)al;
		}
	    }
	}

    }

  free (line);
  fclose(fp);
  
  int Srv; int* Sretval = &Srv;
  gs_marker_stats_01 ( gv ,
  		       ind_filename,
  		       mar_filename,
  		       gen_filename,
		       auxfiles,auxfiles,auxfiles,
		       auxfiles,
  		       Sretval        );
  if (-2 == *Sretval) { *retval = -2; return; }

  *retval = 0; return ;

}

void gs_read_marker_data_03_GV ( char** dta_filename,
				 char** ind_filename,
				 char** mar_filename,
				 char** gen_filename,
				 int*   auxfiles,
				 int*   retval,
				 char** set_name)
{
  gs_read_marker_data_03 ( gs_fdta(*set_name),
			   dta_filename,
			   ind_filename,
			   mar_filename,
			   gen_filename,
			   auxfiles,
			   retval         ) ;

}

double st_time_store[5];

void st_start_timer (int* i, double* time)
{
  st_time_store[*i] = *time;
}

void st_stop_timer (int* i, double* time)
{
  *time = st_time_store[*i];
}

int st_info_level = 0;

extern int INFOLEVEL;

void st_set_info_level (int* i)
{
  st_info_level = *i;          // Selection Tools
  INFOLEVEL     = *i;          // Plabsim
  gs_set_all_info_levels (i); // All data sets
}

void st_get_info_level (int* i)
{
  *i = st_info_level;
}

void st_info ( int level, 
	       char* message )
{
    if ( st_info_level >= level )
    {
        if ( 1 <= level )        Rprintf( "I: %s\n", message );
        else if (  0 == level )  Rprintf( "M: %s\n", message );
        else if ( -1 == level )  Rprintf( "W: %s\n", message );
        else if ( -2 == level )  Rprintf( "E: %s\n", message );
    }
}

void st_r_info(int* level, char** message)
{
  st_info( *level, *message );
}

void gs_cross_init_01 ( gs_varset_type *gv  ,
			int*   retval         )
{
  if ( (NULL == gv->crs) || (NULL == gv->crs_S) )
    {
      gv->NCrs  = gv->NoInd * (gv->NoInd-1) / 2;
      gv->crs_S = realloc ( gv->crs_S, gv->NCrs * sizeof( gs_cross_type ) );
      gv->crs   = realloc ( gv->crs,   gv->NCrs * sizeof( gs_cross_type* ) );
      if ( (NULL == gv->crs) || (NULL == gv->crs_S) ) ERR_M;

      gs_UNSIGNED i,j,k;
      k = 0;
      
      for ( i=0; i < (gv->NoInd-1) ; i++ )  
	for ( j = i+1; j < gv->NoInd ; j++ )
	  {
	    (gv->crs_S[k]).p1 = i;
	    (gv->crs_S[k]).p2 = j;
	    gv->crs[k] = &(gv->crs_S[k]);
	    k++;
	  }
      
      for ( k = 0; k <  gv->NCrs ; k++)
	{
	  gv->crs[k]->gd = 0;
	  gv->crs[k]->mu = 0;
	  gv->crs[k]->mi = 0;
	  gv->crs[k]->ma = 0;
	  gv->crs[k]->va = 0;
	  gv->crs[k]->es = 0;
	}
      
      *retval = 0; return;
    }
}


void gs_calc_rd ( gs_varset_type *gv ,
		  gs_UNSIGNED     k   )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_UNSIGNED  n = 0; // non missing markers
  gs_FLOAT    di = 0; // sum of distances

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];
      
      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
      
      if ( !(  ( -1 == a ) || ( -1 == b ) || 
	       ( -1 == c ) || ( -1 == d )    ) ) 
	{
	  n++;
	  
	  if ( !( ( (c==a)&&(d==b)) ||
		  ( (c==b)&&(d==a))    )  ) 
	    {
	      if ( ( (c==a)||(c==b)) ||
		   ( (d==a)||(d==b))    )  
		{
		  di += 0.5;
		}
	      else 
		{
		  di += 1;
		}
	    }
	}
    }

  (gv->crs_S[k]).gd = (gs_FLOAT) di / (gs_FLOAT)n;
  
}

void gs_calc_mrd ( gs_varset_type *gv ,
	           gs_UNSIGNED     k   )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_UNSIGNED  n = 0; // non missing markers
  gs_FLOAT    di = 0; // sum of distances

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];
      
      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
      
      if ( !(  ( -1 == a ) || ( -1 == b ) || 
	       ( -1 == c ) || ( -1 == d )    ) ) 
	{
	  n++;
	  
	  if ( !( ( (c==a)&&(d==b)) ||
		  ( (c==b)&&(d==a))    )  ) 
	    {
	      if ( ( (c==a)||(c==b)) ||
		   ( (d==a)||(d==b))    )  
		{
		  di += 0.5;
		}
	      else 
		{
		  di += 2;
		}
	    }
	}
    }

  (gv->crs_S[k]).gd = (gs_FLOAT) sqrt( di / (gs_FLOAT)n / 2 );
  
}

void gs_calc_euc ( gs_varset_type *gv ,
	           gs_UNSIGNED     k   )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_UNSIGNED  n = 0; // non missing markers
  gs_FLOAT    di = 0; // sum of distances

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];
      
      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
      
      if ( !(  ( -1 == a ) || ( -1 == b ) || 
	       ( -1 == c ) || ( -1 == d )    ) ) 
	{
	  n++;
	  
	  if ( !( ( (c==a)&&(d==b)) ||
		  ( (c==b)&&(d==a))    )  ) 
	    {
	      if ( ( (c==a)||(c==b)) ||
		   ( (d==a)||(d==b))    )  
		{
		  di += 0.5;
		}
	      else 
		{
		  di += 2;
		}
	    }
	}
    }

  (gv->crs_S[k]).gd = (gs_FLOAT) sqrt( di );
  
}

void gs_calc_mu ( gs_varset_type *gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT   *mu = &((gv->crs_S[k]).mu); 

  *mu = gv->u[0];
 
  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);
      *mu  +=  (gs_FLOAT)count / (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1+m]);

    }
}

void gs_cross_eval_mu_01 ( gs_varset_type* gv  ,
			   int* retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

  #pragma omp parallel for  
  for ( k = 0; k <  gv->NCrs ; k++)
    {
      gs_calc_mu (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_mu_01_GV ( int*    retval  ,
		              char**  set_name   )
{
  gs_cross_eval_mu_01 ( gs_fdta(*set_name) ,
		        retval              );
}

void gs_calc_ma ( gs_varset_type *gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT   *ma = &((gv->crs_S[k]).ma); 

  *ma = gv->u[0];

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);

      if ( ( 0 < count ) && (0 < gv->u[1+m]) )
	*ma  +=  (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1+m]);
    }
}

void gs_cross_eval_ma_01 ( gs_varset_type* gv  ,
			   int* retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

  #pragma omp parallel for  
  for ( k = 0; k <  gv->NCrs ; k++)
    {
      gs_calc_ma (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_ma_01_GV ( int*    retval  ,
		              char**  set_name   )
{
  gs_cross_eval_ma_01 ( gs_fdta(*set_name) ,
		        retval              );
}
void gs_calc_mi ( gs_varset_type *gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT   *mi = &((gv->crs_S[k]).mi); 

  *mi = gv->u[0];

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);

      if ( ( 0 < count ) && (0 > gv->u[1+m]) )
	*mi  +=  (gs_FLOAT)2 * (gs_FLOAT)(gv->u[1+m]);
    }
}

void gs_cross_eval_mi_01 ( gs_varset_type* gv  ,
			   int* retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

  #pragma omp parallel for  
  for ( k = 0; k <  gv->NCrs ; k++)
    {
      gs_calc_mi (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_mi_01_GV ( int*    retval  ,
		              char**  set_name   )
{
  gs_cross_eval_mi_01 ( gs_fdta(*set_name) ,
		        retval              );
}

void gs_cross_eval_gd_01 ( gs_varset_type *gv  ,
			   char ** dist ,
			   int*   retval         )
{

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available, run 'read.marker.data' first");
    *retval = -2; return ;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  gs_UNSIGNED k;

  if ( 0 == strcmp("mrd",*dist) ) 
    {
      #pragma omp parallel for  
      for ( k = 0; k <  gv->NCrs ; k++)  gs_calc_mrd (gv,k);
    }
  else if ( 0 == strcmp("rd",*dist) ) 
    {
      #pragma omp parallel for  
      for ( k = 0; k <  gv->NCrs ; k++)  gs_calc_rd (gv,k);
    }
  else
    {
      #pragma omp parallel for  
      for ( k = 0; k <  gv->NCrs ; k++)  gs_calc_euc (gv,k);
    }

  *retval = 0; return;
}

void gs_cross_eval_gd_01_GV ( char**  dist    ,
			      int*    retval  ,
			      char**  set_name   )
{
  gs_cross_eval_gd_01 ( gs_fdta(*set_name) ,
			dist,
			retval              );
}

void gs_calc_genv_01 (gs_varset_type* gv)

{
  gs_UNSIGNED ii, jj, idx;
  int a,b, e, count;

  gv->genv = realloc (gv->genv, gv->NoMar * gv->NoInd * sizeof (gs_FLOAT));

  for (ii=0; ii < gv->NoMar; ii++)
    {
      e = gv->EffAll[ 1+ii ];

      for (jj=0; jj < gv->NoInd; jj++)
	{
	  idx = ii*gv->NoInd + jj;
	  a = gv->All1[idx] ;
	  b = gv->All2[idx] ;
	  count = (int)(e==a) + (int)(e==b);
	  gv->genv[idx] = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1+ii]) /2;
	}
    }
}


void gs_calc_va_linked ( gs_varset_type* gv,
			 gs_UNSIGNED k      )
{
  gs_UNSIGNED c, ii, jj, m, idx, idx_i, idx_j;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT q, ggu, ggv, hhu, hhv;

  gs_FLOAT* va = &((gv->crs_S[k]).va); 

  *va = 0;

  for ( c=0; c < gv->NChrom; c++ )
    {
       for  ( ii = 0 ; ii < gv->chrom[c].NLocC ; ii ++)
	 for  ( jj = 0 ; jj < gv->chrom[c].NLocC ; jj ++)
	   {
	     idx = ii * gv->chrom[c].NLocC + jj;
	     q = (gv->chrom[c].q) [idx]  ;

	     m = ii + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggu   = gv->genv[idx_i];
	     hhu   = gv->genv[idx_j];  

	     m = jj + gv->chrom[c].first_loc_i;
	     idx_i = m*gv->NoInd + i;
	     idx_j = m*gv->NoInd + j;
	     ggv   = gv->genv[idx_i];
	     hhv   = gv->genv[idx_j];  
	     
	     *va += (q/2 - 0.25) * (ggu*ggv + hhu*hhv - ggu*hhv - hhu*ggv);
	   }
    }

  *va = *va * 4;

}


void gs_calc_va ( gs_varset_type* gv,
		  gs_UNSIGNED k )
{
  gs_UNSIGNED m, idx_i, idx_j;
  int a,b,c,d, e, count;

  gs_UNSIGNED i = (gv->crs_S[k]).p1 ;
  gs_UNSIGNED j = (gv->crs_S[k]).p2 ;  

  gs_FLOAT  EZ2 = 0;
  gs_FLOAT  EZ  = 0;

  gs_FLOAT* va = &((gv->crs_S[k]).va); 

  *va = 0;

  for (m=0; m < gv->NoMar; m++)
    {
      idx_i = m*gv->NoInd + i;
      a = gv->All1[idx_i];
      b = gv->All2[idx_i];

      idx_j = m*gv->NoInd + j;
      c = gv->All1[idx_j];
      d = gv->All2[idx_j];
	  
      e = gv->EffAll[ 1+m ];

      count = (int)(e==a) + (int)(e==b) + (int)(e==c) + (int)(e==d);

      EZ2  = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1+m]) * (gs_FLOAT)(gv->u[1+m]) / 4;
      EZ   = (gs_FLOAT)count * (gs_FLOAT)(gv->u[1+m])  / 4 ;

      *va += ( EZ2 - EZ*EZ ) * 4; // *4 both gametes = one entrire dh line

    }
}

void gs_cross_eval_va_01 ( gs_varset_type* gv  ,
			   char** pop_type     ,
			   int*   t            ,
			   char** map_function ,
			   int*   retval          )
{

  if ( NULL == gv->u ) {
    gs_info(gv,-2,"No estimated effects available");
    *retval = -2; return;
  }

  int Srv=0; int* Sretval = &Srv;
  gs_cross_init_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return ; }

  char*empty=""; char ** Sout_filename = &empty;
  int auf = 0; int*Sauxfiles = &auf;

  gs_UNSIGNED k;

  if ( 0 == strcmp("unlinked",*pop_type) )
    {
       #pragma omp parallel for  
       for ( k = 0; k <  gv->NCrs ; k++)
	 {
	   gs_calc_va (gv,k);
	 }
    }

  else if ( ( 0 == strcmp(*pop_type,"DH")) ||
	    ( 0 == strcmp(*pop_type,"SSD"))  )
    {

      gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_rf_01 ( gv, map_function, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_q_01  ( gv, pop_type, t, Sout_filename, Sauxfiles, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_calc_genv_01 (gv);

       #pragma omp parallel for  
       for ( k = 0; k <  gv->NCrs ; k++)
	 {
	   gs_calc_va_linked (gv,k);
	 }
    }

  else {
    gs_info(gv,-2,"Unknown population type");
    *retval = -2; return;
  }


  *retval = 0; return;
}

void gs_cross_eval_va_01_GV ( char**  pop_type    ,
			      int*    t           ,
			      char**  map_function ,
			      int*    retval      ,
			      char**  set_name      )
{
  gs_cross_eval_va_01 ( gs_fdta(*set_name) ,
			pop_type           ,
			t                  ,
			map_function       ,
		        retval               );
}

void gs_calc_es ( gs_varset_type* gv ,
		  gs_UNSIGNED k      ,  
		  double* i            )
{
  gs_FLOAT* mu = &((gv->crs_S[k]).mu); 
  gs_FLOAT* va = &((gv->crs_S[k]).va); 
  gs_FLOAT* es = &((gv->crs_S[k]).es); 

  *es = *mu + *i * sqrt(*va) ;

}

void gs_cross_eval_es_01 ( gs_varset_type* gv  ,
			   double * i          ,
			   int* retval          )
{
  if ( (NULL == gv->crs) || (NULL == gv->crs_S) 
       || ( 0 == gv->crs[0]->va ) ) {
    gs_info(gv,-2,"No estimated variances available");
    *retval = -2; return;
  }

  sprintf( gv->msg, "Standardized selection differential %5.2f", (float)*i );
  gs_info(gv,0,gv->msg);

  gs_UNSIGNED k;

  #pragma omp parallel for  
  for ( k = 0; k <  gv->NCrs ; k++)
    {
      gs_calc_es (gv,k,i);
    }

  *retval = 0; return;
}

void gs_cross_eval_es_01_GV ( double* i       ,
			      int*    retval  ,
			      char**  set_name   )
{
  gs_cross_eval_es_01 ( gs_fdta(*set_name) ,
			i                  ,
		        retval               );
}

int gs_cross_cmp_gd(const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->gd < (*ib)->gd )  return 1;
  if ( (*ia)->gd > (*ib)->gd )  return -1;
  return 0;

} 

int gs_cross_cmp_mu (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->mu < (*ib)->mu )  return 1;
  if ( (*ia)->mu > (*ib)->mu )  return -1;
  return 0;

} 

int gs_cross_cmp_ma (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->ma < (*ib)->ma )  return 1;
  if ( (*ia)->ma > (*ib)->ma )  return -1;
  return 0;
} 

int gs_cross_cmp_mi (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->ma < (*ib)->ma )  return -1;
  if ( (*ia)->ma > (*ib)->ma )  return 1;
  return 0;
} 

int gs_cross_cmp_va (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->va < (*ib)->va )  return 1;
  if ( (*ia)->va > (*ib)->va )  return -1;
  return 0;

} 

int gs_cross_cmp_es (const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->es < (*ib)->es )  return 1;
  if ( (*ia)->es > (*ib)->es )  return -1;
  return 0;

} 

int gs_cross_cmp_index(const void *a, const void *b) 
{ 

  const gs_cross_type **ia = (const gs_cross_type **)a;
  const gs_cross_type **ib = (const gs_cross_type **)b;

  if ( (*ia)->p1 > (*ib)->p1 )  return 1;
  else {
    if ( (*ia)->p1 < (*ib)->p1 )  return -1;
    else {
      if ( (*ia)->p2 > (*ib)->p2 )  return 1;
      else {
	if ( (*ia)->p1 < (*ib)->p1 )  return -1;
	else return 0;
      }
    }
  }

} 



void gs_cross_info_01 (  gs_varset_type *gv   ,
			 int   * bestn        ,
		         char ** sortby       ,
			 char**  out_filename ,
			 int*    retval         )
{
  gs_UNSIGNED i,m;
  FILE *outfp;
   
  if ( NULL == gv->crs ) {
      gs_info(gv,-2, "No cross information available");
      *retval = -2; return ;
  }

  if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

  if ( 0 == strcmp(*sortby,"index") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_index           );

  if ( 0 == strcmp(*sortby,"gd") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_gd           );

  if ( 0 == strcmp(*sortby,"mu") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_mu           );

  if ( 0 == strcmp(*sortby,"mi") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_mi           );

  if ( 0 == strcmp(*sortby,"ma") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_ma           );

  if ( 0 == strcmp(*sortby,"va") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_va           );

  if ( 0 == strcmp(*sortby,"es") )
    qsort ( gv->crs, 
	    gv->NCrs, 
	    sizeof (gs_cross_type *), 
	    gs_cross_cmp_es           );

  if (0 == *bestn) 
    m = gv->NCrs;
  else
    m = ( (gs_UNSIGNED)(*bestn) < gv->NCrs) ? (gs_UNSIGNED)(*bestn) : gv->NCrs;

  fprintf(outfp,"P1No P2No P1Name P2Name gd mu mi ma va es\n");
  for ( i = 0; i < m ; i++)
    fprintf(outfp,"%lu %lu %s %s %f %f %f %f %f %f\n",
	    1+(gv->crs[i])->p1,
	    1+(gv->crs[i])->p2,
	    gv->IndName[ (gv->crs[i])->p1 ],
	    gv->IndName[ (gv->crs[i])->p2 ],
	    (double)(gv->crs[i])->gd,
	    (double)(gv->crs[i])->mu,
	    (double)(gv->crs[i])->mi,
	    (double)(gv->crs[i])->ma,
	    (double)(gv->crs[i])->va,
	    (double)(gv->crs[i])->es
	    );
  
  fclose(outfp);  

  *retval = 0; 
  return ;
}

void gs_cross_info_01_GV (int   * bestn,
			  char ** sortby,
			  char**  out_filename,
			  int*    retval,
			  char**  set_name         )
{
  gs_cross_info_01 ( gs_fdta(*set_name),
		     bestn,
		     sortby,
		     out_filename,
		     retval         ) ;
}


extern void ps_get_simpop ( char ** pop_name   ,
			    gs_UNSIGNED * NoInd,
			    gs_UNSIGNED * NoMar,
			    char (**IndName) [gs_MAXNME] ,
			    char (**MarName) [gs_MAXNME] ,
			    int *(* All1)    ,
			    int *(* All2)    ,
			    gs_UNSIGNED       *  NLoc ,
			    gs_mappoint_type  ** map_S,
			    gs_mappoint_type *** map  ,
			    gs_UNSIGNED       ** mdl  ,
			    int*    bg,
			    int* Sretval                  );

void st_get_simpop (  gs_varset_type *gv   ,
	              char ** pop_name,
		      int*    bg,
		      int*    retval         )
{
  gs_reset (gv);

  gs_UNSIGNED * NoInd = &(gv->NoInd);
  gs_UNSIGNED * NoMar = &(gv->NoMar);

  char   (**IndName) [gs_MAXNME]   = &(gv->IndName); 
  char   (**MarName) [gs_MAXNME]   = &(gv->MarName); 

  int *(* All1)  = &(gv->All1); 
  int *(* All2)  = &(gv->All2); 

  gs_UNSIGNED       *  NLoc  =  &(gv->NLoc);   
  gs_mappoint_type  ** map_S =  &(gv->map_S);  
  gs_mappoint_type *** map   =  &(gv->map);    
  gs_UNSIGNED       ** mdl   =  &(gv->mdl);      

 int Srv=0; int* Sretval = &Srv;
 ps_get_simpop ( pop_name,
		 NoInd,
		 NoMar,
		 IndName,
		 MarName,
		 All1,
		 All2,
                 NLoc, 
		 map_S,
		 map,  
		 mdl,  
		 bg,
		 Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }

  /* Determine marker statistics */
  {
    int Srv=0; int* Sretval = &Srv; char*empty="";
    char ** ind_filename = &empty; char ** mar_filename = &empty;
    char ** gen_filename = &empty; int auf = 0; int*Sauxfiles = &auf;
    gs_marker_stats_01 ( gv ,
			 ind_filename,
			 mar_filename,
			 gen_filename,
			 Sauxfiles,Sauxfiles,Sauxfiles,
			 Sauxfiles,
			 Sretval        );
    if (-2 == *Sretval) { *retval = -2; return; }
  }


  *retval = 1;   return;
}

void st_get_simpop_GV (  char ** pop_name,
			 int*    bg,
			 int*    retval,
			 char** set_name         )
{
  st_get_simpop ( gs_fdta(*set_name),
		  pop_name,
		  bg,
		  retval   );     
}

void il_get_population_size (  gs_varset_type *gv ,
			       int*   NoLines       )
{
  *NoLines = (int) gv->NoInd;
}
				
void il_get_population_size_GV ( int*   NoLines,
				 char** set_name         )
{
  il_get_population_size ( gs_fdta(*set_name),
			   NoLines           ) ;
}


void il_calc_segment_length ( gs_varset_type *gv )
{
  gs_UNSIGNED l,c,m;
  gs_UNSIGNED i1, i2, i3;
  gs_FLOAT    p1, p2, p3;

  gv->dimu = gv->NLoc;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );

  for ( c=0,l=0 ; c < gv->NChrom ; c ++ ) 
    {
      m=0;
      i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
      i3 =  gv->chrom[c].first_loc_i + m + 1 ; p3 =  gv->map[i3]->pos;
      gv->u[l] = (p3-p2)/2 ; 
      l++;

      for ( m=1 ; m < (-1 + gv->chrom[c].NLocC) ; m++,l++)
	{
	  i1 =  gv->chrom[c].first_loc_i + m - 1 ; p1 =  gv->map[i1]->pos;
	  i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
	  i3 =  gv->chrom[c].first_loc_i + m + 1 ; p3 =  gv->map[i3]->pos;
	  gv->u[l] = (p2-p1)/2 + (p3-p2)/2; 
	}

      i1 =  gv->chrom[c].first_loc_i + m - 1 ; p1 =  gv->map[i1]->pos;
      i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
      gv->u[l] = (p2-p1)/2 ; 
      l++;
    }

}

void il_mark_loci  ( gs_varset_type *gv ,
		     int*    chrom      ,
		     double* lower      ,
		     double* upper      ,
		     int*    exclude       )
{

  gs_UNSIGNED l,c,m;
  gs_UNSIGNED i2;
  gs_FLOAT    p2;

  int hit = 0;

  if ( 0 == *chrom ) return ;

  for ( c=0,l=0 ; c < gv->NChrom ; c ++ ) 
    {
      for ( m=0 ; m <  gv->chrom[c].NLocC ; m++,l++)
	{
	  i2 =  gv->chrom[c].first_loc_i + m     ; p2 =  gv->map[i2]->pos;
	  
	  hit = ( ( c == (gs_UNSIGNED) (-1 + *chrom )) &&
		  ( (0 == *upper) || 
		    ( ( *lower <= p2) && ( p2 <= *upper ) ) ) );

	  if (1 == *exclude) hit = !hit;
		  
	  if ( !hit ) gv->u[l] = 0;
	}
    }
}


double il_calc_tot (gs_varset_type *gv)
{
  gs_UNSIGNED l,m;
  gs_FLOAT    p2;

  int cc, s_c, e_c;
  s_c = 0;
  e_c = (int)gv->NChrom;

  double tot = 0;

  for ( cc = s_c ; cc < e_c ; cc ++ ) 
    {
      for ( m=0 ; m < gv->chrom[cc].NLocC ; m++)
	{
	  l =  gv->chrom[cc].first_loc_i + m ;
	  p2 = gv->map[l]->pos;
	  tot += gv->u[l];
	}
    }
  tot = tot * 2;
  return (tot);
}

void il_eval_lines  (  gs_varset_type *gv ,
		       int*   allele      ,
		       int*   chrom       ,
		       double* lower      ,
		       double* upper      ,
		       int*   exclude     ,
		       double* abs        ,
		       double* tot        ,
		       double* nseg       ,
		       double* lseg       ,
		       int*   retval        )
{
  gs_UNSIGNED i,l,m,idx;
  gs_UNSIGNED l1,l2,idx1,idx2;
  gs_FLOAT    p2;

  int cc, s_c, e_c;

  if ( ( 0 == gv->NoInd) || ( 0 == gv->NoMar) ) {
    gs_info(gv,-2,"No marker data available");
    *retval = -2; return ;
  }

  gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  s_c = 0;
  e_c = (int)gv->NChrom;

  if ( ( *chrom <0 ) || ( *chrom > (int)gv->NChrom )) {
    gs_info(gv,-2,"Incorrect chromosome number");
    *retval = -2; 
    return; 
    }

  for ( i=0 ; i < gv->NoInd; i ++ ) 
    {
      
      il_calc_segment_length ( gv );

      il_mark_loci  ( gv           ,
		      &(chrom[i])  ,
		      &(lower[i])  ,
		      &(upper[i])  ,
		      exclude  );

      tot[i] = il_calc_tot (gv);

      abs[i] = 0;

      // Allele frequencies
      for ( cc = s_c ; cc < e_c ; cc ++ ) 
	{
	  for ( m=0 ; m < gv->chrom[cc].NLocC ; m++)
	    {
	      l =  gv->chrom[cc].first_loc_i + m ;
	      p2 = gv->map[l]->pos;
	      idx = l*gv->NoInd + i;
	      abs[i] +=
		(gv->All1[idx] == *allele) * gv->u[l]  +
		(gv->All2[idx] == *allele) * gv->u[l];

	    }
	}

      // Segments
      nseg[i] = 0;
      for ( cc = s_c ; cc < e_c ; cc ++ ) 
	{
	  m = 0;
	  l2 = gv->chrom[cc].first_loc_i + m ;
	  idx2 = l2*gv->NoInd + i;
	  if  (0 != gv->u[l2]) 
	    {
	       if  (*allele == gv->All1[idx2]) nseg[i] ++ ;
	       if  (*allele == gv->All2[idx2]) nseg[i] ++ ;
	    }

	  for ( m=1 ; m < gv->chrom[cc].NLocC ; m++)
	    {
	      l2 = gv->chrom[cc].first_loc_i + m ;
	      l1 = l2-1;

	      idx1 = l1*gv->NoInd + i;
	      idx2 = l2*gv->NoInd + i;

	      if  ( (0 != gv->u[l1]) && ((0 != gv->u[l2])))
		{
		  if  ( (gv->All1[idx1] != gv->All1[idx2])  &&
			(*allele == gv->All1[idx2])             )
		    nseg[i] ++ ;

		  if  ( (gv->All2[idx1] != gv->All2[idx2])  && 
			(*allele == gv->All2[idx2])             )
		    nseg[i] ++ ;
		} 
	      else if ( (0 == gv->u[l1]) && ((0 != gv->u[l2])))
		{
		  if  (*allele == gv->All1[idx2]) nseg[i] ++ ;
		  if  (*allele == gv->All2[idx2]) nseg[i] ++ ;
		}
	    }
	}
      if ( 0 != nseg[i] )
	lseg[i] = abs[i] / nseg[i];

    } // for i = 0 .. ind

  *retval=0;
}
				
void il_eval_lines_GV ( int*   allele  ,
			int*   chrom   ,
			double* lower  ,
			double* upper  ,
			int* exclude   ,
			double* abs    ,
			double* tot    , 
			double* nseg   ,
			double* lseg   , 
			int*   retval  ,
			char** set_name  )
{
  il_eval_lines ( gs_fdta(*set_name),
		  allele            ,
		  chrom             ,
		  lower             ,
		  upper             ,
		  exclude           ,
		  abs               ,
		  tot               ,
		  nseg              ,
		  lseg              ,
		  retval              ) ;
}

void il_eval_library (  gs_varset_type *gv ,
			double* cov    ,
			double* dep    ,
			double* seglib  ,
			int   * allele ,
			int*    retval       )
{
  gs_UNSIGNED i,j,m,cc,l1,l2,idx,idx1,idx2,cnt;
  double tot;

  gs_chrom_stats_01 ( gv, Sout_filename, Sauxfiles, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  il_calc_segment_length ( gv );

  /* Genome coverage */

  tot = 0;
  *cov = 0;
  for (i=0; i < gv->NoMar; i++)
    {
      tot += gv->u[i] ;
      for (j=0; j < gv->NoInd; j++)
	{
	  idx = i*gv->NoInd + j;
	  if ( (*allele == gv->All1[idx] ) ||
	       (*allele == gv->All2[idx] )    )
	    {
	      *cov += gv->u[i];
	      break;
	    }
	}
      }
  *cov = *cov/tot;

  /* Depth of coverage */

  gs_UNSIGNED *depth =  malloc( gv->NoMar * sizeof(gs_UNSIGNED) );
  if ( NULL == depth ) ERR_M;

  for (i=0; i < gv->NoMar; i++)
    {
      depth[i] = 0;

      for (j=0; j < gv->NoInd; j++)
	{
	  idx = i*gv->NoInd + j;
	  if ( (*allele == gv->All1[idx] ) ||
	       (*allele == gv->All2[idx] )    )
	    {
	      depth[i] ++;
	    }
	}
      }

  *dep = 0;
  cnt  = 0;
  for (i=0; i < gv->NoMar; i++)
    {
      if ( depth[i] > 0 )
	{
	  *dep += (double) depth[i];
	  cnt ++;
	}
    }
  *dep = *dep / (double) cnt;
  
  free(depth);


  /* Number of segments in the library */
  
  *seglib = (double)gv->NChrom ;

  for ( cc = 0 ; cc < gv->NChrom ; cc ++ ) 
    for ( m=1 ; m < gv->chrom[cc].NLocC ; m++)
      for (i=0; i < gv->NoInd; i++)  
	{
	  l2 = gv->chrom[cc].first_loc_i + m ;
	  l1 = l2-1;

	  idx1 = l1*gv->NoInd + i;
	  idx2 = l2*gv->NoInd + i;

	  if (  (gv->All1[idx1] != gv->All1[idx2]) ||
		(gv->All2[idx1] != gv->All2[idx2])    )
	    {
	      (*seglib)++ ;
	      break;
	    }
	}
} 

				
void il_eval_library_GV ( double* cov    ,
			  double* dep    ,
			  double* seglib , 
			  int*    allele ,
			  int*    retval , 
			  char**  set_name  )
{
  il_eval_library ( gs_fdta(*set_name),
		    cov    ,
		    dep    ,
		    seglib ,
		    allele ,
		    retval) ;
}

//##############################################################################
//# BEGIN: Ridge regression after transformation                         2015-01
//##############################################################################

void gs_build_GI_01 ( gs_varset_type* gv ,
		      int*   retval        )
{
  // Inverse of the G matrix: G =ZZ'
  
  gv->t9  = realloc (gv->t9, gv->NoInd * gv->NoInd *sizeof(double) );
  gv->GI  = realloc (gv->GI, gv->NoInd * gv->NoInd *sizeof(double) );
  if (( NULL == gv->t9 )||( NULL == gv->GI )) ERR_M;

  gs_mtm_mul(gv->NoInd, gv->NoMar,  // ZZ'
  	     gv->NoInd, gv->NoMar,
  	     gv->Z, gv->Z, gv->t9);

  { 
    double d;
    gv->u0  = realloc ( gv->u0, gv->NoInd * sizeof(gs_UNSIGNED) );

    gs_ludcmp ( gv->t9, gv->NoInd, gv->u0, &d);
    gs_luinv  ( gv->t9, gv->NoInd, gv->u0 ,gv->GI );
   
  }

  *retval = 0; return ;
}

void gs_mmet_coeff_01 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  // Constant lambda
  // Calculates GI vor the back transformation

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->y ) {
      gs_info(gv,-2,
        "No performance data available, run 'read.performance.data' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  /* Allocate memory for the A matrix */

  gv->dimA = 1 + gv->NoInd;
  gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
  if ( NULL == gv->A ) ERR_M;

  for ( i = 0; i < gv->dimA ; i++) 
    for ( j = 0; j < gv->dimA ; j++) 
	gv->A [ i*gv->dimA + j ] =  0;

  /* Init part 1*/
  
  gv->A [0] = (double) gv->NoInd;

  /* Init part 2 and 3*/

  for ( j = 0; j < gv->NoInd; j++)  
    {
      gv->A [     0*gv->dimA + (1+j) ] = 1;
      gv->A [ (j+1)*gv->dimA + 0     ] = 1;
    }

  /* Init part 4*/

  gs_build_GI_01(gv,Sretval); // Used also in solve for the back transfomation
  if (-2 == *Sretval) { *retval = -2; return; }

  for ( i = 0; i < gv->NoInd; i++)
    for ( j = 0; j < gv->NoInd; j++)
	gv->A [ (1+i)*gv->dimA + (1+j) ] = 
                                      gv->GI[ i*gv->NoInd + j] * gv->ladi[0];

  for ( i = 0; i < gv->NoInd; i++)
    gv->A [ (1+i)*gv->dimA + (1+i) ] += 1;


  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mmet_coeff_01_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_mmet_coeff_01 (gs_fdta(*set_name),
		    out_filename,
		    auxfiles, 
		    retval        ) ;
}

void gs_mmet_coeff_02 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  // Coefficient matrix for gs_lambda_emstep_tc_01. Constant shrinkage
  // factor. For use iterations, re-uses a stored inv(G) matrix, 
  // calculated by gs_build_GI_01

  gs_UNSIGNED i,j, idx_A;
  FILE *outfp;

  if ( NULL == gv->y ) {
      gs_info(gv,-2,
        "No performance data available, run 'read.performance.data' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  /* Allocate memory for the A matrix */

  gv->dimA = 1 + gv->NoInd;
  gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
  if ( NULL == gv->A ) ERR_M;
 
  for ( i = 0; i < gv->dimA ; i++) 
    for ( j = 0; j < gv->dimA ; j++) 
	gv->A [ i*gv->dimA + j ] =  0;

  /* Init part 1*/
  
  gv->A [0] = (double) gv->NoInd;

  /* Init part 2 and 3*/

  for ( j = 0; j < gv->NoInd; j++)  
    {
      gv->A [     0*gv->dimA + (1+j) ] = 1;
      gv->A [ (j+1)*gv->dimA + 0     ] = 1;
    }

  /* Init part 4*/
  
  // Build A 

  for ( i = 0; i < gv->NoInd; i++)
    for ( j = 0; j < gv->NoInd; j++)
	gv->A [ (1+i)*gv->dimA + (1+j) ] = 
	                gv->ladi[0] * gv->GI[ i*gv->NoInd + j];

  for ( i = 0; i < gv->NoInd; i++)
    gv->A [ (1+i)*gv->dimA + (1+i) ] += 1;


  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mmet_coeff_02_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_mmet_coeff_02 (gs_fdta(*set_name),
		    out_filename,
		    auxfiles, 
		    retval        ) ;
}

void gs_mmet_coeff_03 ( gs_varset_type* gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{

  // MME for unequal shrikage factors, Error variance
  // must be available in gv->sigma[0]. From an rmlct run.
  // Provides GI for the back transformation

  gs_UNSIGNED i,j,m, idx_A;
  FILE *outfp;

  if ( NULL == gv->y ) {
      gs_info(gv,-2,
        "No performance data available, run 'read.performance.data' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
      gs_info(gv,-2,
        "No Z matrix available, run 'build.Z' first");
      *retval = -2; return ;
  }

  if ( NULL == gv->ladi ) {
      gs_info(gv,-2,
        "No shrinkage factor defined");
      *retval = -2; return ;
  }

  /* Allocate memory for the A matrix */

  gv->dimA = 1 + gv->NoInd;
  gv->A =  realloc( gv->A, gv->dimA * gv->dimA * sizeof(double) );
  if ( NULL == gv->A ) ERR_M;

  for ( i = 0; i < gv->dimA ; i++) 
    for ( j = 0; j < gv->dimA ; j++) 
	gv->A [ i*gv->dimA + j ] =  0;

  /* Init part 1*/
  
  gv->A [0] = (double) gv->NoInd;

  /* Init part 2 and 3*/

  for ( j = 0; j < gv->NoInd; j++)  
    {
      gv->A [     0*gv->dimA + (1+j) ] = 1;
      gv->A [ (j+1)*gv->dimA + 0     ] = 1;
    }

  /* Init part 4*/

  double *WZ;
  WZ = malloc ( gv->NoInd * gv->NoMar * sizeof(double) );
  if ( NULL == WZ ) ERR_M;

  double *GS;
  GS = malloc ( gv->NoInd * gv->NoInd *sizeof(double) );
  if ( NULL == GS ) ERR_M;

  gv->GI = realloc (gv->GI, gv->NoInd * gv->NoInd *sizeof(double) );
  if ( NULL == gv->GI ) ERR_M;

  for ( i = 0; i < gv->NoInd ; i++) 
    for ( m = 0; m < gv->NoMar ; m++)
      WZ[  i*gv->NoMar + m ] = gv->Z[ i*gv->NoMar + m ] * 
	(gv->sigma[0]/gv->ladi[m]);                      // this is var_i


  gs_mtm_mul(gv->NoInd, gv->NoMar,  // GS = ZWZ'
  	     gv->NoInd, gv->NoMar,
  	     gv->Z, WZ, GS);

  { // GI = inv (GS)
    double d;
    gv->u0  = realloc ( gv->u0, gv->NoInd * sizeof(gs_UNSIGNED) );

    gs_ludcmp ( GS, gv->NoInd, gv->u0, &d);
    gs_luinv  ( GS, gv->NoInd, gv->u0 ,gv->GI );
   
  }
  
  // Build A 

  for ( i = 0; i < gv->NoInd; i++)
    for ( j = 0; j < gv->NoInd; j++)
	gv->A [ (1+i)*gv->dimA + (1+j) ] = gv->GI[ i*gv->NoInd + j];

  for ( i = 0; i < gv->NoInd; i++)              // divide by var_e
    for ( j = 0; j < gv->NoInd; j++)
      gv->A [ (1+i)*gv->dimA + (1+j) ] = 
	             gv->A [ (1+i)*gv->dimA + (1+j) ] * gv->sigma[0];

  for ( i = 0; i < gv->NoInd; i++)              // add 1 to diagonal
    gv->A [ (1+i)*gv->dimA + (1+i) ] += 1;

  free(WZ);
  free(GS);

  /* Write A matrix to output file */
  if (1 == *auxfiles)
    {
      if ( NULL == (outfp = fopen(*out_filename,"w")) ) ERR_F;

      for ( i = 0; i < gv->dimA ; i++) {
	  for ( j = 0; j < gv->dimA ; j++) {
	      idx_A   = i*gv->dimA + j ;
	      fprintf(outfp,"%f ",gv->A [idx_A] );
	  }
	  fprintf(outfp,"\n");
      }
      
      fclose(outfp);  
    }

  *retval = 0; return ;
}

void gs_mmet_coeff_03_GV ( char** out_filename,
			  int*   auxfiles, 
			  int*   retval,
			  char** set_name                )
{
  gs_mmet_coeff_03 (gs_fdta(*set_name),
		    out_filename,
		    auxfiles, 
		    retval        ) ;
}



void gs_mmet_rhs_01 ( gs_varset_type *gv ,
		     char** out_filename,
		     int*   auxfiles, 
		     int*   retval        )
{
  gs_UNSIGNED i,j, dim_b;
  FILE *outfp;

  if ( NULL == gv->y ) {
    gs_info(gv,-2,
      "No performance data available, run 'read.performance.data' first");
    *retval = -2; return ;
  }

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available, run 'build.Z' first");
    *retval = -2; return ;
  }

  /* Allocate memory for b  (Ax=b) */
  dim_b = 1 + gv->NoInd;
  gv->b =  realloc( gv->b, dim_b * sizeof(double) );
  if ( NULL == gv->b ) ERR_M;

  /* Init part 1*/
  
  gv->b[0] = 0;
  for ( j = 0; j < gv->NoInd; j++) 
    gv->b[0] += gv->y[j];

  /* Init part 2*/
  
  for ( j = 0; j < gv->NoInd; j++)  
    gv->b [ 1+j ] = gv->y[j];

  /* Write b vector to output file */
  if (1 == *auxfiles) {
  
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      for ( i = 0; i < dim_b ; i++)
	fprintf(outfp,"%f \n",gv->b [i] );
      
      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mmet_rhs_01_GV ( char** out_filename,
	 		int*   auxfiles, 
	 		int*   retval,
	 		char** set_name    )
{
  gs_mmet_rhs_01 ( gs_fdta(*set_name) ,
		  out_filename,
		  auxfiles, 
		  retval        ) ;
}

void gs_mmet_solve_01 ( gs_varset_type *gv ,
		      char** out_filename,
		      int*   auxfiles, 
		      int*   retval        )
{
  gs_UNSIGNED i;
  FILE *outfp;

  int Srv=0; int* Sretval = &Srv;

  if ( NULL == gv->A ) {
      gs_info(gv,-2, "No coefficient matrix available");
      *retval = -2; return ;
  }

  if ( NULL == gv->b )  {
      gs_info(gv,-2, "No right hand side available");
      *retval = -2; return ;
  }

  /* Allocate memory for the solution vector of the transformed mode*/
  gv->a =  realloc( gv->a, (1+gv->NoInd) * sizeof(double) );
  if ( NULL == gv->a ) ERR_M;

  /* Solution vector for the transformed data*/
  gs_solve_Axb ( gv->dimA, gv->dimA, gv->A, gv->a, gv->b, Sretval );
  if (*Sretval==-2) {*retval = -2; return ;}

  /* Allocate memory for the solution vector of marker effects*/
  gv->dimu = 1 + gv->NoMar;
  gv->u =  realloc( gv->u, gv->dimu * sizeof(double) );
  if ( NULL == gv->u ) ERR_M;

  double *tmp1;
  tmp1 = malloc ( gv->NoInd * gv->NoMar * sizeof(double) );
  if ( NULL == tmp1 ) ERR_M;

  // Calculate  inv(ZZ') without shrinkage factor
  gs_build_GI_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_tmm_mul(gv->NoInd, gv->NoMar,  // Z' GSI
  	     gv->NoInd, gv->NoInd,
  	     gv->Z, gv->GI, tmp1);

  gs_mm_mul    ( gv->NoMar, gv->NoInd,  // (Z' GSI ) a
  	         gv->NoInd,         1,
  	         tmp1, 1+gv->a, 1+gv->u);

  gv->u[0] = gv->a[0];

  free(gv->A); gv->A =NULL;
  free(tmp1);

  /* Write solution vector to output file */
  if (1 == *auxfiles) {
   
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"effect\n");
      for ( i = 0; i < gv->dimu ; i++)
	fprintf(outfp,"%s.%i %f \n",
		gv->EffNme[i],
		gv->EffAll[i],
		gv->u [i]);

      fclose(outfp);  
  }

  *retval = 0; return ;
}

void gs_mmet_solve_01_GV ( char** out_filename ,
			  int*   auxfiles     , 
			  int*   retval       ,
			  char** set_name      )
{
  gs_mmet_solve_01 ( gs_fdta(*set_name) ,
		    out_filename,
		    auxfiles, 
	            retval        ) ;
}


void gs_lambda_emstep_tc_01 ( gs_varset_type *gv     ,
			      char**  out_filename   ,
			      int*    auxfiles       , 
			      int*    retval         )
{

  gs_UNSIGNED i,j,m;

  if (  NULL== gv->ladi ) {
    gs_info(gv,-2,"No lambda defined");
    *retval = -2; return ;
  }

  int Srv; int* Sretval = &Srv;

  gs_mmet_coeff_02 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_mmet_rhs_01 ( gv, nixp, neinp, Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  /* Allocate memory for the solution vector */
  gv->a =  realloc( gv->a, gv->dimA * sizeof(double) );
  if ( NULL == gv->a ) ERR_M;

  /* Initialize the solution vector */
  for ( i = 0; i < gv->dimA ; i++) gv->a [i] =  gv->b[i];

  { // LU decomposition 
    double d;

    gv->u0  = realloc ( gv->u0, gv->dimA * sizeof(gs_UNSIGNED) );

    /* LU decomposition */
     gs_ludcmp ( gv->A, gv->dimA, gv->u0, &d);

    /* Solutions for a */
     gs_lubksb ( gv->A, gv->dimA, gv->u0, gv->a );

    /* Inverse */
     gv->t4 = realloc (  gv->t4, gv->dimA*gv->dimA * sizeof(double) );

    gs_luinv ( gv->A, gv->dimA, gv->u0 ,gv->t4 );

  }

  /* Allocate memory for the variance components */
  gv->sigma =  realloc( gv->sigma, gv->NoMar * sizeof(double) );
  if ( NULL == gv->sigma ) ERR_M;

  gv->t1 = realloc (gv->t1, sizeof(double));

  gv->t2 = realloc (gv->t2, 
		    ( (gv->NoInd > gv->NoMar) ?  gv->NoInd : gv->NoMar    ) 
		    * sizeof(double)            );

  gv->t3 = realloc (gv->t3, 
		    ( (gv->NoInd > gv->NoMar) ?  gv->NoInd : gv->NoMar    ) 
		    * sizeof(double)            );

  for ( i=0 ; i<gv->NoInd ; i++) gv->t2[i] = 1;

  gs_tmm_mul    (gv->NoInd,1 ,           // y'y
	         gv->NoInd,1 ,	     
	         gv->y       ,
	         gv->y       ,
	         gv->sigma ) ;

  gs_tmm_mul    (gv->NoInd,1 ,           // 1'y
	         gv->NoInd,1 ,	     
	         gv->t2          ,
	         gv->y       ,
	         gv->t1          ) ;  

  gs_tmm_mul    (gv->NoInd, 1 ,	         // a'y
	         gv->NoInd, 1 ,
	         1+(gv->a)    ,
	         gv->y        ,
	         gv->t3           ) ;
  
  // Residual error is estimated 

  gv->sigma[0] -= gv->a[0] *gv->t1[0];  // y'y - b0 1' y
  gv->sigma[0] -= gv->t3[0];            //     - a'Z'y         
  gv->sigma[0]  = gv->sigma[0] / (double)gv->NoInd;
  gv->sigma_e   = gv->sigma[0];

  gs_tmm_mul    (gv->NoInd, 1,           // a'GI	     
	         gv->NoInd, gv->NoInd ,
	         1+(gv->a)            ,
	         gv->GI               ,
	         gv->t3                ) ;

  gs_mtm_mul    (1,         gv->NoInd ,  // a' GI a
	         1,         gv->NoInd , 
	         gv->t3               ,
	         1+(gv->a)            ,
	         gv->t2                ) ;

  /* Diagonal elements of the inverse coefficient matrix */

  gv->t5 = realloc ( gv->t5, gv->NoInd * gv->NoInd * sizeof(double) );
  for ( j=1 ; j < gv->dimA ; j++ ) 
    memcpy ( & gv->t5 [ (j-1) * gv->NoInd ] ,
	     & gv->t4 [  j    * gv->dimA   + 1],
	     (gv->NoInd) * sizeof(double) );

  gv->t6 = realloc ( gv->t6, gv->NoInd * gv->NoInd * sizeof(double) );
  gs_mm_mul  (gv->NoInd, gv->NoInd ,          
	      gv->NoInd, gv->NoInd ,
	      gv->GI               ,
	      gv->t5               ,
	      gv->t6                ) ;

  double tr = 0; 
  for ( i=0; i < gv->NoInd; i++ ) tr +=gv->t6 [  gv->NoInd *i + i ];

  gv->sigma[1] = (gv->t2[0] + gv->sigma[0] * tr ) / (double) gv->NoInd;

  // Update the shrinkage factor
  for (i=0; i <gv->NoMar; i++)
    gv->ladi[i] = gv->sigma[0] / gv->sigma[1] ;

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;

      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda var \n");
      fprintf(outfp,"err %f %f \n", 0.0, gv->sigma[0] );
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f %f\n",
		gv->MarName[m],
		gv->ladi[m] ,
		gv->sigma[1+m]
		);

      fclose(outfp);  
    } 

  *retval = 0; return ;

}

void gs_lambda_emstep_tc_01_GV ( char**  out_filename  ,
			         int*    auxfiles      , 
			         int*    retval        ,
			         char** set_name       )

{
  gs_lambda_emstep_tc_01 (  gs_fdta(*set_name),
			    out_filename,
			    auxfiles, 
			    retval         ) ;
}

void gs_lambda_rmlct_01 ( gs_varset_type* gv ,
		 	  int*    maxiter,
			  double* precision,
			  double* hsq,
			  char**  out_filename,
			  int*    auxfiles, 
			  int*    retval         )
{

  gs_UNSIGNED m;
  int i;
  double qu, maxqu=0;
  int Srv; int* Sretval = &Srv;

  if ( NULL == gv->Z ) {
    gs_info(gv,-2,"No Z matrix available");
    *retval = -2; return ;
  }

  // Save variance components to determine convergence
  double vc_save[2] = {0,0};

  // Initialize Vector of shrinkage factors lambda
  gv->ladi = realloc(gv->ladi,gv->NoMar*sizeof(double));
  if (NULL==gv->ladi) ERR_M;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    gv->ladi[m] = ( 1 / *hsq - 1 ) * (double)gv->NoMar;

  // Calculate  inv(ZZ')
  gs_build_GI_01(gv,Sretval);
  if (-2 == *Sretval) { *retval = -2; return; }
      
  for (i=0; i < *maxiter; i ++) 
    {

      if (i>0) memcpy ( vc_save, gv->sigma,  2 * sizeof(double) ); 

      gs_lambda_emstep_tc_01 ( gv,nixp,neinp,Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      maxqu = 0;
      for ( m = 0 ; m < 2 ; m++) {
	qu = fabs ( (vc_save[m] - gv->sigma[m]) / gv->sigma[m] );
	if (qu>maxqu) maxqu=qu;
      }
      if (maxqu < *precision) break; 
  }
  
  sprintf(gv->msg, "Number of EM iterations: %i, precision: %f",i,maxqu);
  gs_info(gv,1,gv->msg);
  sprintf(gv->msg, "(Greatest change of a var. comp. in the last iteration: %6.3f%%)",
	  maxqu*100);
  gs_info(gv,1,gv->msg);
    
  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      FILE *outfp;
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}


void gs_lambda_rmlct_01_GV ( int*    maxiter        ,
			     double* precision      ,
			     double* hsq            ,
			     char**  out_filename   ,
			     int*    auxfiles       , 
			     int*    retval         ,
			     char**  set_name         )
{
  gs_lambda_rmlct_01 ( gs_fdta(*set_name),
		       maxiter,
		       precision,
		       hsq,
		       out_filename,
		       auxfiles, 
		       retval         ) ;
}

void gs_lambda_rmlat_01 (  gs_varset_type *gv     ,
			   double* alpha          ,
			   int*    maxiter        ,
			   double* precision      ,
			   double* hsq            ,
			   char**  out_filename   ,
			   int*    auxfiles       , 
			   int*    retval         )
{
  gs_UNSIGNED m;
  FILE *outfp;
  int Srv; int* Sretval = &Srv;

  gs_single_marker_aov_01 ( gv, alpha, nixp, neinp, Sretval );
  if (-2 == *Sretval) { *retval = -2; return; }

  gs_lambda_rmlct_01 ( gv,
		       maxiter,precision,hsq,
		      out_filename,auxfiles,Sretval ) ;
  if (-2 == *Sretval) { *retval = -2; return; }

  long double Svar = 0.0;
  for ( m = 0 ; m < gv->NoMar ; m++)  
    Svar += gv->SMvar[m] ;

  Svar = Svar / (long double) gv->NoMar;

  for ( m = 0 ; m < gv->NoMar ; m++) 
    gv->ladi[m] =  ( Svar / gv->SMvar[m]  ) * gv->ladi[m];

  /* Write solution vector to output file */
  if (1 == *auxfiles)
    {
      if ( (outfp = fopen(*out_filename,"w")) == NULL ) ERR_F;

      fprintf(outfp,"lambda \n");
      for ( m = 0; m < gv->NoMar ; m++)
	fprintf(outfp,"%s %f\n",
		gv->MarName[m],
		gv->ladi[m]
		);

      fclose(outfp);  
    }  
  *retval = 0; return ;

}

void gs_lambda_rmlat_01_GV ( double* alpha,
			     int*    maxiter       ,
			     double* precision     ,
			     double* hsq           ,
			     char**  out_filename  ,
			     int*    auxfiles      , 
			     int*    retval        ,
			     char** set_name        )
{
  gs_lambda_rmlat_01 ( gs_fdta(*set_name),
		       alpha             ,
		       maxiter           ,
		       precision         ,
		       hsq               ,
		       out_filename      ,
		       auxfiles          , 
		       retval             ) ;
}


void gs_esteff_rr_01 (  gs_varset_type *gv ,
		        char**  methode,
                        int*    maxiter,
		        double* precision,
		        double* hsq,
			double* alpha,
		        char**  out_filename,
		        int*    auxfiles, 
		        int*    retval         )
{
  int Srv; int* Sretval = &Srv;


  if (strcmp (*methode, "rrwet") == 0) 
    {

      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_const_02 ( gv, hsq, nixp, neinp, retval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

    } 

  else if (strcmp (*methode, "rmlct") == 0) 
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmlct_01 ( gv, maxiter, precision, hsq,
			   nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    } 

  else if (strcmp (*methode, "rmlat") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmlat_01 ( gv, alpha, maxiter, precision, hsq,
			   nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_coeff_03 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mmet_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rmlc") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmlc_01 ( gv,
			  maxiter,precision,hsq,
		      nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rmla") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmla_01 ( gv,
			  alpha,
			  maxiter,precision,hsq,
			  nixp,neinp,Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rmlv") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_rmlv_01 ( gv,
			  maxiter,precision,hsq,
			  nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rrwe") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_const_02 ( gv, hsq, nixp, neinp, retval );
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
      
      gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else if (strcmp (*methode, "rrwa") == 0)
    {
      gs_build_Z_01 ( gv , nixp, neinp, Sretval );
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_lambda_aov_01 ( gv, hsq, alpha, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_coeff_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_rhs_01 ( gv, nixp, neinp, Sretval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }

      gs_mme_solve_01 ( gv, out_filename, auxfiles, retval ) ;
      if (-2 == *Sretval) { *retval = -2; return; }
    }

  else 
    { 
      gs_info(gv,-2,"Unknown shrinkage methode");
      *retval = -2;  return;
     }

  *retval = 0; return ;
}

void gs_esteff_rr_01_GV (  char**  methode,
			   int*    maxiter,
		           double* precision,
			   double* hsq,
			   double* alpha,
		           char**  out_filename,
		           int*    auxfiles, 
		           int*    retval,
		           char** set_name         )
{
  gs_esteff_rr_01 ( gs_fdta(*set_name),
                    methode,
		    maxiter,   
	            precision, 
		    hsq,
		    alpha,
	            out_filename,
	            auxfiles, 
	            retval         ) ;
}

void gs_cross_validation_02 ( gs_varset_type *gv   ,
			      char ** scheme       ,
			      int*    n_estimation ,
			      int*    n_runs       , 
			      double* hsq          , 
			      double* alpha        ,
			      int*    maxiter      ,
			      double* precision    ,
			      int*    gibbsburn    ,
			      int*    gibbsiter    ,
			      char**  out_filename_r , 
			      char**  out_filename_u , 
			      int*    auxfiles     , 
			      int*    retval       )
 { 
   if ( RAND_MAX < (int) gv->NoInd )  {
       gs_info(gv,-2,"More individuals than RAND_MAX");
       *retval = -2; return ;
     }

   if ( *n_estimation > (int) gv->NoInd )  {
       gs_info(gv,-2,"Training is set larger than base population");
       *retval = -2; return ;
     }


   int gi,gb;
   gi=*gibbsiter; gb=*gibbsburn;

   int ntr = 0;
   #pragma omp parallel
   {
     ntr = omp_get_num_threads();
   }
   if (0==ntr) ntr=1;

#if defined(use_openblas) && defined(MT)
  int eins = 1;
  openblas_set_num_threads(eins); 
  goto_set_num_threads(eins); 
  omp_set_num_threads (ntr);
#endif

#if defined(use_mkl) && defined(MT)
  int eins = 1;
  mkl_set_num_threads(eins); 
  omp_set_num_threads (ntr);
#endif


   if (!rgs) { srand( (unsigned) time(NULL) ); rgs=1; }

   int CVNoRuns = *n_runs;  

   /* Results */
   double* CVcor = malloc ( CVNoRuns * sizeof(double) );
   if (NULL == CVcor) ERR_M;
   for (gs_UNSIGNED c = 0; c < (gs_UNSIGNED) CVNoRuns ; c++) CVcor[c] = 0;

   gs_FLOAT* CVest = malloc ( CVNoRuns * (1+gv->NoMar) *sizeof(gs_FLOAT) );


     /* --- once for each thread */
     #pragma omp parallel
     {
       int Srv; int* Sretval = &Srv;
 
       gs_varset_type ES; gs_varset_type *es = &ES;
       gs_init (es);

       gs_varset_type VS; gs_varset_type *vs = &VS;
       gs_init (vs);
       
       gs_UNSIGNED i;
       double cor;       
       double sy=0, syy=0, sx=0, sxx=0, sxy=0, k;

       #pragma omp for
       for (int r=0 ; r < CVNoRuns; r ++)  
	 /* --- once for each CVRun */
	 {

	   cor = 0;

	   // Build estimation and validation set
	   gs_build_esvs_01(gv,es,vs,n_estimation, 
			 nixp,nixp,nixp,nixp,neinp, Sretval);
	   if (-2 == *Sretval) { *retval = -2;  }


	   // Estimate effects in estimation set 
	   gs_esteff_rr_01 ( es,scheme,maxiter,precision,hsq,alpha,
			     nixp,neinp,Sretval);
	   if (-2 == *Sretval) { *retval = -2;  }

	   // Save effect estimates
	   for (i =0; i < es->dimu; i++) CVest[(r*es->dimu)+i] = es->u [i];

	   // Estimate breeding values in validation set 
	   gs_estimate_bv_01 (es,vs,nixp,neinp,Sretval);
	   if (-2 == *Sretval) { *retval = -2;}

	   // Correlation estimated - observed 

	   sy=0, syy=0, sx=0, sxx=0, sxy=0;
	   k = (double) vs->NoInd;
	   
	   for ( i = 0; i < vs->NoInd; i++ )  {
	       sx  += vs->est[i] ;
	       sxx += vs->est[i] * vs->est[i] ;
	       sy  += vs->y[i] ;
	       syy += vs->y[i] * vs->y[i] ;
	       sxy += vs->est[i] * vs->y[i] ;
	     }
	   cor =    (      k*sxy - sx*sy            ) /
	       sqrt ( (k*sxx-sx*sx) * (k*syy-sy*sy) ) ;

	   CVcor[r] = cor;
	   

	 } // for r
       
       gs_reset(es);
       gs_reset(vs);

     } // parallel


#if defined(use_openblas) && defined(MT)
  openblas_set_num_threads(ntr); 
  goto_set_num_threads(ntr); 
#endif

#if defined(use_mkl) && defined(MT)
  mkl_set_num_threads(ntr); 
#endif


     gs_UNSIGNED i,r;
     gs_FLOAT sum;

     /* Allocate memory for the solution vector */

     /* Allocate memory for the solution vector */
     int Srv; int* Sretval = &Srv;

     char* meth_ = "rrwet";   char **meth = &meth_;
    
     gs_esteff_rr_01 ( gv,meth,maxiter,precision,hsq,alpha,
		       nixp,neinp,Sretval);

     /* Average effects estimated in CV */
     for ( i=0 ; i < gv->dimu; i++)
       {
	 sum = 0;
	 for ( r=0; r < (gs_UNSIGNED)CVNoRuns ; r++)
	   sum += CVest[(r*gv->dimu)+i] ;
	 gv->u[i] = (double)( sum / (double)CVNoRuns );
       }

   /* Write results to output file */
   if (1 == *auxfiles)
     {
       FILE *fp;
       if ( (fp = fopen(*out_filename_r,"w")) == NULL ) ERR_F;
       fprintf(fp,"cor \n");
       for (int c = 0; c < CVNoRuns ; c++)
       	 fprintf(fp,"%i %f\n",c, CVcor[c] );
       fclose(fp);  

       if ( (fp = fopen(*out_filename_u,"w")) == NULL ) ERR_F;

       fprintf(fp,"effect\n");
       for ( i = 0; i < gv->dimu ; i++)
	 fprintf(fp,"%s.%i %f \n",
		 gv->EffNme[i],
		 gv->EffAll[i],
		 gv->u [i]);
       fclose(fp);  
     }
   
   free(CVcor);
   free(CVest);

   *retval = 0; return;

 } 


void gs_cross_validation_02_GV ( char ** scheme       ,
				 int*    n_estimation ,
				 int*    n_runs       , 
				 double* hsq          , 
				 double* alpha        ,
				 int*    maxiter      ,
				 double* precision    ,
				 int*    gibbsburn,
				 int*    gibbsiter,
				 char**  out_filename_r , 
				 char**  out_filename_u , 
				 int*    auxfiles     , 
				 int*    retval       ,
				 char** set_name      )
 {
   gs_cross_validation_02 ( gs_fdta(*set_name),
			    scheme       ,
			    n_estimation ,
			    n_runs       , 
			    hsq          , 
			    alpha        ,
			    maxiter      ,
			    precision    ,
			    gibbsburn    ,
			    gibbsiter    ,
			    out_filename_r , 
			    out_filename_u , 
			    auxfiles     , 
			    retval           ) ;
 }

//##############################################################################
//# END: Ridge regression after transformation                           2015-01
//##############################################################################
