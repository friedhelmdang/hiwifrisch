/////////////////////////////////////////////////////////////////////////
// SeletionTools main source file genomic-selection.c                   //
// (c) Matthias Frisch                                                   //
///////////////////////////////////////////////////////////////////////////


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
