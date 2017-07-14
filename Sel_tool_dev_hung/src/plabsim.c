/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PLABSIM VERSION 3 */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

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

/* not used

#define EXPYEAR1 "2017"
#define EXPYEAR2 "2015"
#define EXPYEAR3 "2016"

*/

#define AUTHOR "$Author: frisch-m $"
#define DATE "$Date: 2015/11/29 00:00:00 $"
#define REVISION "$Revision: 15.2 $"

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* FUNCTION DECLARATIONS */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void info( int level, char*info, ... );
void remove_genotype_population( char **PopNames );
void remove_evaluate_population( char **PopNames );
void flush_population( char**PopName, int mode );

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* RANDOM NUMBERS */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int rng_initialized = 0;

#ifdef use_gsl

   #include <gsl/gsl_rng.h> 
   #include <gsl/gsl_randist.h>

   gsl_rng *rng[MAXTHREADS];

   const gsl_rng_type **rng_default; // pointer to array element
   const gsl_rng_type **all_types;   // array

   char rng_def_nme [] ="taus";
   char * pp_rng_def_nme = rng_def_nme;
   char ** p_rng_def_nme = &pp_rng_def_nme;

   void rng_init ();

   void rng_info()
   {
     info(0,"Random number generator type: '%s'\n",gsl_rng_name(rng[0]));
   }

   void rng_choose ( char ** name)
   {
     for ( rng_default=all_types; *rng_default != 0; rng_default++)
       {
	 if( 0== strcmp ( (*rng_default)->name, *name) ) break;
       }
     if (*rng_default == 0) rng_choose( p_rng_def_nme);
     rng_init();
     rng_info();
   }


   void rng_list()
   {
     const gsl_rng_type**t;
     
     Rprintf("Available generators:\n");
     for( t = all_types; *t != 0; t++)
       {
       if (0==strcmp(gsl_rng_name(rng[0]),(*t)->name))
	 {
	   info(0,"--> %s <--\n",(*t)->name);
	 }
       else
	 {
	   info(0,"    %s\n",(*t)->name);
	 }
     }
   }

#endif


#ifndef use_gsl

   #define IA 16807
   #define IM 2147483647
   #define AM (1.0/IM)
   #define IQ 127773
   #define IR 2836
   #define NTAB 32
   #define NDIV (1+(IM-1)/NTAB)
   #define EPS 1.2e-7
   #define RNMX (1.0-EPS)

   typedef struct {
     long seed       ;
     long iyy        ;
     long ivv[NTAB]  ;
   } gsl_rng; 

   gsl_rng* rng[MAXTHREADS];



  double gsl_rng_uniform (gsl_rng* rng)
  { 
    int j;
    long k;
  
    long* idum = &(rng->seed);
    long* iy   = &(rng->iyy );
    long* iv   = &(rng->ivv[0]);
  
    double temp;
  
    if ( (*idum) <= 0 || !(*iy) ) 
      {
	if ( -(*idum) < 1 ) (*idum)=1;
	else *idum = -(*idum);
	for ( j = NTAB+7; j >= 0; j-- ) 
	  {
	    k = (*idum) / IQ;
	    *idum = IA * (*idum - k * IQ) - IR * k;
	    if ( *idum < 0) *idum += IM;
	    if (j < NTAB) iv[j] = *idum;
	  }
	(*iy) = iv[0];
      }
    k = (*idum) / IQ;
    *idum = IA * ( *idum - k * IQ ) - IR * k;
    if (*idum < 0) *idum += IM;
    j = (*iy) / NDIV;
    (*iy) = iv[j];
    iv[j] = *idum;
    if ( ( temp = AM * (*iy) ) > RNMX )  return RNMX;
    else return temp;
  }

  #undef IA
  #undef IM
  #undef AM
  #undef IQ
  #undef IR
  #undef NTAB
  #undef NDIV
  #undef EPS
  #undef RNMX

  double gammln(double xx)
  { 
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
  
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x+0.5) * log(tmp);
    ser = 1.000000000190015;
    for ( j = 0; j <= 5; j++ ) ser += cof[j]/++y;
    return -tmp + log(2.5066282746310005*ser/x);
  }


  #define rPI 3.141592654

  unsigned int gsl_ran_poisson (gsl_rng* rng, double xm)
  {
    static double sq,alxm,g,oldm=(-1.0);
    double em,t,y;
    
    if(xm < 12.0)
      {
	if(xm != oldm)
	  {
	    oldm = xm;
	    g = exp(-xm);
	  }
	em = -1;
	t = 1.0;
	do
	  {
	    ++em;
	    t *= gsl_rng_uniform (rng);
	  } while (t>g);
      }
    else
      {
	
	if(xm != oldm)
	  {
	    oldm = xm;
	    sq = sqrt(2.0*xm);
	    alxm = log(xm);
	    g = xm*alxm-gammln(xm+1.0);
	  }
	do
	  {
	    do
	      {
		y = tan(rPI* gsl_rng_uniform (rng) );
		em = sq*y+xm;
	      } while(em < 0.0);
	    
	    em = floor(em);
	    t = 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
	  } while ( gsl_rng_uniform (rng) > t);
      }
    return round(em);
  }

  #undef rPI


  double gsl_rng_uniform_pos (gsl_rng* rng)
  {
    double r = gsl_rng_uniform (rng) ;

    if ( ( r == (double)0 ) || ( r == (double)1) ) return gsl_rng_uniform (rng);
    else return r;
  }

  unsigned long int gsl_rng_uniform_int (gsl_rng* rng, unsigned long int n)
  {
    double r = gsl_rng_uniform (rng) ;

    if ( r == (double)1 )  return ( gsl_rng_uniform_int (rng,n) );
    else return ( floor ( r * n) );
  }

  gsl_rng* gsl_rng_alloc()
  {
    return (gsl_rng*) malloc ( sizeof (gsl_rng) );
  }

  void gsl_rng_set (gsl_rng * rng, unsigned long int s)
  {
    rng->seed = -(long) s;
    rng->iyy = 0;
  }

  void gsl_rng_free (gsl_rng * r)
  {
    free(r);
  }

  void rng_choose(char**name_rng)
  {
    char** dummy = name_rng; dummy++; // avoid compiler warning

    info(0, "Random number generator: Numerical Recipies rand1() \n" );
   }

   void rng_info()
   {
     info(0, "Random number generator type: Numerical Recipies rand1() \n");
   }

   void rng_list()
   {
     info(0, "Available generators: Numerical Recipies rand1() \n" );
   }

  #define p_rng_default

#endif


void rng_init () {
 
  int i;
  unsigned long int seed;

  srand( (unsigned)time(NULL) );

  #ifdef use_gsl
    if ( 0 == rng_initialized )
      {
	all_types = gsl_rng_types_setup();
	
	for ( rng_default=all_types; *rng_default != 0; rng_default++)
	  {
	    if( 0== strcmp ( (*rng_default)->name, rng_def_nme) ) break;
	  }
      }
    const gsl_rng_type * p_rng_default = *rng_default;
  #endif

  if ( 1 == rng_initialized)
    {
      for (i =0; i < MAXTHREADS; i++ )
	gsl_rng_free (rng[i]);
      info(1, "Re-initializing random number generator \n" );
    }
  else
    {
      info(1, "Initializing random number generator \n" );
    }

  for (i =0; i < MAXTHREADS; i++ )
    {
      seed = (unsigned long) rand();
      rng[i] = gsl_rng_alloc (p_rng_default);
      gsl_rng_set (rng[i],seed);
    }

  rng_initialized = 1;

}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* PLABSIM MAIN  */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#define MAXNAME 256
#define MAXLINE 100
#define CPOS double
#define POPULATION long int

#define CCTIME __TIME__
#define CCDATE __DATE__

/* #define xCalloc(...)Calloc(__VA_ARGS__) */
/* #define xRealloc(...)Realloc(__VA_ARGS__) */
/* #define xFree(...)Free(__VA_ARGS__) */

/* Alternative memory allocation without R wrappers */

#define xCalloc(n,t)(t*)malloc((size_t)(n)*sizeof(t))
#define xRealloc(p,n,t)(t*)realloc((void*)(p),(size_t)(n)*sizeof(t))
#define xFree(p)(free((void*)(p)),(p)= NULL)

/* Type definitions */

typedef struct
{
    CPOS Pos;
    int  chrom;
    char *name;
    char *class;
    int  counter;
} 
TMappoint;

typedef struct
{
    double    ChromLen;
    int       NrMappoints;
    TMappoint *Mappoint;
} 
TChromDescr;

typedef struct
{
  CPOS           Pos;
  unsigned short All;
}
TLoc;

typedef struct
{
  int  NoLociAlloc;
  int  NoLociUsed;
  TLoc *Loc;
}
THom;

typedef struct
{
  THom *Hom;
}
TChrom;

typedef struct
{
  TChrom *Chrom;
  float  *GValue;
  int    *genotype;
  char   *info;
}
TInd;

typedef struct
{
  char       Name[ MAXNAME ];
  POPULATION NoIndsAlloc;
  POPULATION NoInds;
  TInd *Ind;
}
TPop;

typedef struct
{
    int  sizeRows, sizeCols;
    int  nPop, nMarker;
    int  *nAll, *nInd;
    int  *nAll_off, *nInd_off;
    int  *data;
}
TPopMatrix;

struct TIntLinkList
{
    int number;
    struct TIntLinkList *child, *parent;
    struct TIntLinkList *next,  *prev;
};

typedef struct TIntLinkList IntNode;

typedef struct SPopListEl
{
    TPop* Pop;
    struct SPopListEl *Prev;
    struct SPopListEl *Next;
}
TPopListEl;

typedef struct
{
    TMappoint *MP;
    int  NrAll;
    int  *All;
}
TGen;

typedef struct
{
    float Eff;
    int   NrLoci;
    TGen  *Loc;
}
TEff;

typedef struct
{
    float GValue_avg;
    float weight;
    int   NrEff;
    int   Pos;
    TEff  *Eff;
    char  Name[ 3*MAXNAME ];
    int   MPs;
    TMappoint**MPset;
}
TEffects;


typedef struct
{
    int  betrag;
    int  *member;
}
TSet;

/* Variable declarations  */

int     NoChroms = 0;
int     NoHoms = 0;
int     NoLociInit = 10;
IntNode *startNode = NULL;
double  COFreq = 1;

TPopListEl *PopListEl = NULL;
int        bg, missing, doubleHetero;
int        dimLA, dimLB, locusA, locusB, nrGeno;
int        **table, *genotype;


TPopListEl  *PopListStart = NULL;
TChromDescr *PopDescr = NULL;
int         cMP = 0;
int         cMPhom = 0;
int         INFOLEVEL = 0;

TEffects*Effects = NULL;
int NrEffects = 0;
TSet*AllSet;
unsigned long int MDP = 0;

/* Generic functions */

void info( int level, char*info, ... )
{
    if ( INFOLEVEL >= level )
    {
        va_list ap;
        char*p, *sval, Info[ 200 ] = "";
        double dval;
        POPULATION ival;
        int pos = 0;

        va_start( ap, info );

        for ( p = info;*p;p++ )
        {
            if ( *p != '%' )
            {
                Info[ pos++ ] = *p;
                Info[ pos ] = '\0';
                continue;
            }

            switch ( *++p )
            {

            case'd':
                {
                    char tmp[ 30 ] = "";
                    ival = va_arg( ap, POPULATION );

                    sprintf( tmp, "%ld", ival );
                    strcat( Info, tmp );
                    pos = strlen( Info );
                    break;
                }

            case'f':
                {
                    char tmp[ 30 ] = "";
                    dval = va_arg( ap, double );

                    sprintf( tmp, "%f", dval );
                    strcat( Info, tmp );
                    pos = strlen( Info );
                    break;
                }

            case's':

                for ( sval = va_arg( ap, char* );*sval;sval++ )
                {
                    Info[ pos++ ] = *sval;
                    Info[ pos ] = '\0';
                }

                break;

            default:
                Info[ pos++ ] = *p;
                Info[ pos ] = '\0';
                break;
            }
        }

        va_end( ap );

        if ( 0 <= level )
            Rprintf( "M: %s", Info );
        else if ( -1 == level )
            Rprintf( "W: %s", Info );
        else if ( -2 == level )
            Rprintf( "E: %s", Info );
    }
}


void r_info(int*level, char**message)
{
  info( *level, "%s\n", *message );
}

void plabsim_info_cctime_get( char**message )
{
    sprintf( *message, "%s - %s", CCTIME, CCDATE );
}

 
void print_version()
{
    info( 0, "%s (%s): Rev. %s \n", AUTHOR, DATE, REVISION );
    info( 0, "compiled: %s - %s\n", CCTIME, CCDATE );
}

void print_version_2( char **message )
{
  char * co_date_complete = DATE;
  char co_date[20];
  strncpy(co_date,co_date_complete+7,10);
  co_date[10] = '\0';
  info( 0, "plabsim c %s\n", co_date );
  info( 0, "plabsim r %s\n", *message );
}

FILE *fileopen ( char* path,
		 char* mode)
{
    FILE*fp;

    if ( ( fp = fopen( path, mode ) ) == NULL )
        info( -2, "%Can't open the file %s\n", path );

    return fp;
}

void fileclose( FILE*fp )
{
    if ( fclose( fp ) != 0 )
        info( -2, "Can't close the file stream\n");
}


void set_info_level( int*level )
{
  INFOLEVEL = *level;
  if ( INFOLEVEL < -2 )  INFOLEVEL = -2;
  if ( INFOLEVEL > 10 )  INFOLEVEL = 10;
  info(2,"Info level set to %d\n",INFOLEVEL);
}

void get_info_level (int* i)
{
  *i = INFOLEVEL;
}

/* Simulations functions */

void find_list_element ( char        *Name, 
			 TPopListEl  **PopListEl )
{
    *PopListEl = PopListStart;

    while ( ( *PopListEl != NULL ) && 
	    ( strcmp( ( *PopListEl ) ->Pop->Name, Name ) != 0 ) )
      {
	( *PopListEl ) = ( *PopListEl ) ->Next;
      }
}


void remove_individuum ( TInd *DelInd )
{
    TChrom  VDelChrom;
    TChrom *DelChrom = &VDelChrom;
    int chrom, hom;

    for ( chrom = 0; chrom < NoChroms; chrom++ ) 
    {
      DelChrom = &(DelInd->Chrom[chrom]);
      for ( hom = 0; hom < NoHoms; hom++ ) xFree (DelChrom->Hom[hom].Loc);
      xFree (DelChrom->Hom);
    }

    xFree (DelInd->Chrom);

    if ( NULL != DelInd->GValue )   xFree (DelInd->GValue);
    if ( NULL != DelInd->genotype ) xFree (DelInd->genotype);
    if ( NULL != DelInd->info )     xFree (DelInd->info);
}


void new_individuum ( TInd *NewInd )
{
    int    chrom, hom;
    THom   VNewHom;
    THom   *NewHom = &VNewHom;
    TChrom VNewChrom;
    TChrom *NewChrom = &VNewChrom;

    NewInd->Chrom    = (TChrom*) xCalloc( NoChroms, TChrom );
    NewInd->GValue   = NULL;
    NewInd->genotype = NULL;
    NewInd->info     = NULL;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        NewChrom = &( NewInd->Chrom[ chrom ] );
        NewChrom->Hom = ( THom* ) xCalloc( NoHoms, THom );

        for ( hom = 0;hom < NoHoms;hom++ )
        {
            NewHom = &( NewChrom->Hom[ hom ] );
            NewHom->NoLociAlloc = NoLociInit;
            NewHom->NoLociUsed = 0;
            NewHom->Loc = ( TLoc* ) xCalloc( NoLociInit, TLoc );
        }
    }
}


TInd* get_individual ( TPop* Pop, POPULATION nrInd )
{
    return &( Pop->Ind[ nrInd ] );
}


void set_nolociinit ( int* NewLociInit)
{
  if ( 0 < *NewLociInit )  NoLociInit = *NewLociInit;
}


void remove_map ()
{
  int chrom, i, j = 0;

  if ( NULL != PopDescr ) 
  {
    for ( chrom = 0;chrom < NoChroms;chrom++ ) {
      if ( PopDescr[ chrom ].NrMappoints != 0 ) {
	j = 1;
	for ( i = 0; i < PopDescr[chrom].NrMappoints; i++ )
	  {
	    xFree( PopDescr[ chrom ].Mappoint[ i ].name );
	    xFree( PopDescr[ chrom ].Mappoint[ i ].class );
	  }
	xFree (PopDescr[chrom].Mappoint);
	PopDescr[chrom].NrMappoints = 0;
      }
    }
    
    cMP = 0;
    cMPhom = 0;
      
    if ( 1 == j ) info( 1, "Linkage map removed\n" );
  }
}


void set_genome_par ( int*    _NoChroms, 
		      int*    _NoHoms, 
		      double* _ChromLen )
{
  time_t ltime;
  char   tme[64];

  time(&ltime);strcpy(tme,ctime(&ltime));

  /* Not used
  if (
      (strstr(tme,EXPYEAR1)==NULL)&&
      (strstr(tme,EXPYEAR2)==NULL)&&
      (strstr(tme,EXPYEAR3)==NULL)
      ) 
    {
      info(-2,"Compilation date: %s - %s\n", CCTIME, CCDATE );
      info(-2,"Software update required: %s\n",tme);
      return ;
    }
  */
  
  if ( NULL != PopListStart ) {
    info( -2, "First remove all populations\n" );
    return;
  }

  if ( ( PopDescr != NULL ) && ( cMP > 0 ) ) {
    info( -2, "First remove previously loaded linkage map\n" );
    return;
  }

  int chrom;

  NoChroms = *_NoChroms;
  NoHoms   = *_NoHoms;

  if ( NULL == PopDescr )
    PopDescr = (TChromDescr*) xCalloc( NoChroms, TChromDescr );
  else
    PopDescr = (TChromDescr*) xRealloc( PopDescr, NoChroms, TChromDescr );
  
  for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
      PopDescr[ chrom ].ChromLen = ( double ) _ChromLen[ chrom ];
      PopDescr[ chrom ].NrMappoints = 0;
      PopDescr[ chrom ].Mappoint = NULL;
    }

  info( 0, "%d chromosomes defined\n", NoChroms);

}

void optimize_population ( char** PopNames )
{
  int         chrom, hom;
  POPULATION  ind;
  TPopListEl* PopListEl;
  TPop*       Pop = NULL;
  char*       PopName;

  PopName = strtok( *PopNames, "[' ','\n']" );
  
  while ( NULL != PopName )
    {
      find_list_element( PopName, &PopListEl );
      
      if ( PopListEl == NULL ) {
	info( -2, "Can't find population %s\n", PopName );
	return;
      }
      
      Pop = PopListEl->Pop;
      
      if ( NULL != Pop->Ind )
	{
	  for ( ind = 0; ind < Pop->NoInds; ind++ )
	    for ( chrom = 0; chrom < NoChroms; chrom++ )
	      for ( hom = 0; hom < NoHoms; hom++ )
		{
		  THom*Hom = &( ( *Pop ).Ind[ind].Chrom[chrom].Hom[hom] );
		  
		  if ( ( Hom->NoLociUsed + NoLociInit ) 
		       <= Hom->NoLociAlloc                )
		    {
		      Hom->NoLociAlloc = ( Hom->NoLociUsed / NoLociInit ) * 
			                 NoLociInit;
		      
		      if ( 0 == Hom->NoLociAlloc )
			Hom->NoLociAlloc += NoLociInit;
		      
		      if ( Hom->NoLociAlloc < Hom->NoLociUsed )
			Hom->NoLociAlloc += NoLociInit;
		      
		      Hom->Loc = ( TLoc* ) 
			xRealloc( Hom->Loc, Hom->NoLociAlloc, TLoc );
		    }
		}
	  
	  if ( ( Pop->NoInds ) != ( Pop->NoIndsAlloc ) )
	    {
	      info( -1, "NoIndsAlloc is greater than NoInds\n" );
	      info( -1, "use resize to remove unused individuals\n" );
	      
	      for ( ind = Pop->NoInds;ind < Pop->NoIndsAlloc;ind++ )
		{
		  for ( chrom = 0;chrom < NoChroms;chrom++ )
		    {
		      for ( hom = 0;hom < NoHoms;hom++ )
			{
			  THom*Hom = &( ( *Pop ).Ind[ ind ].Chrom[ chrom ].Hom[ hom ] );
			  Hom->NoLociAlloc = NoLociInit;
			  Hom->Loc = ( TLoc* ) xRealloc( Hom->Loc, Hom->NoLociAlloc, TLoc );
			}
		    }
		}
	    }
	}
      
      PopName = strtok( NULL, "[' ','\n']" );
    }
}


void int_resize_population( char**PopName, POPULATION newSize )
{
  POPULATION ind;
  TPopListEl*PopListEl;
  TPop*Pop = NULL;
  
  find_list_element( *PopName, &PopListEl );

  if ( PopListEl == NULL ) {
    info( -2, "Can't find population %s\n", *PopName );
    return;
  }

  if ( newSize < 0 ){
    info( -2, "Population size must be a positive number \n" );
    return;
  }
  
  Pop = PopListEl->Pop;

  if ( newSize > Pop->NoIndsAlloc )
    {
      Pop->Ind = ( TInd* ) xRealloc( Pop->Ind, newSize, TInd );
      
      for ( ind = Pop->NoIndsAlloc; 
	    ind < newSize; 
	    ind++ )
	new_individuum( &( Pop->Ind[ ind ] ) );
      
      Pop->NoIndsAlloc = newSize;
    }
  
  if ( newSize < Pop->NoIndsAlloc )
    {
      if ( Pop->NoInds > newSize )
	Pop->NoInds = newSize;

      for ( ind = newSize; 
	    ind < Pop->NoIndsAlloc;
	    ind++ )
	remove_individuum( &( Pop->Ind[ ind ] ) );
      
      Pop->Ind = ( TInd* ) xRealloc( Pop->Ind, newSize, TInd );
      Pop->NoIndsAlloc = newSize;
    }

  info( 1, "Population %s resized to %d individuals\n",	*PopName, newSize );
}


void resize_population ( char** PopName, char** cnewSize )
{
    int_resize_population( PopName, atol( *cnewSize ) );
}

TPop* new_population ( char** Name, POPULATION* NoInds )
{
  POPULATION  ind;
  TPop*       NewPop = NULL;
  TPopListEl* NewPopListEl;
  
  if ( 2 > strlen( *Name ) ) {
    info( -2, "Population name is too short\n" );
    return NULL;
  }

  if ( 0 > *NoInds ) {
    info( -2, "Incorrect number of individuals\n" );
    return NULL;
  }
  
  find_list_element( *Name, &NewPopListEl );
  
  if ( NewPopListEl != NULL )  
    NewPop = NewPopListEl->Pop;
  
  if ( NULL != NewPop ) {
    info( -2, "Population already exists\n" );
    return NULL;
  }

  NewPop = (TPop*) xCalloc( 1, TPop );
  NewPopListEl = (TPopListEl*) xCalloc( 1, TPopListEl );

  strcpy( NewPop->Name, *Name );
  
  NewPop->NoIndsAlloc = *NoInds;
  NewPop->NoInds      = 0;
  NewPopListEl->Pop   = NewPop;
  NewPopListEl->Prev  = NULL;
  NewPopListEl->Next  = PopListStart;

  if ( NewPopListEl->Next != NULL ) 
    NewPopListEl->Next->Prev = NewPopListEl;

  PopListStart = NewPopListEl;

  if ( NewPop->NoIndsAlloc > 0 )
    {
      NewPop->Ind = (TInd*) xCalloc( NewPop->NoIndsAlloc, TInd );
      for ( ind = 0; ind < NewPop->NoIndsAlloc; ind++ )
        {
	  new_individuum( &( NewPop->Ind[ ind ] ) );
        }
    }
  else
    {
      NewPop->Ind = NULL;
    }

    info( 1, "Population %s with %d individuals generated\n", *Name, *NoInds );
    return ( NewPop );
}

void remove_population( char** PopNames )
{
  char*       token;
  TPopListEl* DelPopListEl;
  char*       _PopNames;
  
  _PopNames = ( char* ) xCalloc( strlen( *PopNames ) + 2, char );
  
  strcpy( _PopNames, *PopNames );
  token = strtok( _PopNames, "[' ','\n']" );
  
  while ( NULL != token )
    {
      find_list_element( token, &DelPopListEl );
      if ( DelPopListEl != NULL )
	{
	  POPULATION ind;
	  for ( ind = 0; 
		ind < DelPopListEl->Pop->NoIndsAlloc; 
		ind++ )
	    {
	      remove_individuum( &( DelPopListEl->Pop->Ind[ ind ] ) );
	    }
	  
	  if ( ( DelPopListEl->Prev == NULL ) && 
	       ( DelPopListEl->Next == NULL )     )
	    {
	      PopListStart = NULL;
	    }
	  else
            {
	      if ( DelPopListEl->Prev == NULL )
                {
		  DelPopListEl->Next->Prev = NULL;
		  PopListStart = DelPopListEl->Next;
                }
	      else
                {
		  if ( DelPopListEl->Next == NULL )
                    {
		      DelPopListEl->Prev->Next = NULL;
                    }
		  else
                    {
		      DelPopListEl->Next->Prev = DelPopListEl->Prev;
		      DelPopListEl->Prev->Next = DelPopListEl->Next;
                    }
                }
            }
	  
	  if ( NULL != DelPopListEl->Pop->Ind )
	    xFree( DelPopListEl->Pop->Ind );

	  xFree( DelPopListEl->Pop );
	  xFree( DelPopListEl );

	  info( 1, "Population %s deleted\n", token );
        }
      
      token = strtok( NULL, "[' ','\n']" );
    }
  xFree( _PopNames );
}


TPop *exist_population ( char*      PopName,
			 POPULATION NoInds,
			 int        empty )
{
  TPopListEl *PopListEl;

  find_list_element( PopName, &PopListEl );

  if ( PopListEl == NULL )
    {
      return new_population( &PopName, &NoInds );
    }

  if ( PopListEl->Pop->NoIndsAlloc != NoInds )
    {
      int_resize_population( &PopName, NoInds );
    }

  if ( 1 == empty )
    {
      remove_genotype_population( &PopName );
      remove_evaluate_population( &PopName );
    }
  
  return PopListEl->Pop;
}

int insertLocus( unsigned short All, 
		 CPOS  Pos, 
		 THom* Hom )
{
  if ( ( 0 == Hom->NoLociUsed ) ||
       ( ( Hom->Loc[ Hom->NoLociUsed - 1 ].Pos < Pos ) &&
	 ( Hom->Loc[ Hom->NoLociUsed - 1 ].All != All ) ) )
    {
      TLoc * Loc;
      
      if ( Hom->NoLociUsed == Hom->NoLociAlloc )
        {
	  Hom->NoLociAlloc += NoLociInit;
	  Hom->Loc = ( TLoc* ) xRealloc( Hom->Loc, NoLociInit + Hom->NoLociAlloc, TLoc );
        }
      
      Loc = &( Hom->Loc[ Hom->NoLociUsed ] );
      Loc->Pos = Pos;
      Loc->All = All;
      Hom->NoLociUsed++;
      return 0;
    }

  return 1;
}


void init_population( char** Name,
		      char** cCount,
		      int*_Ind,
		      int*_Chrom,
		      int*_Hom,
		      double*_Pos,
		      int*_All        )
{
  int min_chrom = 0, 
    max_chrom = 0,
    min_hom   = 0, 
    max_hom    = 0;
  POPULATION row, 
    Count = atol( *cCount );
  POPULATION max_Ind = 0, 
    min_Ind = 0;
  TPopListEl* PopListEl;
  TPop* Pop;
  THom* Hom;
  
  find_list_element( *Name, &PopListEl );
  
  if ( PopListEl != NULL ) {
    info( -2, "First remove population %s\n", *Name );
    return;
  }

  for ( row = 0; row < Count; row++ )
    {
      max_chrom = ( *( _Chrom + row ) - 1 > max_chrom ) ? *( _Chrom + row ) - 1 : max_chrom;
      min_chrom = ( *( _Chrom + row ) - 1 < min_chrom ) ? *( _Chrom + row ) - 1 : min_chrom;
      max_hom = ( *( _Hom + row ) - 1 > max_hom ) ? *( _Hom + row ) - 1 : max_hom;
      min_hom = ( *( _Hom + row ) - 1 < min_hom ) ? *( _Hom + row ) - 1 : min_hom;
      max_Ind = ( *( _Ind + row ) - 1 > max_Ind ) ? *( _Ind + row ) - 1 : max_Ind;
      min_Ind = ( *( _Ind + row ) - 1 < min_Ind ) ? *( _Ind + row ) - 1 : min_Ind;
    }
  
  if ( ( min_chrom < 0 ) || ( max_chrom >= NoChroms ) ) {
    info( -2, "Error in the chromosome column\n" );
    return;
  }
 
  if ( ( min_hom < 0 ) || ( max_hom >= NoHoms ) ) {
    info( -2, "Error in the homologue column\n" );
    return;
  }
  
  if ( min_Ind < 0 ) {
    info( -2, "Error in the individual column\n" );
    return;
  }
  
  Pop = exist_population( *Name, max_Ind + 1, 1 );
  if ( NULL == Pop )  return;
  
  Pop->NoInds = max_Ind + 1;
  
  for ( row = 0; row < Count; row++ )
    {
      TInd*Ind = get_individual( Pop, *( _Ind + row ) - 1 );
      Hom = &( Ind->Chrom[*(_Chrom+row)-1].Hom[*(_Hom+row)-1] );
      insertLocus( _All[ row ], _Pos[ row ], Hom );
      
      if ( 0 == _All[ row ] )
	info( -1, "Don't use 0 to code alleles\n" );
    }
}


void get_allele_number( char** Name, 
			int*Count  )
{
    int         chrom, hom;
    POPULATION  ind;
    TPopListEl* PopListEl;
    TPop*       Pop;
    TInd*       Ind;
    TChrom*     Chrom;
    THom*       Hom;

    find_list_element( *Name, &PopListEl );

    if ( PopListEl != NULL )
      {
        *Count = 0;
        Pop = PopListEl->Pop;
        for ( ind = 0; ind < Pop->NoInds; ind++ )
	  {
            Ind = get_individual( PopListEl->Pop, ind );
            for ( chrom = 0;chrom < NoChroms;chrom++ )
	      {
                Chrom = &( Ind->Chrom[ chrom ] );
                for ( hom = 0; hom < NoHoms; hom++ )
		  {
                    Hom = &( Chrom->Hom[ hom ] );
                    *Count += Hom->NoLociUsed;
		  }
	      }
	  }
      }
}


void get_population( 
		    char**  PopName,
		    int*    _Ind,
		    int*    _Chrom,
		    int*    _Hom,
		    double* _Pos,
		    int*    _All      )
{
  int chrom, hom, loc;
  POPULATION ind;
  POPULATION row = 0;
  TPopListEl*PopListEl;
  TPop*Pop;
  TInd*Ind;
  TChrom*Chrom;
  THom*Hom;
  TLoc*Loc;
  
  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl ) {
    info( -2, "Can't find population %s\n", *PopName );
    return;
  }
  
  Pop = PopListEl->Pop;
  
  for ( ind = 0;ind < Pop->NoInds;ind++ )
    {
      Ind = get_individual( PopListEl->Pop, ind );
      for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
	  Chrom = &( Ind->Chrom[chrom] );
	  for ( hom = 0;hom < NoHoms;hom++ )
            {
	      Hom = &( Chrom->Hom[hom] );
	      for ( loc = 0; loc < Hom->NoLociUsed; loc++ )
                {
		  Loc = &( Hom->Loc[ loc ] );
		  _Ind[ row ] = ind + 1;
		  _Chrom[ row ] = chrom + 1;
		  _Hom[ row ] = hom + 1;
		  _Pos[ row ] = ( double ) Loc->Pos;
		  _All[ row++ ] = Loc->All;
                }
            }
        }
    }
}


void list_populations_parameter ( int* lPopName, int* NoPops )
{
  TPopListEl * PopListEl;
  PopListEl = PopListStart;

  while ( PopListEl != NULL )
    {
      (*NoPops) ++ ;
      (*lPopName) += strlen( PopListEl->Pop->Name ) + 1;
      PopListEl = PopListEl->Next;
    }
  ( *lPopName ) += 10;
}


void list_populations ( char **all_pops, int *inds )
{
  TPopListEl * PopListEl;
  PopListEl = PopListStart;
  
  strcpy( *all_pops, "" );
  
  while ( PopListEl != NULL )
    {
      *inds = PopListEl->Pop->NoInds;
      inds = inds + 1;
      strcat( *all_pops, PopListEl->Pop->Name );
      strcat( *all_pops, " " );
      PopListEl = PopListEl->Next;
    }
}


void set_co_freq( double*cofreq )
{
  if ( 0 <= *cofreq )
    {
      COFreq = *cofreq;
      info( 1, "Crossoverfrequency set to %f\n", *cofreq );
    }
  else
    info( -1, "Recombination Frequency must be a positive number\n" );
}


int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}


void meiosis( THom*   PgHom, 
	      THom*   PHom, 
	      double  Length, 
	      double* significant,
	      int     r)
{
  unsigned int NuCO = gsl_ran_poisson(rng[r], COFreq * Length );
  CPOS posCO [NuCO + 2] ;
  unsigned int i, aHom;
  int aPos[] = {0, 0};
  CPOS pos;

  if ( -1 < *significant )
    {
      for ( i = 1;i <= NuCO;i++ )
        {
	  posCO[i] = fround( gsl_rng_uniform_pos(rng[r]) * Length, *significant );
	  posCO[i] = ( *(posCO+i) > Length ) ? Length : *(posCO+i);
        }
    }
  else
    {
      for ( i = 1;i < NuCO + 1;i++ )
        {
	  posCO[i] = gsl_rng_uniform_pos(rng[r]) * Length;
        }
    }

  posCO[ 0 ] = 0;
  posCO[ NuCO + 1 ] = Length;

  switch ( NuCO )
    {

    case 0:

    case 1:
      break;

    case 2:
      if ( ( posCO[ 1 ] > posCO[ 2 ] ) )
        {
	  pos = posCO[ 1 ];
	  posCO[ 1 ] = posCO[ 2 ];
	  posCO[ 2 ] = pos;
        }
      break;

    case 3:
      if ( ( posCO[ 1 ] < posCO[ 2 ] ) && ( posCO[ 2 ] < posCO[ 3 ] ) )
	break;

    default:
      // R_rsort( posCO + 1, NuCO );
      qsort ( posCO+1 , NuCO, sizeof (double), compare_doubles);

    }

  aHom = ( int ) gsl_rng_uniform_int( rng[r], 2 );

    for ( i = 0; i <= NuCO; i++ )
      {
        while ( ( aPos[ aHom ] + 1 < PHom[ aHom ].NoLociUsed ) &&
                ( PHom[ aHom ].Loc[ aPos[ aHom ] ].Pos < posCO[ i ] ) &&
                ( PHom[ aHom ].Loc[ aPos[ aHom ] + 1 ].Pos <= posCO[ i ] ) )
	  {
            aPos[ aHom ] ++ ;
	  }
        while ( ( aPos[ aHom ] < PHom[ aHom ].NoLociUsed ) &&
                ( PHom[ aHom ].Loc[ aPos[ aHom ] ].Pos < posCO[ i + 1 ] ) )
	  {
            pos = PHom[ aHom ].Loc[ aPos[ aHom ] ].Pos;
            pos = ( pos > posCO[ i ] ) ? pos : posCO[ i ];
            insertLocus( PHom[ aHom ].Loc[ aPos[ aHom ] ].All, pos, PgHom );
            aPos[ aHom ] ++;
	  }
        if ( aPos[ aHom ] > 0 )
	  aPos[ aHom ] -- ;
        aHom = ( aHom == 1 ) ? 0 : 1;
      }
}

/* qsort int comparison function */ 
int intcmp(const void *a, const void *b) 
{ 
    const int *ia = (const int *)a; 
    const int *ib = (const int *)b;
    return *ia  - *ib; 
} 

/* qsort unsigned short comparison function */ 
int unsignedshortcmp(const void *a, const void *b) 
{ 
    const unsigned short *ia = (const unsigned short *)a; 
    const unsigned short *ib = (const unsigned short *)b;
    return *ia  - *ib; 
} 


void cross ( 
	   char**  NamePg,  
	   char**  NameP1,
	   char**  NameP2,
	   double* Self,
	   char**  cNoPg,
	   double* significant,
	   int*    mode,
	   int*    crossing_scheme )
{
  POPULATION ind, 
             NoPg = atol( *cNoPg );

  TPopListEl *PopListEl;

  TPop *PgPop = NULL, 
       *P1Pop = NULL, 
       *P2Pop = NULL;


  if ( 2 != NoHoms ) {
    info( -2, "Implemented for diploids only\n" );
    return;
  }

  if ( 1 > NoPg ) {
    info( -2, "Number of progenies must be positive\n" );
    return;
  }

  find_list_element( *NamePg, &PopListEl );
  if ( PopListEl != NULL )  PgPop = PopListEl->Pop;

  find_list_element( *NameP1, &PopListEl );
  if ( PopListEl != NULL )  P1Pop = PopListEl->Pop;

  find_list_element( *NameP2, &PopListEl );
  if ( PopListEl != NULL )  P2Pop = PopListEl->Pop;


  if ( ( P1Pop == NULL ) ||  ( P2Pop == NULL ) ) {
    info( -2, "Parent population not defined\n" );
    return;
  }
  
  if ( ( P1Pop == PgPop ) || ( P2Pop == PgPop ) ) {
    info( -2, "Progeny population incorrect\n" );
    return;
  }

  if ( ( 0 == P1Pop->NoInds ) || ( 0 == P2Pop->NoInds ) ){
    info( -2, "Parental population of size 0\n" );
    return;
  }

  PgPop = exist_population( *NamePg, NoPg, 1 );
  if ( NULL == PgPop )  return;

  if ( -1 == NoPg )  
    PgPop->NoInds = PgPop->NoIndsAlloc;
  else 
    PgPop->NoInds = NoPg;

  #pragma omp parallel  
  {
    int chrom;
    double r;
    POPULATION ind1, ind2; 
    TInd *P1Ind, *P2Ind, *PgInd;

    int i = omp_get_thread_num();

  #pragma omp for 
  for ( ind = 0; ind < PgPop->NoInds; ind++ )
  {

    if ( ( 0 == *crossing_scheme ) || ( 1 == *crossing_scheme ) ) {
      P1Ind = get_individual( P1Pop, 
			      ind1 = (POPULATION)gsl_rng_uniform_int ( rng[i], 
                                     (unsigned long int)P1Pop->NoInds ) );
    }
    else {
      P1Ind = get_individual( P1Pop, 
			      ind1 = ind % P1Pop->NoInds );
    }

    if ( ( 0 == *crossing_scheme ) || ( 2 == *crossing_scheme ) ) {
      P2Ind = get_individual( P2Pop,
			      ind2 = (POPULATION) gsl_rng_uniform_int ( rng[i], 
                                     (unsigned long int)P2Pop->NoInds ) );
    }
    else {
      P2Ind = get_individual( P2Pop, 
			      ind2 = ind % P2Pop->NoInds );
    }

    if ( ( 0 == *mode ) && 
	 ( P1Ind == P2Ind ) && 
	 ( 1 < P1Pop->NoInds ) )
    {
      while ( P1Ind == P2Ind )
	{
	  P2Ind = get_individual( P2Pop,
				  ind2 = (POPULATION)gsl_rng_uniform_int( rng[i],
	                          (unsigned long int)P2Pop->NoInds) );
	}
    }

    r = gsl_rng_uniform_pos(rng[i]);
    PgInd = get_individual( PgPop, ind );

    if ( *Self / 2 > r )
      {
	P2Ind = P1Ind;
	info( 2, "  %s (%d) selfed -> %s (%d)\n", *NameP1, ind1, *NamePg, ind );
      }
    else
      {
	if ( *Self > r )
	  {
	    P1Ind = P2Ind;
	    info( 2, "  %s (%d) selfed -> %s (%d)\n", *NameP2, ind2, *NamePg, ind );
	  }
	else
	  {
	    info( 2, "  %s (%d) x %s (%d) -> %s (%d)\n",
		  *NameP1,
		  ind1,
		  *NameP2,
		  ind2,
		  *NamePg,
		  ind );
	  }
      }
    
    if ( 2 == NoHoms )
      {
	for ( chrom = 0;chrom < NoChroms;chrom++ )
	  {
	    PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed = 0;
	    PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed = 0;
	    
	    meiosis( PgInd->Chrom[ chrom ].Hom,
		     ( *P1Ind ).Chrom[ chrom ].Hom,
		     PopDescr[ chrom ].ChromLen,
		     significant,
		     i);
	    
	    meiosis( &( PgInd->Chrom[ chrom ].Hom[ 1 ] ),
		     ( *P2Ind ).Chrom[ chrom ].Hom,
		     PopDescr[ chrom ].ChromLen,
		     significant,
		     i);
	  }
      }
  }
  } //parallel
}


void ps_set_num_threads (int * activate)
{
  if ( *activate > MAXTHREADS)  *activate = MAXTHREADS ;
  if ( *activate < 1 )          *activate = 1 ;
  omp_set_num_threads (*activate);
}

/*---------------- */
/*vvvvvvvvvvvvvvvv */

void ssd_mating( char**  NamePg,
		 char**  NameP1,
		 char**  cNoPg,
		 int*    maxcycles,
		 double* significant )
{
  POPULATION ind, NoPg = atol( *cNoPg ), nPg = atol( *cNoPg ), prog;
  int chrom;

  TPopListEl *PopListEl;
  TPop       *PgPop = NULL, 
             *P1Pop = NULL;
  TInd       *PgInd, 
             *P1Ind;
  THom       *tmpHomA, 
             *tmpHomB;
  int cycle, 
      nocycles = *maxcycles, 
      hom, 
      loc;

  if ( 2 != NoHoms ){
    info( -2, "only for diploid\n" );
    return;
  }
  
  if ( 1 > NoPg ){
    info( -2, "Number of progenies must be a positive number\n" );
    return;
  }

    find_list_element( *NamePg, &PopListEl );

    if ( PopListEl != NULL )
        PgPop = PopListEl->Pop;

    find_list_element( *NameP1, &PopListEl );

    if ( PopListEl != NULL )
        P1Pop = PopListEl->Pop;

    if ( P1Pop == NULL )
    {
        info( -2, "Can't find parent population: %s\n", *NameP1 );
	return;
    }

    if ( P1Pop == PgPop )
    {
        info(
            -2,
            "Progeny population name (%s) is identical to parent population name (%s)\n",
            *NamePg, *NameP1
        );
	return;
    }

    NoPg = P1Pop->NoInds * NoPg;

    tmpHomA = ( THom* ) xCalloc( NoHoms, THom );

    for ( hom = 0;hom < NoHoms;hom++ )
    {
        tmpHomA[ hom ].NoLociAlloc = NoLociInit;
        tmpHomA[ hom ].NoLociUsed = 0;
        tmpHomA[ hom ].Loc = ( TLoc* ) xCalloc( NoLociInit, TLoc );
    }

    PgPop = exist_population( *NamePg, NoPg, 1 );
    if ( NULL == PgPop )  return;

    PgPop->NoInds = NoPg;

    for ( ind = 0; ind < P1Pop->NoInds; ind++ )
      {
        P1Ind = get_individual( P1Pop, ind );
        for ( prog = 0;prog < nPg;prog++ )
	  {
            PgInd = get_individual( PgPop, ind * nPg + prog );
            for ( chrom = 0;chrom < NoChroms;chrom++ )
	      {
                int diff = 0;
                cycle = nocycles;
                tmpHomB = P1Ind->Chrom[ chrom ].Hom;
                do
		  {
                    diff = 0;
                    for ( hom = 0;hom < NoHoms;hom++ )
		      {
                        if ( tmpHomA[ hom ].NoLociAlloc < tmpHomB[ hom ].NoLociAlloc )
			  {
                            tmpHomA[ hom ].NoLociAlloc = tmpHomB[ hom ].NoLociAlloc;
                            tmpHomA[ hom ].Loc = ( TLoc* ) xRealloc( tmpHomA[ hom ].Loc,
								     tmpHomA[ hom ].NoLociAlloc,
								     TLoc );
                        }
                        tmpHomA[ hom ].NoLociUsed = tmpHomB[ hom ].NoLociUsed;
                        for ( loc = 0;loc < tmpHomA[ hom ].NoLociUsed;loc++ )
			  {
                            tmpHomA[ hom ].Loc[ loc ].All = tmpHomB[ hom ].Loc[ loc ].All;
                            tmpHomA[ hom ].Loc[ loc ].Pos = tmpHomB[ hom ].Loc[ loc ].Pos;
			  }
		      }
                    PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed = 0;
                    PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed = 0;

                    meiosis( PgInd->Chrom[ chrom ].Hom,
                             tmpHomA,
                             PopDescr[ chrom ].ChromLen,
                             significant,
			     0 );

                    meiosis( &( PgInd->Chrom[ chrom ].Hom[ 1 ] ),
                             tmpHomA,
                             PopDescr[ chrom ].ChromLen,
                             significant,
			     0 );

                    if ( PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed !=
                            PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed )
		      {
                        diff++;
		      }
                    else
		      {
                        int hom = 0;
                        do
			  {
                            if ( ( PgInd->Chrom[ chrom ].Hom[ 0 ].Loc[ hom ].All !=
				   PgInd->Chrom[ chrom ].Hom[ 1 ].Loc[ hom ].All ) ||
				 ( PgInd->Chrom[ chrom ].Hom[ 0 ].Loc[ hom ].Pos !=
				   PgInd->Chrom[ chrom ].Hom[ 1 ].Loc[ hom ].Pos ) )
			      {
                                diff++;
			      }
			  }
                        while ( ( !diff ) && ( ++hom < PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed ) );
                    }
                    tmpHomB = PgInd->Chrom[ chrom ].Hom;
		  }
                while ( ( 0 != --cycle ) && ( diff ) );
	      }
	  }
      }
    
    for ( hom = 0;hom < NoHoms;hom++ )
      {
        xFree( tmpHomA[ hom ].Loc );
      }
    
    xFree( tmpHomA );
}


void dh( char**  NamePg,
	 char**  NameP1,
	 char**  cNoPg,
	 double* significant )
{
    POPULATION ind, 
               NoPg = atol( *cNoPg );

    int chrom;

    TPopListEl*PopListEl;

    TPop *PgPop = NULL, 
         *P1Pop = NULL;

    TInd *PgInd, 
         *P1Ind;

    if ( 2 != NoHoms ) {
      info( -2, "Implemented for diploids only\n" );
      return;
    }

    if ( 1 > NoPg ) {
      info( -2, "Number of progenies must be positive\n" );
      return;
    }

    find_list_element( *NamePg, &PopListEl );
    if ( PopListEl != NULL ) PgPop = PopListEl->Pop;

    find_list_element( *NameP1, &PopListEl );
    if ( PopListEl != NULL ) P1Pop = PopListEl->Pop;


    if ( P1Pop == NULL )  {
      info( -2, "Parent population not defined\n" );
      return;
    }

    if ( P1Pop == PgPop ) {
      info( -2, "Progeny population incorrect\n" );
      return;
    }

    if  ( 0 == P1Pop->NoInds ) {
      info( -2, "Parental population of size 0\n" );
      return;
    }

    PgPop = exist_population( *NamePg, NoPg, 1 );
    if ( NULL == PgPop )  return;



    PgPop->NoInds = NoPg;

    for ( ind = 0;ind < PgPop->NoInds;ind++ )
    {
        P1Ind = get_individual( P1Pop,
           (POPULATION)gsl_rng_uniform_int(rng[0], (unsigned long int)P1Pop->NoInds ) );
        PgInd = get_individual( PgPop, ind );

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            int loc, hom;

            PgInd->Chrom[ chrom ].Hom[ 0 ].NoLociUsed = 0;
            PgInd->Chrom[ chrom ].Hom[ 1 ].NoLociUsed = 0;

            meiosis( PgInd->Chrom[ chrom ].Hom,
                     P1Ind->Chrom[ chrom ].Hom,
                     PopDescr[ chrom ].ChromLen,
                     significant,
		     0 );

            for ( hom = 1;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( PgInd->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( PgInd->Chrom[ chrom ] ).Hom[ 0 ];
                PgInd->Chrom[ chrom ].Hom[ hom ].NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    insertLocus( P2Hom->Loc[ loc ].All, P2Hom->Loc[ loc ].Pos, P1Hom );
                }
            }
        }
    }
}


void define_map( fName )
char**fName;
{
    FILE*fp = NULL;
    int i = 0, chrom, failed = 0;
    double pos;
    int counter = 0;
    char lociname[ MAXNAME ], lociclass[ MAXNAME ];

    if ( cMP != 0 ){
      info( -2, "First remove previously loaded linkage map\n" );
      return;
    }

    if ( NULL != ( fp = fileopen( *fName, "r" ) ) )
    {
        for ( i = 0;;i++ )
        {
            if ( fscanf( fp, "%d %lf %s %s[\n,\t]", &chrom, &pos, lociname, lociclass ) == EOF )
            {
                break;
            }

            else
            {
                if ( ( chrom > NoChroms ) || ( 1 > chrom ) )
                {
                    info( -1, "Wrong chromosome (%d)\n", chrom );
                    failed++;
                    continue;
                }

                if ( !( ( pos >= 0 ) && ( pos <= PopDescr[ chrom - 1 ].ChromLen ) ) )
                {
                    info( -1, "Chromosome %d: Locusposition outside of the chromosome (0 <= %f <= %f)\n",
                          chrom, pos, PopDescr[ chrom - 1 ].ChromLen );
                    failed++;
                    continue;
                }

                if ( PopDescr[ chrom - 1 ].NrMappoints > 0 )
                {
                    if ( PopDescr[ chrom - 1 ].Mappoint[ PopDescr[ chrom - 1 ].NrMappoints - 1 ].Pos > pos )
                    {
                        info( -1, "Not sorted! (%d Element)\n", i );
                        failed++;
                        continue;
                    }
                }

                PopDescr[ chrom - 1 ].NrMappoints++;
                cMP++;

                if ( 1 == PopDescr[ chrom - 1 ].NrMappoints )
                    PopDescr[ chrom - 1 ].Mappoint = ( TMappoint* ) xCalloc( 1, TMappoint );
                else
                    PopDescr[ chrom - 1 ].Mappoint = ( TMappoint* ) xRealloc
                                                     ( PopDescr[ chrom - 1 ].Mappoint, PopDescr[ chrom - 1 ].NrMappoints, TMappoint );

                {
                    TMappoint*tmpPoint = &PopDescr[ chrom - 1 ].Mappoint[ PopDescr[ chrom - 1 ].NrMappoints - 1 ];
                    ( *tmpPoint ).chrom = chrom;
                    ( *tmpPoint ).Pos = pos;
                    ( *tmpPoint ).name = ( char* ) xCalloc( strlen( lociname ) + 1, char );
                    ( *tmpPoint ).class = ( char* ) xCalloc( strlen( lociclass ) + 1, char );
                    ( *tmpPoint ).counter = counter++;
                    strcpy( ( *tmpPoint ).name, lociname );
                    strcpy( ( *tmpPoint ).class, lociclass );
                }
            }
        }

        info( 0, "Map file %s loaded\n", *fName );
        info( 0, "%d loci defined - %d errors \n", i, failed );
        fileclose( fp );
    }

    cMPhom = NoHoms * cMP;
}
 
void define_linkage_map(
    int*dim,
    int*chrom,
    double*pos,
    char**names,

    char**class
)
{
    int i;
    int chromosome;
    int counter = 0;
    TMappoint*tmpPoint;

    if ( cMP != 0 ){
      info( -2, "First remove previously loaded linkage map\n" );
      return;
    }

    for ( i = 0;i < *dim;i++ )
    {
        chromosome = *( chrom + i ) - 1;

        if ( ( chromosome >= NoChroms ) || ( 0 > chromosome ) )
        {
            info( -1, "Wrong chromosome (%d)\n", chromosome );
        }

        else
        {
            if ( !( ( *( pos + i ) >= 0 ) && ( *( pos + i ) <= PopDescr[ chromosome ].ChromLen ) ) )
            {
                info( -1, "Chromosome %d: Locusposition outside of the chromosome (0 <= %f <= %f)\n",
                      chromosome + 1, *( pos + i ), PopDescr[ chromosome ].ChromLen );
            }

            else
            {
                PopDescr[ chromosome ].NrMappoints++;
                cMP++;

                if ( 1 == PopDescr[ chromosome ].NrMappoints )
                    PopDescr[ chromosome ].Mappoint = ( TMappoint* ) xCalloc( 1, TMappoint );
                else
                    PopDescr[ chromosome ].Mappoint = ( TMappoint* ) xRealloc
                                                      ( PopDescr[ chromosome ].Mappoint, PopDescr[ chromosome ].NrMappoints, TMappoint );

                tmpPoint = &PopDescr[ chromosome ].Mappoint[ PopDescr[ chromosome ].NrMappoints - 1 ];

                tmpPoint->chrom = chromosome + 1;

                ( *tmpPoint ).Pos = *( pos + i );

                ( *tmpPoint ).name = ( char* ) xCalloc( strlen( *( names + i ) ) + 1, char );

                ( *tmpPoint ).class = ( char* ) xCalloc( strlen( *( class + i ) ) + 1, char );

                ( *tmpPoint ).counter = counter++;

                strcpy( ( *tmpPoint ).name, *( names + i ) );

                strcpy( ( *tmpPoint ).class, *( class + i ) );
            }
        }
    }
}



int mccmp( const void*i, const void*j )
{
    int c0, c1;
    CPOS p0, p1;

    c0 = ( (TMappoint*) i )->chrom;
    c1 = ( (TMappoint*) j )->chrom;
    p0 = ( (TMappoint*) i )->Pos;
    p1 = ( (TMappoint*) j )->Pos;

    if ( c0 > c1 )  return 1;
    if ( c0 < c1 )  return -1;
    if ( p0 > p1 )  return 1;
    if ( p0 < p1 )  return -1;
    return 0;
}


void generate_map_file( fName, Description )
char**fName;
char**Description;

{
    int status = 0, i;
    char*token;
    int number = 1;
    char name[ 3 * MAXNAME ], classname[ 3 * MAXNAME ];
    TMappoint*newMPs = NULL;
    int cMPs = 0;
    FILE*fp = fileopen( *fName, "w" );

    token = strtok( *Description, "[' ',',']" );

    while ( NULL != token )
    {
        switch ( status )
        {

        case 0:

            if ( 1 == sscanf( token, "*%d", &number ) )
            {
                status = 1;
            }

            else
            {
                strcpy( name, token );
                number = 1;
                status = 2;
            }

            break;

        case 1:
            strcpy( name, token );
            status = 2;
            break;

        case 2:
            strcpy( classname, token );
            status = 3;
            break;

        case 3:

            if ( 0 == strcmp( token, "random" ) )
            {
                int lNR = 0;
                double genomelength = 0;
                status = 0;

                for ( i = 0;i < NoChroms;i++ )
                    genomelength += PopDescr[ i ].ChromLen;

                for ( i = 0;i < number;i++ )
                {
                    int chrom = -1;
                    CPOS pos = - gsl_rng_uniform_pos(rng[0]) * genomelength;

                    while ( pos < 0 )
                        pos += PopDescr[ ++chrom ].ChromLen;

                    {
                        char tmpstring[ 4 * MAXNAME ];

                        if ( 0 == cMPs++ )
                            newMPs = ( TMappoint* ) xCalloc( 1, TMappoint );
                        else
                            newMPs = ( TMappoint* ) xRealloc( newMPs, cMPs, TMappoint );

                        newMPs[ cMPs - 1 ].chrom = chrom;

                        newMPs[ cMPs - 1 ].Pos = pos;

                        newMPs[ cMPs - 1 ].counter = cMPs;

                        if ( lNR != -1 )
                        {
                            sprintf( tmpstring, "%s%d", name, lNR++ );
                        }

                        else
                            sprintf( tmpstring, "%s", name );

                        newMPs[ cMPs - 1 ].name = ( char* ) xCalloc( strlen( tmpstring ) + 1, char );

                        strcpy( newMPs[ cMPs - 1 ].name, tmpstring );

                        newMPs[ cMPs - 1 ].class = ( char* ) xCalloc( strlen( classname ) + 1, char );

                        strcpy( newMPs[ cMPs - 1 ].class, classname );
                    }

                    ;
                }

                number = 1;
            }

            else if ( ( 0 == strcmp( token, "terminal" ) ) && ( number % NoChroms == 0 ) )
            {
                int chrom, j, lNR = 0;
                CPOS pos, delta;

                status = 0;

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    delta = ( number / NoChroms > 1 ) ? PopDescr[ chrom ].ChromLen / ( number / NoChroms - 1 ) : 0;

                    for ( j = 0;j < number / NoChroms;j++ )
                    {
                        pos = delta * j;


                        {
                            char tmpstring[ 4 * MAXNAME ];

                            if ( 0 == cMPs++ )
                                newMPs = ( TMappoint* ) xCalloc( 1, TMappoint );
                            else
                                newMPs = ( TMappoint* ) xRealloc( newMPs, cMPs, TMappoint );

                            newMPs[ cMPs - 1 ].chrom = chrom;

                            newMPs[ cMPs - 1 ].Pos = pos;

                            newMPs[ cMPs - 1 ].counter = cMPs;

                            if ( lNR != -1 )
                            {
                                sprintf( tmpstring, "%s%d", name, lNR++ );
                            }

                            else
                                sprintf( tmpstring, "%s", name );

                            newMPs[ cMPs - 1 ].name = ( char* ) xCalloc( strlen( tmpstring ) + 1, char );

                            strcpy( newMPs[ cMPs - 1 ].name, tmpstring );

                            newMPs[ cMPs - 1 ].class = ( char* ) xCalloc( strlen( classname ) + 1, char );

                            strcpy( newMPs[ cMPs - 1 ].class, classname );
                        }

                        ;
                    }
                }

                number = 1;
            }

            else if ( ( 0 == strcmp( token, "nonterminal" ) ) && ( number % NoChroms == 0 ) )
            {
                int chrom, j, lNR = 0;
                CPOS pos, delta;

                status = 0;

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    delta = PopDescr[ chrom ].ChromLen / ( number / NoChroms + 1 );

                    for ( j = 1;j <= number / NoChroms;j++ )
                    {
                        pos = delta * j;


                        {
                            char tmpstring[ 4 * MAXNAME ];

                            if ( 0 == cMPs++ )
                                newMPs = ( TMappoint* ) xCalloc( 1, TMappoint );
                            else
                                newMPs = ( TMappoint* ) xRealloc( newMPs, cMPs, TMappoint );

                            newMPs[ cMPs - 1 ].chrom = chrom;

                            newMPs[ cMPs - 1 ].Pos = pos;

                            newMPs[ cMPs - 1 ].counter = cMPs;

                            if ( lNR != -1 )
                            {
                                sprintf( tmpstring, "%s%d", name, lNR++ );
                            }

                            else
                                sprintf( tmpstring, "%s", name );

                            newMPs[ cMPs - 1 ].name = ( char* ) xCalloc( strlen( tmpstring ) + 1, char );

                            strcpy( newMPs[ cMPs - 1 ].name, tmpstring );

                            newMPs[ cMPs - 1 ].class = ( char* ) xCalloc( strlen( classname ) + 1, char );

                            strcpy( newMPs[ cMPs - 1 ].class, classname );
                        }

                        ;
                    }
                }

                number = 1;

            }

            else
            {
                int chrom = 0;
                CPOS pos = 0;

                if ( ( 2 == sscanf( token, "%d/%lf", &chrom, &pos ) ) && ( 1 == number ) )
                {
                    int lNR = -1;
                    status = 0;
                    chrom--;

                    {
                        char tmpstring[ 4 * MAXNAME ];

                        if ( 0 == cMPs++ )
                            newMPs = ( TMappoint* ) xCalloc( 1, TMappoint );
                        else
                            newMPs = ( TMappoint* ) xRealloc( newMPs, cMPs, TMappoint );

                        newMPs[ cMPs - 1 ].chrom = chrom;

                        newMPs[ cMPs - 1 ].Pos = pos;

                        newMPs[ cMPs - 1 ].counter = cMPs;

                        if ( lNR != -1 )
                        {
                            sprintf( tmpstring, "%s%d", name, lNR++ );
                        }

                        else
                            sprintf( tmpstring, "%s", name );

                        newMPs[ cMPs - 1 ].name = ( char* ) xCalloc( strlen( tmpstring ) + 1, char );

                        strcpy( newMPs[ cMPs - 1 ].name, tmpstring );

                        newMPs[ cMPs - 1 ].class = ( char* ) xCalloc( strlen( classname ) + 1, char );

                        strcpy( newMPs[ cMPs - 1 ].class, classname );
                    }

                    ;
                }

                else {
                    info( -2, "Error in input file %s\n", token );
		    return;
		}

            }
        }

        token = strtok( NULL, "[',',' ']" );
    }

    qsort( newMPs, cMPs, sizeof( TMappoint ), mccmp );

    for ( i = 0;i < cMPs;i++ )
    {
        fprintf( fp, "%d\t%f\t%s\t%s\n", newMPs[ i ].chrom + 1, newMPs[ i ].Pos,
                 newMPs[ i ].name, newMPs[ i ].class );
        xFree( newMPs[ i ].name );
        xFree( newMPs[ i ].class );
    }

    xFree( newMPs );
    fileclose( fp );
}


void generate_effect_file( fName, Description )
char**fName;
char**Description;

{

    int status = 0;
    char*token, classname[ MAXNAME * 3 ], all0[ MAXNAME * 3 ], all1[ MAXNAME * 3 ];

    FILE*fp = fileopen( *fName, "w" );

    fprintf( fp, "# Automatisch erzeugt\n" );
    fprintf( fp, "# %s\n", *Description );

    token = strtok( *Description, "[' ',',']" );

    while ( NULL != token )
    {
        switch ( status )
        {

        case 0:
            {
                int tmpi;

                if ( 1 == sscanf( token, "%d", &tmpi ) )
                {
                    fprintf( fp, "%d\n", tmpi );

                }

                else
                {
                    strcpy( classname, token );
                    strcpy( all0, "" );
                    strcpy( all1, "" );
                    status = 2;
                }

                break;
            }

        case 2:
            strcpy( all0, token );
            status = 4;
            break;

        case 3:
            {
                float mean, var;
                int chrom, i;

                if ( 2 == sscanf( token, "%f/%f", &mean, &var ) )
                {
                    for ( chrom = 0;chrom < NoChroms;chrom++ )
                    {
                        for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
                        {
                            if ( NULL != strstr( PopDescr[ chrom ].Mappoint[ i ].class, classname ) )
                            {
                                if ( 0 == strcmp( all1, "" ) )
                                    fprintf( fp, "%s\t%s\t%f\n", PopDescr[ chrom ].Mappoint[ i ].name, all0, norm_rand() * sqrt( var ) + mean );
                                else
                                    fprintf( fp, "%s\t%s\t%s\t%s\t%f\n",
                                             PopDescr[ chrom ].Mappoint[ i ].name, all0,
                                             PopDescr[ chrom ].Mappoint[ i ].name, all1,
                                             norm_rand() * sqrt( var ) + mean );
                            }
                        }
                    }

                    status = 0;
                }

                else 
		  {
		    info( -2, "Error\n" );
		    return;
		  }
                break;
            }

        case 4:

            if ( 0 == strcmp( token, "normal" ) )
                status = 3;
            else if ( 0 == strcmp( token, "uniform" ) )
                status = 6;
            else
            {
                strcpy( all1, token );
                status = 5;
            }

            break;

        case 5:

            if ( 0 == strcmp( token, "normal" ) )
                status = 3;
            else if ( 0 == strcmp( token, "uniform" ) )
                status = 6;
            else {
                info( -2, "Fehler Verteilungsart erwartet: %s\n", token );
		return;
	    }

            break;

        case 6:
            {
                float eff;
                int chrom, i;

                if ( 1 == sscanf( token, "%f", &eff ) )
                {
                    for ( chrom = 0;chrom < NoChroms;chrom++ )
                    {
                        for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
                        {
                            if ( NULL != strstr( PopDescr[ chrom ].Mappoint[ i ].class, classname ) )
                            {
                                if ( 0 == strcmp( all1, "" ) )
                                    fprintf( fp, "%s\t%s\t%f\n", PopDescr[ chrom ].Mappoint[ i ].name, all0, eff );
                                else
                                    fprintf( fp, "%s\t%s\t%s\t%s\t%f\n",
                                             PopDescr[ chrom ].Mappoint[ i ].name, all0,
                                             PopDescr[ chrom ].Mappoint[ i ].name, all1,
                                             eff );
                            }
                        }
                    }

                    status = 0;
                }

                else {
                    info( -2, "Error\n" );
		    return;
		}
                break;
            }
        }

        token = strtok( NULL, "[',',' ']" );
    }

    fileclose( fp );
}


unsigned short get_allele( Hom, Pos )
THom*Hom;
CPOS Pos;
{
    TLoc*start = Hom->Loc, *end = &( Hom->Loc[ Hom->NoLociUsed - 1 ] ), *pos = start;

    while ( start <= end )
      {
        pos = ( end - start ) / 2 + start;
        if ( Pos < pos->Pos )
	  {
            end = pos - 1;
	  }
        else
	  {
            if ( ( pos + 1 ) ->Pos <= Pos )
	      {
                start = pos + 1;
	      }
            else
	      {
                return pos->All;
	      }
	  }
      }

    return pos->All;
}


IntNode* insertIntLinkList( 
			   IntNode* Genotype,
			   int*     Phaenotyp,
			   int      counter,
			   int*     nrGt )
{
  IntNode*nodepos = Genotype;
  
  if ( NULL == Genotype )
    {
      Genotype = ( IntNode* ) xCalloc( 1, IntNode );
      Genotype->child = NULL;
      Genotype->parent = NULL;
      Genotype->prev = NULL;
      Genotype->next = NULL;

      if ( 0 == counter )
        {
	  Genotype->number = 1;
	  ( *nrGt ) ++;
        }
      else
        {
	  Genotype->number = Phaenotyp[ 0 ];
	  Genotype->child = insertIntLinkList( Genotype->child, &Phaenotyp[ 1 ], counter - 1, nrGt );
	  Genotype->child->parent = Genotype;
        }

        return Genotype;
    }

  if ( 0 == counter )
    {
      ( Genotype->number ) ++;
      return Genotype;
    }
  do
    {
      if ( nodepos->number != Phaenotyp[ 0 ] )
        {
	  if ( nodepos->number > Phaenotyp[ 0 ] )
            {
	      
	      IntNode * tmp;
	      tmp = ( IntNode* ) xCalloc( 1, IntNode );
	      tmp->child = NULL;
	      tmp->parent = nodepos->parent;
	      tmp->prev = NULL;
	      tmp->next = nodepos;
	      nodepos->prev = tmp;
	      
	      if ( nodepos->parent != NULL )
		nodepos->parent->child = tmp;
	      
	      tmp->number = Phaenotyp[ 0 ];
	      
	      Genotype = nodepos = tmp;
            }
	  else
            {
	      if ( ( NULL != nodepos->next ) && ( nodepos->next->number <= Phaenotyp[ 0 ] ) )
                {
		  nodepos = nodepos->next;
                }
	      else
                {
		  IntNode*tmp;
		  tmp = ( IntNode* ) xCalloc( 1, IntNode );
		  tmp->child = NULL;
		  tmp->parent = nodepos->parent;
		  tmp->prev = nodepos;
		  tmp->next = nodepos->next;
		  tmp->number = Phaenotyp[ 0 ];
		  nodepos->next = tmp;
		  
		  if ( tmp->next != NULL )
                    {
		      tmp->next->prev = tmp;
                    }
		  
		  nodepos = tmp;
                }
            }
        }
    }

  while ( nodepos->number < Phaenotyp[ 0 ] );

  nodepos->child = insertIntLinkList( nodepos->child,
				      &Phaenotyp[ 1 ],
				      counter - 1,
				      nrGt );
  
  nodepos->child->parent = nodepos;
  
  return Genotype;

}

IntNode* go_down ( IntNode *Node )
{
  while ( Node->child != NULL ) Node = Node->child;
  return Node;
}

void print_path( IntNode*Node,
		 int*_nrLoci,
		 int*_nrGt,
		 int*m )
{
    int i = *_nrLoci;

    do {
      m[ ( i-- ) * ( *_nrGt ) ] = Node->number;
    }
    while ( NULL != ( Node = Node->parent ) );
}


IntNode* find_next_leaf( IntNode* Node)
{
    IntNode*tmp;

    while ( NULL != Node->parent )
      {
        tmp = Node;
        Node = Node->parent;
        xFree( tmp );
	
        if ( Node->next != NULL )
	  {
            tmp = Node->next;
            xFree( Node );
            return ( go_down( tmp ) );
	  }
      }

    xFree( Node );
    return ( NULL );
}


void print_genotype ( int *_nrLoci,
		      int *_nrGt,
		      int *_matrix,
		      int *mode,
		      char **Loci,
		      char **lName,
		      char **separator,
		      char **tail  )
{
    int     tmp_counter = 0;
    int     hom;
    IntNode *Node = startNode;
    char    *pch;
    char    temp[100];

    pch = strtok( *Loci, " " );
    strcpy( *lName, "" );

    while ( pch != NULL )
    {
      if ( 1 == *mode ) {
	for ( hom = 1;hom < NoHoms + 1;hom++ ) {
	  sprintf( temp, "%s - %d%s", pch, hom, *separator );
	  strcat( *lName, temp );
	}
      }
      else {
	sprintf( temp, "%s%s", pch, *separator );
	strcat( *lName, temp );
      }
      pch = strtok( NULL, " " );
    }

    strcat( *lName, *tail );

    if ( NULL != ( Node = go_down( Node ) ) ) {
      do {
	print_path( Node, _nrLoci, _nrGt, &_matrix[ tmp_counter++ ] );
      }
      while ( NULL != ( Node = find_next_leaf( Node ) ) );
    }
}


SEXP genotype_evaluate( SEXP _PopName,
			SEXP _Loci,
			SEXP _mode )
{
    TPopListEl*PopListEl = NULL;
    int mode = 0;
    char*PopName;

    if ( 0 == cMP ){
      info( -2, "No linkage map defined\n" );
      return R_NilValue;
    }

    if ( !isInteger( _mode ) || ( length( _mode ) != 1 ) ) 
      {
        info( -2, "Value of mode is not an integer\n" );
	return R_NilValue;
      }
    else
      {
        mode = INTEGER( _mode ) [ 0 ];
      }

    PopName = ( char* ) xCalloc( strlen( CHAR( STRING_ELT( _PopName, 0 ) ) ) + 2, char );

    if ( !isString( _PopName ) || ( length( _PopName ) != 1 ) )
      {
        info( -2, "Name of population is not a single string\n" );
	return R_NilValue;
      }
    else
      {
        strcpy( PopName, CHAR( STRING_ELT( _PopName, 0 ) ) );
      }

    if ( !isString( _Loci ) )
      {
        info( -2, "Name of Loci is not a string\n" );
	return R_NilValue;
      }

    find_list_element( PopName, &PopListEl );

    if ( NULL == PopListEl )
      {
	info( -2, "Population in not defined: %s\n", PopName );
	return R_NilValue ;
      }

    else
    {

        typedef struct
        {
            int nr_marker;
            TMappoint*Mappoint;
        }

        TEvChrom;

        TEvChrom*EvChrom = ( TEvChrom* ) xCalloc( NoChroms, TEvChrom );
        int chrom, i, j, l, m;

        int nrLoci = 0;
        int*pnrLoci = &nrLoci;
        int nrGt = 0;
        int*pnrGt = &nrGt;

        SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;
        int*locus_idx = ( int* ) xCalloc( length( _Loci ), int );

        POPULATION ind;
        CPOS Pos;

        startNode = NULL;


        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            EvChrom[ chrom ].nr_marker = 0;
            EvChrom[ chrom ].Mappoint = NULL;
        }

        for ( l = 0;l < length( _Loci );l++ )
        {
            locus_idx[ l ] = 0;

            for ( chrom = 0;( ( chrom < NoChroms ) && ( 0 == locus_idx[ l ] ) );chrom++ )
            {
                for ( j = 0;( ( j < PopDescr[ chrom ].NrMappoints ) && ( 0 == locus_idx[ l ] ) );j++ )
                {
                    if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ j ].name, CHAR( STRING_ELT( _Loci, l ) ) ) )
                    {
                        locus_idx[ l ] = 1;
                        nrLoci += NoHoms;
                        EvChrom[ chrom ].nr_marker++;

                        if ( EvChrom[ chrom ].nr_marker == 1 )
                        {
                            EvChrom[ chrom ].Mappoint = ( TMappoint* ) xCalloc( 1, TMappoint );
                        }

                        else
                        {
                            EvChrom[ chrom ].Mappoint = ( TMappoint* ) xRealloc( EvChrom[ chrom ].Mappoint,
                                                        EvChrom[ chrom ].nr_marker, TMappoint );
                        }

                        EvChrom[ chrom ].Mappoint[ EvChrom[ chrom ].nr_marker - 1 ] = PopDescr[ chrom ].Mappoint[ j ];
                    }
                }
            }
        }

        if ( 0 == nrLoci )
	  {
            info( -2, "Loci not found\n" );
	    return R_NilValue;
	  }

        if ( 0 == PopListEl->Pop->NoInds )
	  {
            info( -2, "Empty population: %s\n", *PopName );
	    return R_NilValue;
	  }

        {
            TChrom*Chrom;
            TInd*Ind;
            IntNode*Node;
            int*alleles = ( int* ) xCalloc( NoHoms, int );
            int*IndGT = ( int* ) xCalloc( nrLoci + 1, int );
            int hom, tmp_counter = 0;

            for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
                Ind = get_individual( PopListEl->Pop, ind );
                m = 0;

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    for ( l = 0;l < EvChrom[ chrom ].nr_marker;l++ )
                    {
                        if ( NULL == Ind->genotype )
                        {
                            Chrom = &( Ind->Chrom[ chrom ] );
                            Pos = EvChrom[ chrom ].Mappoint[ l ].Pos;

                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                alleles[ hom ] = get_allele( &( Chrom->Hom[ hom ] ), Pos );
                            }
                        }

                        else
                        {
                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                alleles[ hom ] = Ind->genotype[ 2 * EvChrom[ chrom ].Mappoint[ l ].counter + hom ];
                            }
                        }

                        switch ( mode )
                        {

                        case 0:
                            break;

                        case 1:

                            if ( EvChrom[ chrom ].nr_marker == 1 )
                                R_isort( &alleles[ 0 ], NoHoms );

                            break;

                        default:
                            R_isort( &alleles[ 0 ], NoHoms );

                            break;
                        }

                        for ( hom = 0;hom < NoHoms;hom++ )
                        {
                            IndGT[ m * NoHoms + hom ] = alleles[ hom ];
                        }

                        m++;
                    }
                }

                startNode = insertIntLinkList( startNode, IndGT, nrLoci, pnrGt );
            }

            {
                int*_matrix = ( int* ) xCalloc( ( nrLoci + 1 ) * nrGt, int );
                int offset;
                char tmp[ 1024 ];
                int N = 0;

                Node = startNode;

                if ( NULL != ( Node = go_down( Node ) ) )
                {
                    do
                    {
                        print_path( Node, pnrLoci, pnrGt, &_matrix[ tmp_counter++ ] );
                    }

                    while ( NULL != ( Node = find_next_leaf( Node ) ) );
                }

                PROTECT( matrix = allocMatrix( INTSXP, nrGt, nrLoci + 2 ) );
                PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
                INTEGER( matrixDim ) [ 0 ] = nrGt;
                INTEGER( matrixDim ) [ 1 ] = nrLoci + 2;
                setAttrib( matrix, R_DimSymbol, matrixDim );
                PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
                PROTECT( xAxis = allocVector( STRSXP, nrLoci + 2 ) );
                PROTECT( yAxis = allocVector( STRSXP, nrGt ) );

                for ( i = offset = 0;i < nrLoci + 1;i++ )
                {
                    for ( j = 0;j < nrGt;j++, offset++ )
                    {
                        INTEGER( matrix ) [ offset ] = _matrix[ offset ];

                        if ( i == nrLoci )
                        {
                            N += _matrix[ offset ];
                        }
                    }
                }

                for ( j = 0;j < nrGt;j++, offset++ )
                {
                    INTEGER( matrix ) [ offset ] = N;
                }

                for ( i = 0, offset = 0;i < length( _Loci );i++ )
                {
                    if ( locus_idx[ i ] )
                    {
                        for ( hom = 1;hom <= NoHoms;hom++ )
                        {
                            sprintf( tmp, "%s-%d ", CHAR( STRING_ELT( _Loci, i ) ), hom );
                            SET_STRING_ELT( xAxis, offset++, mkChar( tmp ) );
                        }
                    }
                }

                SET_STRING_ELT( xAxis, offset++, mkChar( "count" ) );
                SET_STRING_ELT( xAxis, offset++, mkChar( "frequency" ) );

                for ( i = 0;i < nrGt;i++ )
                {
                    sprintf( tmp, "%d ", i + 1 );
                    SET_STRING_ELT( yAxis, i, mkChar( tmp ) );
                }

                SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
                SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
                setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
                UNPROTECT( 5 );
                xFree( PopName );
                xFree( EvChrom );
                xFree( locus_idx );
                xFree( alleles );
                xFree( IndGT );
                xFree( _matrix );
                return ( matrix );
            }
        }
    }
    return(R_NilValue);
}


void evaluate_genotype( PopName,
                        Loci,
                        mode,
                        nrGt,
                        nrLoci,
                        llName,
                        separator,
                        tail )
char**PopName  ;
char**Loci;
int*mode;
int*nrGt;
int*nrLoci;
int*llName;
char**separator;
char**tail;
{

    typedef struct
    {
        int nr_marker;
        TMappoint*Mappoint;
    }

    TEvChrom;

    TEvChrom*EvChrom = ( TEvChrom* ) xCalloc( NoChroms, TEvChrom );
    CPOS Pos;
    int j, chrom;
    TPopListEl*PopListEl = NULL;
    startNode = NULL;
    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Can't find population %s!", *PopName );
	return;
    }

    if ( 0 == cMP ) {
        info( -2, "No linkage map defined\n" );
	return;
    }


    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        EvChrom[ chrom ].nr_marker = 0;
        ( &EvChrom[ chrom ] ) ->Mappoint = NULL;
    }

    {
        char*pch;
        *llName += strlen( *tail ) + 5;
        pch = strtok( *Loci, " " );

        while ( pch != NULL )
        {

            for ( chrom = 0;chrom < NoChroms;chrom++ )
            {
                if ( PopDescr[ chrom ].NrMappoints > 0 )
                {
                    for ( j = 0;j < PopDescr[ chrom ].NrMappoints;j++ )
                    {
                        if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ j ].name, pch ) )
                        {
                            ( *nrLoci ) += NoHoms;
                            EvChrom[ chrom ].nr_marker++;
                            *llName += ( ( strlen( pch ) + 6 + strlen( *separator ) ) * NoHoms );

                            if ( EvChrom[ chrom ].nr_marker == 1 )
                            {
                                EvChrom[ chrom ].Mappoint = ( TMappoint* ) xCalloc( 1, TMappoint );
                            }

                            else
                            {
                                EvChrom[ chrom ].Mappoint = ( TMappoint* ) xRealloc( EvChrom[ chrom ].Mappoint,
                                                            EvChrom[ chrom ].nr_marker,
                                                            TMappoint );
                            }

                            EvChrom[ chrom ].Mappoint[ EvChrom[ chrom ].nr_marker - 1 ] = ( PopDescr[ chrom ].Mappoint[ j ] );
                            break;
                        }
                    }
                }
            }

            pch = strtok( NULL, " " );
        }
    }

    strcpy( *Loci, "" );

    if ( 0 < *nrLoci )
    {
        int k, l, m;
        int*all = ( int* ) xCalloc( NoHoms, int );
        int*IndGT = ( int* ) xCalloc( nrLoci + 1, int );
        POPULATION ind;
        IndGT[ *nrLoci ] = -1;

        if ( 0 == PopListEl->Pop->NoInds ){
	  info( -2, "Empty population %s!", *PopName );
	  return;
	}

        for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
        {
            TInd*Ind = get_individual( PopListEl->Pop, ind );
            m = -1;

            for ( k = 0;k < NoChroms;k++ )
            {
                if ( EvChrom[ k ].nr_marker > 0 )
                {
                    for ( l = 0;l < EvChrom[ k ].nr_marker;l++ )
                    {
                        TChrom*Chrom = &( Ind->Chrom[ k ] );
                        m++;

                        if ( 0 == ind )
                        {
                            if ( 0 != strlen( *Loci ) )
                                strcat( *Loci, " " );

                            strcat( *Loci, EvChrom[ k ].Mappoint[ l ].name );
                        }

                        Pos = ( &( &EvChrom[ k ] ) ->Mappoint[ l ] ) ->Pos;

                        if ( NULL == Ind->genotype )
                        {
                            for ( j = 0;j < NoHoms;j++ )
                            {
                                all[ j ] = get_allele( &( Chrom->Hom[ j ] ), Pos );
                            }
                        }

                        else
                        {
                            for ( j = 0;j < NoHoms;j++ )
                            {
                                all[ j ] = Ind->genotype[ 2 * EvChrom[ k ].Mappoint[ l ].counter + j ];
                            }
                        }

                        switch ( *mode )
                        {

                        case 0:
                            break;

                        case 1:

                            if ( EvChrom[ k ].nr_marker == 1 )
                                R_isort( &all[ 0 ], NoHoms );

                            break;

                        default:
                            R_isort( &all[ 0 ], NoHoms );

                            break;
                        }

                        for ( j = 0;j < NoHoms;j++ )
                        {
                            IndGT[ m * NoHoms + j ] = all[ j ];
                        }
                    }
                }
            }

            startNode = insertIntLinkList( startNode, IndGT, *nrLoci, nrGt );
        }


        xFree( all );
        xFree( IndGT );

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            if ( 0 < EvChrom[ chrom ].nr_marker )
            {
                xFree( ( &EvChrom[ chrom ] ) ->Mappoint );
            }
        }
    }

    else
    {
        info( -2, "Loci not found\n" );
	return;
    }

    xFree( EvChrom );
}


SEXP  allele_evaluate( _PopName, _Loci, _missing, _bg )
SEXP _PopName;
SEXP _Loci;
SEXP _missing;
SEXP _bg;

{
    char*PopName;
    int bg = 0;
    int missing = 0;
    TPopListEl*PopListEl = NULL;

    if ( 0 == cMP ){
        info( -2, "No linkage map defined\n" );
	return R_NilValue;
    }

    PopName = ( char* ) xCalloc( strlen( CHAR( STRING_ELT( _PopName, 0 ) ) ) + 2, char );

    if ( !isString( _PopName ) || ( length( _PopName ) != 1 ) ) {
        info( -2, "Name of population is not a single string\n" );
	return R_NilValue;
    }
    else
        strcpy( PopName, CHAR( STRING_ELT( _PopName, 0 ) ) );

    if ( !isString( _Loci ) ){
        info( -2, "Name of Loci is not a string\n" );
	return R_NilValue;
    }

    if ( !isInteger( _bg ) || ( length( _bg ) != 1 ) ) {
        info( -2, "Value of background allele is not an integer\n" );
	return R_NilValue;
    }
    else
        bg = INTEGER( _bg ) [ 0 ];

    if ( !isInteger( _missing ) || ( length( _missing ) != 1 ) ) {
        info( -2, "Value of missing is not an integer\n" );
	return R_NilValue;
    }
    else
        missing = INTEGER( _missing ) [ 0 ];

    find_list_element( PopName, &PopListEl );

    if ( NULL == PopListEl )
    {
        info( -2, "Can't find population %s\n", PopName );
        xFree( PopName );
	return R_NilValue;
    }

    else
    {

        int*locus_idx = ( int* ) xCalloc( length( _Loci ), int );

        int nrGt = 0;
        int*pnrGt = &nrGt;

        typedef struct
        {
            int nr_marker;
            TMappoint*Mappoint;
        }

        TEvChrom;

        TEvChrom*EvChrom = ( TEvChrom* ) xCalloc( NoChroms, TEvChrom );
        int j, l, chrom;
        int nrLoci = 0;
        int*pnrLoci = &nrLoci;

        startNode = NULL;


        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            EvChrom[ chrom ].nr_marker = 0;
            EvChrom[ chrom ].Mappoint = NULL;
        }


        for ( l = 0;l < length( _Loci );l++ )
        {
            locus_idx[ l ] = 0;

            for ( chrom = 0;( ( chrom < NoChroms ) && ( 0 == locus_idx[ l ] ) );chrom++ )
            {
                for ( j = 0;( ( j < PopDescr[ chrom ].NrMappoints ) && ( 0 == locus_idx[ l ] ) );j++ )
                {
                    if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ j ].name, CHAR( STRING_ELT( _Loci, l ) ) ) )
                    {
                        locus_idx[ l ] = 1;
                        nrLoci++;
                        EvChrom[ chrom ].nr_marker++;

                        if ( EvChrom[ chrom ].nr_marker == 1 )
                        {
                            EvChrom[ chrom ].Mappoint = ( TMappoint* ) xCalloc( 1, TMappoint );
                        }

                        else
                        {
                            EvChrom[ chrom ].Mappoint = ( TMappoint* ) xRealloc( EvChrom[ chrom ].Mappoint,
                                                        EvChrom[ chrom ].nr_marker, TMappoint );
                        }

                        EvChrom[ chrom ].Mappoint[ EvChrom[ chrom ].nr_marker - 1 ] = PopDescr[ chrom ].Mappoint[ j ];
                    }
                }
            }
        }


        if ( 0 == nrLoci )
        {
            info( -2, "Loci not found\n" );
	    return R_NilValue;
        } else

        if ( 0 == PopListEl->Pop->NoInds )
        {
            info( -2, "Empty population %s\n", *PopName );
	    return R_NilValue;
        } else

        {
            CPOS Pos;
            int chrom, hom, l, temp, i;
            int*all = ( int* ) xCalloc( nrLoci, int );
            POPULATION ind;
            TInd*Ind;
            TChrom*Chrom;

            for ( i = 0;i < nrLoci;i++ )
            {
                all[ i ] = 0;
            }

            for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
                Ind = get_individual( PopListEl->Pop, ind );
                temp = 0;

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    for ( l = 0;l < EvChrom[ chrom ].nr_marker;l++ )
                    {
                        if ( NULL == Ind->genotype )
                        {
                            Chrom = &( Ind->Chrom[ chrom ] );
                            Pos = EvChrom[ chrom ].Mappoint[ l ].Pos;

                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                if ( 0 == ( all[ temp ] = get_allele( &( Chrom->Hom[ hom ] ), Pos ) ) )
                                    info( -2, "Allele 0 not allowed\n" );

                                if ( ( 1 == missing ) && ( all[ temp ] == bg ) )
                                {
                                    all[ temp ] = 0;
                                }

                                else
                                {
                                    startNode = insertIntLinkList( startNode, all, nrLoci, pnrGt );
                                }
                            }
                        }

                        else
                        {
                            for ( hom = 0;hom < NoHoms;hom++ )
                            {
                                if ( 0 == ( all[ temp ] = Ind->genotype[ 2 * EvChrom[ chrom ].Mappoint[ l ].counter + hom ] ) )
                                    info( -2, "Allele 0 not allowed\n" );

                                if ( ( 1 == missing ) && ( all[ temp ] == bg ) )
                                {
                                    all[ temp ] = 0;
                                }

                                else
                                {
                                    startNode = insertIntLinkList( startNode, all, nrLoci, pnrGt );
                                }
                            }
                        }

                        all[ temp ] = 0;
                        temp++;
                    }
                }
            }


            for ( chrom = 0;chrom < NoChroms;chrom++ )
            {
                if ( 0 < EvChrom[ chrom ].nr_marker )
                {
                    xFree( ( &EvChrom[ chrom ] ) ->Mappoint );
                }
            }

            {
                int*_matrix = ( int* ) xCalloc( ( nrLoci + 1 ) * nrGt, int );
                int offset;
                char tmp[ 1024 ];
                int N = 0;
                IntNode*Node;
                int tmp_counter = 0;
                SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;


                Node = startNode;

                if ( NULL != ( Node = go_down( Node ) ) )
                {
                    do
                    {
                        print_path( Node, pnrLoci, pnrGt, &_matrix[ tmp_counter++ ] );
                    }

                    while ( NULL != ( Node = find_next_leaf( Node ) ) );
                }

                PROTECT( matrix = allocMatrix( INTSXP, nrGt, nrLoci + 2 ) );
                PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
                INTEGER( matrixDim ) [ 0 ] = nrGt;
                INTEGER( matrixDim ) [ 1 ] = nrLoci + 2;
                setAttrib( matrix, R_DimSymbol, matrixDim );
                PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
                PROTECT( xAxis = allocVector( STRSXP, nrLoci + 2 ) );
                PROTECT( yAxis = allocVector( STRSXP, nrGt ) );

                for ( i = offset = 0;i < nrLoci + 1;i++ )
                {
                    for ( j = 0;j < nrGt;j++, offset++ )
                    {
                        INTEGER( matrix ) [ offset ] = _matrix[ offset ];

                        if ( i == nrLoci )
                        {
                            N += _matrix[ offset ];
                        }
                    }
                }

                for ( j = 0;j < nrGt;j++, offset++ )
                {
                    INTEGER( matrix ) [ offset ] = N;
                }

                for ( i = 0, offset = 0;i < length( _Loci );i++ )
                {
                    if ( locus_idx[ i ] )
                    {
                        sprintf( tmp, "%s", CHAR( STRING_ELT( _Loci, i ) ) );
                        SET_STRING_ELT( xAxis, offset++, mkChar( tmp ) );
                    }
                }

                SET_STRING_ELT( xAxis, offset++, mkChar( "count" ) );
                SET_STRING_ELT( xAxis, offset++, mkChar( "frequency" ) );

                for ( i = 0;i < nrGt;i++ )
                {
                    sprintf( tmp, "%d ", i + 1 );
                    SET_STRING_ELT( yAxis, i, mkChar( tmp ) );
                }

                SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
                SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
                setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
                UNPROTECT( 5 );
                xFree( PopName );
                xFree( locus_idx );
                xFree( EvChrom );
                xFree( all );
                xFree( _matrix );
                return ( matrix );
            }
        }
    }
    return R_NilValue;
}



void evaluate_allele( PopName, Loci, nrGt, nrLoci, llName, separator, tail, missing, bg )
char**PopName;
char**Loci;
int*nrGt;
int*nrLoci;
int*llName;
char**separator;
char**tail;
int*missing;
int*bg;

{

    typedef struct
    {
        int nr_marker;
        TMappoint*Mappoint;
    }

    TEvChrom;

    TEvChrom*EvChrom;
    CPOS Pos;
    int chrom, j;
    TPopListEl*PopListEl;

    EvChrom = ( TEvChrom* ) xCalloc( NoChroms, TEvChrom );
    startNode = NULL;
    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Can't find population %s\n", *PopName );
	return;
    }

    if ( 0 == cMP ){
        info( -2, "No linkage map defined\n" );
	return;
    }


    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        EvChrom[ chrom ].nr_marker = 0;
        ( &EvChrom[ chrom ] ) ->Mappoint = NULL;
    }

    {
        char*pch;
        *llName += strlen( *tail ) + 5;
        pch = strtok( *Loci, " " );

        while ( pch != NULL )
        {



            for ( chrom = 0;chrom < NoChroms;chrom++ )
            {
                if ( PopDescr[ chrom ].NrMappoints > 0 )
                {
                    for ( j = 0;j < PopDescr[ chrom ].NrMappoints;j++ )
                    {
                        if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ j ].name, pch ) )
                        {
                            ( *nrLoci ) ++;
                            *llName += ( strlen( pch ) + strlen( *separator ) );
                            EvChrom[ chrom ].nr_marker++;

                            if ( 1 == EvChrom[ chrom ].nr_marker )
                                ( &EvChrom[ chrom ] ) ->Mappoint = ( TMappoint* ) xCalloc( 1, TMappoint );
                            else
                                ( &EvChrom[ chrom ] ) ->Mappoint = ( TMappoint* )
                                                                   xRealloc( ( &EvChrom[ chrom ] ) ->Mappoint, EvChrom[ chrom ].nr_marker, TMappoint );

                            ( &EvChrom[ chrom ] ) ->Mappoint[ EvChrom[ chrom ].nr_marker - 1 ] = ( PopDescr[ chrom ].Mappoint[ j ] );

                            break;
                        }
                    }
                }
            }


            pch = strtok( NULL, " " );
        }
    }

    strcpy( *Loci, "" );

    if ( 0 < *nrLoci )
    {
        int chrom, hom, l, temp, i;
        int*all;
        POPULATION ind;

        all = ( int* ) xCalloc( *nrLoci, int );

        for ( i = 0;i < *nrLoci;i++ )
            all[ i ] = 0;

        for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
        {
            temp = 0;

            for ( chrom = 0;chrom < NoChroms;chrom++ )
            {
                if ( EvChrom[ chrom ].nr_marker > 0 )
                {
                    for ( l = 0;l < EvChrom[ chrom ].nr_marker;l++ )
                    {
                        TInd*Ind = get_individual( PopListEl->Pop, ind );
                        TChrom*Chrom = &( Ind->Chrom[ chrom ] );

                        if ( 0 == ind )
                        {
                            if ( 0 != strlen( *Loci ) )
                                strcat( *Loci, " " );

                            strcat( *Loci, ( &( &EvChrom[ chrom ] ) ->Mappoint[ l ] ) ->name );
                        }

                        Pos = ( &( &EvChrom[ chrom ] ) ->Mappoint[ l ] ) ->Pos;

                        for ( hom = 0;hom < NoHoms;hom++ )
                        {
                            if ( 0 == ( all[ temp ] = get_allele( &( Chrom->Hom[ hom ] ), Pos ) ) )
                                info( -2, "Allele 0 not allowed\n" );

                            if ( ( 1 == *missing ) && ( all[ temp ] == *bg ) )
                            {
                                all[ temp ] = 0;
                            }

                            else
                            {
                                startNode = insertIntLinkList( startNode, all, *nrLoci, nrGt );
                            }
                        }

                        all[ temp ] = 0;
                        temp++;
                    }
                }
            }
        }


        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            if ( 0 < EvChrom[ chrom ].nr_marker )
            {
                xFree( ( &EvChrom[ chrom ] ) ->Mappoint );
            }
        }

        xFree( all );
    }

    else
    {
        info( -2, "Loci not found.\n" );
	return;
    }

    xFree( EvChrom );
}



void evaluate_locus( PopName, Loci, _matrix, noLoci, bg )
char**PopName;
char**Loci;
int*_matrix;
int*noLoci;
int*bg;

{
    TPopListEl*PopListEl;
    char*pch;
    int chrom, i = 0, j;
    POPULATION noHomo, noHetero, noMissing, ind;

    if ( 2 != NoHoms ){
        info( -2, "This function is only defined for diploid organisms!\n" );
	return;
    }


    if ( 0 == cMP ){
        info( -2, "No linkage map defined\n" );
 	return;
   }


    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Can't find population %s\n", *PopName );
	return;
    }

    pch = strtok( *Loci, " " );

    while ( pch != NULL )
    {
        noHomo = 0;
        noHetero = 0;
        noMissing = 0;

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            if ( PopDescr[ chrom ].NrMappoints > 0 )
            {
                for ( j = 0;j < PopDescr[ chrom ].NrMappoints;j++ )
                {
                    if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ j ].name, pch ) )
                    {
                        for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
                        {
                            TInd*Ind = get_individual( PopListEl->Pop, ind );

                            if ( NULL == Ind->genotype )
                            {
                                info( -2, "first genotype_population (chrom=%d, j=%d, ind=%d, pch=%s)",
                                      chrom, j, ind, pch );
				return;
                            }

                            if ( Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter ] !=
                                    Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter + 1 ] )
                            {
                                noHetero++;
                            }

                            else
                            {
                                if ( ( Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter ] == *bg ) ||
                                        ( Ind->genotype[ 2 * PopDescr[ chrom ].Mappoint[ j ].counter + 1 ] == *bg ) )
                                {
                                    noMissing++;
                                }

                                else
                                {
                                    noHomo++;
                                }
                            }
                        }
                    }
                }
            }
        }

        _matrix[ i ] = noHomo;
        _matrix[ i + ( *noLoci ) ] = noHetero;
        _matrix[ i + 2 * ( *noLoci ) ] = noMissing;
        i++;
        pch = strtok( NULL, " " );
    }
}



SEXP evaluate_haplotype( _PopName, _Loci, _missing, _bg )
SEXP _PopName;
SEXP _Loci;
SEXP _missing;
SEXP _bg;

{
    int bg = 0, missing = 0;
    SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;


    if ( !isString( _PopName ) || ( length( _PopName ) != 1 ) )
      {
        info( -2, "Name of population is not a single string\n" );
	return  R_NilValue;
      }


    if ( !isString( _Loci ) || ( length( _Loci ) != 1 ) )
      {
        info( -2, "Name of Loci is not a single string\n" );
	return R_NilValue;
      }

    if ( !isInteger( _bg ) || ( length( _bg ) != 1 ) )
      {
        info( -2, "Value of background allele is not an integer\n" );
	return  R_NilValue;
      }
    else
      {
        bg = INTEGER( _bg ) [ 0 ];
      }

    if ( !isInteger( _missing ) || ( length( _missing ) != 1 ) )
      {
        info( -2, "Value of missing variable is not an integer\n" );
	return R_NilValue;
      }
    else
      {
        missing = INTEGER( _missing ) [ 0 ];
      }

    {
        TPopListEl*PopListEl;
        POPULATION ind;
        TMappoint*Mappoint = NULL;

        int nLoci = 0;
        int nGt = 0;
        int*pGt = &nGt;

        char*PopName = ( char* ) xCalloc( strlen( CHAR( STRING_ELT( _PopName, 0 ) ) ) + 2, char );
        char*Loci = ( char* ) xCalloc( strlen( CHAR( STRING_ELT( _Loci, 0 ) ) ) + 2, char );

        strcpy( PopName, CHAR( STRING_ELT( _PopName, 0 ) ) );
        strcpy( Loci, CHAR( STRING_ELT( _Loci, 0 ) ) );

        find_list_element( PopName, &PopListEl );

        if ( !PopListEl ){
            info( -2, "Can't find population %s\n", PopName );
	    return R_NilValue;
	}

        if ( 1 > PopListEl->Pop->NoInds ){
            info( -2, "Empty population %s\n", PopName );
	    return R_NilValue;
	}

        {
            int lastLoci;
            int chrom, loc;
            char*locus;

            locus = strtok( Loci, " " );

            while ( locus != NULL )
            {
                lastLoci = nLoci;

                for ( chrom = 0;chrom < NoChroms;chrom++ )
                {
                    if ( PopDescr[ chrom ].NrMappoints > 0 )
                    {
                        for ( loc = 0;loc < PopDescr[ chrom ].NrMappoints;loc++ )
                        {
                            if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ loc ].name, locus ) )
                            {
                                if ( 0 == nLoci )
                                {
                                    Mappoint = ( TMappoint* ) xCalloc( 1, TMappoint );
                                }

                                else
                                {
                                    Mappoint = ( TMappoint* ) xRealloc( Mappoint, nLoci + 1, TMappoint );
                                }

                                Mappoint[ nLoci ] = ( PopDescr[ chrom ].Mappoint[ loc ] );
                                nLoci++;
                                break;
                            }
                        }
                    }
                }

                if ( lastLoci == nLoci )
                {
                    info( -2, "Cant find locus %s\n ", locus );
		    return R_NilValue;
                }

                locus = strtok( NULL, " " );
            }
        }


        {
            TInd*Ind;
            int loc, hom, nrmissing;
            int*hap = ( int* ) xCalloc( nLoci, int );


            startNode = NULL;

            for ( loc = 0;loc < nLoci;loc++ )
            {
                hap[ loc ] = 0;
            }

            startNode = insertIntLinkList( startNode, hap, nLoci, pGt );

            for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
                nrmissing = 0;
                Ind = get_individual( PopListEl->Pop, ind );

                if ( NULL == Ind->genotype ) {
                    info( -2, "First genotype_population\n" );
		    return R_NilValue;
		}

                for ( hom = 0;hom < NoHoms;hom++ )
                {
                    for ( loc = 0;loc < nLoci;loc++ )
                    {
                        hap[ loc ] = Ind->genotype[ Mappoint[ loc ].counter * 2 + hom ];

                        if ( ( 1 == missing ) && ( bg == hap[ loc ] ) )
                        {
                            nrmissing++;
                        }
                    }

                    if ( 0 == nrmissing )
                    {
                        startNode = insertIntLinkList( startNode, hap, nLoci, pGt );
                    }
                }
            }

            xFree( hap );
        }

        {
            int y_dim = nGt;
            int x_dim = nLoci + 1;
            IntNode*Node = startNode;
            int tmpcounter = 0;
            char*locus;

            PROTECT( matrix = allocMatrix( INTSXP, y_dim, x_dim ) );
            PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
            INTEGER( matrixDim ) [ 0 ] = y_dim;
            INTEGER( matrixDim ) [ 1 ] = x_dim;
            setAttrib( matrix, R_DimSymbol, matrixDim );
            PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
            PROTECT( xAxis = allocVector( STRSXP, x_dim ) );
            PROTECT( yAxis = allocVector( STRSXP, y_dim ) );

            if ( ( NULL != Node ) && ( NULL != ( Node = go_down( Node ) ) ) )
            {
                int * path = ( int* ) xCalloc( nGt * x_dim, int );
                int ttt = 0;
                int iii;

                do
                {
                    print_path( Node, &nLoci, &nGt, &path[ ttt++ ] );
                }

                while ( NULL != ( Node = find_next_leaf( Node ) ) );

                for ( iii = 0;iii < x_dim;iii++ )
                {
                    for ( ttt = 0;ttt < y_dim;ttt++, tmpcounter++ )
                    {
                        INTEGER( matrix ) [ tmpcounter ] = path[ tmpcounter ];
                    }
                }

                xFree( path );
            }

            strcpy( Loci, CHAR( STRING_ELT( _Loci, 0 ) ) );
            locus = strtok( Loci, " " );
            {
                int col = 0;

                while ( locus != NULL )
                {
                    SET_STRING_ELT( xAxis, col++, mkChar( locus ) );
                    locus = strtok( NULL, " " );
                }

                SET_STRING_ELT( xAxis, col, mkChar( "count" ) );
            }

            SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
            setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
            UNPROTECT( 5 );
        }

        {
            int x;
            int**table = NULL;
            int dimLA = 0, dimLB = 0;
            int lx, ly, xPos, yPos;

            for ( x = 0;x < nGt;x++ )
            {
                xPos = 0;
                yPos = 0;

                if ( ( 0 == dimLA ) && ( 0 == dimLB ) )
                {
                    dimLA = 1;
                    dimLB = 1;
                    table = ( int** ) xCalloc( dimLA + 1, int );
                    table[ 0 ] = ( int* ) xCalloc( dimLB + 1, int );
                    table[ 1 ] = ( int* ) xCalloc( dimLB + 1, int );
                    table[ 0 ][ 0 ] = INTEGER( matrix ) [ x + 2 * nGt ];
                    table[ 1 ][ 0 ] = INTEGER( matrix ) [ x ];
                    table[ 0 ][ 1 ] = INTEGER( matrix ) [ x + nGt ];
                    table[ 1 ][ 1 ] = INTEGER( matrix ) [ x + 2 * nGt ];
                }

                else
                {
                    for ( lx = 1;lx <= dimLB;lx++ )
                    {
                        if ( table[ 0 ][ lx ] == INTEGER( matrix ) [ x + nGt ] )
                        {
                            xPos = lx;
                        }
                    }

                    if ( 0 == xPos )
                    {
                        dimLB++;

                        for ( ly = 0;ly <= dimLA;ly++ )
                        {
                            table[ ly ] = ( int* ) xRealloc( table[ ly ], dimLB + 1, int );
                            table[ ly ][ dimLB ] = 0;
                        }

                        lx = 1;

                        while ( table[ 0 ][ lx ] < INTEGER( matrix ) [ x + nGt ] )
                        {
                            lx++;

                            if ( lx >= dimLB )
                                break;
                        }

                        xPos = lx;

                        if ( lx != dimLB )
                        {
                            for ( ly = 0;ly <= dimLA;ly++ )
                            {
                                for ( lx = dimLB;lx > xPos;lx-- )
                                {
                                    table[ ly ][ lx ] = table[ ly ][ lx - 1 ];
                                }
                            }

                            for ( ly = 0;ly <= dimLA;ly++ )
                            {
                                table[ ly ][ xPos ] = 0;
                            }
                        }

                        table[ 0 ][ xPos ] = INTEGER( matrix ) [ x + nGt ];
                    }

                    for ( ly = 1;ly <= dimLA;ly++ )
                    {
                        if ( table[ ly ][ 0 ] == INTEGER( matrix ) [ x ] )
                        {
                            yPos = ly;
                        }
                    }

                    if ( 0 == yPos )
                    {
                        dimLA++;
                        table = ( int** ) xRealloc( table, dimLA + 1, int* );
                        table[ dimLA ] = ( int* ) xCalloc( dimLB + 1, int );

                        table[ dimLA ][ 0 ] = INTEGER( matrix ) [ x ];

                        for ( lx = 1;lx <= dimLB;lx++ )
                        {
                            table[ dimLA ][ lx ] = 0;
                        }

                        yPos = dimLA;
                    }

                    if ( ( xPos != 0 ) && ( yPos != 0 ) )
                    {
                        table[ 0 ][ 0 ] += INTEGER( matrix ) [ x + 2 * nGt ];
                        table[ yPos ][ xPos ] = INTEGER( matrix ) [ x + 2 * nGt ];
                    }
                }
            }

            for ( ly = 0;ly <= dimLA;ly++ )
            {
                xFree( table[ ly ] );
            }

            xFree( table );
        }

        xFree( PopName );
        xFree( Loci );
    }

    return ( matrix );
}


void haplotype_init( PopName, _doubleHetero, _missing, _bg )
char**PopName;
int*_missing;
int*_doubleHetero;
int*_bg;

{
    PopListEl = NULL;
    find_list_element( *PopName, &PopListEl );

    if ( !PopListEl ) {
      info( -2, "Can't find population %s\n", *PopName );
      return ;
    }

    if ( 1 > PopListEl->Pop->NoInds ){
      info( -2, "Population %s is empty\n", *PopName );
      return ;
    }

    bg = *_bg;

    missing = *_missing;

    doubleHetero = *_doubleHetero;
}



void haplotype_next( _locusA, _locusB, _noTI, _noI, _nodH, _dLocA, _dLocB )
int*_locusA;
int*_locusB;
int*_noTI;
int*_noI;
int*_nodH;
int*_dLocA;
int*_dLocB;

{
    int nGt = 0;
    int*pGt = &nGt;
    TInd*Ind;
    POPULATION ind;
    int nLoci = 2, hap[ 2 ], dH = 0;

    dimLA = 0, dimLB = 0;
    *_noI = 0;
    locusA = *_locusA;
    locusA *= 2;
    locusB = *_locusB;
    locusB *= 2;

    startNode = NULL;
    *_noTI = PopListEl->Pop->NoInds;

    for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
    {
        Ind = get_individual( PopListEl->Pop, ind );

        if ( NULL == Ind->genotype ){
            info( -2, "first genotype_population\n" );
	    return ;
	}

        if ( doubleHetero &&
                ( Ind->genotype[ locusA ] != Ind->genotype[ locusA + 1 ] ) &&
                ( Ind->genotype[ locusB ] != Ind->genotype[ locusB + 1 ] ) )
        {
            dH++;
        }

        else
        {
            if ( missing ||
                    ( ( Ind->genotype[ locusA ] != bg ) &&
                      ( Ind->genotype[ locusA + 1 ] != bg ) &&
                      ( Ind->genotype[ locusB ] != bg ) &&
                      ( Ind->genotype[ locusB + 1 ] != bg ) ) )
            {
                hap[ 0 ] = Ind->genotype[ locusA ];
                hap[ 1 ] = Ind->genotype[ locusB ];
                startNode = insertIntLinkList( startNode, hap, 2, pGt );
                hap[ 0 ] = Ind->genotype[ locusA + 1 ];
                hap[ 1 ] = Ind->genotype[ locusB + 1 ];
                ( *_noI ) ++;
                nrGeno += 4;

                if ( 4 == nrGeno )
		  genotype = ( int* ) xCalloc( nrGeno, int );
                else
		  genotype = ( int* ) xRealloc( genotype, nrGeno, int );

                genotype[ nrGeno - 4 ] = Ind->genotype[ locusA ];
                genotype[ nrGeno - 3 ] = Ind->genotype[ locusB ];
                genotype[ nrGeno - 2 ] = Ind->genotype[ locusA + 1 ];
                genotype[ nrGeno - 1 ] = Ind->genotype[ locusB + 1 ];

                startNode = insertIntLinkList( startNode, hap, 2, pGt );
            }
        }
    }


    {
      IntNode*Node = startNode;

      if ( ( NULL != Node ) && ( NULL != ( Node = go_down( Node ) ) ) )
        {
	  int ttt = 0;
	  int x, xPos = 0, yPos = 0;
	  int*path = ( int* ) xCalloc( nGt * 3, int );
	  int lx, ly;
	  
	  for ( x = 0;x < nGt*3;x++ )
            {
	      path[ x ] = 0;
            }
	  
	  do
	    {
	      print_path( Node, &nLoci, &nGt, &path[ ttt++ ] );
	    }
	  while ( NULL != ( Node = find_next_leaf( Node ) ) );

	  for ( x = 0;x < nGt;x++ )
            {
	      if ( ( 0 == dimLA ) && ( 0 == dimLB ) )
                {
		  dimLA = 1;
		  dimLB = 1;
		  table = ( int** ) xCalloc( dimLA + 1, int );
		  table[ 0 ] = ( int* ) xCalloc( dimLB + 1, int );
		  table[ 1 ] = ( int* ) xCalloc( dimLB + 1, int );
		  table[ 0 ][ 0 ] = path[ x + 2 * nGt ];
		  table[ 1 ][ 0 ] = path[ x ];
		  table[ 0 ][ 1 ] = path[ x + nGt ];
		  table[ 1 ][ 1 ] = path[ x + 2 * nGt ];
                }
	      else
                {
		  for ( lx = 1, xPos = 0;lx <= dimLB;lx++ )
                    {
		      if ( table[ 0 ][ lx ] == path[ x + nGt ] )
                        {
			  xPos = lx;
                        }
                    }
		  if ( 0 == xPos )
                    {
		      dimLB++;
		      
		      for ( ly = 0;ly <= dimLA;ly++ )
                        {
			  table[ ly ] = ( int* ) xRealloc( table[ ly ], dimLB + 1, int );
			  table[ ly ][ dimLB ] = 0;
                        }
		      
		      lx = 1;
		      
		      while ( table[ 0 ][ lx ] < path[ x + nGt ] )
                        {
			  lx++;
			  
			  if ( lx >= dimLB )
			    break;
                        }
		      
		      xPos = lx;
		      
		      if ( lx != dimLB )
                        {
			  for ( ly = 0;ly <= dimLA;ly++ )
                            {
			      for ( lx = dimLB;lx > xPos;lx-- )
                                {
				  table[ ly ][ lx ] = table[ ly ][ lx - 1 ];
                                }
                            }
			  
			  for ( ly = 0;ly <= dimLA;ly++ )
                            {
			      table[ ly ][ xPos ] = 0;
                            }
                        }
		      
		      table[ 0 ][ xPos ] = path[ x + nGt ];
                    }
		  
		  for ( ly = 1, yPos = 0;ly <= dimLA;ly++ )
                    {
		      if ( table[ ly ][ 0 ] == path[ x ] )
                        {
			  yPos = ly;
                        }
                    }
		  
		  if ( 0 == yPos )
                    {
		      dimLA++;
		      table = ( int** ) xRealloc( table, dimLA + 1, int* );
		      table[ dimLA ] = ( int* ) xCalloc( dimLB + 1, int );
		      
		      table[ dimLA ][ 0 ] = path[ x ];
		      
		      for ( lx = 1;lx <= dimLB;lx++ )
                        {
			  table[ dimLA ][ lx ] = 0;
                        }
		      
		      yPos = dimLA;
                    }
		  
		  if ( ( xPos != 0 ) && ( yPos != 0 ) )
                    {
		      table[ 0 ][ 0 ] += path[ x + 2 * nGt ];
		      table[ yPos ][ xPos ] = path[ x + 2 * nGt ];
                    }
                }
            }
	  
	  xFree( path );
        }
    }
    
    *_nodH = dH;
    *_dLocA = dimLA;
    *_dLocB = dimLB;

    {
        int maxAllelA = 0;
        int maxAllelB = 0;

        int*transAllelA;
        int*transAllelB;

        int i;

        for ( i = 1;i <= dimLA;i++ )
        {
            if ( maxAllelA < table[ i ][ 0 ] )
            {
                maxAllelA = table[ i ][ 0 ];
            }
        }

        for ( i = 1;i <= dimLB;i++ )
        {
            if ( maxAllelB < table[ 0 ][ i ] )
            {
                maxAllelB = table[ 0 ][ i ];
            }
        }

        transAllelA = ( int* ) xCalloc( maxAllelA + 1, int );
        transAllelB = ( int* ) xCalloc( maxAllelB + 1, int );

        for ( i = 0;i <= maxAllelA;i++ )
        {
            transAllelA[ i ] = 0;
        }

        for ( i = 0;i <= maxAllelB;i++ )
        {
            transAllelB[ i ] = 0;
        }

        for ( i = 1;i <= dimLA;i++ )
        {
            transAllelA[ table[ i ][ 0 ] ] = i - 1;
        }

        for ( i = 1;i <= dimLB;i++ )
        {
            transAllelB[ table[ 0 ][ i ] ] = ( i - 1 ) * dimLA;
        }


        for ( i = 0;i < nrGeno;i += 2 )
        {
            genotype[ i ] = transAllelA[ genotype[ i ] ];
            genotype[ i + 1 ] = transAllelB[ genotype[ i + 1 ] ];
        }

        xFree( transAllelA );
        xFree( transAllelB );




    }
}



SEXP haplotype_return_table()

{
  SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;

  PROTECT( matrix = allocMatrix( INTSXP, dimLA, dimLB ) );
  PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
  INTEGER( matrixDim ) [ 0 ] = dimLA;
  INTEGER( matrixDim ) [ 1 ] = dimLB;
  setAttrib( matrix, R_DimSymbol, matrixDim );
  PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
  PROTECT( xAxis = allocVector( STRSXP, dimLB ) );
  PROTECT( yAxis = allocVector( STRSXP, dimLA ) );
  
  {
    int tmpcounter, x, y;
    char tmp[ 200 ];
    
    for ( x = 1, tmpcounter = 0;x <= dimLB;x++ )
      {
	sprintf( tmp, "locusB.%d ", table[ 0 ][ x ] );
	SET_STRING_ELT( xAxis, x - 1, mkChar( tmp ) );
	
	for ( y = 1;y <= dimLA;y++, tmpcounter++ )
	  {
	    if ( x == 1 )
	      {
		sprintf( tmp, "locusA.%d ", table[ y ][ 0 ] );
		SET_STRING_ELT( yAxis, y - 1, mkChar( tmp ) );
	      }
	    
	    INTEGER( matrix ) [ tmpcounter ] = table[ y ][ x ];
	  }
      }
  }
  
  SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
  SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
  setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
  
  UNPROTECT( 5 );
  return ( matrix );
}


SEXP haplotype_return_genotype()
  
{
  SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;
  
  PROTECT( matrix = allocMatrix( INTSXP, nrGeno / 2, 2 ) );
  PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
  INTEGER( matrixDim ) [ 0 ] = nrGeno / 2;
  INTEGER( matrixDim ) [ 1 ] = 2;
  setAttrib( matrix, R_DimSymbol, matrixDim );
  PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
  PROTECT( xAxis = allocVector( STRSXP, 2 ) );
  PROTECT( yAxis = allocVector( STRSXP, nrGeno / 2 ) );
  
  {
    int tmpcounter, x, y;
    char tmp[ 200 ];
    
    SET_STRING_ELT( xAxis, 0, mkChar( "LocusA" ) );
    SET_STRING_ELT( xAxis, 1, mkChar( "LocusB" ) );
    
    for ( x = 1, tmpcounter = 0;x <= ( nrGeno / 4 );x++ )
      {
	for ( y = 1;y <= 2;y++, tmpcounter++ )
	  {
	    sprintf( tmp, "Ind%d.%d", x, y );
	    SET_STRING_ELT( yAxis, 2 * ( x - 1 ) + ( y - 1 ), mkChar( tmp ) );
	    INTEGER( matrix ) [ tmpcounter ] = genotype[ 4 * ( x - 1 ) + 2 * ( y - 1 ) ];
	    INTEGER( matrix ) [ tmpcounter + nrGeno / 2 ] = genotype[ 4 * ( x - 1 ) + 2 * ( y - 1 ) + 1 ];
	  }
      }
  }
  
  SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
  SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
  setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
  
  UNPROTECT( 5 );
  return ( matrix );
}


void haplotype_calc_d ( double *_dd,
		        double *_rr,
		        double *_dprime )

{
  int x, y, xx, yy, sumxx, sumyy, sumTotal;
  double D, DD = 0, RA = 1, RB = 1, dPrime = 0;
  double p, q;

  RA = 1;
  RB = 1;

  for ( x = 1, sumTotal = 0; x <= dimLB; x++ ) {
    for ( y = 1; y <= dimLA; y++ ) {
      sumTotal += table[ y ][ x ];
    }
  }

  for ( x = 1;x <= dimLB;x++ ) {

    for ( yy = 1, sumyy = 0;yy <= dimLA;yy++ ) {
      sumyy += table[ yy ][ x ];
    }
      
    for ( y = 1; y <= dimLA; y++ ) {

      for ( xx = 1, sumxx = 0;xx <= dimLB;xx++ ) {
	sumxx += table[ y ][ xx ];
      }

      D = ( ( (double) table[y][x] - ( (double) sumxx * (double) sumyy ) /
	      ( (double) sumTotal ) ) / ( (double) sumTotal ) );

      DD += D * D;
      
      if ( 1 == x ) {
	RA -= (double) ( ( (double) sumxx / ( (double) sumTotal ) ) *
			 ( (double) sumxx / ( (double) sumTotal ) ) );
      }
      
      if ( 1 == y ) {
	RB -= (double) ( ( (double) sumyy / ( (double) sumTotal ) ) *
			 ( (double) sumyy / ( (double) sumTotal ) ) );
      }
      
      p = (double) sumxx / ( (double) sumTotal );
      q = (double) sumyy / ( ( double ) sumTotal );

      if ( D < 0 ) {

	  if ( (p*q) < ( (1-p)*(1-q) ) ) {
	    dPrime += ( p * q * D * ( -1 ) / ( p * q ) );
	  }
	  else {
	    dPrime += ( p * q * D * ( -1 ) / ( ( 1 - p ) * ( 1 - q ) ) );
	  }

      }
      else {

          if ( ( (1-p)*q ) < (p*(1-q) ) ) {
	    dPrime += ( p * q * D / ( ( 1 - p ) * q ) );
	  }
	  else {
	    dPrime += ( p * q * D / ( p * ( 1 - q ) ) );
	  }

      }
    }
  }

  *_dd = DD;
  *_rr = DD / ( RA * RB );
  *_dprime = dPrime;
}



void haplotype_free()

{
    int ly;

    for ( ly = 0;ly <= dimLA;ly++ )
    {
        xFree( table[ ly ] );
    }

    xFree( table );
    xFree( genotype );
    nrGeno = 0;
}



void concat_population( NameP1, NameP2 )
char**NameP1;
char**NameP2;
{
    TInd tmpInd;
    POPULATION ind = 0, j, oldP1Alloc;
    TPopListEl*PopListEl;
    TPop*P1Pop = NULL;
    TPop*P2Pop = NULL;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl ){
      info( -2, "Can't find population %s!", *NameP2 );
      return ;
    }

    remove_evaluate_population( NameP2 );

    remove_genotype_population( NameP2 );

    P2Pop = PopListEl->Pop;

    find_list_element( *NameP1, &PopListEl );

    if ( NULL == PopListEl )
    {
        P1Pop = new_population( NameP1, &ind );
	if ( NULL == P1Pop ) return ;
    }

    else
    {
        P1Pop = PopListEl->Pop;
        remove_evaluate_population( NameP1 );
        remove_genotype_population( NameP1 );
    }

    info( 1, "Population %s concatenated to population %s\n", *NameP2, *NameP1 );

    ind = P1Pop->NoInds, j = 0, oldP1Alloc = P1Pop->NoIndsAlloc;
    P1Pop->NoInds += P2Pop->NoInds;
    P1Pop->NoIndsAlloc += P2Pop->NoIndsAlloc;

    if ( NULL == P1Pop->Ind )
    {
        P1Pop->Ind = ( TInd* ) xCalloc( P1Pop->NoIndsAlloc, TInd );
    }

    else
    {
        P1Pop->Ind = ( TInd* ) xRealloc( P1Pop->Ind, P1Pop->NoIndsAlloc, TInd );
    }

    for ( ;j < P2Pop->NoIndsAlloc;ind++, j++ )
    {
        tmpInd = P1Pop->Ind[ ind ];
        P1Pop->Ind[ ind ] = P2Pop->Ind[ j ];

        if ( ind < oldP1Alloc )
            P1Pop->Ind[ P1Pop->NoIndsAlloc - j - 1 ] = tmpInd;
    }

    P2Pop->NoIndsAlloc = 0;
    remove_population( NameP2 );
}


void sample_population( NameP1, NameP2, cNoI, rep )
char**NameP1;
char**NameP2;
char**cNoI;
int*rep;

{
    int chrom, hom, loc, r;
    POPULATION ind, NoI = atol( *cNoI ), *samples = NULL;
    TInd*P1Ind, *P2Ind;
    TPopListEl*PopListEl;

    TPop*P1Pop = exist_population( *NameP1, NoI, 1 ), *P2Pop = NULL;
    if ( NULL == P1Pop )  return;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl )
      {
        info( -2, "Can't find population %s \n", *NameP2 );
	return ;
      }

    P2Pop = PopListEl->Pop;

    if ( 0 == P2Pop->NoInds ){
        info( -2, "Population >%s< is empty\n", *NameP2 );
	return;
      }

    if ( -1 == NoI )
        P1Pop->NoInds = P1Pop->NoIndsAlloc;
    else
        P1Pop->NoInds = NoI;

    if ( ( 0 == *rep ) && ( P2Pop->NoInds < P1Pop->NoInds ) )
      {
        info( -2, "sampling only with repetition possible\n" );
	return;        
      }

    for ( ind = 0;ind < P1Pop->NoInds;ind++ )
    {

      r = (int) gsl_rng_uniform_int( rng[0],  (unsigned long int)P2Pop->NoInds );


        if ( 0 == *rep )
        {
            int j = 0;
            POPULATION ind_i;

            do
            {
                j = 0;

                for ( ind_i = 0;ind_i < ind;ind_i++ )
                {
                    if ( r == *( samples + ind_i ) )
                    {
                        j = 1;
                        break;
                    }
                }

                if ( 0 == j )
                {
                    if ( 0 == ind )
                        samples = ( POPULATION* ) xCalloc( 1, POPULATION );
                    else
                        samples = ( POPULATION* ) xRealloc( samples, ind + 1, POPULATION );

                    *( samples + ind ) = r;
                }

                else
		  r = (int) gsl_rng_uniform_int( rng[0], 
						 (unsigned long int)P2Pop->NoInds );
            }

            while ( j == 1 );
        }

        P1Ind = get_individual( P1Pop, ind );
        P2Ind = get_individual( P2Pop, r );

        if ( NULL != P2Ind->genotype )
        {
            int i;
            P1Ind->genotype = ( int* ) xCalloc( cMPhom, int );

            for ( i = 0;i < cMPhom;i++ )
            {
                P1Ind->genotype[ i ] = P2Ind->genotype[ i ];
            }
        }

        else
        {
            P1Ind->genotype = NULL;
        }

        if ( NULL != P2Ind->GValue )
        {
            int i;
            P1Ind->GValue = ( float* ) xCalloc( NrEffects + 1, float );

            for ( i = 0;i < NrEffects + 1;i++ )
            {
                P1Ind->GValue[ i ] = P2Ind->GValue[ i ];
            }
        }

        else
        {
            P1Ind->GValue = NULL;
        }

        if ( NULL != P2Ind->info )
        {
            P1Ind->info = ( char* ) xCalloc( MAXNAME, char );
            strcpy( P1Ind->info, P2Ind->info );
        }

        else
        {
            P1Ind->info = NULL;
        }

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( P1Ind->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( P2Ind->Chrom[ chrom ] ).Hom[ hom ];
                P1Hom->NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    insertLocus( P2Hom->Loc[ loc ].All, P2Hom->Loc[ loc ].Pos, P1Hom );
                }
            }
        }

    }
}




void copy_population( NameP1, NameP2, cstart, cn )
char**NameP1;
char**NameP2;
char**cstart;
char**cn;

{
    int chrom, hom, loc;
    POPULATION ind, start = atol( *cstart ), n = atol( *cn );
    TInd*P1Ind, *P2Ind;
    TPop*P1Pop = NULL, *P2Pop = NULL;
    TPopListEl*PopListEl;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl )
      {
        info( -2, "Can't find population %s !\n", *NameP2 );
	return;
      }

    P2Pop = PopListEl->Pop;

    if ( 0 == P2Pop->NoInds )
      {
        info( -2, "Population >%s< is empty\n", *NameP2 );
	return;
      }

    start--;

    if ( start < 0 )
        start = 0;

    if ( n < 0 )
        n = P2Pop->NoInds;

    if ( start + n > P2Pop->NoInds )
        n = P2Pop->NoInds - start;

    P1Pop = exist_population( *NameP1, n, 1 );
    if (NULL == P1Pop) return;

    P1Pop->NoInds = n;

    for ( ind = 0;ind < n;ind++ )
    {
        P1Ind = get_individual( P1Pop, ind );
        P2Ind = get_individual( P2Pop, start + ind );


        if ( NULL != P2Ind->genotype )
        {
            int i;
            P1Ind->genotype = ( int* ) xCalloc( cMPhom, int );

            for ( i = 0;i < cMPhom;i++ )
            {
                P1Ind->genotype[ i ] = P2Ind->genotype[ i ];
            }
        }

        else
        {
            P1Ind->genotype = NULL;
        }

        if ( NULL != P2Ind->GValue )
        {
            int i;
            P1Ind->GValue = ( float* ) xCalloc( NrEffects + 1, float );

            for ( i = 0;i < NrEffects + 1;i++ )
            {
                P1Ind->GValue[ i ] = P2Ind->GValue[ i ];
            }
        }

        else
        {
            P1Ind->GValue = NULL;
        }

        if ( NULL != P2Ind->info )
        {
            P1Ind->info = ( char* ) xCalloc( MAXNAME, char );
            strcpy( P1Ind->info, P2Ind->info );
        }

        else
        {
            P1Ind->info = NULL;
        }

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( P1Ind->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( P2Ind->Chrom[ chrom ] ).Hom[ hom ];
                P1Hom->NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    insertLocus( P2Hom->Loc[ loc ].All, P2Hom->Loc[ loc ].Pos, P1Hom );
                }
            }
        }

    }

    info( 1, "Copied %d individuals (Index %d to %d) from population %s to population %s.\n",
          n, start + 1, start + n, *NameP2, *NameP1 );
}

void append_population( NameP1, NameP2 )
char**NameP1;
char**NameP2;

{
    int chrom, hom, loc;
    POPULATION P1Inds = 0, ind;
    TInd*P1Ind, *P2Ind;
    TPopListEl*PopListEl;
    TPop*P1Pop = NULL, *P2Pop = NULL;

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Can't find population %s \n", *NameP2 );
	return;
    }

    P2Pop = PopListEl->Pop;

    if ( 0 == P2Pop->NoInds ) {
        info( -2, "Population >%s< is empty\n", *NameP2 );
	return;
    }

    find_list_element( *NameP1, &PopListEl );

    if ( NULL != PopListEl )
        P1Inds = PopListEl->Pop->NoInds;

    P1Pop = exist_population( *NameP1, P2Pop->NoInds + P1Inds, 0 );
    if ( NULL == P1Pop )  return;

    P1Pop->NoInds = P1Pop->NoIndsAlloc;

    P1Ind = get_individual( P1Pop, 0 );

    P2Ind = get_individual( P2Pop, 0 );

    if ( 0 < P1Inds )
    {
        if ( ( NULL != P1Ind->GValue ) && ( NULL == P2Ind->GValue ) )
            remove_evaluate_population( NameP1 );

        if ( ( NULL == P1Ind->GValue ) && ( NULL != P2Ind->GValue ) )
            remove_evaluate_population( NameP2 );

        if ( ( NULL != P1Ind->genotype ) && ( NULL == P2Ind->genotype ) )
            remove_genotype_population( NameP1 );

        if ( ( NULL == P1Ind->genotype ) && ( NULL != P2Ind->genotype ) )
            remove_genotype_population( NameP2 );
    }

    for ( ind = 0;ind < P2Pop->NoInds;ind++ )
    {
        P1Ind = get_individual( P1Pop, ind + P1Inds );
        P2Ind = get_individual( P2Pop, ind );


        if ( NULL != P2Ind->genotype )
        {
            int i;
            P1Ind->genotype = ( int* ) xCalloc( cMPhom, int );

            for ( i = 0;i < cMPhom;i++ )
            {
                P1Ind->genotype[ i ] = P2Ind->genotype[ i ];
            }
        }

        else
        {
            P1Ind->genotype = NULL;
        }

        if ( NULL != P2Ind->GValue )
        {
            int i;
            P1Ind->GValue = ( float* ) xCalloc( NrEffects + 1, float );

            for ( i = 0;i < NrEffects + 1;i++ )
            {
                P1Ind->GValue[ i ] = P2Ind->GValue[ i ];
            }
        }

        else
        {
            P1Ind->GValue = NULL;
        }

        if ( NULL != P2Ind->info )
        {
            P1Ind->info = ( char* ) xCalloc( MAXNAME, char );
            strcpy( P1Ind->info, P2Ind->info );
        }

        else
        {
            P1Ind->info = NULL;
        }

        for ( chrom = 0;chrom < NoChroms;chrom++ )
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                THom*P1Hom = &( P1Ind->Chrom[ chrom ] ).Hom[ hom ];
                THom*P2Hom = &( P2Ind->Chrom[ chrom ] ).Hom[ hom ];
                P1Hom->NoLociUsed = 0;

                for ( loc = 0;loc < P2Hom->NoLociUsed;loc++ )
                {
                    insertLocus( P2Hom->Loc[ loc ].All, P2Hom->Loc[ loc ].Pos, P1Hom );
                }
            }
        }

    }

    info( 1, "Population %s (%d individuals) appended to population %s.\n",
          *NameP2, P2Pop->NoInds, *NameP1 );
}


void select_all_n_best ( char   **NamePop,
		         char   **indexbest,
		         int    *n,
		         double *minscore,
		         double *maxscore  )

{
    TPopListEl *PopListEl;
    TPop       *Pop = NULL;
    TInd       *Ind0, *Ind1;
    POPULATION ind = 0;
    int        _n = *n;

    find_list_element( *NamePop, &PopListEl );

    if ( PopListEl != NULL )
        Pop = PopListEl->Pop;
    else
        info( -2, "Can't find population %s\n", *NamePop );

    Ind0 = get_individual( Pop, ind );

    if ( NULL == Ind0->GValue ){
        info( -2, "No Effects calculated\n" );
	return;
    }

    *maxscore = Ind0->GValue[0];

    for ( ind = 1; ind < Pop->NoInds; ind++ )
      {
        Ind1 = get_individual( Pop, ind );
	
        if ( 0 != ( Ind0->GValue[ 0 ] - Ind1->GValue[ 0 ] ) )
	  if ( 0 == --_n ) break;
	
        Ind0 = Ind1;
      }

    *minscore = Ind0->GValue[0];
    sprintf( *indexbest, "%ld", ind );
    info( 1, 
	  "%d individuals selected, (score %f .. %f)\n", 
	  ind, 
	  *minscore, 
	  *maxscore );
}

void divide_population ( char **NameP1,
			 char **NameP2,
			 char **cNoI   )
{
  POPULATION NoI = atol( *cNoI );
  TPopListEl *PopListEl;
  TPop *P1Pop, *P2Pop = NULL;

  find_list_element( *NameP2, &PopListEl );

  if ( PopListEl != NULL ) 
    {
      P2Pop = PopListEl->Pop;
      P1Pop = exist_population( *NameP1, NoI, 1 );
      if ( NULL == P1Pop )  return;

      if ( NoI > P2Pop->NoInds )
	{
	  info( -1, 
		"population.divide - population size adjusted %d -> %d \n", 
		NoI, 
		P2Pop->NoInds );

	  NoI = P2Pop->NoInds;
	}
      
      {
        TInd tmpInd;
        POPULATION ind, ind_j;

        P1Pop->NoInds = NoI;

        for ( ind = 0; 
	      ind < P1Pop->NoInds; 
	      ind++ )
        {
            tmpInd = P1Pop->Ind[ ind ];
            P1Pop->Ind[ ind ] = P2Pop->Ind[ ind ];
            P2Pop->Ind[ ind ] = tmpInd;
        }

        for ( ind = NoI, ind_j = 0;
	      ind < P2Pop->NoInds;
	      ind++, ind_j++ )
        {
            tmpInd = P2Pop->Ind[ ind ];
            P2Pop->Ind[ ind ] = P2Pop->Ind[ ind_j ];
            P2Pop->Ind[ ind_j ] = tmpInd;
        }

        P2Pop->NoInds -= NoI;
      }
    }   else
    {
      info( -2, "Population %s is not defined.", *NameP2 );
      return;
    }

}

void rename_population( OldName, NewName )
char**OldName;
char**NewName;

{
    TPopListEl*PopListEl1, *PopListEl2;
    TPop*Pop = NULL;

    find_list_element( *OldName, &PopListEl1 );

    if ( PopListEl1 == NULL ){
        info( -2, "Can't find population %s\n", *OldName );
	return;
    }

    find_list_element( *NewName, &PopListEl2 );

    if ( PopListEl2 != NULL ){
        info( -2, "Populationname >%s< is already used\n", *NewName );
	return;
    }

    if ( strlen( *NewName ) >= MAXNAME ){
        info( -2, "New populationname is too long\n" );
	return;
    }

    Pop = PopListEl1->Pop;

    strcpy( Pop->Name, *NewName );
}


void swap_population_name( NameP1, NameP2 )
char**NameP1;
char**NameP2;

{
    TPopListEl*PopListEl;
    find_list_element( *NameP1, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Can't find population %s\n", *NameP1 );
	return;
    }

    find_list_element( *NameP2, &PopListEl );

    if ( NULL == PopListEl ){
        info( -2, "Can't find population  %s\n", *NameP2 );
	return;
    }

    if ( 0 == strcmp( *NameP1, *NameP2 ) ){
        info( -2, "You have to specify two different populations\n" );
	return;
    }

    {

        {
            TPopListEl*PopListEl;
            char tmpName[ MAXNAME ];
            TPop*P1Pop = NULL, *P2Pop = NULL;

            find_list_element( *NameP1, &PopListEl );

            if ( PopListEl != NULL )
                P1Pop = PopListEl->Pop;

            if ( NULL == P1Pop ){
                info( -2, "Can't find population %s\n", *NameP1 );
		return;
	    }

            find_list_element( *NameP2, &PopListEl );

            if ( PopListEl != NULL )
                P2Pop = PopListEl->Pop;

            if ( NULL == P2Pop ){
                info( -2, "Can't find population  %s\n", *NameP2 );
		return;
	    }

            strcpy( tmpName, P1Pop->Name );

            strcpy( P1Pop->Name, P2Pop->Name );

            strcpy( P2Pop->Name, tmpName );
        }

    }
}

void get_population_size( PopNames, cNoInds, cNoIndsAlloc )
char**PopNames;
char**cNoInds;
char**cNoIndsAlloc;

{
    TPopListEl*PopListEl;
    POPULATION NoInds = 0, NoIndsAlloc = 0;
    char*PopName;

    PopName = strtok( *PopNames, " " );

    while ( NULL != PopName )
    {
        find_list_element( PopName, &PopListEl );

        if ( NULL == PopListEl ){
	  info( -2, "Can't find population %s\n", PopName );
	  return;
	}

        NoInds += PopListEl->Pop->NoInds;

        NoIndsAlloc += PopListEl->Pop->NoIndsAlloc;

        PopName = strtok( NULL, " " );
    }

    sprintf( *cNoInds, "%ld", NoInds );
    sprintf( *cNoIndsAlloc, "%ld", NoIndsAlloc );
}


void get_genome_par( int*_NoChroms, int*_NoHoms )

{
    *_NoChroms = NoChroms;
    *_NoHoms = NoHoms;
}

void get_genome_par_b( double*_ChromLen )

{
    int chrom;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        _ChromLen[ chrom ] = PopDescr[ chrom ].ChromLen;
    }

}


void get_map_a(
    int*_NoMappoints,
    int*_maxname,
    char**_classlgth,
    char**_namelgth
)
{
    int chrom;
    int point;
    unsigned long int class_lgth = 0;
    unsigned long int name_lgth = 0;

    *_NoMappoints = 0;
    *_maxname = MAXNAME;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        *_NoMappoints += PopDescr[ chrom ].NrMappoints;

        for ( point = 0;point < PopDescr[ chrom ].NrMappoints;point++ )
        {
            name_lgth += ( 1 + strlen( PopDescr[ chrom ].Mappoint[ point ].name ) );
            class_lgth += ( 1 + strlen( PopDescr[ chrom ].Mappoint[ point ].class ) );
        }
    }

    sprintf( *_classlgth, "%lu", class_lgth );
    sprintf( *_namelgth, "%lu", name_lgth );
}

void get_map_b(
    double*_chrom,
    double*_pos,
    char**_name,
    char**_class
)
{
    int chrom;
    int point;
    long counter = 0;
    char*p = *_name;
    char*q = *_class;
    int len;



    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        for ( point = 0;point < PopDescr[ chrom ].NrMappoints;point++ )
        {
            _chrom[ counter ] = chrom + 1;
            _pos[ counter ] = PopDescr[ chrom ].Mappoint[ point ].Pos;
            len = strlen( PopDescr[ chrom ].Mappoint[ point ].name );

            if ( len > MAXNAME )
                len = MAXNAME;

            strncpy( p, PopDescr[ chrom ].Mappoint[ point ].name, len );

            p += ( len + 1 );

            len = strlen( PopDescr[ chrom ].Mappoint[ point ].class );

            if ( len > MAXNAME )
                len = MAXNAME;

            strncpy( q, PopDescr[ chrom ].Mappoint[ point ].class, len );

            q += ( len + 1 );

            counter++;
        }
    }
}



TMappoint*find_mappoint( char*pch )
{
    int chrom, i;

    for ( chrom = 0;chrom < NoChroms;chrom++ )
    {
        if ( PopDescr[ chrom ].NrMappoints > 0 )
        {
            for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
            {
                if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ i ].name, pch ) )
                    return & ( PopDescr[ chrom ].Mappoint[ i ] );
            }
        }
    }

    return NULL;
}


void insertmember( Set, newmember )
TSet*Set;
int newmember;
{
    int i;

    for ( i = 0;i < Set->betrag;i++ )
    {
        if ( Set->member[ i ] == newmember )
            return ;
    }

    if ( 0 == Set->betrag++ )
    {
        Set->member = ( int* ) xCalloc( Set->betrag, int );
    }

    else
    {
        Set->member = ( int* ) xRealloc( Set->member, Set->betrag, int );
    }

    Set->member[ Set->betrag - 1 ] = newmember;
}


TEffects* find_Effect( char*Name )
{
    int i;

    for ( i = 0; i < NrEffects; i++ )
    {
        if ( 0 == strcmp( Name, Effects[ i ].Name ) )
            return & ( Effects[ i ] );
    }

    return NULL;
}


void list_eff_parameter( int* lEffName, int* nrEffects )
{
    int i;
    *nrEffects = NrEffects;

    for ( i = 0; i < NrEffects; i++ )
    {
        *lEffName += strlen( Effects[ i ].Name ) + 1;
    }

    ( *lEffName ) += 10;
}



void list_eff( char**all_eff, double*weigth_eff )
{
    int i;
    strcpy( *all_eff, "" );

    for ( i = 0;i < NrEffects;i++ )
    {
        strcat( *all_eff, Effects[ i ].Name );
        strcat( *all_eff, " " );
        weigth_eff[ i ] = Effects[ i ].weight;
    }
}



void get_population_gvalue( char   **PopName,
			    double *gvalue,
			    char   **EffName )
{
  int         eff = 0;
  POPULATION  ind;
  TPopListEl  *PopListEl;
  TEffects    *Effect;
  TInd        *Ind;

  find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl )
      {
        info( -2, "Population %s not defined\n", *PopName );
	return;
      }

    if ( NULL != EffName )
      {
        if ( NULL == ( Effect = find_Effect( *EffName ) ) )
	  {
            info( -2, "Effect %s not defined\n", *EffName );
	    return;
	  }
	
        eff = Effect->Pos;
      }

    for ( ind = 0; 
	  ind < PopListEl->Pop->NoInds;
	  ind++ )
    {
        Ind = get_individual( PopListEl->Pop, ind );

        if ( NULL == Ind->GValue )
	  {
            info( -2, "Population not evaluated\n" );
	    return;
	  }

        *( gvalue + ind ) = (double) Ind->GValue[eff];
    }
}


void set_population_gvalue ( char   **PopName,
			     double *gvalue )

{
  POPULATION  ind;
  TPopListEl *PopListEl;
  TInd       *Ind;

  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl )
    {
      info( -2, "Population %s not defined\n", *PopName );
      return;
    }
  
  for ( ind = 0; 
	ind < PopListEl->Pop->NoInds;
	ind++ )
  {
      Ind = ( TInd* ) get_individual( PopListEl->Pop, ind );
      
      if ( NULL == Ind->GValue )
	Ind->GValue = ( float* ) xCalloc( 1, float );
      
      Ind->GValue[0] = ( float ) * ( gvalue + ind );
  }
}


void set_population_info ( char **PopName,
			   char **_ind,
			   char **_info )
{
  POPULATION  ind = atol( *_ind ) - 1;
  TPopListEl *PopListEl;
  TInd       *Ind;
  
  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl ) 
      {
	info( -2, "Population %s not defined\n", *PopName );
	return;
      }
  
  Ind = ( TInd* ) get_individual( PopListEl->Pop, ind );
  
  if ( NULL == Ind->info )
    Ind->info = ( char* ) xCalloc( MAXNAME, char );
  
  if ( strlen( *_info ) >= MAXNAME )
    {
      info( -2, "Info string is too long\n" );
      return;
    }
  
  strcpy( Ind->info, *_info );
}


void get_population_info ( char **PopName,
			   char **_ind,
			   char **_info )
{
    POPULATION  ind = atol( *_ind ) - 1;
    TPopListEl  *PopListEl;
    TInd        *Ind;

    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl ) 
      {
	info( -2, "Population %s not defined\n", *PopName );
	return;
      }

    Ind = ( TInd* ) get_individual( PopListEl->Pop, ind );
    
    if ( NULL == Ind->info )
      strcpy( *_info, "none" );
    else
      strcpy( *_info, Ind->info );
}



void set_eff_weight( char **Name, double *weight )
{
    TEffects*Effect = find_Effect( *Name );

    if ( NULL == Effect ) {
        info( -2, "Can't find effect %s\n", *Name );
	return;
    }

    Effect->weight = ( float ) * weight;
}

void effect_weight_set_all ( double *weight )
{
    int i;

    for ( i = 0;i < NrEffects;i++ )
    {
        Effects[ i ].weight = *weight;
    }
}


void Insert_ALL_Loc( TEff*Eff, TMappoint*MP, int All )
{
    int i;

    for ( i = 0;i < Eff->NrLoci;i++ )
    {
        if ( MP == Eff->Loc[ i ].MP )
        {
            if ( Eff->Loc[ i ].NrAll++ < NoHoms )
            {
                Eff->Loc[ i ].All[ Eff->Loc[ i ].NrAll - 1 ] = All;
            }

            else
            {
                Eff->Loc[ i ].NrAll--;
                info( -2, "Error in the effect-file - too many alleles\n" );
            }

            return ;
        }
    }

    if ( 0 == Eff->NrLoci++ )
        Eff->Loc = ( TGen* ) xCalloc( 1, TGen );
    else
        Eff->Loc = ( TGen* ) xRealloc( Eff->Loc, Eff->NrLoci, TGen );

    Eff->Loc[ Eff->NrLoci - 1 ].MP = MP;
    Eff->Loc[ Eff->NrLoci - 1 ].NrAll = 1;
    Eff->Loc[ Eff->NrLoci - 1 ].All = ( int* ) xCalloc( NoHoms, int );
    Eff->Loc[ Eff->NrLoci - 1 ].All[ 0 ] = All;

    for ( i = 1;i < NoHoms;i++ )
        Eff->Loc[ Eff->NrLoci - 1 ].All[ i ] = -1;
}


void load_effmap_mod ( char**effectName,
		       char**effectSpec,
		       double*effectWeight )

{
    char line[ 10 * MAXLINE ], *token;
    int i, j, k, tmpi, status, maxLoc = 0;
    float tmpf;
    TMappoint*tmpMP = NULL;
    TEffects*Effect;

    if ( NULL != find_Effect( *effectName ) ){
        info( -2, "This effect is already loaded (same name)\n" );
	return;
    }

    if ( 0 == NrEffects++ )
        Effects = ( TEffects* ) xCalloc( 1, TEffects );
    else
        Effects = ( TEffects* ) xRealloc( Effects, NrEffects, TEffects );

    Effect             = &Effects[ NrEffects - 1 ];
    Effect->GValue_avg = 0;
    Effect->weight     = 1;
    Effect->NrEff      = 0;
    Effect->Pos        = NrEffects;
    strcpy( Effect->Name, *effectName );
    Effect->weight = ( float ) * effectWeight;
    Effect->GValue_avg = 0;

    for ( i = 0; i < 1; i++ )
    {
        strcpy( line, *effectSpec );
        status = 0;
        token = strtok( line, "[ \t]" );

        while ( ( NULL != token ) && ( status != -1 ) && ( status != 1 ) )
        {
            switch ( status )
            {

            case 0:

                if ( 0 == strcmp( line, "\n" ) )
                    status = -1;
                else
                {
                    status = 2;

                    if ( 0 == ( Effect->NrEff ) ++ )
                        Effect->Eff = ( TEff* ) xCalloc( 1, TEff );
                    else
                        Effect->Eff = ( TEff* ) xRealloc( Effect->Eff, Effect->NrEff, TEff );

                    Effect->Eff[ ( Effect->NrEff ) - 1 ].Eff = 0;

                    Effect->Eff[ ( Effect->NrEff ) - 1 ].NrLoci = 0;
                }

                break;

            case 2:

                if ( NULL != ( tmpMP = find_mappoint( token ) ) )
                {
                    status = 3;
                }

                else
                {
                    status = -1;
                    info( -1, "Marker %s not found.\nline %d\n", token, i );
                }

                break;

            case 3:

                if ( 1 == sscanf( token, "%d", &tmpi ) )
                {
                    status = 4;
                    Insert_ALL_Loc( &Effect->Eff[ Effect->NrEff - 1 ], tmpMP, tmpi );
                }

                else
                {
                    status = -1;
                    info( -1, "Integer expected\nline %d\n", i );
                }

                break;

            case 4:

                if ( 1 == sscanf( token, "%f", &tmpf ) )
                {
                    status = 1;
                    Effect->Eff[ Effect->NrEff - 1 ].Eff = tmpf;

                    if ( maxLoc < Effect->Eff[ Effect->NrEff - 1 ].NrLoci )
                    {
                        maxLoc = Effect->Eff[ Effect->NrEff - 1 ].NrLoci;
                    }
                }

                else
                {
                    status = 2;
                }

                break;
            }

            if ( 2 != status )
                token = strtok( NULL, "[ \t]" );
        }
    }

    Effect->MPs = 0;

    for ( i = 0;i < Effect->NrEff;i++ )
    {
        for ( j = 0;j < Effect->Eff[ i ].NrLoci;j++ )
        {
            for ( k = 0;k < Effect->MPs;k++ )
            {
                if ( *( Effect->MPset + k ) == ( Effect->Eff[ i ].Loc[ j ].MP ) )
                    break;
            }

            if ( k == Effect->MPs )
            {
                Effect->MPs++;

                if ( 1 == Effect->MPs )
                    Effect->MPset = ( TMappoint** ) xCalloc( 1, TMappoint* );
                else
                    Effect->MPset = ( TMappoint** ) xRealloc( Effect->MPset, Effect->MPs, TMappoint* );

                Effect->MPset[ Effect->MPs - 1 ] = Effect->Eff[ i ].Loc[ j ].MP;
            }
        }
    }
}

void load_effmap( file, name )
char**file;
char**name;
{
    FILE*fp;
    char line[ 10 * MAXLINE ], *token;
    int i, j, k, tmpi, status, maxLoc = 0;
    float tmpf;
    TMappoint*tmpMP = NULL;
    TEffects*Effect;

    if ( NULL != find_Effect( *name ) ){
        info( -2, "This effect is already loaded (same name)\n" );
	return;
    }

    if ( 0 == NrEffects++ )
        Effects = ( TEffects* ) xCalloc( 1, TEffects );
    else
        Effects = ( TEffects* ) xRealloc( Effects, NrEffects, TEffects );

    Effect = &Effects[ NrEffects - 1 ];

    Effect->GValue_avg = 0;

    Effect->weight = 1;

    Effect->NrEff = 0;

    Effect->Pos = NrEffects;

    Effect->Eff = NULL;

    Effect->MPset = NULL;

    strcpy( Effect->Name, *name );

    fp = fileopen( *file, "r" );

    for ( i = 0;;i++ )
    {
        if ( fgets( line, 10 * MAXLINE, fp ) == NULL )
            break;

        if ( (unsigned)( j = strcspn( line, "#" ) ) != strlen( line ) )
        {
            if ( j == 0 )
                strcpy( line, "" );
            else
            {
                line[ j++ ] = '\n';
                line[ j ] = '\0';
            }
        }

        status = 0;
        token = strtok( line, "[ \t]" );

        while ( ( NULL != token ) && ( status != -1 ) && ( status != 1 ) )
        {
            switch ( status )
            {

            case 0:

                if ( 0 == strcmp( token, "Weigth:" ) )
                {
                    status = 5;
                }

                else if ( 1 == sscanf( token, "%f", &Effect->GValue_avg ) )
                {
                    status = 1;
                }

                else
                {
                    if ( 0 == strcmp( line, "\n" ) )
                        status = -1;
                    else
                    {
                        status = 2;

                        if ( 0 == ( Effect->NrEff ) ++ )
                            Effect->Eff = ( TEff* ) xCalloc( 1, TEff );
                        else
                            Effect->Eff = ( TEff* ) xRealloc( Effect->Eff, Effect->NrEff, TEff );

                        Effect->Eff[ ( Effect->NrEff ) - 1 ].Eff = 0;

                        Effect->Eff[ ( Effect->NrEff ) - 1 ].NrLoci = 0;
                    }
                }

                break;

            case 2:

                if ( NULL != ( tmpMP = find_mappoint( token ) ) )
                {
                    status = 3;
                }

                else
                {
                    status = -1;
                    info( -1, "Marker %s not found.\nline %d\n", token, i );
                }

                break;

            case 3:

                if ( 1 == sscanf( token, "%d", &tmpi ) )
                {
                    status = 4;
                    Insert_ALL_Loc( &Effect->Eff[ Effect->NrEff - 1 ], tmpMP, tmpi );
                }

                else
                {
                    status = -1;
                    info( -1, "Integer expected\nline %d\n", i );
                }

                break;

            case 4:

                if ( 1 == sscanf( token, "%f", &tmpf ) )
                {
                    status = 1;
                    Effect->Eff[ Effect->NrEff - 1 ].Eff = tmpf;

                    if ( maxLoc < Effect->Eff[ Effect->NrEff - 1 ].NrLoci )
                    {
                        maxLoc = Effect->Eff[ Effect->NrEff - 1 ].NrLoci;
                    }
                }

                else
                {
                    status = 2;
                }

                break;

            case 5:

                if ( 1 == sscanf( token, "%f", &Effect->weight ) )
                {
                    status = 1;
                }

                else
                {
                    status = -1;
                    info( -1, "Error: weight" );
                }

                break;
            }

            if ( 2 != status )
                token = strtok( NULL, "[ \t]" );
        }
    }

    fileclose( fp );

    Effect->MPs = 0;

    for ( i = 0;i < Effect->NrEff;i++ )
    {
        for ( j = 0;j < Effect->Eff[ i ].NrLoci;j++ )
        {
            for ( k = 0;k < Effect->MPs;k++ )
            {
                if ( *( Effect->MPset + k ) == ( Effect->Eff[ i ].Loc[ j ].MP ) )
                    break;
            }

            if ( k == Effect->MPs )
            {
                Effect->MPs++;

                if ( 1 == Effect->MPs )
                    Effect->MPset = ( TMappoint** ) xCalloc( 1, TMappoint* );
                else
                    Effect->MPset = ( TMappoint** ) xRealloc( Effect->MPset, Effect->MPs, TMappoint* );

                Effect->MPset[ Effect->MPs - 1 ] = Effect->Eff[ i ].Loc[ j ].MP;
            }
        }
    }


    {
        int l, m;

        for ( i = 0;i < maxLoc;i++ )
        {
            for ( k = 0;k < NoHoms;k++ )
            {
                int oldEff = Effect->NrEff;
                int n;

                for ( n = k + 1;n < NoHoms;n++ )
                {
                    for ( j = 0;j < oldEff;j++ )
                    {
                        if ( ( i < Effect->Eff[ j ].NrLoci ) && ( k < Effect->Eff[ j ].Loc[ i ].NrAll ) &&
                                ( Effect->Eff[ j ].Loc[ i ].All[ k ] != Effect->Eff[ j ].Loc[ i ].All[ n ] ) )
                        {
                            Effect->NrEff++;
                            Effect->Eff = ( TEff* ) xRealloc( Effect->Eff, Effect->NrEff, TEff );
                            Effect->Eff[ Effect->NrEff - 1 ].Eff = Effect->Eff[ j ].Eff;
                            Effect->Eff[ Effect->NrEff - 1 ].NrLoci = Effect->Eff[ j ].NrLoci;
                            Effect->Eff[ Effect->NrEff - 1 ].Loc =
                                ( TGen* ) xCalloc( Effect->Eff[ j ].NrLoci, TGen );

                            for ( l = 0;l < Effect->Eff[ j ].NrLoci;l++ )
                            {
                                Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].MP =
                                    Effect->Eff[ j ].Loc[ l ].MP;
                                Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].NrAll =
                                    Effect->Eff[ j ].Loc[ l ].NrAll;
                                Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All =
                                    ( int* ) xCalloc( NoHoms, int );

                                for ( m = 0;m < NoHoms;m++ )
                                {
                                    if ( ( l == i ) && ( m == k ) )
                                    {
                                        Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All[ m ] = Effect->Eff[ j ].Loc[ l ].All[ n ];
                                    }

                                    else if ( ( l == i ) && ( m == n ) )
                                    {
                                        Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All[ m ] = Effect->Eff[ j ].Loc[ l ].All[ k ];
                                    }

                                    else
                                    {
                                        Effect->Eff[ Effect->NrEff - 1 ].Loc[ l ].All[ m ] = Effect->Eff[ j ].Loc[ l ].All[ m ];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void permutate_effmap()
{}



void effmap_remove( effectName )
char**effectName;
{
    if ( NULL == find_Effect( *effectName ) )
    {
        info( -2, "Couldn't find effect: %s\n", *effectName );
	return;
    }

    else
    {
        int i, j, k;

        for ( i = NrEffects - 1;i >= 0;i-- )
        {
            if ( 0 == strcmp( *effectName, Effects[ i ].Name ) )
            {


                for ( j = 0;j < Effects[ i ].NrEff;j++ )
                {
                    for ( k = 0;k < Effects[ i ].Eff[ j ].NrLoci;k++ )
                    {
                        if ( NULL != Effects[ i ].Eff[ j ].Loc[ k ].All )
                        {
                            xFree( Effects[ i ].Eff[ j ].Loc[ k ].All );
                        }
                    }

                    xFree( Effects[ i ].Eff[ j ].Loc );
                }

                if ( NULL != Effects[ i ].Eff )
                {
                    xFree( Effects[ i ].Eff );
                }

                if ( NULL != Effects[ i ].MPset )
                {
                    xFree( Effects[ i ].MPset );
                }


                for ( j = i + 1;j < NrEffects;j++ )
                {
                    Effects[ j - 1 ] = Effects[ j ];
                }

                NrEffects--;

                if ( 0 == NrEffects )
                {
                    xFree( Effects );
                }

                else
                {
                    Effects = ( TEffects* ) xRealloc( Effects, NrEffects, TEffects );
                }

                info( 1, "Effect %s removed\n", *effectName );
                break;
            }
        }
    }
}


void effmap_remove_all()
{
    int i, j, k;

    for ( i = 0;i < NrEffects;i++ )
    {
        for ( j = 0; j < Effects[i].NrEff; j++ )
        {
            for ( k = 0;k < Effects[ i ].Eff[ j ].NrLoci;k++ )
            {
                if ( NULL != Effects[i].Eff[j].Loc[k].All )
                {
                    xFree( Effects[i].Eff[j].Loc[k].All );
                }
            }

            xFree( Effects[ i ].Eff[ j ].Loc );
        }

        if ( NULL != Effects[i].Eff )   xFree( Effects[i].Eff );
        if ( NULL != Effects[i].MPset ) xFree( Effects[i].MPset );

    }

    if ( NULL != Effects ) xFree( Effects );

    if ( 0 < NrEffects )  
      info( 1, "%d effects removed\n", NrEffects );

    NrEffects = 0;
}



void remove_genotype_population ( char **PopNames )
{
    char*PopName;
    PopName = strtok( *PopNames, "[' ','\n']" );

    while ( NULL != PopName )
    {
        TPopListEl * PopListEl;
        find_list_element( PopName, &PopListEl );

        if ( NULL == PopListEl )
            info( -1, "Can't find population %s\n", PopName );
        else
        {
            POPULATION ind;

            for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
                if ( NULL != PopListEl->Pop->Ind[ ind ].genotype )
                    xFree( PopListEl->Pop->Ind[ ind ].genotype );

                PopListEl->Pop->Ind[ ind ].genotype = NULL;
            }

            info( 1, "Genotype removed for population %s\n", PopName );
        }

        PopName = strtok( NULL, "[' ','\n']" );
    }
}


void genotype_population( char **PopNames)
{
    char *PopName;

    TInd *Ind0;
    TPopListEl *PopListEl;

    if ( 0 == cMP )
      {
        info( -2, "No linkage map defined\n" );
	return;
      }

    PopName = strtok( *PopNames, "[' ','\n']" );

    while ( NULL != PopName )
    {

      find_list_element( PopName, &PopListEl );

      if ( NULL == PopListEl )
	info( -1, "Can't find population %s\n", PopName );
      else
        {
	  Ind0 = get_individual( PopListEl->Pop, 0 );

	  if ( Ind0->genotype != NULL ) {
	    info( 1, "Population %s already genotyped - genotyping skipped\n", PopName );
	  }  else {

	    #pragma omp parallel
	    {

	      POPULATION ind;
	      int  chrom, hom;

	      TInd   *Ind;
	      TChrom *Chrom;
	      int    i, m ;

	    #pragma omp for  
	    for ( ind = 0; 
		  ind < PopListEl->Pop->NoInds; 
		  ind++ )
	    {
	      m = 0;

	      for ( chrom = 0; chrom < NoChroms; chrom++ )
		{
		  for ( i = 0; 
			i < PopDescr[chrom].NrMappoints;
			i++, m++ )
		    {
		      Ind = get_individual( PopListEl->Pop, ind );
		      Chrom = &( Ind->Chrom[ chrom ] );
		      
		      if ( NULL == Ind->genotype )
			Ind->genotype = ( int* ) xCalloc( cMPhom, int );
		      else
			Ind->genotype = ( int* ) xRealloc( Ind->genotype, 
							   cMPhom, 
							   int );
		      
		      for ( hom = 0; hom < NoHoms; hom++ )
			{
			  Ind->genotype[m*NoHoms+hom] = 
			    get_allele( &( Chrom->Hom[ hom ] ),
					PopDescr[chrom].Mappoint[i].Pos );
			}
		    }
		}
	    }
	    } //parallel
	    
	    info( 1, "Population %s genotyped\n", PopName );
	  }
        }

      PopName = strtok( NULL, "[' ','\n']" );
    }
}


void evaluate_population ( char**PopNames,
			   char**EffName   )
{
  char       *PopName;
  int        eff = -1;

  TPopListEl *PopListEl;
  TInd       *Ind;
  
  if ( 0 == NrEffects ) 
    {
      info( -2, "No effects defined\n" );
      return;
    }

  PopName = strtok( *PopNames, "[' ','\n']" );

  while ( NULL != PopName )
    {
      find_list_element( PopName, &PopListEl );

      if ( NULL == PopListEl )
	info( -1, "Can't find population %s\n", PopName );
      else
        {
	  Ind = ( TInd* ) get_individual( PopListEl->Pop, 0 );
	  
	  if ( NULL == Ind->genotype )
	    {
	      info( -2, "First genotype population %s\n", PopName );
	      return;
	    }
	  
	  if ( NULL != EffName )
            {
	      TEffects * Effect;
	      
	      if ( NULL == ( Effect = find_Effect( *EffName ) ) ){
		info( -2, "Can't find the effect %s\n", *EffName );
	      	return;
	      }

	      eff = Effect->Pos;
            }


	  #pragma omp parallel 
	  {

	    POPULATION ind;
	    TInd       *Ind;
	    int        i, j, hom, loc, status;


	  #pragma omp for 
	  for ( ind = 0;
		ind < PopListEl->Pop->NoInds;
		ind++ )
          {
	    Ind = ( TInd* ) get_individual( PopListEl->Pop, ind );

	    if ( NULL != Ind->GValue )
	      {
		Ind->GValue = (float*) xRealloc( Ind->GValue, 
						 NrEffects + 1, 
						 float );
	      }
	    else
	      {
		Ind->GValue = (float*) xCalloc( NrEffects + 1, float );
	      }

	    Ind->GValue[0] = 0;
	      
	    for ( i = 1; 
		  i < NrEffects + 1; 
		  i++ )
            {
	      Ind->GValue[i] = Effects[i-1].GValue_avg;

	      for ( j = 0; 
		    j < Effects[i-1].NrEff;
		    j++ )
		{
		  status = 0;
		  loc = 0;

		  while ( ( 0 == status ) 
			  && 
			  ( loc < Effects[i-1].Eff[j].NrLoci ) )
		    {
		      hom = 0;
		      
		    while ( ( 0 == status ) && ( hom < NoHoms ) )
		    {
		      if ( ( -1 < Effects[i-1].Eff[j].Loc[loc].All[hom] ) &&
			   ( Effects[i-1].Eff[j].Loc[loc].All[hom] !=
			     Ind->genotype[Effects[i-1].Eff[j].Loc[loc].MP->counter * NoHoms + hom ] ) )
	              {
			status = 1;
		      }
		      
		      hom++;
		    }
		      
		    loc++;
		  }
		  
		  if ( 0 == status )
		    Ind->GValue[i] += Effects[i-1].Eff[j].Eff;
		}
		
		Ind->GValue[0] += Ind->GValue[i] * Effects[i-1].weight;
	      }
	      
	      if ( eff > -1 )
		Ind->GValue[0] = Ind->GValue[eff];
            }
	  } // parallel
	  
	  info( 1, "Population %s evaluated; selection index: %s\n",
		PopName, ( -1 == eff ) ? "weigthed sum" : *EffName );
        }

        PopName = strtok( NULL, "[' ','\n']" );
    }
}

void remove_evaluate_population ( char **PopNames )
{
  char* PopName;

  PopName = strtok( *PopNames, "[' ','\n']" );

  while ( NULL != PopName )
    {
      TPopListEl * PopListEl;

      find_list_element( PopName, &PopListEl );
      
      if ( NULL == PopListEl )
	info( -1, "Can't find population %s\n", PopName );
      else
        {
	  POPULATION ind;
	  
	  for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
            {
	      if ( NULL != PopListEl->Pop->Ind[ ind ].GValue )
                {
		  xFree( PopListEl->Pop->Ind[ ind ].GValue );
		  PopListEl->Pop->Ind[ ind ].GValue = NULL;
                }
            }
	  
	  info( 1, "Evaluation removed for population %s\n", PopName );
        }
      
      PopName = strtok( NULL, "[' ','\n']" );
    }
}


int indcmp( const void*i, const void*j )
{
    if ( *(float*) ( ( (TInd*) i ) ->GValue ) 
	 > 
	 *(float*) ( ( (TInd*) j ) ->GValue ) )
    {
        return -1;
    }

    if ( *(float*) ( ( (TInd*) i ) ->GValue ) 
	 < 
	 *(float*) ( ( (TInd*) j ) ->GValue ) )
    {
        return 1;
    }

    return 0;
}

int indcmpinv( const void*i, const void*j )
{
    if ( *( float* ) ( ( ( TInd* ) i ) ->GValue ) 
	 < 
	 *( float* ) ( ( ( TInd* ) j ) ->GValue ) )
    {
        return -1;
    }

    if ( *( float* ) ( ( ( TInd* ) i ) ->GValue ) 
	 > 
	 *( float* ) ( ( ( TInd* ) j ) ->GValue ) )
    {
        return 1;
    }

    return 0;
}

void sort_population ( char **PopNames,
		       int  *decreasing )
{
  TPopListEl *PopListEl;
  TInd       *Ind;
  char       *PopName;

  PopName = strtok( *PopNames, " " );

  while ( NULL != PopName )
    {
        find_list_element( PopName, &PopListEl );

        if ( NULL == PopListEl ) {
            info( -2, "Can't find population\n" );
	    return;
	}

        if ( 0 == PopListEl->Pop->NoInds ) {
            info( -2, "Population %s is empty\n", PopName );
	    return;
	}

        Ind = get_individual( PopListEl->Pop, 0 );

        if ( NULL == Ind->GValue ) {
            info( -2, " First evaluate_population\n" );
	    return;
	}

        if ( *decreasing ) 
	  {
	    qsort( PopListEl->Pop->Ind, 
		   PopListEl->Pop->NoInds, 
		   sizeof( TInd ), 
		   indcmp );
	  } 
	else {
	  qsort( PopListEl->Pop->Ind, 
		 PopListEl->Pop->NoInds, 
		 sizeof( TInd ), 
		 indcmpinv );
        }

        PopName = strtok( NULL, " " );
    }
}


int getPos( int all, TSet Set )
{
    int i;

    for ( i = 0;i < Set.betrag;i++ )
    {
        if ( Set.member[ i ] == all )
        {
            return i;
        }
    }

    info( -2, "Fehler in getPos()\n" );
    return 0;
}

SEXP return_population ( SEXP _PopNames, SEXP _bg)
{
  int i_bg = 0;
  SEXP matrix, xAxis, yAxis, matrixDim, matrixDimNames;

  if ( !isString( _PopNames ) || ( length( _PopNames ) != 1 ) )
    {
      info( -2, "Name of populations is not a single string.\n" );
      return R_NilValue;
    }

  if ( !isInteger( _bg ) || ( length( _bg ) != 1 ) )
    {
      info( -2, "Value of background allele is not an integer.\n" );
      return R_NilValue;
    }
  else
    {
        i_bg = INTEGER( _bg ) [ 0 ];
    }

  {
    char * PopNames = 
      (char*) xCalloc( strlen( CHAR( STRING_ELT(_PopNames,0) ) ) + 2, char );
    char tmp[200];
    int* bg = &i_bg;
    unsigned long int col = 0, row = 0;

    char* PopName;
    POPULATION ind, NoInds = 0;
    TPopListEl* PopListEl;
    int chrom, hom,  marker, cmarker, allele;
    unsigned long int NoAll = 0, offset;

    strcpy( PopNames, CHAR( STRING_ELT( _PopNames, 0 ) ) );

    AllSet = ( TSet* ) xCalloc( cMP, TSet );

    for ( marker = 0;marker < cMP;marker++ )
      {
	AllSet[ marker ].member = NULL;
	AllSet[ marker ].betrag = 0;
      }

    PopName = strtok( PopNames, " " );

    while ( NULL != PopName )
      {
	find_list_element( PopName, &PopListEl );
	
	if ( !PopListEl ){
	  info( -2, "Can't find population %s\n", PopName );
	  return R_NilValue;
	}

	if ( 1 > PopListEl->Pop->NoInds ){
	  info( -2, "Empty population %s\n", PopName );
	  return R_NilValue;
	}
	
	
	NoInds += PopListEl->Pop->NoInds;
	
	for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ ) {
	  TInd*Ind = get_individual( PopListEl->Pop, ind );

	  for ( marker = 0; marker < cMP; marker++ ) {
	    for ( hom = 0; hom < NoHoms; hom++ ) {
	      if ( *bg != Ind->genotype[ NoHoms * marker + hom ] )
		insertmember( ( AllSet + marker ), 
			      Ind->genotype[ NoHoms * marker + hom ] );
	    }
	  }
	}

	PopName = strtok( NULL, " " );
      }

    for ( marker = 0; marker < cMP; marker++ )
      {
	if ( 1 < AllSet[ marker ].betrag )
	  {
	    qsort( AllSet[marker].member, 
		   AllSet[marker].betrag, 
		   sizeof(int), 
		   intcmp );
	  }
	NoAll += AllSet[ marker ].betrag;
      }

    PROTECT( matrix = allocMatrix( INTSXP, NoAll, NoInds ) );
    PROTECT( matrixDim = allocVector( INTSXP, 2 ) );
    INTEGER( matrixDim ) [ 0 ] = NoAll;
    INTEGER( matrixDim ) [ 1 ] = NoInds;
    setAttrib( matrix, R_DimSymbol, matrixDim );
    PROTECT( matrixDimNames = allocVector( VECSXP, 2 ) );
    PROTECT( xAxis = allocVector( STRSXP, NoInds ) );
    PROTECT( yAxis = allocVector( STRSXP, NoAll ) );
    for ( offset = 0; offset < NoAll*NoInds; offset++ )
      INTEGER( matrix ) [ offset ] = 0;
    strcpy( PopNames, CHAR( STRING_ELT( _PopNames, 0 ) ) );
    
    PopName = strtok( PopNames, " " );
    offset = 0;

    while ( NULL != PopName )
      {
	find_list_element( PopName, &PopListEl );
	
	if ( !PopListEl ){
	  info( -2, "Can't find population %s\n", PopName );
	  return R_NilValue;
	}

	for ( ind = 0; ind < PopListEl->Pop->NoInds; ind++ ) {
	  TInd*Ind = get_individual( PopListEl->Pop, ind );
	  sprintf( tmp, "%s.%ld ", PopListEl->Pop->Name, ind + 1 );
	  SET_STRING_ELT( xAxis, col++, mkChar( tmp ) );
	    
	  for ( marker = 0; marker < cMP; marker++ ) {
	    if ( NULL == Ind->genotype ){
	      info( -2, "First genotype_population\n" );
	      return R_NilValue;
	    }

	    for ( hom = 0;hom < NoHoms;hom++ ) {
	      if ( *bg != Ind->genotype[ NoHoms * marker + hom ] ) {
        	INTEGER( matrix ) [ offset + 
				    getPos( Ind->genotype[marker*NoHoms+hom], 
					    AllSet[marker] )
				    ] = 1;
	      }
	    }
	    offset += AllSet[ marker ].betrag;
	  }
	}
	PopName = strtok( NULL, " " );
      }

    for ( chrom = 0, marker = 0, row = 0; chrom < NoChroms; chrom++ ) {

      for ( cmarker = 0; 
	    cmarker < PopDescr[chrom].NrMappoints; 
	    cmarker++,marker++ )
	{
	  for ( allele = 0;allele < AllSet[ marker ].betrag;allele++ )
	    {
	      sprintf( tmp, "%s.%d ", 
		       PopDescr[chrom].Mappoint[cmarker].name, 
		       AllSet[marker].member[allele] );
	      SET_STRING_ELT( yAxis, row++, mkChar( tmp ) );
	    }
	}
    }
    for ( marker = 0;marker < cMP;marker++ )
      {
	if ( NULL != AllSet[ marker ].member )
	  xFree( AllSet[ marker ].member );
      }
    
    xFree( AllSet );
    xFree( PopNames );
  }

  SET_VECTOR_ELT( matrixDimNames, 1, xAxis );
  SET_VECTOR_ELT( matrixDimNames, 0, yAxis );
  setAttrib( matrix, R_DimNamesSymbol, matrixDimNames );
  UNPROTECT( 5 );
  return ( matrix );
}



void evaluate_genome ( char**  PopName,
		       int*    all,
		       double* ratio,
		       double* length,
		       double* total,
		       int*    blocks   )
{
  TPopListEl*PopListEl;
  POPULATION ind;
  int chrom, hom, loc;
  double gLength = 0;
  TInd*Ind;

  find_list_element( *PopName, &PopListEl );

  if ( NULL == PopListEl ){
    info( -2, "Can't find Population %s\n", *PopName );
    return;
  }

  for ( chrom = 0;chrom < NoChroms;chrom++ )
    gLength += PopDescr[ chrom ].ChromLen * NoHoms;

  for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ ) {
    Ind = get_individual( PopListEl->Pop, ind);
    for ( chrom = 0; chrom < NoChroms; chrom++ ) {
      for ( hom = 0; hom < NoHoms; hom++ ) {
	for ( loc = 0; loc < Ind->Chrom[chrom].Hom[hom].NoLociUsed; loc++ ) {
	  if ( *all == Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc ].All )
	    {
	      if ( loc + 1 < Ind->Chrom[ chrom ].Hom[ hom ].NoLociUsed )
		{
		  length[ind] += Ind->Chrom[chrom].Hom[hom].Loc[loc+1].Pos -
		                 Ind->Chrom[chrom].Hom[hom].Loc[loc].Pos;
		}
	      else
		{
		  length[ind] += PopDescr[chrom].ChromLen - 
                                 Ind->Chrom[chrom].Hom[hom].Loc[loc].Pos;
		}
	      blocks[ ind ] ++;
	    }
	}
      }
    }
    ratio[ ind ] = length[ ind ] / gLength;
    total[ ind ] = gLength;
  }
}

void evaluate_genome_region ( char**  PopName,
			      int*    all,
			      int*    chrom,
			      double* b,
			      double* e,
			      double* ratio,
			      double* length,
			      double* total,
			      int*    blocks )
{
  TPopListEl*PopListEl;
  POPULATION ind;
  int hom, loc;
  double gLength = 0;
  double begin, end;
  TInd*Ind;
  
  find_list_element( *PopName, &PopListEl );
  
  if ( NULL == PopListEl ){
    info( -2, "Can't find Population %s\n", *PopName );
    return;
  }

  gLength += ( *e - *b ) * NoHoms;
  
  for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
    {
      Ind = get_individual( PopListEl->Pop, ind );
      
      for ( hom = 0;hom < NoHoms;hom++ )
	{
	  for ( loc = 0; loc < Ind->Chrom[*chrom].Hom[hom].NoLociUsed; loc++ )
	    {
	      if ( *all == Ind->Chrom[*chrom].Hom[hom].Loc[loc].All )
		{
		  begin = Ind->Chrom[*chrom].Hom[hom].Loc[loc].Pos;
		  
		  if ( loc + 1 < Ind->Chrom[*chrom].Hom[hom].NoLociUsed )
		    {
		      end = Ind->Chrom[*chrom].Hom[hom].Loc[loc+1].Pos;
		    }
		  else
		    {
		      end = PopDescr[ *chrom ].ChromLen;
		    }
		  
		  if ( ( end >= *b ) && ( begin <= *e ) )
		    {
		      if ( begin < *b ) begin = *b;
		      if ( end   > *e ) end   = *e;
		      length[ ind ] += end - begin;
		      blocks[ ind ] ++;
		    }
		}
	    }
	}
      ratio[ ind ] = length[ ind ] / gLength;
      total[ ind ] = gLength;
    }
}


void reset_mdp()
{
  MDP = 0;
  info(2, "MDP counter reset\n" );
}

void set_mdp( char** cMDP )
{
  MDP = atol( *cMDP );
}

void evaluate_mdp ( char** progeny,
		    char** parent1,
		    char** parent2,
		    char** effectfile )
{
  TEffects   *Effect = NULL;
  POPULATION oldMDP = MDP;
  int i;

  if ( 0 == NrEffects ) {
    info( -2, "No Effectfiles loaded\n" );
    return;
  }

  if ( NULL == effectfile ) {
    info( -2, "No effect specified\n" );
    return;
  }

  if ( ( NULL != effectfile ) 
       && ( NULL == ( Effect = find_Effect(*effectfile) ) ) ) {
    info( -2, "Can't find the effect %s\n", *effectfile );
    return;
  }

  {
    TPopListEl* PopListEl;
    POPULATION  ind;
    int* counter = ( int* ) xCalloc( Effect->MPs, int );
    int* all = ( int* ) xCalloc( Effect->MPs, int );
    int  hom;
    TInd*Ind;

    for ( i = 0; i < Effect->MPs; i++ ) {
      counter[ i ] = 0;
      all[ i ] = -1;
    }

    find_list_element( *parent1, &PopListEl );

    if ( NULL == PopListEl ) {
      info( -2, "Can't find population %s\n", *parent1 );
      return;
    }

    for ( ind = 0; 
	  ind < PopListEl->Pop->NoInds; 
	  ind++ ) {
      Ind = get_individual( PopListEl->Pop, ind );
      for ( i = 0; i < Effect->MPs; i++ ) {
	if ( counter[ i ] < 2 ) {
	  for ( hom = 0;hom < NoHoms;hom++ ) {
	    if ( all[ i ] 
		 != Ind->genotype[ NoHoms * Effect->MPset[i] ->counter + hom ] )
      {
	all[ i ] = Ind->genotype[ NoHoms * Effect->MPset[i] ->counter + hom ];
	counter[ i ] ++;
      }
	  }
	}
      }
    }

    for ( i = 0; i < Effect->MPs; i++ ) {
      if ( 1 == counter[ i ] ) {
	counter[ i ] = 0;
	all[ i ] = -1;
      }
    }

    find_list_element( *parent2, &PopListEl );

    if ( NULL == PopListEl ){
      info( -2, "Can't find population %s\n", *parent2 );
      return;
    }

    for ( ind = 0; ind < PopListEl->Pop->NoInds; ind++ ) {
      Ind = get_individual( PopListEl->Pop, ind );
      for ( i = 0; i < Effect->MPs; i++ ) {
	if ( counter[ i ] < 2 ) {
	  for ( hom = 0; hom < NoHoms; hom++ ) {
	    if ( all[ i ] 
		 != Ind->genotype[ NoHoms * Effect->MPset[ i ] ->counter + hom ] )
     {
       all[ i ] = Ind->genotype[ NoHoms * Effect->MPset[ i ] ->counter + hom ];
       counter[ i ] ++;
     }
	  }
	}
      }
    }
    
    find_list_element( *progeny, &PopListEl );

    if ( NULL == PopListEl ){
      info( -2, "Can't find population %s\n", *progeny );
      return;
    }

    for ( i = 0; i < Effect->MPs; i++ ) {
      if ( counter[ i ] == 2 ) {
	MDP += PopListEl->Pop->NoInds;
      }
    }

    xFree( counter );
    xFree( all );
  }

  info( 1, "%d marker data points required in this step (total %d)\n",
	MDP - oldMDP, MDP );
}


void get_mdp ( char** cmdp )
{
  sprintf( *cmdp, "%lu", MDP );
}

typedef struct {
    int oldPos;
    TMappoint*Mappoint;
} TGenPop;

TGenPop* PGenPop;

int genpopcmp(  const void* i, const void*j  )
{
    int chrom0, chrom1;
    CPOS pos0, pos1;
    chrom0 = ( ( TGenPop* ) i ) ->Mappoint->chrom;
    chrom1 = ( ( TGenPop* ) j ) ->Mappoint->chrom;
    pos0 = ( ( TGenPop* ) i ) ->Mappoint->Pos;
    pos1 = ( ( TGenPop* ) j ) ->Mappoint->Pos;

    if ( chrom0 > chrom1 )
        return 1;

    if ( chrom0 < chrom1 )
        return -1;

    if ( pos0 > pos1 )
        return 1;

    if ( pos0 < pos1 )
        return -1;

    return 0;
}

void generate_population(
    int*  no_pop,
    int*  no_ind,
    int*  no_mar,
    int*  no_all,
    int*  no_pop_u,
    int*  pop_u,
    int*  no_ind_u,
    int*  ind_u,
    int*  no_mar_u,
    int*  mar_u,
    double* data,
    double* freq ,
    char**  pop_list,
    char**  mar_list,
    int*    al_list,
    int*    bg,
    int*    bc )
{
  int i, j, k, l;
  double*ptr_data, *ptr_dataC;
  double count;
  int rowpos_d;
  unsigned long int colpos;
  int nrow_d, nrow_f, ncol_d;
  int*coffset = ( int* ) xCalloc( *no_pop, int );
  int*roff_d = ( int* ) xCalloc( *no_mar, int );
  int*roff_f = ( int* ) xCalloc( *no_mar_u, int );
  
  double freq___ = *freq ; freq___++; 
  /* freq is unused but compiler warnings and change of the
     function call is ommitted by this nonsense assignment*/
  
  coffset[ 0 ] = 0;
  for ( k = 1; k < *no_pop; k++ )
    {
      coffset[ k ] = coffset[ k - 1 ] + no_ind[ k - 1 ];
    }
  
  roff_d[ 0 ] = 0;
  for ( i = 1; i < *no_mar; i++ )
    {
      roff_d[ i ] = roff_d[ i - 1 ] + no_all[ i - 1 ];
    }

  roff_f[ 0 ] = 0;
  for ( i = 1;i < *no_mar_u;i++ )
    {
      roff_f[ i ] = roff_f[ i - 1 ] + no_all[ mar_u[ i - 1 ] - 1 ];
    }
  
    nrow_d = roff_d[ *no_mar - 1 ] + no_all[ *no_mar - 1 ];
    nrow_f = roff_f[ *no_mar_u - 1 ] + no_all[ mar_u[ *no_mar_u - 1 ] - 1 ];
    ncol_d = coffset[ *no_pop - 1 ] + no_ind[ *no_pop - 1 ];

    {
      double*dataC;
      int sum_ind_u = 0;

      dataC = ( double* ) xCalloc( nrow_d * ncol_d, double );

      for ( k = 0;k < *no_pop_u;k++ )
        {
	  for ( l = 0;l < no_ind_u[ k ];l++ )
            {
	      colpos = ( coffset[ pop_u[ k ] - 1 ] + ( ind_u[ l + sum_ind_u ] - 1 ) ) * nrow_d;
	      for ( i = 0;i < *no_mar_u;i++ )
                {
		  count = 0;
		  for ( j = 0;j < no_all[ mar_u[ i ] - 1 ];j++ )
                    {
		      rowpos_d = ( roff_d[ mar_u[ i ] - 1 ] + j );
		      ptr_data = &data[ colpos + rowpos_d ];
		      if ( !ISNA( *ptr_data ) && ( *ptr_data > 0 ) )
                        {
			  *ptr_data = 1;
			  count++;
                        }
                    }
		  for ( j = 0;j < no_all[ mar_u[ i ] - 1 ];j++ )
                    {
		      rowpos_d = ( roff_d[ mar_u[ i ] - 1 ] + j );
		      ptr_data = &data[ colpos + rowpos_d ];
		      ptr_dataC = &dataC[ colpos + rowpos_d ];
		      if ( count > 0 )
			( *ptr_dataC ) = ( *ptr_data ) / ( double ) count;
		      else
			( *ptr_dataC ) = ( *ptr_data );
                    }
                }
            }
	  sum_ind_u += no_ind_u[ k ];
        }
      
      {
	int chrom, hom, tokencounter = 0;
	char*token;
	
	if ( 0 == cMP ){
	  info( -2, "No linkage map defined\n" );
	  return;
	}

	PGenPop = ( TGenPop* ) xCalloc( *no_mar, TGenPop );
	token = strtok( *mar_list, " " );
	
	while ( NULL != token )
	  {
	    PGenPop[ tokencounter ].Mappoint = NULL;
	    for ( chrom = 0;chrom < NoChroms;chrom++ )
	      {
		for ( i = 0;i < PopDescr[ chrom ].NrMappoints;i++ )
		  {
		    if ( 0 == strcmp( PopDescr[ chrom ].Mappoint[ i ].name, token ) )
		      {
			PGenPop[ tokencounter ].oldPos = tokencounter;
			PGenPop[ tokencounter ].Mappoint = &( PopDescr[ chrom ].Mappoint[ i ] );
		      }
		  }
	      }

	    if ( NULL == PGenPop[ tokencounter ].Mappoint ){
	      info( -2, "Can't find marker %s!", token );
	      return;
	    }

	    tokencounter++;
	    
	    token = strtok( NULL, " " );
	  }
	
	qsort( PGenPop, *no_mar, sizeof( TGenPop ), genpopcmp );

	{
	  long int offset = 0;
	  token = strtok( *pop_list, " " );
	  
	  for ( k = 0;k < *no_pop_u;k++ ){
	    int nrow = 0;
	    int ind;
	    int*_ind, *_chrom, *_hom, *_all;
	    double*_pos;

	    _ind = ( int* ) xCalloc( 1, int );
	    _chrom = ( int* ) xCalloc( 1, int );
	    _hom = ( int* ) xCalloc( 1, int );
	    _pos = ( double* ) xCalloc( 1, double );
	    _all = ( int* ) xCalloc( 1, int );
	      
	    for ( ind = 0;ind < no_ind_u[ k ];ind++ ) {
	      for ( chrom = 0;chrom < NoChroms;chrom++ ) {
		int tmpcounter = 0;
		CPOS lastpos = 0;
		for ( i = 0;i < *no_mar;i++ ) {
		  int *candidates;
		  int nr_cand = NoHoms, c, act_cand = 0;
		  candidates = ( int* ) xCalloc( NoHoms, int );
		  if ( PGenPop[ i ].Mappoint->chrom - 1 == chrom ) {
		    for ( c = 0;c < NoHoms;c++ ){
		      candidates[ c ] = *bg;
		    }
		    for ( j = 0;j < no_all[ PGenPop[ i ].oldPos ];j++ ) {
		      while ( dataC[ offset + ind * nrow_d + 
                              roff_d[ PGenPop[ i ].oldPos ] + j ] > 0 ) {
			if ( nr_cand == act_cand ){
			  ++nr_cand;
			  candidates = ( int* ) xRealloc( candidates, 
							  nr_cand, 
							  int );
			}
			candidates [ act_cand++ ] 
			  = ( int ) al_list[ roff_d[ PGenPop[ i ].oldPos ] + j ];
			dataC[ offset + ind * nrow_d + roff_d[ PGenPop[ i ].oldPos ] + j ] 
			  -= ( float ) 1 / NoHoms;
		      }
		    }  
			
		    qsort( candidates, nr_cand, sizeof(int), intcmp );
			
		    for ( hom = 0;hom < NoHoms;hom++ )
		      {
			nrow++;
			    
			_ind = ( int* ) xRealloc( _ind, nrow, int );
			_chrom = ( int* ) xRealloc( _chrom, nrow, int );
			_hom = ( int* ) xRealloc( _hom, nrow, int );
			_pos = ( double* ) xRealloc( _pos, nrow, double );
			_all = ( int* ) xRealloc( _all, nrow, int );
			
			*( _ind + nrow - 1 ) = ind + 1;
			*( _chrom + nrow - 1 ) = chrom + 1;
			*( _hom + nrow - 1 ) = hom + 1;
			  
			
			if ( 0 == tmpcounter ) {
			  * ( _pos + nrow - 1 ) = 0;
			} else {
			  *( _pos + nrow - 1 ) 
			    = ( PGenPop[ i ].Mappoint->Pos + lastpos ) / 2;
			}
			int sel_all;
			if ( 1 == *bc ) {
			  sel_all = hom;
			} else {
			  sel_all = (int) gsl_rng_uniform_int( rng[0], 
							       (unsigned long int)nr_cand );
			}
			nr_cand--;
			*( _all + nrow - 1 ) = candidates[ sel_all ];
			candidates[ sel_all ] = candidates[ nr_cand ];
		      }
		    tmpcounter++;
		    lastpos = PGenPop[ i ].Mappoint->Pos;
			
		    if ( nr_cand > 0 ){
		      info( -1,
			    "More than two alleles. Marker %s Pop: %s Ind: %d\n",
			    PGenPop[ i ].Mappoint->name, token,
			    ind + 1 );
		    }
		  }
		    
		  xFree( candidates );
		}
		  
		if ( 0 == tmpcounter ){
		  if ( 0 == ind ){
		    info( -1, "No marker defined for chromosome %d.\n", chrom + 1 );
		  }
		  for ( hom = 0;hom < NoHoms;hom++ ){
		    nrow++;
		    _ind = ( int* ) xRealloc( _ind, nrow, int );
		    _chrom = ( int* ) xRealloc( _chrom, nrow, int );
		    _hom = ( int* ) xRealloc( _hom, nrow, int );
		    _pos = ( double* ) xRealloc( _pos, nrow, double );
		    _all = ( int* ) xRealloc( _all, nrow, int );
		    *( _ind + nrow - 1 ) = ind + 1;
		    *( _chrom + nrow - 1 ) = chrom + 1;
		    *( _hom + nrow - 1 ) = hom + 1;
		    *( _pos + nrow - 1 ) = 0;
		    *( _all + nrow - 1 ) = *bg;
		  }
		}
		
	      }
	    }
	    
	    offset += ( ind ) * nrow_d;
	    {
	      char cnrow[ 50 ];
	      char*prow = cnrow;
	      
	      sprintf( cnrow, "%d", nrow );
	      init_population( &token, 
			       &prow, 
			       &( _ind[ 0 ] ), 
			       &( _chrom[ 0 ] ), 
			       &( _hom[ 0 ] ), 
			       &( _pos[ 0 ] ), 
			       &( _all[ 0 ] ) );
	    }
	    
	    token = strtok( NULL, " " );
	    nrow = 0;
	    xFree( _ind );
	    xFree( _chrom );
	    xFree( _hom );
	    xFree( _pos );
	    xFree( _all );
	  }
	}
	
	xFree( PGenPop );
      }
      
      xFree( dataC );
    }

    xFree( coffset );
    xFree( roff_d );
    xFree( roff_f );
}


SEXP splitdt( SEXP in_string )
{
  char *ma; 
  char *a;
  int i;
  int n = length( in_string );
  SEXP out_string;
  
  PROTECT( out_string = allocMatrix( STRSXP, n, 2 ) );
  for ( i = 0;i < n;i++ )
    {
      ma = (char*) CHAR( STRING_ELT( in_string, i ) );
      if ( NULL == ( a = strstr( ma, "." ) ) )
	a = strstr( ma, "_" );
      a++;
      SET_STRING_ELT( out_string, n + i, mkChar( a ) );
      a--;
      *a = '\0';
      SET_STRING_ELT( out_string, i, mkChar( ma ) );
      *a = '.';
    }
  UNPROTECT( 1 );
  return ( out_string );
}

void pq_genotype(
    int*nomar,
    int*noind,
    int*mgt,
    int*k
)
{
    int i, j;
    int*p_m = mgt;
    int*p_k = k;

    for ( j = 0;j < *nomar;j++ )
        for ( i = 0;i < *noind;i++ )
        {
            if ( ( *p_m == 0 ) && ( *( 1 + p_m ) == 1 ) )
                * p_k = 0;
            else if ( ( *p_m == 1 ) && ( *( 1 + p_m ) == 1 ) )
                * p_k = 1;
            else if ( ( *p_m == 1 ) && ( *( 1 + p_m ) == 0 ) )
                * p_k = 2;

            p_k++;

            p_m++;

            p_m++;
        }
}

/**************************************************************************
New code for MABC 2010 BEGIN
**************************************************************************/

void evaluate_ld ( 
		  char** PopName,
		  int* all,
		  int* chr,
		  double* pos,
		  double* ld
		  )

/* length of the linkage drag */
{
    TPopListEl*PopListEl;
    POPULATION ind;
    int chrom, hom, loc;
    double gLength = 0;
    TInd*Ind;
    double  upper, lower;

    find_list_element( *PopName, &PopListEl );

    if ( NULL == PopListEl )
      {
        info( -2, "Can't find Population %s\n", *PopName );
	return;
      }

    for ( chrom = 0;chrom < NoChroms;chrom++ )
        gLength += PopDescr[ chrom ].ChromLen * NoHoms;

    for ( ind = 0;ind < PopListEl->Pop->NoInds;ind++ )
    {
        Ind = get_individual( PopListEl->Pop, ind );

        chrom = *chr - 1;
        {
            for ( hom = 0;hom < NoHoms;hom++ )
            {
                for ( loc = 0;loc < Ind->Chrom[ chrom ].Hom[ hom ].NoLociUsed;loc++ )
                {
		  if ( *all == Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc ].All ) 
                    {
                        if ( loc + 1 < Ind->Chrom[ chrom ].Hom[ hom ].NoLociUsed )
                        {
			  upper = Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc + 1 ].Pos;
                        }
                        else
                        {
			  upper = PopDescr[ chrom ].ChromLen;
                        }
			lower = Ind->Chrom[ chrom ].Hom[ hom ].Loc[ loc ].Pos;
			if ( (lower <= *pos) && (*pos < upper) )
			  {
			    if (ld[ ind ] == 0) {         /* first homologue  */
			      ld[ ind ] = upper - lower;
			    } else {                      /* second homologue */
			      ld[ ind ] = ld[ ind ] / 2 + (upper - lower) / 2;
 			    }
			  }
                    }
                }
            }
        }
	/* ld[ ind ] = length[ ind ] ; */
    }
}


/**************************************************************************
New code for MABC 2010 END
**************************************************************************/



/**************************************************************************
Genetic distances BEGIN
**************************************************************************/

#include <R_ext/Random.h> 


static int AllowZeroFrequencies= 0;


void set_allow_zeros(int*allow){
  AllowZeroFrequencies= *allow;
}


void get_allow_zeros(int*allow){
  *allow= AllowZeroFrequencies;
}


SEXP splitdot(SEXP in_string){

  char*ma,*a;
  int i;
  int n= length(in_string);
  SEXP out_string;

  PROTECT(out_string= allocMatrix(STRSXP,n,2));
  for(i= 0;i<n;i++){
    ma= (char*) CHAR(STRING_ELT(in_string,i));
    if(NULL==(a= strstr(ma,".")))
      a= strstr(ma,"_");
    a++;
    SET_STRING_ELT(out_string,n+i,mkChar(a));
    a--;*a= '\0';
    SET_STRING_ELT(out_string,i,mkChar(ma));
    *a= '.';
  }
  UNPROTECT(1);
  return(out_string);
}


void allele_freq(
		 int*no_pop,
		 int*no_ind,
		 int*no_mar,
		 int*no_all,
		 int*no_ind_u,
		 int*ind_u,
		 double*data,
		 double*freq
		 )
{
  int coff[*no_pop];
  int roff[*no_mar];
  int i,j,k,l,rowpos,colpos,nrow;
  double*ptr_data,*ptr_freq,*ptr_dataC;
  double count;

  roff[0]= coff[0]= 0;
  for ( k = 1; k<*no_pop; k++ ) coff[k] = coff[k-1]+no_ind[k-1];
  for ( i = 1; i<*no_mar; i++ ) roff[i] = roff[i-1]+no_all[i-1];
  nrow = roff[*no_mar-1] + no_all[*no_mar-1];

  {

    int ncol_d= coff[*no_pop-1]+no_ind[*no_pop-1];
    double*dataC;
    dataC= (double*)Calloc(nrow*ncol_d,double);
    for(k= 0;k<*no_pop;k++){
      for(l= 0;l<no_ind_u[k];l++){
	colpos= (coff[k]+(ind_u[coff[k]+l]-1))*nrow;
	for(i= 0;i<*no_mar;i++){
	  count= 0;
	  for(j= 0;j<no_all[i];j++){
	    rowpos= roff[i]+j;
	    ptr_data= &data[colpos+rowpos];
	    if(*ptr_data> 0){
	      *ptr_data= 1;
	      count+= 1;
	    }
	  }
	  for(j= 0;j<no_all[i];j++){
	    rowpos= roff[i]+j;
	    ptr_data= &data[colpos+rowpos];
	    ptr_dataC= &dataC[colpos+rowpos];

	    if(AllowZeroFrequencies){*ptr_dataC= *ptr_data;}
	    else{
	      if(count> 0){
		if(*ptr_data> 0)(*ptr_dataC)= (*ptr_data)/(double)count;
		else(*ptr_dataC)= 0;
	      }
	      else{*ptr_dataC= -1;}
	    }
	  }
	}
      }
    }

    for(i= 0;i<*no_mar;i++)
      for(j= 0;j<no_all[i];j++){
	rowpos= roff[i]+j;
	for(k= 0;k<*no_pop;k++){
	  ptr_freq= &freq[k*nrow+rowpos];
	  *ptr_freq= -1;
	  count= 0;
	  for(l= 0;l<no_ind_u[k];l++){
	    colpos= (coff[k]+(ind_u[coff[k]+l]-1))*nrow;
	    ptr_dataC= &dataC[colpos+rowpos];
	    if(*ptr_dataC!=-1){
	      count+= 1;
	      if(*ptr_freq==-1)(*ptr_freq)= 0;
	      (*ptr_freq)+= (*ptr_dataC);
	    }
	  }
	  if(count> 0){
	    if(*ptr_freq> 0)*ptr_freq= *ptr_freq/(double)count;
	    else*ptr_freq= 0;
	  }
	  else*ptr_freq= -1;
	}
      }

    Free(dataC);
  }

}


void correct_missing(
		     int*no_pop,
		     int*no_ind,
		     int*no_mar,
		     int*no_all,
		     int*no_ind_u,
		     int*ind_u,
		     double*data,
		     double*freq
		     )
{

  double freq___ = *freq ; freq___++; 
  /* freq is unused but compiler warnings and change of the
     function call is ommitted by this nonsense assignment*/

  int coff[*no_pop];
  int roff[*no_mar];
  int i,j,k,l,rowpos,colpos,nrow;
  double*ptr_data,*ptr_dataC;
  
  double count;

  roff[0]= coff[0]= 0;
  for(k= 1;k<*no_pop;k++)coff[k]= coff[k-1]+no_ind[k-1];
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];

  {
    int ncol_d = coff[*no_pop-1]+no_ind[*no_pop-1];
    double*dataC;
    dataC= (double*)Calloc(nrow*ncol_d,double);

    for(k= 0;k<*no_pop;k++){
      for(l= 0;l<no_ind_u[k];l++){
	colpos= (coff[k]+(ind_u[coff[k]+l]-1))*nrow;
	for(i= 0;i<*no_mar;i++){
	  count= 0;
	  for(j= 0;j<no_all[i];j++){
	    rowpos= roff[i]+j;
	    ptr_data= &data[colpos+rowpos];
	    if(*ptr_data> 0){
	      *ptr_data= 1;
	      count+= 1;
	    }
	  }
	  for(j= 0;j<no_all[i];j++){
	    rowpos= roff[i]+j;
	    ptr_data= &data[colpos+rowpos];
	    ptr_dataC= &dataC[colpos+rowpos];
	    
	    if(AllowZeroFrequencies){*ptr_dataC= *ptr_data;}
	    else{
	      if(count> 0){
		if(*ptr_data> 0)(*ptr_dataC)= (*ptr_data)/(double)count;
		else(*ptr_dataC)= 0;
	      }
	      else{*ptr_dataC= -1;}
	    }
	  }
	}
      }
    }

    for(k= 0;k<*no_pop;k++){
      for(l= 0;l<no_ind_u[k];l++){
	colpos= (coff[k]+(ind_u[coff[k]+l]-1))*nrow;
	for(i= 0;i<*no_mar;i++){
	  for(j= 0;j<no_all[i];j++){
	    rowpos= roff[i]+j;
	    ptr_data= &data[colpos+rowpos];
	    ptr_dataC= &dataC[colpos+rowpos];
	    (*ptr_data)= (*ptr_dataC);
	  }
	}
      }
    }
 
   Free(dataC);
  }
  
}

void rd(
	int*no_pop,
	int*no_mar,
	int*no_all,
	int*no_mar_u,
	int*mar_u,
	double*freq,
	double*dist
	)
{
  
  long i,j,k1,k2,m,observed,nrow,rowpos;
  double t1,t2;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int roff[*no_mar];
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;

  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      t2= 0;m= 0;
      for(i= 0;i<*no_mar_u;i++){
	observed= 1;
	t1= 0;
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1))
	    t1+= pow((*ptr_freq1-*ptr_freq2),2);
	  else observed= 0;
	}
	if(observed){
	  t2+= sqrt(t1/2);
	  m+= 1;
	}
      }
      (*ptr_dist)= t2/(double)m;
      ptr_dist++;
    }
}


void mrd(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,m,observed,nrow,rowpos;
  double t1;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int roff[*no_mar];
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      t1= 0;m= 0;
      for(i= 0;i<*no_mar_u;i++){
	observed= 1;
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1))
	    t1+= pow((*ptr_freq1-*ptr_freq2),2);
	  else observed= 0;
	}
	if(observed)m+= 1;
      }
      (*ptr_dist)= sqrt(t1/2/(double)m);
      ptr_dist++;
    }
}


void euc(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,m,observed,nrow,rowpos;
  double t1;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int roff[*no_mar];
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;

  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      t1= 0;m= 0;
      for(i= 0;i<*no_mar_u;i++){
	observed= 1;
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1))
	    t1+= pow((*ptr_freq1-*ptr_freq2),2);
	  else observed= 0;
	}
	if(observed)m+= 1;
      }
      (*ptr_dist)= sqrt(t1);
      ptr_dist++;
    }
}


void jac(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,nrow,rowpos;
  double zaehler,nenner,fa,fb;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int roff[*no_mar];
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      zaehler= nenner= 0;
      for(i= 0;i<*no_mar_u;i++){
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1)){
	    fa= (double)(*ptr_freq1> 0?1:0);
	    fb= (double)(*ptr_freq2> 0?1:0);
	    zaehler+= fa*fb;
	    nenner+= fa+fb-fa*fb;
	  }
	}
      }
      if(nenner> 0)
	(*ptr_dist)= zaehler/nenner;
      else
	(*ptr_dist)= 0;
      ptr_dist++;
    }
}


void dic(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,nrow,rowpos;
  double zaehler,nenner,fa,fb;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int roff[*no_mar];
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      zaehler= nenner= 0;
      for(i= 0;i<*no_mar_u;i++){
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1)){
	    fa= (double)(*ptr_freq1> 0?1:0);
	    fb= (double)(*ptr_freq2> 0?1:0);
	    zaehler+= fa*fb;
	    nenner+= fa+fb;
	  }
	}
      }
      if(nenner> 0)
	(*ptr_dist)= 2*zaehler/nenner;
      else
	(*ptr_dist)= 0;
      ptr_dist++;
    }
}


void sma(
	 int*no_pop,
	 int*no_mar,
	 int*no_all,
	 int*no_mar_u,
	 int*mar_u,
	 double*freq,
	 double*dist
	 )
{
  long i,j,k1,k2,nrow,rowpos;
  double zaehler,nenner,fa,fb;
  double*ptr_dist,*ptr_freq1,*ptr_freq2;
  int roff[*no_mar];
  
  roff[0]= 0;
  for(i= 1;i<*no_mar;i++)roff[i]= roff[i-1]+no_all[i-1];
  nrow= roff[*no_mar-1]+no_all[*no_mar-1];
  
  ptr_dist= dist;
  
  for(k1= 0;k1<(*no_pop-1);k1++)
    for(k2= k1+1;k2<(*no_pop);k2++){
      zaehler= nenner= 0;
      for(i= 0;i<*no_mar_u;i++){
	for(j= 0;j<no_all[mar_u[i]-1];j++){
	  rowpos= (roff[mar_u[i]-1]+j);
	  ptr_freq1= &freq[k1*nrow+rowpos];
	  ptr_freq2= &freq[k2*nrow+rowpos];
	  if((*ptr_freq1!=-1)&&(*ptr_freq2!=-1)){
	    fa= (double)(*ptr_freq1> 0?1:0);
	    fb= (double)(*ptr_freq2> 0?1:0);
	    zaehler+= fa*fb+(1-fa)*(1-fb);
	    nenner+= 1;
	  }
	}
      }
      if(nenner> 0)
	(*ptr_dist)= zaehler/nenner;
      else
	(*ptr_dist)= 0;
      ptr_dist++;
    }
}

void gen_dist(
	      int*measure,
	      int*no_pop,
	      int*no_ind,
	      int*no_mar,
	      int*no_all,
	      int*no_pop_u,
	      int*pop_u,
	      int*no_ind_u,
	      int*ind_u,
	      int*no_mar_u,
	      int*mar_u,
	      double*data,
	      double*freq,
	      double*dist,
	      int*s,
	      double*sdev,
	      char**par
	      )
{

  void(*distance)();

  switch(*measure){
    case 1:  distance = rd;  break;
    case 2:  distance = mrd; break;
    case 3:  distance = euc; break;
    case 11: distance = dic; break;
    case 12: distance = jac; break;
    case 13: distance = sma; break;
    default: distance = euc;
    }

  allele_freq(no_pop,no_ind,no_mar,no_all,no_ind_u,ind_u,data,freq);
  (*distance)(no_pop,no_mar,no_all,no_mar_u,mar_u,freq,dist);

  if(*s==1){

    int i,j,k;
    double f= (double)(*no_mar-1)/(double)(*no_mar);
    int n_dist= (pow(*no_pop,2)-*no_pop)/2;
    int no_mar_u[1];
    *no_mar_u= *no_mar-1;

    {
      double dist_j[n_dist],sum_dist_sq[n_dist],sum_dist[n_dist];
      int mar_u[*no_mar_u];
      
      for(j= 0;j<n_dist;j++)
	sum_dist_sq[j]= sum_dist[j]= 0;
      
      for(i= 0;i<*no_mar;i++){
	for(j= k= 0;j<*no_mar;j++)
	  if(i!=j){
	    mar_u[k]= j+1;
	    k++;
	  }
	for(j= 0;j<n_dist;j++)dist_j[j]= 0;

	(*distance)(no_pop,no_mar,no_all,no_mar_u,mar_u,freq,dist_j);

	for(j= 0;j<n_dist;j++){
	  sum_dist_sq[j]+= pow(dist_j[j],2);
	  sum_dist[j]+= dist_j[j];
	}

      }

      for(j= 0;j<n_dist;j++)
	sdev[j] = sqrt(f*(sum_dist_sq[j]-pow(sum_dist[j],2) /
			  (double)(*no_mar)));
      
    }
  }

  if(*s> 1){

    int i,j,k,l,lfdnr,n_dist,sumall= 0,sumind= 0;
    double f= (double)1/(double)(*s);
    for (i = 0; i<(*no_pop_u); i++) sumind += no_ind[pop_u[i]-1];
    n_dist= (pow(*no_pop_u,2)-*no_pop_u)/2;

    {

      int ind_u_boot[sumind],mar_u_boot[*no_mar_u];
      double dist_b[n_dist],sum_dist_sq[n_dist],sum_dist[n_dist];
      
      srand((unsigned)time(NULL));
      
      for(j= 0;j<n_dist;j++)
	sum_dist_sq[j]= sum_dist[j]= 0;
      
      for(i= 0;i<*s;i++){
	if(strstr(*par,"m")){
	  sumall= 0;
	  for(k= 0;k<*no_mar_u;k++){
	    mar_u_boot[k] = 1 + 
	      (int)floor(*no_mar_u*((double)rand())/(RAND_MAX+0.000001));
	    sumall+= no_all[mar_u_boot[k]-1];
	  }
	}
	else{
	  sumall= 0;
	  for(k= 0;k<*no_mar_u;k++){
	    mar_u_boot[k]= mar_u[k];
	    sumall+= no_all[mar_u_boot[k]-1];
	  }
	}
	
	if(strstr(*par,"i")){
	  lfdnr= 0;
	  for(k= 0;k<*no_pop;k++){
	    for(l= 0;l<no_ind[k];l++){
	      ind_u_boot[lfdnr]= 1+
		(int)floor(no_ind[k]*((double)rand())/(RAND_MAX+0.000001));
	      lfdnr++;
	    }
	  }
	}
	else{
	  lfdnr= 0;
	  for(k= 0;k<*no_pop;k++){
	    for(l= 0;l<no_ind[k];l++){
	      ind_u_boot[lfdnr]= ind_u[lfdnr];
	      lfdnr++;
	    }
	  }
	}

	for(j= 0;j<n_dist;j++)dist_b[j]= 0;

	allele_freq(no_pop,no_ind,no_mar,no_all,no_ind_u,ind_u_boot,data,freq);
	(*distance)(no_pop,no_mar,no_all,no_mar_u,mar_u_boot,freq,dist_b);

	for(j= 0;j<n_dist;j++){
	  sum_dist_sq[j]+= pow(dist_b[j],2);
	  sum_dist[j]+= dist_b[j];
	}

      }

      for(j= 0;j<n_dist;j++)
	sdev[j]= sqrt(f*(sum_dist_sq[j]-pow(sum_dist[j],2)/(double)(*s)));
      
    }
  }
}

/**************************************************************************
Genetic distances END
**************************************************************************/

/**************************************************************************
Interface to SelectionTools BEGIN
**************************************************************************/
#define gs_MAXSTR 4096
#define gs_MAXNME 128

typedef unsigned long int  gs_UNSIGNED ;
typedef long double        gs_FLOAT ;

typedef struct
{
    int    chrom;
    double pos;
    char   name [gs_MAXSTR];
    char   class[gs_MAXSTR];
} gs_mappoint_type;


void ps_get_simpop ( char ** pop_name   ,
		     gs_UNSIGNED * NoInd,
		     gs_UNSIGNED * NoMar,
		     char (**IndName) [gs_MAXNME]  ,
		     char (**MarName) [gs_MAXNME]  ,
		     int *(* All1)  ,
		     int *(* All2)  ,
		     gs_UNSIGNED       *  NLoc ,
		     gs_mappoint_type  ** map_S,
		     gs_mappoint_type *** map  ,
		     gs_UNSIGNED       ** mdl  ,
		     int* bg,
		     int* Sretval         )
{

  gs_UNSIGNED i,m;
  int c,l;
  char tmp_str[gs_MAXNME];

  /* Find population */

  TPopListEl *PopListEl;
  TPop *Pop = NULL; 

  find_list_element( *pop_name, &PopListEl );

  if ( PopListEl == NULL ) 
    {
      info( -2, "SIM Population is not defined\n" );
      *Sretval = -2;  return;
    }

  Pop = PopListEl->Pop;
  
  if ( 1 > Pop->NoInds )
    {
      info( -2, "SIM Population is empty\n " );
      *Sretval = -2;  return;
    }

  /* Indidivuals */

  *NoInd  = Pop->NoInds;

  *IndName = realloc ( *IndName,
		       *NoInd * sizeof(char[gs_MAXNME]));

  for (i = 0 ; i< *NoInd; i++) 
    {
      sprintf(tmp_str,"%s.%lu",Pop->Name,1+i) ;
      strcpy( (*IndName)[i] ,tmp_str  );
    }

  /* Markers */ 

  *NoMar = cMP;
  *MarName = realloc( *MarName,
		      *NoMar * sizeof(char[gs_MAXNME]));

  if ( NULL == PopDescr ) 
    {
      info(-2,"SIM No linkage map loaded\n");
      *Sretval = -2;  return;
    }

  m = 0;
  for ( c=0; c < NoChroms; c++ ) 
    {
      for ( l = 0; l < PopDescr[c].NrMappoints; l++ )
	{
	  if ( gs_MAXNME < strlen(PopDescr[c].Mappoint[l].name) ) {
	    info(-2,"SIM Marker name too long\n");
	    *Sretval = -2;  return;
	  }
	  strcpy ( (*MarName)[m] , PopDescr[c].Mappoint[l].name );
	  m++;  
	}
    }

  /* Allocate memory for the genotypic data */

  *All1 =  realloc( *All1, *NoMar * *NoInd * sizeof(int) );
  *All2 =  realloc( *All2, *NoMar * *NoInd * sizeof(int) );

  /* "Genotype population" to a ST data set */

  #pragma omp parallel
  {

    POPULATION ind;
    int  chrom;

    TInd   *Ind;
    TChrom *Chrom;
    int    i, m ;

    gs_UNSIGNED idx;

    #pragma omp for  
    for ( ind = 0; ind < Pop->NoInds; ind++ )
      {
	Ind = get_individual( PopListEl->Pop, ind );

	for ( chrom=0,m=0 ; chrom < NoChroms; chrom++ )
	  {
	    Chrom = &( Ind->Chrom[chrom] );

	    for ( i = 0; i < PopDescr[chrom].NrMappoints; i++,m++ )
	      {

		idx = m * *NoInd + ind;

		(*All1)[idx] =  get_allele( &( Chrom->Hom[0] ),
		                PopDescr[chrom].Mappoint[i].Pos );

		(*All2)[idx] =  get_allele( &( Chrom->Hom[1] ),
  				PopDescr[chrom].Mappoint[i].Pos );

		if ( *bg == (*All1)[idx] ) (*All1)[idx] = -1 ; 
		if ( *bg == (*All2)[idx] ) (*All2)[idx] = -1 ; 

	      }
	  }
      }
  } //parallel

  /* Linkage map */

  *NLoc  = cMP;
  *map_S = realloc ( *map_S, *NLoc*sizeof(gs_mappoint_type) );
  *map   = realloc ( *map,   *NLoc*sizeof(gs_mappoint_type*) );
  *mdl   = realloc ( *mdl,   *NLoc*sizeof(gs_UNSIGNED) );

  for ( c=0,m=0 ; c < NoChroms; c++ ) 
    {
      for ( l = 0; l < PopDescr[c].NrMappoints; l++,m++ )
	{
	  (*map_S)[m].chrom =  (int)    PopDescr[c].Mappoint[l].chrom;
	  (*map_S)[m].pos   =  (double) PopDescr[c].Mappoint[l].Pos*100;
	  strcpy ( (*map_S)[m].name,  PopDescr[c].Mappoint[l].name );
          strcpy ( (*map_S)[m].class, PopDescr[c].Mappoint[l].class );
	}
    }

  /* Array of pointers to mappoints and index in the marker matrix*/
  for ( m=0 ; m < *NLoc ; m++ )
    {
      (*map)[m] = &((*map_S)[m]);
      (*mdl)[m] = m;
    }

  *Sretval = 0;  return;

}

/**************************************************************************
Interface to SelectionTools END
**************************************************************************/
