#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>

#include"mkl.h"
#include <time.h>
#include <math.h>



/*----------------Creating Idendity Matrix----------------*/
void create_idendity_matrix(int n, double * I){
	int i, j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i == j){
				I[i + n*j] = 1;
			}
			else{
				I[i + n*j] = 0;
			}
		}
		
	}
}

/*----------------Generalized Inverse----------------*/
//Solving Moore-Pearson Pseudoinverse with SVD
void inverse(double *a, int nrowA, double *inv) {
	
	create_idendity_matrix(nrowA, inv);	
	lapack_int info, n = nrowA;
	lapack_int *rank = malloc(sizeof(lapack_int));
	double *a_copy = malloc(sizeof(double) *nrowA*nrowA);
	double *s = malloc(sizeof(double) *nrowA);

	memcpy(a_copy, a, sizeof(double)*nrowA*nrowA);
	double rcond = -1;

	info =  LAPACKE_dgelss(LAPACK_COL_MAJOR, n, n, n, a_copy, n, inv, n, s, rcond, rank);
	free(rank);
	free(s);
	free(a_copy);
	if(info != 0) printf("Error in function inverse \n");
}

/*----------------Determinant----------------*/
double solve_det_new(double *a, int n) {
	//Calculate determinante of Matrix with LR decomposition
	int i;
	double det = 1;
	lapack_int info, col = n;
	lapack_int *ipiv = malloc(sizeof(lapack_int) *n);
	double *a_copy = malloc(sizeof(double) *n*n);

	memcpy(a_copy, a, sizeof(double)*n*n);

	info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, col, col, a_copy ,col , ipiv);
	for(i = 0; i < n; i++) {
		det = det * a_copy[i*n + i];
	}	

	free(ipiv);
	free(a_copy);
	if(info != 0) printf("Error in fucntion solve_det_new \n");
	
	return det;
}



/*----------------Sclarproduct with self--------------------------*/
double scalar_own(double *x_temp, int *nrecords) {
	int i;
	double scalar = 0;
	
	for (i = 0; i  < *nrecords; i++) {
	 scalar = scalar + x_temp[i]*x_temp[i];
	}
	return scalar;

}

/*----------------Sclarproduct with self--------------------------*/
double sample_var(int *nrecords, double *e) {
	
	int i;
	double vare = 0;
	for(i = 0; i < *nrecords; i++) {
		vare = vare + e[i]*e[i];
		
	}
	
	return vare/rchisq(*nrecords-2);
}

void calc_residuals(double *e, int *nrecords, int *nmarkers, double *x, double *y, double *g, double mu) {
	
	int i, j;
	double sum = 0;
	
	for(i = 0; i < *nrecords; i++) {
		sum=0;
		
		for(j = 0; j < *nmarkers; j++) {
			sum = sum + x[j**nrecords + i] * g[j];
			
		}
		e[i] = y[i] - sum - mu;
	}
}

double sample_mean(int *nrecords, int *nmarkers, double vare, double *x, double * y, double *g) {
	
	int i, k;	
	double mean = 0;
	double sumy = 0;
	double tempx[*nmarkers];

	for(i = 0; i < *nrecords; i++) {
		sumy = sumy + y[i];
	}

	for(k = 0; k < *nmarkers; k++) {
		tempx[k] = 0;
		for(i = 0; i < *nrecords; i++) {
		tempx[k] = tempx[k] + x[i + k**nrecords];
		}
	}
	
	for(i = 0; i < *nmarkers; i++) {
		mean = mean + tempx[i]*g[i];
	}
	
	return rnorm((sumy - mean)/ *nrecords, sqrt(vare/ *nrecords));
}



double calc_lh_now(double *x_temp, double gvar, double *Ival, double *ycorr, int *nrecords, int a) {
	
	int i, j;
	double V[*nrecords**nrecords], V_inv[*nrecords**nrecords];
	double calc_temp[*nrecords], calc_temp2[*nrecords];	
	double calc_temp3 = 0;
	double result;

	if (a == 0) {
		for (i = 0; i < *nrecords**nrecords; i++) {
			V[i] = Ival[i];
		}

	}
	
	else {
		for (i = 0; i < *nrecords; i++) {
			calc_temp[i] = x_temp[i] * gvar;
		}

		for (i = 0; i < *nrecords; i++) {
			for (j = 0; j < *nrecords; j++) {
				V[j+i**nrecords] = calc_temp[j]*x_temp[i];
			}	
		}

		for (i = 0; i < *nrecords**nrecords; i++) {
			V[i] = V[i] + Ival[i];	
		}
	}

	inverse(V, *nrecords, V_inv);

	for (i = 0; i < *nrecords; i++) {
		calc_temp2[i] = 0;
		for (j = 0; j < *nrecords; j++) {
			calc_temp2[i] = calc_temp2[i] + ycorr[j]*V_inv[j + i**nrecords];
		}
	}

	for (i = 0; i < *nrecords; i++) {
		calc_temp3 = calc_temp3 + calc_temp2[i]*ycorr[i];
		
	}

	result = (1/(2*pow(M_PI,0.5**nrecords)*sqrt(solve_det_new(V,*nrecords)))*exp(-0.5*calc_temp3));
	
	return result;
}	       


double calc_meanval(double *x_temp, double *x, double *y, double *gtemp, double mu, 
		double vare, double gvar, int *nrecords, int *nmarkers){
	
	int i, j;
	
	double mult_temp1 = 0;	
	double mult_temp2[*nmarkers];
	double mult_temp3 = 0;
	double mult_temp4 = 0;
	double mult_temp5 = 0;
	double mult_temp6;
	double bracket1, bracket2;

	//-----first bracket-----
	//init
	for (i = 0; i < *nmarkers; i++) {
		mult_temp2[i] = 0;
	}	

	for (i = 0; i < *nrecords; i++) {
		mult_temp1 = mult_temp1 + x_temp[i] * y[i];
		mult_temp4 = mult_temp4 + x_temp[i];
	}
	
	mult_temp4 = mult_temp4*mu;

	for (i = 0; i < *nmarkers; i++) {
		for (j = 0; j < *nrecords; j++) {
			mult_temp2[i] = mult_temp2[i] + x_temp[j]*x[j + i**nrecords];
		}
	}
	
	for (i = 0; i < *nmarkers; i++) {
		mult_temp3 = mult_temp3 + mult_temp2[i]*gtemp[i];	
	}

	
	bracket1 = mult_temp1 - mult_temp3 - mult_temp4;
	
	//-----second bracket-----
	for (i = 0; i < *nrecords; i++) {
		mult_temp5 = mult_temp5 + x_temp[i]*x_temp[i];
	}

	mult_temp6 = vare/gvar;

	bracket2 = mult_temp5+mult_temp6;
	return bracket1/bracket2;
}




/*------------------------BayesB--------------------------*/

void baysB(int *nrecords, int *nmarkers, int *numit, int *numMHCycles, 
	double *propSegs ,double *chi_parameter, double *x, double *y,	double *gS, double *gvarS) {
	
	//-----Decleration and Initialization-----
	int n, i, j, k;
	double vare, gvar_new, LH0, LH1, LH2, mu, alpha, meanval;
 	double g[*nmarkers], gvar[*nmarkers], gtemp[*nmarkers],	Ival[*nrecords**nrecords], 
		ycorr[*nrecords], e[*nrecords], x_temp[*nrecords];

	//Init e, gvar, mu, g
	vare = 0;
	mu = 0.1;
	
	for (i = 0; i < *nrecords**nrecords; i++) {
		Ival[i] = 0;
	}

	for (i = 0; i < *nmarkers; i++){
		gvar[i] = 0.1;
		g[i] = 0.01;
	}
	
	for (i = 0; i < *nrecords; i++) {
		e[i] = 0;
	}
	

	//-----Beginning the iteration-----/
	for(n = 0; n < *numit; n ++) {	
		//seed = seed + n;

		//calculating residuals
		calc_residuals(e, nrecords, nmarkers, x, y, g, mu);

		//sample vare
		vare = sample_var(nrecords, e);
		//seed = seed + 1;
		
		//sample mean
		mu = sample_mean(nrecords, nmarkers, vare, x,  y, g);
		//seed = seed + 1;

		//sample gvar and g with Metropolis Hasting algorithm
			for(j = 0; j < *nmarkers; j++){
				
				for(i = 0; i < *nmarkers; i++) {
					gtemp[i] = g[i];
				}
				
				gtemp[j] = 0;
				for (i = 0; i < *nrecords**nrecords; i++) {
					Ival[i] = 0;
				}
	
				for (i = 0; i < *nrecords; i++) {
					ycorr[i] = 0;
				}

				for(i = 0; i < *nrecords; i++) {
					//init
					ycorr[i] = y[i] - mu;
					Ival[i**nrecords + i] = vare;
					x_temp[i] = x[i + j**nrecords];
					for (k = 0; k < *nmarkers; k++) {
						ycorr[i] = ycorr[i] - x[i+k**nrecords]*gtemp[k];	
					}
				}

				//calclulating likelihood with gvar
				LH1 = calc_lh_now(x_temp, gvar[j], Ival, ycorr, nrecords, 1);
				LH0 = calc_lh_now(x_temp, gvar[j], Ival, ycorr, nrecords, 0);

				for(i = 0; i < *numMHCycles; i++) {
					//seed = seed + 1;
					if (runif(0, 1) < *propSegs) {

						gvar_new = 0.002/rchisq(*chi_parameter);

						//seed = seed + 1;

						LH2 = calc_lh_now(x_temp, gvar_new, Ival, ycorr, nrecords, 2);

						if (LH2/LH1 < 1) {
							alpha = LH2/LH1;
						}

						else {
							alpha = 1;
						}

						//seed = seed + 1;
						if (runif(0, 1) < alpha) {
							gvar[j] = gvar_new;
							LH1 = LH2;
						}
					}

					else {
						if (LH0/LH1 < 1) {
							alpha = LH0/LH1;
						}

						else {
							alpha = 1;
						}

						
						if (runif(0, 1) < alpha) {
							gvar[j] = 0;
							LH1 = LH0;					
						}
					}

				}


				if (gvar[j] > 0) {
					meanval = calc_meanval(x_temp, x, y, gtemp, mu, vare, 
								gvar[j], nrecords, nmarkers);
					//seed = seed + 1;
					g[j] = rnorm(meanval, sqrt((vare)/(scalar_own(x_temp, nrecords)+(vare)/gvar[j])));

				}

				else {
					g[j] = 0;
				
				}

				gS[j**numit + n] = g[j];
				gvarS[j**numit + n] = gvar[j];	
			}
				printf("\r Running BayesB... %.1f %%", (double)(n+1)/(*numit)*100);
			
		}
		printf("\nEffects has been estimated!");
	}

