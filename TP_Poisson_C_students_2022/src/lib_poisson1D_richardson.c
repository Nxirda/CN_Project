/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

const double PI = 3.14159265358979323846;

void eig_poisson1D(double* eigval, int *la)
{
  for(int i = 0; i < (*la); i++)
  {
    eigval[i] = 2-(2*cos((PI*(i+1))/((*la)+1)));
  }
}

double eigmax_poisson1D(int *la)
{
  double result = 2-(2*cos((PI*(1))/((*la)+1)));
  return result;
}

double eigmin_poisson1D(int *la){
  double result = 2-(2*cos((PI*((*la)-1))/((*la)+1)));
  return result;
}

double richardson_alpha_opt(int *la){
  //double* eigval = (double*)malloc(sizeof(double)*(*la));
  //eig_poisson1D(eigval, la);
  double eig_min = eigmin_poisson1D(la);
  double eig_max = eigmax_poisson1D(la);
  double result = 2/(eig_min + eig_max);
  return result;
}

//By default we can try to validate the model using a value between 0 and 1 as eigen values
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  //Tol = epsilon in the PDF
  //Formula : X[K+1] = X[K] + alpha(b-Ax[k]) 
  int info = 1;
  int *ipiv = (int*) calloc(*la, sizeof(int));
  double* Y = (double*) calloc((*la), sizeof(double));
  cblas_dcopy((*la), RHS, 1, Y, 1);
  
  //R0 = ||b - Ax||/||b||
  double residual = (cblas_dnrm2((*la), RHS, 1));
  //Update b value here
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1, AB, *lab, X, 1, 1, Y, 1);

  double norm_res = (cblas_dnrm2((*la), Y, 1));
  resvec[0] = norm_res/residual;

  while(cblas_dnrm2((*maxit), resvec, 1) > (*tol) && (*nbite) < (*maxit))
  {
    cblas_daxpy((*la), (*alpha_rich), Y, 1, X, 1);
    cblas_dcopy((*la), RHS, 1, Y, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -(*alpha_rich), AB, *lab, X, 1, *alpha_rich, Y, 1);

    norm_res = (cblas_dnrm2((*la), Y , 1));
    (*nbite) ++;
    resvec[(*nbite)] = norm_res/residual;

  } 
  printf("\n%d\n", *nbite);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{

}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{

}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  //On aura besoin de faire une Factorisation LU de M au début 
  //On aura besoin de faire une DGBTRS en début de chaque itération pour calculer 
}

