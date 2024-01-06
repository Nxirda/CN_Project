/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

//Constants
#define TRF 0
#define TRI 1
#define SV  2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr = 0;
  int jj = 0;

  int nbpoints = 10;
  int la = nbpoints -2;
 
  int info = 1;
  int NRHS = 1;
  int IMPLEM = 0;

  double T0 = -5.0;
  double T1 = 5.0;
  double **AAB;


  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }


  printf("--------- Poisson 1D ---------\n\n");
  double* RHS   = (double *) malloc(sizeof(double)*la);
  double* EX_SOL= (double *) malloc(sizeof(double)*la);
  double* X     = (double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  //Number of diagonals in ku and kl, kv is needed for LU decomposition
  int kv=1;
  int ku=1;
  int kl=1;

  //Size of the matrix on the x axis (if col major)
  int lab=kv+kl+ku+1;

  //Actual Matrix
  double *AB = (double *) malloc(sizeof(double)*lab*la);

  //Store AB as a GB matrix
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  //Méthode de validation : erreur avant/erreur arrière
  printf("Solution with LAPACK\n");
  int *ipiv = (int *) calloc(la, sizeof(int));

  /* LU Factorization */
  if(IMPLEM == TRF)
  {
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if(IMPLEM == TRI)
  {
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  /* Solution (Triangular) */
  if(IMPLEM == TRI || IMPLEM == TRF)
  {
    if (info==0)
    {
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0)
      {
        printf("\n INFO DGBTRS = %d\n",info);
      }
      else
      {
        printf("\n INFO = %d\n",info);
      } 
    }
  }
  /*  */

  /* It can also be solved with dgbsv */
  if(IMPLEM == SV)
  {
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
  }

  write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  double relres = relative_forward_error(RHS, EX_SOL, &la);
  // TODO : Compute relative norm of the residual
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);
  printf("\n\n--------- End -----------\n");
}
