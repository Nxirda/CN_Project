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
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;

  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  //double temp, relres;
  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS     =1;
  nbpoints =10;
  la       =nbpoints-2;
  T0       =-5.0;  
  T1       =5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS   = (double *) malloc(sizeof(double)*la);
  EX_SOL= (double *) malloc(sizeof(double)*la);
  X     = (double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  //Number of diagonals in ku and kl, kv is needed for LU decomposition
  kv=1;
  ku=1;
  kl=1;

  //Size of the matrix on the x axis (if col major)
  lab=kv+kl+ku+1;

  //Actual Matrix
  AB = (double *) malloc(sizeof(double)*lab*la);

  //Store AB as a GB matrix
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  //Store result of A*X in RHS
  //cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1, AB, lab, X, 1, 1, RHS, 1);

  //Write the result of dgbmv
  //write_vec(RHS, &la, "RHS.dat");
  //Méthode de validation : erreur avant/erreur arrière

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

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

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB1_Facto.dat");

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
  // TODO : use dgbsv

  write_GB_operator_rowMajor_poisson1D(AB, &la, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  // TODO : Compute relative norm of the residual
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
