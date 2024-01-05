/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

//
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv)
{
  for(int i = 0; i < (*la); i++)
  {
    AB[(*kv) +1 +i*(*lab)] = 2.0;

    //Lower diagonal
    if(i != 0)
      AB[((*kv) + i*(*lab))] = -1;

    //Upper diagonal
    if(i != (*la) -1)
      AB[(*kv) +2 + i*(*lab)] = -1;
  } 
}

//
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv)
{
  int ii, jj, kk;
  for (jj=0; jj<(*la); jj++)
  {
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0; ii< *kv; ii++)
      {
	      AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

//
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
  RHS[0] = *BC0;
  RHS[(*la)-1] = *BC1;
  for(int i = 1; i < ((*la)-1); i++)
  {
    RHS[i] = 0.0;
  }
}  

//
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1)
{
  double DELTA = (*BC1) - (*BC0);
  for(int i = 0; i < (*la); i++)
  {
    EX_SOL[i] = (*BC0) + X[i] * DELTA;
  }
}  

//
void set_grid_points_1D(double* X, int* la)
{
  double h = 1.0/(1.0 * ((*la)+1));
  for(int i = 0; i < *la; i++)
  {
    X[i] = (i+1)*h;
  }
}

//USE LAPACK/BLAS HERE
double relative_forward_error(double* x, double* y, int* la)
{
  double x_norm = cblas_dnrm2((*la), x, 1);
  double y_norm = cblas_dnrm2((*la), y, 1);
  
  double result = y_norm - x_norm;

  result = result/y_norm;
  return result;
}

//
int indexABCol(int i, int j, int *lab)
{
  return j*(*lab)+i;
}

//
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  double m;

  if((*la) != (*n) || (*kl) > 1 || (*ku) > 1)
  {
    (*info) = -1;
  }
  else
  {
    (*info) = 0;

    for(int k = 0; k < (*n)-1; k++)
    {
      int i = k+1;
      //printf("A(k,k) is %f\n", AB[indexABCol(2,k,lab)]);
      printf("Iteration [%d]\n", k);
      for(; i < (*n); i++)
      {
        //A(i,k) /= A(k,k)
        //printf("A(i,k) is %f\n", AB[indexABCol(3,k,lab)]);
        AB[indexABCol(3,k, lab)] /= AB[indexABCol(2,k,lab)];
      }
      for (i = k+1; i < (*n); i++) {
        for (int j = k+1; j < (*n); j++) {

            //A(i,j) -= A(i,k) * A(k,j)
            //printf("A(i,j) is %f\n", AB[indexABCol(i,j,lab)]);
            //printf("A(i,k) is %f\n", AB[indexABCol(i,k,lab)]);
            //printf("A(k,j) is %f\n", AB[indexABCol(k,j,lab)]);
            AB[indexABCol(1,j,lab)] -= AB[indexABCol(3,k,lab)] * AB[indexABCol(2,j,lab)];
        }
      } 
    }
  }
  return *info;
}
