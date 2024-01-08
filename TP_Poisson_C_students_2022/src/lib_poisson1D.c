/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

//
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv)
{
  int src = 0;

  //Create KV columns of zeros
  if(*kv >= 1)
  {
    src = (*kv);
    for(int k = 0; k < (*la); k++)
    {
      for(int i = 0; i < (*kv) ; i++)
      {
        AB[indexABCol(i,k,lab)] = 0.0;
      }
    }
  }

  //Create actual poisson 1D matrix
  for(int k = 0; k < (*la); k++)
  {
    AB[indexABCol(src+1,k,lab)] = 2.0;

    //Lower diagonal
    if(k != 0)
    {
      AB[indexABCol(src,k,lab)] = -1.0;
    }else
    {
      AB[indexABCol(src,k,lab)] = -0.0;
    }
    
    //Upper diagonal
    if(k != (*la)-1)
    {
      AB[indexABCol(src+2,k,lab)] = -1.0;
    }
    else
    {
      AB[indexABCol(src+2,k,lab)] = 0.0;
    }
  }
}

//
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv)
{
  int ii = 0;
  int jj = 0;
  int kk = 0;
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

//
double relative_forward_error(double* x, double* y, int* la)
{
  double x_norm = cblas_dnrm2((*la), x, 1);
  double y_norm = cblas_dnrm2((*la), y, 1);
  
  double result = x_norm - y_norm;

  result = result/x_norm;
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
  if((*la) != (*n) || (*kl) > 1 || (*ku) > 1)
  {
    (*info) = -1;
    return -1;
  }
  else
  {
    //Compute base for upper diagonal
    AB[indexABCol(3,0,lab)] /= AB[indexABCol(2,0,lab)];
    //ipiv[0] = AB[indexABCol(2,0,lab)];

    for(int k = 1; k < (*n); k++)
    {
      //Compute diag values
      AB[indexABCol(2,k,lab)] -= AB[indexABCol(1,k,lab)] * AB[indexABCol(3,k-1,lab)];

      ipiv[k-1] = k;
      //Compute upper diag value
      AB[indexABCol(3,k,lab)] /= AB[indexABCol(2,k,lab)];
    }
    (*info) = 0;
  }
  return *info;
}