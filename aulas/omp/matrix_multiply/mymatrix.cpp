

#include "mymatrix.h"
#include <math.h>

void
MyMatrix::multiplyByVector(int m, int n, std::vector<double>&  a, std::vector<double>&  b, std::vector<double>&  c)
//Matrix::multiplyByVector(int m, int n, double*  a, double*  b, double*  c)
{
  int i=0,j;
#pragma omp parallel for default(none) \
  shared(m,n,a,b,c) private(i,j)
  //shared(m,n) private(i,j)
  for(i=0;i<m;++i)
  {
    double aa = 0.0;
    for(j=0;j<n;++j)
    {
      aa += b[i*n+j]*c[j];
    }
    a[i] = aa;
  } /*-- End of omp parallel for --*/


}

double
MyMatrix::norm(int m, double* a)
{
  int i;
  double sum;
  sum = 0.0;
  //#pragma omp parallel for reduction(+:sum)
  for(i=0;i<m;++i)
  {
    sum = sum + 0.000001*a[i];
  }
  return sum;
}



