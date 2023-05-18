
#ifdef _OPENMP
#include <omp.h>
#endif

//#include <stdio.h>
#include <iostream>
//#include "mymatrix.h"
//#include <ctime>
//#include <vector>

int main (int argc, char *argv[])
{
  int i,j,m=1000,n=120000;
//  Timer m_timer;	

  std::cout <<"Please give m and n:"<<std::endl;
  scanf_s("%d %d",&m,&n);


  double* a = new double[m];
  double* b = new double[m*n];
  double* c = new double[n];
//  std::vector<double> a(m);
//  std::vector<double> b(m*n);
//  std::vector<double> c(n);

  std::cout <<"Initializing matrix B and vector c"<<std::endl;
  for(j=0;j<n;++j)
    c[j] = 2.0;

  for(i=0;i<m;++i)
    for(j=0;j<n;++j)
      b[i*n+j] = i;


//  Matrix A;

  std::cout <<"Executing multiply function for m = "<<m<<" n = "<<n<<std::endl;
//   m_timer.measure();
//   A.multiplyByVector(m,n,a,b,c);
//   std::cout << "Time: " << m_timer.measure() << " s" << std::endl;


  // shared(a) All threads can read from and write to vector a.
  return 0;
}
