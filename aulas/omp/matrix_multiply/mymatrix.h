

#include <vector>
#include <math.h>



class MyMatrix
{

public:

  MyMatrix() {};
  ~MyMatrix() {};

  void multiplyByVector(int m, int n, std::vector<double>&  a, std::vector<double>&  b, std::vector<double>&  c);
  //  void multiplyByVector(int m, int n, double*  a, double*  b, double*  c);
  double norm(int m, double* a);

};
