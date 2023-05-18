#include <omp.h>
#include <stdio.h>

#include <tchar.h>
#include <iostream>

// EXEMPLO 8
// Exemplo com Redução
int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "TEST OMP!" << std::endl;
	int i, n;
	float a[100], b[100], sum;
	n = 100; /* Some initializations */
	
	for (i=0; i < n; i++)
		a[i] = b[i] = i * 1.0;
	
	sum = 0.0;
	
	#pragma omp parallel for reduction(+:sum)
	for (i=0; i < n; i++)
		sum = sum + (a[i] * b[i]);
	
	printf(" Sum = %f\n",sum);

	return 0;
}