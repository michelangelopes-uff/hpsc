#include <omp.h>
#include <stdio.h>

#include <tchar.h>
#include <iostream>

// EXEMPLO 6
// Example of the private clause - Each thread has a local copy of
// variables i and a.
int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "TEST OMP!" << std::endl;
	int i,a;
	int n = 6;
	
	#pragma omp parallel for private(i,a)
	for (i=0; i<n; i++)
	{
		a = i+1;
		printf("Thread %d has a value of a = %d for i = %d\n",
		omp_get_thread_num(),a,i);
	} // End of parallel region 

	// This clause makes the sequentially last value of variable a accessible
	// outside the parallel loop.
	#pragma omp parallel for private(i) lastprivate(a)
	for (i=0; i<n; i++)
	{
		a = i+1;
		printf("Thread %d has a value of a = %d for i = %d\n",
		omp_get_thread_num(),a,i);
	} // End of parallel for 
	
	printf("Value of a after parallel for: a = %d\n",a);

	return 0;
}