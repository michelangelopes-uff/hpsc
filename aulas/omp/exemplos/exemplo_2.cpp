#include <omp.h>
#include <stdio.h>

#include <tchar.h>
#include <iostream>

// EXEMPLO 2
int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "TEST OMP!" << std::endl;
	
	#pragma omp parallel
	{
		printf("The parallel region is executed by thread %d\n",
		omp_get_thread_num());
		
		if ( omp_get_thread_num() == 1 ) {
			printf(" Thread %d does things differently\n",
			omp_get_thread_num());
		}
	} // End of parallel region
	
	return 0;
}