#include <omp.h>
#include <stdio.h>

#include <tchar.h>
#include <iostream>

// EXEMPLO 1
int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "TEST OMP!" << std::endl;
	omp_set_num_threads( 2 );
	int th_id, nthreads;
	
	#pragma omp parallel private(th_id)
	{
		th_id = omp_get_thread_num();
		printf("Hello World from thread %d\n", th_id);
		#pragma omp barrier
		
		if ( th_id == 0 )
		{
			nthreads = omp_get_num_threads();
			printf("There are %d threads\n",nthreads);
		}
	}
	
	return 0;
}