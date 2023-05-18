#include <omp.h>
#include <stdio.h>

#include <tchar.h>
#include <iostream>

// EXEMPLO 4
// The number of sections controls, and limits, the amount of parallelism.
// If there are n of these code blocks, at most n threads can execute 
// in parallel.
void funcA()
{
  	printf("In funcA: this section is executed by thread %d\n",
    omp_get_thread_num());
}

void funcB()
{
  	printf("In funcB: this section is executed by thread %d\n",
    omp_get_thread_num());
}

int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "TEST OMP!" << std::endl;
	
	#pragma omp parallel
	{
		#pragma omp sections
		{
			#pragma omp section
				(void) funcA();
			#pragma omp section
				(void) funcB();
		} // End of sections block
	} // End of parallel region

	return 0;
}