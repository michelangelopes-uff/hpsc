#include <omp.h>
#include <stdio.h>

#include <tchar.h>
#include <iostream>

// EXEMPLO 5
// Example of the single construct - Only one thread initializes the
// shared variable a.
int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "TEST OMP!" << std::endl;
	int i,a;
	int n = 12;
	int *b;
	
	if ( (b=(int *)malloc(n*sizeof(int))) == NULL )
		perror("memory allocation for b");
	
	#pragma omp parallel shared(a,b) private(i)
	{
		#pragma omp single
		{
			a = 10;
			printf("Single construct executed by thread %d\n",
			omp_get_thread_num());
		}
		/* A barrier is automatically inserted here */
		
		#pragma omp for
		for (i=0; i<n; i++)
			b[i] = a;
	} /*-- End of parallel region --*/

	printf("After the parallel region:\n");
	
	for (i=0; i<n; i++)
		printf("b[%d] = %d\n",i,b[i]);
	
	free(b);

	return 0;
}