#include <omp.h>
#include <stdio.h>

#include <tchar.h>
#include <iostream>

// EXEMPLO 3
// Definir as variaveis shared e private no loop.
// Em seguida, testar a 'clause' nowait, para quebar o sincronismo.
int _tmain(int argc, _TCHAR* argv[])
{
	std::cout << "TEST OMP!" << std::endl;
	int i;
	int n = 12;
	
	#pragma omp parallel
	#pragma omp for
	for(i = 0; i < n; i++)
	{
		printf("The parallel region is executed by thread %d in the loop iteration %d\n",
		omp_get_thread_num(), i);
	} // End of parallel region
	
	return 0;
}