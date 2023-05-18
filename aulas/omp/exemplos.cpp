
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


//-------------------------------------------------------------------------


// EXEMPLO 2
// int _tmain(int argc, _TCHAR* argv[])
// {
//   std::cout << "TEST OMP!" << std::endl;
// #pragma omp parallel
//   {
//     printf("The parallel region is executed by thread %d\n",
//       omp_get_thread_num());
//     if ( omp_get_thread_num() == 1 ) {
//       printf(" Thread %d does things differently\n",
//         omp_get_thread_num());
//     }
//   } // End of parallel region
//   return 0;
// }


//-------------------------------------------------------------------------


// EXEMPLO 3
// Definir as variaveis shared e private no loop.
// Em seguida, testar a 'clause' nowait, para quebar o sincronismo.
//int _tmain(int argc, _TCHAR* argv[])
//{
//  std::cout << "TEST OMP!" << std::endl;
//  int i;
//  int n = 12;
//#pragma omp parallel
//#pragma omp for
//  for(i = 0; i < n; i++)
//  {
//    printf("The parallel region is executed by thread %d in the loop iteration %d\n",
//      omp_get_thread_num(), i);
//  } // End of parallel region
//  return 0;
//}


//-------------------------------------------------------------------------


// EXEMPLO 4
// The number of sections controls, and limits, the amount of parallelism.
// If there are n of these code blocks, at most n threads can execute 
// in parallel.
// void funcA()
// {
//   printf("In funcA: this section is executed by thread %d\n",
//     omp_get_thread_num());
// }
// void funcB()
// {
//   printf("In funcB: this section is executed by thread %d\n",
//     omp_get_thread_num());
// }
// int _tmain(int argc, _TCHAR* argv[])
// {
//   std::cout << "TEST OMP!" << std::endl;
//   #pragma omp parallel
//   {
//     #pragma omp sections
//     {
//       #pragma omp section
//         (void) funcA();
//       #pragma omp section
//         (void) funcB();
//     } // End of sections block
//   } // End of parallel region
//   return 0;
// }


//-------------------------------------------------------------------------


// EXEMPLO 5
// Example of the single construct - Only one thread initializes the
// shared variable a.
// int _tmain(int argc, _TCHAR* argv[])
// {
//   std::cout << "TEST OMP!" << std::endl;
//   int i,a;
//   int n = 12;
//   int *b;
//   if ( (b=(int *)malloc(n*sizeof(int))) == NULL )
//     perror("memory allocation for b");
//   #pragma omp parallel shared(a,b) private(i)
//   {
//     #pragma omp single
//     {
//       a = 10;
//       printf("Single construct executed by thread %d\n",
//         omp_get_thread_num());
//     }
//     /* A barrier is automatically inserted here */
//     #pragma omp for
//     for (i=0; i<n; i++)
//       b[i] = a;
//   } /*-- End of parallel region --*/
//   printf("After the parallel region:\n");
//   for (i=0; i<n; i++)
//     printf("b[%d] = %d\n",i,b[i]);
//   free(b);
//   return 0;
// }


//-------------------------------------------------------------------------


// EXEMPLO 6
// Example of the private clause - Each thread has a local copy of
// variables i and a.
// int _tmain(int argc, _TCHAR* argv[])
// {
//   std::cout << "TEST OMP!" << std::endl;
//   int i,a;
//   int n = 6;
//   #pragma omp parallel for private(i,a)
//   for (i=0; i<n; i++)
//   {
//     a = i+1;
//     printf("Thread %d has a value of a = %d for i = %d\n",
//       omp_get_thread_num(),a,i);
//   } // End of parallel region 
// 
//   // This clause makes the sequentially last value of variable a accessible
//   // outside the parallel loop.
//   #pragma omp parallel for private(i) lastprivate(a)
//   for (i=0; i<n; i++)
//   {
//     a = i+1;
//     printf("Thread %d has a value of a = %d for i = %d\n",
//       omp_get_thread_num(),a,i);
//   } // End of parallel for 
//   printf("Value of a after parallel for: a = %d\n",a);
// 
//   return 0;
// }


//-------------------------------------------------------------------------


// EXEMPLO 7
// Exemplo de funções openMP.
// int _tmain(int argc, _TCHAR* argv[])
// {
//   std::cout << "TEST OMP!" << std::endl;
//   int nthreads, tid, procs, maxt,
//     inpar, dynamic, nested;
//   // Start parallel region
//   #pragma omp parallel private(nthreads, tid)
//   {
//     /* Obtain thread number */
//     tid = omp_get_thread_num();
//     /* Only master thread does this */
//     if (tid == 0) 
//     {
//       printf("Thread %d getting info.\n", tid);
//       /* Get environment information */
//       procs = omp_get_num_procs();
//       nthreads = omp_get_num_threads();
//       maxt = omp_get_max_threads();
//       inpar = omp_in_parallel();
//       dynamic = omp_get_dynamic();
//       nested = omp_get_nested();
// 
//       /* Print environment information */
//       printf("Number of processors = %d\n", procs);
//       printf("Number of threads = %d\n", nthreads);
//       printf("Max threads = %d\n", maxt);
//       printf("In parallel? = %d\n", inpar);
//       printf("Dynamic threads? = %d\n", dynamic);
//       printf("Nested parallelism? = %d\n", nested);
//     }
//   } // End of parallel region
// 
//   return 0;
// }


//-------------------------------------------------------------------------


// EXEMPLO 8
// Exemplo com Redução
// int _tmain(int argc, _TCHAR* argv[])
// {
//   std::cout << "TEST OMP!" << std::endl;
//   int i, n;
//   float a[100], b[100], sum;
//   n = 100; /* Some initializations */
//   for (i=0; i < n; i++)
//     a[i] = b[i] = i * 1.0;
//   sum = 0.0;
//   #pragma omp parallel for reduction(+:sum)
//   for (i=0; i < n; i++)
//     sum = sum + (a[i] * b[i]);
//   printf(" Sum = %f\n",sum);
//   return 0;
// }


