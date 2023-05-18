#include "ptime.h"
#include <omp.h>


PTime::PTime( )
: start( )
{ 
start = omp_get_wtime();
time = start;
}

PTime::~PTime( )
{ 
}

double 
PTime::measure( )
{
double old_time = time;
time = omp_get_wtime();
const double seconds = time - old_time;
return seconds;
}

double 
PTime::appendTime( )
{
time = omp_get_wtime();
const double seconds = time - start;
return seconds;
}

double 
PTime::measureTotal( )
{
time = omp_get_wtime();
const double seconds = time - start;
return seconds;
}

