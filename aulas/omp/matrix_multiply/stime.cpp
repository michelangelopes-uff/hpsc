#include "stime.h"


STime::STime()
: m_start( clock() )
, m_CPS( double( CLOCKS_PER_SEC ) )
{ 
  m_time = m_start;
}

STime::~STime()
{ 
}

double 
STime::measure()
{
  clock_t old_time = m_time;
  m_time = clock();
  const double seconds = static_cast< double >( m_time - old_time ) / m_CPS;
  return seconds;
}

double 
STime::appendTime()
{
  m_time = clock();
  const double seconds = static_cast< double >( m_time - m_start ) / m_CPS;
  return seconds;
}

double 
STime::measureTotal()
{
  m_time = clock( );
  const double seconds = static_cast< double >( m_time - m_start ) / m_CPS;
  return seconds;
}