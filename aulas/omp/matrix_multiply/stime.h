#include <ctime>
/*!
* @brief This class measures the processing time.
* @author Andre Pereira
*/
class STime
{

private: 
  clock_t m_start;
  clock_t m_time;
  const double m_CPS;

public:
  STime();
  ~STime();
  double measure();
  double appendTime();
  double measureTotal();

};
