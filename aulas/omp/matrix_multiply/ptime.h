#include <ctime>
/*!
* @brief This class measures the processing time.
* @author Andre Pereira
*/
class PTime
{

private: 
  double start;
  double time;

public:
  PTime();
  ~PTime();
  double measure();
  double appendTime();
  double measureTotal();

};
