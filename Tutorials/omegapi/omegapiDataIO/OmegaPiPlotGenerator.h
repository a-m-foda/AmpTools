#if !(defined OMEGAPIPLOTGENERATOR)
#define OMEGAPIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class OmegaPiPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kOmegaMass = 0, kNumHists};
//mamplitudep_real = 0, mamplitudep_imag = 1, mamplitudem_real = 2, mamplitudem_imag = 3, minetnsity = 4, tamplitudep_real = 5, tamplitudep_imag = 6, tamplitudem_real = 7, tamplitudem_imag = 8, tinetnsity = 9, eamplitudep_real = 10, eamplitudep_imag = 11, eamplitudem_real = 12, eamplitudem_imag = 13, einetnsity = 15, kNumHists};

  OmegaPiPlotGenerator( const FitResults& results );
  OmegaPiPlotGenerator( );
    
  void projectEvent( Kinematics* kin );
 
private:
  
  void createHistograms( );
 
};

#endif

