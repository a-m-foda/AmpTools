//July 26th 2019, Based on DOI: 10.1016/0550-3213(84)90382-1
#if !defined(OMEGAPIANGAMP)
#define OMEGAPIANGAMP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"

#ifdef GPU_ACCELERATION
void GPUomegapiAngAmp_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
			    GDouble m_1p, GDouble m_w_1p, GDouble m_n_1p, GDouble m_1m, GDouble m_w_1m, GDouble m_n_1m,
                            GDouble m_0m, GDouble m_w_0m, GDouble m_n_0m, GDouble m_ds_ratio, GDouble m_phi0_1p,
                            GDouble m_theta_1p, GDouble m_phip_1p, GDouble m_phim_1p, GDouble m_psi_1p,
			    GDouble m_phi0_1m, GDouble m_theta_1m, GDouble m_phip_1m, GDouble m_phim_1m,
                            GDouble m_psi_1m, GDouble m_phi0_0m, GDouble m_theta_0m, bool useCutoff, GDouble polAngle, GDouble polFraction );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class omegapiAngAmp : public UserAmplitude< omegapiAngAmp >
{

public:
  
  omegapiAngAmp() : UserAmplitude< omegapiAngAmp >() { }
  omegapiAngAmp( const vector< string >& args );
  
  string name() const { return "omegapiAngAmp"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;

 void updatePar( const AmpParameter& par );

#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION

private:
  bool useCutoff;
  AmpParameter m_1p;
  AmpParameter m_w_1p;
  AmpParameter m_n_1p;
  AmpParameter m_phi0_1p;
  AmpParameter m_phip_1p;
  AmpParameter m_phim_1p;
  AmpParameter m_theta_1p;
  AmpParameter m_psi_1p;

  AmpParameter m_1m;
  AmpParameter m_w_1m;
  AmpParameter m_n_1m;
  AmpParameter m_phi0_1m;
  AmpParameter m_phip_1m;
  AmpParameter m_phim_1m;
  AmpParameter m_theta_1m;
  AmpParameter m_psi_1m;

  AmpParameter m_0m;
  AmpParameter m_w_0m;
  AmpParameter m_n_0m;
  AmpParameter m_phi0_0m;
  AmpParameter m_theta_0m;

  AmpParameter m_ds_ratio;

  double polAngle, polFraction;
  
  TH1D *totalFlux_vs_E;
  TH1D *polFlux_vs_E;
  TH1D *polFrac_vs_E;

};

#endif
