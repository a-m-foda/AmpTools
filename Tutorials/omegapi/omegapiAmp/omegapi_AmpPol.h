//Dec 14th 2019, Based on model by Adam Szczepaniak & Vincent Mathieu
#if !defined(OMEGAPI_AMPPOL)
#define OMEGAPI_AMPPOL

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
void GPUomegapi_AmpPol_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
			     int term_sign, GDouble hel_c_0_m_1, GDouble hel_c_1_m_1, GDouble hel_c_1_p_0, GDouble hel_c_1_p_2, GDouble hel_c_2_m_1, 
			     GDouble hel_c_2_p_2, GDouble vertex_0_m_m_0, GDouble vertex_0_p_m_0, GDouble vertex_1_m_m_m1, GDouble vertex_1_m_m_0, 
			     GDouble vertex_1_m_m_p1, GDouble vertex_1_p_m_m1, GDouble vertex_1_p_m_0, GDouble vertex_1_p_m_p1, GDouble vertex_2_m_m_m2, 
			     GDouble vertex_2_m_m_m1, GDouble vertex_2_m_m_0, GDouble vertex_2_m_m_p1, GDouble vertex_2_m_m_p2, GDouble vertex_2_p_m_m2, 
			     GDouble vertex_2_p_m_m1, GDouble vertex_2_p_m_0, GDouble vertex_2_p_m_p1, GDouble vertex_2_p_m_p2, GDouble dalitz_alpha, 
			     GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction);

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class omegapi_AmpPol : public UserAmplitude< omegapi_AmpPol >
{

public:
  
  omegapi_AmpPol() : UserAmplitude< omegapi_AmpPol >() { }
  omegapi_AmpPol( const vector< string >& args );
  ~omegapi_AmpPol(){}
  
  string name() const { return "omegapi_AmpPol"; }
  
complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
//complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  
  // **********************
  // The following lines are optional and can be used to precalcualte
  // user-defined data that the amplitudes depend on.
  
  // Use this for indexing a user-defined data array and notifying
  // the framework of the number of user-defined variables.
  
  enum UserVars { uv_cosTheta = 0, uv_Phi = 1, uv_cosThetaH = 2, uv_PhiH = 3, uv_prod_angle = 4, uv_Pgamma = 5, uv_dalitz_z = 6, uv_dalitz_sin3theta = 7, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }
  
  // This function needs to be defined -- see comments and discussion
  // in the .cc file.
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
  
  // This is an optional addition if the calcAmplitude routine
  // can run with only the user-defined data and not the original
  // four-vectors.  It is used to optimize memory usage in GPU
  // based fits.
  bool needsUserVarsOnly() const { return true; }
  // **  end of optional lines **
  
 void updatePar( const AmpParameter& par );

#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION

private:
	AmpParameter term_sign;

        AmpParameter hel_c_0_m_1;

	AmpParameter hel_c_1_m_1;
        AmpParameter hel_c_1_p_0;
        AmpParameter hel_c_1_p_2;

        AmpParameter hel_c_2_m_1;
        AmpParameter hel_c_2_p_2;

	AmpParameter vertex_0_m_m_0;
	AmpParameter vertex_0_p_m_0;

	AmpParameter vertex_1_m_m_m1;
	AmpParameter vertex_1_m_m_0;
	AmpParameter vertex_1_m_m_p1;

	AmpParameter vertex_1_p_m_m1;
	AmpParameter vertex_1_p_m_0;
	AmpParameter vertex_1_p_m_p1;

	AmpParameter vertex_2_m_m_m2;
	AmpParameter vertex_2_m_m_m1;
	AmpParameter vertex_2_m_m_0;
	AmpParameter vertex_2_m_m_p1;
	AmpParameter vertex_2_m_m_p2;

	AmpParameter vertex_2_p_m_m2;
	AmpParameter vertex_2_p_m_m1;
	AmpParameter vertex_2_p_m_0;
	AmpParameter vertex_2_p_m_p1;
	AmpParameter vertex_2_p_m_p2;

// 	AmpParameter vertex_0_m_m_0;
// 	AmpParameter vertex_0_m_p_0;
// 
// 	AmpParameter vertex_0_p_m_0;
// 	AmpParameter vertex_0_p_p_0;
// 
// 	AmpParameter vertex_1_m_m_m2;
// 	AmpParameter vertex_1_m_m_m1;
// 	AmpParameter vertex_1_m_m_0;
// 	AmpParameter vertex_1_m_m_p1;
// 	AmpParameter vertex_1_m_m_p2;
// 
// 	AmpParameter vertex_1_m_p_m2;
// 	AmpParameter vertex_1_m_p_m1;
// 	AmpParameter vertex_1_m_p_0;
// 	AmpParameter vertex_1_m_p_p1;
// 	AmpParameter vertex_1_m_p_p2;
// 
// 	AmpParameter vertex_1_p_m_m2;
// 	AmpParameter vertex_1_p_m_m1;
// 	AmpParameter vertex_1_p_m_0;
// 	AmpParameter vertex_1_p_m_p1;
// 	AmpParameter vertex_1_p_m_p2;
// 
// 	AmpParameter vertex_1_p_p_m2;
// 	AmpParameter vertex_1_p_p_m1;
// 	AmpParameter vertex_1_p_p_0;
// 	AmpParameter vertex_1_p_p_p1;
// 	AmpParameter vertex_1_p_p_p2;
// 
// 	AmpParameter vertex_2_m_m_m2;
// 	AmpParameter vertex_2_m_m_m1;
// 	AmpParameter vertex_2_m_m_0;
// 	AmpParameter vertex_2_m_m_p1;
// 	AmpParameter vertex_2_m_m_p2;
// 
// 	AmpParameter vertex_2_m_p_m2;
// 	AmpParameter vertex_2_m_p_m1;
// 	AmpParameter vertex_2_m_p_0;
// 	AmpParameter vertex_2_m_p_p1;
// 	AmpParameter vertex_2_m_p_p2;
// 
// 	AmpParameter vertex_2_p_m_m2;
// 	AmpParameter vertex_2_p_m_m1;
// 	AmpParameter vertex_2_p_m_0;
// 	AmpParameter vertex_2_p_m_p1;
// 	AmpParameter vertex_2_p_m_p2;
// 
// 	AmpParameter vertex_2_p_p_m2;
// 	AmpParameter vertex_2_p_p_m1;
// 	AmpParameter vertex_2_p_p_0;
// 	AmpParameter vertex_2_p_p_p1;
// 	AmpParameter vertex_2_p_p_p2;

	AmpParameter dalitz_alpha;
	AmpParameter dalitz_beta;
	AmpParameter dalitz_gamma;
	AmpParameter dalitz_delta;

	double polAngle, polFraction;
  
  TH1D *totalFlux_vs_E;
  TH1D *polFlux_vs_E;
  TH1D *polFrac_vs_E;

};

#endif
