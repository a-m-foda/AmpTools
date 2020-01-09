//Dec 14th 2019


#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "breakupMomentum.cuh"
#include "barrierFactor.cuh"

///////////////////////////////////////////////////
#include "cuda.h"

#include "GPUUtils/lorentzBoost.cuh"
#include "GPUUtils/threeVector.cuh"
#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

// Macro to ease definition of loops
#define LOOP(INDEX,START,END,INC) for (int INDEX=START;INDEX<=END;INDEX+=INC)
////////////////////////////////////////////////////////////////////////////////////////
 __device__ WCUComplex CZero = { 0, 0 };
 __device__ WCUComplex COne  = { 1, 0 };
 __device__ WCUComplex ic  = { 0, 1 };

 __device__ GDouble DegToRad = PI/180.0;
///////////////////////////////////////////////////////////////////////////////
__global__ void
GPUomegapi_AmpPol_kernel( GPU_AMP_PROTO, int term_sign, GDouble hel_c_0_m_1, GDouble hel_c_1_m_1, GDouble hel_c_1_p_0, GDouble hel_c_1_p_2, GDouble hel_c_2_m_1, GDouble hel_c_2_p_2, GDouble vertex_0_m_m_0, GDouble vertex_1_m_m_m1, GDouble vertex_1_m_m_0, GDouble vertex_1_m_m_p1, GDouble vertex_1_p_m_m1, GDouble vertex_1_p_m_0, GDouble vertex_1_p_m_p1, GDouble vertex_2_m_m_m2, GDouble vertex_2_m_m_m1, GDouble vertex_2_m_m_0, GDouble vertex_2_m_m_p1, GDouble vertex_2_m_m_p2, GDouble vertex_2_p_m_m2, GDouble vertex_2_p_m_m1, GDouble vertex_2_p_m_0, GDouble vertex_2_p_m_p1, GDouble vertex_2_p_m_p2, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction)
{
	int iEvent = GPU_THIS_EVENT;

 GDouble hel_c[3][2][3] = {
			    { {0.0, hel_c_0_m_1, 0.0}, {0.0, 0.0, 0.0} },
			    { {0.0, hel_c_1_m_1, 0.0}, {hel_c_1_p_0, 0.0, hel_c_1_p_2} },
			    { {0.0, hel_c_2_m_1, 0.0}, {0.0, 0.0, hel_c_2_p_2} }
			      };//hel_c_spin_p_odd = 0 & hel_c_spin_m_even = 0 

  GDouble vertex[3][2][2][5]= {//symmetry eq(7) assuming tau_e = +1
			      { { {0.0, 0.0, vertex_0_m_m_0, 0.0, 0.0}, {0.0, 0.0, vertex_0_m_m_0, 0.0, 0.0} },
			      { {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0} } },

			      { { {0.0, vertex_1_m_m_m1, vertex_1_m_m_0, vertex_1_m_m_p1, 0.0}, 
				{0.0, vertex_1_m_m_p1, -vertex_1_m_m_0, vertex_1_m_m_m1, 0.0} }, 
			      { {0.0, vertex_1_p_m_m1, vertex_1_p_m_0, vertex_1_p_m_p1, 0.0}, 
				{0.0, -vertex_1_p_m_p1, vertex_1_p_m_0, -vertex_1_p_m_m1, 0.0} } },

			      { { {vertex_2_m_m_m2, vertex_2_m_m_m1, vertex_2_m_m_0, vertex_2_m_m_p1, vertex_2_m_m_p2},
			      {vertex_2_m_m_p2, -vertex_2_m_m_p1, vertex_2_m_m_0, -vertex_2_m_m_m1, vertex_2_m_m_m2} },
			      { {vertex_2_p_m_m2, vertex_2_p_m_m1, vertex_2_p_m_0, vertex_2_p_m_p1, vertex_2_p_m_p2}, 
			      {-vertex_2_p_m_p2, vertex_2_p_m_p1, -vertex_2_p_m_0, vertex_2_p_m_m1, -vertex_2_p_m_m2} } }
			      };// vertex_spin_parity_photon-helicity_proton-helicity
    

   GDouble cosTheta =  GPU_UVARS(0);
   GDouble Phi =  GPU_UVARS(1);
   GDouble cosThetaH = GPU_UVARS(2);
   GDouble PhiH = GPU_UVARS(3);
   GDouble prod_angle = GPU_UVARS(4);
   GDouble polfrac = GPU_UVARS(5);
   GDouble dalitz_z = GPU_UVARS(6);
   GDouble dalitz_sin3theta = GPU_UVARS(7);


   GDouble Pgamma_corr = polfrac;
   if (term_sign == -1) Pgamma_corr = -polfrac;

///////////////////////////////////////////////////////////////////////////////////////////

   WCUComplex amplitudep = CZero;
   WCUComplex amplitudem = CZero;

   GDouble G = G_SQRT( 1 + 2 * dalitz_alpha * dalitz_z + 2 * dalitz_beta * G_POW(dalitz_z, (double) 3/2) * dalitz_sin3theta
			 + 2 * dalitz_gamma * G_POW(dalitz_z, (double) 2.0) + 2 * dalitz_delta * G_POW(dalitz_z, (double) 5/2) * dalitz_sin3theta );

   for (int i = 0; i <= 1; i++)//spin J_i
   {
      for (int eta = 0; eta <= 1; eta++)//parity eta_i
      {
	WCUComplex term = CZero;
	WCUComplex termp = CZero;
	WCUComplex termm = CZero;
	    
	  for (int Lambda = -i; Lambda <= i; Lambda++)//spin projection of omegapi
	  {
	      for (int lambda = -1; lambda <= 1; lambda++)//omega helicity
	      {
		  //WCUComplex hel_amp = CZero;
                  GDouble hel_amp = 0.0;

		  for(int l = 0; l <= 1; l++)//partial waves (l). Should probably extend to higher orders.
		  {if ( (eta == 0 && l% 2 == 0) || (eta == 1 && l%2 != 0) ) continue;
		   
		  hel_amp += hel_c[i][eta][l] * clebsch(l, 0, 1, lambda, i, lambda);
		  //CPU --> clebschGordan(j1,j2,m1,m2,J,M) || GPU --> clebsch(j1,m1,j2,m2,J,M);
		  }

		  term += wignerD( i, Lambda, lambda, cosTheta, Phi ) * hel_amp * wignerD( 1, lambda, 0, cosThetaH, PhiH ) * G;
		}//loop over lambda
		
		termm = term * vertex[i][eta][0][Lambda];
		termp = term * vertex[i][eta][1][Lambda];

	  }//loop over Lambda
	  
	  amplitudep += termp;
	  amplitudem += termm;
      }//loop over parity
  }//loop over spin

 WCUComplex pol_factorm = { G_COS(  -2 * prod_angle ) , G_SIN( -2 * prod_angle ) };
 WCUComplex pol_factorp = { G_COS(  2 * prod_angle ) , G_SIN( 2 * prod_angle ) };
 WCUComplex amplitude = ( 1/G_SQRT(2.0) ) * ( G_SQRT(1.0+Pgamma_corr) )* (pol_factorp * amplitudep - pol_factorm * amplitudem);

  pcDevAmp[iEvent] = amplitude;
//  pcDevAmp[iEvent] = COne;
}

void
GPUomegapi_AmpPol_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
			int term_sign, GDouble hel_c_0_m_1, GDouble hel_c_1_m_1, GDouble hel_c_1_p_0, GDouble hel_c_1_p_2, GDouble hel_c_2_m_1, GDouble hel_c_2_p_2, GDouble vertex_0_m_m_0, GDouble vertex_1_m_m_m1, GDouble vertex_1_m_m_0, GDouble vertex_1_m_m_p1, GDouble vertex_1_p_m_m1, GDouble vertex_1_p_m_0, GDouble vertex_1_p_m_p1, GDouble vertex_2_m_m_m2, GDouble vertex_2_m_m_m1, GDouble vertex_2_m_m_0, GDouble vertex_2_m_m_p1, GDouble vertex_2_m_m_p2, GDouble vertex_2_p_m_m2, GDouble vertex_2_p_m_m1, GDouble vertex_2_p_m_0, GDouble vertex_2_p_m_p1, GDouble vertex_2_p_m_p2, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction)

{  

  GPUomegapi_AmpPol_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, term_sign, hel_c_0_m_1, hel_c_1_m_1, hel_c_1_p_0, hel_c_1_p_2, hel_c_2_m_1, hel_c_2_p_2, vertex_0_m_m_0, vertex_1_m_m_m1, vertex_1_m_m_0, vertex_1_m_m_p1, vertex_1_p_m_m1, vertex_1_p_m_0, vertex_1_p_m_p1, vertex_2_m_m_m2, vertex_2_m_m_m1, vertex_2_m_m_0, vertex_2_m_m_p1, vertex_2_m_m_p2, vertex_2_p_m_m2, vertex_2_p_m_m1, vertex_2_p_m_0, vertex_2_p_m_p1, vertex_2_p_m_p2, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta, polAngle, polFraction);

}
