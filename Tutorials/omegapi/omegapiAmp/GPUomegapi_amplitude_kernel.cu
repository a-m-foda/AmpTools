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
GPUomegapi_amplitude_kernel( GPU_AMP_PROTO, int sign, int lambda_gamma, int spin, int parity, int helicity, GDouble c_0, GDouble c_1, GDouble c_2, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction)
{
	int iEvent = GPU_THIS_EVENT;

   GDouble cosTheta =  GPU_UVARS(0);
   GDouble Phi =  GPU_UVARS(1);
   GDouble cosThetaH = GPU_UVARS(2);
   GDouble PhiH = GPU_UVARS(3);
   GDouble prod_angle = GPU_UVARS(4);
   GDouble polfrac = GPU_UVARS(5);
   GDouble dalitz_z = GPU_UVARS(6);
   GDouble dalitz_sin3theta = GPU_UVARS(7);

   int hel_c[3] = { c_0, c_1, c_2};

///////////////////////////////////////////////////////////////////////////////////////////

   WCUComplex amplitude = CZero;

   GDouble G = G_SQRT( 1 + (2 * dalitz_alpha * dalitz_z) + (2 * dalitz_beta * G_POW(dalitz_z, (double) 3/2) * dalitz_sin3theta) + (2 * dalitz_gamma * G_POW(dalitz_z, (double) 2.0)) + (2 * dalitz_delta * G_POW(dalitz_z, (double) 5/2) * dalitz_sin3theta) );

	      for (int lambda = -1; lambda <= 1; lambda++)//omega helicity
	      {
		  //WCUComplex hel_amp = CZero;
                  GDouble hel_amp = 0.0;

		  for(int l = 0; l <= 2; l++)//partial waves (l). Should probably extend to higher orders.
		  {//if ( (eta == 0 && l% 2 == 0) || (eta == 1 && l%2 != 0) ) continue;
		   
		  hel_amp += hel_c[l] * clebsch(l, 0, 1, lambda, spin, lambda);
		  //CPU --> clebschGordan(j1,j2,m1,m2,J,M) || GPU --> clebsch(j1,m1,j2,m2,J,M);
		  }

		  amplitude += wignerD( spin, helicity, lambda, cosTheta, Phi ) * hel_amp * wignerD( 1, lambda, 0, cosThetaH, PhiH ) * G;
		}//loop over lambda
		
 WCUComplex prefactor = { G_COS(  lambda_gamma * prod_angle ) , G_SIN( lambda_gamma * prod_angle ) };

 if (sign == -1 && lambda_gamma == -1){amplitude *= -1*G_SQRT( ( 1 + (sign * polfrac) )/2 ) * prefactor;}
 else{amplitude *= G_SQRT( ( 1 + (sign * polfrac) )/2 ) * prefactor;}

  pcDevAmp[iEvent] = amplitude;
//  pcDevAmp[iEvent] = COne;
}

void
GPUomegapi_amplitude_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
			int sign, int lambda_gamma, int spin, int parity, int helicity, GDouble c_0, GDouble c_1, GDouble c_2, GDouble dalitz_alpha, GDouble dalitz_beta, GDouble dalitz_gamma, GDouble dalitz_delta, GDouble polAngle, GDouble polFraction)

{  

  GPUomegapi_amplitude_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, sign, lambda_gamma, spin, parity, helicity, c_0, c_1, c_2, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta, polAngle, polFraction);

}
