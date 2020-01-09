//Dec 14th 2019, Based on model by Adam Szczepaniak & Vincent Mathieu
#include <ctime>
#include <stdlib.h>
#include <stdio.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
//#include "UTILITIES/CobremsGeneration.hh"
//#include "UTILITIES/BeamProperties.h"

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/AmpParameter.h"
#include "omegapi_AmpPol.h"
#include "barrierFactor.h"
#include "clebschGordan.h"
#include "wignerD.h"
#include "breakupMomentum.h"
#include "omegapiAngles.h"

#include <cmath>
#include <complex>
#include <vector>
#include "TMath.h"

TLorentzVector targetopi2(0,0,0,0.938);

omegapi_AmpPol::omegapi_AmpPol( const vector< string >& args ):
  UserAmplitude< omegapi_AmpPol >( args )
{
	assert( args.size() == (6+17+4+2) || args.size() == (6+17+4+3) );
	
	if(args.size() == (6+17+4+3)){
		polAngle  = atof(args[6+17+4+1].c_str() ); // azimuthal angle of the photon polarization vector in the lab measured in degrees.
		polFraction = AmpParameter( args[6+17+4+2] ); // polarization fraction
		std::cout << "Fixed polarization fraction =" << polFraction << " and pol.angle= " << polAngle << " degrees." << std::endl;
	}
/*     else if (args.size() == (6+17+4+2)){
		// BeamProperties configuration file
		TString beamConfigFile = args[6+17+4+1].c_str();
		BeamProperties beamProp(beamConfigFile);
		polFrac_vs_E = (TH1D*)beamProp.GetPolFrac();
		polAngle = beamProp.GetPolAngle();
		std::cout << "Polarisation angle of " << polAngle << " from BeamProperties." << std::endl;
		if(polAngle == -1)
			std::cout << "This is an amorphous run. Set beam polarisation to 0." << std::endl;
		for(Int_t i=0; i<polFrac_vs_E->GetXaxis()->GetNbins()+2; i++){
			//cout << polFrac_vs_E->GetBinContent(i) << endl;
		}
	}*/
	else
	assert(0);

    //term_sign
    term_sign = AmpParameter(args[0]);
    registerParameter(term_sign);

    //hel_c_spin_parity_partial-wave(l)
    hel_c_0_m_1 = AmpParameter(args[1]);

    hel_c_1_m_1 = AmpParameter(args[2]);
    hel_c_1_p_0 = AmpParameter(args[3]);
    hel_c_1_p_2 = AmpParameter(args[4]);
        
    hel_c_2_m_1 = AmpParameter(args[5]);
    hel_c_2_p_2 = AmpParameter(args[6]);
        
    registerParameter(hel_c_0_m_1);
    
    registerParameter(hel_c_1_m_1);
    registerParameter(hel_c_1_p_0);
    registerParameter(hel_c_1_p_2);

    registerParameter(hel_c_2_m_1);
    registerParameter(hel_c_2_p_2);

   //Implement vertex symmetries
    vertex_0_m_m_0 = AmpParameter(args[6+1]);

    vertex_1_m_m_m1 = AmpParameter(args[6+2]);
    vertex_1_m_m_0 = AmpParameter(args[6+3]);
    vertex_1_m_m_p1 = AmpParameter(args[6+4]);

    vertex_1_p_m_m1 = AmpParameter(args[6+5]);
    vertex_1_p_m_0 = AmpParameter(args[6+6]);
    vertex_1_p_m_p1 = AmpParameter(args[6+7]);

    vertex_2_m_m_m2 = AmpParameter(args[6+8]);
    vertex_2_m_m_m1 = AmpParameter(args[6+9]);
    vertex_2_m_m_0 = AmpParameter(args[6+10]);
    vertex_2_m_m_p1 = AmpParameter(args[6+11]);
    vertex_2_m_m_p2 = AmpParameter(args[6+12]);

    vertex_2_p_m_m2 = AmpParameter(args[6+13]);
    vertex_2_p_m_m1 = AmpParameter(args[6+14]);
    vertex_2_p_m_0 = AmpParameter(args[6+15]);
    vertex_2_p_m_p1 = AmpParameter(args[6+16]);
    vertex_2_p_m_p2 = AmpParameter(args[6+17]);

    registerParameter(vertex_0_m_m_0);

    registerParameter(vertex_1_m_m_m1);
    registerParameter(vertex_1_m_m_0);
    registerParameter(vertex_1_m_m_p1);

    registerParameter(vertex_1_p_m_m1);
    registerParameter(vertex_1_p_m_0);
    registerParameter(vertex_1_p_m_p1);

    registerParameter(vertex_2_m_m_m2);
    registerParameter(vertex_2_m_m_m1);
    registerParameter(vertex_2_m_m_0);
    registerParameter(vertex_2_m_m_p1);
    registerParameter(vertex_2_m_m_p2);

    registerParameter(vertex_2_p_m_m2);
    registerParameter(vertex_2_p_m_m1);
    registerParameter(vertex_2_p_m_0);
    registerParameter(vertex_2_p_m_p1);
    registerParameter(vertex_2_p_m_p2);
    
   //Dalitz Parameters
   dalitz_alpha  = AmpParameter(args[6+17+1]);
   dalitz_beta  = AmpParameter(args[6+17+2]);
   dalitz_gamma  = AmpParameter(args[6+17+3]);
   dalitz_delta  = AmpParameter(args[6+17+4]);

   registerParameter(dalitz_alpha);
   registerParameter(dalitz_beta);
   registerParameter(dalitz_gamma);
   registerParameter(dalitz_delta);

}
////////////////////////////////////////////////// User Vars //////////////////////////////////
void
omegapi_AmpPol::calcUserVars( GDouble** pKin, GDouble* userVars ) const 
{

  TLorentzVector beam  (pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]); 
  TLorentzVector recoil(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);

  TLorentzVector rhos_pip(pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0]);
  TLorentzVector rhos_pim(pKin[5][1], pKin[5][2], pKin[5][3], pKin[5][0]);
  TLorentzVector rho = rhos_pip + rhos_pim;

  TLorentzVector omegas_pi(pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0]);
  TLorentzVector omega = rho + omegas_pi;

  TLorentzVector Xs_pi(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);

  TLorentzVector X = omega + Xs_pi;

    //////////////////////// Boost Particles and Get Angles//////////////////////////////////

  //Helicity coordinate system
  TLorentzVector Gammap = beam + targetopi2;
 
// polarization BeamProperties
	GDouble Pgamma=polFraction;//fixed beam polarization fraction
	if(polAngle == -1)
	Pgamma = 0.;//if beam is amorphous set polarization fraction to 0
	else if(polFrac_vs_E!=NULL){
	int bin = polFrac_vs_E->GetXaxis()->FindBin(beam.E());

	if (bin == 0 || bin > polFrac_vs_E->GetXaxis()->GetNbins()){
		Pgamma = 0.;
	}
	else
	 Pgamma = polFrac_vs_E->GetBinContent(bin);
	}

  //Calculate decay angles in helicity frame
  vector <double> locthetaphi = getomegapiAngles(polAngle, omega, X, beam, Gammap);

  vector <double> locthetaphih = getomegapiAngles(rhos_pip, omega, X, Gammap, rhos_pim);

  userVars[uv_cosTheta] = TMath::Cos(locthetaphi[0]);
  userVars[uv_Phi] = locthetaphi[1];

  userVars[uv_cosThetaH] = TMath::Cos(locthetaphih[0]);
  userVars[uv_PhiH] = locthetaphih[1];

  userVars[uv_prod_angle] = locthetaphi[2];

  userVars[uv_Pgamma] = Pgamma;
  
///////////////////////////////////////////// Dalitz Parameters ///////////////////////////////
  double dalitz_s = rho.M();//s=M2(pip pim)
  double dalitz_t = (rhos_pip+omegas_pi).M2();//t=M2(pip pi0)
  double dalitz_u = (rhos_pim+omegas_pi).M2();//u=M2(pim pi0)
  double m3pi = (2*139.57018)+134.9766;
  double dalitz_d = 2*omega.M()*( omega.M() - m3pi);
  double dalitz_sc = (1/3)*( omega.M2() - rhos_pip.M2() - rhos_pim.M2() - omegas_pi.M2());
  double dalitzx = sqrt(3)*(dalitz_t - dalitz_u)/dalitz_d;
  double dalitzy = 3*(dalitz_sc - dalitz_s)/dalitz_d;
  double dalitz_z = dalitzx*dalitzx + dalitzy*dalitzy;
  double dalitz_sin3theta = TMath::Sin(3 *  TMath::ASin( (dalitzy/sqrt(dalitz_z) )) );
  
  userVars[uv_dalitz_z] = dalitz_z;
  userVars[uv_dalitz_sin3theta] = dalitz_sin3theta;
  
}

////////////////////////////////////////////////// Amplitude Calculation //////////////////////////////////

complex< GDouble >
omegapi_AmpPol::calcAmplitude( GDouble** pKin, GDouble* userVars ) const
{

 GDouble hel_c[3][2][3] = {
			    { {0.0, hel_c_0_m_1, 0.0}, {0.0, 0.0, 0.0} },
			    { {0.0, hel_c_1_m_1, 0.0}, {hel_c_1_p_0, 0.0, hel_c_1_p_2} },
			    { {0.0, hel_c_2_m_1, 0.0}, {0.0, 0.0, hel_c_2_p_2} }
			      };//hel_c_spin_p_odd = 0 & hel_c_spin_m_even = 0 
    
  GDouble vertex[3][2][2][5]= {//symmetry eq(7) assuming tau_e = +1
			      { { {0.0, 0.0, vertex_0_m_m_0, 0.0, 0.0}, {0.0, 0.0, vertex_0_m_m_0, 0.0, 0.0} },
			      { {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0} } },
	  
			      { { {0.0, vertex_1_m_m_m1, vertex_1_m_m_0, vertex_1_m_m_p1, 0.0}, {0.0, vertex_1_m_m_p1, -vertex_1_m_m_0, vertex_1_m_m_m1, 0.0} }, 
			      { {0.0, vertex_1_p_m_m1, vertex_1_p_m_0, vertex_1_p_m_p1, 0.0}, {0.0, -vertex_1_p_m_p1, vertex_1_p_m_0, -vertex_1_p_m_m1, 0.0} } },

			      { { {vertex_2_m_m_m2, vertex_2_m_m_m1, vertex_2_m_m_0, vertex_2_m_m_p1, vertex_2_m_m_p2},
			      {vertex_2_m_m_p2, -vertex_2_m_m_p1, vertex_2_m_m_0, -vertex_2_m_m_m1, vertex_2_m_m_m2} },
			      { {vertex_2_p_m_m2, vertex_2_p_m_m1, vertex_2_p_m_0, vertex_2_p_m_p1, vertex_2_p_m_p2}, 
			      {-vertex_2_p_m_p2, vertex_2_p_m_p1, -vertex_2_p_m_0, vertex_2_p_m_m1, -vertex_2_p_m_m2} } }
			      };// vertex_spin_parity_photon-helicity_proton-helicity
    
  complex <GDouble> COne(1,0);  
  complex <GDouble> CZero(0,0);  

   GDouble cosTheta = userVars[uv_cosTheta];
   GDouble Phi = userVars[uv_Phi];
   GDouble cosThetaH = userVars[uv_cosThetaH];
   GDouble PhiH = userVars[uv_PhiH];
   GDouble prod_angle = userVars[uv_prod_angle];
   GDouble polfrac = userVars[uv_Pgamma];
   GDouble dalitz_z = userVars[uv_dalitz_z];
   GDouble dalitz_sin3theta = userVars[uv_dalitz_sin3theta];

   GDouble Pgamma_corr = polfrac;
   if (term_sign == -1) Pgamma_corr = -polfrac;

   GDouble G = sqrt(1 + 2 * dalitz_alpha * dalitz_z + 2 * dalitz_beta * pow(dalitz_z,3/2) * dalitz_sin3theta
			 + 2 * dalitz_gamma * pow(dalitz_z,2) + 2 * dalitz_delta * pow(dalitz_z,5/2) * dalitz_sin3theta );

   complex <GDouble> amplitudep = CZero;
   complex <GDouble> amplitudem = CZero;

   for (int i = 0; i <= 2; i++)//spin J_i
   {
     
      for (int eta = 0; eta <= 1; eta++)//parity eta_i
      {
	complex <GDouble> term = CZero;
	complex <GDouble> termp = CZero;//lambda_gamma = 1
	complex <GDouble> termm = CZero;//lambda_gamma = -1
	    
	  for (int Lambda = -i; Lambda <= i; Lambda++)//spin projection of omega pi
	  {
	      for (int lambda = -1; lambda <= 1; lambda++)//omega helicity
	      {
		  GDouble hel_amp = 0.0;

		  for(int l = 0; l <= 2; l++)//partial waves (l). Should probably extend to higher orders.
		  {if ( (eta == 0 && l% 2 == 0) || (eta == 1 && l%2 != 0) ) continue;
		  
		  hel_amp += hel_c[i][eta][l] * clebschGordan(l, 1, 0, lambda, i, lambda);
		  }//loop over l
		  
		  term += wignerD( i, Lambda, lambda, cosTheta, Phi ) * hel_amp * wignerD( 1, lambda, 0, cosThetaH, PhiH ) * G;
		}//loop over lambda
		
		termm = term * vertex[i][eta][0][Lambda];
		termp = term * vertex[i][eta][1][Lambda];
	  }//loop over Lambda
	  
	  amplitudep += termp;
	  amplitudem += termm;
      }//loop over parity
  }//loop over spin

 complex <GDouble> pol_factorm (cos(  -2 * prod_angle ), sin( -2 * prod_angle ));
 complex <GDouble> pol_factorp (cos(  2 * prod_angle ), sin( 2 * prod_angle ));

 complex <GDouble> amplitude = (1/sqrt(2)) * (sqrt(1+Pgamma_corr))* (pol_factorp * amplitudep - pol_factorm * amplitudem);

return amplitude;
}

void omegapi_AmpPol::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates  
}

#ifdef GPU_ACCELERATION
void
omegapi_AmpPol::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
    
  GPUomegapi_AmpPol_exec( dimGrid, dimBlock, GPU_AMP_ARGS, 
			  term_sign, hel_c_0_m_1, hel_c_1_m_1, hel_c_1_p_0, hel_c_1_p_2, hel_c_2_m_1, hel_c_2_p_2, vertex_0_m_m_0, 
			  vertex_1_m_m_m1, vertex_1_m_m_0, vertex_1_m_m_p1, vertex_1_p_m_m1, vertex_1_p_m_0, vertex_1_p_m_p1, vertex_2_m_m_m2, 
			  vertex_2_m_m_m1, vertex_2_m_m_0, vertex_2_m_m_p1, vertex_2_m_m_p2, vertex_2_p_m_m2, vertex_2_p_m_m1, vertex_2_p_m_0, 
			  vertex_2_p_m_p1, vertex_2_p_m_p2, dalitz_alpha, dalitz_beta, dalitz_gamma, dalitz_delta, polAngle, polFraction);
}
#endif //GPU_ACCELERATION

