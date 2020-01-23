//#include <ctime>
//#include <stdlib.h>
//#include <stdio.h>

//#include <cassert>
//#include <iostream>
//#include <string>
//#include <sstream>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

//#include "clebschGordan.h"
//#include "wignerD.h"
//#include "breakupMomentum.h"
//#include "omegapiAngles.h"

//#include <cmath>
//#include <complex>
//#include <vector>
//#include "TMath.h"

#include "omegapiDataIO/OmegaPiPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

/* Constructor to display FitResults */
OmegaPiPlotGenerator::OmegaPiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
	createHistograms();
}

/* Constructor for event generator (no FitResult) */
OmegaPiPlotGenerator::OmegaPiPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void OmegaPiPlotGenerator::createHistograms( ) {
  cout << " calls to bookHistogram go here" << endl;
  
   bookHistogram( kOmegaMass, new Histogram1D( 200, 0.2, 0.8, "MOmega", "Invariant Mass of #pi^{0} #gamma") );
/*   bookHistogram( mamplitudep_real, new Histogram1D( 200, 0.8, 3.0, "Re(Amp+) vs M", "Invariant Mass of #omega #pi^{0}") );
   bookHistogram( mamplitudep_imag, new Histogram1D( 200, 0.8, 3.0, "Im(Amp+) vs M", "Invariant Mass of #omega #pi^{0}") );
   bookHistogram( mamplitudem_real, new Histogram1D( 200, 0.8, 3.0, "Re(Amp-) vs M", "Invariant Mass of #omega #pi^{0}") );
   bookHistogram( mamplitudem_imag, new Histogram1D( 200, 0.8, 3.0, "Im(Amp-) vs M", "Invariant Mass of #omega #pi^{0}") );
   bookHistogram( minetnsity, new Histogram1D( 200, 0.8, 3.0, "Intensity vs M", "Invariant Mass of #omega #pi^{0}") );

   bookHistogram( tamplitudep_real, new Histogram1D( 200, 0.0, 4.0, "Re(Amp+) vs Mand-t", "Mand-t") );
   bookHistogram( tamplitudep_imag, new Histogram1D( 200, 0.0, 4.0, "Im(Amp+) vs Mand-t", "Mand-t") );
   bookHistogram( tamplitudem_real, new Histogram1D( 200, 0.0, 4.0, "Re(Amp-) vs Mand-t", "Mand-t") );
   bookHistogram( tamplitudem_imag, new Histogram1D( 200, 0.0, 4.0, "Im(Amp-) vs Mand-t", "Mand-t") );
   bookHistogram( tinetnsity, new Histogram1D( 200, 0.0, 4.0, "Intensity vs Mand-t", "Mand-t") );

   bookHistogram( eamplitudep_real, new Histogram1D( 200, 6.0, 12.0, "Re(Amp+) vs Beam E", "Beam E") );
   bookHistogram( eamplitudep_imag, new Histogram1D( 200, 6.0, 12.0, "Im(Amp+) vs Beam E", "Beam E") );
   bookHistogram( eamplitudem_real, new Histogram1D( 200, 6.0, 12.0, "Re(Amp-) vs Beam E", "Beam E") );
   bookHistogram( eamplitudem_imag, new Histogram1D( 200, 6.0, 12.0, "Im(Amp-) vs Beam E", "Beam E") );
   bookHistogram( einetnsity, new Histogram1D( 200, 6.0, 12.0, "Intensity vs Beam E", "Beam E") );
*/
}

void
OmegaPiPlotGenerator::projectEvent( Kinematics* kin ){

   cout << "project event" << endl;
   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector recoil = kin->particle( 1 );
   TLorentzVector Xs_pi = kin->particle( 2 );//bachelor pi0
   TLorentzVector omegas_pi = kin->particle( 3 );//omega's pi0
   TLorentzVector rhos_pip = kin->particle( 4 );//pi-
   TLorentzVector rhos_pim = kin->particle( 5 );//pi+

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TLorentzVector rho = rhos_pip + rhos_pim;
  TLorentzVector omega = rho + omegas_pi;
  TLorentzVector X = omega + Xs_pi;

  TLorentzVector target(0,0,0,0.938);
//  double polAngle = 0.0;

//  double b1_mass = X.M();
//  double Mandt = fabs((target-recoil).M2());
//  double beamE = beam.E();
  
/*    //////////////////////// Boost Particles and Get Angles//////////////////////////////////

  //Helicity coordinate system
  TLorentzVector Gammap = beam + target;
 
  //Calculate decay angles in helicity frame
  vector <double> locthetaphi = getomegapiAngles(polAngle, omega, X, beam, Gammap);

  vector <double> locthetaphih = getomegapiAngles(rhos_pip, omega, X, Gammap, rhos_pim);

  
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
  
        double Pgamma = 0.4;
        //double polAngle = 0.0;

	double hel_c_0_m_1 = 1.0;

	double hel_c_1_m_1 = 0.0;
        double hel_c_1_p_0 = 0.0;
        double hel_c_1_p_2 = 0.0;

        double hel_c_2_m_1 = 0.0;
        double hel_c_2_p_2 = 0.0;

	double vertex_0_m_m_0 = 1.0;
	double vertex_0_p_m_0 = 1.0;

	double vertex_1_m_m_m1 = 1.0;
	double vertex_1_m_m_0 = 1.0;
	double vertex_1_m_m_p1 = 1.0;

	double vertex_1_p_m_m1 = 1.0;
	double vertex_1_p_m_0 = 1.0;
	double vertex_1_p_m_p1 = 1.0;

	double vertex_2_m_m_m2 = 1.0;
	double vertex_2_m_m_m1 = 1.0;
	double vertex_2_m_m_0 = 1.0;
	double vertex_2_m_m_p1 = 1.0;
	double vertex_2_m_m_p2 = 1.0;

	double vertex_2_p_m_m2 = 1.0;
	double vertex_2_p_m_m1 = 1.0;
	double vertex_2_p_m_0 = 1.0;
	double vertex_2_p_m_p1 = 1.0;
	double vertex_2_p_m_p2 = 1.0;

	double dalitz_alpha = 94.0;
	double dalitz_beta = 0.0;
	double dalitz_gamma = 0.0;
	double dalitz_delta = 0.0;

 GDouble hel_c[3][2][3] = {
			    { {0.0, hel_c_0_m_1, 0.0}, {0.0, 0.0, 0.0} },
			    { {0.0, hel_c_1_m_1, 0.0}, {hel_c_1_p_0, 0.0, hel_c_1_p_2} },
			    { {0.0, hel_c_2_m_1, 0.0}, {0.0, 0.0, hel_c_2_p_2} }
			      };//hel_c_spin_p_odd = 0 & hel_c_spin_m_even = 0 

    
  GDouble vertex[3][2][2][5]= {//symmetry eq(7) assuming tau_e = +1
			      { { {0.0, 0.0, vertex_0_m_m_0, 0.0, 0.0}, {0.0, 0.0, vertex_0_m_m_0, 0.0, 0.0} },
			      { {0.0, 0.0, vertex_0_p_m_0, 0.0, 0.0}, {0.0, 0.0, -vertex_0_p_m_0, 0.0, 0.0} } },

			      { { {0.0, vertex_1_m_m_m1, vertex_1_m_m_0, vertex_1_m_m_p1, 0.0}, {0.0, vertex_1_m_m_p1, -vertex_1_m_m_0, vertex_1_m_m_m1, 0.0} }, 
			      { {0.0, vertex_1_p_m_m1, vertex_1_p_m_0, vertex_1_p_m_p1, 0.0}, {0.0, -vertex_1_p_m_p1, vertex_1_p_m_0, -vertex_1_p_m_m1, 0.0} } },

			      { { {vertex_2_m_m_m2, vertex_2_m_m_m1, vertex_2_m_m_0, vertex_2_m_m_p1, vertex_2_m_m_p2},
			      {vertex_2_m_m_p2, -vertex_2_m_m_p1, vertex_2_m_m_0, -vertex_2_m_m_m1, vertex_2_m_m_m2} },
			      { {vertex_2_p_m_m2, vertex_2_p_m_m1, vertex_2_p_m_0, vertex_2_p_m_p1, vertex_2_p_m_p2}, 
			      {-vertex_2_p_m_p2, vertex_2_p_m_p1, -vertex_2_p_m_0, vertex_2_p_m_m1, -vertex_2_p_m_m2} } }
			      };// vertex_spin_parity_photon-helicity_proton-helicity
    
  complex <GDouble> COne(1,0);  
  complex <GDouble> CZero(0,0);  

   GDouble cosTheta = TMath::Cos(locthetaphi[0]);
   GDouble Phi = locthetaphi[1];
   GDouble cosThetaH = TMath::Cos(locthetaphih[0]);
   GDouble PhiH = locthetaphih[1];
   GDouble prod_angle = locthetaphi[2];
   GDouble polfrac = Pgamma;

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

 complex <GDouble> amplitude1 = (1/sqrt(2)) * (sqrt(1+polfrac))* (pol_factorp * amplitudep - pol_factorm * amplitudem);
 complex <GDouble> amplitude2 = (1/sqrt(2)) * (sqrt(1-polfrac))* (pol_factorp * amplitudep - pol_factorm * amplitudem);

  GDouble amplitude1_real = amplitude1.real();//real(amplitude);
  GDouble amplitude1_imag = amplitude1.imag();//imag(amplitude);
  GDouble amplitude2_real = amplitude2.real();//real(amplitude);
  GDouble amplitude2_imag = amplitude2.imag();//imag(amplitude);
  GDouble intensity = amplitude1_real*amplitude1_real + amplitude1_imag*amplitude1_imag + amplitude2_real*amplitude2_real + amplitude2_imag*amplitude2_imag;//amplitude1.norm() + amplitude2.norm();
*/
   cout << "calls to fillHistogram go here" << endl;
   fillHistogram( kOmegaMass, 1.0 );
/*   fillHistogram( mamplitudep_real, b1_mass);//, amplitude1_real);
   fillHistogram( mamplitudep_imag, b1_mass);//, amplitude1_imag);
   fillHistogram( mamplitudem_real, b1_mass);//, amplitude2_real);
   fillHistogram( mamplitudem_imag, b1_mass);//, amplitude2_imag);
   fillHistogram( minetnsity, b1_mass);//, intensity);

   fillHistogram( tamplitudep_real, Mandt);//, amplitude1_real);
   fillHistogram( tamplitudep_imag, Mandt);//, amplitude1_imag);
   fillHistogram( tamplitudem_real, Mandt);//, amplitude2_real);
   fillHistogram( tamplitudem_imag, Mandt);//, amplitude2_imag);
   fillHistogram( tinetnsity, Mandt);//, intensity);

   fillHistogram( eamplitudep_real, beamE);//, amplitude1_real);
   fillHistogram( eamplitudep_imag, beamE);//, amplitude1_imag);
   fillHistogram( eamplitudem_real, beamE);//, amplitude2_real);
   fillHistogram( eamplitudem_imag, beamE);//, amplitude2_imag);
   fillHistogram( einetnsity, beamE);//, intensity);
*/
}

