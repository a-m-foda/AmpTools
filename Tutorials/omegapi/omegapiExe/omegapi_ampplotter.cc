
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include "omegapiDataIO/ROOTDataReader.h"
#include "omegapiDataIO/ROOTDataReaderBootstrap.h"
#include "omegapiDataIO/ROOTDataReaderWithTCut.h"
#include "omegapiDataIO/ROOTDataReaderTEM.h"
#include "omegapiAmp/omegapiAngAmp.h"
#include "omegapiAmp/omegapi_AmpPol.h"
#include "omegapiAmp/omegapiAngles.h"

// #include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/Kinematics.h"

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){
	
  // set default parameters
  
//   bool useMinos = false;

  string configfile;
  string amp_name;
//   string seedfile;
  
  // parse command line
  
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-c"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-n"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  amp_name = argv[++i]; }
//     if (arg == "-s"){
//       if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
//       else  seedfile = argv[++i]; }
//     if (arg == "-n") useMinos = true;
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
//       cout << "   -n \t\t\t\t\t use MINOS instead of MIGRAD" << endl;
      cout << "   -c <file>\t\t\t\t config file" << endl;
//       cout << "   -s <output file>\t\t\t for seeding next fit based on this fit (optional)" << endl;
      exit(1);}
  }
  
  if (configfile.size() == 0){
    cout << "No config file specified" << endl;
    exit(1);
  }

  ConfigFileParser parser(configfile);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  ReactionInfo* reaction = cfgInfo->reactionList()[0];
  string reactionName(reaction->reactionName());
    
    
  AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
  AmpToolsInterface::registerAmplitude( omegapi_AmpPol() );
  
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
  
  AmpToolsInterface ati( cfgInfo );
  
  TH1D *IntensityM = new TH1D("IntensityM", "Intensity Vs M(#omega #pi); M(#omega #pi); Intensity", 200, 0.0, 3.0);
  TH1D *IntensityMandt = new TH1D("IntensityMandt", "Intensity Vs Mand-t; Mand-t; Intensity", 200, 0.0, 4.0);
  TH1D *IntensityBeamE = new TH1D("IntensityBeamE", "Intensity Vs Beam E; Beam E; Intensity", 200, 5.0, 12.0);
  TH1D *IntensityCosTheta = new TH1D("IntensityCosTheta", "Intensity Vs Cos #theta; Cos #theta; Intensity", 200, -1.0, 1.0);
  TH1D *IntensityPhi = new TH1D("IntensityPhi", "Intensity Vs #phi; #phi ; Intensity", 200, -3.15, 3.15);
  TH1D *IntensityCosThetaH = new TH1D("IntensityCosThetaH", "Intensity Vs Cos #theta_{H}; Cos #theta_{H}; Intensity", 200, -1.0, 1.0);
  TH1D *IntensityPhiH = new TH1D("IntensityPhiH", "Intensity Vs #phi_{H}; #phi_{H} ; Intensity", 200, -3.15, 3.15);
  TH1D *IntensityProdAng = new TH1D("IntensityProdAng", "Intensity Vs #Phi; #Phi ; Intensity", 200, -3.15, 3.15);

  TH1D *ReAmpM = new TH1D("ReAmpM", "Re(Amp) Vs M(#omega #pi); M(#omega #pi); Re(Amp)", 200, 0.0, 3.0);
  TH1D *ReAmpMandt = new TH1D("ReAmpMandt", "Re(Amp) Vs Mand-t; Mand-t; Re(Amp)", 200, 0.0, 4.0);
  TH1D *ReAmpBeamE = new TH1D("ReAmpBeamE", "Re(Amp) Vs Beam E; Beam E; Re(Amp)", 200, 5.0, 12.0);
  TH1D *ReAmpCosTheta = new TH1D("ReAmpCosTheta", "Re(Amp) Vs Cos #theta; Cos #theta; Re(Amp)", 200, -1.0, 1.0);
  TH1D *ReAmpPhi = new TH1D("ReAmpPhi", "Re(Amp) Vs #phi; #phi ; Re(Amp)", 200, -3.15, 3.15);
  TH1D *ReAmpCosThetaH = new TH1D("ReAmpCosThetaH", "Re(Amp) Vs Cos #theta_{H}; Cos #theta_{H}; Re(Amp)", 200, -1.0, 1.0);
  TH1D *ReAmpPhiH = new TH1D("ReAmpPhiH", "Re(Amp) Vs #phi_{H}; #phi_{H} ; Re(Amp)", 200, -3.15, 3.15);
  TH1D *ReAmpProdAng = new TH1D("ReAmpProdAng", "Re(Amp) Vs #Phi; #Phi ; Re(Amp)", 200, -3.15, 3.15);

  TH1D *ImAmpM = new TH1D("ImAmpM", "Im(Amp) Vs M(#omega #pi); M(#omega #pi); Im(Amp)", 200, 0.0, 3.0);
  TH1D *ImAmpMandt = new TH1D("ImAmpMandt", "Im(Amp) Vs Mand-t; Mand-t; Im(Amp)", 200, 0.0, 4.0);
  TH1D *ImAmpBeamE = new TH1D("ImAmpBeamE", "Im(Amp) Vs Beam E; Beam E; Im(Amp)", 200, 5.0, 12.0);
  TH1D *ImAmpCosTheta = new TH1D("ImAmpCosTheta", "Im(Amp) Vs Cos #theta; Cos #theta; Im(Amp)", 200, -1.0, 1.0);
  TH1D *ImAmpPhi = new TH1D("ImAmpPhi", "Im(Amp) Vs #phi; #phi ; Im(Amp)", 200, -3.15, 3.15);
  TH1D *ImAmpCosThetaH = new TH1D("ImAmpCosThetaH", "Im(Amp) Vs Cos #theta_{H}; Cos #theta_{H}; Im(Amp)", 200, -1.0, 1.0);
  TH1D *ImAmpPhiH = new TH1D("ImAmpPhiH", "Im(Amp) Vs #phi_{H}; #phi_{H} ; Im(Amp)", 200, -3.15, 3.15);
  TH1D *ImAmpProdAng = new TH1D("ImAmpProdAng", "Im(Amp) Vs #Phi; #Phi ; Im(Amp)", 200, -3.15, 3.15);

  ati.loadEvents(ati.dataReader( reactionName ));
  ati.processEvents(reactionName);
/** Retrieve the intensity calculated for the specified event.
   *  This requires a prior call to processEvents.
   *
   *  \see clearEvents
   *  \see loadEvent
   *  \see loadEvents
   *  \see processEvents
   */
    int NumofEvents = ati.numEvents();
    double intensity = 0.0;
    complex<double> amp (0.0,0.0);
    double amp_real = 0.0;
    double amp_imag = 0.0;

    string ampname1 = "pone";
    string ampname2 = "mone";

    for ( int i= 0; i < NumofEvents; i++)
    {
    Kinematics* kin = ati.kinematics(i);

   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector recoil = kin->particle( 1 );
   TLorentzVector Xs_pi = kin->particle( 2 );//bachelor pi0
   TLorentzVector omegas_pi = kin->particle( 3 );//omega's pi0
   TLorentzVector rhos_pip = kin->particle( 4 );//pi+
   TLorentzVector rhos_pim = kin->particle( 5 );//pi-

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      TLorentzVector rho = rhos_pip + rhos_pim;
      TLorentzVector omega = rho + omegas_pi;
      TLorentzVector X = omega + Xs_pi;
    
      TLorentzVector target(0,0,0,0.938);
    
      TLorentzVector Gammap = beam + target;

      double b1_mass = X.M();
      double Mandt = fabs((target-recoil).M2());
      double beamE = beam.E();

  double polAngle = 90.0;
	
  vector <double> locthetaphi = getomegapiAngles(polAngle, omega, X, beam, Gammap);

  vector <double> locthetaphih = getomegapiAngles(rhos_pip, omega, X, Gammap, rhos_pim);

  double cosTheta = TMath::Cos(locthetaphi[0]);
  double Phi = locthetaphi[1];
 
  double cosThetaH = TMath::Cos(locthetaphih[0]);
  double PhiH = locthetaphih[1];
 
  double prod_angle = locthetaphi[2];
             
    intensity = ati.intensity(i);//int iEvent, unsigned int iDataSet = 0) const;i

    IntensityM->Fill(b1_mass, intensity);
    IntensityMandt->Fill(Mandt, intensity);
    IntensityBeamE->Fill(beamE, intensity);
    IntensityCosTheta->Fill(cosTheta, intensity);
    IntensityPhi->Fill(Phi, intensity);
    IntensityCosThetaH->Fill(cosThetaH, intensity);
    IntensityPhiH->Fill(PhiH, intensity);
    IntensityProdAng->Fill(prod_angle, intensity);

    
//    amp  = ati.decayAmplitude(i, "pone");//int iEvent, string ampName, unsigned int iDataSet = 0) const;

//    amp_real = std::real(intensity);
//    amp_imag = std::imag(intensity);
    
//    cout << "amp = " << std::real(intensity) << " +i " << std::imag(intensity) << endl;
//.real() << " +i " << intensity.imag() << endl;
    }    

     char STRING[256];
     printf(STRING,"Intesnity for partial wave %01s", amp_name);
     TCanvas *c1 = new TCanvas(STRING, STRING, 800, 1200);
//	c1->Divide(2,4);
      c1->Divide(2,2);

        c1->cd(1);
        IntensityM->Draw("e1");

/*        c1->cd(2);
        IntensityMandt->Draw("e1");

        c1->cd(3);
        IntensityBeamE->Draw("e1");

	c1->cd(4);
	IntensityCosTheta->Draw("e1");
*/
//        c1->cd(5);
        c1->cd(2);
        IntensityPhi->Draw("e1");

//        c1->cd(6);
	c1->cd(3);        
        IntensityCosThetaH->Draw("e1");

//        c1->cd(7);
//        IntensityPhiH->Draw("e1");

//        c1->cd(8);
        c1->cd(4);
        IntensityProdAng->Draw("e1");

//	char*[256] STRING;
	sprintf(STRING,"fitresults_%01s.C", amp_name);
	c1->SaveAs(STRING);
        sprintf(STRING,"fitresults_%01s.pdf", amp_name);
        c1->SaveAs(STRING);
   
/*     TCanvas *c2 = new TCanvas("ReAmp", "ReAmp", 800, 1200);
	c2->Divide(2,4);

        c2->cd(1);
        ReAmpM->Draw("e1");

        c2->cd(2);
        ReAmpMandt->Draw("e1");

        c2->cd(3);
        ReAmpBeamE->Draw("e1");

	c2->cd(4);
	ReAmpCosTheta->Draw("e1");

        c2->cd(5);
        ReAmpPhi->Draw("e1");

        c2->cd(6);
        ReAmpCosThetaH->Draw("e1");

        c2->cd(7);
        ReAmpPhiH->Draw("e1");

        c2->cd(8);
        ReAmpProdAng->Draw("e1");

	c2->SaveAs("int.pdf");
        c2->SaveAs("int.C");
   
     TCanvas *c3 = new TCanvas("ImAmp", "ImAmp", 800, 1200);
	c3->Divide(2,4);

        c3->cd(1);
        ImAmpM->Draw("e1");

        c3->cd(2);
        ImAmpMandt->Draw("e1");

        c3->cd(3);
        ImAmpBeamE->Draw("e1");

	c3->cd(4);
	ImAmpCosTheta->Draw("e1");

        c3->cd(5);
        ImAmpPhi->Draw("e1");

        c3->cd(6);
        ImAmpCosThetaH->Draw("e1");

        c3->cd(7);
        ImAmpPhiH->Draw("e1");

        c3->cd(8);
        ImAmpProdAng->Draw("e1");

	c3->SaveAs("int.pdf)");
        c3->SaveAs("int.C)");
*/      
  return 0;
}



