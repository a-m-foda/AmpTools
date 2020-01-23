#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

//#include "omegapiDataIO/OmegaPiPlotGenerator.h"
#include "omegapiDataIO/ROOTDataReader.h"
#include "omegapiDataIO/ROOTDataReaderBootstrap.h"
#include "omegapiDataIO/ROOTDataReaderWithTCut.h"
#include "omegapiDataIO/ROOTDataReaderTEM.h"
#include "omegapiAmp/omegapiAngAmp.h"
#include "omegapiAmp/omegapi_AmpPol.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

// typedef OmegaPiPlotGenerator PlotGen;

void atiSetup(){
  
  AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
  AmpToolsInterface::registerAmplitude( omegapi_AmpPol() );
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
}

using namespace std;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Viewing Results Using AmpPlotter and writing root histograms *** " << endl << endl;

  if (argc < 2){
    cout << "Usage:" << endl << endl;
    cout << "\tomegapi_plotter <results file name> -o <output file name>" << endl << endl;
    return 0;
  }

  const int NParam = 28;
  const int NumTBins = 10;
  const int NumEBins = 10;
  const int NumMBins = 10;
  
  bool showGui = false;
  string outName = "omegapi_plot.root";

  for (int i = 1; i < argc; i++){

    string arg(argv[i]);

    if (arg == "-g"){
      showGui = true;
    }
    if (arg == "-o"){
      outName = argv[++i];
    }

/*    if (arg == "-t"){  
                    if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                    else  NumTBins = atoi( argv[++i] ); }
                    
    if (arg == "-e"){  
                    if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                    else  NumEBins = atoi( argv[++i] ); }
                    
    if (arg == "-m"){  
                    if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
                    else  NumMBins = atoi( argv[++i] ); }
*/                  
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -o <file>\t output file path" << endl;
      cout << "\t -g <file>\t show GUI" << endl;
      exit(1);
    }
  }


    // ************************
    // parse the command line parameters
    // ************************

  cout << "Output file name    = " << outName << endl << endl;
  TFile* plotfile = new TFile( outName.c_str(), "recreate");
  TH1::AddDirectory(kFALSE);

  ofstream outfile;
  outfile.open( "omegapi_fitPars.txt" );

  char STRING[256];

  double zero[NumMBins] = {0.0};
    
  double m_ll = 0.8;
  double m_ul = 2.0;
    
  double t_ll = 0.1;
  double t_ul = 1.0;
    
  double e_ll = 6.0;
  double e_ul = 11.6;
    
  double MBins[NumMBins];

    for(int k=0; k<NumMBins; k++){

        MBins[k] = m_ll + (k*((m_ul-m_ll)/NumMBins));
    
    }//Mass bin Loop

  double param[NParam][NumTBins][NumEBins][NumMBins];
  double parame[NParam][NumTBins][NumEBins][NumMBins];

// parameters to check
  vector< string > pars;

  pars.push_back("hel_c_0_m_1");
  pars.push_back("hel_c_1_m_1");
  pars.push_back("hel_c_1_p_0");
  pars.push_back("hel_c_1_p_2");
  pars.push_back("hel_c_2_m_1");
  pars.push_back("hel_c_2_p_2");
  pars.push_back("vertex_0_m_m_0");
  pars.push_back("vertex_0_p_m_0");
  pars.push_back("vertex_1_m_m_m1");
  pars.push_back("vertex_1_m_m_0");
  pars.push_back("vertex_1_m_m_p1");
  pars.push_back("vertex_1_p_m_m1");
  pars.push_back("vertex_1_p_m_0");
  pars.push_back("vertex_1_p_m_p1");
  pars.push_back("vertex_2_m_m_m2");
  pars.push_back("vertex_2_m_m_m1");
  pars.push_back("vertex_2_m_m_0");
  pars.push_back("vertex_2_m_m_p1");
  pars.push_back("vertex_2_m_m_p2");
  pars.push_back("vertex_2_p_m_m2");
  pars.push_back("vertex_2_p_m_m1");
  pars.push_back("vertex_2_p_m_0");
  pars.push_back("vertex_2_p_m_p1");
  pars.push_back("vertex_2_p_m_p2");
  pars.push_back("dalitz_alpha");
  pars.push_back("dalitz_beta");
  pars.push_back("dalitz_gamma");
  pars.push_back("dalitz_delta");

  TMultiGraph *mg[NParam][NumTBins];
  TGraphErrors *gr[NParam][NumTBins][NumEBins];

  for(int i=0; i < NumTBins; i++)
  {//loop over Mand T bin
      
        for(int j=0; j < NumEBins; j++)
        {//loop over Beam E bin
            
            for(int k=0; k < NumMBins; k++)
            {//loop over Inv Mass bin
        
            sprintf(STRING,"t_%01i_e_%01i_m_%01i/t_%01i_e_%01i_m_%01i.fit", i, j, k, i, j, k);
	    std::string resultsName = STRING;
            cout << "Fit results file name    = " << resultsName << endl;

                // ************************
                // load the results and display the configuration info
                // ************************

            //cout << "Loading Fit results" << endl;
            FitResults results( resultsName );
            if( !results.valid() ){
                
                cout << "Invalid fit results in file:  " << resultsName << endl;
                continue;//exit( 1 );
            }
            //cout << "Fit results loaded" << endl;
                // ************************
                // set up the plot generator
                // ************************

            atiSetup();

            //cout << " Initialized ati" << endl;

	          // ************************
                // retrieve parameters for plotting and asymmetry
                // ************************

            //cout << "Checking Parameters" << endl;
            // file for writing parameters (later switch to putting in ROOT file)

            for(unsigned int n = 0; n<pars.size(); n++) 
            {
                double parValue = results.parValue( pars[i] );
                double parError = results.parError( pars[i] );
                outfile << parValue << "\t" << parError << "\t";
                
            //     param[id][mandt][beame][mass];
                 param[n][i][j][k] = parValue;
                 parame[n][i][j][k] = parError;           

            }

            // covariance matrix
            vector< vector< double > > covMatrix;
            covMatrix = results.errorMatrix();

            outfile << endl;

            resultsName.clear();

            }//Inv Mass loop
            
            for( int  n= 0; n < NParam; n++ ) 
            {
                
            gr[n][i][j] = new TGraphErrors( NumMBins, MBins, param[n][i][j], zero, parame[n][i][j] );
            gr[n][i][j]->SetMarkerColor(j+1);
            gr[n][i][j]->SetLineColor(j+1);
            gr[n][i][j]->SetMarkerStyle(2);
	    
	    mg[n][i] = new TMultiGraph();
            mg[n][i]->Add( gr[n][i][j] ); 

            }//parameter
         
           }//Beam E loop
        
    TCanvas* can = new TCanvas( "can", "Amplitude Analysis Plots", 800, 1400);
    can->Divide( 3, 5 );
    
    can->cd( 1 );
    mg[0][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; C^{0^{-}}_{1}");
    mg[0][i]->Draw( "ALP" );
    mg[0][i]->Write();

    can->cd( 2 );
    mg[1][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; C^{0^{-}}_{1}");
    mg[1][i]->Draw( "ALP" );
    mg[1][i]->Write();
    
    can->cd( 3 );
    mg[2][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; C^{0^{-}}_{1}");
    mg[2][i]->Draw( "ALP" );
    mg[2][i]->Write();

    can->cd( 4 );
    mg[3][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; C^{0^{-}}_{1}");
    mg[3][i]->Draw( "ALP" );
    mg[3][i]->Write();
    
    can->cd( 5 );
    mg[4][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; C^{0^{-}}_{1}");
    mg[4][i]->Draw( "ALP" );
    mg[4][i]->Write();
    
    can->cd( 6 );
    mg[5][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; C^{0^{-}}_{1}");
    mg[5][i]->Draw( "ALP" );
    mg[5][i]->Write();
    
    can->cd( 7 );
    mg[6][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[6][i]->Draw( "ALP" );
    mg[6][i]->Write();
    
    can->cd( 8 );
    mg[7][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[7][i]->Draw( "ALP" );
    mg[7][i]->Write();
    
    can->cd( 9 );
    mg[8][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[8][i]->Draw( "ALP" );
    mg[8][i]->Write();
    
    can->cd( 10 );
    mg[9][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[9][i]->Draw( "ALP" );
    mg[9][i]->Write();
    
    can->cd( 11 );
    mg[10][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[10][i]->Draw( "ALP" );
    mg[10][i]->Write();
    
    can->cd( 12 );
    mg[11][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[11][i]->Draw( "ALP" );
    mg[11][i]->Write();
    
    can->cd( 13 );
    mg[12][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[12][i]->Draw( "ALP" );
    mg[12][i]->Write();
    
    can->cd( 14 );
    mg[13][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[13][i]->Draw( "ALP" );
    mg[13][i]->Write();

    sprintf(STRING,"fitresults_t_%01i.C(", i);
    can->SaveAs(STRING);
    sprintf(STRING,"fitresults_t_%01i.pdf(", i);
    can->SaveAs(STRING);

    TCanvas* can2 = new TCanvas( "can2", "Amplitude Analysis Plots cont.", 800, 1400);
    can2->Divide( 3, 5 );

    
    can2->cd( 1 );
    mg[14][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[14][i]->Draw( "ALP" );
    mg[14][i]->Write();
    
    can2->cd( 2 );
    mg[15][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[15][i]->Draw( "ALP" );
    mg[15][i]->Write();
    
    can2->cd( 3 );
    mg[16][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[16][i]->Draw( "ALP" );
    mg[16][i]->Write();
    
    can2->cd( 4 );
    mg[17][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[17][i]->Draw( "ALP" );
    mg[17][i]->Write();
    
    can2->cd( 5 );
    mg[18][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[18][i]->Draw( "ALP" );
    mg[18][i]->Write();
    
    can2->cd( 6 );
    mg[19][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[19][i]->Draw( "ALP" );
    mg[19][i]->Write();
    
    can2->cd( 7 );
    mg[20][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[20][i]->Draw( "ALP" );
    mg[20][i]->Write();
    
    can2->cd( 8 );
    mg[21][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[21][i]->Draw( "ALP" );
    mg[21][i]->Write();
    
    can2->cd( 9 );
    mg[22][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[22][i]->Draw( "ALP" );
    mg[22][i]->Write();
    
    can2->cd( 10 );
    mg[23][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; V^{0^{-}}_{-1,0}");
    mg[23][i]->Draw( "ALP" );
    mg[23][i]->Write();
    
    can2->cd( 11 );
    mg[24][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; #alpha");
    mg[24][i]->Draw( "ALP" );
    mg[24][i]->Write();
    
    can2->cd( 12 );
    mg[25][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; #beta");
    mg[25][i]->Draw( "ALP" );
    mg[25][i]->Write();
    
    can2->cd( 13 );
    mg[26][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; #gamma");
    mg[26][i]->Draw( "ALP" );
    mg[26][i]->Write();
    
    can2->cd( 14 );
    mg[27][i]->SetTitle("Mand t = ; M(#omega #pi^{0}) GeV/c^{2}; #delta");
    mg[27][i]->Draw( "ALP" );
    mg[27][i]->Write();

    sprintf(STRING,"fitresults_t_%01i.C)", i);
    can2->SaveAs(STRING);
    sprintf(STRING,"fitresults_t_%01i.pdf)", i);
    can2->SaveAs(STRING);

  }//Mand T loop


        
        
  
  plotfile->Close();

    // ************************
    // start the GUI
    // ************************

  if(showGui) {

	  cout << ">> Plot generator ready, starting GUI..." << endl;
	  
	  int dummy_argc = 0;
	  char* dummy_argv[] = {};  
	  TApplication app( "app", &dummy_argc, dummy_argv );
	  
	  gStyle->SetFillColor(10);
	  gStyle->SetCanvasColor(10);
	  gStyle->SetPadColor(10);
	  gStyle->SetFillStyle(1001);
	  gStyle->SetPalette(1);
	  gStyle->SetFrameFillColor(10);
	  gStyle->SetFrameFillStyle(1001);
   
     cout << " Initialized App " << endl;     
// 	  PlotFactory factory( plotGen );	
     cout << " Created Plot Factory " << endl;
	  //PlotterMainWindow mainFrame( gClient->GetRoot(), factory );
     cout << " Main frame created " << endl;
	  
	  app.Run();
     cout << " App running" << endl;
  }
    
  return 0;

}


