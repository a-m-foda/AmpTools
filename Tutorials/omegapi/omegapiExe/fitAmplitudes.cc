
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>

#include "omegapiDataIO/ROOTDataReader.h"
#include "omegapiDataIO/ROOTDataReaderBootstrap.h"
#include "omegapiDataIO/ROOTDataReaderWithTCut.h"
#include "omegapiDataIO/ROOTDataReaderTEM.h"
#include "omegapiAmp/omegapiAngAmp.h"
#include "omegapiAmp/omegapi_AmpPol.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){
	
  // set default parameters
  
  bool useMinos = false;

  string configfile;
  string seedfile;
  
  // parse command line
  
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-c"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-s"){
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  seedfile = argv[++i]; }
    if (arg == "-n") useMinos = true;
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "   -n \t\t\t\t\t use MINOS instead of MIGRAD" << endl;
      cout << "   -c <file>\t\t\t\t config file" << endl;
      cout << "   -s <output file>\t\t\t for seeding next fit based on this fit (optional)" << endl;
      exit(1);}
  }
  
  if (configfile.size() == 0){
    cout << "No config file specified" << endl;
    exit(1);
  }

  ConfigFileParser parser(configfile);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  AmpToolsInterface::registerAmplitude( omegapiAngAmp() );
  AmpToolsInterface::registerAmplitude( omegapi_AmpPol() );
  
  AmpToolsInterface::registerDataReader( ROOTDataReader() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderBootstrap() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderWithTCut() );
  AmpToolsInterface::registerDataReader( ROOTDataReaderTEM() );
  
  AmpToolsInterface ati( cfgInfo );
  
  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ati.likelihood() << endl;
  
  MinuitMinimizationManager* fitManager = ati.minuitMinimizationManager();
  
  if( useMinos ){
    
    fitManager->minosMinimization();
  }
  else{
    
    fitManager->migradMinimization();
  }
  
  bool fitFailed =
    ( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 );
  
  if( fitFailed ){
    cout << "ERROR: fit failed use results with caution..." << endl;
  }
  
  cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ati.likelihood() << endl;
  
  ati.finalizeFit();
  
  if( seedfile.size() != 0 && !fitFailed ){
    
    ati.fitResults()->writeSeed( seedfile );
  }
  
  return 0;
}


