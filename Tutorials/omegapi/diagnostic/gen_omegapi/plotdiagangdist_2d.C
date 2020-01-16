#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLine.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TEllipse.h>

void plotdiagangdist_2d() {
  char STRING[256];
    sprintf(STRING,"gen_omegapi_diagnostic.root");


 TFile *file = TFile::Open(STRING);

// set the paper & margin sizes
gStyle->SetPaperSize(20,26);
gStyle->SetPadTopMargin(0.12);
gStyle->SetPadRightMargin(0.05);
gStyle->SetPadBottomMargin(0.16);
gStyle->SetPadLeftMargin(0.12);

// use large Times-Roman fonts
//gStyle->SetTextFont(132);
gStyle->SetTextSize(0.08);
//gStyle->SetLabelFont(132,"x");
//gStyle->SetLabelFont(132,"y");
//gStyle->SetLabelFont(132,"z");
gStyle->SetLabelSize(0.05,"x");
gStyle->SetTitleSize(0.07,"x");
gStyle->SetLabelSize(0.05,"y");
gStyle->SetTitleSize(0.07,"y");
gStyle->SetLabelSize(0.05,"z");
gStyle->SetTitleSize(0.07,"z");
   gStyle->SetTitleSize(0.08,"t");
//   gStyle->SetTitleSize(0.065,"xy");
   gStyle->SetTitleOffset(0.8,"y");
   gROOT->ForceStyle();

   TH1F *M = (TH1F*)file->Get("M");
   TH1F *M_isobar = (TH1F*)file->Get("M_isobar");
   TH2F *M_Dalitzxy = (TH2F*)file->Get("M_dalitzxy");

   TH2F *M_CosTheta = (TH2F*)file->Get("M_CosTheta");
   TH2F *M_Phi = (TH2F*)file->Get("M_Phi");
   TH2F *M_CosThetaH = (TH2F*)file->Get("M_CosThetaH");
   TH2F *M_PhiH = (TH2F*)file->Get("M_PhiH");
   TH2F *M_Phi_Prod = (TH2F*)file->Get("M_Phi_Prod");

  ////////////////////////////////////////////////////////////////////////////////////////////////////

      TCanvas *q5= new TCanvas("Diagnostic Hist.","Diagnostic Hist.");
      q5->Divide(2,4);

   q5->cd(1);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M->GetXaxis()->SetRangeUser(0.8, 2.);
     M->SetTitle("M(#pi^+ #pi^- #pi^0 #pi^0) Gev/c^2; Events");
     M->SetStats(false);
     M->Draw("E1");


     /////////////////////////////////////////////////////////////////     
     q5->cd(2);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M_isobar->GetXaxis()->SetRangeUser(0.6, 1.0);
     M_isobar->SetTitle("M(#pi^+ #pi^- #pi^0) Gev/c^2; Events");
     M_isobar->SetStats(false);
     M_isobar->Draw("E1");


     /////////////////////////////////////////////////////////////////     
     q5->cd(3);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M_Dalitzxy->GetXaxis()->SetRangeUser(-1.55, 1.55);
     M_Dalitzxy->GetYaxis()->SetRangeUser(-2.55, 0.05);
     M_Dalitzxy->SetTitle("Dalitz X; Dalitz Y");
     M_Dalitzxy->SetStats(false);
     M_Dalitzxy->Draw("colz");


     /////////////////////////////////////////////////////////////////     
     q5->cd(4);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M_Phi_Prod->GetXaxis()->SetRangeUser(0.6, 2.0);
     M_Phi_Prod->SetStats(false);
     M_Phi_Prod->Draw("colz");


     /////////////////////////////////////////////////////////////////     
     q5->cd(5);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M_CosTheta->GetXaxis()->SetRangeUser(0.6, 2.0);
     M_CosTheta->SetStats(false);
     M_CosTheta->Draw("colz");


     /////////////////////////////////////////////////////////////////     
     q5->cd(6);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M_Phi->GetXaxis()->SetRangeUser(0.6, 2.0);
     M_Phi->SetStats(false);
     M_Phi->Draw("colz");


     /////////////////////////////////////////////////////////////////     
     q5->cd(7);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M_CosThetaH->GetXaxis()->SetRangeUser(0.6, 2.0);
     M_CosThetaH->SetStats(false);
     M_CosThetaH->Draw("colz");


     /////////////////////////////////////////////////////////////////     
     q5->cd(8);
     gStyle->SetFuncWidth(2);
     gStyle->SetOptFit();

     M_PhiH->GetXaxis()->SetRangeUser(0.6, 2.0);
     M_PhiH->SetStats(false);
     M_PhiH->Draw("colz");

     q5->Print("gen_amp_angulardist.C");
     q5->Print("gen_amp_angulardist.pdf");

}
