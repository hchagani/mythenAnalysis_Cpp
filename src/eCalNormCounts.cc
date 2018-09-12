// Plot Normalised Counts vs. Energy for no, Mn, Cu & Zr targets
// Plot Standard, High Gain & Fast settings for each target
// Counts normalised to unity at inflection points
//
// Requires energy calibration and fit parameter files from eCalFit
// Assumes data at each setting taken in following order: shutter closed; no target; Mn target; Cu target; Zr target
// If shutter closed data does not exist use <firstFileIndex> - 1
//
// Hassan Chagani 2016-10-21 12:53
//
// gcc v4.4.7
// root v5.99/01
//
// To compile $: g++ -o ../bin/eCalNormCounts eCalNormCounts.cc `root-config --cflags --glibs`
// To run $: ../bin/eCalNormCounts <firstStandardFileIndex> <firstHighGainFileIndex> <firstFastFileIndex> <moduleID>

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

bool CheckFile(const char * fileName);
inline double GetMedian(const int sortCounts[], const int & nStrips);
void SetHistProperties(TH2I * hist);

#ifndef __CINT__
int main(int argc, char *argv[]) {

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // Provide guidance and exit if incorrect number of arguments
  if(argc != 5) {

    cout << "Insufficient arguments!" << endl << endl
         << "Usage:" << endl
         << "./eCalNormCounts <firstStandardFileIndex> <firstHighGainFileIndex> <firstFastFileIndex> <moduleID>" << endl;

  return 1;

  }

  // Assign command line arguments to constants and variables
  const int nSettings = 3;   // Standard, high gain & fast
  int fileIndices[nSettings];
  for(int i = 0 ; i < nSettings ; i++) fileIndices[i] = atoi(argv[i+1]) + 1;
  char moduleID[10];
  strcpy(moduleID,argv[4]);

  // Fit parameter variables
  const int nFits = 5;   // no, Mn, Cu, Zr & (Mn - no)
  double maximum[nSettings][nFits],steepness[nSettings][nFits],midPoint[nSettings][nFits];
  double offsetMnCor[nSettings];   // Offset from (Mn target - no target) fit

  // Calibration variables
  double offsetCal[nSettings],gainCal[nSettings];

  // PART 1 : FIT PARAMETER AND CALIBRATION FILES HANDLING

  // Handling of fit parameter and calibration file names
  // Check to see if file exists
  // If so read in data
  // Otherwise exit
  for(int i = 0 ; i < nSettings ; i++) {

    char settings[20];
    if(i == 0) strcpy(settings,"standard");
    else if(i == 1) strcpy(settings,"highgain");
    else if(i == 2) strcpy(settings,"fast");

    // Fit parameter file first
    char fitParameterFileName[100];
    sprintf(fitParameterFileName,"/localhome/local/mythen/analysis/main/settings/%s_%s_fitParameters.txt",moduleID,settings);
    if(CheckFile(fitParameterFileName)) {

      ifstream fitParameterFile;
      fitParameterFile.open(fitParameterFileName);

      // Read data from fit parameter file line-by-line
      for(int j = 0 ; j < nFits ; j++) {

        double skipRange;   // Throwaway variable to read in fit limits
        fitParameterFile >> skipRange >> skipRange >> maximum[i][j] >> steepness[i][j] >> midPoint[i][j];
        if(j == 4) fitParameterFile >> offsetMnCor[i];

      }

      fitParameterFile.close();   // Tidy up

    } else {

      cout << fitParameterFileName << " does not exist. Please create one with:" << endl
           << "./eCalFit <firstFileIndex> <settings> <moduleID>" << endl;

      return 1;

    }

    // Calibration file last
    char calFileName[100];
    sprintf(calFileName,"/localhome/local/mythen/analysis/main/settings/newcal_%s.sn%s",settings,moduleID);
    if(CheckFile(calFileName)) {

      ifstream calFile;
      calFile.open(calFileName);

      // Read data from calibration file
      calFile >> offsetCal[i] >> gainCal[i];

      calFile.close();   // Tidy up

    } else {

      cout << calFileName << " does not exist. Please create one with:" << endl
           << "./eCalFit <firstFileIndex> <settings> <moduleID> fit" << endl;

      return 1;

    }

  }

  // PART 2 : READ DATA FROM FILES

  // Arrays for data
  const int nThresh = 651;
  const int offsetThresh = 200;
  const int nStrips = 1280;
  const int nTargets = nFits - 1;   // Remove (Mn target - no target) from configurations
  double energy[nSettings][nThresh];
  double normCounts0T[nSettings][nThresh],normCountsMn[nSettings][nThresh],normCountsCu[nSettings][nThresh],normCountsZr[nSettings][nThresh];
  double normCountsMnSub0T[nSettings][nThresh];

  // Data filename string handling
  ifstream dataFile;
  char fileNamePrefix[100],fileNamePostfix[5];
  char currentFileName[100];
  strcpy(fileNamePrefix,"/localhome/local/Desktop/OutputDir/run_S");
  strcpy(fileNamePostfix,".raw");

  // Loop over all settings
  for(int i = 0 ; i < nSettings ; i++) {

    // Loop over all targets
    for(int j = 0 ; j < nTargets ; j++) {

      // Loop over all files where k is vthreshold
      for(int k = 0 ; k < nThresh ; k++) {

        strcpy(currentFileName,fileNamePrefix);
        char fileID[10];
        sprintf(fileID,"%d_%d",k+offsetThresh,j+fileIndices[i]);
        strcat(currentFileName,fileID);
        strcat(currentFileName,fileNamePostfix);

        dataFile.open(currentFileName);
        if(k%100 == 0) cout << "File: " << currentFileName << endl;

        int skipStrip;   // Throwaway variable to read in strip number
        int sortCounts[nStrips];   // Array for sorted counts

        // Loop over data in file where l is channel number
        for(int l = 0 ; l < nStrips ; l++) dataFile >> skipStrip >> sortCounts[l];

        // Fill energy array once per setting
        if(j == 0) energy[i][k] = (k + offsetThresh - offsetCal[i]) / gainCal[i];

        sort(sortCounts,sortCounts+nStrips);   // Sort counts array in ascending order

        // Inflection point occurs at half of maximum counts from S-curve fit
        if(j == 0) normCounts0T[i][k] = GetMedian(sortCounts,nStrips) / (maximum[i][j] / 2.);
        else if(j == 1) {
          normCountsMn[i][k] = GetMedian(sortCounts,nStrips) / (maximum[i][j] / 2.);
          normCountsMnSub0T[i][k] = GetMedian(sortCounts,nStrips) / (maximum[i][j+3] / 2.);
        }
        else if(j == 2) normCountsCu[i][k] = GetMedian(sortCounts,nStrips) / (maximum[i][j] / 2.);
        else if(j == 3) normCountsZr[i][k] = GetMedian(sortCounts,nStrips) / (maximum[i][j] / 2.);

        dataFile.close();   // Tidy up

      }

    }

  }

  // PART 3 : PLOT DATA

  TCanvas *c1 = new TCanvas("c1","",2970,2100);
  c1->Divide(2,2);

  // Empty histograms with axes limits for no, Mn, Cu & Zr targets normalised counts vs. energy
  TH2I *h0T = new TH2I("h0T","",10,0,12,10,0,3);
  SetHistProperties(h0T);
  char h0TTitle[100];
  sprintf(h0TTitle,"no Target : %s",moduleID);
  h0T->SetTitle(h0TTitle);
  TH2I *hMn = new TH2I("hMn","",10,0,12,10,0,7);
  SetHistProperties(hMn);
  char hMnTitle[100];
  sprintf(hMnTitle,"Mn Target : %s",moduleID);
  hMn->SetTitle(hMnTitle);
  TH2I *hCu = new TH2I("hCu","",10,0,12,10,0,3);
  SetHistProperties(hCu);
  char hCuTitle[100];
  sprintf(hCuTitle,"Cu Target : %s",moduleID);
  hCu->SetTitle(hCuTitle);
  TH2I *hZr = new TH2I("hZr","",10,0,21,10,0,5);
  SetHistProperties(hZr);
  char hZrTitle[100];
  sprintf(hZrTitle,"Zr Target : %s",moduleID);
  hZr->SetTitle(hZrTitle);

  // Declare and fill graphs with normalised counts vs. energy for no, Mn, Cu & Zr targets
  TGraph *gr[nSettings][nTargets];
  TGraph *grMnSub0T[nSettings];   // Use inflection point from (Mn - no) target fit

  for(int i = 0 ; i < nSettings ; i++) {

    gr[i][0] = new TGraph(nThresh,&energy[i][0],&normCounts0T[i][0]);
    gr[i][1] = new TGraph(nThresh,&energy[i][0],&normCountsMn[i][0]);
    gr[i][2] = new TGraph(nThresh,&energy[i][0],&normCountsCu[i][0]);
    gr[i][3] = new TGraph(nThresh,&energy[i][0],&normCountsZr[i][0]);
    grMnSub0T[i] = new TGraph(nThresh,&energy[i][0],&normCountsMnSub0T[i][0]);

  }

  c1->cd(1)->SetGrid();
  h0T->Draw();
  for(int i = 0 ; i < nSettings ; i++) {   // no target
    gr[i][0]->SetMarkerStyle(5);
    gr[i][0]->SetMarkerColor(1+i);
    gr[i][0]->Draw("p");
  }
  c1->cd(2)->SetGrid();
  hMn->Draw();
  for(int i = 0 ; i < nSettings ; i++) {   // Mn target
    gr[i][1]->SetMarkerStyle(5);
    gr[i][1]->SetMarkerColor(4+i);
    gr[i][1]->Draw("p");
    grMnSub0T[i]->SetMarkerStyle(5);
    grMnSub0T[i]->SetMarkerColor(1+i);
    grMnSub0T[i]->Draw("p");
  }
  c1->cd(3)->SetGrid();
  hCu->Draw();
  for(int i = 0 ; i < nSettings ; i++) {   // Cu target
    gr[i][2]->SetMarkerStyle(5);
    gr[i][2]->SetMarkerColor(1+i);
    gr[i][2]->Draw("p");
  }
  c1->cd(4)->SetGrid();
  hZr->Draw();
  for(int i = 0 ; i < nSettings ; i++) {   // Zr target
    gr[i][3]->SetMarkerStyle(5);
    gr[i][3]->SetMarkerColor(1+i);
    gr[i][3]->Draw("p");
  }

  // Legend for no, Cu & Zr targets graphs
  TLegend *leg = new TLegend(0.15,0.15,0.35,0.35);
  leg->SetTextSize(0.035);
  leg->SetHeader("Settings");
  leg->SetFillColor(10);
  leg->AddEntry(gr[1][0],"High Gain","p");
  leg->AddEntry(gr[0][0],"Standard","p");
  leg->AddEntry(gr[2][0],"Fast","p");

  // Legend for Mn target graphs
  TLegend *legMnSub0T = new TLegend(0.45,0.65,0.9,0.9);
  legMnSub0T->SetTextSize(0.035);
  legMnSub0T->SetFillColor(10);
  legMnSub0T->SetHeader("Inflection point from:");
  legMnSub0T->SetNColumns(2);
  legMnSub0T->AddEntry((TObject*)0,"Mn Target only","");
  legMnSub0T->AddEntry((TObject*)0,"(Mn - no) Target","");
  legMnSub0T->AddEntry(gr[1][1],"High Gain","p");
  legMnSub0T->AddEntry(grMnSub0T[1],"High Gain","p");
  legMnSub0T->AddEntry(gr[0][1],"Standard","p");
  legMnSub0T->AddEntry(grMnSub0T[0],"Standard","p");
  legMnSub0T->AddEntry(gr[2][1],"Fast","p");
  legMnSub0T->AddEntry(grMnSub0T[2],"Fast","p");

  // Draw legends
  c1->cd(1);
  leg->Draw();
  c1->cd(2);
  legMnSub0T->Draw();
  c1->cd(3);
  leg->Draw();
  c1->cd(4);
  leg->Draw();

  // Print canvas to file in pdf format
  char outputFileName[100];
  sprintf(outputFileName,"/localhome/local/mythen/analysis/main/plots/%s_eCalNormCounts.pdf",moduleID);
  c1->Print(outputFileName);

  // Tidy up
  delete c1;
  c1 = 0;
  delete h0T;
  h0T = 0;
  delete hMn;
  hMn = 0;
  delete hCu;
  hCu = 0;
  delete hZr;
  hZr = 0;
  for(int i = 0 ; i < nSettings ; i++) {
    for(int j = 0 ; j < nTargets ; j++) {
      delete gr[i][j];
      gr[i][j] = 0;
    }
    delete grMnSub0T[i];
    grMnSub0T[i] = 0;
  }
  delete leg;
  leg = 0;
  delete legMnSub0T;
  legMnSub0T = 0;

  return 0;

}
#endif

// Check to determine if file exists by opening it
bool CheckFile(const char * fileName) {

  ifstream fileCheck;
  fileCheck.open(fileName);

  return fileCheck;

}

// Get median from sorted array of counts
inline double GetMedian(const int sortCounts[], const int & nStrips) {

  return ((double)sortCounts[nStrips/2 - 1] + (double)sortCounts[nStrips/2]) / 2.;

}

// Set empty histogram properties
void SetHistProperties(TH2I * hist) {

  hist->GetXaxis()->SetTitle("Energy (keV)");
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitle("Normalised Counts");
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);

}
