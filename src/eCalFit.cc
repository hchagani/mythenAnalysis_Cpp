// Energy Calibration
// Canvas 1 : Plot Median Counts vs. vthreshold for no, Mn, Cu & Zr targets
//            Fit S-curve to distributions
// Canvas 2 : Plot Median Counts vs. vthreshold for all and (Mn - no) targets
//            Plot inflection points from S-curve fits to Mn, Cu & Zr targets
//            Linear fit to inflection points
//            Mn only (bottom left) ; (Mn - no) (bottom right)
// Canvas 3 : Plot Mean Counts vs. vthreshold for Cu & Zr targets (top)
//            Plot Standard Deviation from Mean vs. vthreshold (bottom)
//
// Does not produce inflection point plots for energy calibration if fit argument is not given
// Allows user to create fit parameter file with intial guesses if fit argument is not given
// Performs fit using values from fit parameter file if fit argument is given
// Overwrites values in fit parameter file if fit argument is given
// Produces calibration file
//
// Assumes data taken in following order: shutter closed; no target; Mn target; Cu target; Zr target
//
// Hassan Chagani 2016-10-24 09:55
//
// gcc v4.4.7
// root v.5.99/01
//
// To compile:
// $: g++ -o ../bin/eCalFit eCalFit.cc `root-config --cflags -- glibs`
// To run:
// $: ../bin/eCalFit <firstFileIndex> <settings> <moduleID>
// or
// $: ../bin/eCalFit <firstFileIndex> <settings> <moduleID> fit

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
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TMath.h"
#include "TLatex.h"

using namespace std;

bool CheckFile(const char * fileName);
inline int GetMedian(const int sortCounts[], const int & nStrips);
double GetMean(const int sortCounts[], const int & nStrips);
double GetStandDev(const int sortCounts[], const double & meanCounts, const int & nStrips);
void SetHistProperties(TH2I * hist);
void SetHistStandDevProperties(TH2I * hist);
void SetGraphProperties(TGraphErrors * graph);
void SetFitParameterOptions(const char * fileName);
void SetFitTextLabel(TPaveText * pText, const TF1 * fit);
void Set4ParFitTextLabel(TPaveText * pText, const TF1 * fit);
inline void SetChiSqString(char * eqString, const TF1 * fit);
inline void SetEqString(char * eqString, const TF1 * fit);
inline void Set4ParEqString(char * eqString, const TF1 * fit);
inline void SetCalEqString(char * eqString, const TF1 * fit);
void SetNoTargetEnergy(char * energyString, const TF1 * calFit, const TF1 * sCurveFit);

#ifndef __CINT__
int main(int argc, char *argv[]) {

  // Provide guidance and exit if incorrect number of arguments
  if(argc != 4 && argc != 5) {

    cout << "Insufficient arguments!" << endl << endl
         << "Usage:" << endl
         << "./eCalFit <firstFileIndex> <settings> <moduleID>" << endl
         << "./eCalFit <firstFileIndex> <settings> <moduleID> fit" << endl;

    return 1;

  }

  bool fit = false;   // Fit argument does not exist
  if(argc == 5) {   // Provide guidance and exit if fit argument is incorrectly specified

    if(strcmp(argv[4],"fit")) {

      cout << "If you would like to perform fits, issue the following arguments:" << endl
           << "./eCalFit <firstFileIndex> <settings> <moduleID> fit" << endl;

      return 1;

    } else fit = true;   // Otherwise fit argument exists

  }

  // Assign command line arguments to constants and variables
  const int firstFileIndex = atoi(argv[1]);
  char settings[20],settingsTitle[20];
  if(strcmp(argv[2],"standard") == 0 || strcmp(argv[2],"Standard") == 0) {
    strcpy(settings,"standard");
    strcpy(settingsTitle,"Standard");
  } else if(strcmp(argv[2],"highgain") == 0 || strcmp(argv[2],"Highgain") == 0 || strcmp(argv[2],"HighGain") == 0) {
    strcpy(settings,"highgain");
    strcpy(settingsTitle,"High Gain");
  } else if(strcmp(argv[2],"fast") == 0 || strcmp(argv[2],"Fast") == 0) {
    strcpy(settings,"fast");
    strcpy(settingsTitle,"Fast");
  } else {
    cout << "Settings argument <" << argv[2] << "> incorrect. Exiting..." << endl;
    return 1;
  }
  char moduleID[10];
  strcpy(moduleID,argv[3]);

  // Handling of fit parameter file name and declaration of variables to store data from it
  char fitParameterFileName[100];
  sprintf(fitParameterFileName,"/localhome/local/mythen/analysis/main/settings/%s_%s_fitParameters.txt",moduleID,settings);
  bool fitParameterFileExists = false;
  const int noOfFits = 5;
  double minRange[noOfFits],maxRange[noOfFits],maximum[noOfFits],midPoint[noOfFits];
  double offsetMnSubNo;
  double skipSteepness;   // Steepness unnecessary for initial fit

  // If fit argument was given by user check to see if fit parameter file exists
  if(fit == true) {

    fitParameterFileExists = CheckFile(fitParameterFileName);

    if(fitParameterFileExists == false) {   // Provide guidance and exit if fit parameter file does not exist

      cout << "Fit parameter file does not exist. Please create one with:" << endl
           << "./eCalFit <firstFileIndex> <settings> <moduleID>" << endl;

      return 1;

    } else {   // Otherwise read in data from file

      ifstream fitParameterFile;
      fitParameterFile.open(fitParameterFileName);

      for(int i = 0 ; i < noOfFits ; i++) {

        if(i == noOfFits - 1) fitParameterFile >> minRange[i] >> maxRange[i] >> maximum[i] >> skipSteepness >> midPoint[i] >> offsetMnSubNo;
        else fitParameterFile >> minRange[i] >> maxRange[i] >> maximum[i] >> skipSteepness >> midPoint[i];

      }

      fitParameterFile.close();   // Tidy up

    }

  }

  // Ask user to check that command line arguments are correct
  // Exit if they are not
  string correctSettings;
  cout << "Settings: " << endl
       << firstFileIndex << "\t\tFirst File Index" << endl
       << settings << "\tSettings" << endl
       << moduleID << "\tModule ID" << endl << endl
       << "Is this correct? [Y/n] ";
  getline(cin,correctSettings);

  if(correctSettings == "N" || correctSettings == "n" || correctSettings == "No" || correctSettings == "NO" || correctSettings == "no") {

    cout << "Exiting..." << endl;

    return 1;

  }

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  // PART 1 : READ DATA FROM FILES

  // Filename string declarations and initial handling
  ifstream dataFile;

  char fileNamePrefix[100],fileNamePostfix[5];
  char currentFileName[100];
  strcpy(fileNamePrefix,"/localhome/local/Desktop/OutputDir/run_S");
  strcpy(fileNamePostfix,".raw");

  // Arrays for data
  const int nThresh = 651;
  const int offsetThresh = 200;
  const int nStrips = 1280;
  const int nConfig = 5;
  int vthreshold[nThresh];
  int medianCounts[nConfig][nThresh];
  double meanCounts[nConfig-3][nThresh],meanCountsError[nConfig-3][nThresh];
  double standDevMeanCounts[nConfig-3][nThresh];
  int countsNoTarget[nThresh][nStrips];
  int sortCountsMn[nStrips];
  int medianCountsMn[nThresh];

  // Loop over all files where i is configuration
  // i = 0 : shutter closed
  // i = 1 : shutter open, no target
  // i = 2 : shutter open, Mn target
  // i = 3 : shutter open, Cu target
  // i = 4 : shutter open, Zr target
  for(int i = 0 ; i < nConfig ; i++) {

    // Loop over all files where j is vthreshold
    for(int j = 0 ; j < nThresh ; j++) {

      strcpy(currentFileName,fileNamePrefix);
      char fileID[10];
      sprintf(fileID,"%d_%d",j+offsetThresh,i+firstFileIndex);
      strcat(currentFileName,fileID);
      strcat(currentFileName,fileNamePostfix);

      dataFile.open(currentFileName);
      if(j%10 == 0) cout << "File: " << currentFileName << endl;

      int skipStrip;   // Throwaway variable to read in strip number
      int sortCounts[nStrips];   // Array for sorted counts

      // Loop over data in file where k is channel number
      for(int k = 0 ; k < nStrips; k++) {

        dataFile >> skipStrip >> sortCounts[k];

        if(i == 1) countsNoTarget[j][k] = sortCounts[k];   // Shutter open, no target
        else if(i == 2) sortCountsMn[k] = sortCounts[k] - countsNoTarget[j][k];   // Shutter open, Mn target

      }

      if(i == 0) vthreshold[j] = j + offsetThresh;   // Fill vthreshold array only once

      sort(sortCounts,sortCounts+nStrips);   // Sort counts array in ascending order
      medianCounts[i][j] = GetMedian(sortCounts,nStrips);   // Fill medianCounts array

      if(i == 2) {   // Subtract no target counts from Mn target counts and fill medianCounts array

        sort(sortCountsMn,sortCountsMn+nStrips);
        medianCountsMn[j] = GetMedian(sortCountsMn,nStrips);

      }

      if(i == 3 || i == 4) {   // Calculate mean, standard deviation and error for Cu & Zr targets

        meanCounts[i-3][j] = GetMean(sortCounts,nStrips);
        standDevMeanCounts[i-3][j] = GetStandDev(sortCounts,meanCounts[i-3][j],nStrips);
        meanCountsError[i-3][j] = meanCounts[i-3][j] / TMath::Sqrt((double)nStrips);

      }

      dataFile.close();

    }

  }

  // PART 2 : PLOT DATA

  TCanvas *c1 = new TCanvas("c1","",2970,2100);
  c1->Divide(2,2);
  TCanvas *c2 = new TCanvas("c2","",2970,2100);
  c2->Divide(2,2);
  TCanvas *c3 = new TCanvas("c3","",2970,2100);
  c3->Divide(2,2);

  // Declare and assign axes limits for no, Mn, Cu & Zr target median counts vs. vthreshold plots
  // Helps readability
  double h1XMin,h1XMax,h1YMax;
  const double h1YMin = 0.0;   // y-axis minimum always 0
  double h2XMin,h2XMax,h2YMax;
  const double h2YMin = h1YMin;
  double h3XMin,h3XMax,h3YMax;
  const double h3YMin = h1YMin;
  double h4XMin,h4XMax,h4YMax;
  const double h4YMin = h1YMin;
  if(strcmp(settings,"fast") == 0) {   // Axes limits for fast settings

    h1XMin = offsetThresh-0.5+350;
    h1XMax = nThresh+offsetThresh+0.5-80;
    h1YMax = *max_element(&medianCounts[1][250],&medianCounts[1][650])/2000.;
    if(strcmp(moduleID,"2e8") == 0) h1YMax = *max_element(&medianCounts[1][250],&medianCounts[1][650])/1000.;

    h2XMin = h1XMin;
    h2XMax = h1XMax;
    h2YMax = *max_element(&medianCounts[2][250],&medianCounts[2][650])/1000.;
    if(strcmp(moduleID,"2e8") == 0) h2YMax = *max_element(&medianCounts[2][250],&medianCounts[2][650])/500.;

    h3XMin = h1XMin;
    h3XMax = h1XMax;
    h3YMax = *max_element(&medianCounts[3][250],&medianCounts[3][650])/500.;
    if(strcmp(moduleID,"2e8") == 0) h3YMax = *max_element(&medianCounts[3][250],&medianCounts[3][650])/200.;

    h4XMin = h1XMin;
    h4XMax = h1XMax;
    h4YMax = *max_element(&medianCounts[4][250],&medianCounts[4][650])/1000.;
    if(strcmp(moduleID,"2e8") == 0) h4YMax = *max_element(&medianCounts[4][250],&medianCounts[4][650])/500.;

  } else {   // Axes limits for standard and high gain settings

    h1XMin = offsetThresh-0.5+250;
    h1XMax = nThresh+offsetThresh+0.5-100;
    h1YMax = *max_element(&medianCounts[1][250],&medianCounts[1][500])*1.5;

    h2XMin = h1XMin;
    h2XMax = h1XMax;
    h2YMax = *max_element(&medianCounts[2][250],&medianCounts[2][500])*3.;

    h3XMin = h1XMin;
    h3XMax = h1XMax;
    h3YMax = *max_element(&medianCounts[3][250],&medianCounts[3][500])*1.5;

    h4XMin = h1XMin;
    h4XMax = h1XMax;
    h4YMax = *max_element(&medianCounts[4][250],&medianCounts[4][500])*1.5;

  }

  // Empty histograms with axes limits for no, Mn, Cu & Zr target median counts vs. vthreshold plots
  TH2I *h1 = new TH2I("h1","",10,h1XMin,h1XMax,10,h1YMin,h1YMax);
  SetHistProperties(h1);
  char h1Title[100];
  sprintf(h1Title,"%s Settings : no target : %s",settingsTitle,moduleID);
  h1->SetTitle(h1Title);
  TH2I *h2 = new TH2I("h2","",10,h2XMin,h2XMax,10,h2YMin,h2YMax);
  SetHistProperties(h2);
  char h2Title[100];
  sprintf(h2Title,"%s Settings : Mn target : %s",settingsTitle,moduleID);
  h2->SetTitle(h2Title);
  TH2I *h3 = new TH2I("h3","",10,h3XMin,h3XMax,10,h3YMin,h3YMax);
  SetHistProperties(h3);
  char h3Title[100];
  sprintf(h3Title,"%s Settings : Cu target : %s",settingsTitle,moduleID);
  h3->SetTitle(h3Title);
  TH2I *h4 = new TH2I("h4","",10,h4XMin,h4XMax,10,h4YMin,h4YMax);
  SetHistProperties(h4);
  char h4Title[100];
  sprintf(h4Title,"%s Settings : Zr target : %s",settingsTitle,moduleID);
  h4->SetTitle(h4Title);

  // Empty histograms with axes limits for Cu & Zr target mean counts vs. vthreshold plots
  TH2I *h3Mean = new TH2I("h3Mean","",10,h3XMin,h3XMax,10,h3YMin,h3YMax);
  SetHistProperties(h3Mean);
  char h3MeanTitle[100];
  sprintf(h3MeanTitle,"%s Settings : Cu target : %s : Mean counts",settingsTitle,moduleID);
  h3Mean->SetTitle(h3MeanTitle);
  TH2I *h4Mean = new TH2I("h4Mean","",10,h4XMin,h4XMax,10,h4YMin,h4YMax);
  SetHistProperties(h4Mean);
  char h4MeanTitle[100];
  sprintf(h4MeanTitle,"%s Settings : Zr target : %s : Mean counts",settingsTitle,moduleID);
  h4Mean->SetTitle(h4MeanTitle);
  TH2I *h3StandDev = new TH2I("h3StandDev","",10,h3XMin,h3XMax,10,*min_element(&meanCounts[0][(int)h3XMin-offsetThresh],&meanCounts[0][(int)h3XMax-offsetThresh]),
                                                                  *max_element(&meanCounts[0][(int)h3XMin-offsetThresh],&meanCounts[0][(int)h3XMax-offsetThresh]));
  SetHistStandDevProperties(h3StandDev);
  char h3StandDevTitle[100];
  sprintf(h3StandDevTitle,"%s Settings : Cu target : %s : Standard deviation from mean",settingsTitle,moduleID);
  h3StandDev->SetTitle(h3StandDevTitle);
  TH2I *h4StandDev = new TH2I("h4StandDev","",10,h4XMin,h4XMax,10,*min_element(&meanCounts[1][(int)h4XMin-offsetThresh],&meanCounts[1][(int)h4XMax-offsetThresh]),
                                                                  *max_element(&meanCounts[1][(int)h4XMin-offsetThresh],&meanCounts[1][(int)h4XMax-offsetThresh]));
  SetHistStandDevProperties(h4StandDev);
  char h4StandDevTitle[100];
  sprintf(h4StandDevTitle,"%s Settings : Zr target : %s : Standard deviation from mean",settingsTitle,moduleID);
  h4StandDev->SetTitle(h4StandDevTitle);

  // Empty histogram for all target median counts vs. vthreshold plots
  // Determine y-axis maximum by finding maximum counts within range looping over each target
  double maxCounts = *max_element(&medianCounts[1][250],&medianCounts[1][500]);
  for(int i = 2 ; i < 5 ; i++) {

    if(maxCounts < *max_element(&medianCounts[i][250],&medianCounts[i][500])) maxCounts = *max_element(&medianCounts[i][250],&medianCounts[i][500]);

  }
  double hAllYMax;   // y-axis maximum is variable
  // Other axes limits are constants from first plot
  const double hAllXMin = h1XMin;
  const double hAllXMax = h1XMax;
  const double hAllYMin = h1YMin;
  if(strcmp(settings,"fast") == 0) hAllYMax = maxCounts*6.;   // y-axis maximum for fast settings
  else hAllYMax = maxCounts*3.;   // y-axis maximum for standard and high gain settings
  TH2I *hAll = new TH2I("hAll","",10,hAllXMin,hAllXMax,10,hAllYMin,hAllYMax);
  SetHistProperties(hAll);
  char hAllTitle[100];
  sprintf(hAllTitle,"%s Settings : %s",settingsTitle,moduleID);
  hAll->SetTitle(hAllTitle);

  // Empty histogram for Mn target - no target counts
  double h3SubYMax;   // y-axis maximum is variable
  // Other axes limits are constants from first plot
  const double h3SubXMin = h1XMin;
  const double h3SubXMax = h1XMax;
  const double h3SubYMin = h1YMin;
  if(strcmp(settings,"fast") == 0) h3SubYMax = *max_element(&medianCountsMn[250],&medianCountsMn[500])*9.;   // y-axis maximum for fast settings
//  else h3SubYMax = *max_element(&medianCountsMn[250],&medianCountsMn[500])*3.;   // y-axis maximum for standard and high gain settings
//  TH2I *h3Sub = new TH2I("h3Sub","",10,h3SubXMin,h3SubXMax,10,h3SubYMin,h3SubYMax);
  TH2I *h3Sub = new TH2I("h3Sub","",10,h3SubXMin,h3SubXMax,10,h3SubYMin,h2YMax-h1YMax);
  SetHistProperties(h3Sub);
  char h3SubTitle[100];
  sprintf(h3SubTitle,"%s Settings : Mn target - no target : %s",settingsTitle,moduleID);
  h3Sub->SetTitle(h3SubTitle);

  // Declare and fill graphs with median counts for no, Mn, Cu & Zr targets
  TGraph *gr0 = new TGraph(nThresh,vthreshold,&medianCounts[0][0]);
  gr0->SetMarkerStyle(20);
  gr0->SetMarkerSize(2);
  TGraph *gr1 = new TGraph(nThresh,vthreshold,&medianCounts[1][0]);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(2);
  gr1->SetMarkerColor(2);
  TGraph *gr2 = new TGraph(nThresh,vthreshold,&medianCounts[2][0]);
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerSize(2);
  gr2->SetMarkerColor(3);
  TGraph *gr3 = new TGraph(nThresh,vthreshold,&medianCounts[3][0]);
  gr3->SetMarkerStyle(23);
  gr3->SetMarkerSize(2);
  gr3->SetMarkerColor(4);
  TGraph *gr4 = new TGraph(nThresh,vthreshold,&medianCounts[4][0]);
  gr4->SetMarkerStyle(33);
  gr4->SetMarkerSize(2);
  gr4->SetMarkerColor(5);

  // Clone graphs of median counts for no, Mn, Cu & Zr targets to change marker style for displaying them on same plot
  TGraph *gr1Fit = (TGraph*)gr1->Clone();
  gr1Fit->SetMarkerStyle(25);
  TGraph *gr2Fit = (TGraph*)gr2->Clone();
  gr2Fit->SetMarkerStyle(26);
  TGraph *gr3Fit = (TGraph*)gr3->Clone();
  gr3Fit->SetMarkerStyle(32);
  TGraph *gr4Fit = (TGraph*)gr4->Clone();
  gr4Fit->SetMarkerStyle(27);

  TF1 *fitCounts[4];   // Array of S-curve fits for no, Mn, Cu & Zr targets
  if(fit == true) {

    for(int i = 0 ; i < 4 ; i++) {

      char fitCountsName[10];
      sprintf(fitCountsName,"f%dFit",i+1);
      fitCounts[i] = new TF1(fitCountsName,"[0]/(1 + TMath::Exp(-[1]*(x - [2])))",minRange[i],maxRange[i]);
      fitCounts[i]->SetParameter(0,maximum[i]);   // From fit parameter file
      fitCounts[i]->SetParameter(2,midPoint[i]);   // From fit parameter file

    }

  }

  // Declare and fill graph with Mn target - no target counts
  TGraph *gr3Sub = new TGraph(nThresh,vthreshold,medianCountsMn);
  gr3Sub->SetMarkerStyle(26);
  gr3Sub->SetMarkerSize(2);
  gr3Sub->SetMarkerColor(3);

  // S-curve fit for Mn target - no target counts
  TF1 *fit3Sub = new TF1("f3SubFit","[3]+[0]/(1 + TMath::Exp(-[1]*(x - [2])))");
  if(fit == true) {

    fit3Sub->SetRange(minRange[4],maxRange[4]);
    fit3Sub->SetParameter(0,maximum[4]);
    fit3Sub->SetParameter(2,midPoint[4]);
    fit3Sub->SetParameter(3,offsetMnSubNo);

  }

  // Mean counts graphs
  double vthresholdDouble[nThresh];
  double vthresholdError[nThresh];
  for(int i = 0 ; i < nThresh ; i++) {
    vthresholdDouble[i] = (double)vthreshold[i];
    vthresholdError[i] = 0.0;
  }
  TGraphErrors *gr3Mean = new TGraphErrors(nThresh,vthresholdDouble,&meanCounts[0][0],vthresholdError,&meanCountsError[0][0]);
  gr3Mean->SetMarkerStyle(23);
  gr3Mean->SetMarkerSize(2);
  gr3Mean->SetMarkerColor(4);
  TGraphErrors *gr4Mean = new TGraphErrors(nThresh,vthresholdDouble,&meanCounts[1][0],vthresholdError,&meanCountsError[1][0]);
  gr4Mean->SetMarkerStyle(33);
  gr4Mean->SetMarkerSize(2);
  gr4Mean->SetMarkerColor(5);
  TGraph *gr3StandDev = new TGraph((int)h3XMax-(int)h3XMin,&vthresholdDouble[(int)h3XMin-offsetThresh],&standDevMeanCounts[0][(int)h3XMin-offsetThresh]);
  gr3StandDev->SetMarkerStyle(23);
  gr3StandDev->SetMarkerSize(2);
  gr3StandDev->SetMarkerColor(4);
  TGraph *gr4StandDev = new TGraph((int)h4XMax-(int)h4XMin,&vthresholdDouble[(int)h4XMin-offsetThresh],&standDevMeanCounts[1][(int)h4XMin-offsetThresh]);
  gr4StandDev->SetMarkerStyle(33);
  gr4StandDev->SetMarkerSize(2);
  gr4StandDev->SetMarkerColor(5);

  // Legend for plot with all graphs
  TLegend *leg = new TLegend(0.1,0.5,0.5,0.9);
  leg->SetFillColor(10);
  leg->AddEntry(gr0,"Shutter closed","p");
  leg->AddEntry(gr1,"Shutter open, no target","p");
  leg->AddEntry(gr2,"Shutter open, Mn target","p");
  leg->AddEntry(gr3,"Shutter open, Cu target","p");
  leg->AddEntry(gr4,"Shutter open, Zr target","p");

  // Text boxes for fit \chi^2/NDF and equation
  TPaveText *p1Fit = new TPaveText(h1XMin,0.7*h1YMax,h1XMin+(h1XMax-h1XMin)/2.,h1YMax);
  TPaveText *p2Fit = new TPaveText(h2XMin,0.7*h2YMax,h2XMin+(h2XMax-h2XMin)/2.,h2YMax);
  TPaveText *p3Fit = new TPaveText(h3XMin,0.7*h3YMax,h3XMin+(h3XMax-h3XMin)/2.,h3YMax);
  TPaveText *p4Fit = new TPaveText(h4XMin,0.7*h4YMax,h4XMin+(h4XMax-h4XMin)/2.,h4YMax);
//  TPaveText *p3SubFit = new TPaveText(h3SubXMin,0.7*h3SubYMax,h3SubXMin+(h3SubXMax-h3SubXMin)/2.,h3SubYMax);
  TPaveText *p3SubFit = new TPaveText(h3SubXMin,0.7*(h2YMax-h1YMax),h3SubXMin+(h3SubXMax-h3SubXMin)/2.,h2YMax-h1YMax);

  const double eMn = 5.90;   // Mn K\alpha_1 line
  const double eCu = 8.05;   // Cu K\alpha_1 line
  const double eZr = 15.77;   // Zr K\alpha_1 line
  double xEnergy[2][3] = { { eMn , eCu , eZr } , { 0.0 , 0.0 , 0.0 } };   // No errors in K\alpha_1 line energies
  double iPoints[2][3],iPointsMnCor[2][3];   // Arrays of inflection points (parameter 2 in S-curve fits)

  // Draw and fit to no, Mn, Cu & Zr target median counts graphs
  // Extract inflection points and associated errors where necessary
  c1->cd(1)->SetGrid();
  h1->Draw();
  gr1Fit->Draw("p");
  if(fit == true) {
    gr1Fit->Fit("f1Fit","R");
    SetFitTextLabel(p1Fit,fitCounts[0]);
    p1Fit->Draw();
  }
  c1->cd(2)->SetGrid();
  h2->Draw();
  gr2Fit->Draw("p");
  if(fit == true) {
    gr2Fit->Fit("f2Fit","R");
    iPoints[0][0] = fitCounts[1]->GetParameter(2);
    iPoints[1][0] = fitCounts[1]->GetParError(2);
    SetFitTextLabel(p2Fit,fitCounts[1]);
    p2Fit->Draw();
  }
  c1->cd(3)->SetGrid();
  h3->Draw();
  gr3Fit->Draw("p");
  if(fit == true) {
    gr3Fit->Fit("f3Fit","R");
    iPoints[0][1] = fitCounts[2]->GetParameter(2);
    iPoints[1][1] = fitCounts[2]->GetParError(2);
    iPointsMnCor[0][1] = iPoints[0][1];
    iPointsMnCor[1][1] = iPoints[1][1];
    SetFitTextLabel(p3Fit,fitCounts[2]);
    p3Fit->Draw();
  }
  c1->cd(4)->SetGrid();
  h4->Draw();
  gr4Fit->Draw("p");
  if(fit == true) {
    gr4Fit->Fit("f4Fit","R");
    iPoints[0][2] = fitCounts[3]->GetParameter(2);
    iPoints[1][2] = fitCounts[3]->GetParError(2);
    iPointsMnCor[0][2] = iPoints[0][2];
    iPointsMnCor[1][2] = iPoints[1][2];
    SetFitTextLabel(p4Fit,fitCounts[3]);
    p4Fit->Draw();
  }
  char outputFileNameP1[100];
  sprintf(outputFileNameP1,"/localhome/local/mythen/analysis/main/plots/%s_%s_eCal01.pdf",moduleID,settings);
  c1->Print(outputFileNameP1);

  // Declare and fill graph for intial energy calibration
  TGraphErrors *grCal = new TGraphErrors(3,&xEnergy[0][0],&iPoints[0][0],&xEnergy[1][0],&iPoints[1][0]);
  SetGraphProperties(grCal);
  char grCalTitle[100];
  sprintf(grCalTitle,"Energy Calibration : %s Settings : %s",settingsTitle,moduleID);
  grCal->SetTitle(grCalTitle);
  grCal->SetMarkerStyle(5);
  grCal->SetMarkerSize(3);

  TF1 *fCal3E = new TF1("fCal3E","[0]*x + [1]",xEnergy[0][0],xEnergy[0][2]);   // 3-point (Mn, Cu & Zr) initial calibration fit
  TF1 *fCal2E = new TF1("fCal2E","[0]*x + [1]",xEnergy[0][1],xEnergy[0][2]);   // 2-point (Cu & Zr) initial calibration fit
  fCal2E->SetLineColor(4);
  TF1 *fCalMnCor = new TF1("fCalMnCor","[0]*x + [1]",xEnergy[0][0],xEnergy[0][2]);   // 3-point (Mn, Cu & Zr) Mn corrected calibration fit

  // Character arrays for calibration fit equations
  char eqCal3E[100],eqCal2E[100];
  char eqCalMnCor[100];
  char noTargetE[100];

  // Draw all graphs on same plot
  c2->cd(1)->SetGrid();
  hAll->Draw();
  gr0->Draw("p");
  gr1->Draw("p");
  gr2->Draw("p");
  gr3->Draw("p");
  gr4->Draw("p");
  leg->Draw();
  // Draw and fit to Mn target - no target graph
  c2->cd(2)->SetGrid();
  h3Sub->Draw();
  gr3Sub->Draw("p");
  if(fit == true) {
    gr3Sub->Fit("f3SubFit","R");
    iPointsMnCor[0][0] = fit3Sub->GetParameter(2);   // Mn inflection point (corrected)
    iPointsMnCor[1][0] = fit3Sub->GetParError(2);   // Mn inflection point error (corrected)
    Set4ParFitTextLabel(p3SubFit,fit3Sub);
    p3SubFit->Draw();
  }
  if(fit == true) {   // Only draw on pad if fitting
    // Initial calibration
    c2->cd(3)->SetGrid();
    grCal->Draw("ap");
    grCal->Fit("fCal3E","R");
    SetCalEqString(eqCal3E,fCal3E);
    grCal->Fit("fCal2E","R+");
    SetCalEqString(eqCal2E,fCal2E);
    TText *tCal3E = new TText(6.,iPoints[0][2]+15.,eqCal3E);
    tCal3E->SetTextSize(0.035);
    tCal3E->Draw();
    TText *tCal2E = new TText(6.,iPoints[0][2]+5.,eqCal2E);
    tCal2E->SetTextSize(0.035);
    tCal2E->SetTextColor(4);
    tCal2E->Draw();
  }
  // Declare and fill Mn corrected calibration graph
  TGraphErrors *grCalMnCor = new TGraphErrors(3,&xEnergy[0][0],&iPointsMnCor[0][0],&xEnergy[1][0],&iPointsMnCor[1][0]);
  SetGraphProperties(grCalMnCor);
  char grCalMnCorTitle[100];
  sprintf(grCalMnCorTitle,"Mn Corrected Energy Calibration : %s Settings : %s",settingsTitle,moduleID);
  grCalMnCor->SetTitle(grCalMnCorTitle);
  grCalMnCor->SetMarkerStyle(5);
  grCalMnCor->SetMarkerSize(3);
  if(fit == true) {   // Only draw on pad if fitting
    c2->cd(4)->SetGrid();
    grCalMnCor->Draw("ap");
    grCalMnCor->Fit("fCalMnCor","R");
    SetCalEqString(eqCalMnCor,fCalMnCor);
    TText *tCalMnCor = new TText(6.,iPointsMnCor[0][2]+15.,eqCalMnCor);
    tCalMnCor->SetTextSize(0.035);
    tCalMnCor->Draw();
    SetNoTargetEnergy(noTargetE,fCalMnCor,fitCounts[0]);
    TLatex *tNoTargetE = new TLatex(6.,iPointsMnCor[0][2]+5.,noTargetE);
    tNoTargetE->SetTextSize(0.035);
    tNoTargetE->Draw();
  }
  char outputFileNameP2[100];
  sprintf(outputFileNameP2,"/localhome/local/mythen/analysis/main/plots/%s_%s_eCal02.pdf",moduleID,settings);
  c2->Print(outputFileNameP2);

  c3->cd(1)->SetGrid();
  h3Mean->Draw();
  gr3Mean->Draw("p");
  c3->cd(2)->SetGrid();
  h4Mean->Draw();
  gr4Mean->Draw("p");
  c3->cd(3)->SetGrid();
  c3->cd(3)->SetLogy();
  h3StandDev->Draw();
  gr3StandDev->Draw("p");
  c3->cd(4)->SetGrid();
  c3->cd(4)->SetLogy();
  h4StandDev->Draw();
  gr4StandDev->Draw("p");
  char outputFileNameP3[100];
  sprintf(outputFileNameP3,"/localhome/local/mythen/analysis/main/plots/%s_%s_eCalMean.pdf",moduleID,settings);
  c3->Print(outputFileNameP3);

  // Ask user to write fit parameters to file if fit argument was not given
  if(fit == false) {

    string fitParOpt;
    cout << "Please take a look at output file. Would you like to set fit parameters now? [Y/n] ";
    getline(cin,fitParOpt);

    if(fitParOpt == "N" || fitParOpt == "n" || fitParOpt == "No" || fitParOpt == "NO" || fitParOpt == "no") {

      cout << "Exiting..." << endl;

      return 0;

    } else SetFitParameterOptions(fitParameterFileName);

  }

  // Write gain and offset to calibration file if fits are performed
  if(fit == true) {

    ofstream calFile;
    char calFileName[100];
    sprintf(calFileName,"/localhome/local/mythen/analysis/main/settings/newcal_%s.sn%s",settings,moduleID);
    calFile.open(calFileName);

    // Write offset then gain to file
    // Extract parameters from Mn corrected calibration fit
    calFile << fCalMnCor->GetParameter(1) << " " << fCalMnCor->GetParameter(0) << endl;

    calFile.close();   // Tidy up

    // Re-use object to write to fit parameter file
    calFile.open(fitParameterFileName);

    for(int i = 0 ; i < 4 ; i++) calFile << minRange[i] << " " << maxRange[i] << " " << fitCounts[i]->GetParameter(0) << " " << fitCounts[i]->GetParameter(1) << " "
                                         << fitCounts[i]->GetParameter(2) << endl;

    calFile << minRange[4] << " " << maxRange[4] << " " << fit3Sub->GetParameter(0) << " " << fit3Sub->GetParameter(1) << " " << fit3Sub->GetParameter(2) << " "
            << fit3Sub->GetParameter(3) << endl;

    calFile.close();

  }

  // Tidy up
  delete c1;
  c1 = 0;
  delete c2;
  c2 = 0;
  delete c3;
  c3 = 0;
  delete h1;
  h1 = 0;
  delete h2;
  h2 = 0;
  delete h3;
  h3 = 0;
  delete h4;
  h4 = 0;
  delete h3Mean;
  h3Mean = 0;
  delete h4Mean;
  h4Mean = 0;
  delete h3StandDev;
  h3StandDev = 0;
  delete h4StandDev;
  h4StandDev = 0;
  delete hAll;
  hAll = 0;
  delete h3Sub;
  h3Sub = 0;
  delete gr0;
  gr0 = 0;
  delete gr1;
  gr1 = 0;
  delete gr2;
  gr2 = 0;
  delete gr3;
  gr3 = 0;
  delete gr4;
  gr4 = 0;
  delete gr1Fit;
  gr1Fit = 0;
  delete gr2Fit;
  gr2Fit = 0;
  delete gr3Fit;
  gr3Fit = 0;
  delete gr4Fit;
  gr4Fit = 0;
  delete gr3Sub;
  gr3Sub = 0;
  delete gr3Mean;
  gr3Mean = 0;
  delete gr4Mean;
  gr4Mean = 0;
  delete gr3StandDev;
  gr3StandDev = 0;
  delete gr4StandDev;
  gr4StandDev = 0;
  delete grCal;
  grCal = 0;
  delete grCalMnCor;
  grCalMnCor = 0;
  for(int i = 0 ; i < 4 ; i++) {
    delete fitCounts[i];
    fitCounts[i] = 0;
  }
  delete fit3Sub;
  fit3Sub = 0;
  delete fCal3E;
  fCal3E = 0;
  delete fCal2E;
  fCal2E = 0;
  delete fCalMnCor;
  fCalMnCor = 0;
  delete leg;
  leg = 0;
  delete p1Fit;
  p1Fit = 0;
  delete p2Fit;
  p2Fit = 0;
  delete p3Fit;
  p3Fit = 0;
  delete p4Fit;
  p4Fit = 0;
  delete p3SubFit;
  p3SubFit = 0;

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
inline int GetMedian(const int sortCounts[], const int & nStrips) {

  return (sortCounts[nStrips/2 - 1] + sortCounts[nStrips/2]) / 2;

}

// Get mean from sorted array of counts
double GetMean(const int sortCounts[], const int & nStrips) {

  double meanCounts = 0.0;
  for(int i = 0; i < nStrips ; i++) meanCounts += ((double)sortCounts[i] / (double)nStrips);
  return meanCounts;

}

// Get standard deviation from mean
double GetStandDev(const int sortCounts[], const double & meanCounts, const int & nStrips) {

  double standDev = 0.0;
  for(int i = 0 ; i < nStrips ; i++) standDev += ((((double)sortCounts[i] - meanCounts) * ((double)sortCounts[i] - meanCounts)) / (double)nStrips);
  return TMath::Sqrt(standDev);

}

// Set empty histogram properties
void SetHistProperties(TH2I * hist) {

  hist->GetXaxis()->SetTitle("vthreshold (DACu)");
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitle("Number of Counts / N");
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);

}

// Set empty histogram properties for standard deviation
void SetHistStandDevProperties(TH2I * hist) {

  hist->GetXaxis()->SetTitle("vthreshold (DACu)");
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitle("Standard Deviation from Mean / #sigma");
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);

}

// Set graph properties for calibration fits
void SetGraphProperties(TGraphErrors * graph) {

  graph->GetXaxis()->SetTitle("X-Ray Energy (keV) / E");
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetTitleSize(0.03);
  graph->GetXaxis()->SetLabelSize(0.03);
  graph->GetYaxis()->SetTitle("vthreshold (DACu) / V");
  graph->GetYaxis()->SetTitleOffset(1.4);
  graph->GetYaxis()->CenterTitle();
  graph->GetYaxis()->SetTitleSize(0.03);
  graph->GetYaxis()->SetLabelSize(0.03);

}

// Write fit parameter options to file from user input
void SetFitParameterOptions(const char * fileName) {

  int rangeMin = 0;
  int rangeMax = 0;
  int maximum = 0;
  int midPoint = 0;
  int offset = 0;
  const int steepness = 0;   // Not used to set initial fit parameters so user is not asked for it

  ofstream fitParFile;
  fitParFile.open(fileName);

  // Loop over all targets and finally (Mn target - no target)
  for(int i = 0 ; i < 5 ; i++) {

    if(i == 0) cout << "No target" << endl;
    if(i == 1) cout << "Mn target" << endl;
    if(i == 2) cout << "Cu target" << endl;
    if(i == 3) cout << "Zr target" << endl;
    if(i == 4) cout << "Mn target - no target" << endl;

    // Ask user for data
    cout << "Minimum range: ";
    cin >> rangeMin;
    cout << "Maximum range: ";
    cin >> rangeMax;
    cout << "Maximum: ";
    cin >> maximum;
    cout << "Midpoint: ";
    cin >> midPoint;
    if(i == 4) {
      cout << "Offset: ";
      cin >> offset;
    }

    // Display user's entries and ask user to confirm that settings are correct
    string acceptSettings;
    cout << endl << "Minimum range: " << rangeMin << endl
         << "Maximum range: " << rangeMax << endl
         << "Maximum: " << maximum << endl
         << "Midpoint: " << midPoint << endl;
    if(i == 4) cout << "Offset: " << offset << endl;
    cout << "Are these settings acceptable? [Y/n] ";
    cin.get();
    getline(cin,acceptSettings);

    // Decrement counter if settings are incorrect so user has another opportunity to enter them
    if(acceptSettings == "N" || acceptSettings == "n" || acceptSettings == "No" || acceptSettings == "NO" || acceptSettings == "no") i--;
    else {   // Write to file if settings are correct

      if(i == 4) fitParFile << rangeMin << " " << rangeMax << " " << maximum << " " << steepness << " " << midPoint << " " << offset << endl;   // (MN target - no target)
      else fitParFile << rangeMin << " " << rangeMax << " " << maximum << " " << steepness << " " << midPoint << endl;   // Individual targets

    }

  }

  fitParFile.close();   // Tidy up

}

// Write chi^2/NDF and equation from fit to text box
void SetFitTextLabel(TPaveText * pText, const TF1 * fit) {

  char fitChiSq[100],fitEq[100];
  pText->SetFillColor(10);
  SetChiSqString(fitChiSq,fit);
  SetEqString(fitEq,fit);
  pText->AddText(fitChiSq);
  pText->AddText(fitEq);

}

// Write chi^2/NDF and equation from fit to text box for (Mn target - no target)
void Set4ParFitTextLabel(TPaveText * pText, const TF1 * fit) {

  char fitChiSq[100],fitEq[100];
  pText->SetFillColor(10);
  SetChiSqString(fitChiSq,fit);
  Set4ParEqString(fitEq,fit);
  pText->AddText(fitChiSq);
  pText->AddText(fitEq);  

}

// Fill character array with fit's chi^2/NDF value
inline void SetChiSqString(char * eqString, const TF1 * fit) {

  sprintf(eqString,"#Chi^{2} / NDF = %.3E",fit->GetChisquare()/(double)fit->GetNDF());

}

// Fill character array with fit equation
inline void SetEqString(char * eqString, const TF1 * fit) {

  sprintf(eqString,"N(x) = #frac{%.0f}{1 + e^{-%.3f(x - %.2f)}}",fit->GetParameter(0),fit->GetParameter(1),fit->GetParameter(2));

}

// Fill character array with fit equation for (Mn target - no target)
inline void Set4ParEqString(char * eqString, const TF1 * fit) {

  sprintf(eqString,"N(x) = %.0f + #frac{%.0f}{1 + e^{-%.3f(x - %.2f)}}",fit->GetParameter(3),fit->GetParameter(0),fit->GetParameter(1),fit->GetParameter(2));

}

// Fill character array with calibration fit equation
inline void SetCalEqString(char * eqString, const TF1 * fit) {

  sprintf(eqString,"V = %.4fE + %.3f",fit->GetParameter(0),fit->GetParameter(1));

}

// Fill character array with threshold energy of no target inflection point
void SetNoTargetEnergy(char * energyString, const TF1 * calFit, const TF1 * sCurveFit) {

  double gain = calFit->GetParameter(0);
  double gainError = calFit->GetParError(0);
  double offset = calFit->GetParameter(1);
  double offsetError = calFit->GetParError(1);
  double inflectionPoint = sCurveFit->GetParameter(2);
  double inflectionPointError = sCurveFit->GetParError(2);

  double energy0T = (inflectionPoint - offset) / gain;

  double errVSubO = TMath::Sqrt((inflectionPointError * inflectionPointError) + (offsetError * offsetError));
  double energy0TError = energy0T * TMath::Sqrt(((errVSubO / (inflectionPoint - offset)) * (errVSubO / (inflectionPoint - offset))) + ((gainError / gain) * (gainError / gain)));

  sprintf(energyString,"no Target E_{t} = %.2f #pm %.2f keV",energy0T,energy0TError);

}
