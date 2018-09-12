// Plot untrimmed (top left) and trimmed threshold scan (bottom left)
// Plot Counts / Number of Channels from trimmed threshold scan data (bottom right)
// Plot Trimbits distribution from noise file (top right)
//
// Hassan Chagani 2016-10-25 11:31
//
// gcc v4.4.7
// root v5.99/01
//
// To compile:
// $: g++ -o ../bin/threshScan threshScan.cc `root-config --cflags --glibs`
// To run:
// $: ../bin/threshScan <untrimmedFileIndex> <trimmedFileIndex> <settings> <moduleID>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

void HistSettings2D(TH2I * hist);
void HistSettings1D(TH1I * hist);
void TrimBits(char * fileName, TH1I * const hist, const int &);

#ifndef __CINT__
int main(int argc, char *argv[]) {

  if(argc != 5) {

    cout << "Insufficient arguments!" << endl << endl
         << "Usage:" << endl
         << "./threshScan <untrimmedFileIndex> <trimmedFileIndex> <settings> <moduleID>" << endl;

    return 1;

  }

  // Assign command line arguments to constants and variables
  const int untrimmedFileIndex = atoi(argv[1]);
  const int trimmedFileIndex = atoi(argv[2]);
  char settings[20],settingsTitle[20];
  if(strcmp(argv[3],"standard") == 0 || strcmp(argv[3],"Standard") == 0) {
    strcpy(settings,"standard");
    strcpy(settingsTitle,"Standard");
  } else if(strcmp(argv[3],"highgain") == 0 || strcmp(argv[3],"Highgain") == 0 || strcmp(argv[3],"HighGain") == 0) {
    strcpy(settings,"highgain");
    strcpy(settingsTitle,"High Gain");
  } else if(strcmp(argv[3],"fast") == 0 || strcmp(argv[3],"Fast") == 0) {
    strcpy(settings,"fast");
    strcpy(settingsTitle,"Fast");
  } else {
    cout << "Settings argument <" << argv[3] << "> incorrect. Exiting..." << endl;
    return 1;
  }
  char moduleID[10];
  strcpy(moduleID,argv[4]);

  // Check command line settings are correct
  // Exit if not
  string correctSettings;
  cout << "Settings: " << endl
       << untrimmedFileIndex << "\t\tUntrimmed File Index" << endl
       << trimmedFileIndex << "\t\tTrimmed File Index" << endl
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

  // PART 1 : READ IN DATA FROM FILES

  // Filename string declarations and initial handling
  ifstream dataFile;

  char fileNamePrefix[100],fileNamePostfix[5];
  char currentFileName[100];
  strcpy(fileNamePrefix,"/localhome/local/Desktop/OutputDir/run_S");
  strcpy(fileNamePostfix,".raw");

  // Arrays for data
  const int nThresh = 401;
  const int offsetThresh = 500;
  const int nStrips = 1280;
  int strip[nStrips],vthreshold[nThresh];
  int untrimCounts[nThresh][nStrips],trimCounts[nThresh][nStrips];

  // Loop over all files where i is vthreshold
  for(int i = 0 ; i < nThresh ; i++) {

    // Filename string handling for untrimmed data
    strcpy(currentFileName,fileNamePrefix);
    char fileID[10];
    sprintf(fileID,"%d_%d",i+offsetThresh,untrimmedFileIndex);
    strcat(currentFileName,fileID);
    strcat(currentFileName,fileNamePostfix);

    dataFile.open(currentFileName);
    if(i%10 == 0) cout << "File: " << currentFileName << " ";

    int skipStrip;   // Throwaway variable to read in strip number

    if(i + offsetThresh == 500) {   // Fill both strip and untrimCounts array from initial file

      for(int j = 0 ; j < nStrips ; j++) dataFile >> strip[j] >> untrimCounts[i][j];

    } else {   // Only fill untrimCounts array from subsequent files

      for(int j = 0 ; j < nStrips ; j++) dataFile >> skipStrip >> untrimCounts[i][j];

    }

    vthreshold[i] = i + offsetThresh;   // Fill vthreshold array

    // Close untrimmed data file and open trimmed data file
    dataFile.close();

    strcpy(currentFileName,fileNamePrefix);
    sprintf(fileID,"%d_%d",i+offsetThresh,trimmedFileIndex);
    strcat(currentFileName,fileID);
    strcat(currentFileName,fileNamePostfix);

    dataFile.open(currentFileName);
    if(i%10 == 0) cout << currentFileName << endl;

    // Only fill trimCounts array from data files
    // strip array is already filled
    for(int j = 0 ; j < nStrips ; j++) dataFile >> skipStrip >> trimCounts[i][j];

    dataFile.close();   // Tidy up

  }

  // PART 2 : PLOT DATA

  TCanvas *c1 = new TCanvas("c1","",2970,2100);
  c1->Divide(2,2);

  // Untrimmed threshold scan histogram
  TH2I *h1 = new TH2I("h1","",nStrips,-0.5,nStrips-0.5,nThresh,offsetThresh-0.5,nThresh+offsetThresh+0.5);
  h1->GetXaxis()->SetTitle("Channel Number");
  h1->GetYaxis()->SetTitle("vthreshold (DACu)");
  HistSettings2D(h1);

  // Clone untrimmed threshold scan histogram for trimmed threshold scan
  TH2I *h2 = (TH2I*)h1->Clone("h2");

  // Histogram for counts / number of channels to investigate vthreshold
  TH1I *h3 = new TH1I("h3","",nThresh,offsetThresh-0.5,nThresh+offsetThresh+0.5);
  h3->GetXaxis()->SetTitle("vthreshold (DACu)");
  h3->GetYaxis()->SetTitle("#frac{Summed Counts from All Channels}{Number of Channels}");
  HistSettings1D(h3);

  char h1Title[100],h2Title[100];
  sprintf(h1Title,"Untrimmed Threshold Scan : %s Settings : %s",settingsTitle,moduleID);
  h1->SetTitle(h1Title);
  sprintf(h2Title,"Trimmed Threshold Scan : %s Settings : %s",settingsTitle,moduleID);
  h2->SetTitle(h2Title);
  h3->SetTitle("Summed Counts from All Channels After Trimming");

  for(int i = 0 ; i < nThresh ; i++) {   // Loop over thresholds

    for(int j = 0 ; j < nStrips; j++) {   // Loop over strips

      // Fill histograms with associated weights
      h1->Fill(strip[j],vthreshold[i],untrimCounts[i][j]);
      h2->Fill(strip[j],vthreshold[i],trimCounts[i][j]);
      h3->Fill(vthreshold[i],trimCounts[i][j]);

    }

  }

  // Create histogram for trimbits and call function to read data and fill it
  char trimBitsFileName[100];
  sprintf(trimBitsFileName,"/localhome/local/slsDetectorsPackage/settingsdir/mythen/%s/noise.sn%s",settings,moduleID);
  TH1I *h4 = new TH1I("h4","",64,-0.5,63.5);
  char h4Title[100];
  sprintf(h4Title,"Distribution of Trimbits : %s Settings : %s",settingsTitle,moduleID);
  h4->SetTitle(h4Title);
  TrimBits(trimBitsFileName,h4,nStrips);

  c1->cd(1)->SetGrid();
  h1->Draw("colz");   // Untrimmed threshold scan
  c1->cd(3)->SetGrid();
  h2->Draw("colz");   // Trimmed threshold scan
  c1->cd(4)->SetGrid();
  c1->cd(4)->SetLogy();
  h3->Scale(1./(float)nStrips);   // Counts / number of channels
  h3->Draw("hist");
  c1->cd(2)->SetGrid();
  h4->Draw("hist");   // Trimbits

  // Print canvas to file in pdf format
  char outputFileName[100];
  sprintf(outputFileName,"/localhome/local/mythen/analysis/main/plots/%s_%s_trimming.pdf",moduleID,settings);
  c1->Print(outputFileName);

  // Tidy up
  delete c1;
  c1 = 0;
  delete h1;
  h1 = 0;
  delete h2;
  h2 = 0;
  delete h3;
  h3 = 0;
  delete h4;
  h4 = 0;

  return 0;

}
#endif

// Set empty 2D histogram properties
void HistSettings2D(TH2I * hist) {

  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetZaxis()->SetLabelSize(0.03);

}

// Set empty 1D histogram properties
void HistSettings1D(TH1I * hist) {

  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);

}

// Read in data from trimbits file
// Fill histogram with data
void TrimBits(char * fileName, TH1I * const hist, const int & nStrips) {

  ifstream trimBitsDataFile;
  trimBitsDataFile.open(fileName);
  cout << "TrimBits File: " << fileName << endl;

  // Two variables for skipping data entry
  string skipString;
  int skipInt;
  int trimBits[nStrips];   // Array for data

  // Loop over data file
  // i goes over individual values in data file
  // j counts lines in data file
  // k corrects for lines where all data is skipped
  for(int i = 0 , j = 0 , k = 6 ; j < 1296 ; i++ , j++) {

    if(j%129-6 == 0) k++;   // Increment k when outBuffEnable line reached

    // First 6 lines and outBuffEnable lines to be skipped completely
    if(j < 6 || j%129-6 == 0) trimBitsDataFile >> skipString >> skipInt;
    else {

      // Trimbits are first integer in remaining lines
      // Values from five integers that follow to be thrown away
      trimBitsDataFile >> trimBits[j-k] >> skipInt >> skipInt >> skipInt >> skipInt >> skipInt;

    }

  }

  hist->GetXaxis()->SetTitle("Trimbits");
  hist->GetYaxis()->SetTitle("Number of Channels");
  HistSettings1D(hist);
  
  // Fill histogram with trimbits from array
  for(int i = 0 ; i < nStrips ; i++) hist->Fill(trimBits[i]);

}
