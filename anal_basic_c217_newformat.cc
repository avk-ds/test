/*
  Raj Comments:
  Need to chnage trigger layer
  also change layer Orddering, now 11 RPC data was routed to 10th layer DFE

*/


/*
  g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -o anal_basic_c217_newformat EveTree.C StraightLineFit.cc anal_basic_c217_newformat.cc `root-config --cflags --libs` -lMinuit


*/


#define ISALIGN


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <new>
#include<climits>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "EveTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"

#include "TStyle.h"
#include "TAttFill.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TPostScript.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"

#include "StraightLineFit.h"
using namespace std;

bool fDebug = false;
const double cval = 29.979; // velocity of light in cm/ns


int layerOrder[12]={0,1,2,3,4,5,6,7,8,9,11,10};
int layer4stripNumbering[32] = {16,17,2,3,4,5,6,7,18,19,10,11,12,13,14,15,20,21,18,19,20,21,22,23,22,23,26,27,28,29,30,31};                                                                                   



const int  nlayer =12; //Maximum RPC layers
const int  nstrip =32; // Strips in RPC chamber
const int  firstXstrip =0;  //1st strip in X
const int  lastXstrip =31; // last strip in X

const int layfirst = 0;
const int laylast = 11;

const int  firstYstrip =0;  //1st strip in Y
const int  lastYstrip =31; // last strip in Y


const double  accprng =1.0; //Acceptance region with respect tot strip centre
const double  effirng =0.25; //Additional region for efficeincy calculation


const double pival=acos(-1);
const int nmxhits=4;
const int nmxusedhits=3;


const double stripwidth = 3.0; // cm
const double stripgap = 0.2; // cm
const double layergap = 16.0; // cm

// const int trigly1 =5; //0;extpo
// const int trigly2 =7;
// const int trigly3 =8;
// const int trigly4 =9; //11;

const int trigly1 =6; //0;extpo
const int trigly2 =7;
const int trigly3 =8;
const int trigly4 =11; //11;

int trigger_layers[4] = {trigly1, trigly2, trigly3, trigly4};

const float xyPosDev=7.0;
const double  mxchisq =2.0;   //Maximum value of Normalised chi^2 (position fit);

const int nmnhits =4; //6; //Minimum number of layers required for position and time measurements


#ifdef ISALIGN
const int nmxiter = 1;
#else
const int nmxiter = 5;
#endif

const double layerzpos[nlayer]={0.0, 16.25, 32.05, 48.3, 64.2, 80.3, 96.4, 112.1, 128.2, 144.4, 160.3, 176.05};

float errxco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
float erryco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};



double xposerrsq[nmxhits][nlayer] = {
  {0.0745597, 0.0601179, 0.0501015, 0.0394726, 0.0350947, 0.055039, 0.0390593, 0.0372429, 0.0398186, 0.0568854, 0.04426, 0.0647753},
  {0.0459906, 0.0543877, 0.0371751, 0.0559855, 0.0453876, 0.0360321, 0.0656442, 0.0641278, 0.043064, 0.053713, 0.0444662, 0.058869},
  {0.162488, 0.109485, 0.0785132, 0.0824725, 0.0803276, 0.0997201, 0.0863911, 0.0762118, 0.0899414, 0.0973464, 0.0703962, 0.0933884},
  {0.957458, 0.751783, 0.666736, 0.637636, 0.647876, 0.757628, 0.472752, 0.623146, 0.766528, 0.807635, 0.673675, 0.884292}
};

double yposerrsq[nmxhits][nlayer] = {
  {0.0459599, 0.0382816, 0.0381181, 0.0500019, 0.0393566, 0.0400298, 0.0522574, 0.0332281, 0.0377729, 0.0391096, 0.0295895, 0.0408621},
  {0.033369, 0.0253477, 0.0479287, 0.042458, 0.0278329, 0.0277059, 0.0528327, 0.038016, 0.0188255, 0.0222188, 0.0273548, 0.0275486},
  {0.092609, 0.112, 0.0768894, 0.0688024, 0.118428, 0.108849, 0.0775983, 0.0854793, 0.0995298, 0.110006, 0.071602, 0.0937076},
  {0.689821, 0.759425, 0.598038, 0.607062, 0.748, 0.710522, 0.621418, 0.675569, 0.720422, 0.719581, 0.702963, 0.877633}
};

// float xoff[12]={0.0};
// float yoff[12]={0.0};

Double_t gausX(Double_t* x, Double_t* par){
  return par[0]*(TMath::Gaus(x[0], par[1], par[2], kFALSE)); //kTRUE));
}


void GetXposInStrip(double* ext, double* off, double* pos) {
  for (int ij=0; ij<nlayer; ij++) { 
    int istr = int(ext[ij]+off[ij]);
    if (ext[ij]+off[ij]<0.0) {
      pos[ij] = ext[ij]+off[ij] - istr + 0.5;
    } else {
      pos[ij] = ext[ij]+off[ij] - istr - 0.5;
    }
    //is it corrected
    //istr = int(ext[ij]);
    //  pos[ij] = ext[ij] - istr - 0.5;
  }
}
double GetAverage(TH2* hist) {

  double sum = 0.0;
  int totalBins = 0;

  for (int xBin = 1; xBin <= hist->GetNbinsX(); ++xBin) {
    for (int yBin = 1; yBin <= hist->GetNbinsY(); ++yBin) {
      double content = hist->GetBinContent(xBin, yBin);
      sum += content;
      if(content>0.)totalBins++;
    }
  }

  return (totalBins > 0) ? sum / totalBins : 0.0;
}





int main() {

  gStyle->SetPalette(1,0);
  gStyle->SetFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatStyle(1001);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);

  gStyle->SetStatFont(22);        // Times New Roman
  gStyle->SetTextFont(22);        // Times New Roman
  gStyle->SetTitleFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelSize(0.06, "XYZ"); // Times New Roman  
  gStyle->SetNdivisions(606, "XYZ");

  gStyle->SetOptTitle(0);
  gStyle->SetFuncWidth(1);
  gStyle->SetFuncColor(2);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(101);
  gStyle->SetOptLogy(0);
  gStyle->SetStatW(.18);
  gStyle->SetStatH(.08);
  gStyle->SetPadTopMargin(.02); //0.09
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.02);
  gStyle->SetPadRightMargin(0.02);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.12);
  latex.SetTextFont(42);
  latex.SetTextAlign(1); //(31); // align right

  char name[100];
  char title[100];
  
  char rootfiles[100];
  char outfilx[100];
  char outfile[100];
  char datafile[200];
  char infile[200];


  

  sprintf(rootfiles, "test_basic.log");
  int len = strlen(rootfiles);
  strncpy(outfilx, rootfiles, len-4);
  outfilx[len-4]='\0';
  sprintf(outfilx, "%s", outfilx);
  len = strlen(outfilx);
  outfilx[len]='\0';

#ifdef ISALIGN
  /*
    sprintf(outfile, "fileOut_INORUN_20241206_180901/Alignment.root", outfilx);
    TFile* AlignfileIn = TFile::Open(outfile);
    TTree* AlignTree = (TTree*)AlignfileIn->Get("AlignTree");
    AlignTree->SetBranchAddress("xoff",xoff);
    AlignTree->SetBranchAddress("yoff",yoff);

    AlignTree->GetEntry(0);
    cout<<"double xoff[12] = {";
    for (int ij=0; ij<nlayer; ij++) {
    cout<<xoff[ij];
    if(ij<nlayer-1)cout<<", ";

    }
    cout<<"};"<<endl;

    cout<<"double yoff[12] = {";
    for (int ij=0; ij<nlayer; ij++) {
    cout<<yoff[ij];
    if(ij<nlayer-1)cout<<", ";
    }
    cout<<"};"<<endl;
  */

#else 

 

  sprintf(outfile, "Alignment.root", outfilx);
  TFile* AlignfileOut = new TFile(outfile, "recreate");
  TTree* AlignTree = new TTree("AlignTree","AlignTree");
  AlignTree->Branch("xoff",xoff,"xoff[12]/F");
  AlignTree->Branch("yoff",yoff,"yoff[12]/F");


#endif

  double xoff[nlayer] = {-0.150392, -0.140311, -0.024575, 0.0918862, -0.0481311, -0.147595, 0.0374531, -0.01373, 0.038743, 0.0377367, 0, -0.0490905};
  double yoff[nlayer] = {0, 0.0583675, 0, 0.149097, 0, -0.0762682, -0.108328, 0.246808, -0.0333007, 0.0331146, 0, -0.000868276};

  
  double xxerr[nlayer], yyerr[nlayer];
  for ( int ix=0; ix<nlayer; ix++) {xxerr[ix] = yyerr[ix] = errxco[ix]*errxco[ix];}





  sprintf(outfile, "%s.root", outfilx);
  TFile* fileOut = new TFile(outfile, "recreate");
  
  sprintf(outfile, "%s.ps", outfilx);
  TPostScript ps(outfile,111);  
  ps.Range(20,30); //ps.Range(10,20);

  sprintf(outfile, "%s.txt", outfilx);
  ofstream file_out(outfile);
  


  TH1F* occu_x[nlayer];
  TH1F* occu_y[nlayer];
  TH2F* raw_occu[nlayer];
  
  TH1F* xlayer_mult[nlayer];
  TH1F* ylayer_mult[nlayer];

  for (int ij=0; ij<nlayer; ij++) {
    sprintf(title, "occu_x%i", ij);
    occu_x[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "occu_y%i", ij);
    occu_y[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
    
    sprintf(title, "raw_occu_l%i", ij);
    raw_occu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "xlayer_mult_l%i", ij);
    xlayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
    
    sprintf(title, "ylayer_mult_l%i", ij);
    ylayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

  }

  TH1F* xlayer_reso[nlayer][2*nmxiter];
  TH1F* ylayer_reso[nlayer][2*nmxiter];

  TH1F* xlayer_reso_mul[nlayer][2*nmxiter][nmxhits];
  TH1F* ylayer_reso_mul[nlayer][2*nmxiter][nmxhits];


  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "xlayer_reso_l%i_i%i", ij, jk);
      xlayer_reso[ij][jk]=new TH1F(title, title, 150, -6.0, 6.0);
      sprintf(title, "ylayer_reso_l%i_i%i", ij, jk);
      ylayer_reso[ij][jk]=new TH1F(title, title, 150, -6.0, 6.0);

      for (int kl=0; kl<nmxhits; kl++) {
        sprintf(title, "xlayer_reso_l%i_i%i_mul%i", ij, jk, kl+1);
        xlayer_reso_mul[ij][jk][kl]=new TH1F(title, title, 150, -6., 6.);
        sprintf(title, "ylayer_reso_l%i_i%i_mul%i", ij, jk, kl+1);
        ylayer_reso_mul[ij][jk][kl]=new TH1F(title, title, 150, -6., 6.);
        
      } 

    }
  }

  TH2F* totalentry[nlayer][2*nmxiter];
  TH2F* triggereffi_x[nlayer][2*nmxiter];
  TH2F* triggereffi_y[nlayer][2*nmxiter];
  
  TH2F* effi_x[nlayer][2*nmxiter];
  TH2F* effi_y[nlayer][2*nmxiter];
  
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "totalentry_l%i_i%i", ij,jk);
    totalentry[ij][jk]=new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "triggereffi_x_l%i_i%i", ij,jk);
    triggereffi_x[ij][jk]=new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    sprintf(title, "triggereffi_y_l%i_i%i", ij,jk);
    triggereffi_y[ij][jk]=new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "effi_x_l%i_i%i", ij,jk);
    effi_x[ij][jk]=new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    sprintf(title, "effi_y_l%i_i%i", ij,jk);
    effi_y[ij][jk]=new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5); 
    }
  }

  TH2D* strp_xmul[nlayer][2*nmxiter];
  TH2D* strp_ymul[nlayer][2*nmxiter];

  for (int ij=0; ij<nlayer; ij++) { 
    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "strp_xmul_l%i_i%i", ij,jk);
    strp_xmul[ij][jk] = new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);      
    sprintf(title, "strp_ymul_l%i_i%i", ij,jk);
    strp_ymul[ij][jk] = new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);
    }
  }
  

  int nentrymx;

  int xhits[nlayer],yhits[nlayer];
  double xrms[nlayer]={0};
  double yrms[nlayer]={0};

  int nTotalp=0;
  int firstiter = true;
  int occulyr=12;



  //*********************************************************************************************
  //               Iteration Starts
  //*********************************************************************************************
  int iiterrs = 0;
  const int ntcormx = 2;

  for (int iiter=0; iiter<nmxiter; iiter++) {

    file_out <<"iiter: "<< iiter<<endl;
    //file_out <<"iiter: "<< iiter<<endl;

    

    for(int ntcor = 0; ntcor < ntcormx; ntcor++){

      iiterrs = nmxiter * ntcor + iiter;//group all removing the layer form fit first for diff iiter
      
      int lstr = 0;
      int lend = nlayer;

      
      file_out<<lstr<<" "<<lend<<endl;

      if(ntcor==1){lend=lstr + 1;}
      
      for(int laye = lstr; laye < lend; laye++){
	
	file_out<<"Removing layer "<<laye<<" from fit "<<endl;
	if(fDebug)cout<<"Removing layer "<<laye<<" from fit "<<endl;
    
	occulyr = (ntcor==0) ? laye : nlayer;
	ifstream file_db;
	file_db.open(rootfiles);  
      
	while(!(file_db.eof())) {
	  file_db >> datafile>>nentrymx;
	  if(fDebug)cout<<"datafile "<<datafile<<endl;
    
	  if (strstr(datafile,"#")) continue;
	  if(file_db.eof()) break;
    
	  sprintf(infile, "newformat/%s", datafile);
	  TFile *fileIn = new TFile(infile, "read");
	  TTree *event_tree= (TTree*)fileIn->Get("evetree");
    
	  EveTree *event = new EveTree(event_tree);
	  event->Loop();
    
	  int nentry = event_tree->GetEntries();        
	  //nentry = 1000;
	  nentry = min(nentry,nentrymx);
	  if(fDebug)cout <<infile<<" has "<< nentry<<" events "<<endl;
    

	  for(int ievt=0;ievt<nentry;ievt++) { 
	    //file_out<<"ievt: "<<ievt<<endl;
	    fileIn->cd();
	    event_tree->GetEntry(ievt);
     
	    vector<int> xptsall[nlayer];
	    vector<int> yptsall[nlayer];
	    vector<int> xpts[nlayer];
	    vector<int> ypts[nlayer];
      
      
	    for(int jk=0;jk<nlayer;jk++) {
	      xptsall[jk].clear(); yptsall[jk].clear();
	      xpts[jk].clear(); ypts[jk].clear();
	    }




      
	    for(int jk=0;jk<nlayer;jk++) {
	      for(int kl=0; kl<nstrip; kl++) {
		//cout<<ievt<<" "<<jk<<" "<<kl<<endl;
		if(event->xLayer[layerOrder[jk]]->TestBitNumber(kl)) {
		  xptsall[jk].push_back(kl);
		}                
		if(event->yLayer[layerOrder[jk]]->TestBitNumber(kl)) {
		  // if(jk==3)yptsall[jk].push_back(layer4stripNumbering[kl]);
		  // else yptsall[jk].push_back(kl);
		  yptsall[jk].push_back(kl);


		}

	      }
	    }

	    bool istrigger = true;

	    for (int trig = 0; trig < 4; trig++) {
	      if (xptsall[trigger_layers[trig]].size() == 0) {
		istrigger = false; // At least one layer is missing hits
		break;
	      }
	    }


	    //if (!istrigger) continue; 
      
      



	    if(firstiter){      
	      nTotalp++;
	

	
	      //Fill Raw Occupancy
	      for(int jk=0;jk<nlayer;jk++) {
		for (int kl=0; kl< xptsall[jk].size(); kl++) {
		  occu_x[jk]->Fill(xptsall[jk][kl]); 
		  for (int lm=0; lm< yptsall[jk].size(); lm++) {
		    raw_occu[jk]->Fill(xptsall[jk][kl], yptsall[jk][lm]);
		  }
		}
	      }    
	      for(int jk=0;jk<nlayer;jk++) {
		for (int kl=0; kl< yptsall[jk].size(); kl++) {
		  occu_y[jk]->Fill(yptsall[jk][kl]);
		}
	      }
      
	      //Fill Multiplicity
	      for(int jk=0;jk<nlayer;jk++) {
		xlayer_mult[jk]->Fill(xptsall[jk].size());
		ylayer_mult[jk]->Fill(yptsall[jk].size());
	      }

	    }
      

	    ////////////////////////////////////////////
	    //
	    //  First clean up noisey layers
	    //
	    ////////////////////////////////////////////
      

	    for (int iz=0; iz<nlayer; iz++) {

	      for (int ix=0; ix<xptsall[iz].size(); ix++) {
		xpts[iz].push_back(xptsall[iz][ix]);	  
	      }
	      xhits[iz] = xpts[iz].size();
	      //file_out<<"xhits[iz] "<<iz<<" "<<xhits[iz]<<endl;
	      for (int ix=0; ix<yptsall[iz].size(); ix++) {
		ypts[iz].push_back(yptsall[iz][ix]);	  
	      }
	      yhits[iz] = ypts[iz].size();

	    }
     

	    double Xpos[nlayer] = {0.0}; 
	    bool Xusedpos[nlayer];

      
      
	    for (int ij=0;ij<nlayer;ij++) {
	      if (xhits[ij]<=0 || xhits[ij]>nmxhits) {
		Xpos[ij]= -100;
	      } else {
		for (int ix=0; ix<xhits[ij]; ix++) {
		  if (ix<xhits[ij]-1 && abs(xpts[ij][ix]-xpts[ij][ix+1])>1) { Xpos[ij]=-100; break;}
		  Xpos[ij] +=xpts[ij][ix];
		}
		if (Xpos[ij]>=0.0) {
		  Xpos[ij]  = Xpos[ij]/xhits[ij] + 0.5 - xoff[ij];
		  xxerr[ij] = xposerrsq[xhits[ij]-1][ij];
	    
		}
	      }

      
	      //Sort out hits, which can be used for fit
	      for(int ij=0;ij<nlayer;ij++){
		Xusedpos[ij] = (Xpos[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false;
	      }

      
	      if(fDebug)cout<<"Event No.: "<<ievt<<" Layer: "<<ij<<" "<<xptsall[ij].size()<<" "<<xhits[ij]<<endl;
	      for (int ix=0; ix<xhits[ij]; ix++) {
		if(fDebug)cout<<xpts[ij][ix]<<" ";
	      }
	      if(fDebug)cout<< Xpos[ij]<<" "<<Xusedpos[ij]<<endl;
	    }

	    int Nx=0;
	    int nxfail = 0;
	    double xchi2 = 0;
	    double Xdev[nlayer];for (int ij=0; ij<nlayer; ij++) { Xdev[ij] = 100;}
      
	    double zval[nlayer],xext[nlayer],xexter[nlayer],xposinstr[nlayer];
	    double xslope = 0.;
	    double xinters = 0.;
	    //for (int ix=0; ix<nlayer; ix++) { zval[ix]=ix;}
	    for (int ix=0; ix<nlayer; ix++) { zval[ix]=layerzpos[ix];}
      
	    for (int ix=0; ix<nlayer; ix++) { xext[ix]= xexter[ix] =xposinstr[ix] =  100;}

	    if(fDebug)cout<<"Fitting occulyr: "<<occulyr<<endl;
	    StraightLineFit xposfit(1,zval,Xpos,xxerr,Xusedpos,occulyr,occulyr,layfirst,laylast,xyPosDev);
	    xposfit.GetParameters(nxfail,xinters,xslope);
	    xposfit.GetChisqure(Nx,xchi2);
	    xposfit.GetFitValues(xext,Xdev,xexter);

	    GetXposInStrip(xext, xoff, xposinstr);

	    for (int ij=0;ij<nlayer;ij++) {
	      //if(occulyr>=nlayer) cout<<"xpos: "<<Xpos[ij]<<" "<<zval[ij]<<" "<<occulyr<<" "<<xxerr[ij]<<" "<<xext[ij]<<" "<<Xdev[ij]<<" "<<xexter[ij]<<" "<<xposinstr[ij]<<" "<<Xusedpos[ij]<<" "<<xyPosDev<<endl;
	    }
	    //if(occulyr>=nlayer) cout<<"Fit Results:: Nx: "<<Nx<<" xchi2: "<<xchi2<<" "<<nxfail<<endl;
	   
	   
	   for (int ij=0; ij<nlayer; ij++) {
	      if(fDebug)cout<<zval[ij]<<" "<<xext[ij]<<" "<<Xpos[ij]<<" "<<Xdev[ij]<<" "<<xext[ij]-Xpos[ij]<<endl;
	    }
	    //Fill Resolutions
	   //cout<<"ievt: "<<ievt<<" "<<Nx<<" "<<nmnhits<<" "<<xchi2<<" "<<nxfail<<endl;
	   if (Nx>= nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {
	     //cout<<"ievt: "<<ievt<<" "<<Nx<<" "<<nmnhits<<" "<<xchi2<<" "<<nxfail<<endl;
	     if(occulyr>=nlayer){
	       for (int ij=0; ij<nlayer; ij++) {
		 //cout<<ij<<" "<<xexter[ij]<<" "<<xposinstr[ij]<<" "<<xhits[ij]<<endl;
		 if (xexter[ij]<0.5) {
		   strp_xmul[ij][iiterrs]->Fill(xposinstr[ij], xhits[ij]);
		   //if(ij==0)cout<<"ievt:"<<ievt<<" "<<ij<<" "<<iiterrs<<" "<<xposinstr[ij]<<" "<<xhits[ij]<<endl;
		 }
		 
		  if (abs(Xdev[ij])<6.0 && Xusedpos[ij]) {
		    xlayer_reso[ij][iiterrs]->Fill(Xdev[ij]);
		    int ixx=min(xhits[ij], nmxhits);
		    if (ixx>0 && Xpos[ij]>1.0 && Xpos[ij]<nstrip-1) { 
		      xlayer_reso_mul[ij][iiterrs][ixx-1]->Fill(Xdev[ij]);
		    }//if (ixx>0 && Xpos[ij]>1.0 && Xpos[ij]<nstrip-1) { 
		  }
		}
	      }
	      else{
		if (xexter[laye]<0.5) {strp_xmul[laye][iiterrs]->Fill(xposinstr[laye], xhits[laye]);}
		if (abs(Xdev[laye])<6.0 && Xusedpos[laye]) {
		  xlayer_reso[laye][iiterrs]->Fill(Xdev[laye]);
		  int ixx=min(xhits[laye], nmxhits);
		  if (ixx>0 && Xpos[laye]>1.0 && Xpos[laye]<nstrip-1) { 
		    xlayer_reso_mul[laye][iiterrs][ixx-1]->Fill(Xdev[laye]);
		  }//if (ixx>0 && Xpos[laye]>1.0 && Xpos[laye]<nstrip-1) { 
		}
	      }

	    } 
      

      

	    double Ypos[nlayer] = {0.0}; 
	    bool Yusedpos[nlayer];

      
      
	    for (int ij=0;ij<nlayer;ij++) {
	      if (yhits[ij]<=0 || yhits[ij]>nmxhits) {
		Ypos[ij]= -100;
	      } else {
		for (int iy=0; iy<yhits[ij]; iy++) {
		  if (iy<yhits[ij]-1 && abs(ypts[ij][iy]-ypts[ij][iy+1])>1) { Ypos[ij]=-100; break;}
		  Ypos[ij] +=ypts[ij][iy];
		}
		if (Ypos[ij]>=0.0) {
		  Ypos[ij]  = Ypos[ij]/yhits[ij] + 0.5 - yoff[ij];
		  yyerr[ij] = yposerrsq[yhits[ij]-1][ij];
	    
		}
	      }

      
	      //Sort out hits, which can be used for fit
	      for(int ij=0;ij<nlayer;ij++){
		Yusedpos[ij] = (Ypos[ij]>-100 && yhits[ij]<=nmxusedhits) ? true : false;
	      }


      
	      if(fDebug)cout<<"Event No.: "<<ievt<<" Layer: "<<ij<<" "<<yptsall[ij].size()<<" "<<yhits[ij]<<endl;
	      for (int iy=0; iy<yhits[ij]; iy++) {
		if(fDebug)cout<<ypts[ij][iy]<<" ";
	      }
	      if(fDebug)cout<< Ypos[ij]<<" "<<Yusedpos[ij]<<endl;
	    }

      
	    int Ny=0;
	    int nyfail = 0;
	    double ychi2 = 0;
	    double Ydev[nlayer];for (int ij=0; ij<nlayer; ij++) { Ydev[ij] = 100;}
      
	    double yext[nlayer],yexter[nlayer],yposinstr[nlayer];
	    double yslope = 0.;
	    double yinters = 0.;
	    //for (int iy=0; iy<nlayer; iy++) { zval[iy]=iy;}
	    for (int ix=0; ix<nlayer; ix++) { zval[ix]=layerzpos[ix];}

	    for (int iy=0; iy<nlayer; iy++) { yext[iy]= yexter[iy] =yposinstr[iy] =  100;}

	    if(fDebug)cout<<"Fitting occulyr: "<<occulyr<<endl;
	    StraightLineFit yposfit(1,zval,Ypos,yyerr,Yusedpos,occulyr,occulyr,0,11,xyPosDev);
	    yposfit.GetParameters(nyfail,yinters,yslope);
	    yposfit.GetChisqure(Ny,ychi2);
	    yposfit.GetFitValues(yext,Ydev,yexter);
	    GetXposInStrip(yext, yoff, yposinstr);

	    for (int ij=0;ij<nlayer;ij++) {
	      //file_out<<"ypos: "<<Ypos[ij]<<" "<<zval[ij]<<" "<<occulyr<<" "<<yyerr[ij]<<" "<<yext[ij]<<" "<<Ydev[ij]<<" "<<yexter[ij]<<" "<<yposinstr[ij]<<endl;
	    }
	    


	    //file_out<<"Fit Results:: Ny: "<<Ny<<" ychi2: "<<ychi2<<" "<<nyfail<<endl;
	  
	  for (int ij=0; ij<nlayer; ij++) {
	      if(fDebug)cout<<zval[ij]<<" "<<yext[ij]<<" "<<Ypos[ij]<<" "<<Ydev[ij]<<" "<<yext[ij]-Ypos[ij]<<endl;
	    }
	    //Fill Resolutions
      
	    if (Ny>= nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {

	      if(occulyr>=nlayer){
		for (int ij=0; ij<nlayer; ij++) {

		  if (yexter[ij]<0.5) {strp_ymul[ij][iiterrs]->Fill(yposinstr[ij], yhits[ij]);}

		  if (abs(Ydev[ij])<6.0 && Yusedpos[ij]) {
		    ylayer_reso[ij][iiterrs]->Fill(Ydev[ij]);
		    int iyy = min(yhits[ij],nmxhits); 
		    if (iyy>0 && Ypos[ij]>1.0 && Ypos[ij]<nstrip-1 ) { 
		      ylayer_reso_mul[ij][iiterrs][iyy-1]->Fill(Ydev[ij]);
		    }// ylayer_reso_mul[ij][iiterrs][iyy-1]->Fill(Ydev[ij]);
		
		  }
		}
	      }
	      else{
		if (yexter[laye]<0.5) {strp_ymul[laye][iiterrs]->Fill(yposinstr[laye], yhits[laye]);}

		if (abs(Ydev[laye])<6.0 && Yusedpos[laye]) {
		  ylayer_reso[laye][iiterrs]->Fill(Ydev[laye]);
		  int iyy = min(yhits[laye],nmxhits); 
		  if (iyy>0 && Ypos[laye]>1.0 && Ypos[laye]<nstrip-1 ) { 
		    ylayer_reso_mul[laye][iiterrs][iyy-1]->Fill(Ydev[laye]);
		  }// ylayer_reso_mul[laye][iiterrs][iyy-1]->Fill(Ydev[laye]);
		}
	      }

	    }
      

	    //Efficiency and Multiplicity
	    if (Ny>=nmnhits /*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {//4Nov
	      if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {

		if (occulyr>=nlayer) {
		  for (int ij=0; ij<nlayer; ij++) {
		    //cout<<"layer 4 "<<ij<<" "<<xexter[ij]<<" "<<yexter[ij]<<endl;
		    if ( yexter[ij]<0.2 &&  xexter[ij]<0.2) {
        
		      bool  isfiducial = (int(xext[ij]+xoff[ij])>1 && int(xext[ij]+xoff[ij])<31 && int(yext[ij]+yoff[ij])>1 && int(yext[ij]+yoff[ij])<31) ? 1 : 0;
		      isfiducial = true;
		
		      double gap=accprng+effirng; //find_max_gap(ij);	
		
		      if (abs(xposinstr[ij])+xexter[ij]<accprng && abs(yposinstr[ij])+yexter[ij]<accprng) {
			//if (abs(xext[ij])+xexter[ij]<accprng && abs(yext[ij])+yexter[ij]<accprng) {
		    
			totalentry[ij][iiterrs]->Fill(xext[ij],yext[ij]);
		 
			if (isfiducial) {
			  if (xptsall[ij].size()>0) { triggereffi_x[ij][iiterrs]->Fill(xext[ij],yext[ij]);}
			  if (yptsall[ij].size()>0) { triggereffi_y[ij][iiterrs]->Fill(xext[ij],yext[ij]);}
                
			  for (int ix=0; ix<xptsall[ij].size(); ix++) {
			    if (abs(xext[ij]- xptsall[ij][ix] + xoff[ij])<gap) {
			      effi_x[ij][iiterrs]->Fill(xext[ij],yext[ij]);
			      break;
			    }
			  }
                      
			  for (int iy=0; iy<yptsall[ij].size(); iy++) {
			    if (abs(yext[ij]- yptsall[ij][iy] + yoff[ij])<gap) {
			      effi_y[ij][iiterrs]->Fill(xext[ij],yext[ij]);
			      break;
			    }
			  }
		    
		    


			}
		      }
		
		    }// if ( yexter[ij]<0.2 &&  xexter[ij]<0.2) {

		  }//for (int ij=0; ij<nlayer; ij++) {
		}//if (occulyr>=nlayer) {
		else{
		    if ( yexter[occulyr]<0.2 &&  xexter[occulyr]<0.2) {
        
		      bool  isfiducial = (int(xext[occulyr]+xoff[occulyr])>1 && int(xext[occulyr]+xoff[occulyr])<31 && int(yext[occulyr]+yoff[occulyr])>1 && int(yext[occulyr]+yoff[occulyr])<31) ? 1 : 0;
		      isfiducial = true;
		
		      double gap=accprng+effirng; //find_max_gap(occulyr);	
		
		      if (abs(xposinstr[occulyr])+xexter[occulyr]<accprng && abs(yposinstr[occulyr])+yexter[occulyr]<accprng) {
			//if (abs(xext[occulyr])+xexter[occulyr]<accprng && abs(yext[occulyr])+yexter[occulyr]<accprng) {
		    
			totalentry[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
		 
			if (isfiducial) {
			  if (xptsall[occulyr].size()>0) { triggereffi_x[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);}
			  if (yptsall[occulyr].size()>0) { triggereffi_y[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);}
                
			  for (int ix=0; ix<xptsall[occulyr].size(); ix++) {
			    if (abs(xext[occulyr]- xptsall[occulyr][ix] + xoff[occulyr])<gap) {
			      effi_x[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
			      break;
			    }
			  }
                      
			  for (int iy=0; iy<yptsall[occulyr].size(); iy++) {
			    if (abs(yext[occulyr]- yptsall[occulyr][iy] + yoff[occulyr])<gap) {
			      effi_y[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
			      break;
			    }
			  }
		    
		    


			}
		      }
		
		    }// if ( yexter[occulyr]<0.2 &&  xexter[occulyr]<0.2) {

		}
	      }
	    }



      
	  }// for(int ievt=0;ievt<nentry;ievt++) { 
    
	  fileIn->Close();
	  delete event;
	  delete fileIn;  

	}// while(!(file_db.eof())) {



	file_db.close();
	fileOut->cd();

	/////////////////////////////////////
	//
	//  Plotting
	//
	////////////////////////////////////

	if(firstiter){

	  //Plot Occupancy
	  double scale  = 100./max(1, nTotalp);
	  gStyle->SetOptLogy(0);
	  gStyle->SetPadLeftMargin(0.12);
	  gStyle->SetPadBottomMargin(0.14);
	  gStyle->SetOptStat(0);
	  ps.NewPage();
  
	  TCanvas *c1 = new TCanvas ("c1","Occupancy", 500, 700);
	  c1->Divide(3,4);
	  for (int ij=0; ij<nlayer; ij++) {
	    c1->cd(ij+1);
	    occu_x[ij]->Scale(scale);
	    occu_x[ij]->GetXaxis()->SetTitle(occu_x[ij]->GetTitle());
	    occu_x[ij]->GetXaxis()->SetTitleSize(0.07);
	    occu_x[ij]->GetXaxis()->CenterTitle();
	    occu_x[ij]->Draw();
	    occu_x[ij]->SetMarkerStyle(24);
	    occu_x[ij]->SetMarkerSize(0.8);
	  }
	  c1->Update();
          
	  ps.NewPage();
	  for (int ij=0; ij<nlayer; ij++) {
	    c1->cd(ij+1);
	    occu_y[ij]->Scale(scale);
	    occu_y[ij]->GetXaxis()->SetTitle(occu_y[ij]->GetTitle());
	    occu_y[ij]->GetXaxis()->SetTitleSize(0.07);
	    occu_y[ij]->GetXaxis()->CenterTitle();
	    occu_y[ij]->Draw();
	    occu_y[ij]->SetMarkerStyle(24);
	    occu_y[ij]->SetMarkerSize(0.8);

	  }
	  c1->Update();
	  ps.NewPage();
	  if (c1) { delete c1; c1=0;}



  
	  gStyle->SetOptStat(1100);
	  gStyle->SetOptTitle(1);
	  gStyle->SetOptLogy(1);
	  gStyle->SetStatW(.40); //40);
	  gStyle->SetStatH(.20); //30);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetPadLeftMargin(0.12);
	  gStyle->SetPadBottomMargin(0.06);
	  gStyle->SetPadTopMargin(0.09); //(0.03);
	  gStyle->SetPadRightMargin(0.01);

	  TCanvas *c2 = new TCanvas ("c2","Multiplicity", 500, 700);
	  c2->Divide(3,4);
	  gStyle->SetStatY(.99); gStyle->SetStatTextColor(3);//green
	  for (int ij=0; ij<nlayer; ij++) {
	    c2->cd(ij+1);
	    xlayer_mult[ij]->Scale(1./nTotalp);
	    xlayer_mult[ij]->SetLineColor(3); xlayer_mult[ij]->Draw();
	  }
	  c2->Update();
	  gStyle->SetStatY(.79); gStyle->SetStatTextColor(4);//blue
	  for (int ij=0; ij<nlayer; ij++) {
	    c2->cd(ij+1);
	    ylayer_mult[ij]->Scale(1./nTotalp);
	    ylayer_mult[ij]->SetLineColor(4); ylayer_mult[ij]->Draw("sames");
	  }
          
	  c2->Update();
	  ps.NewPage();
	  if (c2) { delete c2; c2=0;}

  
	  firstiter = false;
	}


	ps.NewPage();

	//Fill other Things for all occulyr and iiter

	#ifndef ISALIGN
	
	int istr = (occulyr >= nlayer) ? 0 : occulyr;
	int iend = (occulyr >= nlayer) ? nlayer : occulyr+1;
    
	for (int iocc = istr; iocc<iend; iocc++) {
	  gStyle->SetOptStat(1110);
	  gStyle->SetOptFit(101);
	  gStyle->SetOptLogy(1);
          
	  gStyle->SetPadBottomMargin(0.11);
	  gStyle->SetPadTopMargin(0.08);
	  gStyle->SetPadLeftMargin(0.07);
	  gStyle->SetPadRightMargin(0.02);
	  gStyle->SetOptTitle(1);
	  gStyle->SetTitleFontSize(0.07);
	  gStyle->SetStatW(.20);
	  gStyle->SetStatH(.12);
	  gStyle->SetStatY(.99);
	  gStyle->SetStatX(.99);
          
          
	  TCanvas *c3=new TCanvas ("c3","Residue",500,700);
	  c3->Divide(2,2);
          
	  c3->cd(1);
	  //xlayer_reso[iocc][iiterrs]->Draw("hist");          
    

	  double alw= xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmin();
	  double ahg= xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmax();
	  TF1* fitfx = new TF1("fitfx", gausX, alw, ahg, 3);
          
	  double  parx[3]={xlayer_reso[iocc][iiterrs]->GetMaximum(), xlayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*xlayer_reso[iocc][iiterrs]->GetRMS())};
	  fitfx->SetParameters(parx);
	  //	    fitfx->SetParLimits(0, 0.17*parx[0], 1.5*parx[0]);
	  fitfx->SetParLimits(1, parx[1]-1., parx[1]+1.);
	  fitfx->SetParLimits(2, 0.12, 1.5*parx[2]);
          

	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("X-residues (pitch)");
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
          
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	  xlayer_reso[iocc][iiterrs]->Fit(fitfx, "WW:BMRQ");

	  c3->cd(2);
	  //xlayer_reso[iocc][iiter]->Draw("hist");          
    

	  alw= ylayer_reso[iocc][iiterrs]->GetXaxis()->GetXmin();
	  ahg= ylayer_reso[iocc][iiterrs]->GetXaxis()->GetXmax();
	  TF1* fitfy = new TF1("fitfy", gausX, alw, ahg, 3);
          
	  double  pary[3]={ylayer_reso[iocc][iiterrs]->GetMaximum(), ylayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*ylayer_reso[iocc][iiterrs]->GetRMS())};
	  fitfy->SetParameters(pary);
	  //	    fitfy->SetParLimits(0, 0.17*parx[0], 1.5*parx[0]);
	  fitfy->SetParLimits(1, pary[1]-1., pary[1]+1.);
	  fitfy->SetParLimits(2, 0.12, 1.5*pary[2]);
          

	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("X-residues (pitch)");
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
          
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
	  ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
	  ylayer_reso[iocc][iiterrs]->Fit(fitfy, "WW:BMRQ");
    
  
	  c3->Update();
	  ps.NewPage();
    
	  if (iiter<nmxiter-1) {
	    if (fabs(xlayer_reso[iocc][iiterrs]->GetMean()-
		     fitfx->GetParameter(1)) 
		< xlayer_reso[iocc][iiterrs]->GetRMS()) {
	      xoff[iocc] += fitfx->GetParameter(1);
	    }
            
	    if (fabs(ylayer_reso[iocc][iiterrs]->GetMean()-
		     fitfy->GetParameter(1)) 
		< ylayer_reso[iocc][iiterrs]->GetRMS()) {
	      yoff[iocc] += fitfy->GetParameter(1);
	    }
	  }
          
	  xrms[iocc] = fitfx->GetParameter(2);
	  yrms[iocc] = fitfy->GetParameter(2);	
          
	  file_out <<"lay "<< iocc<<" it "<< iiter <<" X-sh "<< fitfx->GetParameter(1)<<" "<<fitfx->GetParameter(2)<<" "<<xoff[iocc]<<" Y-sh "<< fitfy->GetParameter(1)<<" "<<fitfy->GetParameter(2)<<" "<<yoff[iocc]<<endl;

    
	  delete fitfx; fitfx=0;
	  delete fitfy; fitfy=0;

	  if(c3) {delete c3; c3=0;}


	}// for (int iocc = istr; iocc<iend; iocc++) {

  
  

      file_out<<endl;
      
      


#endif
      }
      }
      }
      
    TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
    c4c->Divide(3,4);
  
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    latex.SetTextSize(0.08);

  
    //Fill other Things for all occulyr and iiter
 
    for (int iiter = 0; iiter<2*nmxiter; iiter++) {

      for (int jkl=0; jkl <2; jkl++) { 
  
	for (int klm =0; klm<nmxhits; klm++) {
	  TH1F* histz[nlayer];
	  ps.NewPage(); 

	  for (int ij=0; ij<nlayer; ij++) {
	    c4c->cd(ij+1);
	    switch(jkl) {
	    case 0 : histz[ij] = (TH1F*)xlayer_reso_mul[ij][iiter][klm]->Clone(); break;
	    case 1 : histz[ij] = (TH1F*)ylayer_reso_mul[ij][iiter][klm]->Clone(); break;
        
	    default : histz[ij] = (TH1F*)ylayer_reso_mul[ij][iiter][klm]->Clone(); break;
	    }
	    histz[ij]->GetXaxis()->SetLabelSize(.07);
	    histz[ij]->GetYaxis()->SetLabelSize(.055);
	    if (jkl<=1) {
	      histz[ij]->GetXaxis()->SetTitle("#Delta X (pitch)");
	    } else {
	      histz[ij]->GetXaxis()->SetTitle("#Delta Y (pitch)");
	    }
	    histz[ij]->GetXaxis()->CenterTitle();
	    histz[ij]->GetXaxis()->SetTitleOffset(0.8);
	    histz[ij]->GetXaxis()->SetTitleSize(0.07); 
          
	    TFitResultPtr ptr = histz[ij]->Fit("gaus","SQ");
	    latex.DrawLatex(0.15, 0.70,Form("%g", int(1000*histz[ij]->GetMean())/1000.));
	    latex.DrawLatex(0.15, 0.56,Form("%g", int(1000*histz[ij]->GetRMS())/1000.));
	    int fitStatus = ptr;
	    if (fitStatus==0 && histz[ij]->GetEntries()>3) {
	      latex.DrawLatex(0.75, 0.70,Form("%g", int(1000*ptr->Parameter(1))/1000.));
	      latex.DrawLatex(0.75, 0.56,Form("%g", int(1000*ptr->Parameter(2))/1000.));
	    }
	  }
	  c4c->Update();
	  for (int ij=0; ij<nlayer; ij++) {
	    if (histz[ij]) { delete histz[ij]; histz[ij]=0;}
	  }
	}
      }
    }


    for (int ij=0; ij<nlayer; ij++) {
      for (int iiter = 0; iiter<2*nmxiter; iiter++) {
      triggereffi_x[ij][iiter]->Divide(totalentry[ij][iiter]);
      triggereffi_y[ij][iiter]->Divide(totalentry[ij][iiter]);
      effi_x[ij][iiter]->Divide(totalentry[ij][iiter]);
      effi_y[ij][iiter]->Divide(totalentry[ij][iiter]);
      }
    }



    gStyle->SetOptLogy(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.11);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetTitleFontSize(0.07);
    latex.SetTextSize(0.08);

    TCanvas* c4d = new TCanvas("c4d", "c4d", 700, 900);
    c4d->Divide(3,4);

   
    
    for (int iiter = 0; iiter<2*nmxiter; iiter++) {

    for (int jkl = 0; jkl < 5; jkl++) {
      TH2F* histy[nlayer];
      ps.NewPage();
  
      for (int ij=0; ij<nlayer; ij++) {
	switch(jkl) {
	case 0 : histy[ij] = (TH2F*)totalentry[ij][iiter]->Clone(); break;
	case 1 : histy[ij] = (TH2F*)triggereffi_x[ij][iiter]->Clone(); break;
	case 2 : histy[ij] = (TH2F*)triggereffi_y[ij][iiter]->Clone(); break;
	case 3 : histy[ij] = (TH2F*)effi_x[ij][iiter]->Clone(); break;
	case 4 : histy[ij] = (TH2F*)effi_y[ij][iiter]->Clone(); break;	
	default: histy[ij] = (TH2F*)effi_y[ij][iiter]->Clone(); break;	

	}// switch(jkl) {
      }// for (int ij=0; ij<nlayer; ij++) {

      for (int ij=0; ij<nlayer; ij++) {
	c4d->cd(ij+1);
      
	//      histy[ij]->SetMinimum(xmn);
	//      histy[ij]->SetMaximum(xmx);
	histy[ij]->GetXaxis()->SetRangeUser(-0.5, 32.5);
	histy[ij]->GetYaxis()->SetRangeUser(-0.5, 32.5);
      
	histy[ij]->GetXaxis()->SetLabelSize(.07);
	histy[ij]->GetYaxis()->SetLabelSize(.07);
	histy[ij]->GetZaxis()->SetLabelSize(.06);
	histy[ij]->GetXaxis()->SetLabelOffset(-.001);
	histy[ij]->GetYaxis()->SetLabelOffset(.01);
      
	histy[ij]->GetXaxis()->SetTitle("X-Strip No");
	histy[ij]->GetXaxis()->CenterTitle();
	histy[ij]->GetXaxis()->SetTitleOffset(0.8);
	histy[ij]->GetXaxis()->SetTitleSize(0.07); 
      
	histy[ij]->GetYaxis()->SetTitle("Y-Strip No");
	histy[ij]->GetYaxis()->CenterTitle();
	histy[ij]->GetYaxis()->SetTitleOffset(0.8);
	histy[ij]->GetYaxis()->SetTitleSize(0.07); 
      
	if(jkl==47 || jkl == 48) {histy[ij]->GetZaxis()->SetRangeUser(0.,1.); }
      
	//      histy[ij]->GetXaxis()->SetTitle(histy[ij]->GetTitle());
	//      histy[ij]->GetXaxis()->SetTitleSize(.10);
	//      histy[ij]->GetXaxis()->SetTitleOffset(.5);
	histy[ij]->SetNdivisions(506,"XY");
	histy[ij]->Draw("colz");
      
	if(jkl > 0){
	  double effi = GetAverage(histy[ij]);
	  //cout<<effi<<endl;
	  latex.DrawLatex(0.5,0.5,Form("%.2f",effi));
	}


      }
      c4d->Update();
  



      for (int ij=0; ij<nlayer; ij++) {
	if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
      }
    
    }//  for (int jkl = 0; jkl < 5; jkl++) {
    }
    



    if(c4d){delete c4d;c4d=0;}







    //Fill Multiplicity
    ps.NewPage();
    gStyle->SetOptLogz(1);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetPadLeftMargin(0.09);
    gStyle->SetPadBottomMargin(0.09);
    gStyle->SetPadTopMargin(0.09); //0.03);
      
    TCanvas* c5a = new TCanvas("c5a", "c5a", 700, 900);
    c5a->Divide(3,4);
  

    for (int iiter = 0; iiter<2*nmxiter; iiter++) {

      
      for(int ixy=0;ixy<2;ixy++){
	
	TH2D* histz[nlayer];
	
	
	for (int ij=0; ij<nlayer; ij++) {
	  switch(ixy){
	  case 0 : histz[ij] = (TH2D*)strp_xmul[ij][iiter]->Clone(); break;
	  case 1 : histz[ij] = (TH2D*)strp_ymul[ij][iiter]->Clone(); break;
	  default : histz[ij] = (TH2D*)strp_xmul[ij][iiter]->Clone(); break;
	    
	  }
	}
	
	
      for (int ij = 0; ij < nlayer; ij++) {

	c5a->cd(ij+1);
	int nbinx = histz[ij]->GetNbinsX();
	int nbiny = histz[ij]->GetNbinsY();

	// Create total array to store sum of Y bins for each X bin
	double total[nbinx] = {0};

	// Calculate total content for each X bin
	for (int ix = 1; ix <= nbinx; ix++) {
	  for (int iy = 1; iy <= nbiny; iy++) {
            total[ix - 1] += histz[ij]->GetBinContent(ix, iy);
	  }
	  if (total[ix - 1] < 1.0) total[ix - 1] = 1.0; // Avoid division by zero
	}

	// Normalize each Y-bin projection

	int icol=0;
	for (int iy = 1; iy <= nbiny; iy++) {
	  icol++; if (icol%5==0) icol++;
	  TH1D* proj_hist = histz[ij]->ProjectionX(Form("proj_hist_layer%d_bin%d", ij, iy), iy, iy);

	  // Normalize bin content
	  for (int ix = 1; ix <= nbinx; ix++) {
            double content = proj_hist->GetBinContent(ix);
            double normalized = content / total[ix - 1];
            double error = sqrt(normalized * (1 - normalized) / total[ix - 1]);
            proj_hist->SetBinContent(ix, normalized);
            proj_hist->SetBinError(ix, error);
	    proj_hist->SetLineColor(icol);
	    proj_hist->SetLineWidth(0);
	    proj_hist->SetMarkerColor(icol);
	    proj_hist->SetMarkerSize(0.9);
	    proj_hist->SetMarkerStyle(kFullCircle);
	    proj_hist->GetYaxis()->SetRangeUser(0., 1.0);
	  }

	  // Draw normalized projection
	  if (iy == 1)
            proj_hist->Draw("P");
	  else
            proj_hist->Draw("P same");
	}
      }


      c5a->Update();
      ps.NewPage();
    
      for (int ij=0; ij<nlayer; ij++) {
	if (histz[ij]) { delete histz[ij]; histz[ij]=0;}
      }
     
    }
  
    }

  
















    ps.Close();

    fileOut->cd();
    fileOut->Write();
    fileOut->Close();

#ifdef ISALIGN    
    //  AlignfileIn->Close();
#else
    AlignTree->Fill();
    AlignfileOut->Write();
    AlignfileOut->Close();
#endif









    return 0;
  }
