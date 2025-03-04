/*
test_basic.log

*/

#define C217STRIP       // Strip convention of C217 for time shift
#define ONEMBYONEM      // RPC of 1m x 1m with readout in 32 strips
//#define OLDFORMAT       //old C217 formated data
#define ONETIMEFORMAT   //Number of TDC per layer (need manual intervention to map strip and TDC)
//#define MONTECARLO      //MC simulated event
//#define FASTSIM     //Simple simulation using aperture_mc.C
//#define FULLG4SIM    //Simulation using Geant4 (C217Sim)
//#define ISEFFICIENCY   
#define TIMESLOPE
//#define MADURAIINTM1 //Madurai intermediate stage, where middle layers are occupied by 50%
//#define MADURAIINTM2  //Madurai intermediate stage, where 4th layer with NINO chip and 5,6,7 & 8 
//      were partially filled from file INORUN_20160421_185949 (thoughOn 15/04/2016, around 12 hrs 
//      SG00013 was moved to stacker to mount the NINO boards ) 
//#define MADURAIINTM3 //With 2mx2m trigger INORUN_20160804_153224
//#define MADURAIINTM4 //With 2mx2m trigger and L1,L4,L5,L6,L7,L8,L9 RPC Daq switch only when combined IDaq and RPCdaq (RPC_INORUN_20170425_185300.ire)
//#define NOISE_SIM // Noise file for simulation 
//#define NEWRPCDAQ // timing data changed in new RPCDaq data in madurai (8 channel per X and Y side)
//#define CAUDATA
//#define NEWRPCDAQ1 //all layers are instumented using NINO-RPCDaq and layer 11 instrumented by Anush-RPCDAQ //Data from RPC_evtraw-r231738.rre
//#define CAUDATA1
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

#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TDatime.h"

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
#include "EveTree.h"
//GMA161203 #include "EveTreeG4.h"
#include "StraightLineFit.h"
#include "TMath.h"
#include <vector>

using namespace std;

const int plot_level=90; // Which figure want to see, more level means less plots


const double cval = 29.979; // velocity of light in cm/ns

#ifdef ONETIMEFORMAT
const int  nTDCpLayer=1; // Number of TDC per layer
#elif defined(NEWRPCDAQ) || defined(NEWRPCDAQ1)
const int  nTDCpLayer=8;
#else 
const int  nTDCpLayer=2;
#endif

const int  nlayer =12; //Maximum RPC layers

const int  nstrip =32; // Strips in RPC chamber
const int  nusedstrip =32; // Strips connected with Amplifier/ADC in RPC chamber
const int  lastXstrip =32; // last strip in X
const int  lastYstrip =32; // last strip in Y

bool  isOnlyCom = false;//true;//false;//true;//false;//true; //Combined fit only, no iteration, not proper efficiency
bool isTimeCorOrReso=true; // 10th Aug 2016, must be true always //false; //true; // true: Time offset Correction and 
                                    // false: Time resolution, they should come one by one 
const int  firstlayer =0; //1; //Used first layer for efficiency calculation
const int  lastlayer =11; // 10; //Used last layer for efficiency calculation

const int  layfirst =0; // 1; //Used first layer in track fitting
const int  laylast =11; // 10; //Used last layer in track fitting

const int  firstXstrip =0;  //1st strip in X
const int  firstYstrip =0;  //1st strip in Y

const double  mxchisq =2.0;//15;//2.0;   //Maximum value of Normalised chi^2 (position fit);
const double  mxtimechisq=2.0; // Maximum value of Normalised chi^2 (time fit);
const double  accprng =1.0; //Acceptance region with respect tot strip centre
const double  effirng =0.25; //Additional region for efficeincy calculation
                           //Extrapolation error < 0.2 of strip width

const int xtcorstr=0;    //Starting layer for X-time correction
const int xtcorend=11;   // End layer for X-time correction
const int ytcorstr=0;    //Starting layer for Y-time correction
const int ytcorend=11;   //End layer for Y-time correction

double find_max_gap(int ix) {
  if (ix >=layfirst && ix <=laylast) return 1.5;
  if (ix <layfirst) return 1.5+0.3*(layfirst-ix);
  if (ix >layfirst) return 1.5+0.3*(ix-laylast);
  return 1.5;
}

const char* labels[nstrip]={"S00", "S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10",
			    "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20",
			    "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", "S30",
			    "S31"};

const char* xlabels[2*nlayer]={"X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11",
			       "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", "Y11"};

const double layerzpos[nlayer]={0.0, 16.25, 32.05, 48.3, 64.2, 80.3, 96.4, 112.1, 128.2, 144.4, 160.3, 176.05};

const double pival=acos(-1);
const int nmxhits=4;
const int nmxusedhits=3;

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

const int nmxtimehit=4;
const int nmxusedtimehit=3;

Double_t gausX(Double_t* x, Double_t* par){
  return par[0]*(TMath::Gaus(x[0], par[1], par[2], kFALSE)); //kTRUE));
}

const int nmxchn = 200;
int nchannel =0;
double m_data[nmxchn];
double m_xpos[nmxchn];

void fcnsg(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t flag) { 
  double fval=0;
  double x[2];
  for (int ij=0; ij<nchannel; ij++) {
    x[0] = m_xpos[ij];
    fval += pow( (m_data[ij] - gausX(x, par)), 2.) / TMath::Max(1., m_data[ij]);
  }
  f = fval;
}

Double_t fitspec(Double_t* x, Double_t* par) {
  //  int nx=int (x[0]/nstrip);
  double yy = x[0]; // - nx*nstrip;
  
  double yval = par[0] + par[1] + par[2]*sin(pival*yy/nusedstrip);
  //  double yval = par[2]*pow(sin(pival*yy/32), par[3]);
  //  double yval = par[2]*sin(pival*yy/32);
  return yval;
}

//bias in time resolution due to uncertainties in other layers
double bias_intime_xreso2[12]={0.404097, 0.294606, 0.250296, 0.168577, 0.128642, 0.114688, 0.120965, 0.138543, 0.158922, 0.251805, 0.316103, 0.46363};
double bias_intime_yreso2[12]={0.508908, 0.390572, 0.342497, 0.227587, 0.167038, 0.142859, 0.140968, 0.173794, 0.230109, 0.359566, 0.454661, 0.635356};

double bias_inpos_xreso2[12]={0.0240879, 0.0176354, 0.0150597, 0.0098424, 0.00719886, 0.00583462, 0.00570846, 0.00690567, 0.00941585, 0.0144282, 0.0174868, 0.0240973};
double bias_inpos_yreso2[12]={0.0242652, 0.0176931, 0.0146815, 0.0100788, 0.00727195, 0.00586189, 0.00567934, 0.00681836, 0.00907941, 0.0137591, 0.0167822, 0.0235308};

//Time shift correction based on postion of position of track in a strip
double parabolic(double x, double* par) {
  double yy = par[0]+par[1]*abs(x)+par[2]*x*x;
  //  cout<<" yyy "<< par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<x <<" "<<fabs(x)<<" "<<abs(x)<<" "<<yy<<endl;
  return yy;
}

double strpos_vs_time[2*nmxtimehit][nlayer][3] = { // 0:1x, 1:2x, 3:3x, 4:4x, 5:1y, 6:2y, 7:3y & 8:4y
{{-0.108298, 0.234881, 0.673403}, 
{-0.0428534, 0.178566, 1.14984}, 
{-0.00179806, -0.058392, 1.48385}, 
{-0.0537578, -0.236439, 2.14217}, 
{-0.0253757, -0.187405, 3.09114}, 
{-0.0696471, -0.183559, 3.78667}, 
{0.00250928, -0.276274, 2.84852}, 
{0.00835515, -0.302841, 2.87873}, 
{-0.167504, 0.238886, 1.81847}, 
{-0.0788643, -0.159428, 1.91216}, 
{-0.0894312, 0.0856115, 1.11812}, 
{-0.0779776, 0.280962, 0.816331}}, 

{{-0.214927, 1.42247, -1.16535}, 
{-0.529811, 2.10891, -1.24669}, 
{-0.334925, 1.51751, -0.815884}, 
{-0.23738, 0.973092, 0.138824}, 
{-0.531924, 2.05702, -0.448773}, 
{-0.856707, 2.48331, 0.00234152}, 
{-0.336572, 0.837857, 0.942916}, 
{-0.34803, 1.08262, 0.572813}, 
{-0.489068, 1.92972, -0.656595}, 
{-0.30959, 1.66386, -0.946176}, 
{-0.20136, 1.233, -0.594369}, 
{-0.320393, 1.41313, -0.852566}}, 

 {{-1.13823, 0.415615, -0.0956477}, 
{-0.938263, 0.444087, -0.333743}, 
{-0.691977, 0.240515, -0.2864}, 
{-0.543465, 0.348418, -0.170337}, 
{-0.97317, 0.539603, -0.902608}, 
{-1.39851, 0.387069, -0.574954}, 
{-0.521833, 0.18826, 0.461717}, 
{-0.582583, 0.382221, -0.114176}, 
{-1.00909, -0.0518222, -0.00510963}, 
{-0.798994, 0.286262, -0.0887585}, 
{-0.551053, 0.278699, -0.448895}, 
{-0.607885, 0.123587, -0.00279212}}, 

{{-1.00401, 1.8883, -2.63737}, 
{-0.758471, 0.544713, -0.793395}, 
{-0.685073, -0.094587, -0.216575}, 
{-0.617107, 0.49322, -1.20662}, 
{-0.946502, 0.632696, -0.874812}, 
{-1.19493, 1.00998, -1.33925}, 
{-0.69642, 0.0394355, -0.33321}, 
{-0.640847, 0.168112, -0.245247}, 
{-0.670089, 0.702369, -1.27111}, 
{-0.736155, 0.147291, 0.0974525}, 
{-0.692997, 0.0315732, -0.6247}, 
{-0.607172, -0.669918, 0.388367}},

{{-0.00329939, 0.242462, 1.06257}, 
{-0.0902285, 0.0624302, 2.03717}, 
{-0.0578305, 0.0933354, 1.96958}, 
{-0.0585525, 0.401106, 1.45776}, 
{-0.0667166, -0.319948, 3.5735}, 
{-0.120356, -0.369183, 4.42829}, 
{-0.298342, 0.70353, 1.81377}, 
{-0.150394, -0.433648, 4.82061}, 
{-0.102069, 0.0859978, 2.08871}, 
{-0.124022, 0.0217614, 1.94001}, 
{-0.0913154, -0.114765, 2.79464}, 
{-0.114681, 0.102416, 1.63289}}, 

{{-0.272125, 1.68014, -1.17568}, 
{-0.487944, 2.29154, -1.52053}, 
{-0.323852, 1.38964, -0.240395}, 
{-0.34994, 1.6006, -0.671717}, 
{-0.718347, 2.41741, -0.471764}, 
{-0.721208, 2.50796, -0.40962}, 
{-0.267072, 1.43869, -0.586951}, 
{-0.649484, 2.28938, -0.168097}, 
{-0.95087, 3.55793, -1.35827}, 
{-0.881629, 3.59384, -2.48311}, 
{-0.435896, 1.95779, -1.00669}, 
 {-0.390051, 1.95935, -1.36556}},

{{-1.18847, 0.0757272, 0.117748}, 
{-1.19774, 0.267215, 0.0940863}, 
{-0.844947, 0.544917, -0.412739}, 
{-0.731456, 0.455506, -0.742424}, 
{-1.49598, 0.688469, 0.132483}, 
{-1.38607, 0.569129, -0.348203}, 
{-0.514269, 0.568307, -1.28993}, 
{-1.04816, 0.668811, -0.60223}, 
{-1.29873, 0.645937, -0.852456}, 
{-1.51009, 0.171537, 0.353457}, 
{-0.923922, 0.40255, -0.626108}, 
{-1.04374, 0.076331, -0.0552816}},

{{-0.92741, -0.989598, 1.67875}, 
{-0.670879, -0.485287, 1.01524}, 
{-0.791233, 0.0619398, -0.343801}, 
{-0.781125, -0.409228, -0.243095}, 
{-0.975324, 0.466529, -0.269322}, 
{-1.04694, 0.801549, -0.997753}, 
{-0.556878, 0.961198, -2.14108}, 
{-0.928105, 1.69794, -2.72732}, 
{-1.01487, 0.592797, 0.340055}, 
{-1.00587, -0.0115691, 1.0955}, 
{-0.855476, 0.00859748, -0.137312},
{-0.803551, 0.247394, -1.63056}}
};

double biasinxtime[nlayer]={0};

double biasinytime[nlayer]={0};

const double stripwidth = 3.0; // cm
const double stripgap = 0.2; // cm
const double layergap = 16.0; // cm
const double fact = stripwidth; //layergap;
const double sigmaz = 0.1; //in cm as the smallest division in scale
//  const double sigmapos = 0.8; //in cm as 2.8/sqrt(12). its a uniform distribution
const double sigmarpc = 1.5; //in ns RPC time resolution
const int nmnhits =4;//5;//3;// 4;//3;//4;//9;//5;//4; //9; //Minimum number of layers required for position and time measurements
const int nmnentry = 10; // Statistics requred in a strip for time correction and use in next iteration

bool isTiming = false; //true;//false;//true;//false;//true;//true; // true //Don't do timing fit etc
bool useallpos=false; //true; //While use all strip irrespecive of aligned or not
bool usealltime=false; //While use timing of all strip irrespecive of aligned or not
bool timefit;

const float xyPosDev=7.0; // seven sigma 2.0; //maximum deviation of points from fit line

const int trigly1 =2; //0;
const int trigly2 =4;
const int trigly3 =7;
const int trigly4 =9; //11;


const int nDelayPar=6;
double timecorpar[2][nlayer][nstrip][nDelayPar]={0.0};

Double_t polfunc(Double_t* x, Double_t* par) {
  return par[0]+
    par[1]*x[0]+
    par[2]*x[0]*x[0]+
    par[3]*x[0]*x[0]*x[0]+
    par[4]*x[0]*x[0]*x[0]*x[0];
}

double path_correction_time(int ixy, int il, int ist, double xx) {
  double xval[2];
  double par[nDelayPar];
  xval[0] = xx;
 
  for (int ij=0; ij<nDelayPar; ij++) {
    par[ij] = timecorpar[ixy][il][ist][ij];
  }

  double cor= polfunc(xval, par) - par[nDelayPar-1]; //compared wrt central one

  if (cor > 5.0) cor =  5.0;
  if (cor < -5.0) cor = -5.0; 
          
  return cor;
}

void GetXposInStrip(double* ext, double* off, double* pos) {
  for (int ij=0; ij<nlayer; ij++) { 
    int istr = int(ext[ij]+off[ij]);
    if (ext[ij]+off[ij]<0.0) {
      pos[ij] = ext[ij]+off[ij] - istr + 0.5;
    } else {
      pos[ij] = ext[ij]+off[ij] - istr - 0.5;
    }
  }
}

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -1;
  for (int ix=1; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return 1000;
}

int main() {
  // ----------------- Don't Change Anything Starting from Here ------------       ********************   ---------------------------
  static unsigned int mypow_2[32];
  for (int ij=0; ij<32; ij++) {
    mypow_2[ij] = pow(2, ij);
  }
  
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
  gStyle->SetTitleYOffset(0.05);

  gStyle->SetStatFont(22);        // Times New Roman
  gStyle->SetTextFont(22);        // Times New Roman
  gStyle->SetTitleFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelFont(22,"XYZ"); // Times New Roman
  gStyle->SetLabelSize(0.06, "XYZ"); // Times New Roman  
  gStyle->SetNdivisions(606, "XYZ");
  gStyle->SetPaintTextFormat("6.4f"); 

  //  gStyle->SetStatFontSize(.015);

  //  gStyle->SetPadGridX(3);
  //  gStyle->SetPadGridY(3);
  //  gStyle->SetGridStyle(3);

  gStyle->SetOptTitle(0);
  gStyle->SetFuncWidth(1);
  gStyle->SetFuncColor(2);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(101);
  gStyle->SetOptLogy(0);
  gStyle->SetStatW(.18);
  gStyle->SetStatH(.08);
  gStyle->SetPadTopMargin(.001); //0.09
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.001);
  gStyle->SetPadRightMargin(0.001);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.12);
  latex.SetTextFont(42);
  latex.SetTextAlign(1); //(31); // align right

#ifdef MONTECARLO
  TH2F* tmp_sel_theta_phi[10];
  TH2F* sel_theta_phi[10]; //Added to store theta phi distributions for different sets, mainly to obtain acceptance (+selection) efficiency to compare data and MC
  TH2F* tmp_sel_theta_phi_8[10][8];
  TH2F* sel_theta_phi_8[10][8];
  TH2F* sel_theta_phi_16[10][16];
  TH2F* tmp_sel_theta_phi_16[10][16];
#ifdef FASTSIM
  float thgen, phgen, xpgen, ypgen;
  float xtrue[nlayer];
  float ytrue[nlayer];
  float xposa[nlayer];
  float yposa[nlayer];
  int   xystrp[nlayer];
  int   xstrp[nlayer][3];
  int   ystrp[nlayer][3];
  float xtimea[nlayer];
  float ytimea[nlayer];

  //TH2F* tmp_sel_theta_phi[10];
  //TH2F* sel_theta_phi[10]; //Added to store theta phi distributions for different sets, mainly to obtain acceptance (+selection) efficiency to compare data and MC
  
#else //FULLG4SIM
  
  static const unsigned int ngenmx=50;
  static const unsigned int nlayermx = 12;
  static const unsigned int nsimhtmx= 1000;
  
  UInt_t    irun;                // Run number of these events
  UInt_t    ievt;                //Event number
  UInt_t    ngent;
  Float_t   ievt_wt;                //*GMa
  Int_t      intxn_id;               //*GMa
  
  Int_t   pidin[ngenmx];    //PID of incident particle
  Float_t  momin[ngenmx];    //Energy of incident particle
  Float_t  thein[ngenmx];    //Initial polar angle of incident particle
  Float_t  thegen[ngenmx];    //Initial polar angle of incident particle
  Float_t  phiin[ngenmx];     //Initial azimuthal angle of incident particle 
  Float_t  posxin[ngenmx];   //Initial X-position
  Float_t  posyin[ngenmx];     //Initial Y-position
  Float_t  poszin[ngenmx];     //Initial Z-position
  
  // For Simulation output
  UInt_t nsimhtx;  // Number of Hits in X side
  UInt_t nsimhty;  // Number of Hits in Y side
  
  Int_t nlayert; // Number of Layers
  Int_t nstript; // Number of Layers
  Int_t digiXtime[nlayermx]; // Time for each layer. in multiples of 100 ps. for each layer. Time stamp using the earliest signal. X side.
  Int_t digiYtime[nlayermx]; // Time for each layer. in multiples of 100 ps. for each layer. Time stamp using the earliest signal. Y side.
  ULong64_t xdata[nlayermx]; // 32 strips data bit wise for 12 layers. X side.
  ULong64_t ydata[nlayermx]; // 32 strips data bit wise for 12 layers. X side.
  UInt_t triggerinfoX;
  UInt_t triggerinfoY;  // Number of Hits in X side
#endif
#endif

  int EvtIndx,yyindx,xxindx,BID;

  bool posinuse[8][nlayer][nstrip]; //Reject strips in timing, which are not aligned
  bool timeinuse[8][nlayer][nstrip]; //Reject strips in timing, which are not aligned
  
  for ( int ij=0; ij<8; ij++) {
    for ( int jk=0; jk<nlayer; jk++) {
      for ( int kl=0; kl<nstrip; kl++) {
        posinuse[ij][jk][kl] = timeinuse[ij][jk][kl]=true;
      }
    }
  }
  // To calcualte average X-time and Y-time to synchronise them.
  int ntotxtime=0;
  double totxtime=0;
  int ntotytime=0;
  double totytime=0;
  
  int isequence[nlayer]={0,1,2,3,4,5,6,7,8,9,10,11};
  int nentrymx=-1;
  int ntotal = 0;
  int nseltot=0;
  
  UInt_t nsec=0;
  int isalign=0;
  
  char outfile[100];
  char outfilx[100];
  
  char name[100];
  char title[100];
  char infile[200];
  char datafile[100];
  char rootfiles[100];
  //  cout <<"Give the input file name"<<endl;;
  //  cin>> rootfiles;
  //  sprintf(rootfiles, "test_2479.log");
  sprintf(rootfiles, "test_basic.log");
  if (isOnlyCom) { 
    isalign = 0;
  } else {
    //    cout <<"Do you want to align ? yes/no 1/0" <<endl;;
    //    cout <<"for alignment, it can be any number greater than zero"<<endl;
    //    cin>> isalign;
    isalign =5;
  }
  
  //# of iteration where all layers are included
  // For the time being it is not implemented, but can be used
  // Can be use in other way, first few iteration with combined+individual, then
  // Only combined
  //Total # of iteration all layers + individual layers
//#ifdef MONTECARLO
 // const int nmxiter = 1; 
//#else
  const int nmxiter = (isalign>0) ? 1 : 1;// 0;//1; //Less than 12, otherwise change the pad in canvas
//#endif
  const int nlayerit = nlayer; // (isalign >0) ? nlayer : 1;
  
  int /*ievt,*/ nhits;
  double xslope, xinters, yslope, yinters, timexslope, timeyslope, timex2slope, timey2slope;
  float txslop, tyslop, timexinters, timeyinters, errtimexslope,errtimeyslope,errtimexinters, errtimeyinters;
  double xchi2, ychi2, xt0chi2, yt0chi2;
  int nxstrip, nystrip,Nx,Ny,nxtime,nytime, ntxyla;
  double zen;
  double xtexter[nlayer];
  double ytexter[nlayer]; 

  double Xpos[nlayer], Xdev[nlayer]; 
  bool Xusedpos[nlayer];
  double Ypos[nlayer], Ydev[nlayer];
  bool Yusedpos[nlayer];//=new float[nlayer];

  double Xpos1[nlayer], Xdev1[nlayer]; 
  bool Xusedpos1[nlayer];
  double Ypos1[nlayer], Ydev1[nlayer];
  bool Yusedpos1[nlayer];//=new float[nlayer];

  double Xpos2[nlayer], Xdev2[nlayer]; 
  bool Xusedpos2[nlayer];
  double Ypos2[nlayer], Ydev2[nlayer];
  bool Yusedpos2[nlayer];//=new float[nlayer];

  double xslope1, xinters1, yslope1, yinters1;
  double xslope2, xinters2, yslope2, yinters2;
  double xchi21, ychi21,xchi22, ychi22;
  int Nx1,Ny1,Nx2,Ny2;

  int len = strlen(rootfiles);
  strncpy(outfilx, rootfiles, len-4);
  outfilx[len-4]='\0';
  sprintf(outfilx, "%s_jj0_", outfilx);
  len = strlen(outfilx);
  outfilx[len]='\0';

#if defined(NOISE_SIM)
  //#ifndef MONTECARLO

  int ntot_uncsim=0;
  sprintf(outfile, "%snoise_%i.root", outfilx, isalign);
  TFile* NoisefileOut = new TFile(outfile, "recreate");	      
  
  TH1F* n_multall[2];
  
  TH1F* layer_hitsall[2][nlayer];
  TH1F* layer_multall[2][nlayer];
  TH1F* layer_timeall[2][nlayer];
  
  TH2F* total_vs_indmulall[2][nlayer];
  
  TH2F* mult_correlall[2][nlayer];
  TH2F* hist_correlall[2][nlayer];
  
  //After removal of muon hits
  TH1F* n_mult[2];
  TH1F* layer_hits[2][nlayer];
  TH1F* layer_mult[2][nlayer];
  TH1F* layer_time[2][nlayer];
  
  TH2F* total_vs_indmul[2][nlayer];
  
  TH2F* mult_correl[2][nlayer];
  TH2F* hist_correl[2][nlayer];
  
  for (int ixy=0; ixy<2; ixy++) { 
    sprintf(title, "%s", (ixy==0) ? "x" : "y");
    sprintf(name, "%sn_mult", title);
    n_mult[ixy] = new TH1F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5);
    sprintf(name, "%sn_multall", title);
    n_multall[ixy] = new TH1F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5);
    
    for (int ij = 0; ij<nlayer; ij++) {
      sprintf(name, "%slayer_hitsall_%i", title, ij);
      layer_hitsall[ixy][ij] = new TH1F(name, name, nstrip, -0.5, nstrip-0.5);
      
      sprintf(name, "%slayer_multall_%i", title, ij);
      layer_multall[ixy][ij] = new TH1F(name, name, nstrip+1, -0.5, nstrip+0.5);
      
      sprintf(name, "%slayer_timeall_%i", title, ij);
      layer_timeall[ixy][ij] = new TH1F(name, name, 200, 0.0, 200.);
      
      sprintf(name, "%shist_correlall_l%i", title, ij);
      hist_correlall[ixy][ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(name, "%stotal_vs_indmulall_l%i", title, ij);
      total_vs_indmulall[ixy][ij] = new TH2F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5, nstrip+1, -0.5, nstrip+0.5); 
      
      sprintf(name, "%smult_correlall_l%i", title, ij);
      mult_correlall[ixy][ij] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5, nstrip+1, -0.5, nstrip+0.5); 
      
      sprintf(name, "%slayer_hits_%i", title, ij);
      layer_hits[ixy][ij] = new TH1F(name, name, nstrip, -0.5, nstrip-0.5);
      
      sprintf(name, "%slayer_mult_%i", title, ij);
      layer_mult[ixy][ij] = new TH1F(name, name, nstrip+1, -0.5, nstrip+0.5);
      
      sprintf(name, "%slayer_time_%i", title, ij);
      layer_time[ixy][ij] = new TH1F(name, name, 200, 0.0, 200.);
      
      sprintf(name, "%shist_correl_l%i", title, ij);
      hist_correl[ixy][ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(name, "%stotal_vs_indmul_l%i", title, ij);
      total_vs_indmul[ixy][ij] = new TH2F(name, name, nlayer*nstrip+1, -0.5, nlayer*nstrip+0.5, nstrip+1, -0.5, nstrip+0.5); 
      
      sprintf(name, "%smult_correl_l%i", title, ij);
      mult_correl[ixy][ij] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5, nstrip+1, -0.5, nstrip+0.5); 
    }
  }
  //#endif
#endif

  sprintf(outfile, "%s%i.root", outfilx, isalign);
  TFile* fileOut = new TFile(outfile, "recreate");
  
  //  sprintf(outfile, "cor_%s%i.root", outfilx, isalign);
  // TFile* filecorOut = new TFile(outfile, "recreate");
  Int_t EvNum;
  UShort_t XYhitall[12];
  // UShort_t Yhitall[12];
  ULong_t XYhitcor;
  // ULong_t Yhitcor[12];
  
  sprintf(outfile, "%s%i.ps", outfilx, isalign);
  TPostScript ps(outfile,111);  
  ps.Range(20,30); //ps.Range(10,20);
  
  sprintf(outfile, "%s%i.txt", outfilx, isalign);
  ofstream file_out(outfile);

  sprintf(outfile, "%s%i_str.txt", outfilx, isalign);
  ofstream file_outstr(outfile);

  TTree* T2 = new TTree("T2", "store"); //("Tree Name","Tree Title")

  //  T2->Branch("Evt",&ievt,"ievt/I");
  T2->Branch("nhits", &nhits, "nhits/I");
  T2->Branch("xslope", &xslope, "xslope/D");
  T2->Branch("xinters", &xinters, "xinters/D");
  T2->Branch("yslope", &yslope, "yslope/D");
  T2->Branch("yinters", &yinters, "yinters/D");

  T2->Branch("xchi2", &xchi2, "xchi2/D"); 
  T2->Branch("ychi2", &ychi2, "ychi2/D"); 
  T2->Branch("zen",&zen,"zen/D");


  TTree* T3 = new TTree("T3","store_T3");
  T3->Branch("EvNum",&EvNum,"EvNum/i");
  T3->Branch("XYhitall",XYhitall,"XYhitall[12]/s");
  T3->Branch("XYhitcor",&XYhitcor,"XYhitcor/L");

  if (isTiming) {
    
    T2->Branch("txslop", &txslop, "txslop/F");
    T2->Branch("tyslop", &tyslop, "tyslop/F");
    T2->Branch("xt0chi2", &xt0chi2, "xt0chi2/D");
    T2->Branch("yt0chi2", &yt0chi2, "yt0chi2/D");
    T2->Branch("nxtime", &nxtime, "nxtime/I");
    T2->Branch("nytime", &nytime, "nytime/I");
    T2->Branch("ntxyla", &ntxyla, "ntxyla/I");
    T2->Branch("timexslope",&timexslope,"timexslope/D");
    T2->Branch("timeyslope",&timeyslope,"timeyslope/D");
    T2->Branch("errtimexslope",&errtimexslope,"errtimexslope/F");
    T2->Branch("xtexter", xtexter, "xtexter[12]/D");
    T2->Branch("ytexter", ytexter, "ytexter[12]/D");
  }


#ifdef MONTECARLO
#ifdef FASTSIM
  
  T2->Branch("thgen", &thgen, "thgen/F");
  T2->Branch("phgen", &phgen, "phgen/F"); 
  T2->Branch("xpgen", &xpgen, "xpgen/F");
  T2->Branch("ypgen", &ypgen, "ypgen/F");  
  
  T2->Branch("xystrp", xystrp, "xystrp[12]/I");
  T2->Branch("xstrp", xstrp, "xstrp[12][3]/I");
  T2->Branch("ystrp", ystrp, "ystrp[12][3]/I");
  T2->Branch("xtimea", xtimea, "xtimea[12]/F");  
  T2->Branch("ytimea", ytimea, "ytimea[12]/F"); 

#else  //FULLG4SIM
  T2->Branch("irun",&irun,"irun/i"); //VALGRIND
  T2->Branch("ievt",&ievt,"ievt/i");
  
  T2->Branch("ngent",&ngent,"ngent/i");
  T2->Branch("pidin",pidin,"pidin[ngent]/I");
  T2->Branch("ievt_wt",&ievt_wt,"ievt_wt/F");
  T2->Branch("intxn_id",&intxn_id,"intxn_id/I");
  T2->Branch("momin",momin,"momin[ngent]/F");
  T2->Branch("thein",thein,"thein[ngent]/F");
  T2->Branch("phiin",phiin,"phiin[ngent]/F");
  T2->Branch("posxin",posxin,"posxin[ngent]/F");
  T2->Branch("posyin",posyin,"posyin[ngent]/F");
  T2->Branch("poszin",poszin,"poszin[ngent]/F");
  T2->Branch("nsimhtx", &nsimhtx, "nsimhtx/i");
  T2->Branch("nsimhty", &nsimhty, "nsimhty/i");
  T2->Branch("nlayert", &nlayert, "nlayert/I");
  T2->Branch("xdata",xdata,"xdata[nlayert]/L");
  T2->Branch("ydata",ydata,"ydata[nlayert]/L");
  T2->Branch("triggerinfoX",&triggerinfoX,"triggerinfoX/i");
  T2->Branch("triggerinfoY",&triggerinfoY,"triggerinfoY/i");

  if (isTiming) {
    T2->Branch("digiXtime", digiXtime, "digiXtime[nlayert]/I");
    T2->Branch("digiYtime", digiYtime, "digiYtime[nlayert]/I");
  }
#endif
#endif

  int narray =75;//60;
  const int nselcrit=16;//6;
  TH1F* costhe[nselcrit];
  costhe[0] = new TH1F("costhe_all", "costhe_all", narray, 0., narray);//filling as dn/d(theta)
  costhe[1] = new TH1F("costhe_trig", "costhe_trig", narray, 0., narray);
  costhe[2] = new TH1F("costhe_selecx", "costhe_selecx", narray, 0., narray);
  costhe[3] = new TH1F("costhe_selecy", "costhe_selecy", narray, 0., narray);
  costhe[4] = new TH1F("costhe_accep", "costhe_accep", narray, 0., narray);
  costhe[5] = new TH1F("costhe_accep_H", "costhe_accep_H", 10, 0.5, 1.0);//10Nov Honda also give 0 to 60 deg in same bin filling as dn/d(costheta)

  // To find the multiple scattering in Layer 5 due to lead block
  const int npixel_theta = 16;
  TH2F* pixel_diff_theta[16];; 
  TH1F* theta_diff_1 = new TH1F("theta_diff_1","theta_diff_1",100.,0.,10.);
  TH1F* theta_diff_2 = new TH1F("theta_diff_2","theta_diff_2",100.,0.,10.);
  
  for(int ij=0;ij<npixel_theta;ij++) {
    sprintf(name,"diff_theta_L5_i%i",ij);
    pixel_diff_theta[ij] = new TH2F(name,name,stripwidth*nstrip/2,0.,96.,stripwidth*nstrip/2,0.,96.);
  }
  
  TH1F* pixel_scatang1[nstrip][nstrip] = {{0}};
  TH1F* pixel_scatang[nstrip][nstrip] = {{0}};
  double val_scatang[nstrip][nstrip]={{0.}};
  for(int ij =0;ij< nstrip ;ij++) {
    for(int jk=0;jk<nstrip;jk++) {
      sprintf(name,"pixel_scatang1_x%i_y%i",ij,jk);
      pixel_scatang1[ij][jk] = new TH1F(name,name,100,0.,10.);
      sprintf(name,"pixel_scatang_x%i_y%i",ij,jk);
      pixel_scatang[ij][jk] = new TH1F(name,name,100,0.,10.);
    }
  }
  
  TH2F* pixel_scatmean1= new TH2F("scatang_mean1","scatang_mean1",32,-0.5,31.5,32,-0.5,31.5);
  TH2F* pixel_scatmean= new TH2F("scatang_mean","scatang_mean",32,-0.5,31.5,32,-0.5,31.5);
  // 9th July 2016 : Added to test "n" for different criteria
  const int ntrkselcrit=140;
  TH1F* zenithang[nselcrit][ntrkselcrit];
  TH1F* azimuthang[nselcrit][ntrkselcrit];
  TH1F* zenithang_azimuth[nselcrit][ntrkselcrit];
  TH1F* zenithang_azimuth_8[nselcrit][ntrkselcrit];
  for (int ij=0; ij<nselcrit; ij++) {
    for (int jk=0; jk<ntrkselcrit; jk++) { 
      sprintf(name, "zenithang_%i_%i", ij, jk);
      zenithang[ij][jk] = new TH1F(name, name, 2*narray, 0., narray);
      sprintf(name, "azimuthang_%i_%i", ij, jk);
      azimuthang[ij][jk] = new TH1F(name, name, 36, 0.0, 360.0);
      sprintf(name, "zenithang_azimuth_%i_%i", ij, jk);
      zenithang_azimuth[ij][jk] = new TH1F(name, name, 2*narray, 0., narray);
      sprintf(name, "zenithang_azimuth_8_%i_%i", ij, jk);
      zenithang_azimuth_8[ij][jk] = new TH1F(name, name, 2*narray, 0., narray);
    }
  }
  
  TH1F* phiang[6];
  
  phiang[0] = new TH1F("phiang_all", "phiang_all", 36, -pival, pival);
  phiang[1] = new TH1F("phiang_trig", "phiang_trig", 36, -pival, pival);
  phiang[2] = new TH1F("phiang_selecx", "phiang_selecx", 36, -pival, pival);
  phiang[3] = new TH1F("phiang_selecy", "phiang_selecy", 36, -pival, pival);
  phiang[4] = new TH1F("phiang_accep", "phiang_accep", 36, -pival, pival);

  TH2D* strp_xmul[nlayer][2*nmxiter];
  TH2D* strp_ymul[nlayer][2*nmxiter];

  TH2D* strp_xmulsim[nlayer][2*nmxiter];
  TH2D* strp_ymulsim[nlayer][2*nmxiter];
  

#ifdef ISEFFICIENCY

  TH2D* defefficiency_uncx[nlayer];
  TH2D* defefficiency_uncy[nlayer];

  TH2D* deftriggereffi_x[nlayer];
  TH2D* deftriggereffi_y[nlayer];

  TH1D* inefficiency_xallpixel[nlayer][2*nmxiter];
  TH1D* inefficiency_yallpixel[nlayer][2*nmxiter];

  TH1D* inefficiency_xpixel[nlayer][2*nmxiter];
  TH1D* inefficiency_ypixel[nlayer][2*nmxiter];

  TH2D* inefficiency_uncx[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy[nlayer][2*nmxiter];

  TH2D* inefficiency_corx[nlayer][2*nmxiter];
  TH2D* inefficiency_cory[nlayer][2*nmxiter];


  TH2D* inefficiency_xt[nlayer][2*nmxiter];
  TH2D* inefficiency_yt[nlayer][2*nmxiter];

  TH2D* total_xt[nlayer][2*nmxiter];
  TH2D* total_yt[nlayer][2*nmxiter];

  TH2D* triggereffi_x[nlayer][2*nmxiter];
  TH2D* triggereffi_y[nlayer][2*nmxiter];

  TH2D* triggereffi_xevt[nlayer][2*nmxiter];
  TH2D* triggereffi_yevt[nlayer][2*nmxiter];

  TH2D* totalentry[nlayer][2*nmxiter];

  TH2D* difefficiency_uncx[nlayer][2*nmxiter];
  TH2D* difefficiency_uncy[nlayer][2*nmxiter];

  TH2D* difefficiency_xt[nlayer][2*nmxiter];
  TH2D* difefficiency_yt[nlayer][2*nmxiter];

  TH2D* diftriggereffi_x[nlayer][2*nmxiter];
  TH2D* diftriggereffi_y[nlayer][2*nmxiter];

  TH2D* diftriggereffi_xevt[nlayer][2*nmxiter];
  TH2D* diftriggereffi_yevt[nlayer][2*nmxiter];

  //corr efficiencies for Monte Carlo
  
  TH2D* totalentry_5[nlayer][2*nmxiter];
  TH2D* totalentry_6[nlayer][2*nmxiter];
  TH2D* totalentry_7[nlayer][2*nmxiter];
  TH2D* totalentry_8[nlayer][2*nmxiter];
  TH2D* totalentry_9[nlayer][2*nmxiter];
  TH2D* totalentry_10[nlayer][2*nmxiter];
  TH2D* totalentry_11[nlayer][2*nmxiter];
  
  TH2D* triggereffi_x_5[nlayer][2*nmxiter];
  TH2D* triggereffi_y_5[nlayer][2*nmxiter];
  TH2D* triggereffi_x_6[nlayer][2*nmxiter];
  TH2D* triggereffi_y_6[nlayer][2*nmxiter];
  TH2D* triggereffi_x_7[nlayer][2*nmxiter];
  TH2D* triggereffi_y_7[nlayer][2*nmxiter];
  TH2D* triggereffi_x_8[nlayer][2*nmxiter];
  TH2D* triggereffi_y_8[nlayer][2*nmxiter];
  TH2D* triggereffi_x_9[nlayer][2*nmxiter];
  TH2D* triggereffi_y_9[nlayer][2*nmxiter];
  TH2D* triggereffi_x_10[nlayer][2*nmxiter];
  TH2D* triggereffi_y_10[nlayer][2*nmxiter];
  TH2D* triggereffi_x_11[nlayer][2*nmxiter];
  TH2D* triggereffi_y_11[nlayer][2*nmxiter];

  TH2D* inefficiency_uncx_5[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy_5[nlayer][2*nmxiter];
  TH2D* inefficiency_uncx_6[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy_6[nlayer][2*nmxiter];
  TH2D* inefficiency_uncx_7[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy_7[nlayer][2*nmxiter];
  TH2D* inefficiency_uncx_8[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy_8[nlayer][2*nmxiter];
  TH2D* inefficiency_uncx_9[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy_9[nlayer][2*nmxiter];
  TH2D* inefficiency_uncx_10[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy_10[nlayer][2*nmxiter];
  TH2D* inefficiency_uncx_11[nlayer][2*nmxiter];
  TH2D* inefficiency_uncy_11[nlayer][2*nmxiter];
  
  TH2D* inefficiency_corx_5[nlayer][2*nmxiter];
  TH2D* inefficiency_cory_5[nlayer][2*nmxiter];
  TH2D* inefficiency_corx_6[nlayer][2*nmxiter];
  TH2D* inefficiency_cory_6[nlayer][2*nmxiter];
  TH2D* inefficiency_corx_7[nlayer][2*nmxiter];
  TH2D* inefficiency_cory_7[nlayer][2*nmxiter];
  TH2D* inefficiency_corx_8[nlayer][2*nmxiter];
  TH2D* inefficiency_cory_8[nlayer][2*nmxiter];
  TH2D* inefficiency_corx_9[nlayer][2*nmxiter];
  TH2D* inefficiency_cory_9[nlayer][2*nmxiter];
  TH2D* inefficiency_corx_10[nlayer][2*nmxiter];
  TH2D* inefficiency_cory_10[nlayer][2*nmxiter];
  TH2D* inefficiency_corx_11[nlayer][2*nmxiter];
  TH2D* inefficiency_cory_11[nlayer][2*nmxiter];
#endif

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "strp_xmul_l%i_i%i", ij, jk);
      strp_xmul[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);
      
      sprintf(title, "strp_ymul_l%i_i%i", ij, jk);
      strp_ymul[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);

      sprintf(title, "strp_xmulsim_l%i_i%i", ij, jk);
      strp_xmulsim[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);
      
      sprintf(title, "strp_ymulsim_l%i_i%i", ij, jk);
      strp_ymulsim[ij][jk]=new TH2D(title, title, 60, -0.5, 0.5, 5, -0.5, 4.5);
    }
    
#ifdef ISEFFICIENCY

    sprintf(title, "defefficiency_uncx_l%i", ij);
    defefficiency_uncx[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
    sprintf(title, "defefficiency_uncy_l%i", ij);
    defefficiency_uncy[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
    sprintf(title, "deftriggereffi_x_l%i", ij);
    deftriggereffi_x[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
    sprintf(title, "deftriggereffi_y_l%i", ij);
    deftriggereffi_y[ij]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
    for (int jk=0; jk<2*nmxiter; jk++) {
      sprintf(title, "totalentry_l%i_i%i", ij, jk);
      totalentry[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_l%i_i%i", ij, jk);
      inefficiency_uncx[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_l%i_i%i", ij, jk);
      inefficiency_uncy[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_l%i_i%i", ij, jk);
      inefficiency_corx[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_l%i_i%i", ij, jk);
      inefficiency_cory[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      
      sprintf(title, "totalentry_5_l%i_i%i", ij, jk);
      totalentry_5[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_5_l%i_i%i", ij, jk);
      inefficiency_uncx_5[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(title, "inefficiency_uncy_5_l%i_i%i", ij, jk);
      inefficiency_uncy_5[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_5_l%i_i%i", ij, jk);
      inefficiency_corx_5[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_5_l%i_i%i", ij, jk);
      inefficiency_cory_5[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "totalentry_6_l%i_i%i", ij, jk);
      totalentry_6[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_6_l%i_i%i", ij, jk);
      inefficiency_uncx_6[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_6_l%i_i%i", ij, jk);
      inefficiency_uncy_6[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_6_l%i_i%i", ij, jk);
      inefficiency_corx_6[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_6_l%i_i%i", ij, jk);
      inefficiency_cory_6[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "totalentry_7_l%i_i%i", ij, jk);
      totalentry_7[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_7_l%i_i%i", ij, jk);
      inefficiency_uncx_7[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_7_l%i_i%i", ij, jk);
      inefficiency_uncy_7[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_7_l%i_i%i", ij, jk);
      inefficiency_corx_7[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_7_l%i_i%i", ij, jk);
      inefficiency_cory_7[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "totalentry_8_l%i_i%i", ij, jk);
      totalentry_8[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_8_l%i_i%i", ij, jk);
      inefficiency_uncx_8[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_8_l%i_i%i", ij, jk);
      inefficiency_uncy_8[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_8_l%i_i%i", ij, jk);
      inefficiency_corx_8[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_8_l%i_i%i", ij, jk);
      inefficiency_cory_8[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(title, "totalentry_9_l%i_i%i", ij, jk);
      totalentry_9[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_9_uncx_l%i_i%i", ij, jk);
      inefficiency_uncx_9[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_9_l%i_i%i", ij, jk);
      inefficiency_uncy_9[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_9_l%i_i%i", ij, jk);
      inefficiency_corx_9[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_9_l%i_i%i", ij, jk);
      inefficiency_cory_9[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "totalentry_10_l%i_i%i", ij, jk);
      totalentry_10[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_10_l%i_i%i", ij, jk);
      inefficiency_uncx_10[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_10_l%i_i%i", ij, jk);
      inefficiency_uncy_10[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_10_l%i_i%i", ij, jk);
      inefficiency_corx_10[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_10_l%i_i%i", ij, jk);
      inefficiency_cory_10[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "totalentry_11_l%i_i%i", ij, jk);
      totalentry_11[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncx_11_l%i_i%i", ij, jk);
      inefficiency_uncx_11[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_uncy_11_l%i_i%i", ij, jk);
      inefficiency_uncy_11[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_corx_11_l%i_i%i", ij, jk);
      inefficiency_corx_11[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_cory_11_l%i_i%i", ij, jk);
      inefficiency_cory_11[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_x_5_l%i_i%i", ij, jk);
      triggereffi_x_5[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_5_l%i_i%i", ij, jk);
      triggereffi_y_5[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(title, "triggereffi_x_6_l%i_i%i", ij, jk);
      triggereffi_x_6[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_6_l%i_i%i", ij, jk);
      triggereffi_y_6[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(title, "triggereffi_x_7_l%i_i%i", ij, jk);
      triggereffi_x_7[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_7_l%i_i%i", ij, jk);
      triggereffi_y_7[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(title, "triggereffi_x_8_l%i_i%i", ij, jk);
      triggereffi_x_8[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_8_l%i_i%i", ij, jk);
      triggereffi_y_8[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(title, "triggereffi_x_9_l%i_i%i", ij, jk);
      triggereffi_x_9[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_9_l%i_i%i", ij, jk);
      triggereffi_y_9[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      
      sprintf(title, "triggereffi_x_10_l%i_i%i", ij, jk);
      triggereffi_x_10[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_10_l%i_i%i", ij, jk);
      triggereffi_y_10[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
     
      sprintf(title, "triggereffi_x_11_l%i_i%i", ij, jk);
      triggereffi_x_11[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_11_l%i_i%i", ij, jk);
      triggereffi_y_11[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
      

      sprintf(title, "inefficiency_xpixel_l%i_i%i", ij, jk);
      inefficiency_xpixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "inefficiency_ypixel_l%i_i%i", ij, jk);
      inefficiency_ypixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "inefficiency_xallpixel_l%i_i%i", ij, jk);
      inefficiency_xallpixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "inefficiency_yallpixel_l%i_i%i", ij, jk);
      inefficiency_yallpixel[ij][jk]=new TH1D(title, title, 60, -0.5, 0.5);

      sprintf(title, "inefficiency_xt_l%i_i%i", ij, jk);
      inefficiency_xt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "inefficiency_yt_l%i_i%i", ij, jk);
      inefficiency_yt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "total_xt_l%i_i%i", ij, jk);
      total_xt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "total_yt_l%i_i%i", ij, jk);
      total_yt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_x_l%i_i%i", ij, jk);
      triggereffi_x[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_y_l%i_i%i", ij, jk);
      triggereffi_y[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_xevt_l%i_i%i", ij, jk);
      triggereffi_xevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "triggereffi_yevt_l%i_i%i", ij, jk);
      triggereffi_yevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_uncx_l%i_i%i", ij, jk);
      difefficiency_uncx[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_uncy_l%i_i%i", ij, jk);
      difefficiency_uncy[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_xt_l%i_i%i", ij, jk);
      difefficiency_xt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "difefficiency_yt_l%i_i%i", ij, jk);
      difefficiency_yt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_x_l%i_i%i", ij, jk);
      diftriggereffi_x[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_y_l%i_i%i", ij, jk);
      diftriggereffi_y[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_xevt_l%i_i%i", ij, jk);
      diftriggereffi_xevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

      sprintf(title, "diftriggereffi_yevt_l%i_i%i", ij, jk);
      diftriggereffi_yevt[ij][jk]=new TH2D(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    }
#endif
  }
  
  TH2F* time_layer = new TH2F("time_layer", "time_layer", 2*nTDCpLayer*nlayer+1, -0.5, 2*nTDCpLayer*nlayer+0.5, 140, -0.5, 1399.5);
  
  TH1F* rawtimex = new TH1F("time_rawx", "time_rawx",1800,0.,1800.);
  TH1F* rawtimey = new TH1F("time_rawy", "time_rawy",1800,0.,1800.);
  TH1F* xlay_timediff[11];
  TH1F* ylay_timediff[11];
  TH1F* xlay_timediff_l1[6];
  TH1F* xlay_timediff_l4[5];
  TH1F* ylay_timediff_l1[6];
  TH1F* ylay_timediff_l4[5];
  for(int ij=0;ij<nlayer-6;ij++) {
    sprintf(name,"xlayer_timediff_l1_%i%i",1,ij+4);
    xlay_timediff_l1[ij] = new TH1F(name,name,500,-49.5,49.5);
    sprintf(name,"ylayer_timediff_l1_%i%i",1,ij+4);
    ylay_timediff_l1[ij] = new TH1F(name,name,500,-49.5,49.5);
  }
  for(int ij=0;ij<nlayer-7;ij++) {
    sprintf(name,"xlayer_timediff_l4_%i%i",4,ij+5);
    xlay_timediff_l4[ij] = new TH1F(name,name,500,-49.5,49.5);
    sprintf(name,"ylayer_timediff_l4_%i%i",4,ij+5);
    ylay_timediff_l4[ij] = new TH1F(name,name,500,-49.5,49.5);
  }
  
  for(int ij=0;ij<nlayer-1;ij++){
    sprintf(name,"xlayer_timediff_%i%i",ij,ij+1);
    xlay_timediff[ij] = new TH1F(name,name,500,-49.5,49.5);
    sprintf(name,"ylayer_timediff_%i%i",ij,ij+1);
    ylay_timediff[ij] = new TH1F(name,name,500,-49.5,49.5);
  }
  TH2F* time_layerstrip = new TH2F("time_layerstrip", "time_layerstrip", 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5, 2000, -0.5, 1999.5);
  
  TH1F* xlayer_alloccu=new TH1F("xlayer_alloccu","xlayer_alloccu",nlayer*nstrip,-0.5, nlayer*nstrip-0.5);
  TH1F* ylayer_alloccu=new TH1F("ylayer_alloccu","ylayer_alloccu",nlayer*nstrip,-0.5, nlayer*nstrip-0.5); 
  
  TH1F* xlayer_alloccusel=new TH1F("xlayer_alloccusel","xlayer_alloccusel",nlayer*nstrip,-0.5, nlayer*nstrip-0.5);
  TH1F* ylayer_alloccusel=new TH1F("ylayer_alloccusel","ylayer_alloccusel",nlayer*nstrip,-0.5, nlayer*nstrip-0.5); 
  
  int nbin = 400000;

  TH1F* trigrate = new TH1F("trigrate", "Trigger rate (Hz)", 100, -0.5, 99.5);
  TH1F* muposxrate = new TH1F("muposxrate", "Muon (xpos) rate (Hz)", 100, -0.5, 99.5);
  TH1F* muposyrate = new TH1F("muposyrate", "Muon (ypos) rate (Hz)", 100, -0.5, 99.5);
  TH1F* mutimexrate = new TH1F("mutimexrate", "Muon (xtime) rate (Hz)", 100, -0.5, 99.5);
  TH1F* mutimeyrate = new TH1F("mutimeyrate", "Muon (ytime) rate (Hz)", 100, -0.5, 99.5);

  TH1F* trigratex = new TH1F("trigratex", "Trigger ratex (Hz)", nbin, -0.5, nbin-0.5);  
  TH1F* muposxratex = new TH1F("muposxratex", "Muon (xpos) rate (Hz)", nbin, -0.5, nbin-0.5);
  TH1F* muposyratex = new TH1F("muposyratex", "Muon (ypos) rate (Hz)", nbin, -0.5, nbin-0.5);
  TH1F* mutimexratex = new TH1F("mutimexratex", "Muon (xtime) rate (Hz)", nbin, -0.5, nbin-0.5);
  TH1F* mutimeyratex = new TH1F("mutimeyratex", "Muon (ytime) rate (Hz)", nbin, -0.5, nbin-0.5);

  TH1F* xlayer_occu[nlayer];
  TH1F* ylayer_occu[nlayer];
  
  TH1F* xlayer_seloccu[nlayer];
  TH1F* ylayer_seloccu[nlayer];

  TH1F* xlayer_sel2occu[nlayer];
  TH1F* ylayer_sel2occu[nlayer];

  TH1F* xlayer_noiseoccu[nlayer];
  TH1F* ylayer_noiseoccu[nlayer];
  
  TH2F* raw_occu[nlayer];
  TH2F* raw_seloccu[nlayer];
  TH2F* raw_noiseoccu[nlayer];
  
  for (int ij=0; ij<nlayer; ij++) {
    sprintf(title, "xlayer_occu_l%i", ij);
    xlayer_occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "ylayer_occu_l%i", ij);
    ylayer_occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "raw_occu_l%i", ij);
    raw_occu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
    sprintf(title, "raw_seloccu_l%i", ij);
    raw_seloccu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
   
    sprintf(title, "raw_noiseoccu_l%i", ij);
    raw_noiseoccu[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
    
    sprintf(title, "xlayer_seloccu_l%i", ij);
    xlayer_seloccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
    
    sprintf(title, "ylayer_seloccu_l%i", ij);
    ylayer_seloccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "xlayer_sel2occu_l%i", ij);
    xlayer_sel2occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "ylayer_sel2occu_l%i", ij);
    ylayer_sel2occu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "xlayer_noiseoccu_l%i", ij);
    xlayer_noiseoccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);

    sprintf(title, "ylayer_noiseoccu_l%i", ij);
    ylayer_noiseoccu[ij] = new TH1F(title, title, nstrip, -0.5, nstrip-0.5);
  }
  
  TH1F* xlayer_mult[nlayer];
  TH1F* ylayer_mult[nlayer];
  
  TH1F* xlayer_allmult[nlayer];
  TH1F* ylayer_allmult[nlayer];
  
  TH1F* xlayer_allmumult[nlayer];
  TH1F* ylayer_allmumult[nlayer];

  TH1F* xlayer_allmutimemult[nlayer];
  TH1F* ylayer_allmutimemult[nlayer];

  for (int ij=0; ij<nlayer; ij++) {
    sprintf(title, "xlayer_mult_l%i", ij);
    xlayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
    
    sprintf(title, "ylayer_mult_l%i", ij);
    ylayer_mult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "xlayer_allmult_l%i", ij);
    xlayer_allmult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
    
    sprintf(title, "ylayer_allmult_l%i", ij);
    ylayer_allmult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "xlayer_allmumult_l%i", ij);
    xlayer_allmumult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
    
    sprintf(title, "ylayer_allmumult_l%i", ij);
    ylayer_allmumult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

    sprintf(title, "xlayer_allmutimemult_l%i", ij);
    xlayer_allmutimemult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);
    
    sprintf(title, "ylayer_allmutimemult_l%i", ij);
    ylayer_allmutimemult[ij] = new TH1F(title, title, nstrip+1, -0.5, nstrip+0.5);

  } 
  
  TH1F* xstrip_mult = new TH1F("xstrip_mult", "xstrip_mult", 6*nlayer, -0.5,  6*nlayer-0.5);
  TH1F* ystrip_mult = new TH1F("ystrip_mult", "ystrip_mult", 6*nlayer, -0.5,  6*nlayer-0.5);
  
  TH2F* dir_cxchi = new TH2F("dir_cxchi", "dir_cxchi", 100, 0., 100., 200, -20., 20.);
  TH2F* dir_cychi = new TH2F("dir_cychi", "dir_cychi", 100, 0., 100., 200, -20., 20.);
  
  TH2F* dir_cxy = new TH2F("dir_cxy", "dir_cxy", 100, -2., 2., 100, -2., 2.);
  
  TH1F* timex_shift[nlayer][nmxiter+1];
  TH1F* timey_shift[nlayer][nmxiter+1];
  
  TH2F* timex_2dshift[nlayer][nmxiter+1];
  TH2F* timey_2dshift[nlayer][nmxiter+1];  

  //  TH1F* time_xstshift[nlayer][nstrip][2*nmxiter];
  //  TH1F* time_ystshift[nlayer][nstrip][2*nmxiter];
  
  TH2F* timex_fy[nlayer];
  TH2F* timey_fx[nlayer];

  TH2F* timex_fy2[nlayer];
  TH2F* timey_fx2[nlayer];

  TProfile* timex_pfy[nlayer];
  TProfile* timey_pfx[nlayer];

  TProfile* timex_pfy2[nlayer];
  TProfile* timey_pfx2[nlayer];
  if (isTiming) {
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<=nmxiter; jk++) {
        sprintf(title, "timex_shift_l%i_i%i", ij, jk);
        timex_shift[ij][jk] = new TH1F(title, title, 120, -50.0, 50.0);
        
        sprintf(title, "timey_shift_l%i_i%i", ij, jk);
        timey_shift[ij][jk] = new TH1F(title, title, 120, -50.0, 50.0);
        
        sprintf(title, "timex_2dshift_l%i_i%i", ij, jk);
        timex_2dshift[ij][jk] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 60, -30.0, 30.0);
        
        sprintf(title, "timey_2dshift_l%i_i%i", ij, jk);
        timey_2dshift[ij][jk] = new TH2F(title, title,  nstrip, -0.5, nstrip-0.5, 60, -30.0, 30.0);
      }
      sprintf(title, "timex_fy_l%i", ij);
      timex_fy[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);
      sprintf(title, "timey_fx_l%i", ij);
      timey_fx[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);
      sprintf(title, "timex_fy2_l%i", ij);
      timex_fy2[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);
      sprintf(title, "timey_fx2_l%i", ij);
      timey_fx2[ij] = new TH2F(title, title, nstrip*nstrip, -0.5, nstrip*nstrip-0.5, 120, -12., 12.);
      
      sprintf(title, "timex_pfy_l%i", ij);
      timex_pfy[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
      sprintf(title, "timey_pfx_l%i", ij);
      timey_pfx[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
      sprintf(title, "timex_pfy2_l%i", ij);
      timex_pfy2[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
      sprintf(title, "timey_pfx2_l%i", ij);
      timey_pfx2[ij] = new TProfile(title, title, nstrip*(nstrip+1), -0.5, nstrip*(nstrip+1)-0.5, -12., 12.);
    }
  }

  TH1F* xlayer_reso[nlayer][2*nmxiter];
  TH1F* ylayer_reso[nlayer][2*nmxiter];
  TH2F* xylayer_reso[nlayer][2*nmxiter];
  bool is_xyreso[nlayer][2*nmxiter];

  TH1F* xlayer_reso_mul[nlayer][2*nmxiter][nmxhits];
  TH1F* ylayer_reso_mul[nlayer][2*nmxiter][nmxhits];

  TH1F* xlayer_exterr[nlayer][2*nmxiter];
  TH1F* ylayer_exterr[nlayer][2*nmxiter];

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      double xposmx=(jk<nmxiter) ? 5.9 : 6.0;
      sprintf(title, "xlayer_reso_l%i_i%i", ij, jk);
      xlayer_reso[ij][jk]=new TH1F(title, title, 150, -xposmx, xposmx);
      sprintf(title, "ylayer_reso_l%i_i%i", ij, jk);
      ylayer_reso[ij][jk]=new TH1F(title, title, 150, -xposmx, xposmx);

      sprintf(title, "xylayer_reso_l%i_i%i", ij, jk);
      xylayer_reso[ij][jk]=new TH2F(title, title, 90, -xposmx, xposmx, 90, -xposmx, xposmx);

      sprintf(title, "xlayer_exterr_l%i_i%i", ij, jk);
      xlayer_exterr[ij][jk]=new TH1F(title, title, 120, 0.0, 1.0);
      sprintf(title, "ylayer_exterr_l%i_i%i", ij, jk);
      ylayer_exterr[ij][jk]=new TH1F(title, title, 120, 0.0, 1.0);

      for (int kl=0; kl<nmxhits; kl++) {
        sprintf(title, "xlayer_reso_l%i_i%i_mul%i", ij, jk, kl+1);
        xlayer_reso_mul[ij][jk][kl]=new TH1F(title, title, 150, -xposmx, xposmx);
        sprintf(title, "ylayer_reso_l%i_i%i_mul%i", ij, jk, kl+1);
        ylayer_reso_mul[ij][jk][kl]=new TH1F(title, title, 150, -xposmx, xposmx);
        
      }
    }
  }
  
  const int nstr_posmx=nmxhits; //3; //This could be same as nmxhits;
  TH2F* xstr_xdev[nlayer][nmxiter][nstr_posmx];
  TH2F* ystr_xdev[nlayer][nmxiter][nstr_posmx];
  TH2F* xstr_ydev[nlayer][nmxiter][nstr_posmx];
  TH2F* ystr_ydev[nlayer][nmxiter][nstr_posmx];

  TH2F* xstr_xtdev[nlayer][nmxiter];
  TH2F* ystr_xtdev[nlayer][nmxiter];
  TH2F* xstr_ytdev[nlayer][nmxiter];
  TH2F* ystr_ytdev[nlayer][nmxiter];

  TH2F* nxystr_xtdev[nlayer][nmxiter+6];
  TH2F* xystr_xtdev[nlayer][nmxiter+6]; // + rawytime, rawytime1, timesy & ytimewo, rawytime2 & ryawtime3
  
  TH2F* nxystr_ytdev[nlayer][nmxiter+6];
  TH2F* xystr_ytdev[nlayer][nmxiter+6];

  TH2F* nxypos_xytdev[8];
  for (int ij=0; ij<8; ij++) {
    sprintf(name, "nxypos_xytdev_%i", ij);
    nxypos_xytdev[ij] = new TH2F(name, name, 2*(nstrip+2), -1.5, nstrip+0.5, 2*(nstrip+2), -1.5, nstrip+0.5);
  }
  TH2F* xtdev_ytdev[4];
  for (int ij=0; ij<4; ij++) {
    sprintf(name, "xtdev_ytdev_%i", ij);
    xtdev_ytdev[ij] = new TH2F(name, name, 120, -30., 30., 120, -30., 30.);
  }
  
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nmxiter; jk++) {
      for (int kl=0; kl<nstr_posmx; kl++) {
        sprintf(name, "xstr_xdev_l%i_i%i_mul%i", ij, jk, kl+1);    
        xstr_xdev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
        sprintf(name, "ystr_xdev_l%i_i%i_mul%i", ij, jk, kl+1);
        ystr_xdev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
        sprintf(name, "xstr_ydev_l%i_i%i_mul%i", ij, jk, kl+1);
        xstr_ydev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
        sprintf(name, "ystr_ydev_l%i_i%i_mul%i", ij, jk, kl+1);
        ystr_ydev[ij][jk][kl] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 90, -4.5, 4.5);
      }
    }
  }
  
  if (isTiming) { 
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<nmxiter; jk++) {
        sprintf(name, "xstr_xtdev_l%i_i%i", ij, jk);    
        xstr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
        sprintf(name, "ystr_xtdev_l%i_i%i", ij, jk);    
        ystr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
        sprintf(name, "xstr_ytdev_l%i_i%i", ij, jk);    
        xstr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
        sprintf(name, "ystr_ytdev_l%i_i%i", ij, jk); 
        ystr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, 120, -18., 18.);
      }
      
      for (int jk=0; jk<nmxiter+6; jk++) {
        sprintf(name, "nxystr_xtdev_l%i_i%i", ij, jk);    
        nxystr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
        
        sprintf(name, "nxystr_ytdev_l%i_i%i", ij, jk);    
        nxystr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
        
        sprintf(name, "xystr_xtdev_l%i_i%i", ij, jk);    
        xystr_xtdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
        
        sprintf(name, "xystr_ytdev_l%i_i%i", ij, jk);    
        xystr_ytdev[ij][jk] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);	
      }
    }
  }

  TH2F* shift_pos;
  TH2F* rms_pos;
  if (isalign>0) {
    shift_pos = new TH2F("shift_pos", "shift_pos", 2*nlayer, -0.5, 2*nlayer-0.5, nmxiter, -0.5, nmxiter-0.5);
    rms_pos = new TH2F("rms_pos", "rms_pos", 2*nlayer, -0.5, 2*nlayer-0.5, nmxiter, -0.5, nmxiter-0.5);
  }
  
  TH1F* passed_strip[nmxiter]; 
  for (int ij=0; ij<nmxiter; ij++) {
    sprintf(title, "passed_strip_i%i", ij);
    passed_strip[ij] = new TH1F(title, title, 4*(nlayer+2), -0.5, 4*(nlayer+2)-0.5);
  } 
  
  TH1F* h_chisqx = new TH1F("chisqx", "chisqx", 120, 0.0, 150.0);
  TH1F* h_reduchisqx = new TH1F("reduchisqx", "reduced chisqx", 90, 0.0, 30.0/*45.0*/);
  TH1F* h_chisqy = new TH1F("chi2y", "chisqy", 120, 0.0, 150.0);
  TH1F* h_reduchisqy = new TH1F("reduchisqy", "reduced chisqy", 90, 0.0, 30.0/*45.0*/);
  TH1F* h_xndf = new TH1F("xndf", "xndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_yndf = new TH1F("yndf", "yndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_xprob = new TH1F("xprob", "xprob", 120, 0.0, 1.0);
  TH1F* h_yprob = new TH1F("yprob", "yprob", 120, 0.0, 1.0);

  TH1F* h_tmp1xndf = new TH1F("tmp1xndf", "tmp1xndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_tmp2xndf = new TH1F("tmp2xndf", "tmp2xndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_tmp3xndf = new TH1F("tmp3xndf", "tmp3xndf", nlayer, 0.5, nlayer+0.5);

  TH1F* h_tchisqx = new TH1F("tchisqx", "tchisqx", 120, 0.0, 90.0);
  TH1F* h_treduchisqx = new TH1F("treduchisqx", "reduced tchisqx", 90, 0.0, 30.0);
  TH1F* h_tchisqy = new TH1F("tchi2y", "tchisqy", 120, 0.0, 90.0);
  TH1F* h_treduchisqy = new TH1F("treduchisqy", "reduced tchisqy", 90, 0.0, 30.0);
  TH1F* h_txndf = new TH1F("txndf", "txndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_tyndf = new TH1F("tyndf", "tyndf", nlayer, 0.5, nlayer+0.5);
  TH1F* h_xtprob = new TH1F("xtprob", "xtprob", 120, 0.0, 1.0);
  TH1F* h_ytprob = new TH1F("ytprob", "ytprob", 120, 0.0, 1.0);

  const int nprob=6;
  double probs[nprob+1]={3.5, 5.5, 7.5, 8.5, 9.5, 10.5, 12.5};
  //  double probs[nprob+1]={2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 12.5};

  TH1F* h_xnprob[nprob];
  TH1F* h_ynprob[nprob];  
  TH1F* h_xtnprob[nprob];
  TH1F* h_ytnprob[nprob];  
  
  for (int ij=0; ij<nprob; ij++) {
    sprintf(name, "h_xnprob_%i", ij);
    sprintf(title, "h_xnprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_xnprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);
    
    sprintf(name, "h_ynprob_%i", ij);
    sprintf(title, "h_ynprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_ynprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);
    
    sprintf(name, "h_xtnprob_%i", ij);
    sprintf(title, "h_xtnprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_xtnprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);
    
    sprintf(name, "h_ytnprob_%i", ij);
    sprintf(title, "h_ytnprob_%i_%i", int(probs[ij]+1), int(probs[ij+1]));
    h_ytnprob[ij] = new TH1F(name, title, 120, 0.0, 1.0);
  }

  TH2F* correction_xtime[nlayer];
  TH2F* correction_ytime[nlayer];

  TH2F* fitted_rms_xtime[nlayer];
  TH2F* fitted_rms_ytime[nlayer];

  TH1F* time_xreso[nlayer][2*nmxiter];
  TH1F* time_yreso[nlayer][2*nmxiter];
  TH2F* time_xyreso[nlayer][2*nmxiter];  
  bool istime_xyreso[nlayer][2*nmxiter];

  bool xposEffPass[nlayer][2*nmxiter];
  bool yposEffPass[nlayer][2*nmxiter];

  TH1F* time_xstrreso[nlayer][nstrip][2*nmxiter];
  TH1F* time_ystrreso[nlayer][nstrip][2*nmxiter];

  TH2F* time_mean_reso = new TH2F("time_mean_reso", "time_mean_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);
  TH2F* time_rms_reso = new TH2F("time_rms_reso", "time_rms_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);
  TH2F* time_corrms_reso = new TH2F("time_corrms_reso", "time_corrms_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);
  TH2F* time_exterr_reso = new TH2F("time_exterr_reso", "time_exterr_reso", 2*nlayer, -0.5, 2*nlayer-0.5, 2*nmxiter, -0.5, 2*nmxiter-0.5);
  
  TH1F* dir_cxlay[nlayerit+1][nmxiter][nprob];
  TH1F* dir_cylay[nlayerit+1][nmxiter][nprob];

  TH1F* dir_cx[nlayerit+1][nmxiter];
  TH1F* dir_cy[nlayerit+1][nmxiter];

  TH1F* dir_cx2[nlayerit+1][nmxiter];
  TH1F* dir_cy2[nlayerit+1][nmxiter];

  TH1F* dir_c0x[nlayerit+1][nmxiter];
  TH1F* dir_c0y[nlayerit+1][nmxiter];

  TH1F* time_offset[nmxiter];
  TH1F* shift_time_mnft[nmxiter];
  TH1F* statmean_time[nmxiter];
  TH1F* statrms_time[nmxiter];
  TH1F* statskew_time[nmxiter];
  TH1F* statkurt_time[nmxiter];
 
  TH1F* rms_time[nmxiter];
  TH1F* rms_timeused[nmxiter];

  TH1F* time_offsetx[nmxiter];
  TH1F* shift_time_mnftx[nmxiter];
  TH1F* statmean_timex[nmxiter];
  TH1F* statrms_timex[nmxiter];
  TH1F* statskew_timex[nmxiter];
  TH1F* statkurt_timex[nmxiter];

  TH1F* rms_timex[nmxiter];
  TH1F* rms_timeusedx[nmxiter];

  TH1F* time_offsety[nmxiter];
  TH1F* shift_time_mnfty[nmxiter];
  TH1F* statmean_timey[nmxiter];
  TH1F* statrms_timey[nmxiter];
  TH1F* statskew_timey[nmxiter];
  TH1F* statkurt_timey[nmxiter];

  TH1F* rms_timey[nmxiter];
  TH1F* rms_timeusedy[nmxiter];

  TH1F* time_mulxreso[nlayer][2*nmxiter][nmxtimehit];
  TH1F* time_mulyreso[nlayer][2*nmxiter][nmxtimehit];

  TH1F* xtime_exterr[nlayer][2*nmxiter];
  TH1F* ytime_exterr[nlayer][2*nmxiter];
  
  TH2F* rawhits_corr_xymul[nlayer];
  TH2F* rawhits_xlay_corr_mul[nlayer][nlayer];
  TH2F* rawhits_ylay_corr_mul[nlayer][nlayer];
  
  for (int ij=0; ij<nlayer; ij++) {
    sprintf(name, "rawhits_corr_xymul_l%i", ij);
    rawhits_corr_xymul[ij] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
    for (int jk=ij+1; jk<nlayer; jk++) {
      sprintf(name, "rawhits_xlay_corr_mul_l%i_l%i", ij, jk);
      rawhits_xlay_corr_mul[ij][jk] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
      
      sprintf(name, "rawhits_ylay_corr_mul_l%i_l%i", ij, jk);
      rawhits_ylay_corr_mul[ij][jk] = new TH2F(name, name, nstrip+1, -0.5, nstrip+0.5,  nstrip+1, -0.5, nstrip+0.5);
    }
  }  
  
  if (isTiming) {
    if (isalign>0) {
      for (int ij=0; ij<nlayer; ij++) {
        sprintf(title, "correction_xtime_l%i", ij);
        correction_xtime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);
        
        sprintf(title, "correction_ytime_l%i", ij);
        correction_ytime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);
        
        sprintf(title, "fitted_rms_xtime_l%i", ij);
        fitted_rms_xtime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);
        
        sprintf(title, "fitted_rms_ytime_l%i", ij);
        fitted_rms_ytime[ij] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nmxiter, -0.5, nmxiter-0.5);
      }
    }
    
    for (int kl=0; kl<2*nmxiter; kl++ ){
      double trange = (kl<nmxiter) ? 25. : 20.;//7.5 : 10.;
      double trange2 = (kl<nmxiter) ? 25. : 20.; //7.5 : 10.;
      for (int ij=0; ij<nlayer; ij++) {
        sprintf(title, "time_xreso_l%i_i%i", ij, kl);
        time_xreso[ij][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);
        
        sprintf(title, "time_yreso_l%i_i%i", ij, kl);
        time_yreso[ij][kl] = new TH1F(title, title, 120, -trange2, trange2); //150, -80.0, 70.0); //120, -trange2, trange2);   
        
        sprintf(title, "time_xyreso_l%i_i%i", ij, kl);
        time_xyreso[ij][kl] = new TH2F(title, title, 120, -trange2-1., trange2-1., 120, -trange2-1., trange2-1.);   
        
        for (int lm=0; lm<nmxtimehit; lm++) {
          sprintf(title, "time_mulxreso_l%i_i%i_m%i", ij, kl, lm+1);
          time_mulxreso[ij][kl][lm] = new TH1F(title, title, 120, -trange2, trange2); 
          
          sprintf(title, "time_mulyreso_l%i_i%i_m%i", ij, kl, lm+1);
          time_mulyreso[ij][kl][lm] = new TH1F(title, title, 120, -trange2, trange2); 
        }
        
        sprintf(title, "xtime_exterr_l%i_i%i", ij, kl);
        xtime_exterr[ij][kl]=new TH1F(title, title, 120, 0.0, 2.1);
        sprintf(title, "ytime_exterr_l%i_i%i", ij, kl);
        ytime_exterr[ij][kl]=new TH1F(title, title, 120, 0.0, 2.1);
        
        for (int jk=0; jk<nstrip; jk++) {
          sprintf(title, "time_xstrreso_l%i_s%i_i%i", ij, jk, kl);
          time_xstrreso[ij][jk][kl] = new TH1F(title, title, 120, -trange, trange);
          
          sprintf(title, "time_ystrreso_l%i_s%i_i%i", ij, jk, kl);
          time_ystrreso[ij][jk][kl] = new TH1F(title, title, 120, -trange, trange);   
          
          //	  sprintf(title, "time_xstshift_l%i_s%i_i%i", ij, jk, kl);
          //	  time_xstshift[ij][jk][kl] = new TH1F(title, title, 180, -trange, trange);
          
          //	  sprintf(title, "time_ystshift_l%i_s%i_i%i", ij, jk, kl);
          //	  time_ystshift[ij][jk][kl] = new TH1F(title, title, 180, -trange, trange);
        }
      }
    }
    for (int kl=0; kl<nmxiter; kl++) {
      for (int ij=0; ij<=nlayerit; ij++) {
        for (int lm=0; lm<nprob; lm++) {
          sprintf(name, "dir_cxlay_l%i_i%i_n%i", ij, kl, lm);
          sprintf(title, "dir_cxlay_l%i_i%i_n[%i-%i]", ij, kl, int(probs[lm]+1), int(probs[lm+1]));
          dir_cxlay[ij][kl][lm] = new TH1F(name, title, 180, -3.5, 5.5);
          sprintf(name, "dir_cylay_l%i_i%i_n%i", ij, kl, lm);
          sprintf(title, "dir_cylay_l%i_i%i_n[%i-%i]", ij, kl, int(probs[lm]+1), int(probs[lm+1]));
          
          dir_cylay[ij][kl][lm] = new TH1F(name, title, 180, -3.5, 5.5);
        }
        
        sprintf(title, "dir_cx_l%i_i%i", ij, kl);
        dir_cx[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);
        
        sprintf(title, "dir_cy_l%i_i%i", ij, kl);
        dir_cy[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);
        
        sprintf(title, "dir_cx2_l%i_i%i", ij, kl);
        dir_cx2[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);
        
        sprintf(title, "dir_cy2_l%i_i%i", ij, kl);
        dir_cy2[ij][kl] = new TH1F(title, title, 180, -3.5, 5.5);
        
        sprintf(title, "dir_c0x_l%i_i%i", ij, kl);
        dir_c0x[ij][kl] = new TH1F(title, title, 120, 70., 130.);
        
        sprintf(title, "dir_c0y_l%i_i%i", ij, kl);
        dir_c0y[ij][kl] = new TH1F(title, title, 120, 70., 130.);
      }
    }
    
    if (isalign>0) {
      for (int kl=0; kl<nmxiter; kl++) {
        sprintf(title, "time_offset_i%i", kl);
        time_offset[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);
        
        sprintf(title, "shift_time_mnft_i%i", kl);
        shift_time_mnft[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);
        
        sprintf(title, "statmean_time_i%i", kl);
        statmean_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);
        
        sprintf(title, "statrms_time_i%i", kl);
        statrms_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);
        
        sprintf(title, "statskew_time_i%i", kl);
        statskew_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);
        
        sprintf(title, "statkurt_time_i%i", kl);
        statkurt_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);
        
        sprintf(title, "rms_time_i%i", kl);
        rms_time[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);
        
        sprintf(title, "rms_timeused_i%i", kl);
        rms_timeused[kl] = new TH1F(title, title, 2*nlayer*nstrip, -0.5, 2*nlayer*nstrip-0.5);   
        
        /////////////////////////////////////////////
        sprintf(title, "time_offsetx_i%i", kl);
        time_offsetx[kl] = new TH1F(title, title, 60, -7.0, 8.0);
        
        sprintf(title, "shift_time_mnftx_i%i", kl);
        shift_time_mnftx[kl] = new TH1F(title, title, 60, -0.3, 0.3);
        
        sprintf(title, "statmean_timex_i%i", kl);
        statmean_timex[kl] = new TH1F(title, title, 60, -1.0, 1.0);
        
        sprintf(title, "statrms_timex_i%i", kl);
        statrms_timex[kl] = new TH1F(title, title, 60, 0.6, 3.0);
        
        sprintf(title, "statskew_timex_i%i", kl);
        statskew_timex[kl] = new TH1F(title, title, 60, -1.0, 3.0);
        
        sprintf(title, "statkurt_timex_i%i", kl);
        statkurt_timex[kl] = new TH1F(title, title, 60, -2.0, 6.0);
        
        sprintf(title, "rms_timex_i%i", kl);
        rms_timex[kl] = new TH1F(title, title, 60, 0.6, 2.7);
        
        sprintf(title, "rms_timeusedx_i%i", kl);
        rms_timeusedx[kl] = new TH1F(title, title, 60, 0.6, 2.4);
        
        /////////////////////////////////////////////
        sprintf(title, "time_offsety_i%i", kl);
        time_offsety[kl] = new TH1F(title, title, 60, -7.0, 8.0);
        
        sprintf(title, "shift_time_mnfty_i%i", kl);
        shift_time_mnfty[kl] = new TH1F(title, title, 60, -0.3, 0.3);
        
        sprintf(title, "statmean_timey_i%i", kl);
        statmean_timey[kl] = new TH1F(title, title, 60, -1.0, 1.0);
        
        sprintf(title, "statrms_timey_i%i", kl);
        statrms_timey[kl] = new TH1F(title, title, 60, 0.6, 3.0);
        
        sprintf(title, "statskew_timey_i%i", kl);
        statskew_timey[kl] = new TH1F(title, title, 60, -1.0, 3.0);
        
        sprintf(title, "statkurt_timey_i%i", kl);
        statkurt_timey[kl] = new TH1F(title, title, 60, -2.0, 6.0);
        
        sprintf(title, "rms_timey_i%i", kl);
        rms_timey[kl] = new TH1F(title, title, 60, 0.6, 2.7);
        
        sprintf(title, "rms_timeusedy_i%i", kl);
        rms_timeusedy[kl] = new TH1F(title, title, 60, 0.5, 2.4);
      }
    }
  } //  if (isTiming)
  TH2F* h_xcorhits = new TH2F("xcorhits","xcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_ycorhits = new TH2F("ycorhits","ycorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xycorhits = new TH2F("xycorhits","xycorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xtcorhits = new TH2F("xtcorhits","xtcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_ytcorhits = new TH2F("ytcorhits","ytcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xytcorhits = new TH2F("xytcorhits","xytcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xrawcorhits = new TH2F("xrawcorhits","xrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_yrawcorhits = new TH2F("yrawcorhits","yrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xyrawcorhits = new TH2F("xyrawcorhits","xyrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xtrawcorhits = new TH2F("xtrawcorhits","xtrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_ytrawcorhits = new TH2F("ytrawcorhits","ytrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_xytrawcorhits = new TH2F("xytrawcorhits","xytrawcorhits",nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);

  TH2F* h_xcorstrips[nlayer];
  TH2F* h_ycorstrips[nlayer];
  TH2F* h_xycorstrips[nlayer];
  TH2F* h_yxcorstrips[nlayer]; 

  TH2F* h_raw_xcorstrips[nlayer];
  TH2F* h_raw_ycorstrips[nlayer];
  TH2F* h_raw_xystrpnhits[nlayer];
  TH2F* h_raw_yxstrpnhits[nlayer]; 
  TH2F* h_raw_xstrpnhits[nlayer];
  TH2F* h_raw_ystrpnhits[nlayer];
  
  TH2F* h_xtcorstrips[nlayer];
  TH2F* h_ytcorstrips[nlayer];
  TH2F* h_xytcorstrips[nlayer]; 
  TH2F* h_yxtcorstrips[nlayer]; 

  for (int ij=0; ij<nlayer; ij++) {
    sprintf(name, "xcorstrips_l%i", ij);
    h_xcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "ycorstrips_l%i", ij);
    h_ycorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xycorstrips_l%i", ij);
    h_xycorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "yxcorstrips_l%i", ij);
    h_yxcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

     sprintf(name, "raw_xcorstrips_l%i", ij);
     h_raw_xcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "raw_ycorstrips_l%i", ij);
    h_raw_ycorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
    sprintf(name, "xstrpnhits_l%i", ij);
    h_raw_xstrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "ystrpnhits_l%i", ij);
    h_raw_ystrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xystrpnhits_l%i", ij);
    h_raw_xystrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "yxstrpnhits_l%i", ij);
    h_raw_yxstrpnhits[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    

    sprintf(name, "xtcorstrips_l%i", ij);
    h_xtcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "ytcorstrips_l%i", ij);
    h_ytcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "xytcorstrips_l%i", ij);
    h_xytcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);

    sprintf(name, "yxtcorstrips_l%i", ij);
    h_yxtcorstrips[ij] = new TH2F(name, name, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
    
  }

  TH1F* time_xraw[nlayer][nstrip];
  TH1F* time_yraw[nlayer][nstrip];
  TH2F* time_xraw2d = new TH2F("time_xraw2d", "time_xraw2d", nlayer, -0.5, nlayer-0.5, nstrip, -0.5, nstrip-0.5);
  TH2F* time_yraw2d = new TH2F("time_yraw2d", "time_yraw2d", nlayer, -0.5, nlayer-0.5, nstrip, -0.5, nstrip-0.5);

#ifdef TIMESLOPE
  TH2F* time_xslope[nlayer][nTDCpLayer];
  TH2F* time_yslope[nlayer][nTDCpLayer]; 

  TH2F* time_xslope_pr[nlayer][nstrip][nmxiter];
  TH2F* time_yslope_pr[nlayer][nstrip][nmxiter]; 
#endif
  for (int ij=0; ij <nlayer; ij++) {
    for (int jk=0; jk<nstrip; jk++) {
      sprintf(title, "time_xraw_l%i_s%i", ij, jk);
      time_xraw[ij][jk] = new TH1F(title, title, 120, 750., 1350.);
      sprintf(title, "time_yraw_l%i_s%i", ij, jk);
      time_yraw[ij][jk] = new TH1F(title, title, 120, 750., 1350.);
#ifdef TIMESLOPE
      for (int kl=0; kl<nmxiter; kl++) {
        
        sprintf(title, "time_xslope_pr_l%i_s%i_i%i", ij, jk, kl);
        time_xslope_pr[ij][jk][kl] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 90, -6.0, 6.0);
        sprintf(title, "time_yslope_pr_l%i_s%i_i%i", ij, jk, kl);
        time_yslope_pr[ij][jk][kl] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 90, -6.0, 6.0);
      }
#endif
      
    }
#ifdef TIMESLOPE
    for (int jk=0; jk<nTDCpLayer; jk++) { 
      sprintf(title, "time_xslope_l%i_tdc%i", ij, jk);
      time_xslope[ij][jk] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 130, 50., 180.);
      sprintf(title, "time_yslope_l%i_tdc%i", ij, jk);
      time_yslope[ij][jk] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, 130, 50., 180.);
    }
#endif
  }
  TH1F* timex_correl[nlayer][nlayer][nmxiter+1];
  TH1F* timey_correl[nlayer][nlayer][nmxiter+1];
  
  TH2F* time_both_cormean[nmxiter+1];
  TH2F* time_both_corrms[nmxiter+1];

  const int npixel=16;
  TH1F* timexy_correl[nlayer][npixel+1][nmxiter+1];

  TH2F* timexy_cormean[npixel+1];
  TH2F* timexy_corrms[npixel+1];

  const int ntimecor=2; //how many correlation do we want
  TProfile* indtimexy_prof[nlayer][ntimecor];
  TH1F* indtimexy_correl[nlayer][nstrip][nstrip][ntimecor];
  TH2F* indtimexy_cormean[nlayer][ntimecor];
  TH2F* indtimexy_corrms[nlayer][ntimecor];
  TH2F* indtimexy_fitmean[nlayer][ntimecor];
  TH2F* indtimexy_fitrms[nlayer][ntimecor];

  if (isTiming) { 
    for (int ij=0; ij<nlayer-1; ij++) {
      for (int jk=ij+1; jk<nlayer; jk++) {
        for (int kl=0; kl<nmxiter+1; kl++) {
          sprintf(title, "timex_correl_l%i_l%i_i%i", ij, jk, kl);
          timex_correl[ij][jk][kl] = new TH1F(title, title, 120, -8., 8.);
          
          sprintf(title, "timey_correl_l%i_l%i_i%i", ij, jk, kl);
          timey_correl[ij][jk][kl] = new TH1F(title, title, 120, -8., 8.);
        }
      }
    }
    
    for (int kl=0; kl<nmxiter+1; kl++) {
      sprintf(title, "time_both_cormean_i%i", kl);
      time_both_cormean[kl] = new TH2F(title, title, nlayer+1, -0.5, nlayer+0.5, nlayer, -0.5, nlayer-0.5);
      
      sprintf(title, "time_both_corrms_i%i", kl);
      time_both_corrms[kl] = new TH2F(title, title, nlayer+1, -0.5, nlayer+0.5, nlayer, -0.5, nlayer-0.5);
    }
    
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<npixel+1; jk++) {
        for (int kl=0; kl<nmxiter+1; kl++) {
          sprintf(title, "timexy_correl_l%i_pixel%i_i%i", ij, jk, kl);
          timexy_correl[ij][jk][kl] = new TH1F(title, title, 120, -8., 8.);
        }
      }
    }
    
    for (int kl=0; kl<npixel+1; kl++) {
      sprintf(title, "timexy_cormean_pixel%i", kl);
      timexy_cormean[kl] = new TH2F(title, title, nlayer, -0.5, nlayer-0.5, nmxiter+1, -0.5, nmxiter+0.5);
      sprintf(title, "timexy_corrms_pixel%i", kl);
      timexy_corrms[kl] = new TH2F(title, title, nlayer, -0.5, nlayer-0.5, nmxiter+1, -0.5, nmxiter+0.5);
    }
    
    for (int ij=0; ij<nlayer; ij++) {
      for (int lm=0; lm<ntimecor; lm++) {
        sprintf(title, "indtimexy_prof_l%i_%i", ij, lm);
        indtimexy_prof[ij][lm] = new TProfile(title, title, 2*nstrip, -0.5, 2*nstrip-0.5, -5.0, 15.0);
        
        sprintf(title, "indtimexy_cormean_l%i_%i", ij, lm);
        indtimexy_cormean[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
        sprintf(title, "indtimexy_corrms_l%i_%i", ij, lm);
        indtimexy_corrms[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
        
        sprintf(title, "indtimexy_fitmean_l%i_%i", ij, lm);
        indtimexy_fitmean[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
        sprintf(title, "indtimexy_fitrms_l%i_%i", ij, lm);
        indtimexy_fitrms[ij][lm] = new TH2F(title, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
        
        for (int jk=0; jk<nstrip; jk++) {
          for (int kl=0; kl<nstrip; kl++) {
            sprintf(name, "indtimexy_correl_l%i_x%i_y%i_%i", ij, jk, kl, lm);
            sprintf(title, "indtimexy_correl_l%i_x%i_y%i_%i", ij, jk+8, kl+8, lm);
            indtimexy_correl[ij][jk][kl][lm] = new TH1F(name, title, 120, -5.0, 15.0);
          }
        }
      }
    }
  } // if (isTiming)

  int xhits[nlayer],yhits[nlayer];   //number of hits after noise rejection for position fit
  int xallhits[nlayer][nTDCpLayer],yallhits[nlayer][nTDCpLayer]; // raw : for one hit, return strip number otherwise -10-multiplicity
  
  bool passxtime[nlayer], passytime[nlayer]; // passed through timing criteria in all hitted strips in x/ystr_x/ytdev
  bool passxtmy[nlayer], passytmx[nlayer]; //for X/Y strips, use boundary of Y/X extrapolation
  int istrxtime[nlayer], istrytime[nlayer]; //strip with early timing

  float errxco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
  float erryco[nlayer]={0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675,0.288675};
  
  double xxerr[nlayer], yyerr[nlayer];
  for ( int ix=0; ix<nlayer; ix++) {xxerr[ix] = yyerr[ix] = errxco[ix]*errxco[ix];}
  
  double xrms[nlayer]={0};
  double yrms[nlayer]={0};

  double deltaposxcov[nlayer][nlayer];
  int deltaposxCount[nlayer][nlayer];

  double deltaposycov[nlayer][nlayer];
  int deltaposyCount[nlayer][nlayer];

  double deltatcov[nlayer][nlayer];
  int deltatCount[nlayer][nlayer];

  double deltatcov2[nlayer][nlayer];
  int deltatCount2[nlayer][nlayer];

  double deltatcovy[nlayer][nlayer];
  int deltatCounty[nlayer][nlayer];

  TH2F* h_deltaposxcov = new TH2F("deltaposxcov", "deltaposxcov", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltaposycov = new TH2F("deltapUosycov", "deltaposycov", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltatcov = new TH2F("deltatcov", "deltatcov", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltatcov2 = new TH2F("deltatcov2", "deltatcov2", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_deltatcovy = new TH2F("deltatcovy", "deltatcovy", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  //Pethu layerwise corr effic
 
  double posxcoreff[nlayer][nlayer];
  double posxcoreffCount;

  double posycoreff[nlayer][nlayer];
  double posycoreffCount;
 
  TH2F* h_posxcoreff = new TH2F("posxcoreff", "posxcoreff", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posycoreff = new TH2F("posycoreff", "posycoreff", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posxcoreff_def = new TH2F("posxcoreff_def", "posxcoreff_def", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posycoreff_def = new TH2F("posycoreff_def", "posycoreff_def", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posxcoreff_dif = new TH2F("posxcoreff_dif", "posxcoreff_dif", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* h_posycoreff_dif = new TH2F("posycoreff_dif", "posycoreff_dif", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  
  for(int ij=0;ij<nlayer;ij++){
    for(int jk=0;jk<nlayer;jk++){
      deltatcov[ij][jk]= deltatcov2[ij][jk]=deltatcovy[ij][jk]=0.0;
      deltatCount[ij][jk]=deltatCount2[ij][jk]=deltatCounty[ij][jk]=0;
      deltaposxcov[ij][jk]= deltaposycov[ij][jk]=0.0;
      deltaposxCount[ij][jk]=deltaposyCount[ij][jk]=0;
      posxcoreff[ij][jk]= posycoreff[ij][jk]=0.0;
    }
  }
  posxcoreffCount=0;posycoreffCount=0;
  //correlation of extrapolated position and time delay
  
  TH2F* xpos_xdev[nlayer][nmxtimehit];
  TH2F* xpos_ydev[nlayer][nmxtimehit];  
  TH2F* ypos_xdev[nlayer][nmxtimehit];
  TH2F* ypos_ydev[nlayer][nmxtimehit];  
  TH2F* xpos_xtdev[nlayer][nmxtimehit];
  TH2F* xpos_ytdev[nlayer][nmxtimehit];  
  TH2F* ypos_xtdev[nlayer][nmxtimehit];
  TH2F* ypos_ytdev[nlayer][nmxtimehit];  

  TH2F* xpos_xtdev_str[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_ytdev_str[nlayer][nmxtimehit];

  TH2F* xpos_xtdev_glb[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_ytdev_glb[nlayer][nmxtimehit];

  TH2F* xpos_ytdev_str[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_xtdev_str[nlayer][nmxtimehit];

  TH2F* xpos_ytdev_glb[nlayer][nmxtimehit]; //Without strip position correction
  TH2F* ypos_xtdev_glb[nlayer][nmxtimehit];

  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nmxtimehit; jk++) {
      sprintf(name, "xpos_xdev_l%i_mul%i", ij, jk+1);
      xpos_xdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);

      sprintf(name, "xpos_ydev_l%i_mul%i", ij, jk+1);
      xpos_ydev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);

      sprintf(name, "ypos_xdev_l%i_mul%i", ij, jk+1);
      ypos_xdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);

      sprintf(name, "ypos_ydev_l%i_mul%i", ij, jk+1);
      ypos_ydev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -2.0, 2.0);

      sprintf(name, "xpos_xtdev_l%i_mul%i", ij, jk+1);
      xpos_xtdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "xpos_ytdev_l%i_mul%i", ij, jk+1);
      xpos_ytdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "ypos_xtdev_l%i_mul%i", ij, jk+1);
      ypos_xtdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "ypos_ytdev_l%i_mul%i", ij, jk+1);
      ypos_ytdev[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

     	 sprintf(name, "xpos_xtdev_str_l%i_mul%i", ij, jk+1);
      xpos_xtdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);
      
      sprintf(name, "ypos_ytdev_str_l%i_mul%i", ij, jk+1);
      ypos_ytdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "xpos_xtdev_glb_l%i_mul%i", ij, jk+1);
      xpos_xtdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);
      
      sprintf(name, "ypos_ytdev_glb_l%i_mul%i", ij, jk+1);
      ypos_ytdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);

      sprintf(name, "xpos_ytdev_str_l%i_mul%i", ij, jk+1);
      xpos_ytdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);
      
      sprintf(name, "ypos_xtdev_str_l%i_mul%i", ij, jk+1);
      ypos_xtdev_str[ij][jk] = new TH2F(name, name, 60, -0.5, 0.5, 60, -10.0, 10.0);

      sprintf(name, "xpos_ytdev_glb_l%i_mul%i", ij, jk+1);
      xpos_ytdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);
      
      sprintf(name, "ypos_xtdev_glb_l%i_mul%i", ij, jk+1);
      ypos_xtdev_glb[ij][jk] = new TH2F(name, name, nstrip+2, -1.5, nstrip+0.5, 60, -15.0, 15.0);
    }
  }

  TH1F* pos_xslope[nlayerit+1][nmxiter];
  TH1F* pos_yslope[nlayerit+1][nmxiter];

  TH1F* pos_theta[nlayerit+1][nmxiter];
  TH1F* pos_phi[nlayerit+1][nmxiter];
  
  for (int ij=0; ij<nlayerit+1; ij++) {
    for (int jk=0; jk<nmxiter; jk++) {
      sprintf(name, "pos_xslope_l%i_i%i", ij, jk);
      pos_xslope[ij][jk] = new TH1F(name, name, 120, -3.0, 3.0);

      sprintf(name, "pos_yslope_l%i_i%i", ij, jk);
      pos_yslope[ij][jk] = new TH1F(name, name, 120, -3.0, 3.0);

      sprintf(name, "pos_theta_l%i_i%i", ij, jk);
      pos_theta[ij][jk] = new TH1F(name, name, 120, 0.0, 75.0);

      sprintf(name, "pos_phi_l%i_i%i", ij, jk);
      pos_phi[ij][jk] = new TH1F(name, name, 120, -180.0, 180.0);

    }
  }
     
  TH2F* widthx_timex[nlayer][nmxiter]; //Not yet implemented in data. Width of pulse hegight vs time delay
  TH2F* widthy_timey[nlayer][nmxiter];
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nmxiter; jk++) {
      sprintf(name, "widthx_timex_l%i_i%i", ij, jk);
      widthx_timex[ij][jk] = new TH2F(name, name, 60., -20., 10., 90, 100., 1000.);
      sprintf(name, "widthy_timey_l%i_i%i", ij, jk);
      widthy_timey[ij][jk] = new TH2F(name, name, 60., -20., 10., 90, 100., 1000.);
    }
  }
  
  TH2F* sel_theta_phi_g4 = new TH2F("sel_theta_phi_g4", "sel_theta_phi_g4", 180,-180.,180.0, 90, 0.0, 90.0);
  TH2F* sel_theta_phi_g4_trig = new TH2F("sel_theta_phi_g4_trig", "sel_theta_phi_g4_trig", 180,-180.,180.0, 90, 0.0, 90.0);
  TH2F* sel_theta_mom = new TH2F("sel_theta_mom","sel_theta_mom",240,0.,60.,1000,0.150,3.);
  TH2F* sel_phi_mom = new TH2F("sel_phi_mom","sel_phi_mom",180,-180.,180.0,1000,0.150,3.);
  TH2F* sel_theta_rec_gen = new TH2F("sel_theta_rec_gen","sel_theta_rec_gen",90,0.,90.0,90,0.,90.0);
  
  const int nCount=100000; //0; //Number of events //Time is second
  const int nsetmx=2000; //Maximum number of set
  int iset=0;
  int isetold=-1;
  TH2F* strp_count_set[nsetmx]={0};
  TH1F* strp_xmult_set[nlayer][nsetmx]={0};
  TH1F* strp_ymult_set[nlayer][nsetmx]={0};
  TH2F* raw_occu_set[nlayer][nsetmx]={0};

  TH1F* xlayer_reso_set[nlayer][nsetmx]={0};
  TH1F* ylayer_reso_set[nlayer][nsetmx]={0}; 
  TH1F* time_xreso_set[nlayer][nsetmx]={0};
  TH1F* time_yreso_set[nlayer][nsetmx]={0};  

#ifdef MONTECARLO
#include "effic_pos_time_madurai64_aug252627.txt"
#include "effic_trig_madurai64_aug252627.txt"
#include "time_corr_madurai64_aug252627.txt"
#include "pos_time_inuse_madurai64_aug252627.txt"
  
#elif defined(ONEMBYONEM)
  
  //#include "effic_trig_160707.txt" //Using trigger on 0-1-3-4 & 7-8-10-11
  //#include "effi_postime_160707.txt" //Using trigger on 0-1-3-4 & 7-8-10-11
  
  //#include "pos_time_inuse_160711.txt" //using only 2479 trigger data (otherwise wider time 
                                     //resolution paticularly in few strip of Layer-1 & 2 in X-side)
#include "time_corr_160711.txt"
  
#else
  
   //   //#include "pos_time_inuse_madurai64_160316.txt" //pos_time_inuse_150705.txt" //"pos_time_inuse_150228.txt"
#endif

// 051116 changed ---only in MC
#ifdef MONTECARLO
#ifdef FASTSIM
for(int ij =0;ij<nmxhits;ij++) {
	for(int jk=0;jk<nlayer;jk++) {
		xposerrsq[ij][jk] *= 1.1015;//1.0358; 1.0815;
		yposerrsq[ij][jk] *= 1.1207;//1.0490; 1.105;
	}
} 
#endif
#ifdef FULLG4SIM
for(int ij =0;ij<nmxhits;ij++) {
	for(int jk=0;jk<nlayer;jk++) {
	  xposerrsq[ij][jk] *= 1.;//1.647/1.556;//(1.696/1.556);//1.500;//00.5000;//1.200;//0.900;//1.1015;//1.0358; 1.0815;
	  yposerrsq[ij][jk] *= 1.;//1.64/1.472;//;//(1.705/1.472);//1.500;//0.5000;//1.200;//0.900;//1.1207;//1.0490; 1.105;
	}
} 
#endif
#endif 

#ifdef ISEFFICIENCY
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<nstrip; jk++) {
      for (int kl=0; kl<nstrip; kl++) {
        defefficiency_uncx[ij]->Fill(jk, kl, ineffi_table_uncX[ij][jk][kl]);
        defefficiency_uncy[ij]->Fill(jk, kl, ineffi_table_uncY[ij][jk][kl]);
        deftriggereffi_x[ij]->Fill(jk, kl, effi_trig_X[ij][jk][kl]);
        deftriggereffi_y[ij]->Fill(jk, kl, effi_trig_Y[ij][jk][kl]);
      }
    }
  }

#endif
//-------------------------------------------------------------------------
  
  if (isalign >0) {
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<nstrip; jk++) {
        for (int kl=0; kl<nstrip; kl++) {
        }
      }
    }
  }

  //nonlinear shift in x/y timing with the number in in y/x strip  
  //  double xtoffystr[nlayer][nstrip]={0};
  //  double ytoffxstr[nlayer][nstrip]={0};

  double xtrms[nlayer][nstrip]={0};
  double ytrms[nlayer][nstrip]={0};
  
  double timesx[nlayer]={0};
  double timesy[nlayer]={0};

  int cau_calib1[nlayer] = {0};
  int cau_calib2[nlayer] = {0};
  
  double widthx[nlayer]={0};
  double widthy[nlayer]={0};

  //Position resolution in layers. 150219
  double pos_xrms[nlayer][nmxhits];
  double pos_yrms[nlayer][nmxhits];

  //Timing resolution in layers. 150112
  double time_xrms[nlayer];
  double time_yrms[nlayer];

  for (int ij=0; ij<nlayer; ij++) {
    time_xrms[ij] = time_yrms[ij] = 1.1; //ns
    for (int jk=0; jk<nmxhits; jk++) {
      if (jk<=1) { 
        pos_xrms[ij][jk] = pos_yrms[ij][jk] =0.25;
      } else {
        pos_xrms[ij][jk] = pos_yrms[ij][jk] =0.30;
      }
    }
  }
  
  bool filloccu = true; // For first entries fill occuplancy plot and remove noisy channels
  
  double errcst, errcov, errlin;
  double dist[nlayer], xtime[nlayer],ytime[nlayer],  xtdev[nlayer],ytdev[nlayer];
  double initxtime[nlayer][nTDCpLayer], rawxtime[nlayer], rawxtime0[nlayer], rawytime0[nlayer], rawxtime1[nlayer], rawytime1[nlayer], rawxtime2[nlayer], rawytime2[nlayer];
  double initytime[nlayer][nTDCpLayer], rawytime[nlayer], rawxtime3[nlayer], rawytime3[nlayer]; //150628 : Time without path lenght correction, but with all other correction

  bool xusedtime[nlayer], yusedtime[nlayer];

  double xval, yval; 
  int nxtfail, nytfail, nentry, nTotalp, nTotalt,nTotallx,nTotally;
  double xc0inters,yc0inters,DDx,errcst_tx,errcov_tx,errlin_tx,DDy,errcst_ty,errcov_ty,errlin_ty;
  int firstiter,ntcormx, iiterrs, lstr, lend,occulyr ,isfill,ntrigX,ntrigY;

  double xtext[nlayer]; //, xtexter[nlayer];
  double ytext[nlayer]; //, ytexter[nlayer];

  //*********************************************************************************************
  //               Iteration Starts
  //*********************************************************************************************
      
  nTotalp = nTotalt = 0;
  nTotallx=nTotally=0;
  

  const int timebin=60; //60 second bin
  int tmptimerate=-1;
  int tmpoldtime=-1;
  int tmpoldmuposxrate=-1;
  int tmpoldmuposyrate=-1;
  int tmpoldmutimexrate=-1;
  int tmpoldmutimeyrate=-1;

  int ntimerate = 0;  //counter for number of events in a second
  int nmuposxrate = 0;
  int nmuposyrate = 0;
  int nmutimexrate = 0;
  int nmutimeyrate = 0;

  int nalltimerate = 0; //counter for seconds
  int nallmuposxrate = 0;
  int nallmuposyrate = 0;
  int nallmutimexrate = 0;
  int nallmutimeyrate = 0;  
  
  //Total time of data taken in the file
  double init_time=0,end_time=0,total_time=0;
  double start_time =0.;
  double time_diff = 0.;
  double tmptime_diff = 0.;
  TH1F* Time_diff_evt = new TH1F("event_time_diff","event_time_diff",1000,0.,100000.);

  for (int iiter=0; iiter<nmxiter; iiter++) {
    
    ntotxtime = totxtime = ntotytime = totytime = 0;
    
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<nstrip; jk++) {
        xtrms[ij][jk]= ytrms[ij][jk]=0.0;
      }
    }
    
    firstiter=0; // before fit, calculate time shift in each layer with all hists
    
    //    int ntcormx = 2;
    //Check this condition properly, why should we use this
    ntcormx =(layfirst-firstlayer>0 || lastlayer-laylast>0) ? 1 : 2;
    if (isOnlyCom) ntcormx=1;
    for (int ntcor=0; ntcor< ntcormx; ntcor++) {
      iiterrs=nmxiter*ntcor+iiter;
      
      lstr = max(firstlayer,0);  //Always 0
      //      lstr = max(firstlayer,2);  //GMA140130

      lend = min(lastlayer+1,nlayerit); //nlayerit=nlayer=12     //Always 12
      
      if (ntcor==0) { lend = lstr+1;}
      
      file_out <<"lay1 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<endl;
      
      for (int laye= lstr; laye<lend; laye++) {
        occulyr = (ntcor==1) ? isequence[abs(laye)] : nlayer;
        //	if (ntcor==1 && (occulyr!=8)) continue;
        //        if (ntcor==1 && (occulyr==0||occulyr==1||occulyr==2||occulyr==3||occulyr==4 || occulyr==5||occulyr==7||occulyr==9||occulyr==10||occulyr==11)) continue;
        isfill = (iiter==nmxiter-1 && ntcor==0) ? true : false; //Fill rootuple and other histogramme for the final iteration with all layers
        
        ifstream file_db;
        file_db.open(rootfiles);  
        
        file_out <<"lay2 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<" "<<laye<<" "<<occulyr<<" "<<ytimeshift<<endl;
        cout <<"lay2 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<" "<<laye<<" "<<occulyr<<" "<<ytimeshift<<endl;
        
        int nfile=0; //Initialisation for new histogramme booking
        while(!(file_db.eof())) {
          file_db >> datafile>>nentrymx;
          if (strstr(datafile,"#")) continue;
          if(file_db.eof()) break;
          nfile++;
          int nx4 =0,nx5=0,nx6=0,nx7=0,ny4=0,ny5=0,ny6=0,ny7=0,nxy4=0,nxy5=0,nxy6=0,nxy7=0,nxory4=0,nxory5=0,nxory6=0,nxory7=0, nxytrig=0;
          int ntotxext[12][8] ={0};
          int ntotyext[12][8] ={0};
          int ntxtime[12][8] = {0};
          int ntytime[12][8] ={0};
          
          int ntotxext2[12][8] ={0};
          int ntotyext2[12][8] ={0};
          int ntxtime2[12][8] = {0};
          int ntytime2[12][8] ={0};
          
#ifdef MONTECARLO
// #ifdef FASTSIM	  
          
//           sprintf(infile, "/media/pethuraj/9CFCFF46FCFF196A/Madurai_MC_64_130916/%s", datafile); //Standalone Monte carlo
//           //	 sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_Geant4_Data/%s", datafile);
//           //	  sprintf(infile, "/home/pethuraj/Pictures/RPCStackSim27122016_BackMuon/%s", datafile);
//           //  sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_East_West_Ayymetry_Data/%s",datafile);
//           //  sprintf(infile, "/media/9CFCFF46FCFF196A/Muon_BackScatt_DATA/%s",datafile);
//           TFile *fileIn = new TFile(infile, "read");
//           if (isfill) { 
//             fileOut->cd();
//             for (int ij=0; ij<10; ij++) {
//               sprintf(name, "sel_theta_phi_%i", ij);
//               tmp_sel_theta_phi[ij] = (TH2F*)fileIn->Get(name);
//               if (nfile==1) { 
//                 sel_theta_phi[ij] = new TH2F(name, 
//                                              tmp_sel_theta_phi[ij]->GetTitle(), 
//                                              tmp_sel_theta_phi[ij]->GetNbinsX(), 
//                                              tmp_sel_theta_phi[ij]->GetXaxis()->GetXmin(),  
//                                              tmp_sel_theta_phi[ij]->GetXaxis()->GetXmax(),
//                                              tmp_sel_theta_phi[ij]->GetNbinsY(), 
//                                              tmp_sel_theta_phi[ij]->GetYaxis()->GetXmin(),  
//                                              tmp_sel_theta_phi[ij]->GetYaxis()->GetXmax()); 
//               }
//               sel_theta_phi[ij]->Add(tmp_sel_theta_phi[ij]);
//             }
//             fileIn->cd();
//           }
          
//           TTree *T1= (TTree*)fileIn->Get("T1");
          
//           T1->SetBranchAddress("thgen", &thgen);
//           T1->SetBranchAddress("phgen", &phgen);
//           T1->SetBranchAddress("xpgen", &xpgen);
//           T1->SetBranchAddress("ypgen", &ypgen);
          
//           T1->SetBranchAddress("xystrp", xystrp);
//           T1->SetBranchAddress("xstrp", xstrp);
//           T1->SetBranchAddress("ystrp", ystrp);
//           T1->SetBranchAddress("xtimea", xtimea);
//           T1->SetBranchAddress("ytimea", ytimea);
          
//           nentry = T1->GetEntries();
//           nentry = min(nentry,nentrymx);
          
//           cout <<"MC file "<< datafile<<" entries "<<nentry<<" iiter= "<<iiter<<" Occuly= "<<occulyr<<endl; 
//           for(int iev=0;iev<nentry;iev++) { 
            
//             xslope= xinters= yslope= yinters= timexslope= timeyslope= timex2slope= timey2slope=-100.;
//             xchi2= ychi2= xt0chi2= yt0chi2=-100.;
//             nxstrip= nystrip=Nx=Ny=nxtime=nytime=ntxyla=0;
            
//             for (int ix=0; ix<nlayer; ix++) {
//               for (int iy=0; iy<2*nmxiter; iy++) {
//                 istime_xyreso[ix][iy]=false;
//                 is_xyreso[ix][iy]=false;
//                 xposEffPass[ix][iy] = false;
//                 yposEffPass[ix][iy] = false;
//               }
//             }
            
//             fileIn->cd();
//             T1->GetEntry(iev); 
//             //
//             vector<int> xpts[nlayer];
//             vector<int> ypts[nlayer];
            
//             vector<int> xptsall[nlayer];
//             vector<int> yptsall[nlayer];	    
            
//             int tmpnx=0;
//             int xx[nlayer][nmxusedhits]= {-100};
//             int yy[nlayer][nmxusedhits]= {-100};
//             int xmultcount[nlayer]={0};
//             int ymultcount[nlayer]={0};
//             for (int ij=0; ij<nlayer; ij++) {
//               // Xpos[ij] = Ypos[ij] = -100;
//               for(int xxx=0;xxx<nmxusedhits;xxx++) {
//                 xx[ij][xxx]= -100;
//                 yy[ij][xxx]= -100;
//               }
              
//               xmultcount[ij]=0;
//               ymultcount[ij]=0;
//               int nxmul = xystrp[ij]%10;
//               int nymul = int(xystrp[ij]/10)%10;
//               int ix = int(xystrp[ij]/1000)%100;
//               int iy = int(xystrp[ij]/100000);
//               for(int nml=0;nml<nmxusedhits;nml++) {
//                 xx[ij][nml] = xstrp[ij][nml];if(xx[ij][nml]>-100) {xmultcount[ij]++;} 
//                 yy[ij][nml] = ystrp[ij][nml];if(yy[ij][nml]>-100) {ymultcount[ij]++;} 
                
//               }
//               //cout<<"mult"<<"	"<<xmultcount<<"	"<<ymultcount<<endl;
//               if (nxmul>0) { tmpnx++;}
//               float tmpxpos=0;
//               for (int jk=0; jk<xmultcount[ij]; jk++) {
//                 if (xx[ij][jk]>=nstrip) {
//                   cout <<"Error in MC sample, X-strip id = "<<iev<<" "<<jk<<"+"<<ix<<" "<<xystrp[ij]<<endl;
//                 } else {
//                   //   if(( ij!=5 && ij!=6 && ij!=7  && ij!=8) ||  xx[ij][jk] <29) {
//                   // if(ij!= 4) {
//                   xptsall[ij].push_back(xx[ij][jk]);
//                   xpts[ij].push_back(xx[ij][jk]);
//                   xlayer_occu[ij]->Fill(xptsall[ij][jk]);
//                   tmpxpos +=(xx[ij][jk]);
//                 }
//               }
              
//               tmpxpos = tmpxpos/max(1, xmultcount[ij]) + 0.5 - xoff[ij];
              
//               float tmpypos = 0;
//               for (int jk=0; jk<ymultcount[ij]; jk++) {
//                 if (yy[ij][jk]>=nstrip) { 
//                   cout <<"Error in MC sample, Y-strip id = "<<iev<<" "<<jk<<"+"<<iy<<" "<<xystrp[ij]<<endl;
//                 } else {
//                   // if((ij!=5 && ij!=6 && ij!=7  && ij!=8) ||  yy[ij][jk] <29) {
//                   //  if( ij !=4 ) {
//                   yptsall[ij].push_back(yy[ij][jk]);
//                   ypts[ij].push_back(yy[ij][jk]);
//                   ylayer_occu[ij]->Fill(yptsall[ij][jk]);
//                   tmpypos +=(yy[ij][jk]);
//                   //}
//                   //}
//                 }
//               }
//               tmpypos = tmpypos/max(1, ymultcount[ij]) + 0.5 - yoff[ij];
              
//               xhits[ij] = xpts[ij].size();
//               yhits[ij] = ypts[ij].size();

//               if (tmpxpos >lastXstrip) { xhits[ij] = 0;}
//               if (tmpypos >lastYstrip) { yhits[ij] = 0;}
//               if (nxmul>0 && nxmul<4 && xxerr[ij]<=0.001) cout <<"xij "<<iiterrs<<" "<<ij<<" "<< nxmul<<" "<<nymul<<" "<<xposerrsq[nxmul-1][ij]<<" "<<xxerr[ij]<<endl;
//               if (nymul>0 && nymul<4 && yyerr[ij]<=0.001) cout <<"yij "<<iiterrs<<" "<<ij<<" "<< nxmul<<" "<<nymul<<" "<<yposerrsq[nymul-1][ij]<<" "<<yyerr[ij]<<endl;
//             }
//             h_tmp1xndf->Fill(tmpnx);
// #else	//FULLG4SIM  
// 	    //  sprintf(infile, "/home/pethuraj/Documents/RPCStackSim27122016/%s",datafile);
//             sprintf(infile, "/home/pethuraj/Documents/RPCStackSim20170906/%s",datafile);
            
//             // sprintf(infile, "/media/9CFCFF46FCFF196A/Madurai_Geant4_Data/%s", datafile); //media/pethuraj/9CFCFF46FCFF196A/Madurai_Geant4_Data/
//             TFile *fileIn = new TFile(infile, "read");
//             if (isfill) { 
//               fileOut->cd();
//               for (int ij=0; ij<10; ij++) {
//                 sprintf(name, "sel_theta_phi_%i", ij);
//                 tmp_sel_theta_phi[ij] = (TH2F*)fileIn->Get(name);
//                 if (nfile==1) { 
//                   sel_theta_phi[ij] = new TH2F(name, 
//                                                tmp_sel_theta_phi[ij]->GetTitle(), 
//                                                tmp_sel_theta_phi[ij]->GetNbinsX(), 
//                                                tmp_sel_theta_phi[ij]->GetXaxis()->GetXmin(),  
//                                                tmp_sel_theta_phi[ij]->GetXaxis()->GetXmax(),
//                                                tmp_sel_theta_phi[ij]->GetNbinsY(), 
//                                                tmp_sel_theta_phi[ij]->GetYaxis()->GetXmin(),  
//                                                tmp_sel_theta_phi[ij]->GetYaxis()->GetXmax()); 
//                 }
//                 sel_theta_phi[ij]->Add(tmp_sel_theta_phi[ij]);
                
//                 for(int jk=0;jk<8;jk++) {
//                   sprintf(name, "sel_theta_phi_%i_%i", ij,jk);
//                   tmp_sel_theta_phi_8[ij][jk] = (TH2F*)fileIn->Get(name);
//                   if (nfile==1) { 
//                     sel_theta_phi_8[ij][jk] = new TH2F(name, 
//                                                        tmp_sel_theta_phi_8[ij][jk]->GetTitle(), 
//                                                        tmp_sel_theta_phi_8[ij][jk]->GetNbinsX(), 
//                                                        tmp_sel_theta_phi_8[ij][jk]->GetXaxis()->GetXmin(),  
//                                                        tmp_sel_theta_phi_8[ij][jk]->GetXaxis()->GetXmax(),
//                                                        tmp_sel_theta_phi_8[ij][jk]->GetNbinsY(), 
//                                                        tmp_sel_theta_phi_8[ij][jk]->GetYaxis()->GetXmin(),  
//                                                        tmp_sel_theta_phi_8[ij][jk]->GetYaxis()->GetXmax()); 
//                   }
//                   sel_theta_phi_8[ij][jk]->Add(tmp_sel_theta_phi_8[ij][jk]);
                  
//                 }
                
//                 for(int jk=0;jk<16;jk++) {
//                   sprintf(name, "sel_theta_phi_16_%i_%i", ij,jk);
//                   tmp_sel_theta_phi_16[ij][jk] = (TH2F*)fileIn->Get(name);
//                   if (nfile==1) { 
//                     sel_theta_phi_16[ij][jk] = new TH2F(name, 
//                                                         tmp_sel_theta_phi_16[ij][jk]->GetTitle(), 
//                                                         tmp_sel_theta_phi_16[ij][jk]->GetNbinsX(), 
//                                                         tmp_sel_theta_phi_16[ij][jk]->GetXaxis()->GetXmin(),  
//                                                         tmp_sel_theta_phi_16[ij][jk]->GetXaxis()->GetXmax(),
//                                                         tmp_sel_theta_phi_16[ij][jk]->GetNbinsY(), 
//                                                         tmp_sel_theta_phi_16[ij][jk]->GetYaxis()->GetXmin(),  
//                                                         tmp_sel_theta_phi_16[ij][jk]->GetYaxis()->GetXmax()); 
//                   }
//                   sel_theta_phi_16[ij][jk]->Add(tmp_sel_theta_phi_16[ij][jk]);
//                 }
                
//               }
//               fileIn->cd();
//             }
//             //  cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
     
//             TTree *T1= (TTree*)fileIn->Get("T1");
            
//             T1->SetBranchAddress("irun", &irun);
//             T1->SetBranchAddress("ievt", &ievt);
            
//             T1->SetBranchAddress("ngent", &ngent);
//             T1->SetBranchAddress("pidin", pidin);
//             T1->SetBranchAddress("ievt_wt", &ievt_wt);
//             T1->SetBranchAddress("intxn_id", &intxn_id);
//             T1->SetBranchAddress("momin", momin);
//             T1->SetBranchAddress("thein", thein);
//             T1->SetBranchAddress("thegen", thegen);
//             T1->SetBranchAddress("phiin", phiin);
//             T1->SetBranchAddress("posxin", posxin);
//             T1->SetBranchAddress("posyin", posyin);
//             T1->SetBranchAddress("poszin", poszin);
//             T1->SetBranchAddress("nsimhtx", &nsimhtx);
//             T1->SetBranchAddress("nsimhty", &nsimhty);
//             T1->SetBranchAddress("nlayert", &nlayert);
//             T1->SetBranchAddress("xdata", xdata);
//             T1->SetBranchAddress("ydata", ydata);
//             T1->SetBranchAddress("triggerinfoX", &triggerinfoX);
//             T1->SetBranchAddress("triggerinfoY", &triggerinfoY);
            
//             if (isTiming) {
//               T1->SetBranchAddress("digiXtime", digiXtime);
//               T1->SetBranchAddress("digiYtime", digiYtime);
//             }	  
            
//             nentry = T1->GetEntries();
//             nentry = min(nentry,nentrymx);
//             // if (filloccu) {ntotal += nentry;}
//             cout<<"nentry "<< nentry<<endl;
//             for(int iev=0;iev<nentry;iev++) { 
//               xslope= xinters= yslope= yinters= timexslope= timeyslope= timex2slope= timey2slope=-100.;
//               xchi2= ychi2= xt0chi2= yt0chi2=-100.;
//               nxstrip= nystrip=Nx=Ny=nxtime=nytime=ntxyla=0;
              
//               for (int ix=0; ix<nlayer; ix++) {
//                 for (int iy=0; iy<2*nmxiter; iy++) {
//                   istime_xyreso[ix][iy]=false;
//                   is_xyreso[ix][iy]=false;
//                   xposEffPass[ix][iy] = false;
//                   yposEffPass[ix][iy] = false;
//                 }
//               }
              
//               fileIn->cd();
              
//               //T1->Show(iev);
//               T1->GetEntry(iev); 
//               if (isfill) { 
//                 fileOut->cd();
//                 sel_theta_phi_g4->Fill(180.*(phiin[0]/pival),180.-180.*(thein[0]/pival));
//               }
              
//               bool trginfo = false;
//               trginfo = (triggerinfoX==15 || triggerinfoY==15) ? 1:0;
//               //	if(triggerinfoX!=15 && triggerinfoY!=15) continue;
//               if(!trginfo) continue;
//               ntotal++;
//               if (isfill) { 
//                 fileOut->cd();
//                 sel_theta_phi_g4_trig->Fill(180.*(phiin[0]/pival),180.-180.*(thein[0]/pival));
//               }
//               for (int ij=0; ij<nlayer; ij++) { Xdev[ij] = 100; Xpos[ij]=0.0;}
//               for (int ij=0; ij<nlayer; ij++) { Ydev[ij] = 100; Ypos[ij]=0.0;}
//               vector<int> xpts[nlayer];
//               vector<int> ypts[nlayer];
              
//               vector<int> xptsall[nlayer];
//               vector<int> yptsall[nlayer];	    
              
//               for (int ij=0; ij<nlayer; ij++) { 
//                 for (int jk=0; jk<nstrip; jk++) {
//                   if ((xdata[ij]>>jk)&0x01) { 
//                     xptsall[ij].push_back(jk);
//                   }
//                   if ((ydata[ij]>>jk)&0x01) { 
//                     yptsall[ij].push_back(jk);
//                   }
//                   if (isTiming) {
//                     xtime[ij] = digiXtime[ij]*0.1;
//                     ytime[ij] = digiYtime[ij]*0.1;
//                   }
//                 }
//               }
	    
//               xslope = 0;
//               xinters = 0;
//               yslope = 0;
//               yinters = 0;
//               ntrigX = ntrigY = 0;
//               for (int iz=0; iz<nlayer; iz++) {
//                 xhits[iz] = yhits[iz] =0;
//                 //Don't use it here	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && xptsall[iz].size()>0) {ntrigX++;} //18 01 2012 
//                 if (filloccu) {
//                   xlayer_allmult[iz]->Fill(xptsall[iz].size());
                  
//                 }
//                 for (int ix=0; ix<xptsall[iz].size(); ix++) {
//                   if (filloccu) {
//                     xlayer_alloccu->Fill(nstrip*iz+xptsall[iz][ix]);
//                     xlayer_occu[iz]->Fill(xptsall[iz][ix]);
                    
//                   }
//                   bool failed=false;
//                   //	if (!posinuse[0][iz][xptsall[iz][ix]]) {failed=true;}
//                   //		if(iz==11 && xptsall[iz][ix]>=0 && xptsall[iz][ix]<=7)  {failed=true;}
//                   if (!failed){
//                     xpts[iz].push_back(xptsall[iz][ix]);
//                     if (filloccu) {xlayer_alloccusel->Fill(nstrip*iz+xptsall[iz][ix]);}
//                   }
//                 }
//                 xhits[iz] = xpts[iz].size();
//                 if (filloccu) {
//                   ylayer_allmult[iz]->Fill(yptsall[iz].size());
//                 }
//                 for (int iy=0; iy<yptsall[iz].size(); iy++) {
//                   // cout<<"yptsall"<<"	"<<iz<<"	"<<iy<<"	"<<yptsall[iz][iy]<<endl;
//                   if (filloccu) {
//                     ylayer_occu[iz]->Fill(yptsall[iz][iy]);
//                     ylayer_alloccu->Fill(nstrip*iz+yptsall[iz][iy]);
//                   }
//                   bool failed=false;
//                   //	if (!posinuse[1][iz][yptsall[iz][iy]]) {failed=true;}
//                   if (!failed){
//                     ypts[iz].push_back(yptsall[iz][iy]);
//                     if (filloccu) {ylayer_alloccusel->Fill(nstrip*iz+yptsall[iz][iy]);}
//                   }
//                 }
//                 yhits[iz] = ypts[iz].size();
//               }
              
//               if(filloccu){
//                 for (int ij=0; ij<nlayer; ij++) {
//                   xlayer_mult[ij]->Fill(xpts[ij].size());
//                   ylayer_mult[ij]->Fill(ypts[ij].size());
                  
//                   rawhits_corr_xymul[ij]->Fill(xptsall[ij].size(), yptsall[ij].size());
                  
//                   for (int ix=0; ix<xptsall[ij].size(); ix++) {
//                     for (int iy=0; iy<yptsall[ij].size(); iy++) {
//                       raw_occu[ij]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
//                     }
//                   }
//                   for (int jk=ij+1; jk<nlayer; jk++) {
//                     rawhits_xlay_corr_mul[ij][jk]->Fill(xptsall[ij].size(), xptsall[jk].size());
//                     rawhits_ylay_corr_mul[ij][jk]->Fill(yptsall[ij].size(), yptsall[jk].size());
//                   }
//                 }
//               }
// #endif	 
              
#else
              //	  sprintf(infile, "hrdc217/%s", datafile);
          sprintf(infile, "newformat/%s", datafile);
          TFile *fileIn = new TFile(infile, "read");
          TTree *event_tree=(TTree*)fileIn->Get("evetree");
          EveTree *event = new EveTree(event_tree);
          event->Loop();
          
          nentry = event_tree->GetEntries();
          nentry = min(nentry,nentrymx);
          
          cout <<infile<<" has "<< nentry<<" events "<<endl;
          for(int iev=0; iev<nentry; iev++) {    //ij is event loop variable. while checking put ij<3 or a small number.
            fileIn->cd();
            event_tree->GetEntry(iev);    // while running for entire event file put "ij<numentries".
            if (filloccu) {
              if (ntotal==0) {nsec=event->Evetime->GetSec(); cout<<"nsec"<<"	"<<nsec<<endl;}
              if(iev == 1) { init_time = event->Evetime->GetSec();/*(event->Evetime->AsGMST()*3600.)*/; }
              if(iev == nentry-1) { end_time = event->Evetime->GetSec();/* event->Evetime->AsGMST()*3600.;*/ 
                total_time +=  end_time -init_time ;
              }
              
              int tmptime = tmptimerate = event->Evetime->GetSec(); //-nsec;
              tmptime_diff = event->Evetime->AsGMST()*3600.*1000000.-start_time;
              start_time = event->Evetime->AsGMST()*3600.*1000000.;
              if(isfill) {
                Time_diff_evt->Fill(tmptime_diff);
              }
              
              if (iev>0 && tmptimerate >=tmpoldtime+timebin) {
                file_out <<"totrate "<<datafile<<" "<<iev<<" "<<nalltimerate<<" "<<ntimerate<<" "<<tmptimerate<<endl;
                trigrate->Fill(ntimerate*(1.0/timebin));
                trigratex->Fill(nalltimerate, ntimerate*(1.0/timebin));
                ntimerate = 0;
                nalltimerate++;
                tmpoldtime = tmptimerate;
              }
              ntimerate++;
              
              ntotal++;
              //	      iset = int(tmptime/nCount);
              iset = int(ntotal/nCount);
              if (iset >=nsetmx) { iset = nsetmx-1;}
              if (isetold != iset) {
                fileOut->cd();
                TDatime T0x(tmptime);
                sprintf(name, "strp_count_set_%i", iset);
                sprintf(title, "strp_count_set_%i [%s]", iset, T0x.AsString());
                strp_count_set[iset] = new TH2F(name, title, 2*nlayer, -0.5,  2*nlayer-0.5, nstrip, -0.5, nstrip-0.5);
                for (int ix=0; ix<nlayer; ix++) {
                  sprintf(name, "strp_xmult_set_l%i_%i", ix, iset);
                  sprintf(title, "strp_xmult_set_l%i_%i [%s]", ix, iset, T0x.AsString());
                  strp_xmult_set[ix][iset] = new TH1F(name, title, nstrip+1, -0.5, nstrip+0.5);
                  
                  sprintf(name, "strp_ymult_set_l%i_%i", ix, iset);
                  sprintf(title, "strp_ymult_set_l%i_%i [%s]", ix, iset, T0x.AsString());
                  strp_ymult_set[ix][iset] = new TH1F(name, title, nstrip+1, -0.5, nstrip+0.5);
                  
                  sprintf(name, "raw_occu_set_l%i_%i", ix, iset);
                  sprintf(title, "raw_occu_set_l%i_%i [%s]", ix, iset, T0x.AsString());
                  raw_occu_set[ix][iset] = new TH2F(name, title, nstrip, -0.5, nstrip-0.5, nstrip, -0.5, nstrip-0.5);
                  
                  sprintf(name, "xlayer_reso_set_l%i_i%i", ix, iset);
                  sprintf(title, "xlayer_reso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
                  xlayer_reso_set[ix][iset]=new TH1F(name, title, 60, -1.8, 1.8);
                  
                  sprintf(name, "ylayer_reso_set_l%i_i%i", ix, iset);
                  sprintf(title, "ylayer_reso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
                  ylayer_reso_set[ix][iset]=new TH1F(name, title, 60, -1.8, 1.8);
                  
                  sprintf(name, "time_xreso_set_l%i_i%i", ix, iset);
                  sprintf(title, "time_xreso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
                  time_xreso_set[ix][iset] = new TH1F(name, title, 60, -7.5, 7.5);
                  
                  sprintf(name, "time_yreso_set_l%i_i%i", ix, iset);
                  sprintf(title, "time_yreso_set_l%i_i%i [%s]", ix, iset, T0x.AsString());
                  time_yreso_set[ix][iset] = new TH1F(name, title, 60, -7.5, 7.5);
                  
                }
                isetold = iset;
                fileIn->cd();
              }
            }
            if (iev%20000==1) {
              cout <<"iev "<< iev<<" "<<event->ENum<<" "<<" "<<event->Evetime->GetDate()<<" "<<event->Evetime->GetTime()<<" "<<event->Evetime->GetSec()<<" "<<event->Evetime->GetSec()<<" "<<event->Evetime->GetTime(0, nsec)<<" "<<nsec<<" "<<event->Evetime->GetSec()-nsec<<endl;   
            }
            
            xslope= xinters= yslope= yinters= timexslope= timeyslope= timex2slope= timey2slope=-100.;
            xchi2= ychi2= xt0chi2= yt0chi2=-100.;
            nxstrip= nystrip=Nx=Ny=nxtime=nytime=ntxyla=0;
            
            for (int ix=0; ix<nlayer; ix++) {
              for (int iy=0; iy<2*nmxiter; iy++) {
                istime_xyreso[ix][iy]=false;
                is_xyreso[ix][iy]=false;
                xposEffPass[ix][iy] = false;
                yposEffPass[ix][iy] = false;
              }
            }
            
            vector<int> xpts[nlayer];
            vector<int> ypts[nlayer];
            
            vector<int> xptsall[nlayer]; //For trigger criteria
            vector<int> yptsall[nlayer];	    
            
            fileOut->cd();
            
            vector<int> xptsalltdc[nlayer][nTDCpLayer]; //signal for each TDC
            vector<int> yptsalltdc[nlayer][nTDCpLayer];	    
            for(int jk=0;jk<nlayer;jk++) {
              for (int kl=0; kl<nTDCpLayer; kl++) {
                xallhits[jk][kl] = yallhits[jk][kl] =-1;
                xptsalltdc[jk][kl].clear();  yptsalltdc[jk][kl].clear();
              }
              xpts[jk].clear(); ypts[jk].clear();  xptsall[jk].clear(); yptsall[jk].clear();
              //printf("layer\t");
              for(int kl=0; kl<nstrip; kl++) {
                if(event->xLayer[jk]->TestBitNumber(kl)) {
                  xptsall[jk].push_back(kl);
#ifdef NEWRPCDAQ1
                  xptsalltdc[jk][kl%8].push_back(kl);
#else 
                  xptsalltdc[jk][int((nTDCpLayer*kl)/nstrip)].push_back(kl); //GMA20241126
#endif		  
                }
                if(event->yLayer[jk]->TestBitNumber(kl)) {
                  
                  yptsall[jk].push_back(kl);
#ifdef NEWRPCDAQ1
                  yptsalltdc[jk][kl%8].push_back(kl);
#else 
                  yptsalltdc[jk][int((nTDCpLayer*kl)/nstrip)].push_back(kl);
#endif
                }
              }
            }
            
            // Store total number of hits in X/Y layer irrespective of noise etc.
            for (int iz=0; iz<nlayer; iz++) {
              for (int itdc=0; itdc<nTDCpLayer; itdc++) { 
                xallhits[iz][itdc] = (xptsalltdc[iz][itdc].size()==1) ? xptsalltdc[iz][itdc][0] : -10-xptsalltdc[iz][itdc].size(); //avoid confusion with no hits and only hits at strip# 1
                yallhits[iz][itdc] = (yptsalltdc[iz][itdc].size()==1) ? yptsalltdc[iz][itdc][0] : -10-yptsalltdc[iz][itdc].size();
                
                //	      cout <<"xallpts "<<iz<<" "<<xallhits[iz]<<" "<<xptsall[iz].size()<<" "<<yallhits[iz]<<" "<<yptsall[iz].size()<<endl;
              }
            }
            
            ////////////////////////////////////////////
            //
            //  First clean up noisey layers and then
            //  Strip efficiency and resolutions etc etc
            //
            ////////////////////////////////////////////
            // cout<<iev<<"  "<<"PASSSSSSSSSSSS 2"<<endl;
            xslope = 0;
            xinters = 0;
            yslope = 0;
            yinters = 0;
            ntrigX = ntrigY = 0;
            
            for (int iz=0; iz<nlayer; iz++) {
              xhits[iz] = yhits[iz] =0;
              //Don't use it here	      if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && xptsall[iz].size()>0) {ntrigX++;} //18 01 2012 
              if (filloccu) {
                xlayer_allmult[iz]->Fill(xptsall[iz].size());
#ifndef MONTECARLO
                strp_xmult_set[iz][iset]->Fill(xptsall[iz].size());
#endif		
              }
              for (int ix=0; ix<xptsall[iz].size(); ix++) {
                if (filloccu) {
                  xlayer_alloccu->Fill(nstrip*iz+xptsall[iz][ix]);
                  xlayer_occu[iz]->Fill(xptsall[iz][ix]);
#ifndef MONTECARLO
                  strp_count_set[iset]->Fill(iz, xptsall[iz][ix]); 
#endif
                }
                bool failed=false;
                if (!posinuse[0][iz][xptsall[iz][ix]]) {failed=true;}
                //		if(iz==1 && xptsall[iz][ix]>=25 && xptsall[iz][ix]<=31) {failed=true;}
                
                if (!failed){
                  xpts[iz].push_back(xptsall[iz][ix]);
                  if (filloccu) {xlayer_alloccusel->Fill(nstrip*iz+xptsall[iz][ix]);}
                }
                
              }
              xhits[iz] = xpts[iz].size();
              //	      cout <<"xhits[iz] "<< iz<<" "<<xhits[iz]<<endl;
              // if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && yptsall[iz].size()>0) { ntrigY++;} //18 01 2012
              if (filloccu) {
                ylayer_allmult[iz]->Fill(yptsall[iz].size());
#ifndef MONTECARLO
                strp_ymult_set[iz][iset]->Fill(yptsall[iz].size());
#endif			
              }
              for (int iy=0; iy<yptsall[iz].size(); iy++) {
                if (filloccu) {
                  ylayer_occu[iz]->Fill(yptsall[iz][iy]);
                  ylayer_alloccu->Fill(nstrip*iz+yptsall[iz][iy]);
#ifndef MONTECARLO
                  strp_count_set[iset]->Fill(nlayer+iz, yptsall[iz][iy]); 
#endif		  
                }
                bool failed=false;
                if (!posinuse[1][iz][yptsall[iz][iy]]) {failed=true;}
                if(iz==0 && yptsall[iz][iy]==21) {    //rejection of ystrip = 21 when ystrip 45 is having hit
                  for(int iyy=iy+1;iyy<yptsall[iz].size();iyy++){
                    if(yptsall[iz][iyy]==45) {
                      failed=true;
                      break;
                    }
                  }
                }
                //cout<<"After"<<
                if (!failed){
                  ypts[iz].push_back(yptsall[iz][iy]);
                  if (filloccu) {ylayer_alloccusel->Fill(nstrip*iz+yptsall[iz][iy]);}
                }
              }
              yhits[iz] = ypts[iz].size();
            }
            
            //Should not be only on MC code   if(ntrigX<4) continue;  
            
            if (filloccu) { // (isfill) {}
              h_xrawcorhits->Fill(-1.0, -1.0);
              
              for (int ix=0; ix<nlayer; ix++) {
                if (xptsall[ix].size()==0) continue;
                for (int iy=ix+1; iy<nlayer; iy++) {
                  if (yptsall[iy].size()==0) continue;
                  h_xrawcorhits->Fill(ix, iy);
                  if (rawxtime[ix]>-50 && rawxtime[iy]>-50) {
                    h_xtrawcorhits->Fill(ix, iy);
                  }
                }
                for (int iy=0; iy<nlayer; iy++) {
                  if (yptsall[iy].size()==0) continue;
                  h_xyrawcorhits->Fill(ix, iy);
                  if (rawxtime[ix]>-50 && rawytime[iy]>-50) {
                    h_xytrawcorhits->Fill(ix, iy);
                  }
                }
              }
              
              for (int ix=0; ix<nlayer-1; ix++) {
                if (yptsall[ix].size()==0) continue;
                for (int iy=ix+1; iy<nlayer; iy++) {
                  if (yptsall[iy].size()==0) continue;
                  h_yrawcorhits->Fill(ix, iy);
                  if (rawytime[ix]>-50 && rawytime[iy]>-50) {
                    h_ytrawcorhits->Fill(ix, iy);
                  }
                }
              }
              
              for (int ij=0; ij<nlayer; ij++) {
                xlayer_mult[ij]->Fill(xpts[ij].size());
                ylayer_mult[ij]->Fill(ypts[ij].size());
                
                rawhits_corr_xymul[ij]->Fill(xptsall[ij].size(), yptsall[ij].size());
                
                for (int ix=0; ix<xptsall[ij].size(); ix++) {
                  for (int iy=0; iy<yptsall[ij].size(); iy++) {
                    raw_occu[ij]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
                    raw_occu_set[ij][iset]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
                  }
                }
                for (int jk=ij+1; jk<nlayer; jk++) {
                  rawhits_xlay_corr_mul[ij][jk]->Fill(xptsall[ij].size(), xptsall[jk].size());
                  rawhits_ylay_corr_mul[ij][jk]->Fill(yptsall[ij].size(), yptsall[jk].size());
                }
              }
            }
            
#endif  //end of else of MONTECARLO	     
            for (int ij=0; ij<nlayer; ij++) { Xdev[ij] = 100; Xpos[ij]=0.0;}
            for (int ij=0; ij<nlayer; ij++) { Xdev1[ij] =Xdev2[ij] = 100; Xpos1[ij]= Xpos2[ij]=0.0;}
            
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
              if (isfill && Xpos[ij] >=0) xstrip_mult->Fill(6*ij+xhits[ij]);
              
              //cout<<"Xpos/Ypos"<<"		"<<Xpos[ij]<<"		"<<Ypos[ij]<<endl;
            }
            
            for (int ij=0; ij<nlayer; ij++) { Ydev[ij] = 100; Ypos[ij]=0.0;}
            for (int ij=0; ij<nlayer; ij++) { Ydev1[ij] = Ydev2[ij] = 100; Ypos1[ij]= Ypos2[ij]=0.0;}
            for (int ij=0;ij<nlayer;ij++) {
              if (yhits[ij]<=0 || yhits[ij]>nmxhits) {
                Ypos[ij]=-100;
              } else {
                for (int iy=0; iy<yhits[ij]; iy++) {
                  if (iy<yhits[ij]-1 && abs(ypts[ij][iy]-ypts[ij][iy+1])>1) { Ypos[ij]=-100; break;}
                  Ypos[ij] +=ypts[ij][iy];
                }
                if (Ypos[ij]>=0.0) {
                  Ypos[ij]  = Ypos[ij]/yhits[ij] + 0.5 - yoff[ij];
                  yyerr[ij] = yposerrsq[yhits[ij]-1][ij];//+(20.*yposerrsq[yhits[ij]-1][ij])/100.);
                }
              }
              if (isfill && Ypos[ij] >=0) ystrip_mult->Fill(6*ij+yhits[ij]);
            }
            
            //#endif  //end of else of MONTECARLO	    
            ntrigX=0; ntrigY=0;
            for (int iz=0; iz<nlayer; iz++) {
              if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && xptsall[iz].size()>0) {ntrigX++;} //18 01 2012 
              if ((iz==trigly1 || iz==trigly2 || iz==trigly3 || iz==trigly4) && yptsall[iz].size()>0) {ntrigY++;} //18 01 2012 
              
            }
            
            if(isfill) {
              // filecorOut->cd();
              EvNum = iev;
              XYhitcor = 0;
              int strp_cor[64]={0,1,21,45,1,2,3,4,3,4,45,0,1,2,32,3,4,5,0,37,0,28,0,1,20,0,1,2,58,28,56,57,62,0,30,0,20,23,0,53,54,55,58,0,1,22,0,59,0,1,2,3,28,29,56,57,58,59,10,11,48,49,50,51}; //GMA 20241126
              ULong64_t tmphitcor=0;
              for(int ij=0;ij<nlayer;ij++) {
                XYhitall[ij] = 0;
                XYhitall[ij] += xptsall[ij].size();
                XYhitall[ij]<<=8;
                XYhitall[ij] += yptsall[ij].size();
                int istx,isty,iendx,iendy;
                if(ij==0) { istx = 0;iendx=1; isty=2;iendy=3;}
                else if(ij==1) { istx = 4;iendx=6; isty=7;iendy=7;}
                else if(ij==2) { istx = 8;iendx=9; isty=10;iendy=10;}
                else if(ij==3) { istx = 11;iendx=14; isty=15;iendy=17;}
                else if(ij==4) { istx = 18;iendx=19; isty=20;iendy=21;}
                else if(ij==5) { istx = 23;iendx=24; isty=0;iendy=-1;}
                else if(ij==6) { istx = 25;iendx=28; isty=29;iendy=32;}
                else if(ij==7) { istx = 33;iendx=34; isty=35;iendy=37;}
                else if(ij==8) { istx = 38;iendx=38; isty=39;iendy=42;}
                else if(ij==9) { istx = 43;iendx=45; isty=0;iendy=-1;}
                else if(ij==10) { istx = 46;iendx=47; isty=0;iendy=-1;}
                else if(ij==11) { istx = 48;iendx=57; isty=58;iendy=63;}
                
                ULong64_t tmpxy=0;
                for(int strx=istx;strx<=iendx;strx++) {
                  ULong64_t tmpx =0;
                  for(int ii=0;ii<xptsall[ij].size();ii++) {
                    if(strp_cor[strx]==xptsall[ij][ii]){
                      tmpx = 1;
                      tmpx<<=strx;
                    }
                  }
                  tmpxy+=tmpx;
                }
                
                for(int stry=isty;stry<=iendy;stry++) {
                  ULong64_t tmpy =0;
                  for(int jj=0;jj<yptsall[ij].size();jj++) {
                    if(strp_cor[stry]==yptsall[ij][jj]){
                      tmpy  = 1;
                      tmpy<<=stry;
                    }
                  }
                  tmpxy+=tmpy;
                }
                tmphitcor += tmpxy;
                    
                if (xptsall[ij].size() >0) {
                  for (int ix=0; ix<xptsall[ij].size(); ix++) {
                    h_raw_xcorstrips[ij]->Fill(xptsall[ij][ix], -1.);
                    h_raw_xstrpnhits[ij]->Fill(xptsall[ij][ix], -1.);
                    h_raw_xystrpnhits[ij]->Fill(xptsall[ij][ix], -1.);
                    h_raw_xstrpnhits[ij]->Fill(xptsall[ij][ix],xptsall[ij].size());
                    h_raw_xystrpnhits[ij]->Fill(xptsall[ij][ix],yptsall[ij].size());
                    for(int ixx=0;ixx<xptsall[ij].size();ixx++) { 
                      
                      if(abs(xptsall[ij][ix]-xptsall[ij][ixx])>0) {
                        h_raw_xcorstrips[ij]->Fill(xptsall[ij][ix],xptsall[ij][ixx]);
                      }
                    }
                  }
                }
                
                
                if (yptsall[ij].size() >0) {
                  for (int iy=0; iy<yptsall[ij].size(); iy++) {
                    h_raw_ystrpnhits[ij]->Fill(yptsall[ij][iy], -1.);
                    h_raw_ycorstrips[ij]->Fill(yptsall[ij][iy], -1.);
                    h_raw_yxstrpnhits[ij]->Fill(yptsall[ij][iy], -1.);
                    h_raw_ystrpnhits[ij]->Fill(yptsall[ij][iy],yptsall[ij].size());
                    h_raw_yxstrpnhits[ij]->Fill(yptsall[ij][iy],xptsall[ij].size());		 
                    for(int iyy=0;iyy<yptsall[ij].size();iyy++) { //choose all of them and then use maximum value before plot
                      
                      if(abs(yptsall[ij][iy]-yptsall[ij][iyy])>0) {
                        h_raw_ycorstrips[ij]->Fill(yptsall[ij][iy],yptsall[ij][iyy]);
                      }
                    }
                  }
                }
              }
              XYhitcor=tmphitcor;
              
              T3->Fill();
            }
            
            if(ntrigX<4 && ntrigY<4) {/*cout<<iev<<"  without trigger hits"<<endl;*/ nxytrig++; }
            
            if(xptsall[trigly1].size()==0) { nx4++;}
            if(xptsall[trigly2].size()==0) { nx5++;}
            if(xptsall[trigly3].size()==0) { nx6++;}
            if(xptsall[trigly4].size()==0) { nx7++;}
            if(yptsall[trigly1].size()==0) { ny4++;}
            if(yptsall[trigly2].size()==0) { ny5++;}
            if(yptsall[trigly3].size()==0) { ny6++;}
            if(yptsall[trigly4].size()==0) { ny7++;}
            if(xptsall[trigly1].size()==0 && yptsall[trigly1].size()==0) { nxy4++;}
            if(xptsall[trigly2].size()==0 && yptsall[trigly2].size()==0) { nxy5++;}
            if(xptsall[trigly3].size()==0 && yptsall[trigly3].size()==0) { nxy6++;}
            if(xptsall[trigly4].size()==0 && yptsall[trigly4].size()==0) { nxy7++;}
            if(xptsall[trigly1].size()==0 || yptsall[trigly1].size()==0) { nxory4++;}
            if(xptsall[trigly2].size()==0 || yptsall[trigly2].size()==0) { nxory5++;}
            if(xptsall[trigly3].size()==0 || yptsall[trigly3].size()==0) { nxory6++;}
            if(xptsall[trigly4].size()==0 || yptsall[trigly4].size()==0) { nxory7++;}
            
            //Sort out hits, which can be used for fit
            for (int ij=0;ij<nlayer;ij++) {
              Xusedpos[ij] = (Xpos[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false; //Xpos[ij] : -101;
            }
            //   cout<<xhits[0]<<"  "<<yhits[0]<<"   "<<Xpos[0]<<"   "<<Ypos[0]<<endl; 
            
            int tmpnx2=0;
            for (int ix=0; ix<nlayer; ix++) { if (Xpos[ix]>-90) {tmpnx2++;} }
            h_tmp2xndf->Fill(tmpnx2);
            
            Nx=0;
            int nxfail = 0;
            xchi2 = 0;
            double xresol = 0;
            double zval[nlayer], xext[nlayer], xexter[nlayer], xposinstr[nlayer];
            for (int ix=0; ix<nlayer; ix++) { zval[ix]=layerzpos[ix];}
            for (int ix=0; ix<nlayer; ix++) { xext[ix]= xexter[ix] =xposinstr[ix] =  100;}
            
            StraightLineFit xposfit(1, zval, Xpos,  xxerr, Xusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
            xposfit.GetParameters(nxfail, xinters, xslope);
            //	    xposfit.GetError(errcst, errlin, errcov);
            xposfit.GetChisqure(Nx,xchi2);
            h_tmp3xndf->Fill(Nx);
            xposfit.GetFitValues(xext, Xdev, xexter);
            
            GetXposInStrip(xext, xoff, xposinstr);
            
            //Sort out hits, which can be used for fit
            for (int ij=0;ij<nlayer;ij++) {
              Yusedpos[ij] = (Ypos[ij]>-99 && yhits[ij]<=nmxusedhits) ? true : false; //Ypos[ij] :  -101;
            }
            
            Ny=0;
            int nyfail = 0;
            ychi2 = 0;
            double yresol = 0;
            double yext[nlayer], yexter[nlayer], yposinstr[nlayer];;
            for (int ix=0; ix<nlayer; ix++) { yext[ix]= yexter[ix] =yposinstr[ix] =  100;}
            
            StraightLineFit yposfit(1, zval, Ypos,  yyerr, Yusedpos, occulyr, occulyr, layfirst, laylast, xyPosDev);
            yposfit.GetParameters(nyfail, yinters, yslope);
            //	    yposfit.GetError(errcst, errlin, errcov);
            yposfit.GetChisqure(Ny, ychi2);
            yposfit.GetFitValues(yext, Ydev, yexter);
            GetXposInStrip(yext, yoff, yposinstr);
            
            nhits = 100*Nx + Ny;
            
            for(int ij = 0;ij<nlayer;ij++) { Xpos1[ij] = Xpos2[ij]=Xpos[ij]; Ypos1[ij] = Ypos2[ij] = Ypos[ij];}
            for (int ij=0;ij<nlayer;ij++) {
              if(ij<=5) { Xusedpos1[ij] = (Xpos1[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false; } else { Xusedpos1[ij]=false;}
            }
            
            Nx1=0;
            int nxfail1 = 0;
            xchi21 = 0;
            
            double  xext1[nlayer], xexter1[nlayer], xposinstr1[nlayer];
            
            for (int ix=0; ix<nlayer; ix++) { xext1[ix]= xexter1[ix] = 100;}
            StraightLineFit xposfit1(1, zval, Xpos1,  xxerr, Xusedpos1, occulyr,occulyr,0, 5, xyPosDev);
            xposfit1.GetParameters(nxfail1, xinters1, xslope1);
            //	    xposfit.GetError(errcst, errlin, errcov);
            xposfit1.GetChisqure(Nx1, xchi21);
            
            xposfit1.GetFitValues(xext1, Xdev1, xexter1);
            
            for (int ij=0;ij<nlayer;ij++) {
              if(ij<=5) { Yusedpos1[ij] = (Ypos1[ij]>-99 && yhits[ij]<=nmxusedhits) ? true : false; } else { Yusedpos1[ij]=false;}
            }
            
            Ny1=0;
            int nyfail1 = 0;
            ychi21 = 0;
            
            double yext1[nlayer], yexter1[nlayer], yposinstr1[nlayer];;
            for (int ix=0; ix<nlayer; ix++) { yext1[ix]= yexter1[ix] =  100;}
            
            StraightLineFit yposfit1(1, zval, Ypos1,  yyerr, Yusedpos1, occulyr, occulyr, 0, 5, xyPosDev);
            yposfit1.GetParameters(nyfail1, yinters1, yslope1);
            yposfit1.GetChisqure(Ny1, ychi21);
            yposfit1.GetFitValues(yext1, Ydev1, yexter1);
            
            for (int ij=0;ij<nlayer;ij++) {
              if(ij>5) { Xusedpos2[ij] = (Xpos2[ij]>-100 && xhits[ij]<=nmxusedhits) ? true : false; } else { Xusedpos2[ij]=false;}
            }
            
            Nx2=0;
            int nxfail2 = 0;
            xchi22 = 0;
            
            double  xext2[nlayer], xexter2[nlayer], xposinstr2[nlayer];
            
            for (int ix=0; ix<nlayer; ix++) { xext2[ix]= xexter2[ix] = 100;}
            StraightLineFit xposfit2(1, zval, Xpos2,  xxerr, Xusedpos2, occulyr, occulyr, 6, 11, xyPosDev);
            xposfit2.GetParameters(nxfail2, xinters2, xslope2);
            xposfit2.GetChisqure(Nx2, xchi22);
            
            xposfit2.GetFitValues(xext2, Xdev2, xexter2);
            
            for (int ij=0;ij<nlayer;ij++) {
              if(ij>5) { Yusedpos2[ij] = (Ypos2[ij]>-99 && yhits[ij]<=nmxusedhits) ? true : false; } else { Yusedpos2[ij]=false;}
            }
            
            Ny2=0;
            int nyfail2 = 0;
            ychi22 = 0;
            
            double yext2[nlayer], yexter2[nlayer], yposinstr2[nlayer];;
            for (int ix=0; ix<nlayer; ix++) { yext2[ix]= yexter2[ix] =  100;}
            
            StraightLineFit yposfit2(1, zval, Ypos2,  yyerr, Yusedpos2, occulyr, occulyr, 6, 11, xyPosDev);
            yposfit2.GetParameters(nyfail2, yinters2, yslope2);
            yposfit2.GetChisqure(Ny2, ychi22);
            yposfit2.GetFitValues(yext2, Ydev2, yexter2);
            
            if (Nx>= nmnhits && xchi2/(Nx-2)<mxchisq && nxfail==0) {//4Nov,2011
              if (filloccu) {
                if (iev>0 && tmptimerate >=tmpoldmuposxrate+timebin) {
                  file_out <<"posxrate "<<datafile<<" "<<iev<<" "<<nallmuposxrate<<" "<<nmuposxrate<<endl;
                  muposxrate->Fill(nmuposxrate*(1.0/timebin));
                  muposxratex->Fill(nallmuposxrate, nmuposxrate*(1.0/timebin));
                      nmuposxrate = 0;
                      nallmuposxrate++;
                      tmpoldmuposxrate = tmptimerate;
                }
                nmuposxrate++;
              }
                  
              pos_xslope[occulyr][iiter]->Fill(stripwidth*xslope);
              if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
                  if (abs(Xdev[ij])<1.5 && Xusedpos[ij]) {
                    for (int jk=0; jk<nlayer; jk++) {
                      if (abs(Xdev[jk])<1.5 && Xusedpos[ij]) {
                        deltaposxcov[ij][jk] += Xdev[ij]*Xdev[jk];
                        deltaposxCount[ij][jk] +=1 ;
                      }
                    }
                  }
#ifndef MONTECARLO
                  if (filloccu && Xusedpos[ij]) {xlayer_reso_set[ij][iset]->Fill(Xdev[ij]);}
#endif
                  if (abs(Xdev[ij])<6.0 && Xusedpos[ij]) { 
                    xlayer_reso[ij][iiterrs]->Fill(Xdev[ij]);
                    is_xyreso[ij][iiterrs] = true;
                    xlayer_exterr[ij][iiterrs]->Fill(xexter[ij]);
                    int ixx=min(xhits[ij], nmxhits);
                    if (ixx>0 && Xpos[ij]>1.0 && Xpos[ij]<nstrip-1) { 
                      xlayer_reso_mul[ij][iiterrs][ixx-1]->Fill(Xdev[ij]);}
                    
                  }
                }
                nTotallx++;
                for (int ij=0; ij<nlayer; ij++) {
                  if (abs(Xdev[ij])<1.5 && Xusedpos[ij]) {
                    for (int jk=0; jk<nlayer; jk++) {
                      if (abs(Xdev[jk])<1.5 && Xusedpos[jk]) {
                        h_posxcoreff->Fill(ij,jk) ;
                        posxcoreffCount +=1 ;
                      }
                    }
                  }
                }
              } else {
                //	xstrip_xdev[occulyr][iiterrs]->Fill(xext[occulyr],Xdev[occulyr]);
                if (abs(Xdev[occulyr])<6.0 && Xusedpos[occulyr]) {
                  is_xyreso[occulyr][iiterrs] = true;
                  xlayer_reso[occulyr][iiterrs]->Fill(Xdev[occulyr]);
                  xlayer_exterr[occulyr][iiterrs]->Fill(xexter[occulyr]);
                  
                  int ixx=min(xhits[occulyr], nmxhits);
                  if (ixx>0 && Xpos[occulyr]>1&&Xpos[occulyr]<nstrip-1 ) { 
                    xlayer_reso_mul[occulyr][iiterrs][ixx-1]->Fill(Xdev[occulyr]);}
                }
                
                if (iiter==nmxiter-1) { //Correlated hits
                  double expp = -100;
                  for (int ix=0; ix<xpts[occulyr].size(); ix++) {
                    if (abs(xext[occulyr] - xpts[occulyr][ix])<1) {
                      expp = xpts[occulyr][ix]; break;
                    }
                  }
                  if (expp >-50) {
                    h_xcorstrips[occulyr]->Fill(expp, -1.);
                    h_xycorstrips[occulyr]->Fill(expp, -1.);
                    for (int ix=0; ix<xpts[occulyr].size(); ix++) {
                      if (abs(expp - xpts[occulyr][ix])>0) { //choose all of them and then use maximum value before plot
                        h_xcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
                        //	xstrip_xdev[occulyr][iiterrs]->Fill(xepp[occulyr],Xdev[occulyr]);
                      }
                    }
                    
                    for (int ix=0; ix<ypts[occulyr].size(); ix++) {
                      h_xycorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
                    }
                  }
                } //if (iiter==nmxiter-1) 
              }
            }
            if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {//4Nov
              if (filloccu) { 
                if (iev>0 && tmptimerate >=tmpoldmuposyrate+timebin) {
                  file_out <<"posyrate "<<datafile<<" "<<iev<<" "<<nallmuposyrate<<" "<<nmuposyrate<<endl;
                  muposyrate->Fill(nmuposyrate*(1.0/timebin));
                  muposyratex->Fill(nallmuposyrate, nmuposyrate*(1.0/timebin));
                  nmuposyrate = 0;
                  nallmuposyrate++;
                  tmpoldmuposyrate = tmptimerate;
                }
                nmuposyrate++;
              }
              
              pos_yslope[occulyr][iiter]->Fill(stripwidth*yslope);
              if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
                  if (abs(Ydev[ij])<1.5 && Yusedpos[ij]) {
                    for (int jk=0; jk<nlayer; jk++) {
                      if (abs(Ydev[jk])<1.5 && Yusedpos[ij]) {
                        deltaposycov[ij][jk] += Ydev[ij]*Ydev[jk];
                        deltaposyCount[ij][jk] +=1 ;
                      }
                    }
                  }
#ifndef MONTECARLO		  
                  if (filloccu && Yusedpos[ij]) {ylayer_reso_set[ij][iset]->Fill(Ydev[ij]);}
#endif
                  if (abs(Ydev[ij])<6.0  && Yusedpos[ij]) {
                    ylayer_reso[ij][iiterrs]->Fill(Ydev[ij]);
                    if (is_xyreso[ij][iiterrs]) {xylayer_reso[ij][iiterrs]->Fill(Xdev[ij], Ydev[ij]);}
                    ylayer_exterr[ij][iiterrs]->Fill(yexter[ij]);
                    int iyy = min(yhits[ij],nmxhits); 
                    if (iyy>0 && Ypos[ij]>1.0 && Ypos[ij]<nstrip-1 ) { 
                      ylayer_reso_mul[ij][iiterrs][iyy-1]->Fill(Ydev[ij]);}
                  }
                }
                nTotally++;
                for (int ij=0; ij<nlayer; ij++) {
                  if (abs(Ydev[ij])<1.5 && Yusedpos[ij]) {
                    for (int jk=0; jk<nlayer; jk++) {
                      if (abs(Ydev[jk])<1.5 && Yusedpos[jk]) {
                        h_posycoreff ->Fill(ij,jk) ;
                        posycoreffCount +=1 ;
                      }
                    }
                  }
                }
              } else {
                if (abs(Ydev[occulyr])<6.0 && Yusedpos[occulyr]) {
                  ylayer_reso[occulyr][iiterrs]->Fill(Ydev[occulyr]);
                  if (is_xyreso[occulyr][iiterrs]) {xylayer_reso[occulyr][iiterrs]->Fill(Xdev[occulyr], Ydev[occulyr]);}
                  ylayer_exterr[occulyr][iiterrs]->Fill(yexter[occulyr]);
                  int iyy = min(yhits[occulyr],nmxhits); 
                  if (iyy>0 && Ypos[occulyr]>1.0 && Ypos[occulyr]<nstrip-1 ) { 
                    ylayer_reso_mul[occulyr][iiterrs][iyy-1]->Fill(Ydev[occulyr]);}
                }
                
                if (iiter==nmxiter-1) { //Correlated hits
                  double expp = -100;
                  for (int ix=0; ix<ypts[occulyr].size(); ix++) {
                    if (abs(yext[occulyr] - ypts[occulyr][ix])<1) {
                      expp = ypts[occulyr][ix]; break;
                    }
                  }
                  if (expp >-50) { 
                    h_ycorstrips[occulyr]->Fill(expp, -1.0);
                    h_yxcorstrips[occulyr]->Fill(expp, -1.0);
                    for (int ix=0; ix<ypts[occulyr].size(); ix++) {
                      if (abs(expp - ypts[occulyr][ix])>0) { //choose all of them and then use maximum value before plot
                        h_ycorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
                      }
                    }
                    for (int ix=0; ix<xpts[occulyr].size(); ix++) {
                      h_yxcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
                    }
                  }
                } //if (iiter==nmxiter-1) 
              }
            }
            
            if (nxfail==0 && isfill) {
              h_chisqx->Fill(xchi2);
              if (Nx-2>0) {
                h_reduchisqx->Fill(xchi2/(Nx-2));
                double probx = TMath::Prob(xchi2, Nx-2);
                h_xprob->Fill(probx);
                int ibin = getbinid(Nx, nprob, probs);
                if (ibin>=0 && ibin<nprob) {h_xnprob[ibin]->Fill(probx);}
              }
              h_xndf->Fill(Nx);
            }
            if (nyfail==0 && isfill) {
              h_chisqy->Fill(ychi2);
              if (Ny-2>0) {
                h_reduchisqy->Fill(ychi2/(Ny-2));
                double probx =TMath::Prob(ychi2, Ny-2);
                h_yprob->Fill(probx);
                int ibin = getbinid(Ny, nprob, probs);
                if (ibin>=0 && ibin<nprob) {h_ynprob[ibin]->Fill(probx);}
              }
              h_yndf->Fill(Ny);
            }
            
            double theta11 = acos(sqrt(1./(1+pow(stripwidth*xslope1,2.)+pow(stripwidth*yslope1,2.))));
            double theta1 = (180./pival)*theta11;
            double theta22 = acos(sqrt(1./(1+pow(stripwidth*xslope2,2.)+pow(stripwidth*yslope2,2.))));
            double theta2 = (180./pival)*theta22;
            double theta_diff_val1 = 0.;//fabs(theta1-theta2);
            
            int  xext1l5 = int(stripwidth*xext1[5]); int yext1l5 = int(stripwidth*yext1[5]);
            int xstp = int(xext1[5]); int ystp = int(yext1[5]);
            double zval1[12];
            for(int ij=0;ij<nlayer;ij++) { zval1[ij] = zval[ij]/stripwidth;}
            double R0 = sqrt(pow(xext1[5]-xext1[0],2.)+pow(yext1[5]-yext1[0],2.)+pow(zval1[5]-zval1[0],2.));
            double A0 = (xext1[5]-xext1[0])/R0;
            double B0 = (yext1[5]-yext1[0])/R0;
            double C0 = (zval1[5]-zval1[0])/R0;
            double R1 = sqrt(pow(xext2[11]-xext2[6],2.)+pow(yext2[11]-yext2[6],2.)+pow(zval1[11]-zval1[6],2.));
            double A1 = (xext2[11]-xext2[6])/R1;
            double B1 = (yext2[11]-yext2[6])/R1;
            double C1 = (zval1[11]-zval1[6])/R1;
            double theta_diff_val = 0.;//fabs(acos(A0*A1+B0*B1+C0*C1))*180./TMath::Pi();
            
            if(Nx1>3 && Ny1>3 && Nx2>3 && Ny2>3 && ychi2/(Ny-2)<15 && nyfail==0 && xchi2/(Nx-2)<15 && nxfail==0 && Nx>=6 && Ny>=6 ) {
              if (isfill) { 
                theta_diff_val = abs(acos(A0*A1+B0*B1+C0*C1))*180./TMath::Pi();
                theta_diff_val1 = abs(theta1-theta2);
                //                cout<<"xstp][ystp " <<xstp<<" "<<ystp<<endl;
                if(xstp>=0. && xstp <28.0 && ystp>=0. && ystp <31.0){pixel_scatang[xstp][ystp]->Fill(theta_diff_val,1.0); }
                theta_diff_1->Fill(theta_diff_val1);
                theta_diff_2->Fill(theta_diff_val);
                if(theta_diff_val1>0.0 && theta_diff_val1<=.3) pixel_diff_theta[0]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>0.3 && theta_diff_val1<=.6) pixel_diff_theta[1]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>0.6 && theta_diff_val1<=.9) pixel_diff_theta[2]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>0.9 && theta_diff_val1<=1.2) pixel_diff_theta[3]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>1.2 && theta_diff_val1<=1.5) pixel_diff_theta[4]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>1.5 && theta_diff_val1<=1.8) pixel_diff_theta[5]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>1.8 && theta_diff_val1<=2.1) pixel_diff_theta[6]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>2.1 && theta_diff_val1<=2.4) pixel_diff_theta[7]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>2.4 && theta_diff_val1<=2.7) pixel_diff_theta[8]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>2.7 && theta_diff_val1<=3.0) pixel_diff_theta[9]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>3.0 && theta_diff_val1<=3.3) pixel_diff_theta[10]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>3.3 && theta_diff_val1<=3.6) pixel_diff_theta[11]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>3.6 && theta_diff_val1<=3.9) pixel_diff_theta[12]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>3.9 && theta_diff_val1<=4.2) pixel_diff_theta[13]->Fill(xext1l5,yext1l5);
                if(theta_diff_val1>4.2 && theta_diff_val1<=4.5) pixel_diff_theta[14]->Fill(xext1l5,yext1l5);
                if(theta2<10.0){ if(theta_diff_val>=0.5 && theta_diff_val<=2.5) pixel_diff_theta[15]->Fill(xext1l5,yext1l5); }
              }
            }
            
            double fitthe1 = acos(sqrt(1./(1+pow(stripwidth*xslope,2.)+pow(stripwidth*yslope,2.)))); //10Nov//S= sqrt(dx^2+dy^2+dz^2)) theta= acos(height / S) height = dz;
            double fitthe = (180./pival)*fitthe1; //acos(sqrt(1./(1+pow(stripwidth*xslope,2.)+pow(stripwidth*yslope,2.)))); 
            
#ifdef C217STRIP
            double fitphi = atan2(-xslope, yslope); 
#else		
            double fitphi = atan2(yslope, xslope);  // What is the direction
#endif
            
            if (Ny>=nmnhits /*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) {//4Nov
              if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) {
                if (filloccu) { 
                  for (int iz=0; iz<nlayer; iz++) { 
                    xlayer_allmumult[iz]->Fill(xptsall[iz].size());
                    ylayer_allmumult[iz]->Fill(yptsall[iz].size());
                  }
                }
                
                //       C217
                //       x1 x2 x3 ...............................x30   x31
                //       ^  ^  ^  ^  ^  ^  ^  ^  W  ^  ^  ^  ^  ^  ^  ^  ^
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y1
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y2
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
                //  S<-- |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.--> N
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->.
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y30
                //       |--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|-->y31
                //                               | 
                //                               E
                
                pos_theta[occulyr][iiter]->Fill(fitthe);
                
                pos_phi[occulyr][iiter]->Fill(fitphi*180./pival);
                if (isfill) { 
#ifdef MONTECARLO
#ifdef FULLG4SIM		
                  sel_theta_mom->Fill(fitthe,abs(momin[0]));
                  sel_phi_mom->Fill(fitphi*180./pival,abs(momin[0]));
                  sel_theta_rec_gen->Fill(60.*(thegen[0]),fitthe);
#endif
#endif		  
                  costhe[0]->Fill(fitthe, 1.0);
                  phiang[0]->Fill(fitphi, 1.0);
                  
                  if (xpts[trigly1].size()>0 &&
                      xpts[trigly2].size()>0 &&
                      xpts[trigly3].size()>0 &&
                      xpts[trigly4].size()>0) {
                    costhe[1]->Fill(fitthe, 1.0);
                    phiang[1]->Fill(fitphi, 1.0);
                        
                    if (Xpos[trigly1] >-99 && Xpos[trigly2] >-99 &&
                        Xpos[trigly3] >-99 && Xpos[trigly4] >-99) {
                      costhe[2]->Fill(fitthe, 1.0);
                      phiang[2]->Fill(fitphi, 1.0);
                      
                      if (Ypos[trigly1] >-99 && Ypos[trigly2] >-99 &&
                          Ypos[trigly3] >-99 && Ypos[trigly4] >-99) {
                        costhe[3]->Fill(fitthe, 1.0);
                        phiang[3]->Fill(fitphi, 1.0);
                        
                        if (abs(Xdev[trigly1]) < 1.5 && Xusedpos[trigly1] && abs(Ydev[trigly1]) < 1.5 &&
                            abs(Xdev[trigly2]) < 1.5 && Xusedpos[trigly2] && abs(Ydev[trigly2]) < 1.5 &&
                            abs(Xdev[trigly3]) < 1.5 && Xusedpos[trigly3] && abs(Ydev[trigly3]) < 1.5 &&
                            abs(Xdev[trigly4]) < 1.5 && Xusedpos[trigly4] && abs(Ydev[trigly4]) < 1.5) {
                          costhe[4]->Fill(fitthe, 1.0);
                          phiang[4]->Fill(fitphi, 1.0);
                          costhe[5]->Fill(cos(fitthe1),1.0);//10Nov
                          zen = fitthe;
                        }
                      }
                    }
                  }
                } // if (isfill)
              } // if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0)
            } // if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0)
            
            //9th July 2016 : redundant costheta distributions for "n" check
            if (isfill && nxfail==0 && nyfail==0) { 
              bool passx[ntrkselcrit]={0};
              // if(Nx<12 && Ny<12) { 
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0) {passx[0]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0) {passx[1]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0) {passx[2]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0) {passx[3]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0) {passx[4]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0) {passx[5]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0) {passx[6]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<2.5 && ychi2/(Ny-2)<2.5) {passx[7]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<2.5 && ychi2/(Ny-2)<2.5) {passx[8]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<2.5 && ychi2/(Ny-2)<2.5) {passx[9]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<2.5 && ychi2/(Ny-2)<2.5) {passx[10]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<2.5 && ychi2/(Ny-2)<2.5) {passx[11]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<2.5 && ychi2/(Ny-2)<2.5) {passx[12]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<2.5 && ychi2/(Ny-2)<2.5) {passx[13]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<3.0 && ychi2/(Ny-2)<3.0) {passx[14]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<3.0 && ychi2/(Ny-2)<3.0) {passx[15]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<3.0 && ychi2/(Ny-2)<3.0) {passx[16]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<3.0 && ychi2/(Ny-2)<3.0) {passx[17]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<3.0 && ychi2/(Ny-2)<3.0) {passx[18]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<3.0 && ychi2/(Ny-2)<3.0) {passx[19]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<3.0 && ychi2/(Ny-2)<3.0) {passx[20]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<3.5 && ychi2/(Ny-2)<3.5) {passx[21]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<3.5 && ychi2/(Ny-2)<3.5) {passx[22]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<3.5 && ychi2/(Ny-2)<3.5) {passx[23]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<3.5 && ychi2/(Ny-2)<3.5) {passx[24]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<3.5 && ychi2/(Ny-2)<3.5) {passx[25]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<3.5 && ychi2/(Ny-2)<3.5) {passx[26]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<3.5 && ychi2/(Ny-2)<3.5) {passx[27]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<4.0 && ychi2/(Ny-2)<4.0) {passx[28]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<4.0 && ychi2/(Ny-2)<4.0) {passx[29]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<4.0 && ychi2/(Ny-2)<4.0) {passx[30]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<4.0 && ychi2/(Ny-2)<4.0) {passx[31]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<4.0 && ychi2/(Ny-2)<4.0) {passx[32]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<4.0 && ychi2/(Ny-2)<4.0) {passx[33]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<4.0 && ychi2/(Ny-2)<4.0) {passx[34]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<4.5 && ychi2/(Ny-2)<4.5) {passx[35]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<4.5 && ychi2/(Ny-2)<4.5) {passx[36]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<4.5 && ychi2/(Ny-2)<4.5) {passx[37]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<4.5 && ychi2/(Ny-2)<4.5) {passx[38]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<4.5 && ychi2/(Ny-2)<4.5) {passx[39]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<4.5 && ychi2/(Ny-2)<4.5) {passx[40]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<4.5 && ychi2/(Ny-2)<4.5) {passx[41]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<5.0 && ychi2/(Ny-2)<5.0) {passx[42]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<5.0 && ychi2/(Ny-2)<5.0) {passx[43]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<5.0 && ychi2/(Ny-2)<5.0) {passx[44]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<5.0 && ychi2/(Ny-2)<5.0) {passx[45]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<5.0 && ychi2/(Ny-2)<5.0) {passx[46]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<5.0 && ychi2/(Ny-2)<5.0) {passx[47]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<5.0 && ychi2/(Ny-2)<5.0) {passx[48]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<5.5 && ychi2/(Ny-2)<5.5) {passx[49]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<5.5 && ychi2/(Ny-2)<5.5) {passx[50]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<5.5 && ychi2/(Ny-2)<5.5) {passx[51]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<5.5 && ychi2/(Ny-2)<5.5) {passx[52]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<5.5 && ychi2/(Ny-2)<5.5) {passx[53]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<5.5 && ychi2/(Ny-2)<5.5) {passx[54]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<5.5 && ychi2/(Ny-2)<5.5) {passx[55]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<6.0 && ychi2/(Ny-2)<6.0) {passx[56]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<6.0 && ychi2/(Ny-2)<6.0) {passx[57]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<6.0 && ychi2/(Ny-2)<6.0) {passx[58]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<6.0 && ychi2/(Ny-2)<6.0) {passx[59]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<6.0 && ychi2/(Ny-2)<6.0) {passx[60]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<6.0 && ychi2/(Ny-2)<6.0) {passx[61]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<6.0 && ychi2/(Ny-2)<6.0) {passx[62]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<6.5 && ychi2/(Ny-2)<6.5) {passx[63]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<6.5 && ychi2/(Ny-2)<6.5) {passx[64]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<6.5 && ychi2/(Ny-2)<6.5) {passx[65]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<6.5 && ychi2/(Ny-2)<6.5) {passx[66]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<6.5 && ychi2/(Ny-2)<6.5) {passx[67]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<6.5 && ychi2/(Ny-2)<6.5) {passx[68]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<6.5 && ychi2/(Ny-2)<6.5) {passx[69]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<7.0 && ychi2/(Ny-2)<7.0) {passx[70]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<7.0 && ychi2/(Ny-2)<7.0) {passx[71]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<7.0 && ychi2/(Ny-2)<7.0) {passx[72]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<7.0 && ychi2/(Ny-2)<7.0) {passx[73]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<7.0 && ychi2/(Ny-2)<7.0) {passx[74]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<7.0 && ychi2/(Ny-2)<7.0) {passx[75]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<7.0 && ychi2/(Ny-2)<7.0) {passx[76]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<7.5 && ychi2/(Ny-2)<7.5) {passx[77]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<7.5 && ychi2/(Ny-2)<7.5) {passx[78]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<7.5 && ychi2/(Ny-2)<7.5) {passx[79]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<7.5 && ychi2/(Ny-2)<7.5) {passx[80]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<7.5 && ychi2/(Ny-2)<7.5) {passx[81]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<7.5 && ychi2/(Ny-2)<7.5) {passx[82]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<7.5 && ychi2/(Ny-2)<7.5) {passx[83]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<8.0 && ychi2/(Ny-2)<8.0) {passx[84]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<8.0 && ychi2/(Ny-2)<8.0) {passx[85]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<8.0 && ychi2/(Ny-2)<8.0) {passx[86]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<8.0 && ychi2/(Ny-2)<8.0) {passx[87]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<8.0 && ychi2/(Ny-2)<8.0) {passx[88]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<8.0 && ychi2/(Ny-2)<8.0) {passx[89]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<8.0 && ychi2/(Ny-2)<8.0) {passx[90]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<9.0 && ychi2/(Ny-2)<9.0) {passx[91]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<9.0 && ychi2/(Ny-2)<9.0) {passx[92]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<9.0 && ychi2/(Ny-2)<9.0) {passx[93]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<9.0 && ychi2/(Ny-2)<9.0) {passx[94]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<9.0 && ychi2/(Ny-2)<9.0) {passx[95]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<9.0 && ychi2/(Ny-2)<9.0) {passx[96]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<9.0 && ychi2/(Ny-2)<9.0) {passx[97]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<10.0 && ychi2/(Ny-2)<10.0) {passx[98]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<10.0 && ychi2/(Ny-2)<10.0) {passx[99]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<10.0 && ychi2/(Ny-2)<10.0) {passx[100]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<10.0 && ychi2/(Ny-2)<10.0) {passx[101]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<10.0 && ychi2/(Ny-2)<10.0) {passx[102]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<10.0 && ychi2/(Ny-2)<10.0) {passx[103]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<10.0 && ychi2/(Ny-2)<10.0) {passx[104]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<11.0 && ychi2/(Ny-2)<11.0) {passx[105]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<11.0 && ychi2/(Ny-2)<11.0) {passx[106]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<11.0 && ychi2/(Ny-2)<11.0) {passx[107]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<11.0 && ychi2/(Ny-2)<11.0) {passx[108]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<11.0 && ychi2/(Ny-2)<11.0) {passx[109]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<11.0 && ychi2/(Ny-2)<11.0) {passx[110]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<11.0 && ychi2/(Ny-2)<11.0) {passx[111]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<12.0 && ychi2/(Ny-2)<12.0) {passx[112]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<12.0 && ychi2/(Ny-2)<12.0) {passx[113]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<12.0 && ychi2/(Ny-2)<12.0) {passx[114]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<12.0 && ychi2/(Ny-2)<12.0) {passx[115]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<12.0 && ychi2/(Ny-2)<12.0) {passx[116]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<12.0 && ychi2/(Ny-2)<12.0) {passx[117]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<12.0 && ychi2/(Ny-2)<12.0) {passx[118]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<13.0 && ychi2/(Ny-2)<13.0) {passx[119]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<13.0 && ychi2/(Ny-2)<13.0) {passx[120]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<13.0 && ychi2/(Ny-2)<13.0) {passx[121]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<13.0 && ychi2/(Ny-2)<13.0) {passx[122]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<13.0 && ychi2/(Ny-2)<13.0) {passx[123]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<13.0 && ychi2/(Ny-2)<13.0) {passx[124]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<13.0 && ychi2/(Ny-2)<13.0) {passx[125]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<14.0 && ychi2/(Ny-2)<14.0) {passx[126]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<14.0 && ychi2/(Ny-2)<14.0) {passx[127]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<14.0 && ychi2/(Ny-2)<14.0) {passx[128]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<14.0 && ychi2/(Ny-2)<14.0) {passx[129]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<14.0 && ychi2/(Ny-2)<14.0) {passx[130]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<14.0 && ychi2/(Ny-2)<14.0) {passx[131]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<14.0 && ychi2/(Ny-2)<14.0) {passx[132]=1;}
              
              if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<15.0 && ychi2/(Ny-2)<15.0) {passx[133]=1;}
              if (Nx>=5 && Ny>=5 && xchi2/(Nx-2)<15.0 && ychi2/(Ny-2)<15.0) {passx[134]=1;}
              if (Nx>=6 && Ny>=6 && xchi2/(Nx-2)<15.0 && ychi2/(Ny-2)<15.0) {passx[135]=1;}
              if (Nx>=7 && Ny>=7 && xchi2/(Nx-2)<15.0 && ychi2/(Ny-2)<15.0) {passx[136]=1;}
              if (Nx>=8 && Ny>=8 && xchi2/(Nx-2)<15.0 && ychi2/(Ny-2)<15.0) {passx[137]=1;}
              if (Nx>=9 && Ny>=9 && xchi2/(Nx-2)<15.0 && ychi2/(Ny-2)<15.0) {passx[138]=1;}
              if (Nx>=10 && Ny>=10 && xchi2/(Nx-2)<15.0 && ychi2/(Ny-2)<15.0) {passx[139]=1;}
              
              if(fitphi<0.) { fitphi= fitphi+2*pival;}
              
              for (int jk=0; jk<ntrkselcrit; jk++) { 
                if (passx[jk]) {
                  //zenithang[0][jk]->Fill(fitthe, 1.0);
                  if ((xpts[trigly1].size()>0 &&
                       xpts[trigly2].size()>0 &&
                       xpts[trigly3].size()>0 &&
                       xpts[trigly4].size()>0) &&
                      (ypts[trigly1].size()>0 &&
                       ypts[trigly2].size()>0 &&
                       ypts[trigly3].size()>0 &&
                       ypts[trigly4].size()>0)) {
                    
                    if ((Xpos[trigly1] >-99 && Xpos[trigly2] >-99 &&
                         Xpos[trigly3] >-99 && Xpos[trigly4] >-99) && (Ypos[trigly1] >-99 && Ypos[trigly2] >-99 &&
                                                                       Ypos[trigly3] >-99 && Ypos[trigly4] >-99)) {
                      
                      zenithang[0][jk]->Fill(fitthe, 1.0);
                      azimuthang[0][jk]->Fill(fitphi*180./pival, 1.0);
                      //  if (abs(momin[0]) > 0.25) { 
                      zenithang[1][jk]->Fill(fitthe, 1.0);
                      azimuthang[1][jk]->Fill(fitphi*180./pival, 1.0);
                      //  if (abs(momin[0]) > 0.35) { 
                      zenithang[2][jk]->Fill(fitthe, 1.0);
                      azimuthang[2][jk]->Fill(fitphi*180./pival, 1.0);
                      //   if (abs(momin[0]) > 0.45) { 
                      zenithang[3][jk]->Fill(fitthe, 1.0);
                      azimuthang[3][jk]->Fill(fitphi*180./pival, 1.0);
                      //    if (abs(momin[0]) > 0.55) { 
                      zenithang[4][jk]->Fill(fitthe, 1.0);
                      azimuthang[4][jk]->Fill(fitphi*180./pival, 1.0);
                      //   if (abs(momin[0]) > 0.65) { 
                      zenithang[5][jk]->Fill(fitthe, 1.0);
                      azimuthang[5][jk]->Fill(fitphi*180./pival, 1.0);
                      //}
                      //}
                      //}
                      //}
                      //}
                      
                      //	if(fitphi<0.) { fitphi= fitphi+2*pival;}
                      
                      if(fitphi*180./pival >=0. && fitphi*180./pival <45.0) { zenithang_azimuth_8[0][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=45.0 &&  fitphi*180./pival <90.0) {zenithang_azimuth_8[1][jk]->Fill(fitthe, 1.0); }
                      else if(fitphi*180./pival >=90.0 &&  fitphi*180./pival <135.0) {zenithang_azimuth_8[2][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=135.0 &&  fitphi*180./pival <180.0){ zenithang_azimuth_8[3][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=180.0 &&  fitphi*180./pival <225.0){ zenithang_azimuth_8[4][jk]->Fill(fitthe, 1.0); }
                      else if(fitphi*180./pival >=225.0 &&  fitphi*180./pival <270.0){ zenithang_azimuth_8[5][jk]->Fill(fitthe, 1.0); }
                      else if(fitphi*180./pival >=270.0 &&  fitphi*180./pival <315.0){ zenithang_azimuth_8[6][jk]->Fill(fitthe, 1.0); }
                      else if(fitphi*180./pival >=315.0 &&  fitphi*180./pival <360.0) {zenithang_azimuth_8[7][jk]->Fill(fitthe, 1.0);}
                      
                      // if(fitphi*180./pival >=0. && fitphi*180./pival <45.0
                      if(fitphi*180./pival >=0. && fitphi*180./pival <22.5){ zenithang_azimuth[0][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=22.5 && fitphi*180./pival <45.0) {zenithang_azimuth[1][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=45.0 && fitphi*180./pival <67.5) {zenithang_azimuth[2][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=67.5 && fitphi*180./pival <90.0) {zenithang_azimuth[3][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=90.0 && fitphi*180./pival <112.5) {zenithang_azimuth[4][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=112.5 && fitphi*180./pival <135.0) {zenithang_azimuth[5][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=135.0 && fitphi*180./pival <157.5) {zenithang_azimuth[6][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=157.5 && fitphi*180./pival <180.0) {zenithang_azimuth[7][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=180.0 && fitphi*180./pival <202.5) {zenithang_azimuth[8][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=202.5 && fitphi*180./pival <225.0) {zenithang_azimuth[9][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=225.0 && fitphi*180./pival <247.5) {zenithang_azimuth[10][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=247.5 && fitphi*180./pival <270.0) {zenithang_azimuth[11][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=270.0 && fitphi*180./pival <292.5) {zenithang_azimuth[12][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=292.5 && fitphi*180./pival <315.0){ zenithang_azimuth[13][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=315.0 && fitphi*180./pival <337.5) {zenithang_azimuth[14][jk]->Fill(fitthe, 1.0);}
                      else if(fitphi*180./pival >=337.5 && fitphi*180./pival <360.0) {zenithang_azimuth[15][jk]->Fill(fitthe, 1.0);}
                      
                      if (abs(Xdev[trigly1]) < 1.5 && Xusedpos[trigly1]  &&
                          abs(Xdev[trigly2]) < 1.5 && Xusedpos[trigly2]  &&
                          abs(Xdev[trigly3]) < 1.5 && Xusedpos[trigly3]  &&
                          abs(Xdev[trigly4]) < 1.5 && Xusedpos[trigly4] ) {
                        
                        if (abs(Xdev[trigly1]) < 1.5 && Xusedpos[trigly1] && abs(Ydev[trigly1]) < 1.5 &&
                            abs(Xdev[trigly2]) < 1.5 && Xusedpos[trigly2] && abs(Ydev[trigly2]) < 1.5 &&
                            abs(Xdev[trigly3]) < 1.5 && Xusedpos[trigly3] && abs(Ydev[trigly3]) < 1.5 &&
                            abs(Xdev[trigly4]) < 1.5 && Xusedpos[trigly4] && abs(Ydev[trigly4]) < 1.5) {
                        }
                      }
                    }
                  }
                }
              } // for (int jk=0; jk<ntrkselcrit; jk++) 
            } // if (isfill && nxfail==0 && nyfail==0)
            
            if (isfill) { 
              nTotalp++;
              for (int ix=0; ix<nlayer; ix++) {
                if (Xpos[ix]>=-1 && Xpos[ix]<=nstrip) {
                  for (int iy=ix+1; iy<nlayer; iy++) {
                    if (Xpos[iy]>=-1 && Xpos[iy]<=nstrip) {
                      h_xcorhits->Fill(ix, iy);
                    }
                  }
                  for (int iy=0; iy<nlayer; iy++) {
                    if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip) {
                      h_xycorhits->Fill(ix, iy);
                    }
                  }
                }
              }
              
              for (int ix=0; ix<nlayer-1; ix++) {
                if (Ypos[ix]>=-1 && Ypos[ix]<=nstrip) {
                  for (int iy=ix+1; iy<nlayer; iy++) {
                    if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip) {
                      h_ycorhits->Fill(ix, iy);
                    }
                  }
                }
              }
            }
            
            //15th Oct
            //	    if (isfiducial) {
            if (Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0) { 
              if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
                  if ( xexter[ij]<0.5) {
                    strp_xmul[ij][iiterrs]->Fill(xposinstr[ij], xhits[ij]);
                  }
                }
              } else {
                strp_xmul[occulyr][iiterrs]->Fill(xposinstr[occulyr], xhits[occulyr]);  
              }
            }
            if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0) { 
              if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
                  if ( yexter[ij]<0.5) {
                    strp_ymul[ij][iiterrs]->Fill(yposinstr[ij], yhits[ij]);
                  }
                }
              } else {
                strp_ymul[occulyr][iiterrs]->Fill(yposinstr[occulyr], yhits[occulyr]);  
              }
            }
            
            if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0 &&
                Nx>=nmnhits/*-ntcor*/ && xchi2/(Nx-2)<mxchisq && nxfail==0 ) {//4Nov
              
              if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
                  if ( yexter[ij]<0.2 &&  xexter[ij]<0.2) {
                    //  if(ij<=3 || ij>=9) {
                    bool  isfiducial = (int(xext[ij]+xoff[ij])>1 && int(xext[ij]+xoff[ij])<58 && int(yext[ij]+yoff[ij])>1 && int(yext[ij]+yoff[ij])<61) ? 1 : 0;
                    //  } else {
                    //	 isfiducial = (int(xext[ij]+xoff[ij])>1 && int(xext[ij]+xoff[ij])<27 && int(yext[ij]+yoff[ij])>1 && int(yext[ij]+yoff[ij])<30) ? 1 : 0;
                    //  }
                    double gap=accprng+effirng; //find_max_gap(ij);
                    if (abs(xposinstr[ij])+xexter[ij]<accprng && abs(yposinstr[ij])+yexter[ij]<accprng) {
#ifdef ISEFFICIENCY
                      totalentry[ij][iiterrs]->Fill(xext[ij],yext[ij]);
#endif
                      if (isfiducial) {
                        
                        if (Xusedpos[ij]) { 
                          strp_xmulsim[ij][iiterrs]->Fill(xposinstr[ij], xhits[ij]);
                        } else {
                          strp_xmulsim[ij][iiterrs]->Fill(xposinstr[ij], 0.0);
                        }
                        if (Yusedpos[ij]) { 
                          strp_ymulsim[ij][iiterrs]->Fill(yposinstr[ij], yhits[ij]);
                        } else {
                          strp_ymulsim[ij][iiterrs]->Fill(yposinstr[ij], 0.0);
                        }
#ifdef ISEFFICIENCY			
                        inefficiency_xallpixel[ij][iiterrs]->Fill(xposinstr[ij]);
                        inefficiency_yallpixel[ij][iiterrs]->Fill(yposinstr[ij]);
#endif
                      }
#ifdef ISEFFICIENCY		      
                      if (xptsall[ij].size()>0) { triggereffi_xevt[ij][iiterrs]->Fill(xext[ij],yext[ij]);}
                      if (yptsall[ij].size()>0) { triggereffi_yevt[ij][iiterrs]->Fill(xext[ij],yext[ij]);}
                      
                      for (int ix=0; ix<xptsall[ij].size(); ix++) {
                        if (abs(xext[ij]- xptsall[ij][ix] + xoff[ij])<gap) {
                          triggereffi_x[ij][iiterrs]->Fill(xext[ij],yext[ij]);
                          break;
                        }
                      }
                      
                      for (int iy=0; iy<yptsall[ij].size(); iy++) {
                        if (abs(yext[ij]- yptsall[ij][iy] + yoff[ij])<gap) {
                          triggereffi_y[ij][iiterrs]->Fill(xext[ij],yext[ij]);
                          break;
                        }
                      }
                      if (abs(Xdev[ij])>accprng || (!Xusedpos[ij])) {
                        if (abs(Ydev[ij])>accprng || (!Yusedpos[ij])) {
                          inefficiency_corx[ij][iiterrs]->Fill(xext[ij],yext[ij]);
                        } else {
                          inefficiency_uncx[ij][iiterrs]->Fill(xext[ij],yext[ij]);
                        }
                        xposEffPass[ij][iiterrs]=true;
                      }
                      if (abs(Xdev[ij])<accprng && Xusedpos[ij]) {
                        if (isfiducial) { inefficiency_xpixel[ij][iiterrs]->Fill(xposinstr[ij]);}
                      }
                      if (abs(Ydev[ij])>accprng || (!Yusedpos[ij])) {
                        //inefficiency_uncy[ij][iiterrs]->Fill(xext[ij],yext[ij]);
                        if (abs(Xdev[ij])>accprng || (!Xusedpos[ij])) {
                          inefficiency_cory[ij][iiterrs]->Fill(xext[ij],yext[ij]);
                        } else {
                          inefficiency_uncy[ij][iiterrs]->Fill(xext[ij],yext[ij]);
                        }
                        yposEffPass[ij][iiterrs]=true;
                      }
                      if (abs(Ydev[ij])<accprng && Yusedpos[ij]) {
                        if (isfiducial) { inefficiency_ypixel[ij][iiterrs]->Fill(yposinstr[ij]);}
                      }
#endif
                    } // if (abs(xposinstr[ij])+xexter[ij]<accprng && abs(yposinstr[ij])+yexter[ij]<accprng) 
                  } // if ( yexter[ij]<0.2 &&  xexter[ij]<0.2)
                } // for (int ij=0; ij<nlayer; ij++)
              } else { // if (occulyr>=nlayer)
                if ( yexter[occulyr]<0.2 &&  xexter[occulyr]<0.2) {
                  
                  //position resolution as a function of position and multiplicity
                  if (occulyr<nlayer) { 
                    int ixx = xpts[occulyr].size()-1; if (ixx>nstr_posmx-1) ixx=nstr_posmx-1;
                    int iyy = ypts[occulyr].size()-1; if (iyy>nstr_posmx-1) iyy=nstr_posmx-1;
                    
                    if (ixx>=0 && ixx<nstr_posmx) {
                      xstr_xdev[occulyr][iiter][ixx]->Fill(xext[occulyr], Xdev[occulyr]);
                      ystr_xdev[occulyr][iiter][ixx]->Fill(yext[occulyr], Xdev[occulyr]);
                    }
                    if (iyy>=0 && iyy<nstr_posmx) {
                      xstr_ydev[occulyr][iiter][iyy]->Fill(xext[occulyr], Ydev[occulyr]);
                      ystr_ydev[occulyr][iiter][iyy]->Fill(yext[occulyr], Ydev[occulyr]);
                    }
                  }
                  
                  bool isfiducial = (int(xext[occulyr]+xoff[occulyr])>1 && int(xext[occulyr]+xoff[occulyr])<57 && 
                                     int(yext[occulyr]+yoff[occulyr])>1 && int(yext[occulyr]+yoff[occulyr])<61) ? 1 : 0;
                  double gap=accprng+effirng; // find_max_gap(occulyr);
                      
                  if (abs(xposinstr[occulyr])+xexter[occulyr]<accprng && abs(yposinstr[occulyr])+yexter[occulyr]<accprng) {
#ifdef ISEFFICIENCY
                    totalentry[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
                    if (xptsall[occulyr].size()>0) { triggereffi_xevt[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);}
                    if (yptsall[occulyr].size()>0) { triggereffi_yevt[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);}
                    
                    for (int ix=0; ix<xptsall[occulyr].size(); ix++) {
                      if (abs(xext[occulyr]- xptsall[occulyr][ix] + xoff[occulyr])<gap) {
                        triggereffi_x[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
                        break;
                      }
                    }
                    
                    for (int iy=0; iy<yptsall[occulyr].size(); iy++) {
                      if (abs(yext[occulyr]- yptsall[occulyr][iy] + yoff[occulyr])<gap) {
                        triggereffi_y[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
                        break;
                      }
                    }
#endif
                    if (isfiducial) {
                          
                      if (Xusedpos[occulyr]) { 
                        strp_xmulsim[occulyr][iiterrs]->Fill(xposinstr[occulyr], xhits[occulyr]);
                      } else {
                        strp_xmulsim[occulyr][iiterrs]->Fill(xposinstr[occulyr], 0.0);
                      }
                      if (Yusedpos[occulyr]) { 
                        strp_ymulsim[occulyr][iiterrs]->Fill(yposinstr[occulyr], yhits[occulyr]);
                      } else {
                        strp_ymulsim[occulyr][iiterrs]->Fill(yposinstr[occulyr], 0.0);
                      }
                      
#ifdef ISEFFICIENCY		      
                      inefficiency_xallpixel[occulyr][iiterrs]->Fill(xposinstr[occulyr]);
                      inefficiency_yallpixel[occulyr][iiterrs]->Fill(yposinstr[occulyr]);
#endif
                    }
#ifdef ISEFFICIENCY		    
                    if (abs(Xdev[occulyr])>gap || (!Xusedpos[occulyr])) {
                      if (abs(Ydev[occulyr])>gap || (!Yusedpos[occulyr])) {
                        
                        inefficiency_corx[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
                      } else {
                        inefficiency_uncx[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
                      }
                      xposEffPass[occulyr][iiterrs]=true;
                    }
                    if (abs(Xdev[occulyr])<gap && Xusedpos[occulyr]) {
                      if (isfiducial) { inefficiency_xpixel[occulyr][iiterrs]->Fill(xposinstr[occulyr]);}
                    }
                    if (abs(Ydev[occulyr])>gap || (!Yusedpos[occulyr])) {
                      if (abs(Xdev[occulyr])>gap || (!Xusedpos[occulyr])) {
                        inefficiency_cory[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
                      } else {
                        inefficiency_uncy[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]);
                      }
                      yposEffPass[occulyr][iiterrs]=true;
                    }
                    if (abs(Ydev[occulyr])<gap && Yusedpos[occulyr]) {
                      if (isfiducial) { inefficiency_ypixel[occulyr][iiterrs]->Fill(yposinstr[occulyr]);}
                    }
#endif
                  } //if (abs(dx)+xexter[occulyr]<accprng && abs(dy)+yexter[occulyr]<accprng)
                } //if ( yexter[occulyr]<0.2 &&  xexter[occulyr]<0.2)
              } //else of  if (occulyr>=nlayer)
            } //if (Ny>=nmnhits/*-ntcor*/ && ychi2/(Ny-2)<mxchisq && nyfail==0 ....)
            
            if ( ychi2/(Ny-2)<mxchisq && nyfail==0 && xchi2/(Nx-2)<mxchisq && nxfail==0 ) {
              if (occulyr>=nlayer) {
                for (int ij=0; ij<nlayer; ij++) {
                  if ( yexter[ij]<0.2 &&  xexter[ij]<0.2) { 
                    double gap=accprng+effirng; //find_max_gap(ij);
                    if (abs(xposinstr[ij])+xexter[ij]<accprng && abs(yposinstr[ij])+yexter[ij]<accprng) {
#ifdef ISEFFICIENCY
                      if(Nx>=5 && Ny >= 5) { totalentry_5[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                      if(Nx>=6 && Ny >= 6) { totalentry_6[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                      if(Nx>=7 && Ny >= 7) { totalentry_7[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                      if(Nx>=8 && Ny >= 8) { totalentry_8[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                      if(Nx>=9 && Ny >= 9) { totalentry_9[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                      if(Nx>=10 && Ny >= 10) { totalentry_10[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                      if(Nx>=11 && Ny >= 11) { totalentry_11[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          
                      for (int ix=0; ix<xptsall[ij].size(); ix++) {
                        if (abs(xext[ij]- xptsall[ij][ix] + xoff[ij])<gap) {
                          if(Nx>=5 && Ny >=5 ) { triggereffi_x_5[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=6 && Ny >=6 ) { triggereffi_x_6[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=7 && Ny >=7 ) { triggereffi_x_7[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=8 && Ny >=8 ) { triggereffi_x_8[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=9 && Ny >=9 ) { triggereffi_x_9[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=10 && Ny >=10 ) { triggereffi_x_10[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=11 && Ny >=11 ) { triggereffi_x_11[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          
                          break;
                        }
                      }
                      for (int iy=0; iy<yptsall[ij].size(); iy++) {
                        if (abs(yext[ij]- yptsall[ij][iy] + yoff[ij])<gap) {
                          if(Nx>=5 && Ny >=5 ) { triggereffi_y_5[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=6 && Ny >=6 ) { triggereffi_y_6[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=7 && Ny >=7 ) { triggereffi_y_7[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=8 && Ny >=8 ) { triggereffi_y_8[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=9 && Ny >=9 ) { triggereffi_y_9[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=10 && Ny >=10 ) { triggereffi_y_10[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=11 && Ny >=11 ) { triggereffi_y_11[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          break;
                        }
                      }
                          
                      if (abs(Xdev[ij])>accprng || (!Xusedpos[ij])) {
                        if (abs(Ydev[ij])>accprng || (!Yusedpos[ij])) {
                          if(Nx>=5 && Ny >=5 ) { inefficiency_corx_5[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=6 && Ny >=6 ) { inefficiency_corx_6[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=7 && Ny >=7 ) { inefficiency_corx_7[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=8 && Ny >=8 ) { inefficiency_corx_8[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=9 && Ny >=9 ) { inefficiency_corx_9[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=10 && Ny >=10 ) { inefficiency_corx_10[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=11 && Ny >=11 ) { inefficiency_corx_11[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                        } else {
                          if(Nx>=5 && Ny >=5 ) { inefficiency_uncx_5[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=6 && Ny >=6 ) { inefficiency_uncx_6[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=7 && Ny >=7 ) { inefficiency_uncx_7[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=8 && Ny >=8 ) { inefficiency_uncx_8[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=9 && Ny >=9 ) { inefficiency_uncx_9[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=10 && Ny >=10 ) { inefficiency_uncx_10[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=11 && Ny >=11 ) { inefficiency_uncx_11[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                        }
                      }
                      if (abs(Ydev[ij])>accprng || (!Yusedpos[ij])) {
                        if (abs(Xdev[ij])>accprng || (!Xusedpos[ij])) {
                          if(Nx>=5 && Ny >=5 ) { inefficiency_cory_5[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=6 && Ny >=6 ) { inefficiency_cory_6[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=7 && Ny >=7 ) { inefficiency_cory_7[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=8 && Ny >=8 ) { inefficiency_cory_8[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=9 && Ny >=9 ) { inefficiency_cory_9[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=10 && Ny >=10 ) { inefficiency_cory_10[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=11 && Ny >=11 ) { inefficiency_cory_11[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                        } else {
                          if(Nx>=5 && Ny >=5 ) { inefficiency_uncy_5[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=6 && Ny >=6 ) { inefficiency_uncy_6[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=7 && Ny >=7 ) { inefficiency_uncy_7[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=8 && Ny >=8 ) { inefficiency_uncy_8[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=9 && Ny >=9 ) { inefficiency_uncy_9[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=10 && Ny >=10 ) { inefficiency_uncy_10[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                          if(Nx>=11 && Ny >=11 ) { inefficiency_uncy_11[ij][iiterrs]->Fill(xext[ij],yext[ij]); }
                        }
                      }
#endif
                      
                    }
                  }
                }
              } else { 
                if ( yexter[occulyr]<0.2 &&  xexter[occulyr]<0.2) { 
                  double gap=accprng+effirng; //find_max_gap(ij);
                  if (abs(xposinstr[occulyr])+xexter[occulyr]<accprng && abs(yposinstr[occulyr])+yexter[occulyr]<accprng) {
#ifdef ISEFFICIENCY
                    if(Nx>=5 && Ny >= 5) { totalentry_5[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                    if(Nx>=6 && Ny >= 6) { totalentry_6[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                    if(Nx>=7 && Ny >= 7) { totalentry_7[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                    if(Nx>=8 && Ny >= 8) { totalentry_8[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                    if(Nx>=9 && Ny >= 9) { totalentry_9[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                    if(Nx>=10 && Ny >= 10) { totalentry_10[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                    if(Nx>=11 && Ny >= 11) { totalentry_11[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                    
                    for (int ix=0; ix<xptsall[occulyr].size(); ix++) {
                      if (abs(xext[occulyr]- xptsall[occulyr][ix] + xoff[occulyr])<gap) {
                        if(Nx>=5 && Ny >=5 ) { triggereffi_x_5[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=6 && Ny >=6 ) { triggereffi_x_6[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=7 && Ny >=7 ) { triggereffi_x_7[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=8 && Ny >=8 ) { triggereffi_x_8[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=9 && Ny >=9 ) { triggereffi_x_9[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=10 && Ny >=10 ) { triggereffi_x_10[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=11 && Ny >=11 ) { triggereffi_x_11[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        break;
                      }
                    }
                    for (int iy=0; iy<yptsall[occulyr].size(); iy++) {
                      if (abs(yext[occulyr]- yptsall[occulyr][iy] + yoff[occulyr])<gap) {
                        if(Nx>=5 && Ny >=5 ) { triggereffi_y_5[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=6 && Ny >=6 ) { triggereffi_y_6[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=7 && Ny >=7 ) { triggereffi_y_7[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=8 && Ny >=8 ) { triggereffi_y_8[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=9 && Ny >=9 ) { triggereffi_y_9[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=10 && Ny >=10 ) { triggereffi_y_10[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=11 && Ny >=11 ) { triggereffi_y_11[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        break;
                      }
                    }
                    if (abs(Xdev[occulyr])>gap || (!Xusedpos[occulyr])) {
                      if (abs(Ydev[occulyr])>gap || (!Yusedpos[occulyr])) {
                        if(Nx>=5 && Ny >=5 ) { inefficiency_corx_5[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=6 && Ny >=6 ) { inefficiency_corx_6[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=7 && Ny >=7 ) { inefficiency_corx_7[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=8 && Ny >=8 ) { inefficiency_corx_8[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=9 && Ny >=9 ) { inefficiency_corx_9[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=10 && Ny >=10 ) { inefficiency_corx_10[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=11 && Ny >=11 ) { inefficiency_corx_11[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                      } else {
                        if(Nx>=5 && Ny >=5 ) { inefficiency_uncx_5[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=6 && Ny >=6 ) { inefficiency_uncx_6[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=7 && Ny >=7 ) { inefficiency_uncx_7[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=8 && Ny >=8 ) { inefficiency_uncx_8[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=9 && Ny >=9 ) { inefficiency_uncx_9[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=10 && Ny >=10 ) { inefficiency_uncx_10[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=11 && Ny >=11 ) { inefficiency_uncx_11[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                      }
                    }
                    
                    if (abs(Ydev[occulyr])>gap || (!Yusedpos[occulyr])) {
                      if (abs(Xdev[occulyr])>gap || (!Xusedpos[occulyr])) {
                        if(Nx>=5 && Ny >=5 ) { inefficiency_cory_5[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=6 && Ny >=6 ) { inefficiency_cory_6[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=7 && Ny >=7 ) { inefficiency_cory_7[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=8 && Ny >=8 ) { inefficiency_cory_8[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=9 && Ny >=9 ) { inefficiency_cory_9[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=10 && Ny >=10 ) { inefficiency_cory_10[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=11 && Ny >=11 ) { inefficiency_cory_11[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                      } else {
                        if(Nx>=5 && Ny >=5 ) { inefficiency_uncy_5[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=6 && Ny >=6 ) { inefficiency_uncy_6[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=7 && Ny >=7 ) { inefficiency_uncy_7[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=8 && Ny >=8 ) { inefficiency_uncy_8[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=9 && Ny >=9 ) { inefficiency_uncy_9[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=10 && Ny >=10 ) { inefficiency_uncy_10[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                        if(Nx>=11 && Ny >=11 ) { inefficiency_uncy_11[occulyr][iiterrs]->Fill(xext[occulyr],yext[occulyr]); }
                      }
                    }
#endif
                    
                  }
                }
              }
            }
            //GMA            if ( !isTiming && isfill &&  nxfail==0 &&  nyfail==0 && Nx>nmnhits && Ny>nmnhits && xchi2/(Ny-2) >mxchisq && ychi2/(Nx-2) >mxchisq) T2->Fill();
            if ((!isTiming) && isfill &&  nxfail==0 &&  nyfail==0 && (Nx>11 || Ny>11 || gRandom->Uniform()<1.e-4) ) T2->Fill(); //07/09/2011
            //	    T2->Fill();
            // Now filling T2 contains only proper track fit data
            nxtime=0;
            nytime=0;
            ntxyla=0;
            if (isTiming) {
              //////////////////////////////////////////////
              //                                          //
              //  Timing informations and directionality  //
              //                                          //
              //////////////////////////////////////////////
              
              //Store time informations
#ifndef MONTECARLO
              
              for(int jk=0;jk<nlayer;jk++) {
                for (int itdc=0; itdc<nTDCpLayer; itdc++) { 
                  
                  int jkrd = (itdc==0) ? jk : jk + 32;
                  int jkfl = (itdc==0) ? jk : jk + 12;
                  //x-side tdc data
                  //		cout <<"time "<<iev<<" "<<jk<<" "<<event->tdcdata[jk] <<" "<<event->tdcdata[jk+16]<<endl;
                  int tmpxxtime = event->tdcdata[jkrd];
                  
                  initxtime[jk][itdc]=(tmpxxtime<50.) ? -100.0 : 0.1*tmpxxtime;
                  if(filloccu) {
                    if (jk==0 && itdc==0) {time_layerstrip->Fill(-1., -1.);} //For total count
                    //printf("inside tdcdate read\n");
                    time_layer->Fill(jkfl, tmpxxtime);
                    
                    if(itdc==0 && jk==3) rawtimex->Fill(tmpxxtime);
                    
                    if (xallhits[jk][itdc]>=0 && xallhits[jk][itdc]<nstrip) {
                      time_xraw[jk][xallhits[jk][itdc]]->Fill(tmpxxtime);
                      time_layerstrip->Fill(nstrip*jk+xallhits[jk][itdc],tmpxxtime);
#ifdef TIMESLOPE
                      if (initxtime[jk][itdc]>0 && abs(xslope)<1.0&&abs(yslope)<1.0) { 
                        time_xslope[jk][itdc]->Fill(yext[jk], initxtime[jk][itdc]); 
                      }
#endif
                    }
                  }
                  
                  jkrd +=16;
                  jkfl +=25; // one gap after nTDCpLayer*nlayer
                  
                  //y-side tdc data
                  int tmpyytime = event->tdcdata[jkrd];
                  
                  initytime[jk][itdc]=(tmpyytime<50.0) ? -100.0 : 0.1*tmpyytime;
                  
                  if(filloccu) {
                    time_layer->Fill(jkfl, tmpyytime);
                    if (yallhits[jk][itdc]>=0 && yallhits[jk][itdc]<nstrip) {
                      time_yraw[jk][yallhits[jk][itdc]]->Fill(tmpyytime);
                      time_layerstrip->Fill(nstrip*(jk+nlayer)+yallhits[jk][itdc], tmpyytime);
#ifdef TIMESLOPE
                      if (initytime[jk][itdc]>0 && abs(xslope)<1.0&&abs(yslope)<1.0) { 
                        time_yslope[jk][itdc]->Fill(xext[jk], initytime[jk][itdc]); 
                      }
#endif
                    }
                  }
                } // for(int jk=0;jk<nlayer;jk++)
              } // for (int itdc=0; itdc<nTDCpLayer; itdc++)
              
#endif // ifndef MONTECARLO	  
              for (int ij=0; ij<nlayer; ij++) { 
                dist[ij]= -100.; istrxtime[ij] = istrytime[ij] = -1;
                passxtime[ij] = passytime[ij] = false; 
                passxtmy[ij] = passytmx[ij] = true;
#ifdef MONTECARLO
                rawxtime[ij] = rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] = timesx[ij] = xtime[ij];
                rawytime[ij] = rawytime1[ij] = rawytime2[ij] = rawytime3[ij] = timesy[ij] = ytime[ij];
                xusedtime[ij] = yusedtime[ij] = true;
#else		
                rawxtime[ij] = rawxtime0[ij] = rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] =-100;
                rawytime[ij] = rawytime0[ij] = rawytime1[ij] = rawytime2[ij] = rawytime3[ij] =-100;
                xusedtime[ij] = yusedtime[ij] = false;
#endif
                xtdev[ij] = 100;
                ytdev[ij] = 100;
              }
              
              xval=-100; yval=-100;
              int init=-1;
              double initzpos=0;
              for( int ij=0;ij<nlayer;ij++) {
                if (init<0) {
                  xval = xext[ij];
                  yval = yext[ij];
                  dist[ij] = 0.0;
                  init = ij;
                  initzpos = layerzpos[ij];
                } else {
                  dist[ij] = sqrt( pow((xext[ij] - xval)*stripwidth, 2.) + 
                                   pow((yext[ij] - yval)*stripwidth, 2.) +
                                   pow(layerzpos[ij] - initzpos, 2.));
                }
                
                if(abs(Xdev[ij]) < 2.0 && abs(Ydev[ij]) < 2.0 && 
                   xpts[ij].size()>0 && xpts[ij].size() <=nmxhits && ypts[ij].size()>0 && ypts[ij].size() <=nmxhits) {
                  
                  
                  //calcualte strips where signal come earlier and the value of offset
                  double tshft=1000.0;
                  bool tpass=true;
                  passxtime[ij] = true;
                  passxtmy[ij] = true; // (yext[ij] > firstYstrip && yext[ij] < lastYstrip) ? true : false;
                  
                  int istr = istrxtime[ij] = xpts[ij][0];
#ifndef MONTECARLO
                  for (int ix=0; ix<xpts[ij].size(); ix++) {
                    if ((!timeinuse[0][ij][xpts[ij][ix]]) || 
                        xpts[ij][ix]< firstXstrip || 
                        xpts[ij][ix]> lastXstrip) tpass= passxtime[ij] = false;
                    
                    if (xtoffset[ij][xpts[ij][ix]]<tshft) { 
                      tshft =xtoffset[ij][xpts[ij][ix]]; istr = istrxtime[ij] = xpts[ij][ix];
                    }
                  }
                  
                  timesx[ij] = initxtime[ij][0];
                  
                  rawxtime[ij] = timesx[ij];
                  for(int ijk=0;ijk<nlayer-1;ijk++) {
                    double xtime_diff = rawxtime[ijk]-rawxtime[ijk+1];
                    if(rawxtime[ijk]>0. && rawxtime[ijk+1]>0.) xlay_timediff[ijk]->Fill(xtime_diff);
                  }
                  for(int ijk=0;ijk<nlayer-6;ijk++) {
                    double xtime_diff_l1 = rawxtime[1]-rawxtime[ijk+4];
                    if(rawxtime[1]>0. && rawxtime[ijk+4]>0.)  xlay_timediff_l1[ijk]->Fill(xtime_diff_l1);
                  }
                  for(int ijk=0;ijk<nlayer-7;ijk++) {
                    double xtime_diff_l4 = rawxtime[4]-rawxtime[ijk+5];
                    if(rawxtime[4]>0. && rawxtime[ijk+5]>0.)  xlay_timediff_l4[ijk]->Fill(xtime_diff_l4);
                  }
                  
                  timesx[ij] -=timeoffsetx[ij];
                  
                  rawxtime0[ij] = timesx[ij];
                  
                  timesx[ij] -=slope_path*yext[ij]; //(5./32.)*yext[ij];
                  
                  rawxtime1[ij] = timesx[ij];
                  timesx[ij] -=tshft;
                  rawxtime2[ij] = timesx[ij];
                  
                  istr = int(yext[ij]+0.5);
                  if (istr<0) istr=0;
                  if (istr>=nstrip) istr = nstrip-1;
                  double dx = yext[ij]-istr;
                  
                  // Linear extrapolation using only two points
                  if ((istr==0 && dx<=0.0) || (istr==nstrip-1 && dx>=0.0)) { 
                    timesx[ij] -=xt_slope_cor[ij][istrxtime[ij]][istr];
                  } else if (dx>0) {
                    timesx[ij] -=(1-dx)*xt_slope_cor[ij][istrxtime[ij]][istr]+dx*xt_slope_cor[ij][istrxtime[ij]][istr+1]; 
                  } else {
                    timesx[ij] -=abs(dx)*xt_slope_cor[ij][istrxtime[ij]][istr-1]+(1-abs(dx))*xt_slope_cor[ij][istrxtime[ij]][istr]; 
                  }
                  //  rawxtime[ij] =  rawxtime0[ij] =  rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] = timesx[ij];
                  if ((xpts[ij].size() >0 && xpts[ij].size()<=nmxtimehit) && (usealltime || (ij==occulyr) || (tpass && passxtmy[ij]))) {
                    int isiz = max(0, min(nmxtimehit,int(xpts[ij].size()))-1);
                    xtime[ij] = timesx[ij] - parabolic(xposinstr[ij], strpos_vs_time[isiz][ij]);
                    rawxtime3[ij] = xtime[ij] + slope_path*yext[ij]; //undo the pathlenght correction
                    xusedtime[ij] = (xpts[ij].size()<=nmxusedtimehit) ? true : false; //  xtime[ij] : -101;
                  } // else {rawxtime[ij] = rawxtime0[ij] = rawxtime1[ij] = rawxtime2[ij] = rawxtime3[ij] =-90;}
#endif
                  
                  tshft=1000.0;
                  tpass=true;
                  passytime[ij] = true;
                  passytmx[ij] = true; //(xext[ij] > firstXstrip && xext[ij] < lastXstrip) ? true : false;
                  istr = istrytime[ij] = ypts[ij][0];
#ifndef MONTECARLO
                  for (int iy=0; iy<ypts[ij].size(); iy++) {
                    if ((!timeinuse[1][ij][ypts[ij][iy]]) || 
                        ypts[ij][iy]< firstYstrip || 
                        ypts[ij][iy]> lastYstrip) tpass = passytime[ij] = false;
                    
                    if (ytoffset[ij][ypts[ij][iy]]<tshft) { 
                      tshft =ytoffset[ij][ypts[ij][iy]]; istr = istrytime[ij] = ypts[ij][iy];
                    }
                  }
                  
                  timesy[ij] = initytime[ij][0];
                  
                  rawytime[ij] = timesy[ij];
                  for(int ijk=0;ijk<nlayer-1;ijk++) {
                    double ytime_diff = rawytime[ijk]-rawytime[ijk+1];
                    if(rawytime[ijk]>0. && rawytime[ijk+1]>0.) ylay_timediff[ijk]->Fill(ytime_diff);
                  }
                  
                  for(int ijk=0;ijk<nlayer-6;ijk++) {
                    double ytime_diff_l1 = rawytime[1]-rawytime[ijk+4];
                    if(rawytime[1]>0.&& rawytime[ijk+4]>0.)  ylay_timediff_l1[ijk]->Fill(ytime_diff_l1);
                  }
                  for(int ijk=0;ijk<nlayer-7;ijk++) {
                    double ytime_diff_l4 = rawytime[4]-rawytime[ijk+5];
                    if(rawytime[4]>0.&& rawytime[ijk+5]>0.) ylay_timediff_l4[ijk]->Fill(ytime_diff_l4);
                  }
                  
                  timesy[ij] -=timeoffsety[ij];
                  
                  rawytime0[ij] = timesy[ij];
                  timesy[ij] -=ytimeshift; // shift y-time to match with x-time
#ifdef C217STRIP
                  timesy[ij] -= 5. - slope_path*xext[ij];
#else	
                  timesy[ij] -= slope_path*xext[ij]; // Here both X/Y will follow same convention 5. - slope_path*xext[ij]; //(5./32.)*xext[ij];
#endif
                  rawytime1[ij] = timesy[ij];
                  
                  
                  timesy[ij] -=tshft;
                  rawytime2[ij] = timesy[ij];
                  
                  istr = int(xext[ij]+0.5);
                  if (istr<0) istr=0;
                  if (istr>=nstrip) istr = nstrip-1;
                  
                  double dy = xext[ij]-istr;
                  
                  // Linear extrapolation using only two points
                  if ((istr==0 && dy<=0.0) || (istr==nstrip-1 && dy>=0.0)) { 
                    timesy[ij] -=yt_slope_cor[ij][istrytime[ij]][istr];
                  } else if (dy>0) {
                    timesy[ij] -=(1-dy)*yt_slope_cor[ij][istrytime[ij]][istr]+dy*yt_slope_cor[ij][istrytime[ij]][istr+1]; 
                  } else {
                    timesy[ij] -=abs(dy)*yt_slope_cor[ij][istrytime[ij]][istr-1]+(1-abs(dy))*yt_slope_cor[ij][istrytime[ij]][istr]; 
                  }
                  
                  if ((ypts[ij].size() >0 && ypts[ij].size()<=nmxtimehit) && (usealltime || (ij==occulyr) || (tpass && passytmx[ij]))) { 
                    int isiz = max(0, min(nmxtimehit,int(ypts[ij].size()))-1)+nmxtimehit;
                    ytime[ij] = timesy[ij] - parabolic(yposinstr[ij], strpos_vs_time[isiz][ij]);
#ifdef C217STRIP 
                    rawytime3[ij] = ytime[ij] + 5. - slope_path*xext[ij];
#else
                    rawytime3[ij] = ytime[ij] + slope_path*xext[ij];
#endif
                    yusedtime[ij] = (ypts[ij].size()<=nmxusedtimehit) ?  true : false; //ytime[ij] : -101;
                    if (ypts[ij].size()>nmxusedtimehit && yusedtime[ij]) cout <<"ypts "<<ij<<" "<< ypts[ij].size()<<" "<<nmxusedtimehit<<" "<<int(yusedtime[ij])<<endl;
                  } // else { rawytime[ij] = rawytime0[ij] = rawytime1[ij] = rawytime2[ij] = rawytime3[ij] =-90;}
#endif
                  int ifirst=0;//4;//1;//0; //was used one, while Layer-0 was not active
                  double dtime = sqrt (pow((xext[ij] - xext[ifirst])*stripwidth, 2.) + 
                                       pow((yext[ij] - yext[ifirst])*stripwidth, 2.) +
                                       pow(layerzpos[ij] - layerzpos[ifirst], 2.0))/cval;
                  
                  double dtime2 = sqrt (pow((xext[ij] - xext[nlayer-2])*stripwidth, 2.) +
                                        pow((yext[ij] - yext[nlayer-2])*stripwidth, 2.) + 
                                        pow(layerzpos[ij] - layerzpos[nlayer-2], 2.0))/cval;
                  if (firstiter==0 && ntcor==0) {    
                    if (init==ifirst && xpts[ifirst].size()>0 && xpts[ifirst].size()<=nmxhits && xtime[ifirst]>-50 && xpts[ij].size()>0 && xpts[ij].size()<=nmxhits && xtime[ij]>-50) {
                      
                      if ((usealltime || timeinuse[0][ij][istrxtime[ij]]) &&  xusedtime[ij] && xusedtime[ifirst]) {
                        timex_shift[ij][iiter+1]->Fill(xtime[ij] - xtime[ifirst] + dtime);
                        timex_2dshift[ij][iiter+1]->Fill(xext[ij], xtime[ij] - xtime[ifirst] + dtime);
                        if (iiter==0) { 
                          timex_shift[ij][0]->Fill(rawxtime[ij] - rawxtime[ifirst] + dtime);
                          timex_2dshift[ij][0]->Fill(xext[ij], rawxtime[ij] - rawxtime[ifirst] + dtime);
                        }
                      }
                    }
                    
                    if (init==ifirst && ypts[ifirst].size()>0 && ypts[ifirst].size()<=nmxhits && ytime[ifirst]>-50 && ypts[ij].size()>0 && ypts[ij].size()<=nmxhits && ytime[ij]>-50) {
                      if ((usealltime || timeinuse[1][ij][istrytime[ij]]) &&  yusedtime[ij] && yusedtime[ifirst]) {
                        timey_shift[ij][iiter+1]->Fill(ytime[ij] - ytime[ifirst] + dtime);
                        timey_2dshift[ij][iiter+1]->Fill(yext[ij], ytime[ij] - ytime[ifirst] + dtime);
                        if (iiter==0) { 
                          timey_shift[ij][0]->Fill(rawytime[ij] - rawytime[ifirst] + dtime);
                          timey_2dshift[ij][0]->Fill(yext[ij], rawytime[ij] - rawytime[ifirst] + dtime);
                        }
                      }
                    }
                  }
                  
                  if (isfill) {
                    if (xpts[ij].size()==1 && timesx[ij]>0 && passxtime[ij]) {
                      if (xpts[0].size() ==1 && xpts[0][0] >=13 && xpts[0][0] <=15 &&
                          ypts[0].size() ==1 && ypts[0][0] >=13 && ypts[0][0] <=15 ) {
                        
                        timex_fy[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[0] + dtime);
                        timex_pfy[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[0] + dtime);
                        timex_pfy[ij]->Fill(nstrip*nstrip + ypts[ij][0], timesx[ij] - timesx[0] + dtime);
                        
                      }
                      
                      if (xpts[nlayer-1].size() ==1 && xpts[nlayer-1][0] >=13 && xpts[nlayer-1][0] <=15 &&
                          ypts[nlayer-1].size() ==1 && ypts[nlayer-1][0] >=13 && ypts[nlayer-1][0] <=15 ) {
                        
                        timex_fy2[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[nlayer-1] - dtime2);
                        timex_pfy2[ij]->Fill(xpts[ij][0]*nstrip + ypts[ij][0], timesx[ij] - timesx[nlayer-1] - dtime2);
                        timex_pfy2[ij]->Fill(nstrip*nstrip + ypts[ij][0], timesx[ij] - timesx[nlayer-1] - dtime2);
                      }
                    } // if (xpts[ij].size()==1 && timesx[ij]>0 && passxtime[ij])
                    
                    if (ypts[ij].size()==1 && timesy[ij]>0.0 && passytime[ij] ) {
                      if (xpts[0].size() ==1 && xpts[0][0] >=13 && xpts[0][0] <=15 &&
                          ypts[0].size() ==1 && ypts[0][0] >=13 && ypts[0][0] <=15 ) {
                        
                        timey_fx[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[0] + dtime);
                        timey_pfx[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[0] + dtime);
                        timey_pfx[ij]->Fill(nstrip*nstrip + xpts[ij][0], timesy[ij] - timesy[0] + dtime);
                      }
                      
                      if (xpts[nlayer-1].size() ==1 && xpts[nlayer-1][0] >=13 && xpts[nlayer-1][0] <=15 &&
                          ypts[nlayer-1].size() ==1 && ypts[nlayer-1][0] >=13 && ypts[nlayer-1][0] <=15 ) {
                        
                        timey_fx2[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[nlayer-1] - dtime2);
                        timey_pfx2[ij]->Fill(ypts[ij][0]*nstrip + xpts[ij][0], timesy[ij] - timesy[nlayer-1] - dtime2);
                        timey_pfx2[ij]->Fill(nstrip*nstrip + xpts[ij][0], timesy[ij] - timesy[nlayer-1] - dtime2);
                      }
                    } // if (ypts[ij].size()==1 && timesy[ij]>0.0 && passytime[ij] )
                  } // if (isfill)
                } else { //if(abs(Xdev[ij]) < xyPosDev && abs(Ydev[ij]) < xyPosDev)
                  xusedtime[ij] = yusedtime[ij] = false;
                }
              } //for(int ij=0;ij<nlayer;ij++)
              
              int tmpxtent[nlayer];
              int tmpytent[nlayer];
              
              timexslope = 0;
              xc0inters = 0;
              xt0chi2 = 1.e+6;
              nxtfail = 0;
              
              timeyslope = 0;
              yc0inters = 0;
              yt0chi2 = 1.e+6;
              nytfail = 0;
              
              for (int ix=0; ix<nlayer; ix++) {xtext[ix]= xtexter[ix] = 1000;}
              
              // Xtime fit
              //	      int iTimeSlopeConst = 0; //((isTimeCorOrReso) && ntcor==1) ? 0 : -1; //-1 : 0;
              int iTimeSlopeConst = (iiter<nmxiter-1 && ntcor==1) ? 0 : -1;
              StraightLineFit xtimefit(iTimeSlopeConst, dist, xtime,  timeserrx2, xusedtime, occulyr, occulyr, xtcorstr, xtcorend, float(7.0));
              xtimefit.GetParameters(nxtfail, xc0inters, timexslope);
              //	    xtimefit.GetError(errcst, errlin, errcov);
              xtimefit.GetChisqure(nxtime, xt0chi2);
              xtimefit.GetFitValues(xtext, xtdev, xtexter);
              
              if (nxtfail==0 && isfill) {
                h_tchisqx->Fill(xt0chi2);
                if (nxtime-2>0) {
                  h_treduchisqx->Fill(xt0chi2/(nxtime-2));
                  double probx =TMath::Prob(xt0chi2, nxtime-2);
                  h_xtprob->Fill(probx);
                  int ibin = getbinid(nxtime, nprob, probs);
                  if (ibin>=0 && ibin<nprob) {h_xtnprob[ibin]->Fill(probx);}
                }
                h_txndf->Fill(nxtime);
              }
              
              if ((isfill || occulyr==nlayer) && nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq) {
                for(int ij=0;ij<nlayer;ij++){
                  
                  if (dist[ij]<0 || xtime[ij] <-50) continue;
                  double dt1 = xtdev[ij]-biasinxtime[ij];
                  if (abs(dt1)>3*sigmarpc || (!xusedtime[ij])) continue;
                  
                  if (isfill) { 
                    for(int jk=0;jk<nlayer;jk++){
                      
                      if (dist[jk]<0 || xtime[jk] <-50) continue;
                      //                                      if(jk != ij )continue;
                      double dt2 = xtdev[jk]-biasinxtime[jk];
                      if (abs(dt2) >3*sigmarpc || (!xusedtime[jk])) continue;
                      deltatcov[ij][jk] += dt1*dt2;
                      deltatCount[ij][jk] +=1 ;
                    }
                  }
                  //140923
                  if (occulyr==nlayer) {
                    if (passxtime[ij]) { 
                      ntotxtime++; totxtime +=xtime[ij];
                      for(int jk=ij+1;jk<nlayer;jk++){
                        if (dist[jk]<0 || xtime[jk] <-50 || abs(xtdev[jk]-biasinxtime[jk])>4*sigmarpc || (!xusedtime[jk]) || (!passxtime[jk])) continue;
                        timex_correl[ij][jk][iiter+1]->Fill(xtime[jk]-xtime[ij] - (dist[ij]-dist[jk])/cval);
                        if (iiter==0) { timex_correl[ij][jk][0]->Fill(rawxtime[jk]-rawxtime[ij] - (dist[ij]-dist[jk])/cval);}
                      } // for(jk=ij+1;jk<nlayer;jk++)
                    }
                  }
                }// for(int ij=0;ij<nlayer;ij++)
              } //if ((isfill || occulyr==nlayer) && nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq) 
              if (nxtime>=nmnhits/*-ntcor*/ && nxtfail==0) {//4Nov For resolution 8layers are considered
                
                //Was it twice 150101
                //Or has memory overflow !!!!!!!!!
                if (xt0chi2/(nxtime-2)<mxtimechisq) { 
                  if (filloccu) {
                    if (iev>0 && tmptimerate >=tmpoldmutimexrate+timebin) {
                      file_out <<"timexrate "<<datafile<<" "<<iev<<" "<<nallmutimexrate<<" "<<nmutimexrate<<endl;
                      mutimexrate->Fill(nmutimexrate*(1.0/timebin));
                      mutimexratex->Fill(nallmutimexrate, nmutimexrate*(1.0/timebin));
                      nmutimexrate = 0;
                      nallmutimexrate++;
                      tmpoldmutimexrate = tmptimerate;
                    }
                    nmutimexrate++;
                  }	
                  
                  dir_cx[occulyr][iiter]->Fill(-cval*timexslope);
                  int ibin = getbinid(nxtime, nprob, probs);
                  if (ibin>=0 && ibin<nprob) {dir_cxlay[occulyr][iiter][ibin]->Fill(-cval*timexslope);if(ibin>4 && (-cval*timexslope >2.0)){ cout<<"dir_bias"<<"	"<<iev<<"	"<<ibin<<"	"<<-cval*timexslope<<endl;}}
                  
                  dir_cx2[occulyr][iiter]->Fill(xtimefit.GetSlope2());
                  dir_c0x[occulyr][iiter]->Fill(xc0inters);
                }
                if (isfill) dir_cxchi->Fill(xt0chi2, -1./timexslope/cval);
                
                if (occulyr >=nlayer) {
                  for (int ij=0; ij<nlayer; ij++) {
                    if (dist[ij]<0 || xtime[ij] <-50) continue;
                    if (xusedtime[ij]) { 
                      if (istrxtime[ij]>=0 && istrxtime[ij] <nstrip && passxtmy[ij] ) { 
                        time_xstrreso[ij][istrxtime[ij]][iiterrs]->Fill(xtdev[ij]); 
                      }
                      xtime_exterr[ij][iiterrs]->Fill(xtexter[ij]);   
                    }
                    if (passxtime[ij] && passxtmy[ij]) {
                      if (xusedtime[ij]) {
#ifndef MONTECARLO
                        if (filloccu) { time_xreso_set[ij][iset]->Fill(xtdev[ij]);}
#endif
                        time_xreso[ij][iiterrs]->Fill(xtdev[ij]);
                        istime_xyreso[ij][iiterrs] = true;
                      }  // Filling and saving. Not fitted.
                      int ixx = max(0,min(xhits[ij],nmxtimehit)-1); 
                      time_mulxreso[ij][iiterrs][ixx]->Fill(xtdev[ij]);
                    }
#ifdef ISEFFICIENCY
                    if (xposEffPass[ij][iiterrs]) { 
                      total_xt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]); //11/12/2011
                      if((!xusedtime[ij]) || abs(xtdev[ij])>5*time_xrms[ij]) { 
                        inefficiency_xt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]);
                      }
                    }
#endif
                  }
                } else {
                  if (dist[occulyr] >=0 && xtime[occulyr] >-50.0) {
                    if (passxtime[occulyr]) { 
                      if (abs(xtdev[occulyr]) >8.0) {
                        if (abs(xtdev[occulyr]) <30.0) {
                          if (xtdev[occulyr]>0) {
                            nxypos_xytdev[0]->Fill(xext[occulyr], yext[occulyr]);
                          } else {
                            nxypos_xytdev[1]->Fill(xext[occulyr], yext[occulyr]);
                          }
                        }
                      }
                      
                      if (xusedtime[occulyr]) {
                        xstr_xtdev[occulyr][iiter]->Fill(xext[occulyr], xtdev[occulyr]);
                        widthx_timex[occulyr][iiter]->Fill(xtdev[occulyr], widthx[occulyr]); //160315
                        
                        ystr_xtdev[occulyr][iiter]->Fill(yext[occulyr], xtdev[occulyr]);
                        
                        if (abs(xtdev[occulyr])<20) {
                          nxystr_xtdev[occulyr][iiter]->Fill(xext[occulyr], yext[occulyr], 1.);
                          xystr_xtdev[occulyr][iiter]->Fill(xext[occulyr], yext[occulyr], xtdev[occulyr]+20.0);
                        }
                      }
                      if (iiter==nmxiter-1) {
                        double diff = -30;
                        for (int jkl=0; jkl<6; jkl++) {
                          switch(jkl) {
                          case 0 : diff = rawxtime0[occulyr]-xtext[occulyr]; break;
                          case 1 : diff = rawxtime1[occulyr]-xtext[occulyr]; break;
                          case 2 : diff = rawxtime2[occulyr]-xtext[occulyr]; break;
                          case 3 : diff = timesx[occulyr]-xtext[occulyr]; break;
                          case 4 : diff = rawxtime3[occulyr]-xtext[occulyr]; break;
                          default : diff = xtime[occulyr]-xtext[occulyr]; break;
                          }
                          if (abs(diff)<20 && xusedtime[occulyr]) {
                            nxystr_xtdev[occulyr][nmxiter+jkl]->Fill(xext[occulyr], yext[occulyr], 1.);
                            xystr_xtdev[occulyr][nmxiter+jkl]->Fill(xext[occulyr], yext[occulyr], diff+20.0);
                          }
                        }
                        
                        int nsize=xpts[occulyr].size();
                        if (nsize>=1 && nsize<=nmxtimehit) {
                          xpos_xdev[occulyr][nsize-1]->Fill(xposinstr[occulyr], Xdev[occulyr]);
                          ypos_xdev[occulyr][nsize-1]->Fill(yposinstr[occulyr], Xdev[occulyr]);
                          xpos_xtdev[occulyr][nsize-1]->Fill(xposinstr[occulyr], xtdev[occulyr]);
                          ypos_xtdev[occulyr][nsize-1]->Fill(yposinstr[occulyr], xtdev[occulyr]);
                          
                          xpos_xtdev_str[occulyr][nsize-1]->Fill(xposinstr[occulyr], diff); 
                          xpos_xtdev_glb[occulyr][nsize-1]->Fill(xext[occulyr], diff);
                          ypos_xtdev_str[occulyr][nsize-1]->Fill(yposinstr[occulyr], diff); 
                          ypos_xtdev_glb[occulyr][nsize-1]->Fill(yext[occulyr], diff);
                          
                        }
                      }
                    }
                    if (xusedtime[occulyr]) {
                      if (istrxtime[occulyr]>=0 && istrxtime[occulyr] <nstrip) { 
#ifdef TIMESLOPE
                        time_xslope_pr[occulyr][istrxtime[occulyr]][iiter]->Fill( yext[occulyr], xtdev[occulyr]);
#endif
                        if (passxtmy[occulyr]) {time_xstrreso[occulyr][istrxtime[occulyr]][iiterrs]->Fill(xtdev[occulyr]);}
                      }
                      xtime_exterr[occulyr][iiterrs]->Fill(xtexter[occulyr]);
                    }
                    if (passxtime[occulyr] && passxtmy[occulyr]) {
                      if (xusedtime[occulyr]) {
                        time_xreso[occulyr][iiterrs]->Fill(xtdev[occulyr]);
                        istime_xyreso[occulyr][iiterrs] = true;
                      }
                      int ixx = max(0,min(xhits[occulyr],nmxtimehit)-1); 
                      time_mulxreso[occulyr][iiterrs][ixx]->Fill(xtdev[occulyr]);
                    }
#ifdef ISEFFICIENCY
                    if (xposEffPass[occulyr][iiterrs]) { 
                      total_xt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]); //11/12/2011
                      if((!xusedtime[occulyr]) || abs(xtdev[occulyr])>5*time_xrms[occulyr]) { 
                        inefficiency_xt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]); //150112
                      }
                    }
#endif
                  } //if (dist[occulyr] >=0 && xtime[occulyr] >-50.0)
                  
                  if (iiter==nmxiter-1) { //Correlated hits
                    if (abs(timesx[occulyr] - xtext[occulyr])<10.0) {
                      double expp = -100;
                      for (int ix=0; ix<xpts[occulyr].size(); ix++) {
                        if (abs(xext[occulyr] - xpts[occulyr][ix])<1) {
                          expp = xpts[occulyr][ix]; break;
                        }
                      }
                      
                      if (expp >-50) {
                        h_xtcorstrips[occulyr]->Fill(expp, -1.0);
                        h_xytcorstrips[occulyr]->Fill(expp, -1.0);
                        for (int ix=0; ix<xpts[occulyr].size(); ix++) {
                          if (abs(expp - xpts[occulyr][ix])>2) { //choose all of them and then use maximum value before plot
                            h_xtcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
                          }
                        }
                        
                        for (int ix=0; ix<ypts[occulyr].size(); ix++) {
                          h_xytcorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
                        }
                      }
                    }
                  } //if (iiter==nmxiter-1)  
                } //else of if (occulyr >=nlayer) 
              } // if (nxtime>=nmnhits/*-ntcor*/ && nxtfail==0)
              
              for (int ix=0; ix<nlayer; ix++) { ytext[ix]= ytexter[ix] = 100;}
              
              int tmpntxy = 0;
              StraightLineFit ytimefit(iTimeSlopeConst, dist, ytime, timeserry2, yusedtime, occulyr, occulyr, ytcorstr, ytcorend, float(7.0));
              
              ytimefit.GetParameters(nytfail, yc0inters, timeyslope);
              //	    ytimefit.GetError(errcst, errlin, errcov);
              ytimefit.GetChisqure(nytime, yt0chi2);
              ytimefit.GetFitValues(ytext, ytdev, ytexter);
              
              //150101		  tmpytent[ij]= 0;
              
              if (nytfail==0 && isfill) {
                h_tchisqy->Fill(yt0chi2);
                if (nytime-2>0) {
                  h_treduchisqy->Fill(yt0chi2/(nytime-2));
                  double probx =TMath::Prob(yt0chi2, nytime-2);
                  h_ytprob->Fill(probx);
                  int ibin = getbinid(nytime, nprob, probs);
                  if (ibin>=0 && ibin<nprob) {h_ytnprob[ibin]->Fill(probx);}
                }
                h_tyndf->Fill(nytime);
              }
              
              if ((isfill || occulyr==nlayer) && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq) {
                for(int ij=0;ij<nlayer;ij++) {
                  if (dist[ij]<0 || ytime[ij]<-50) continue;
                  double dt1 = ytdev[ij]-biasinytime[ij];
                  if (abs(dt1)>4*sigmarpc || (!yusedtime[ij])) continue;
                  
                  //140923
                  //		  cout <<"occulyr "<< occulyr<<endl;
                  if (occulyr==nlayer) {
                    if (passytime[ij]) { 
                      ntotytime++; totytime +=ytime[ij];
                      for(int jk=ij+1;jk<nlayer;jk++){
                        if (dist[jk]<0 || ytime[jk]<-50 || abs(ytdev[jk]-biasinytime[jk])>4*sigmarpc || (!yusedtime[jk])  || (!passytime[jk])) continue;
                        timey_correl[ij][jk][iiter+1]->Fill(ytime[jk]-ytime[ij] - (dist[ij]-dist[jk])/cval);
                        if (iiter==0) { timey_correl[ij][jk][0]->Fill(rawytime[jk]-rawytime[ij] - (dist[ij]-dist[jk])/cval);}
                      } // for(jk=ij+1;jk<nlayer;jk++)
                      
                      if (passxtime[ij]) {
                        //X-Y correlation
                        
                        if (nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq && xtime[ij] >-50 &&
                            abs(xtdev[ij]-biasinxtime[ij])<4*sigmarpc && xusedtime[ij]) {
                          int ipix = 4*int((4*istrxtime[ij])/nstrip) + int((4*istrytime[ij])/nstrip);
                          if (ipix<0 || ipix>15) cout <<" ipix "<<ipix<<" "<<istrxtime[ij]<<" "<<istrytime[ij]<<endl;
                          
                          timexy_correl[ij][ipix][iiter+1]->Fill(ytime[ij]-xtime[ij]);
                          timexy_correl[ij][npixel][iiter+1]->Fill(ytime[ij]-xtime[ij]);
                          if (iiter==0) { 
                            timexy_correl[ij][ipix][0]->Fill(rawytime[ij]-rawxtime[ij]);
                            timexy_correl[ij][npixel][0]->Fill(rawytime[ij]-rawxtime[ij]);
                          }
                          //			  cout <<"iiter "<<itter<<endl;
                          if (iiter==0) { //compare timing of all combination of X/Y strips
                            int ixx = istrxtime[ij];
                            int iyy = istrytime[ij];
                            if (ixx>=0 && ixx <nstrip && iyy>=0 && iyy <nstrip) {
                              //			      cout <<"ijxy "<<ij<<" "<<ixx<<" "<<iyy<<" "<<rawytime0[ij]-rawxtime0[ij]<<" "<<rawytime3[ij]-rawxtime3[ij]<<endl;
                              indtimexy_correl[ij][ixx][iyy][0]->Fill(rawytime0[ij]-rawxtime0[ij]);
                              indtimexy_correl[ij][ixx][iyy][1]->Fill(rawytime3[ij]-rawxtime3[ij]);
#ifdef C217STRIP			      
                              indtimexy_prof[ij][0]->Fill(ixx+iyy, rawytime0[ij]-rawxtime0[ij]);
                              indtimexy_prof[ij][1]->Fill(ixx+iyy, rawytime3[ij]-rawxtime3[ij]);
#else
                              indtimexy_prof[ij][0]->Fill(iyy-ixx+nstrip, rawytime0[ij]-rawxtime0[ij]);
                              indtimexy_prof[ij][1]->Fill(iyy-ixx+nstrip, rawytime3[ij]-rawxtime3[ij]);			      
#endif
                            } else {
                              cout <<"wrong istr "<< ij<<" "<<ixx<<" "<<iyy<<endl;
                            }
                          }
                        }
                      } // if (passxtime[ij])
                    } // if (passytime[ij])
                  } //if (occulyr==nlayer) 
                } // ij
              } // if ((isfill || occulyr==nlayer) && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq)
              
              ntxyla = (xtimefit.GetLayerIds()<<nlayer) +  ytimefit.GetLayerIds();
              
              if (isfill) { 
                nTotalt++;
                for (int ix=0; ix<nlayer; ix++) {
                  if (Xpos[ix]>=-1 && Xpos[ix]<=nstrip && abs(xtdev[ix])<5 && xusedtime[ix]) {
                    for (int iy=ix+1; iy<nlayer; iy++) {
                      if (Xpos[iy]>=-1 && Xpos[iy]<=nstrip && abs(xtdev[iy])<5 && xusedtime[iy]) {
                        h_xtcorhits->Fill(ix, iy);
                      }
                    }
                    for (int iy=0; iy<nlayer; iy++) {
                      if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip && abs(ytdev[iy])<5 && yusedtime[iy]) {
                        h_xytcorhits->Fill(ix, iy);
                      }
                    }
                  }
                }
                for (int ix=0; ix<nlayer-1; ix++) {
                  if (Ypos[ix]>=-1 && Ypos[ix]<=nstrip && abs(ytdev[ix])<5 && yusedtime[ix]) {
                    for (int iy=ix+1; iy<nlayer; iy++) {
                      if (Ypos[iy]>=-1 && Ypos[iy]<=nstrip && abs(ytdev[iy])<5  && yusedtime[iy]) {
                        h_ytcorhits->Fill(ix, iy);
                      }
                    }
                  }
                }
              }
              
              if (nytime>=nmnhits/*-ntcor*/ && nytfail==0) {//4Nov
                if (yt0chi2/(nytime-2)<mxtimechisq) { 
                  
                  if (filloccu) {
                    if (iev>0 && tmptimerate >=tmpoldmutimeyrate+timebin) {
                      file_out <<"timeyrate "<<datafile<<" "<<iev<<" "<<nallmutimeyrate<<" "<<nmutimeyrate<<endl;
                      mutimeyrate->Fill(nmutimeyrate*(1.0/timebin));
                      mutimeyratex->Fill(nallmutimeyrate, nmutimeyrate*(1.0/timebin));
                      nmutimeyrate = 0;
                      nallmutimeyrate++;
                      tmpoldmutimeyrate = tmptimerate;
                    }
                    nmutimeyrate++;
                  }
                  
                  dir_cy[occulyr][iiter]->Fill(-cval*timeyslope);
                  int ibin = getbinid(nytime, nprob, probs);
                  if (ibin>=0 && ibin<nprob) {dir_cylay[occulyr][iiter][ibin]->Fill(-cval*timeyslope);}
                  dir_cy2[occulyr][iiter]->Fill(ytimefit.GetSlope2());
                  dir_c0y[occulyr][iiter]->Fill(yc0inters);
                }
                
                if (isfill) dir_cychi->Fill(yt0chi2, -1./timeyslope/cval);	  
                //		cout<<"nxtime "<<nxtime<<" "<< nmnhits<<" "<<nxtfail<<" "<<isfill<<endl;
                if (nxtime>=nmnhits/*-ntcor*/ && nxtfail==0) {
                  if (filloccu) { 
                    for (int iz=0; iz<nlayer; iz++) {
                      xlayer_allmutimemult[iz]->Fill(xptsall[iz].size());
                      ylayer_allmutimemult[iz]->Fill(yptsall[iz].size());
                    }
                  }
                  
                  double xt = (timexslope <0) ? log10(1+abs(1./timexslope/cval)) : -log10(1+abs(1./timexslope/cval));
                  double yt = (timeyslope <0) ? log10(1+abs(1./timeyslope/cval)) : -log10(1+abs(1./timeyslope/cval));
                  
                  dir_cxy->Fill(xt, yt);
                }
                
                if (occulyr >=nlayer) {
                  for (int ij=0; ij<nlayer; ij++) {
                    if (dist[ij]<0 || ytime[ij] <-50) continue;
                    if (yusedtime[ij]) { 
                      if (istrytime[ij]>=0 && istrytime[ij] <nstrip && passytmx[ij]) {
                        time_ystrreso[ij][istrytime[ij]][iiterrs]->Fill(ytdev[ij]);
                      }
                      ytime_exterr[ij][iiterrs]->Fill(ytexter[ij]);  
                    }
                    if (passytime[ij] && passytmx[ij]) {
                      if (yusedtime[ij]) {
#ifndef MONTECARLO
                        if (filloccu) {time_yreso_set[ij][iset]->Fill(ytdev[ij]);}
#endif
                        time_yreso[ij][iiterrs]->Fill(ytdev[ij]);
                        if (istime_xyreso[ij][iiterrs]) { 
                          time_xyreso[ij][iiterrs]->Fill(xtdev[ij], ytdev[ij]);
                        }
                      }
                      int iyy = max(0,min(yhits[ij],nmxtimehit)-1); 
                      time_mulyreso[ij][iiterrs][iyy]->Fill(ytdev[ij]);
                    }
#ifdef ISEFFICIENCY
                    if (yposEffPass[ij][iiterrs]) { 
                      total_yt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]); //11/12/2011
                      if((!yusedtime[ij]) || abs(ytdev[ij])>5*time_yrms[ij]) { 
                        inefficiency_yt[ij][iiterrs]->Fill(istrxtime[ij], istrytime[ij]);
                      }
                    }
#endif
                  }
                } else {
                  if (dist[occulyr]>=0 && ytime[occulyr] >-50) {
                    if (passytime[occulyr]) {
                      if (xtime[occulyr] >-50 && passxtime[occulyr]) { 
                        if ((abs(xtdev[occulyr])>5.0 || abs(ytdev[occulyr])>5.0) && 
                            abs(xtdev[occulyr])<30.0 && Ny>=9 && abs(ytdev[occulyr])<30.0 && Nx>=9) { 
                          double xx = xext[occulyr]*3.0+3.5;
                          double yy = yext[occulyr]*3.0+3.5;
                          
                          if (xext[occulyr]<1 || xext[occulyr]>30) {
                            xtdev_ytdev[0]->Fill(xtdev[occulyr], ytdev[occulyr]);
                          } else if (yext[occulyr]<1 || yext[occulyr]>30) {
                            xtdev_ytdev[1]->Fill(xtdev[occulyr], ytdev[occulyr]);
                          } else if ((abs(xext[occulyr]-5.0)<1.0 || abs(xext[occulyr]-11.0)<1.0 || abs(xext[occulyr]-16.0)<1.0 || abs(xext[occulyr]-22.0)<1.0 || abs(xext[occulyr]-28.0)<1.0) && (abs(yext[occulyr]-5.0)<1.0 || abs(yext[occulyr]-11.0)<1.0 || abs(yext[occulyr]-17.0)<1.0 || abs(yext[occulyr]-23.0)<1.0 || abs(yext[occulyr]-29)<1.0)) {
                            
                            xtdev_ytdev[2]->Fill(xtdev[occulyr], ytdev[occulyr]);
                          } else {
                            xtdev_ytdev[3]->Fill(xtdev[occulyr], ytdev[occulyr]);
                          }
                        }
                      }
                      
                      if (abs(ytdev[occulyr])>8.0) {
                        if (abs(ytdev[occulyr]) <30.0) {
                          if (ytdev[occulyr]>0) {
                            nxypos_xytdev[2]->Fill(xext[occulyr], yext[occulyr]);
                          } else {
                            nxypos_xytdev[3]->Fill(xext[occulyr], yext[occulyr]);
                          }
                          
                          if (xtime[occulyr] >-50 && passxtime[occulyr]) {
                            if (abs(xtdev[occulyr])>8.0 &&abs(xtdev[occulyr])<30.0) {
                              if (xtdev[occulyr]>0) {
                                if (ytdev[occulyr]>0) { 
                                  nxypos_xytdev[4]->Fill(xext[occulyr], yext[occulyr]);
                                } else {
                                  nxypos_xytdev[5]->Fill(xext[occulyr], yext[occulyr]);
                                }
                              } else {
                                if (ytdev[occulyr]>0) { 
                                  nxypos_xytdev[6]->Fill(xext[occulyr], yext[occulyr]);
                                } else {
                                  nxypos_xytdev[7]->Fill(xext[occulyr], yext[occulyr]);
                                }
                              }
                            }
                          }
                        }
                      }
                      if (yusedtime[occulyr]) {
                        widthx_timex[occulyr][iiter]->Fill(ytdev[occulyr], widthy[occulyr]);
                        xstr_ytdev[occulyr][iiter]->Fill(xext[occulyr], ytdev[occulyr]);
                        ystr_ytdev[occulyr][iiter]->Fill(yext[occulyr], ytdev[occulyr]);
                        
                        if (abs(ytdev[occulyr])<20) {
                          nxystr_ytdev[occulyr][iiter]->Fill(xext[occulyr], yext[occulyr], 1.);
                          xystr_ytdev[occulyr][iiter]->Fill(xext[occulyr], yext[occulyr], ytdev[occulyr]+20.0);
                        }
                      }
                      
                      if (iiter==nmxiter-1) {
                        double diff = -30;
                        for (int jkl=0; jkl<6; jkl++) {
                          switch(jkl) {
                          case 0 : diff = rawytime0[occulyr]-ytext[occulyr]; break;
                          case 1 : diff = rawytime1[occulyr]-ytext[occulyr]; break;
                          case 2 : diff = rawytime2[occulyr]-ytext[occulyr]; break;
                          case 3 : diff = timesy[occulyr]-ytext[occulyr]; break;
                          case 4 : diff = rawytime3[occulyr]-ytext[occulyr]; break;
                          default : diff = ytime[occulyr]-ytext[occulyr]; break;
                          }
                          if (abs(diff)<20 && yusedtime[occulyr]) {
                            nxystr_ytdev[occulyr][nmxiter+jkl]->Fill(xext[occulyr], yext[occulyr], 1.);
                            xystr_ytdev[occulyr][nmxiter+jkl]->Fill(xext[occulyr], yext[occulyr], diff+20.0);
                          }
                        }
                        
                        int nsize=ypts[occulyr].size();
                        if (nsize>=1 && nsize<=nmxtimehit) {
                          xpos_ydev[occulyr][nsize-1]->Fill(xposinstr[occulyr], Ydev[occulyr]);
                          ypos_ydev[occulyr][nsize-1]->Fill(yposinstr[occulyr], Ydev[occulyr]);
                          xpos_ytdev[occulyr][nsize-1]->Fill(xposinstr[occulyr], ytdev[occulyr]);
                          ypos_ytdev[occulyr][nsize-1]->Fill(yposinstr[occulyr], ytdev[occulyr]);
                          
                          ypos_ytdev_str[occulyr][nsize-1]->Fill(yposinstr[occulyr], diff);
                          ypos_ytdev_glb[occulyr][nsize-1]->Fill(yext[occulyr], diff);
                          xpos_ytdev_str[occulyr][nsize-1]->Fill(xposinstr[occulyr], diff);
                          xpos_ytdev_glb[occulyr][nsize-1]->Fill(xext[occulyr], diff);
                        }
                      }
                    }
                    if (yusedtime[occulyr]) { 
                      if (istrytime[occulyr]>=0 && istrytime[occulyr] <nstrip) {
#ifdef TIMESLOPE
                        time_yslope_pr[occulyr][istrytime[occulyr]][iiter]->Fill(xext[occulyr], ytdev[occulyr]);
#endif
                        if (passytmx[occulyr]) {time_ystrreso[occulyr][istrytime[occulyr]][iiterrs]->Fill(ytdev[occulyr]);}
                      }
                      ytime_exterr[occulyr][iiterrs]->Fill(ytexter[occulyr]);
                    }
                    if (passytime[occulyr] && passytmx[occulyr]) {
                      if (yusedtime[occulyr]) {
                        time_yreso[occulyr][iiterrs]->Fill(ytdev[occulyr]);
                        if (istime_xyreso[occulyr][iiterrs]) { 
                          time_xyreso[occulyr][iiterrs]->Fill(xtdev[occulyr], ytdev[occulyr]);
                        }
                      }
                      int iyy = max(0,min(yhits[occulyr],nmxtimehit)-1); 
                      time_mulyreso[occulyr][iiterrs][iyy]->Fill(ytdev[occulyr]);
                    }
#ifdef ISEFFICIENCY		    
                    if (yposEffPass[occulyr][iiterrs]) { 
                      total_yt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]); //11/12/2011
                      if((!yusedtime[occulyr]) || abs(ytdev[occulyr])<5*time_yrms[occulyr]) {
                        inefficiency_yt[occulyr][iiterrs]->Fill(istrxtime[occulyr], istrytime[occulyr]);
                      }
                    }
#endif
                  } //if (dist[occulyr] >=0 && xtime[occulyr] >-50.0)
                  if (iiter==nmxiter-1) { //Correlated hits
                    if (abs(timesy[occulyr] - ytext[occulyr])<10.0) {
                      double expp = -100;
                      for (int ix=0; ix<ypts[occulyr].size(); ix++) {
                        if (abs(yext[occulyr] - ypts[occulyr][ix])<1) {
                          expp = ypts[occulyr][ix]; break;
                        }
                      }
                      if (expp >-50) { 
                        h_ytcorstrips[occulyr]->Fill(expp, -1.0);
                        h_yxtcorstrips[occulyr]->Fill(expp, -1.0);
                        for (int ix=0; ix<ypts[occulyr].size(); ix++) {
                          if (abs(expp - ypts[occulyr][ix])>2) {//choose all of them and then use maximum value before plot
                            h_ytcorstrips[occulyr]->Fill(expp, ypts[occulyr][ix]);
                          }
                        }
                        
                        for (int ix=0; ix<xpts[occulyr].size(); ix++) {
                          h_yxtcorstrips[occulyr]->Fill(expp, xpts[occulyr][ix]);
                        }
                      }
                    }
                  } //if (iiter==nmxiter-1) 
                } // else of if (occulyr >=nlayer) 
              } //if (nytime>=nmnhits/*-ntcor*/ && nytfail==0)
            } // if (Nx>=4 && Ny>=4 && xchi2/(Nx-2)<2.0 && ychi2/(Ny-2)<2.0 && nxfail==0 && nyfail==0); Position fit cut
            if (isTiming && isfill &&  nxfail==0 &&  nyfail==0 && nxtime>2 && nytime>2) {
              
              if (isfill && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq) {
                if (nxtfail==0 && nxtime >nmnhits && xt0chi2/(nxtime-2)<mxtimechisq ) {
                  for(int ij=0;ij<nlayer;ij++){
                    
                    if (dist[ij]<0 || xtime[ij] <-50 || (!xusedtime[ij])) continue;
                    double dt1 = xtdev[ij]-biasinxtime[ij];
                    
                    if (abs(dt1)>3*sigmarpc) continue;
                    for(int jk=0;jk<nlayer;jk++){
                      
                      if (dist[jk]<0 || xtime[jk] <-50 || (!xusedtime[jk])) continue;
                      double dt2 = xtdev[jk]-biasinxtime[jk];
                      if (abs(dt2) >3*sigmarpc) continue;
                      deltatcov2[ij][jk] += dt1*dt2;
                      deltatCount2[ij][jk] +=1 ;
                    }
                  }// ij
                }
                
                for(int ij=0;ij<nlayer;ij++) {
                  if (dist[ij]<0 || ytime[ij] <-50 || (!yusedtime[ij])) continue;
                  double dt1 = ytdev[ij]-biasinytime[ij];
                  if (abs(dt1)>3*sigmarpc || (!yusedtime[ij])) continue;
                  for(int jk=0;jk<nlayer;jk++){
                    
                    if (dist[jk]<0 || ytime[jk] <-50 || (!yusedtime[jk]) ) continue;
                    //                                      if(jk != ij )continue;
                    double dt2 = ytdev[jk]-biasinytime[jk];
                    if (abs(dt2) >3*sigmarpc) continue;
                    deltatcovy[ij][jk] += dt1*dt2;
                    deltatCounty[ij][jk] +=1 ;
                  }
                }// ij
              } //if (isfill && nytfail==0 && nytime >nmnhits && yt0chi2/(nytime-2)<mxtimechisq)
              
              txslop = -1./timexslope/cval;
              tyslop = -1./timeyslope/cval;
              
              errtimexslope =sqrt(errlin_tx);
              
              if (nxtime>11 || nytime>11 || gRandom->Uniform()<1.e-5) T2->Fill(); //07/09/2011
            }
            
            if (ntcor==0) {
              passed_strip[iiter]->Fill(Nx, 1.);
              passed_strip[iiter]->Fill(nlayer+2+Ny, 1.);
              passed_strip[iiter]->Fill(2*(nlayer+2)+nxtime, 1.);
              passed_strip[iiter]->Fill(3*(nlayer+2)+nytime, 1.);
            }
            // cout<<iev<<"		"<<"XXXXXXXXXXXXXXXXXXXXXXX"<<endl;
            //--------------------- Don't change anything here --------------
            if((iev%100000)==1) 
              { cout<<"Processed " <<ntotal<<" ("<<iev<<") events so far for iteration# "<< iiter<<" occu "<<occulyr<<endl; }
            
            //#if ndefined(MONTECARLO) && defined(NOISE_SIM)
#if defined(NOISE_SIM)
            //#ifndef MONTECARLO
            //Generate root files for noise and subsequently use for simulation
            if (isfill && nxfail==0 && Nx >5 && xchi2<50 && nyfail==0 && Ny >5 && ychi2<50) {
              NoisefileOut->cd();
              //Copy isolated events
              bool isxnoise = false;
              bool isynoise = false;
              vector<int> xx2ptsall[nlayer];
              vector<int> yy2ptsall[nlayer];
              nseltot++;
              
              if (nxfail==0 && Nx >5 && xchi2<50) {
                isxnoise = true;
                for (int ij=0; ij<nlayer; ij++) {
                  for (int ix=0; ix<xptsall[ij].size(); ix++) {
                    xlayer_seloccu[ij]->Fill(xptsall[ij][ix]);
                    if (xptsall[ij].size()<4) {xlayer_sel2occu[ij]->Fill(xptsall[ij][ix]);}
                    
                    if (abs(xext[ij] - xptsall[ij][ix])>2.0) {
                      xx2ptsall[ij].push_back(xptsall[ij][ix]);
                      xlayer_noiseoccu[ij]->Fill(xptsall[ij][ix]);
                    } 
                  }
                }
              }
              
              if (nyfail==0 && Ny >5 && ychi2<50) {
                isynoise = true;
                for (int ij=0; ij<nlayer; ij++) {
                  for (int iy=0; iy<yptsall[ij].size(); iy++) {
                    ylayer_seloccu[ij]->Fill(yptsall[ij][iy]);
                    if (yptsall[ij].size()<4) {ylayer_sel2occu[ij]->Fill(yptsall[ij][iy]);}
                    
                    if (abs(yext[ij] - yptsall[ij][iy])>2.0) {
                      yy2ptsall[ij].push_back(yptsall[ij][iy]);
                      ylayer_noiseoccu[ij]->Fill(yptsall[ij][iy]);
                    } 
                  }
                }
              }
              
              for (int ij=0; ij<nlayer; ij++) {
                for (int ix=0; ix<xptsall[ij].size(); ix++) {
                  for (int iy=0; iy<yptsall[ij].size(); iy++) {
                    raw_seloccu[ij]->Fill(xptsall[ij][ix], yptsall[ij][iy]);
                  }
                }
                
                for (int ix=0; ix<xx2ptsall[ij].size(); ix++) {
                  for (int iy=0; iy<yy2ptsall[ij].size(); iy++) {
                    raw_noiseoccu[ij]->Fill(xx2ptsall[ij][ix], yy2ptsall[ij][iy]);
                  }
                }
              }
              
              //Store hits after removing muon hits, simple noise
              if (isxnoise && isynoise) { 
                
                int nxmulall = 0; int nymulall = 0;
                //	      int xhitsall[nlayer]={0};
                for (int ij=0; ij<nlayer; ij++) {
                  //		xhits[ij] = xptsall[ij].size();
                  //		yhits[ij] = yptsall[ij].size();
                  //		nmulall += xhitsall[ij] = xptsall[ij].size()  + yptsall[ij].size();
                  nxmulall += xptsall[ij].size();
                  nymulall += yptsall[ij].size();
                }
                n_multall[0]->Fill(nxmulall);
                n_multall[1]->Fill(nymulall);
                
                for (int ij=0; ij<nlayer; ij++) {
                  //for (int jk=0; jk<=xptsall[ij].size(); jk++) {
                  total_vs_indmulall[0][ij]->Fill(nxmulall, xptsall[ij].size()); //check this
                  //}
                  //for (int jk=0; jk<=yptsall[ij].size(); jk++) {
                  total_vs_indmulall[1][ij]->Fill(nymulall,yptsall[ij].size()); //check this
                  //}
                }
                
                for (int ij=0; ij<nlayer; ij++) {
                  layer_multall[0][ij]->Fill(xptsall[ij].size());
                  layer_multall[1][ij]->Fill(yptsall[ij].size());
                  
                  for (int jk=0; jk<xptsall[ij].size(); jk++) { // X
                    layer_hitsall[0][ij]->Fill(xptsall[ij][jk]);
                    hist_correlall[0][ij]->Fill(xptsall[ij][jk], xptsall[ij][jk]);
                    for (int kl=jk+1; kl<xptsall[ij].size(); kl++) {
                      hist_correlall[0][ij]->Fill(xptsall[ij][jk], xptsall[ij][kl]);
                    }
                  }
                  
                  for (int jk=0; jk<yptsall[ij].size(); jk++) { // Y
                    layer_hitsall[1][ij]->Fill(yptsall[ij][jk]);
                    hist_correlall[1][ij]->Fill(yptsall[ij][jk], yptsall[ij][jk]);
                    for (int kl=jk+1; kl<yptsall[ij].size(); kl++) {
                      hist_correlall[1][ij]->Fill(yptsall[ij][jk], yptsall[ij][kl]);
                    }
                  }
                  
                  for (int jk=ij+1; jk<nlayer; jk++) {
                    mult_correlall[0][jk]->Fill(xptsall[ij].size(), xptsall[jk].size());
                    mult_correlall[1][jk]->Fill(yptsall[ij].size(), yptsall[jk].size());
                  }
                  
                  if (xptsall[ij].size()>0) {layer_timeall[0][ij]->Fill(timesx[ij]);}
                  if (yptsall[ij].size()>0) {layer_timeall[1][ij]->Fill(timesy[ij]);}
                }
                
                //Hits excluding muon trajectory
                
                nxmulall = nymulall = 0;
                ntot_uncsim++;
                //	      for (int ij=0; ij<nlayer; ij++) { xhitsall[ij]=0;}
                for (int ij=0; ij<nlayer; ij++) {
                  //		xhits[ij] = xptsall[ij].size();
                  //		yhits[ij] = yptsall[ij].size();
                  //		nmulall += xhitsall[ij] = xptsall[ij].size()  + yptsall[ij].size();
                  nxmulall += xx2ptsall[ij].size();
                  nymulall += yy2ptsall[ij].size();
                }
                
                n_mult[0]->Fill(nxmulall);
                n_mult[1]->Fill(nymulall);
                
                for (int ij=0; ij<nlayer; ij++) {
                  //for (int jk=0; jk<=xx2ptsall[ij].size(); jk++) {
                  total_vs_indmul[0][ij]->Fill(nxmulall, xx2ptsall[ij].size()); //check this
                  //}
                  //for (int jk=0; jk<=yy2ptsall[ij].size(); jk++) {
                  total_vs_indmul[1][ij]->Fill(nymulall, yy2ptsall[ij].size()); //check this
                  //}
                }
                
                for (int ij=0; ij<nlayer; ij++) {
                  layer_mult[0][ij]->Fill(xx2ptsall[ij].size());
                  layer_mult[1][ij]->Fill(yy2ptsall[ij].size());
                  
                  for (int jk=0; jk<xx2ptsall[ij].size(); jk++) { // X&&Y
                    layer_hits[0][ij]->Fill(xx2ptsall[ij][jk]);
                    hist_correl[0][ij]->Fill(xx2ptsall[ij][jk], xx2ptsall[ij][jk]);
                    for (int kl=jk+1; kl<xx2ptsall[ij].size(); kl++) {
                      hist_correl[0][ij]->Fill(xx2ptsall[ij][jk], xx2ptsall[ij][kl]);
                    }
                  }
                  for (int jk=0; jk<yy2ptsall[ij].size(); jk++) { // X&&Y
                    layer_hits[1][ij]->Fill(yy2ptsall[ij][jk]);
                    hist_correl[1][ij]->Fill(yy2ptsall[ij][jk], yy2ptsall[ij][jk]);
                    for (int kl=jk+1; kl<yy2ptsall[ij].size(); kl++) {
                      hist_correl[1][ij]->Fill(yy2ptsall[ij][jk], yy2ptsall[ij][kl]);
                    }
                  }
                  
                  for (int jk=ij+1; jk<nlayer; jk++) {
                    mult_correl[0][jk]->Fill(xx2ptsall[ij].size(), xx2ptsall[jk].size());
                    mult_correl[1][jk]->Fill(yy2ptsall[ij].size(), yy2ptsall[jk].size()); 
                  }
                  
                  if (xx2ptsall[ij].size()>0) {layer_time[0][ij]->Fill(timesx[ij]);}
                  if (yy2ptsall[ij].size()>0) {layer_time[1][ij]->Fill(timesy[ij]);}
                }
              }
            }
            //#endif  // #if ndefined(MONTECARLO) && defined(NOISE_SIM)
#endif
          }  // EVENT loop
          cout<<"missing hits from X trigger layer"<<"   "<<nentry<<"    "<<nx4<<"  "<<nx5<<"   "<<nx6<<"  "<<nx7<<endl;
          cout<<"missing hits from Y trigger layer"<<"   "<<nentry<<"    "<<ny4<<"  "<<ny5<<"   "<<ny6<<"  "<<ny7<<endl;
          cout<<"missing hits from X and Y trigger layer"<<"   "<<nentry<<"    "<<nxy4<<"  "<<nxy5<<"   "<<nxy6<<"  "<<nxy7<<endl;
          cout<<"missing hits from X or Y trigger layer"<<"   "<<nentry<<"    "<<nxory4<<"  "<<nxory5<<"   "<<nxory6<<"  "<<nxory7<<endl;
          cout<<"missing hits from X and Y final trigger"<<"   "<<nentry<<"    "<<nxytrig<<endl;
          cout<<"latched vs tdc counts"<<endl;
          for(int jj=0;jj<nlayer;jj++) {
            for(int kk=0;kk<8;kk++) {
              cout<<jj<<"   "<<kk<<"   "<<ntotxext[jj][kk]<<"    "<<ntxtime[jj][kk]<<"   "<<ntotyext[jj][kk]<<"    "<<ntytime[jj][kk]<<endl;
            }
          }
          
          cout<<"tdc vs latched"<<endl;
          for(int jj=0;jj<nlayer;jj++) {
            for(int kk=0;kk<8;kk++) {
              cout<<jj<<"   "<<kk<<"   "<<ntxtime2[jj][kk]<<"    "<<ntotxext2[jj][kk]<<"   "<<ntytime2[jj][kk]<<"   "<<ntotyext2[jj][kk]<<"    "<<endl;
            }
          }
          
          fileIn->cd();
#ifdef MONTECARLO
          delete T1;
#else
          // #elif DLEVEL
          delete event_tree;
          //	  delete cau_tree;
#endif
          delete fileIn;
        } // while(!(file_db.eof()))
        
        if (ntotal<1) ntotal=1;
        file_db.close();
        
        if(isfill) {
          for(int ij=0;ij<nstrip/2;ij++) {
            for(int jk=0;jk<nstrip/2;jk++) {
              double scatmean1 = pixel_scatang1[ij][jk]->GetMean();
              // double scatmed1 =  Median(pixel_scatang1[ij][jk]);
              pixel_scatmean1->SetBinContent(ij,jk,scatmean1);
              
              double scatmean = pixel_scatang[ij][jk]->GetMean();
              // double scatmed =  Median(pixel_scatang[ij][jk]);
              pixel_scatmean->SetBinContent(ij,jk,scatmean);
            }
          }
        }
        
        //#if ndefined(MONTECARLO) && defined(NOISE_SIM)
#if defined(NOISE_SIM)
        //#ifndef MONTECARLO
        if (isfill) { 
          NoisefileOut->cd(); 
          // Normalise layer multiplity for each individual value of total multiplicity
          for (int ixy=0; ixy<2; ixy++) { 
            for (int ij=0; ij<nlayer; ij++) {
              for (int jk=0; jk<total_vs_indmulall[ixy][ij]->GetNbinsX(); jk++) {
                double ntt = 0.0;
                for (int kl=0; kl<total_vs_indmulall[ixy][ij]->GetNbinsY(); kl++) { 
                  ntt += total_vs_indmulall[ixy][ij]->GetBinContent(jk+1,kl+1);
                }
                if (ntt<1.0) ntt=1.0;
                for (int kl=0; kl<total_vs_indmulall[ixy][ij]->GetNbinsY(); kl++) { 
                  total_vs_indmulall[ixy][ij]->SetBinContent(jk+1, kl+1, total_vs_indmulall[ixy][ij]->GetBinContent(jk+1,kl+1)/ntt); // (max(1, ntot_uncsim))/*ntt*/);
                }
              }
              
              for (int jk=0; jk<total_vs_indmul[ixy][ij]->GetNbinsX(); jk++) {
                double ntt = 0; 
                for (int kl=0; kl<total_vs_indmul[ixy][ij]->GetNbinsY(); kl++) {
                  ntt += total_vs_indmul[ixy][ij]->GetBinContent(jk+1,kl+1);
                }
                if (ntt<1.0) ntt=1.0;
                for (int kl=0; kl<total_vs_indmul[ixy][ij]->GetNbinsY(); kl++) { 
                  total_vs_indmul[ixy][ij]->SetBinContent(jk+1, kl+1, total_vs_indmul[ixy][ij]->GetBinContent(jk+1,kl+1)/ntt); //(max(1, ntot_uncsim))/*ntt*/);
                }
              }
            }
          }
          
          double scal = 1./(max(1, ntot_uncsim));
          double scalunc = 1./(max(1, ntot_uncsim));
          cout <<"scal "<<scal<<endl;
          
          for (int ixy=0; ixy<2; ixy++) {
            n_multall[ixy]->Scale(scal);
            n_mult[ixy]->Scale(scalunc);
            
            for (int ij=0; ij<nlayer; ij++) {
              layer_hitsall[ixy][ij]->Scale(scal); 
              layer_multall[ixy][ij]->Scale(scal);  
              layer_timeall[ixy][ij]->Scale(scal); 
              
              mult_correlall[ixy][ij]->Scale(scal);
              hist_correlall[ixy][ij]->Scale(scal);
              
              layer_hits[ixy][ij]->Scale(scalunc); 
              layer_mult[ixy][ij]->Scale(scalunc);  
              layer_time[ixy][ij]->Scale(scalunc); 
              
              mult_correl[ixy][ij]->Scale(scalunc);
              hist_correl[ixy][ij]->Scale(scalunc);
            }
          }
          NoisefileOut->Write(); 
          NoisefileOut->Close();
        }
#endif
        fileOut->cd();
        //Reject moisy channels
        if (filloccu) {
          filloccu = false;
          
          TF1* fitxy[nlayer]={0};
          TH1F* xy_occu[nlayer]={0};
          
          double xycount[2][nlayer]={0}; //Sum of X and Y strip for all layer
          
          gStyle->SetOptStat(0);
          gStyle->SetStatW(.36);
          gStyle->SetStatH(.20);
          gStyle->SetStatY(.99);
          gStyle->SetStatX(.99);
          gStyle->SetOptFit(0);
          gStyle->SetOptLogy(1);
          gStyle->SetPadBottomMargin(0.11);
          gStyle->SetPadTopMargin(0.07);
          gStyle->SetPadLeftMargin(0.10);
          gStyle->SetPadRightMargin(0.11);
          gStyle->SetPaintTextFormat("g");
          gStyle->SetOptTitle(1);
          
          for (int ixy=0; ixy<8; ixy++) {
            ps.NewPage();
            gStyle->SetPadRightMargin(0.01);
            gStyle->SetOptLogy(0);
            
            TCanvas *c4=new TCanvas ("c4","Strip occupancy",500,700);
            c4->Divide(3,4);
            
            int nxx = (ixy<2) ? max(1, ntotal) : max(1, nseltot);
            double ascl=100./nxx;
            
            for (int il=0; il<nlayer; il++) {
              switch (ixy) { 
              case 0 : xy_occu[il] = (TH1F*)xlayer_occu[il]->Clone(); break;
              case 1 : xy_occu[il] = (TH1F*)ylayer_occu[il]->Clone(); break;
              case 2 : xy_occu[il] = (TH1F*)xlayer_seloccu[il]->Clone(); break;
              case 3 : xy_occu[il] = (TH1F*)ylayer_seloccu[il]->Clone(); break;	
              case 4 : xy_occu[il] = (TH1F*)xlayer_sel2occu[il]->Clone(); break;
              case 5 : xy_occu[il] = (TH1F*)ylayer_sel2occu[il]->Clone(); break;	
              case 6 : xy_occu[il] = (TH1F*)xlayer_noiseoccu[il]->Clone(); break;
              case 7 : xy_occu[il] = (TH1F*)ylayer_noiseoccu[il]->Clone(); break;
              default : xy_occu[il] = (TH1F*)xlayer_occu[il]->Clone(); break;
              }
              
              xy_occu[il]->Scale(ascl);
              c4->cd(il+1);
              
              sprintf(name, "fitxy_%i_%i", il, ixy);
              fitxy[il] = new TF1(name, fitspec, 0, nusedstrip, 4);
              
              //  if(ixy<2) {
              double hgh ;//= 0.5*xy_occu[il]->GetMaximum();
              double parx[4] ={hgh, 0., hgh, 1.};
              
              if(ixy<2) {
                hgh = 0.5*xy_occu[il]->GetMaximum();
                fitxy[il]->SetParLimits(0, 0.0, 2.0*hgh);
                fitxy[il]->SetParLimits(1, -1.0, 1.0);
                fitxy[il]->SetParLimits(2, 0.5*hgh, 4.0*hgh);
                
                fitxy[il]->FixParameter(3, 1.0);
                fitxy[il]->SetParameters(parx);
              }
              
              xy_occu[il]->SetMarkerStyle(24);
              xy_occu[il]->SetMarkerSize(0.8);
              xy_occu[il]->GetXaxis()->SetTitle("Strip No");
              xy_occu[il]->GetXaxis()->SetTitleOffset(.8);
              xy_occu[il]->GetXaxis()->SetTitleSize(.06);
              xy_occu[il]->GetXaxis()->CenterTitle();
              xy_occu[il]->GetYaxis()->SetTitle("Occupancy (%)");
              xy_occu[il]->GetYaxis()->SetTitleOffset(.8);
              xy_occu[il]->GetYaxis()->SetTitleSize(.06);
              xy_occu[il]->GetYaxis()->CenterTitle();
              
              xy_occu[il]->GetXaxis()->SetLabelSize(0.06);
              xy_occu[il]->GetYaxis()->SetLabelSize(0.06);
              xy_occu[il]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
              if(ixy<2) {
                xy_occu[il]->Fit(name, "WW:BRQ");
                fitxy[il]->GetParameters(parx);
                
                double xx[2]={0,0};
                for (int istr=0; istr<nusedstrip; istr++) {
                  xx[0] = xy_occu[il]->GetBinCenter(istr+1);
                  double yy = xy_occu[il]->GetBinContent(istr+1);
                  // 		cout <<"ixy "<< ixy<<" "<<il<<" "<<ntotal<<" "<<istr<<" "<<xx[0]<<" "<<yy<<" "<<xy_occu[il]->GetBinContent(istr+1)<<endl;
                  xycount[ixy][il] +=yy;
                  
                  if (/*ixy==0 &&*/ yy >2.8*fitspec(xx, parx)) { posinuse[ixy][il][istr]=false; cout<<"XXXXXXXXXXXXXXXXXXXXXX"<<endl;}
                }
              } else {
                xy_occu[il]->Draw();
              }
            }
            c4->Update();
            if (c4) { delete c4; c4=0;}
            
            for (int il=0; il<nlayer; il++) {
              if (fitxy[il]) {delete fitxy[il]; fitxy[il]=0;}
              if (xy_occu[il]) {delete xy_occu[il]; xy_occu[il]=0;}
            }
          } // for (int ixy=0; ixy<2; ixy++)
          
          ps.NewPage();
          gStyle->SetOptStat(1110);
          gStyle->SetOptTitle(1);
          gStyle->SetOptLogy(1);
          gStyle->SetStatW(.40); //40);
          gStyle->SetStatH(.24); //30);
          gStyle->SetTitleFontSize(0.07);
          gStyle->SetPadLeftMargin(0.09);
          gStyle->SetPadBottomMargin(0.06);
          gStyle->SetPadTopMargin(0.09); //(0.03);
          gStyle->SetPadRightMargin(0.01);
          TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
          c4c->Divide(3,4);
          
          //  ps.NewPage();
          gStyle->SetOptStat(1100);
          gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            xlayer_allmult[ij]->Scale(1./ntotal);
            xlayer_allmult[ij]->GetXaxis()->SetLabelSize(0.054);
            xlayer_allmult[ij]->GetYaxis()->SetLabelSize(0.054);
            xlayer_allmult[ij]->GetYaxis()->SetRangeUser(1.e-5, 0.8);
            xlayer_allmult[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
            xlayer_allmult[ij]->SetLineColor(1); xlayer_allmult[ij]->Draw();
          }
              
          c4c->Update();
          gStyle->SetStatY(.87); gStyle->SetStatTextColor(2);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            ylayer_allmult[ij]->Scale(1./ntotal);
            ylayer_allmult[ij]->SetLineColor(2); ylayer_allmult[ij]->Draw("sames");
          }
          c4c->Update();
          gStyle->SetStatY(.75); gStyle->SetStatTextColor(3);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            xlayer_mult[ij]->Scale(1./ntotal);
            xlayer_mult[ij]->SetLineColor(3); xlayer_mult[ij]->Draw("sames");
          }
          c4c->Update();
          gStyle->SetStatY(.64); gStyle->SetStatTextColor(4);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            ylayer_mult[ij]->Scale(1./ntotal);
            ylayer_mult[ij]->SetLineColor(4); ylayer_mult[ij]->Draw("sames");
          }
          
          c4c->Update();
          
          ps.NewPage();
          gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            xlayer_allmumult[ij]->Scale(1./ntotal);
            xlayer_allmumult[ij]->GetXaxis()->SetLabelSize(0.054);
            xlayer_allmumult[ij]->GetYaxis()->SetLabelSize(0.054);
            xlayer_allmumult[ij]->GetYaxis()->SetRangeUser(1.e-5, 0.8);
            xlayer_allmumult[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
            xlayer_allmumult[ij]->SetLineColor(1); xlayer_allmumult[ij]->Draw();
          }
          c4c->Update();
          gStyle->SetStatY(.87); gStyle->SetStatTextColor(2);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            ylayer_allmumult[ij]->Scale(1./ntotal);
            ylayer_allmumult[ij]->SetLineColor(2); ylayer_allmumult[ij]->Draw("sames");
          }
          c4c->Update();
          gStyle->SetStatY(.75); gStyle->SetStatTextColor(3);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            xlayer_allmutimemult[ij]->Scale(1./ntotal);
            xlayer_allmutimemult[ij]->SetLineColor(3); xlayer_allmutimemult[ij]->Draw("sames");
          }
          c4c->Update();
          gStyle->SetStatY(.64); gStyle->SetStatTextColor(4);
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            ylayer_allmutimemult[ij]->Scale(1./ntotal);
            ylayer_allmutimemult[ij]->SetLineColor(4); ylayer_allmutimemult[ij]->Draw("sames");
          }
          
          c4c->Update();
          
          ps.NewPage();
          gStyle->SetStatTextColor(1); gStyle->SetStatY(.99);
          
          gStyle->SetOptStat(1110);
          gStyle->SetStatW(.36);
          gStyle->SetStatH(.20);
          gStyle->SetStatY(.99);
          gStyle->SetStatX(.99);
          gStyle->SetOptFit(0);
          gStyle->SetOptLogy(1);
          gStyle->SetPadBottomMargin(0.12);
          gStyle->SetPadTopMargin(0.09);
          gStyle->SetPadLeftMargin(0.12);
          gStyle->SetPadRightMargin(0.01);
          gStyle->SetPaintTextFormat("g");
          gStyle->SetOptTitle(1);
          
#ifndef MONTECARLO
          if (plot_level>90) { 
            TCanvas *c44=new TCanvas ("c44","Strip occupancy",500,700);
            c44->Divide(3,4);
            for (int ix=0; ix<=iset; ix++) {
              if (ix>0) ps.NewPage();
              gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);
              
              for (int il=0; il<nlayer; il++) {
                c44->cd(il+1);
                strp_xmult_set[il][ix]->SetLineColor(1);
                strp_xmult_set[il][ix]->Scale(1./nCount);
                strp_xmult_set[il][ix]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
                strp_xmult_set[il][ix]->GetXaxis()->SetTitle("Multiplicity");
                strp_xmult_set[il][ix]->GetYaxis()->SetTitle("Normalised Entry");
                strp_xmult_set[il][ix]->GetXaxis()->SetTitleOffset(.75);
                strp_xmult_set[il][ix]->GetXaxis()->SetTitleSize(.06);
                strp_xmult_set[il][ix]->GetXaxis()->CenterTitle();
                strp_xmult_set[il][ix]->GetXaxis()->SetLabelSize(.06);
                strp_xmult_set[il][ix]->GetXaxis()->SetLabelOffset(-0.001);
                strp_xmult_set[il][ix]->GetYaxis()->SetTitleOffset(0.9);
                strp_xmult_set[il][ix]->GetYaxis()->SetTitleSize(.06);
                strp_xmult_set[il][ix]->GetYaxis()->CenterTitle();
                strp_xmult_set[il][ix]->GetYaxis()->SetLabelSize(.05);
                strp_xmult_set[il][ix]->GetYaxis()->SetLabelOffset(-0.001);
                strp_xmult_set[il][ix]->Draw();
              }
              c44->Update();
              gStyle->SetStatY(0.80); gStyle->SetStatTextColor(2);
              for (int il=0; il<nlayer; il++) {
                c44->cd(il+1);
                strp_ymult_set[il][ix]->SetLineColor(2);
                strp_ymult_set[il][ix]->Scale(1./nCount);
                strp_ymult_set[il][ix]->Draw("sames");
              }
              c44->Update();
            }
            
            //	  ps.NewPage();
            gStyle->SetOptStat(0);  gStyle->SetTitleFontSize(0.05);
            gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
            for (int ix=0; ix<=iset; ix++) {
              ps.NewPage();
              c44->cd();
              gPad->SetLogy(0); 
              gPad->SetLeftMargin(0.05);
              gPad->SetRightMargin(0.07);
              gStyle->SetPadTopMargin(0.05);
              
              strp_count_set[ix]->Scale(100./nCount);
              strp_count_set[ix]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
              strp_count_set[ix]->GetXaxis()->SetLabelSize(.040);
              strp_count_set[ix]->GetXaxis()->SetLabelOffset(0.003);
              strp_count_set[ix]->GetYaxis()->SetLabelSize(.030);
              strp_count_set[ix]->GetYaxis()->SetLabelOffset(0.003);
              strp_count_set[ix]->SetMarkerSize(0.75);
              strp_count_set[ix]->GetZaxis()->SetLabelSize(.030);
              for (int ij=0; ij<nusedstrip; ij++) {
                strp_count_set[ix]->GetYaxis()->SetBinLabel(ij+1, labels[ij]);
              }
              strp_count_set[ix]->GetYaxis()->LabelsOption("h");
              
              for (int ij=0; ij<2*nlayer; ij++) {
                strp_count_set[ix]->GetXaxis()->SetBinLabel(ij+1, xlabels[ij]);
              }
              strp_count_set[ix]->GetXaxis()->LabelsOption("v");
              strp_count_set[ix]->Draw("colz:text25");
              c44->Update();
            }
            c44->Clear();
          }
#endif
          ps.NewPage();
          gStyle->SetOptTitle(1);
          gPad->SetLeftMargin(0.05);
          gStyle->SetPadRightMargin(0.12);
          gStyle->SetPadTopMargin(0.10);
          gStyle->SetPadBottomMargin(0.12);
          
          gStyle->SetOptLogy(0);
          gStyle->SetOptLogz(0);
          //	  gStyle->SetOptStat(1110);
          //	  gStyle->SetOptFit(100);
          gStyle->SetStatW(.36); //40);
          gStyle->SetStatH(.28); //30);
          gStyle->SetTitleFontSize(0.07);
          
          TCanvas* c4a = new TCanvas("c4a", "c4a", 700, 900);
          c4a->Divide(3,4);
#ifndef MONTECARLO
          if (plot_level>90) { 
            for (int iyy=0; iyy<4; iyy++) {
              if (!isTiming && iyy>=2) continue;
              for (int ix=0; ix<=iset; ix++) {
                TH1F* histxx[nlayer]={0};
                ps.NewPage();
                for (int ij=0; ij<nlayer; ij++) {
                  c4a->cd(ij+1);
                  switch(iyy) { 
                  case 0 : histxx[ij] = (TH1F*)xlayer_reso_set[ij][ix]->Clone(); break;
                  case 1 : histxx[ij] = (TH1F*)ylayer_reso_set[ij][ix]->Clone(); break;
                  case 2 : histxx[ij] = (TH1F*)time_xreso_set[ij][ix]->Clone(); break;
                  case 3 : histxx[ij] = (TH1F*)time_yreso_set[ij][ix]->Clone(); break;
                  default : histxx[ij] = (TH1F*)xlayer_reso_set[ij][ix]->Clone(); break;
                  }
                  histxx[ij]->GetXaxis()->SetLabelSize(0.07);
                  histxx[ij]->GetXaxis()->SetTitle(histxx[ij]->GetName());
                  histxx[ij]->GetXaxis()->CenterTitle();
                  histxx[ij]->GetXaxis()->SetTitleSize(0.07);
                  histxx[ij]->GetXaxis()->SetTitleOffset(.9);
                  
                  histxx[ij]->GetYaxis()->SetLabelSize(0.07);
                  TFitResultPtr ptr = histxx[ij]->Fit("gaus", "SQ");
                  Int_t fitStatus = ptr;
                  if (fitStatus==0) { 
                    latex.DrawLatex(0.60, 0.84,Form("#scale[0.5]{%g/%i}", int(100*ptr->Chi2())/100., ptr->Ndf()));
                    latex.DrawLatex(0.72, 0.76,Form("#scale[0.5]{%g}", int(1000*ptr->Parameter(1))/1000.));
                    latex.DrawLatex(0.72, 0.68,Form("#scale[0.5]{%g}", int(1000*ptr->Parameter(2))/1000.));
                  }
                  latex.DrawLatex(0.72, 0.60, Form("#scale[0.5]{%g}", histxx[ij]->Integral()));
                  latex.DrawLatex(0.72, 0.52,Form("#scale[0.5]{%g}", int(1000*histxx[ij]->GetMean())/1000.));
                  latex.DrawLatex(0.72, 0.44,Form("#scale[0.5]{%g}", int(1000*histxx[ij]->GetRMS())/1000.));
                  
                  
                }
                c4a->Update();
                for (int ij=0; ij<nlayer; ij++) {
                  if (histxx[ij]) { delete histxx[ij]; histxx[ij]=0;}
                }
              }
            }
          }
#endif
          gStyle->SetOptStat(0);
          gStyle->SetOptFit(0);
          
          for (int ij=0; ij<nlayer; ij++) {
            c4a->cd(ij+1);
            raw_occu[ij]->Scale(100./ntotal);
            raw_occu[ij]->GetXaxis()->SetLabelSize(.07);
            raw_occu[ij]->GetXaxis()->SetTitle("X-strip");
            raw_occu[ij]->GetXaxis()->CenterTitle();
            raw_occu[ij]->GetXaxis()->SetTitleSize(0.07);
            raw_occu[ij]->GetXaxis()->SetTitleOffset(.8);
            
            raw_occu[ij]->GetYaxis()->SetTitle("Y-strip");
            raw_occu[ij]->GetYaxis()->CenterTitle();
            raw_occu[ij]->GetYaxis()->SetTitleSize(0.07);
            raw_occu[ij]->GetYaxis()->SetTitleOffset(.78);
            raw_occu[ij]->GetYaxis()->SetLabelSize(.07);
            raw_occu[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
            raw_occu[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
            //	    raw_occu[ij]->GetZaxis()->SetLabelSize(.07);
            raw_occu[ij]->Draw("colz");
          }
          c4a->Update();
          
          ps.NewPage();
          gStyle->SetOptTitle(1);
          gPad->SetLeftMargin(0.05);
          gStyle->SetPadRightMargin(0.12);
          gStyle->SetPadTopMargin(0.10);
          gStyle->SetPadBottomMargin(0.12);
          
          gStyle->SetOptLogy(0);
          gStyle->SetOptLogz(0);
          //	  gStyle->SetOptStat(1110);
          //	  gStyle->SetOptFit(100);
          gStyle->SetStatW(.36); //40);
          gStyle->SetStatH(.28); //30);
          gStyle->SetTitleFontSize(0.07);
          gStyle->SetOptStat(0);
          gStyle->SetOptFit(0);
          for (int ij=0; ij<nlayer; ij++) {
            c4a->cd(ij+1);
            raw_seloccu[ij]->Scale(100./nseltot);
            raw_seloccu[ij]->GetXaxis()->SetLabelSize(.07);
            raw_seloccu[ij]->GetXaxis()->SetTitle("X-strip");
            raw_seloccu[ij]->GetXaxis()->CenterTitle();
            raw_seloccu[ij]->GetXaxis()->SetTitleSize(0.07);
            raw_seloccu[ij]->GetXaxis()->SetTitleOffset(.8);
            
            raw_seloccu[ij]->GetYaxis()->SetTitle("Y-strip");
            raw_seloccu[ij]->GetYaxis()->CenterTitle();
            raw_seloccu[ij]->GetYaxis()->SetTitleSize(0.07);
            raw_seloccu[ij]->GetYaxis()->SetTitleOffset(.78);
            raw_seloccu[ij]->GetYaxis()->SetLabelSize(.07);
            raw_seloccu[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
            raw_seloccu[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
            //	    raw_occu[ij]->GetZaxis()->SetLabelSize(.07);
            raw_seloccu[ij]->Draw("colz");
          }
          c4a->Update();
          
          ps.NewPage();
          gStyle->SetOptTitle(1);
          gPad->SetLeftMargin(0.05);
          gStyle->SetPadRightMargin(0.12);
          gStyle->SetPadTopMargin(0.10);
          gStyle->SetPadBottomMargin(0.12);
          
          gStyle->SetOptLogy(0);
          gStyle->SetOptLogz(0);
          gStyle->SetStatW(.36); //40);
          gStyle->SetStatH(.28); //30);
          gStyle->SetTitleFontSize(0.07);
          gStyle->SetOptStat(0);
          gStyle->SetOptFit(0);
          for (int ij=0; ij<nlayer; ij++) {
            c4a->cd(ij+1);
            raw_noiseoccu[ij]->Scale(100./nseltot);
            raw_noiseoccu[ij]->GetXaxis()->SetLabelSize(.07);
            raw_noiseoccu[ij]->GetXaxis()->SetTitle("X-strip");
            raw_noiseoccu[ij]->GetXaxis()->CenterTitle();
            raw_noiseoccu[ij]->GetXaxis()->SetTitleSize(0.07);
            raw_noiseoccu[ij]->GetXaxis()->SetTitleOffset(.8);
            
            raw_noiseoccu[ij]->GetYaxis()->SetTitle("Y-strip");
            raw_noiseoccu[ij]->GetYaxis()->CenterTitle();
            raw_noiseoccu[ij]->GetYaxis()->SetTitleSize(0.07);
            raw_noiseoccu[ij]->GetYaxis()->SetTitleOffset(.78);
            raw_noiseoccu[ij]->GetYaxis()->SetLabelSize(.07);
            raw_noiseoccu[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
            raw_noiseoccu[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
            //	    raw_occu[ij]->GetZaxis()->SetLabelSize(.07);
            raw_noiseoccu[ij]->Draw("colz");
          }
          c4a->Update();
          
          if (plot_level>90) { 
            for (int ix=0; ix<iset; ix++) {
              ps.NewPage();
              for (int ij=0; ij<nlayer; ij++) {
                c4a->cd(ij+1);
                raw_occu_set[ij][ix]->Scale(100./nCount);
                raw_occu_set[ij][ix]->GetXaxis()->SetLabelSize(.07);
                raw_occu_set[ij][ix]->GetXaxis()->SetTitle("X-strip");
                raw_occu_set[ij][ix]->GetXaxis()->CenterTitle();
                raw_occu_set[ij][ix]->GetXaxis()->SetTitleSize(0.07);
                raw_occu_set[ij][ix]->GetXaxis()->SetTitleOffset(.78);
                
                raw_occu_set[ij][ix]->GetYaxis()->SetTitle("Y-strip");
                raw_occu_set[ij][ix]->GetYaxis()->CenterTitle();
                raw_occu_set[ij][ix]->GetYaxis()->SetTitleSize(0.07);
                raw_occu_set[ij][ix]->GetYaxis()->SetTitleOffset(.7);
                raw_occu_set[ij][ix]->GetYaxis()->SetLabelSize(.07);
                
                raw_occu_set[ij][ix]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
                raw_occu_set[ij][ix]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
                raw_occu_set[ij][ix]->Draw("colz");
              }
              c4a->Update();
            }
          }
          ps.NewPage();
          for (int ij=0; ij<nlayer; ij++) {
            c4a->cd(ij+1);
            rawhits_corr_xymul[ij]->Scale(100./ntotal);
            rawhits_corr_xymul[ij]->GetXaxis()->SetLabelSize(.07);
            rawhits_corr_xymul[ij]->GetYaxis()->SetLabelSize(.07);
            rawhits_corr_xymul[ij]->GetXaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
            rawhits_corr_xymul[ij]->GetYaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
            rawhits_corr_xymul[ij]->Draw("colz");
          }
          c4a->Update();
          if(plot_level>90) {
            int icol=0;
            for (int ij=0; ij<nlayer; ij++) {
              for (int jk=ij+1; jk<nlayer; jk++) {
                icol++;
                if (icol==1) {ps.NewPage();}
                c4a->cd(2*icol-1);
                rawhits_xlay_corr_mul[ij][jk]->Scale(100./ntotal);
                rawhits_xlay_corr_mul[ij][jk]->GetXaxis()->SetLabelSize(.07);
                rawhits_xlay_corr_mul[ij][jk]->GetYaxis()->SetLabelSize(.07);
                rawhits_xlay_corr_mul[ij][jk]->GetXaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
                rawhits_xlay_corr_mul[ij][jk]->GetYaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
                rawhits_xlay_corr_mul[ij][jk]->Draw("colz");
                
                c4a->cd(2*icol);
                rawhits_ylay_corr_mul[ij][jk]->Scale(100./ntotal);
                rawhits_ylay_corr_mul[ij][jk]->GetXaxis()->SetLabelSize(.07);
                rawhits_ylay_corr_mul[ij][jk]->GetYaxis()->SetLabelSize(.07);
                rawhits_ylay_corr_mul[ij][jk]->GetXaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
                rawhits_ylay_corr_mul[ij][jk]->GetYaxis()->SetRangeUser(-0.5, 5.0); //(nusedstrip);
                rawhits_ylay_corr_mul[ij][jk]->Draw("colz");
                if (icol==6) {c4a->Update(); icol=0;}
              }
            }
            if (icol!=0) { c4a->Update();}
          }
#ifdef TIMESLOPE	  
          //#ifndef MONTECARLO
          if (isTiming) {
            
            ps.NewPage();
            gStyle->SetOptTitle(1);
            gStyle->SetTitleOffset(-0.06,"XYZ");
            gStyle->SetTitleFillColor(10);
            gStyle->SetTitleFontSize(0.06);
            gStyle->SetPadTopMargin(0.08);
            gStyle->SetPadRightMargin(0.15);
            gStyle->SetPadLeftMargin(0.12);
            gStyle->SetPadBottomMargin(0.11);
            
            TCanvas* c3a = new TCanvas("c3a", "c3a", 500, 700);
            c3a->Divide(3,4);
            for (int itdc=0; itdc<nTDCpLayer; itdc++) { 
              if (itdc>0) {ps.NewPage();}
              for (int il=0; il<nlayer; il++) {
                c3a->cd(il+1); 
                time_xslope[il][itdc]->GetXaxis()->SetLabelSize(.07);
                time_xslope[il][itdc]->GetYaxis()->SetLabelSize(.07);
                time_xslope[il][itdc]->GetXaxis()->SetTitle("Strip No");
                time_xslope[il][itdc]->GetXaxis()->SetTitleOffset(.85);
                time_xslope[il][itdc]->GetXaxis()->SetTitleSize(.06);
                time_xslope[il][itdc]->GetXaxis()->CenterTitle();
                time_xslope[il][itdc]->GetYaxis()->SetTitle("Measured Time (ns)");
                time_xslope[il][itdc]->GetYaxis()->SetTitleOffset(.8);
                time_xslope[il][itdc]->GetYaxis()->SetTitleSize(.06);
                time_xslope[il][itdc]->GetYaxis()->CenterTitle();
                
                time_xslope[il][itdc]->Draw("colz");
                time_xslope[il][itdc]->ProfileX()->Draw("same");
              }
              c3a->Update();
              
              ps.NewPage();
              for (int il=0; il<nlayer; il++) {
                c3a->cd(il+1); 
                time_yslope[il][itdc]->GetXaxis()->SetLabelSize(.07);
                time_yslope[il][itdc]->GetYaxis()->SetLabelSize(.07);
                time_yslope[il][itdc]->GetXaxis()->SetTitle("Strip No");
                time_yslope[il][itdc]->GetXaxis()->SetTitleOffset(.6);
                time_yslope[il][itdc]->GetXaxis()->SetTitleSize(.06);
                time_yslope[il][itdc]->GetXaxis()->CenterTitle();
                time_yslope[il][itdc]->GetYaxis()->SetTitle("Occupancy (%)");
                time_yslope[il][itdc]->GetYaxis()->SetTitleOffset(.6);
                time_yslope[il][itdc]->GetYaxis()->SetTitleSize(.06);
                time_yslope[il][itdc]->Draw("colz");
                time_yslope[il][itdc]->ProfileX()->Draw("same");
              }
              c3a->Update();
            }
            delete c3a;
            
            gStyle->SetOptTitle(1);
            gStyle->SetTitleOffset(-0.06,"XYZ");
            gStyle->SetTitleFillColor(10);
            gStyle->SetTitleFontSize(0.06);
            gStyle->SetPadTopMargin(0.08);
            gStyle->SetPadRightMargin(0.15);
            gStyle->SetPadLeftMargin(0.12);
            gStyle->SetPadBottomMargin(0.09);
            
            ps.NewPage();
            TCanvas* c7 = new TCanvas("c7", "c7", 700, 900);
            c7->Divide(1,3);
            TH2F* tmp2d = (TH2F*)time_layerstrip->Clone("hits2d");
            c7->cd(1); tmp2d->Draw("colz");
            time_layerstrip->GetYaxis()->SetRangeUser(400., 1800.0);
            c7->cd(2); time_layerstrip->Draw("colz");
            TH1F* tmp1d = (TH1F*)tmp2d->ProjectionX("zeroTime", 0, 1);
            TH1F* tmp1dall= (TH1F*)tmp2d->ProjectionX("allTime", 0, tmp2d->GetNbinsY());
            tmp1d->Divide(tmp1dall);
            //	    tmp1d->Scale(1./max(1., time_layerstrip->GetBinContent(0,0)));
            tmp1d->GetYaxis()->SetTitle("Fraction of missign time");
            c7->cd(3); tmp1d->Draw();
            
            c7->Update();
            
            for (int ij=0; ij<nlayer; ij++) {
              for (int jk=0; jk<nstrip; jk++) {
                for (int kl=0; kl<nstrip; kl++) {
                  for (int lm=0; lm<ntimecor; lm++) {
                    double ent = indtimexy_correl[ij][jk][kl][lm]->GetEntries();
                    if (ent >10) {
                      double mean = indtimexy_correl[ij][jk][kl][lm]->GetMean();
                      double rms = indtimexy_correl[ij][jk][kl][lm]->GetRMS();
                      
                      indtimexy_cormean[ij][lm]->Fill(jk, kl, mean+20);
                      indtimexy_corrms[ij][lm]->Fill(jk, kl, rms);
                      
                      if (indtimexy_correl[ij][jk][kl][lm]->GetEntries()>3) {
                        TFitResultPtr ptr = indtimexy_correl[ij][jk][kl][lm]->Fit("gaus","SQ0"); //,"",istr, iend);
                        Int_t fitStatus = ptr;
                        if (fitStatus==0) { 
                          mean = ptr->Parameter(1);
                          rms = ptr->Parameter(2);
                          indtimexy_fitmean[ij][lm]->Fill(jk, kl, mean+20);
                          indtimexy_fitrms[ij][lm]->Fill(jk, kl, rms);
                        }
                      }
                    }
                  }
                }
              }
            }
            
            //	    ps.NewPage();
            gStyle->SetTitleFontSize(0.07);
            for (int ixx=0; ixx<4; ixx++) {
              for (int lm=0; lm<ntimecor; lm++) {
                TH2F* histxx[nlayer]={0};
                ps.NewPage();
                for (int ij=0; ij<nlayer; ij++) {
                  c4a->cd(ij+1);
                  switch(ixx) { 
                  case 0 : histxx[ij] = (TH2F*)indtimexy_cormean[ij][lm]->Clone(); break;
                  case 1 : histxx[ij] = (TH2F*)indtimexy_corrms[ij][lm]->Clone(); break;
                  case 2 : histxx[ij] = (TH2F*)indtimexy_fitmean[ij][lm]->Clone(); break;
                  case 3 : histxx[ij] = (TH2F*)indtimexy_fitrms[ij][lm]->Clone(); break;
                  default : histxx[ij] = (TH2F*)indtimexy_cormean[ij][lm]->Clone(); break;
                  }
                  
                  histxx[ij]->GetXaxis()->SetLabelSize(.07);
                  histxx[ij]->GetXaxis()->SetTitle("X-strip");
                  histxx[ij]->GetXaxis()->CenterTitle();
                  histxx[ij]->GetXaxis()->SetTitleSize(0.07);
                  histxx[ij]->GetXaxis()->SetTitleOffset(.9);
                  
                  histxx[ij]->GetYaxis()->SetTitle("Y-strip");
                  histxx[ij]->GetYaxis()->CenterTitle();
                  histxx[ij]->GetYaxis()->SetTitleSize(0.07);
                  histxx[ij]->GetYaxis()->SetTitleOffset(.9);
                  histxx[ij]->GetYaxis()->SetLabelSize(.07);
                  
                  double amn = histxx[ij]->GetMinimum();
                  double amx = histxx[ij]->GetMaximum();
                  if (amx >amn+0.5) {
                    switch(ixx) {
                    case 0 : ;
                    case 2 : 
                      histxx[ij]->SetMaximum(min(30.0, amx));
                      histxx[ij]->SetMinimum(max(10.0, amn)); break;
                    case 1 : 
                      histxx[ij]->SetMaximum(min(3.0, amx));
                      histxx[ij]->SetMinimum(max(0.2, amn)); break;
                    case 3 : 
                      histxx[ij]->SetMaximum(min(2.0, amx));
                      histxx[ij]->SetMinimum(max(0.2, amn)); break; 
                    default : break;
                    }
                  }
                  histxx[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
                  histxx[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
                  histxx[ij]->Draw("colz");
                }
                
                c4a->Update();
                for (int ij=0; ij<nlayer; ij++) {
                  if (histxx[ij]) { delete histxx[ij]; histxx[ij]=0;}
                }
              }
            }
            
            ps.NewPage();
            gStyle->SetOptFit(101);
            gStyle->SetStatW(.30);
            gStyle->SetStatH(.20);
            for (int lm=0; lm<ntimecor; lm++) {
              ps.NewPage();
              
              for (int ij=0; ij<nlayer; ij++) {
                c4a->cd(ij+1);
                int istr=1;
                int iend =63;
                int nbin = indtimexy_prof[ij][lm]->GetNbinsX();
                for (int ix=nbin/2; ix>=0; ix--) {
                  if (indtimexy_prof[ij][lm]->GetBinError(ix+1) <1.e-4) {
                    istr = ix+1; break;
                  }
                }
                
                for (int ix=nbin/2; ix<nbin; ix++) {
                  if (indtimexy_prof[ij][lm]->GetBinError(ix+1) <1.e-4) {
                    iend = ix-1; break;
                  }
                }
                
                if (istr<nbin/2 && iend>nbin/2 ) { 
                  indtimexy_prof[ij][lm]->GetXaxis()->SetLabelSize(.07);
                  indtimexy_prof[ij][lm]->GetYaxis()->SetLabelSize(.07);
                  /*TFitResultPtr ptr = */ indtimexy_prof[ij][lm]->Fit("pol2","SQ","",istr, iend);
                }
              }
              c4a->Update();
            }
          } // if (isTiming)
          //#endif //MONTECARLO
#endif //ISTIMING
        } // if (filloccu)
        
        if (occulyr==nlayer) {
          totxtime /=max(1,ntotxtime);
          totytime /=max(1,ntotytime);
          ytimeshift += totytime - totxtime;
          file_out<<"ytimeshift "<< iiter <<" "<< ytimeshift<<" "<<totytime<<" "<<totxtime<<endl;
        }
        
        double ncx = 0;
        double ncy = 0;
        
        if (isTiming) {
          if (firstiter==0 && ntcor==0) {    
            
            firstiter=1;
            
            for (int ij=1; ij<nlayer; ij++) {
              ps.NewPage();
              gStyle->SetOptTitle(1);
              gStyle->SetOptStat(111110);
              gStyle->SetOptFit(101);
              gStyle->SetOptLogy(1);
              gStyle->SetTitleFontSize(0.06);
              gStyle->SetPadBottomMargin(0.10);
              gStyle->SetPadTopMargin(0.08);
              gStyle->SetPadLeftMargin(0.10);
              gStyle->SetPadRightMargin(0.02);
              gStyle->SetStatW(.20);
              gStyle->SetStatH(.14);
              
              TCanvas *c3=new TCanvas ("c3","Time Residual",500,700);
              c3->Divide(2,2);
              TF1* fittx[2]={0};
              TF1* fitty[2]={0};
              
              for (int ix=0; ix<2; ix++) { 
                
                c3->cd(2*ix+1);
                int iyy=(ix==0) ? 0 : iiter+1;
                sprintf(name, "fittx_%i", ix);
                fittx[ix] = new TF1(name, gausX, timex_shift[ij][iyy]->GetXaxis()->GetXmin(), timex_shift[ij][iyy]->GetXaxis()->GetXmax(), 3);
                double  parx[3]={timex_shift[ij][iyy]->GetMaximum(), timex_shift[ij][iyy]->GetMean(), timex_shift[ij][iyy]->GetRMS()};
                fittx[ix]->SetParameters(parx);
                
                timex_shift[ij][iyy]->GetXaxis()->SetTitle("#Deltat (ns)");
                timex_shift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
                timex_shift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
                timex_shift[ij][iyy]->GetXaxis()->CenterTitle();
                timex_shift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
                timex_shift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);
                
                timex_shift[ij][iyy]->Fit(fittx[ix],"Q");
                
                if (isTimeCorOrReso && ix>0 && isalign >0 && iiter<nmxiter-1 && 
                    timex_shift[ij][iyy]->Integral() >nmnentry && 
                    abs(timex_shift[ij][iyy]->GetMean() - fittx[ix]->GetParameter(1)) < 
                    timex_shift[ij][iyy]->GetRMS()) {
                  timeoffsetx[ij] +=timex_shift[ij][iyy]->GetMean(); //fittx[ix]->GetParameter(1);
                }
                
                c3->cd(2*ix+2);
                sprintf(name, "fitty_%i", ix);
                fitty[ix] = new TF1(name, gausX, timey_shift[ij][iyy]->GetXaxis()->GetXmin(), timey_shift[ij][iyy]->GetXaxis()->GetXmax(), 3);
                double pary[3]={timey_shift[ij][iyy]->GetMaximum(), timey_shift[ij][iyy]->GetMean(), timey_shift[ij][iyy]->GetRMS()};
                fitty[ix]->SetParameters(pary);
                timey_shift[ij][iyy]->GetXaxis()->SetTitle("#Deltat (ns)");
                timey_shift[ij][iyy]->GetXaxis()->CenterTitle();
                timey_shift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
                timey_shift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
                timey_shift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
                timey_shift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);
                timey_shift[ij][iyy]->Fit(fitty[ix],"Q");
                
                if (ix>0) { 
                  if (isTimeCorOrReso && isalign>0 && iiter<nmxiter-1 && 
                      timey_shift[ij][iyy]->Integral() >nmnentry && 
                      abs(timey_shift[ij][iyy]->GetMean() - fitty[ix]->GetParameter(1)) < 
                      timey_shift[ij][iyy]->GetRMS()) {
                    timeoffsety[ij] +=timey_shift[ij][iyy]->GetMean(); //fitty[ix]->GetParameter(1);
                    //		    file_out<<"timeoffsety["<<ij<<"]= "<<ntcor<<" "<<timeoffsety[ij]<<" "<<timex_shift[ij][iyy]->GetMean()<<endl;
                  }
                  
                  file_out <<"global time "<<iiter<<" "<<ij<<" "<< parx[0]<<" "<<parx[1]<<" "<<timeoffsetx[ij]<<" "<<parx[2]<< " y "
                           << pary[0]<<" "<<pary[1]<<" "<<timeoffsety[ij]<<" "<<pary[2]<<endl;
                  
                  file_out<<"fit "<<iiter<<" "<<ij<<" "<<fittx[ix]->GetParameter(0)<<" "<<fittx[ix]->GetParameter(1)<<" "<<fittx[ix]->GetParameter(2)<<" "<<fitty[ix]->GetParameter(0)<<" "<<fitty[ix]->GetParameter(1)<<" "<<fitty[ix]->GetParameter(2)<<endl;
                }
              } //  for (int ix=0; ix<2; ix++)
              c3->Update();
              if (c3) { delete c3; c3=0;}
              for (int ix=0; ix<2; ix++) {
                if (fittx[ix]) { delete fittx[ix]; fittx[ix]=0;}
                if (fitty[ix]) { delete fitty[ix]; fitty[ix]=0;}
              }
            } //for (int ij=1; ij<nlayer; ij++)
            
            for (int ij=1; ij<nlayer; ij++) {
              ps.NewPage();
              
              gStyle->SetOptStat(0);
              gStyle->SetOptFit(0);
              gStyle->SetOptLogy(0);
              gStyle->SetTitleFontSize(0.06);
              gStyle->SetPadBottomMargin(0.10);
              gStyle->SetPadTopMargin(0.08);
              gStyle->SetPadLeftMargin(0.08);
              gStyle->SetPadRightMargin(0.12);
              
              TCanvas *c3=new TCanvas ("c3","Time Residual",500,700);
              c3->Divide(2,2);
              
              for (int ix=0; ix<2; ix++) { 
                c3->cd(2*ix+1);
                int iyy=(ix==0) ? 0 : iiter+1;
                timex_2dshift[ij][iyy]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
                timex_2dshift[ij][iyy]->GetXaxis()->SetTitle("Strip No");
                timex_2dshift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
                timex_2dshift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
                timex_2dshift[ij][iyy]->GetXaxis()->CenterTitle();
                timex_2dshift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
                timex_2dshift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);
                
                timex_2dshift[ij][iyy]->GetYaxis()->SetTitle("#Deltat (ns)");
                timex_2dshift[ij][iyy]->GetYaxis()->SetTitleOffset(.7);
                timex_2dshift[ij][iyy]->GetYaxis()->SetTitleSize(.06);
                timex_2dshift[ij][iyy]->GetYaxis()->CenterTitle();
                timex_2dshift[ij][iyy]->GetYaxis()->SetLabelSize(.05);
                timex_2dshift[ij][iyy]->GetZaxis()->SetLabelSize(.04);
                timex_2dshift[ij][iyy]->Draw("colz");
                
                c3->cd(2*ix+2);
                timey_2dshift[ij][iyy]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
                timey_2dshift[ij][iyy]->GetXaxis()->SetTitle("Strip No");
                timey_2dshift[ij][iyy]->GetXaxis()->CenterTitle();
                timey_2dshift[ij][iyy]->GetXaxis()->SetTitleOffset(.7);
                timey_2dshift[ij][iyy]->GetXaxis()->SetTitleSize(.06);
                timey_2dshift[ij][iyy]->GetXaxis()->SetLabelSize(.05);
                timey_2dshift[ij][iyy]->GetXaxis()->SetLabelOffset(-0.001);
                
                timey_2dshift[ij][iyy]->GetYaxis()->SetTitle("#Deltat (ns)");
                timey_2dshift[ij][iyy]->GetYaxis()->CenterTitle();
                timey_2dshift[ij][iyy]->GetYaxis()->SetTitleOffset(.7);
                timey_2dshift[ij][iyy]->GetYaxis()->SetTitleSize(.06);
                timey_2dshift[ij][iyy]->GetYaxis()->SetLabelSize(.05);
                timey_2dshift[ij][iyy]->GetZaxis()->SetLabelSize(.04);
                timey_2dshift[ij][iyy]->Draw("colz");
                
              } //  for (int ix=0; ix<2; ix++)
              c3->Update();
              if (c3) { delete c3; c3=0;}
            } //for (int ij=1; ij<nlayer; ij++)
          } // if (firstiter==0 && ntcor==0) 
          //	  double ncx = 0;
          //	  double ncy = 0;
          
          for (int ij=0; ij<dir_cx[occulyr][iiter]->GetNbinsX(); ij++) {
            if (dir_cx[occulyr][iiter]->GetBinCenter(ij)<0.0) {
              ncx +=dir_cx[occulyr][iiter]->GetBinContent(ij);
              ncy +=dir_cy[occulyr][iiter]->GetBinContent(ij);
            } else { break;}
          }
          file_out <<" occulyr "<< occulyr<<" "<< iiter <<" xc "<<dir_cx[occulyr][iiter]->Integral()<<" "<<dir_cx[occulyr][iiter]->GetBinContent(0)<<" "<<dir_cx[occulyr][iiter]->GetEntries()<<" "<<dir_cx[occulyr][iiter]->GetMean()<<" "<<dir_cx[occulyr][iiter]->GetRMS()<<" yc "<<dir_cy[occulyr][iiter]->Integral()<<" "<<dir_cy[occulyr][iiter]->GetBinContent(0)<<" "<<dir_cy[occulyr][iiter]->GetEntries()<<" "<<dir_cy[occulyr][iiter]->GetMean()<<" "<<dir_cy[occulyr][iiter]->GetRMS()<<" ncx "<<ncx<<" ncy "<<ncy<<endl; 
          
        } //if (isTiming)
        
        
        ps.NewPage();
        gStyle->SetOptStat(111110);
        gStyle->SetOptFit(101);
        gStyle->SetOptLogy(1);
        
        gStyle->SetPadBottomMargin(0.11);
        gStyle->SetPadTopMargin(0.08);
        gStyle->SetPadLeftMargin(0.07);
        gStyle->SetPadRightMargin(0.02);
        gStyle->SetOptTitle(1);
        gStyle->SetTitleFontSize(0.07);
        gStyle->SetStatW(.30);
        gStyle->SetStatH(.14);
        gStyle->SetStatY(.99);
        gStyle->SetStatX(.99);
        TCanvas *c2a=new TCanvas ("c2a","Slope",500,700);
        if (isTiming) {c2a->Divide(2,3);} else {c2a->Divide(2,2);}
        const int nxplot=6;
        //	latex.SetTextSize(0.08);
        TH1F* histax[nxplot]={0};
        for (int ix=0; ix<nxplot; ix++) {
          if (!isTiming && ix>=4) continue;
          switch(ix) {
          case 0 : histax[ix] = (TH1F*)pos_xslope[occulyr][iiter]->Clone(); break;
          case 1 : histax[ix] = (TH1F*)pos_yslope[occulyr][iiter]->Clone(); break;
          case 2 : histax[ix] = (TH1F*)pos_theta[occulyr][iiter]->Clone(); break;
          case 3 : histax[ix] = (TH1F*)pos_phi[occulyr][iiter]->Clone(); break;
          case 4 : histax[ix] = (TH1F*)dir_cx[occulyr][iiter]->Clone(); break;
          case 5 : histax[ix] = (TH1F*)dir_cy[occulyr][iiter]->Clone(); break;
          default : histax[ix] = (TH1F*)pos_xslope[occulyr][iiter]->Clone(); break;
          }
          c2a->cd(ix+1);
          histax[ix]->GetXaxis()->CenterTitle();
          histax[ix]->GetXaxis()->SetTitleOffset(0.7);
          histax[ix]->GetXaxis()->SetTitleSize(0.06);
          histax[ix]->GetXaxis()->SetLabelOffset(-0.01);
          histax[ix]->GetXaxis()->SetLabelSize(0.06);
          histax[ix]->Draw();
          //	  if (ix==4) {latex.DrawLatex(0.16, 0.66,Form("#font[12]#color[2]{%g}", int(100000*ncx/max(1.,dir_cx[occulyr][iiter]->GetEntries()))/1000.));}
          //	  if (ix==5) {latex.DrawLatex(0.16, 0.66,Form("#scale[1.2]#color[2]{%g}", int(100000*ncy/max(1.,dir_cy[occulyr][iiter]->GetEntries()))/1000.));}
          if (ix==4) {latex.DrawLatex(0.12, 0.66,Form("#scale[0.6]{%g\%}", int(1000000*ncx/max(1.,dir_cx[occulyr][iiter]->GetEntries()))/10000.));}
          if (ix==5) {latex.DrawLatex(0.12, 0.66,Form("#scale[0.6]{%g\%}", int(1000000*ncy/max(1.,dir_cy[occulyr][iiter]->GetEntries()))/10000.));}
          
        }
        c2a->Update();
        if (c2a) { delete c2a; c2a=0;}
        for (int ix=0; ix<nxplot; ix++) {
          if (histax[ix]) { delete histax[ix]; histax[ix]=0;}
        }
        if (isTiming) { 
          ps.NewPage();
          TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
          c4c->Divide(3,4);
          TH1F* histbx[2*nprob]={0};
          for (int ix=0; ix<2*nprob; ix++) {
            if (ix<nprob) { 
              histbx[ix] = (TH1F*)dir_cxlay[occulyr][iiter][ix]->Clone();
            } else {
              histbx[ix] = (TH1F*)dir_cylay[occulyr][iiter][ix-nprob]->Clone();
            }
            c4c->cd(ix+1);
            int ncx=0;
            for (int ij=0; ij<histbx[ix]->GetNbinsX(); ij++) {
              if (histbx[ix]->GetBinCenter(ij)<0.0) {
                ncx +=histbx[ix]->GetBinContent(ij);
              } else { break;}
            }
            histbx[ix]->GetXaxis()->CenterTitle();
            histbx[ix]->GetXaxis()->SetTitleOffset(0.8);
            histbx[ix]->GetXaxis()->SetTitleSize(0.07);
            histbx[ix]->GetXaxis()->SetLabelOffset(-0.01);
            histbx[ix]->GetXaxis()->SetLabelSize(0.07);
            histbx[ix]->Draw();
            latex.DrawLatex(0.12, 0.66,Form("#scale[0.8]{%g\%}", int(1.e6*ncx/max(1.,histbx[ix]->GetEntries()))/1.e4));
          }
          c4c->Update();
          if (c4c) {delete c4c; c4c=0;}
          for (int ix=0; ix<2*nprob; ix++) {
            if (histbx[ix]) { delete histbx[ix]; histbx[ix]=0;}
          }
        } // if (isTiming) 
        
        int istr = (occulyr >= nlayer) ? 0 : occulyr;
        int iend = (occulyr >= nlayer) ? nlayer : occulyr+1;
        
        for (int iocc = istr; iocc<iend; iocc++) {
          double absWidth=2.5;//6.0;//3.0;//2.5;//min(2.5,max(1.9, 0.15*(nmxiter -iiter)));
          double widthScale=1.0;//2.0;//1.0;//min(0.70,max(0.40, 0.05*(nmxiter -iiter))); //Initially one can allow very large asymmetric distribution and and gradually make it smaller

          ps.NewPage();
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
          
          TCanvas *c2=new TCanvas ("c2","Residue",500,700);
          if (isTiming) { c2->Divide(2,2);} else {c2->Divide(1,2);}
          
          c2->cd(1);
          
          double alw= (ntcor==0) ? -2.05 : -2.8; //xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmin();
          double ahg= (ntcor==0) ?  2.05 :  2.8; //xlayer_reso[iocc][iiterrs]->GetXaxis()->GetXmax();
          TF1* fitfx = new TF1("fitfx", gausX, alw, ahg, 3);
          
          double  parx[3]={xlayer_reso[iocc][iiterrs]->GetMaximum(), xlayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*xlayer_reso[iocc][iiterrs]->GetRMS())};
          fitfx->SetParameters(parx);
          fitfx->SetParLimits(1, parx[1]-1., parx[1]+1.);
          fitfx->SetParLimits(2, 0.12, 1.5*parx[2]);
          
          xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("X-residues (pitch)");
          xlayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
          xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
          xlayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
          
          xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
          xlayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
          
          xlayer_reso[iocc][iiterrs]->Fit(fitfx, "BRQ");
          
          c2->cd(2);
          
          TF1* fitfy = new TF1("fitfy", gausX, alw, ahg, 3);
          double  pary[3]={ylayer_reso[iocc][iiterrs]->GetMaximum(), ylayer_reso[iocc][iiterrs]->GetMean(), max(0.15,0.7*ylayer_reso[iocc][iiterrs]->GetRMS())};
          fitfy->SetParameters(pary);
          fitfy->SetParLimits(1, pary[1]-1., pary[1]+1.);
          
          ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitle("Y-residues (pitch)");
          ylayer_reso[iocc][iiterrs]->GetXaxis()->CenterTitle();
          ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
          ylayer_reso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
          
          ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
          ylayer_reso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
          ylayer_reso[iocc][iiterrs]->Fit(fitfy, "BRQ");
          
          //time resolution only for plot, but not for any calculation
          alw= (ntcor==0) ? -7.0 :-10.0;
          ahg= (ntcor==0) ?  7.0 : 10.0;
          TF1* fittfx = new TF1("fittfx", gausX, alw, ahg, 3);
          TF1* fittfy = new TF1("fittfy", gausX, alw, ahg, 3);
          if (isTiming) { 
            c2->cd(3);
            
            double  partx[3]={time_xreso[iocc][iiterrs]->GetMaximum(), time_xreso[iocc][iiterrs]->GetMean(), max(0.15,0.7*time_xreso[iocc][iiterrs]->GetRMS())};
            fittfx->SetParameters(partx);
            fittfx->SetParLimits(1, partx[1]-1., partx[1]+1.);
            fittfx->SetParLimits(2, 0.6, 1.5*partx[2]);
            
            time_xreso[iocc][iiterrs]->GetXaxis()->SetTitle("X - #Deltat (ns)");
            time_xreso[iocc][iiterrs]->GetXaxis()->CenterTitle();
            
            time_xreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
            time_xreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
            time_xreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
            time_xreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
            
            time_xreso[iocc][iiterrs]->Fit(fittfx, "BRQ");
            time_xrms[iocc] = fittfx->GetParameter(2);
            c2->cd(4);
            
            double  party[3]={time_yreso[iocc][iiterrs]->GetMaximum(), time_yreso[iocc][iiterrs]->GetMean(), max(0.15,0.7*time_yreso[iocc][iiterrs]->GetRMS())};
            fittfy->SetParameters(party);
            fittfy->SetParLimits(1, party[1]-1., party[1]+1.);
            fittfy->SetParLimits(2, 0.6, 1.5*party[2]);
            
            time_yreso[iocc][iiterrs]->GetXaxis()->SetTitle("Y - #Deltat (ns)");
            time_yreso[iocc][iiterrs]->GetXaxis()->CenterTitle();
            time_yreso[iocc][iiterrs]->GetXaxis()->SetTitleOffset(0.7);
            time_yreso[iocc][iiterrs]->GetXaxis()->SetTitleSize(0.06);
            time_yreso[iocc][iiterrs]->GetXaxis()->SetLabelOffset(-0.01);
            time_yreso[iocc][iiterrs]->GetXaxis()->SetLabelSize(0.06);
            
            time_yreso[iocc][iiterrs]->Fit(fittfy, "RQ");
            time_yrms[iocc] = fittfy->GetParameter(2);
            file_out <<"lay "<< iocc<<" it "<< iiterrs <<" X-sh "<< fitfx->GetParameter(1)<<" "<<fitfx->GetParameter(2)<<" Y-sh "<< fitfy->GetParameter(1)<<" "<<fitfy->GetParameter(2)<<" X-sh "<< fittfx->GetParameter(1)<<" "<<fittfx->GetParameter(2)<<" "<<time_xreso[iocc][iiterrs]->GetRMS()<<" Y-sh "<< fittfy->GetParameter(1)<<" "<<fittfy->GetParameter(2)<<" "<<time_yreso[iocc][iiterrs]->GetRMS()<<endl; 
          } //if (isTiming)
          
          c2->Update();
          if (c2) { delete c2; c2=0;}
          
          if (iiter<nmxiter-1 && ntcor==1) {
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
          
          delete fitfx; fitfx=0;
          delete fitfy; fitfy=0;
          delete fittfx; fittfx=0;
          delete fittfy; fittfy=0;
          
          if (ntcor==1 && iiter<nmxiter-1) {//Don't correct for last iteration //150112  
            bias_inpos_xreso2[iocc] = (xlayer_exterr[iocc][iiterrs]->GetMean())*(xlayer_exterr[iocc][iiterrs]->GetMean());
            bias_inpos_yreso2[iocc] = (ylayer_exterr[iocc][iiterrs]->GetMean())*(ylayer_exterr[iocc][iiterrs]->GetMean());
            
            for (int ix=0; ix<nmxhits; ix++) {
              
              TFitResultPtr ptrx = xlayer_reso_mul[iocc][iiterrs][ix]->Fit("gaus","SQ0");
              TFitResultPtr ptry = ylayer_reso_mul[iocc][iiterrs][ix]->Fit("gaus","SQ0");
              Int_t fitStatus = ptrx;
              pos_xrms[iocc][ix] = (fitStatus==0 && xlayer_reso_mul[iocc][iiterrs][ix]->GetEntries()>3) ? ptrx->Parameter(2) : 0.0;
              fitStatus = ptry;
              pos_yrms[iocc][ix] = (fitStatus==0 && ylayer_reso_mul[iocc][iiterrs][ix]->GetEntries()>3) ? ptry->Parameter(2) : 0.0;
              // 18th Feb 2015, for estimation of resolution
              xposerrsq[ix][iocc] =  pos_xrms[iocc][ix]*pos_xrms[iocc][ix] - bias_inpos_xreso2[iocc]; 
              yposerrsq[ix][iocc] =  pos_yrms[iocc][ix]*pos_yrms[iocc][ix] - bias_inpos_yreso2[iocc];
              
              if (xposerrsq[ix][iocc]<0.02) xposerrsq[ix][iocc]=0.02; //resolution should be more than 0.5 ns
              if (yposerrsq[ix][iocc]<0.02) yposerrsq[ix][iocc]=0.02;
            }
          }
          
          if (isTiming) { 
            if (ntcor==1 && iiter<nmxiter-1) { //Don't correct for last iteration
              
              TProfile* str_tdev[2]={0};
              TH1F* str_tdev1d[2]={0};
              for (int ixy=0; ixy<2; ixy++) {
                double meanx[nstrip]={0};
                double weightx[nstrip]={0};
                double sum=0;
                double sumwt=0;
                if (ixy==0) {
                  file_outstr <<"double xtoffystr[nlayer][nstrip] = {"<<endl;	
                } else {
                  file_outstr <<"double ytoffxstr[nlayer][nstrip] = {"<<endl;
                }
                if (ixy==0) {
                  str_tdev[ixy] = (TProfile*)ystr_xtdev[iocc][iiter]->ProfileX("xxx");
                  str_tdev1d[ixy] = (TH1F*)ystr_xtdev[iocc][iiter]->ProjectionX("xxx1d");
                } else {
                  str_tdev[ixy] = (TProfile*)xstr_ytdev[iocc][iiter]->ProfileX("yyy");
                  str_tdev1d[ixy] = (TH1F*)xstr_ytdev[iocc][iiter]->ProjectionX("yyy1d");
                }
                int nbin = min(nstrip,str_tdev[ixy]->GetNbinsX());
                for (int jk=0; jk<nbin; jk++) {
                  meanx[jk] =str_tdev[ixy]->GetBinContent(jk+1);
                  weightx[jk] =str_tdev1d[ixy]->GetBinContent(jk+1);
                  sumwt +=weightx[jk];
                  sum +=meanx[jk]*weightx[jk];
                }
                sum /=max(1.0,sumwt);
                
                file_outstr<<"//toffstr "<<nbin<<" "<<ixy<<" "<<iocc<<" "<<iiter<<" "<<sum<<endl;
                if (isTimeCorOrReso) { 
                  if (ixy==0) {
                    for (int jk=0; jk<nbin; jk++) {
                      xtoffystr[iocc][jk] +=meanx[jk]-sum;
                      file_outstr<<xtoffystr[iocc][jk]<<", ";
                      if ((jk+1)%8==0) file_outstr<<endl;
                    }
                  } else {
                    for (int jk=0; jk<nbin; jk++) {
                      ytoffxstr[iocc][jk] +=meanx[jk]-sum;
                      file_outstr<<ytoffxstr[iocc][jk]<<", ";
                      if ((jk+1)%8==0) file_outstr<<endl;
                    }
                  }
                } //if (isTimeCorOrReso)  
              } // for (int ixy=0; ixy<2; ixy++)
              file_outstr <<"};"<<endl;
              for (int ixy=0; ixy<2; ixy++) {
                if (str_tdev[ixy]) { delete str_tdev[ixy]; str_tdev[ixy]=0;}
                if (str_tdev1d[ixy]) { delete str_tdev1d[ixy]; str_tdev1d[ixy]=0;}
              }
            }

            TH1F* time_shift[nstrip][2]={0};
            double fitmean[nstrip][2]={0};
            double fitrms[nstrip][2]={0};
            double fitchi[nstrip][2]={0};
            double statmean[nstrip][2]={0};
            
            TF1* fity[nstrip][2]={0}; 
            
            const int nsgpr=3;
            double fitres[nsgpr];
            double parerr[nsgpr];
            double fchisq;
            ps.NewPage();
            gStyle->SetOptTitle(0);
            gStyle->SetOptStat(0);
            gStyle->SetOptFit(0);
            gStyle->SetOptLogy(1);
            gStyle->SetPadTopMargin(.001);
            gStyle->SetPadBottomMargin(0.001); //0.07
            gStyle->SetPadLeftMargin(0.001);
            gStyle->SetPadRightMargin(0.001);
            
            TCanvas* c1 = new TCanvas("c1", "c1", 700, 900);
            c1->Divide(8,8,1.e-6, 1.e-6);
            for (int ij=0; ij<2; ij++) {
              for (int jk=0; jk<nstrip; jk++) {
                c1->cd(nstrip*ij+jk+1);
                
                if (ij==0) {
                  time_shift[jk][ij] = (TH1F*)time_xstrreso[iocc][jk][iiterrs]->Clone();
                } else {
                  time_shift[jk][ij] = (TH1F*)time_ystrreso[iocc][jk][iiterrs]->Clone();
                } 
                
                if (time_shift[jk][ij]->Integral()>2) {
                  time_shift[jk][ij]->GetXaxis()->SetLabelSize(.07);
                  
                  double  par[3]={time_shift[jk][ij]->GetMaximum(),time_shift[jk][ij]->GetMean(), time_shift[jk][ij]->GetRMS()};
                  
                  int nbinx = time_shift[jk][ij]->GetNbinsX();
                  int nlow=0;
                  for (int kl=0; kl<nbinx; kl++) {
                    if (time_shift[jk][ij]->GetBinContent(kl+1) >0) {
                      nlow=kl; break;
                    }
                  }
                  int nhig = 0;
                  for (int kl=nbinx; kl>0; kl--) {
                    if (time_shift[jk][ij]->GetBinContent(kl+1) >0) {
                      nhig=kl; break;
                    }
                  }
                  float amean = 0.5*(nlow + nhig);
                  nlow = int(TMath::Max(0., nlow - 1.0*(amean - nlow)));
                  nhig = int(TMath::Min(nbinx-1., nhig + 1.0*(nhig - amean)));
                  
                  nchannel = 0;
                  for (int kl=nlow; kl<=nhig; kl++) {
                    if (nchannel <nmxchn) {
                      m_data[nchannel] = time_shift[jk][ij]->GetBinContent(kl+1);
                      m_xpos[nchannel] = time_shift[jk][ij]->GetBinCenter(kl+1);
                      nchannel++;
                    }
                  }
                  double alw= time_shift[jk][ij]->GetBinCenter(nlow+1);
                  double ahg = time_shift[jk][ij]->GetBinCenter(nhig+1);
                  
                  time_shift[jk][ij]->Draw();
                  
                  TMinuit *gMinuit = new TMinuit(nsgpr);
                  gMinuit->SetPrintLevel(-1);
                  
                  TString hname[nsgpr] = {"height", "mean", "rms"};
                  
                  int nmx = time_shift[jk][ij]->GetMaximumBin();
                  double hgh = 0.35*(time_shift[jk][ij]->GetBinContent(nmx-1) + 
                                     time_shift[jk][ij]->GetBinContent(nmx) + 
                                     time_shift[jk][ij]->GetBinContent(nmx+1));
                  double strt[nsgpr] = {hgh,time_shift[jk][ij]->GetMean(), max(0.6, min(3.0, 0.9*time_shift[jk][ij]->GetRMS()))};
                      
                  double alow[nsgpr] = {0.5*strt[0], strt[1]-2.0, max(0.4*strt[2],0.5)};
                  double ahig[nsgpr] = {2.0*strt[0], strt[1]+2.0, min(1.3*strt[2]+0.1,3.5)};
                  double step[nsgpr] = {0.5, 0.01, 0.01};
                  
                  gMinuit->SetFCN(fcnsg);
                  
                  double arglist[10];
                  int ierflg = 0;
                  arglist[0] =  1 ;
                  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
                  
                  for (int kl=0; kl<nsgpr; kl++) {
                    gMinuit->mnparm(kl, hname[kl], strt[kl], step[kl], alow[kl], ahig[kl],ierflg);
                  }
                  
                  arglist[0] = 0;
                  //	      gMinuit->mnexcm("MIGRAD", arglist, 0, ierflg);
                  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
                  
                  arglist[0] = 0;
                  gMinuit->mnexcm("IMPROVE", arglist, 0, ierflg);
                  
                  TString chnam;
                  double parv,err,xlo,xup, plerr, mierr, eparab, gcc;
                  int iuit;
                  
                  for (int kl=0; kl<nsgpr; kl++) {
                    gMinuit->mnpout(kl, chnam, parv, err, xlo, xup, iuit);
                    gMinuit->mnerrs(kl, plerr, mierr, eparab, gcc);
                    fitres[kl] = parv;
                    parerr[kl] = err;
                  }
                  double  fedm, errdef;
                  int  nparx, istat, fitndof;
                  gMinuit->mnstat(fchisq, fedm, errdef, fitndof, nparx, istat);
                  
                  if (istat==0) fchisq =100000000.0;
                  
                  sprintf(name, "fity_%i_%i", ij, jk);
                  
                  fity[jk][ij] = new TF1(name, gausX, alw, ahg, 3);
                  fity[jk][ij]->SetParameters(fitres);
                  fity[jk][ij]->SetLineColor(2);
                  fity[jk][ij]->SetLineWidth(1);
                  fity[jk][ij]->Draw("same");
                  
                  fitmean[jk][ij] = fitres[1]; // fity[jk][ij]->GetParameter(1);
                  fitrms[jk][ij] = fitres[2];
                  fitchi[jk][ij] = fchisq;
                  
                  statmean[jk][ij] = time_shift[jk][ij]->GetMean();
                  
                  latex.DrawLatex(0.32, 0.36,Form("%g", int(1000*fitres[1])/1000.));
                  latex.DrawLatex(0.32, 0.26,Form("%g", int(1000*statmean[jk][ij])/1000.));
                  latex.DrawLatex(0.32, 0.16,Form("%g", int(1000*fitres[2])/1000.));
                  latex.DrawLatex(0.32, 0.06,Form("%g", int(1000*time_shift[jk][ij]->GetRMS())/1000.));
                  delete gMinuit; gMinuit=0;
                } // if (time_shift[jk][ij]->GetEntries()>15)
              } // for (int jk=0; jk<nstrip; jk++)
            } //for (int ij=0; ij<2; ij++) 
            
            c1->Update();
            if (c1) { delete c1; c1=0;} 
            
            for (int ij=0; ij<2; ij++) {
              for (int jk=0; jk<nstrip; jk++) {
                if (time_shift[jk][ij]) { delete time_shift[jk][ij]; time_shift[jk][ij]=0;}
                if (fity[jk][ij]) { delete fity[jk][ij]; fity[jk][ij]=0;}
              }
            }
            file_out <<"lay3 "<< iiter<<" "<<ntcor<<" "<<iiterrs<<" "<<lstr<<" "<<lend<<" "<<laye<<" "<<occulyr<<" "<<iocc<<" "<<nlayer<<endl;
#ifdef 	TIMESLOPE    
            if (ntcor==1 && iiter<nmxiter-1) { //Don't correct for last iteration
              
              //a	      ps.NewPage();
              //a	      gStyle->SetOptLogy(0);
              //a	      TCanvas* c1x = new TCanvas("c1x", "c1x", 700, 900);
              //a	      c1x->Divide(8,8,1.e-6, 1.e-6);
              
              TH2F* time_slope_pr[2][nstrip]={0};
              TProfile* time_slope_proj[2][nstrip]={0};
              TH1F* time_slope_1d[2][nstrip]={0};
              TF1* time_slope_prfit[2][nstrip]={0};
              
              for (int ixy=0; ixy<2; ixy++) {
                //	      for (int ixy=1; ixy>=0; ixy--) {
                file_outstr<<"//slope_cor "<<ixy<<" "<<iocc<<" "<<iiter<<endl;
                for (int str=0; str<nstrip; str++) {
                  
                  double meanx[nstrip]={0};
                  double weightx[nstrip]={0};
                  double sum=0;
                  double sumwt=0;
                  
                  if (ixy==0) { 
                    time_slope_pr[ixy][str] = (TH2F*) time_xslope_pr[iocc][str][iiter]->Clone();
                  } else {
                    time_slope_pr[ixy][str] = (TH2F*) time_yslope_pr[iocc][str][iiter]->Clone();
                  }
                  time_slope_pr[ixy][str]->GetXaxis()->SetLabelSize(.08);
                  time_slope_pr[ixy][str]->GetYaxis()->SetLabelSize(.08);
                  //a		  c1x->cd(ixy*nstrip+str+1); 
                  //a		  if (iiter==0) { time_slope_pr[ixy][str]->Draw("colz");} //Fit("pol1");
                  sprintf(name, "xxx_%i_%i", ixy, str);
                  time_slope_proj[ixy][str] =(TProfile*)time_slope_pr[ixy][str]->ProfileX(name);
                  time_slope_proj[ixy][str]->SetMarkerSize(0.4);
                  time_slope_proj[ixy][str]->SetMarkerStyle(24);
                  
                  //a		  time_slope_proj[ixy][str]->Draw((iiter==0)?"same":"");
                  sprintf(name, "xxx1d_%i_%i", ixy, str);
                  time_slope_1d[ixy][str] =(TH1F*)time_slope_pr[ixy][str]->ProjectionX(name);
                  
                  int nbin = min(nstrip,time_slope_1d[ixy][str]->GetNbinsX());
                  const int nmnhit=3;
                  for (int jk=0; jk<nbin; jk++) {
                    meanx[jk] =time_slope_proj[ixy][str]->GetBinContent(jk+1);
                    weightx[jk] =time_slope_1d[ixy][str]->GetBinContent(jk+1);
                    if (weightx[jk]>nmnhit) {
                      sumwt +=weightx[jk];
                      sum +=meanx[jk]*weightx[jk];
                    }
                  }
                  sum /=max(1.0,sumwt);
                  file_outstr<<"//slope_corstr "<<ixy<<" "<<nbin<<" "<<str<<" "<<sum<<" "<<statmean[str][ixy]<<" "<<time_slope_pr[ixy][str]->GetEntries()<<endl;
                  if (isTimeCorOrReso && iiter<nmxiter-1) { 
                    if (ixy==0) {
                      for (int jk=0; jk<nbin; jk++) {
                        if (weightx[jk]>nmnhit) {xt_slope_cor[iocc][str][jk] +=meanx[jk]-sum;} //statmean[str][ixy];} //sum;}
                        file_outstr<<xt_slope_cor[iocc][str][jk]<<", ";
                        if ((jk+1)%8==0) file_outstr<<endl;
                      }
                    } else {
                      for (int jk=0; jk<nbin; jk++) {
                        if (weightx[jk]>nmnhit) {yt_slope_cor[iocc][str][jk] +=meanx[jk]-sum;} //statmean[str][ixy];} // sum;}
                        file_outstr<<yt_slope_cor[iocc][str][jk]<<", ";
                        if ((jk+1)%8==0) file_outstr<<endl;
                      }
                    }
                  }
                  
                  double par[nDelayPar]={0,0,0,0,0,0};
                  if (time_slope_pr[ixy][str]->GetEntries()>120) { 
                    int istr=1;
                    int iend =31;
                    int nbin = time_slope_pr[ixy][str]->GetNbinsX();
                    for (int ix=nbin/2; ix>=0; ix--) {
                      if (time_slope_proj[ixy][str]->GetBinError(ix+1) <1.e-4) {
                        istr = ix+1; break;
                      }
                    }
                    
                    for (int ix=nbin/2; ix<nbin; ix++) {
                      if (time_slope_proj[ixy][str]->GetBinError(ix+1) <1.e-4) {
                        iend = ix-1; break;
                      }
                    }
                    
                    if (istr<nbin/2 && iend>nbin/2 ) { 
                      TFitResultPtr ptr = time_slope_proj[ixy][str]->Fit("pol4","SQ0","",istr, iend);
                      Int_t fitStatus = ptr;
                      if (fitStatus==0) {
                        par[0] = ptr->Parameter(0);
                        par[1] = ptr->Parameter(1);
                        par[2] = ptr->Parameter(2);
                        par[3] = ptr->Parameter(3);
                        par[4] = ptr->Parameter(4);
                        
                        file_outstr <<" ixy "<< ixy<<" "<<iocc<<" "<<str<<" "
                                    <<par[0]<<" " <<par[1]<<" " <<par[2]<<" " <<par[3]<<" " <<par[4]<<" "<<
                          par[0]+par[1]*16+par[2]*256+par[3]*4096+par[4]*65536<<" "<<sum<<endl; //central strip
                      }
                    } else {
                      file_outstr <<" ixx "<< ixy<<" "<<iocc<<" "<<str<<" "<<par[0]<<" " <<par[1]<<" " <<par[2]<<" " <<par[3]<<" " <<par[4]<<" "<<par[5]<<" "<<sum<<endl;
                    }
                  } else {
                    file_outstr <<" iyy "<< ixy<<" "<<iocc<<" "<<str<<" "<<par[0]<<" " <<par[1]<<" " <<par[2]<<" " <<par[3]<<" " <<par[4]<<" "<<par[5]<<" "<<sum<<endl;
                  }
                  
                  time_slope_proj[ixy][str]->SetMarkerSize(0.6);
                  time_slope_proj[ixy][str]->SetMarkerStyle(24);
                  //a 		  time_slope_proj[ixy][str]->Draw("same");
                  sprintf(name, "pedfun_%i",str);
                  time_slope_prfit[ixy][str] = new TF1(name, polfunc,0.,31.,nDelayPar-1);
                  time_slope_prfit[ixy][str]->SetParameters(par);
                  time_slope_prfit[ixy][str]->SetLineColor(2);
                  time_slope_prfit[ixy][str]->SetLineWidth(1);
                } // for (int str=0; str=nstrip; str++)
              } // for (int ixy=0; ixy<2; ixy++)
              
              for (int ixy=0; ixy<2; ixy++) {
                for (int str=0; str<nstrip; str++) {
                  if (time_slope_pr[ixy][str]) {
                    delete time_slope_pr[ixy][str]; time_slope_pr[ixy][str]=0;
                  }
                  if (time_slope_proj[ixy][str]) {
                    delete time_slope_proj[ixy][str]; time_slope_proj[ixy][str]=0;
                  }
                  if (time_slope_prfit[ixy][str]) {
                    delete time_slope_prfit[ixy][str]; time_slope_prfit[ixy][str]=0;
                  }
                }
              } // for (int ixy=0; ixy<2; ixy++)
            } // if (ntcor==1) 
#endif
            
            if (ntcor==1 && iiter<nmxiter-1) { //Don't correct for last iteration //150112
              bias_intime_xreso2[iocc] = (xtime_exterr[iocc][iiterrs]->GetMean())*(xtime_exterr[iocc][iiterrs]->GetMean());
              bias_intime_yreso2[iocc] = (ytime_exterr[iocc][iiterrs]->GetMean())*(ytime_exterr[iocc][iiterrs]->GetMean());
              
              // 18th Feb 2015, for estimation of resolution
              //	      if (!isTimeCorOrReso) { 
              //		timeserrx2[iocc] =  time_xrms[iocc]*time_xrms[iocc] - bias_intime_xreso2[iocc]; 
              //		timeserry2[iocc] =  time_yrms[iocc]*time_yrms[iocc] - bias_intime_yreso2[iocc];
              
              timeserrx2[iocc] =  (time_xreso[iocc][iiterrs]->GetRMS())*time_xreso[iocc][iiterrs]->GetRMS() - bias_intime_xreso2[iocc]; 
              timeserry2[iocc] =  (time_yreso[iocc][iiterrs]->GetRMS())*time_yreso[iocc][iiterrs]->GetRMS() - bias_intime_yreso2[iocc];
              if (timeserrx2[iocc]<0.25) timeserrx2[iocc]=0.25; //resolution should be more than 0.5 ns
              if (timeserry2[iocc]<0.25) timeserry2[iocc]=0.25;
            }
            
            if (iiter<nmxiter && iocc<nlayer) {
              for (int jk=0; jk<nstrip; jk++) {
                if (isalign >0 ) {
                  double ncont= time_xstrreso[iocc][jk][iiterrs]->Integral();
                  double statmean = time_xstrreso[iocc][jk][iiterrs]->GetMean();
                  double statrms = time_xstrreso[iocc][jk][iiterrs]->GetRMS();
                  double diff = abs(time_xstrreso[iocc][jk][iiterrs]->GetMean()-fitmean[jk][0]);
                  double width = sqrt(max(0.25, fitrms[jk][0]*fitrms[jk][0] -  bias_intime_xreso2[iocc] )); //time_xstrreso[iocc][jk][iiterrs]->GetRMS();
                  
                  if (ncont >nmnentry && iocc>=xtcorstr && iocc<=xtcorend && ntcor==1) {
                    if (iiter<nmxiter-1) { //Do not update for last iteration
                      timeinuse[0][iocc][jk] = true;
                      if (width<absWidth && diff<widthScale*width) {
                        if (isTimeCorOrReso) { 
                          xtoffset[iocc][jk] +=fitmean[jk][0];//statmean; // time_xstrreso[iocc][jk][iiterrs]->GetMean(); // fitmean[jk][0];
                        }
                      } else {
                        timeinuse[0][iocc][jk] = false;
                      }
                    }
                    correction_xtime[iocc]->Fill(jk, iiter, statmean);
                    fitted_rms_xtime[iocc]->Fill(jk, iiter, fitrms[jk][0]);
                    
                    shift_time_mnft[iiter]->Fill(nstrip*iocc+jk, statmean-fitmean[jk][0]);
                    statmean_time[iiter]->Fill(nstrip*iocc+jk, statmean);
                    statrms_time[iiter]->Fill(nstrip*iocc+jk, statrms);
                    statskew_time[iiter]->Fill(nstrip*iocc+jk, time_xstrreso[iocc][jk][iiterrs]->GetSkewness());
                    statkurt_time[iiter]->Fill(nstrip*iocc+jk, time_xstrreso[iocc][jk][iiterrs]->GetKurtosis());
                    
                    time_offset[iiter]->Fill(nstrip*iocc+jk, xtoffset[iocc][jk]);
                    rms_time[iiter]->Fill(nstrip*iocc+jk, fitrms[jk][0]);
                    if (timeinuse[0][iocc][jk]) rms_timeused[iiter]->Fill(nstrip*iocc+jk, fitrms[jk][0]);
                    
                    //////////////////////////////////
                    shift_time_mnftx[iiter]->Fill(statmean-fitmean[jk][0]);
                    statmean_timex[iiter]->Fill(statmean);
                    statrms_timex[iiter]->Fill(statrms);
                    statskew_timex[iiter]->Fill(time_xstrreso[iocc][jk][iiterrs]->GetSkewness());
                    statkurt_timex[iiter]->Fill(time_xstrreso[iocc][jk][iiterrs]->GetKurtosis());
                    
                    if (abs(xtoffset[iocc][jk])>1.e-6) time_offsetx[iiter]->Fill(xtoffset[iocc][jk]);
                    rms_timex[iiter]->Fill(fitrms[jk][0]);
                    if (timeinuse[0][iocc][jk]) rms_timeusedx[iiter]->Fill(fitrms[jk][0]);
                  }
                  
                  ncont=time_ystrreso[iocc][jk][iiterrs]->Integral();
                  statmean = time_ystrreso[iocc][jk][iiterrs]->GetMean();
                  statrms = time_ystrreso[iocc][jk][iiterrs]->GetRMS();
                  diff = abs(time_ystrreso[iocc][jk][iiterrs]->GetMean()-fitmean[jk][1]);
                  width = sqrt(max(0.25, fitrms[jk][1]*fitrms[jk][1] -  bias_intime_yreso2[iocc]));
                  
                  if (ncont>nmnentry && iocc>=ytcorstr && iocc<=ytcorend && ntcor==1) {
                    if (iiter<nmxiter-1) { //Do not update for last iteration
                      timeinuse[1][iocc][jk] = true;
                      if (width<absWidth && diff<widthScale*width) {
                        if (isTimeCorOrReso) { 
                          ytoffset[iocc][jk] +=fitmean[jk][1];//statmean; // time_ystrreso[iocc][jk][iiterrs]->GetMean(); // fitmean[jk][1];
                          //			} else if (diff<0.50*width) {
                          //			  ytoffset[iocc][jk] +=0.5*(statmean+fitmean[jk][1]); 
                        }
                      } else {
                        timeinuse[1][iocc][jk] = false;
                        //			  ytoffset[iocc][jk] +=statmean;  
                      }
                    }
                    
                    correction_ytime[iocc]->Fill(jk, iiter, statmean);
                    fitted_rms_ytime[iocc]->Fill(jk, iiter, fitrms[jk][1]);
                    
                    shift_time_mnft[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, statmean-fitmean[jk][1]);
                    statmean_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, statmean);
                    statrms_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, statrms);
                    statskew_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, time_ystrreso[iocc][jk][iiterrs]->GetSkewness());
                    statkurt_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, time_ystrreso[iocc][jk][iiterrs]->GetKurtosis());
                    
                    time_offset[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, ytoffset[iocc][jk]);
                    rms_time[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, fitrms[jk][1]);
                    if (timeinuse[1][iocc][jk]) rms_timeused[iiter]->Fill(nstrip*nlayer+nstrip*iocc+jk, fitrms[jk][1]);
                    
                    shift_time_mnfty[iiter]->Fill(statmean-fitmean[jk][1]);
                    statmean_timey[iiter]->Fill(statmean);
                    statrms_timey[iiter]->Fill(statrms);
                    statskew_timey[iiter]->Fill(time_ystrreso[iocc][jk][iiterrs]->GetSkewness());
                    statkurt_timey[iiter]->Fill(time_ystrreso[iocc][jk][iiterrs]->GetKurtosis());
                    
                    if (abs(ytoffset[iocc][jk])>1.e-6) time_offsety[iiter]->Fill(ytoffset[iocc][jk]);
                    rms_timey[iiter]->Fill(fitrms[jk][1]);
                    if (timeinuse[1][iocc][jk]) rms_timeusedy[iiter]->Fill(fitrms[jk][1]);
                  }
                } //if (isalign >0 ) 
                if (ntcor==1 && iiter<nmxiter-1) { //Don't update with last iteration
                  
                  xtrms[iocc][jk] =fitrms[jk][0];
                  ytrms[iocc][jk] =fitrms[jk][1];
                }
              } //for (int jk=0; jk<nstrip; jk++)
              
              time_mean_reso->Fill(iocc, iiterrs, time_xreso[iocc][iiterrs]->GetMean());
              time_rms_reso->Fill(iocc, iiterrs, time_xreso[iocc][iiterrs]->GetRMS());
              time_mean_reso->Fill(nlayer+iocc, iiterrs, time_yreso[iocc][iiterrs]->GetMean());
              time_rms_reso->Fill(nlayer+iocc, iiterrs, time_yreso[iocc][iiterrs]->GetRMS());
              
              time_corrms_reso->Fill(iocc, iiterrs, sqrt(timeserrx2[iocc]));
              time_corrms_reso->Fill(nlayer+iocc, iiterrs, sqrt(timeserry2[iocc]));
              
              time_exterr_reso->Fill(iocc, iiterrs, sqrt(bias_intime_xreso2[iocc]));
              time_exterr_reso->Fill(nlayer+iocc, iiterrs, sqrt(bias_intime_yreso2[iocc]));
              
              file_out <<"time+pos " <<iiterrs<<" iocc "<<iocc<<" "
                       << setw(6)<<xoff[iocc]<<"+-"<<setw(6)<<xrms[iocc]<<" "
                       << setw(6)<<yoff[iocc]<<"+-"<<setw(6)<<yrms[iocc]<<" X "
                       << setw(6)<<time_xrms[iocc]<<" "<<xtime_exterr[iocc][iiterrs]->GetMean()<<" "
                       << setw(6)<<time_xreso[iocc][iiterrs]->GetMean()<<"+-"
                       << setw(6)<<time_xreso[iocc][iiterrs]->GetRMS()<<" "<<sqrt(timeserrx2[iocc])<<" Y "
                       << setw(6)<<time_yrms[iocc]<<" "<<ytime_exterr[iocc][iiterrs]->GetMean()<<" "
                       << setw(6)<<time_yreso[iocc][iiterrs]->GetMean()<<"+-"
                       << setw(6)<<time_yreso[iocc][iiterrs]->GetRMS()<<" "<<sqrt(timeserry2[iocc])<<endl;
              
              file_out <<"// X-rms "<<iocc<<" "<< iiter<<endl;
              for (int jk=0; jk<nstrip; jk++) {
                file_out <<xtrms[iocc][jk]<<", ";
                if ((jk+1)%8==0) file_out<<endl;
              }
              
              file_out <<"// Y-rms "<<iocc<<" "<< iiter<<endl;
              for (int jk=0; jk<nstrip; jk++) {
                file_out <<ytrms[iocc][jk]<<", ";
                if ((jk+1)%8==0) file_out<<endl;
              }
              
              file_out <<"// X-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;
              for (int jk=0; jk<nstrip; jk++) {
                file_out <<xtoffset[iocc][jk]<<", ";
                if ((jk+1)%8==0) file_out<<endl;
              }
                  
              file_out <<"// Y-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;
              
              for (int jk=0; jk<nstrip; jk++) {
                file_out <<ytoffset[iocc][jk]<<", ";
                if ((jk+1)%8==0) file_out<<endl;
              }
              
              cout <<"// X-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;
              for (int jk=0; jk<nstrip; jk++) {
                if (abs(xtoffset[iocc][jk])>40.0) cout <<" xtoffset "<<jk <<" "<< xtoffset[iocc][jk]<<" "<< fitmean[jk][0]<<endl;
              }
              for (int jk=0; jk<nstrip; jk++) {
                cout <<xtoffset[iocc][jk]<<", ";
                if ((jk+1)%8==0) cout<<endl;
              }
              cout <<"// Y-str "<<iocc<<" "<< iiter<<" "<<iiterrs<<endl;
              for (int jk=0; jk<nstrip; jk++) {
                if (abs(ytoffset[iocc][jk])>40.0) cout <<" ytoffset "<<jk <<" "<< ytoffset[iocc][jk]<<" "<< fitmean[jk][1]<<endl;
              }
              for (int jk=0; jk<nstrip; jk++) {
                cout <<ytoffset[iocc][jk]<<", ";
                if ((jk+1)%8==0) cout<<endl;
              }
            } //if (isalign >0 && iiter<nmxiter-1 && iocc<nlayer)
          } //if (isTiming)
        } // for (int iocc = istr; iocc<iend; iocc++)
      } // for (int laye=0; laye<nlayerit; laye++) 
    } // or (int ntcor=0; ntcor< (iiter<nmxiter-1) ? 2 : 1; ntcor++)
    
    if (isalign>0) {
      for (int ij=0; ij<nlayer; ij++) {
        shift_pos->Fill(ij, iiter, xoff[ij]);
        shift_pos->Fill(nlayer+ij, iiter, yoff[ij]);
        
        rms_pos->Fill(ij, iiter, xrms[ij]);
        rms_pos->Fill(nlayer+ij, iiter, yrms[ij]);
      }
    }
  } // for (int iiter=0; iiter<nmxiter; iiter++)
  if (isalign>0) {
    fileOut->cd();
    file_out <<endl;
    file_out<<"double xoff[nlayer] = {";
    for (int ij=0; ij<nlayer; ij++) {
      if ( ij<nlayer-1) {
        file_out<<xoff[ij]<<", ";
      } else {
        file_out<<xoff[ij]<<"};"<<endl;
      }
    }
    
    file_out<<"double yoff[nlayer] = {";
    for (int ij=0; ij<nlayer; ij++) {
      if ( ij<nlayer-1) {
        file_out<<yoff[ij]<<", ";
      } else {
        file_out<<yoff[ij]<<"};"<<endl;
      }
    }
    file_out<<endl;
    
    file_out<<"double xposerrsq[nmxhits][nlayer] = {";
    for (int ix=0; ix<nmxhits; ix++) {
      file_out<<"{"; 
      for (int ij=0; ij<nlayer; ij++) {
        if ( ij<nlayer-1) {
          file_out<<xposerrsq[ix][ij]<<", ";
        } else if (ix<nmxhits-1) {
          file_out<<xposerrsq[ix][ij]<<"},"<<endl;
        } else {
          file_out<<xposerrsq[ix][ij]<<"}};"<<endl;
        }
      }
    }
    
    file_out<<"double yposerrsq[nmxhits][nlayer] = {";
    for (int ix=0; ix<nmxhits; ix++) {
      file_out<<"{"; 
      for (int ij=0; ij<nlayer; ij++) {
        if ( ij<nlayer-1) {
          file_out<<yposerrsq[ix][ij]<<", ";
        } else if (ix<nmxhits-1) {
          file_out<<yposerrsq[ix][ij]<<"},"<<endl;
        } else {
          file_out<<yposerrsq[ix][ij]<<"}};"<<endl;
        }
      }
    }
    
    if (isTiming) {
      file_out<<"double timeoffsetx[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
        if (ij<nlayer-1) {
          file_out<<timeoffsetx[ij]<<", ";
        } else {
          file_out<<timeoffsetx[ij]<<"};"<<endl;
        }
      }
      
      file_out<<"double timeoffsety[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
        if ( ij<nlayer-1) {
          file_out<<timeoffsety[ij]<<", ";
        } else {
          file_out<<timeoffsety[ij]<<"};"<<endl;
        }
      }
      file_out<<endl;
      
      file_out<<"double timeserrx2[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
        if ( ij<nlayer-1) {
          file_out<<timeserrx2[ij]<<", ";
        } else {
          file_out<<timeserrx2[ij]<<"};"<<endl;
        }
      }
      file_out<<endl;
      
      file_out<<"double timeserry2[nlayer] = {";
      for (int ij=0; ij<nlayer; ij++) {
        if ( ij<nlayer-1) {
          file_out<<timeserry2[ij]<<", ";
        } else {
              file_out<<timeserry2[ij]<<"};"<<endl;
        }
      }
      file_out<<endl;
      
      file_out <<"double xtoffset[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
        file_out <<"// X-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
          if (ij<nlayer-1 || jk<nstrip-1) {
            file_out <<xtoffset[ij][jk]<<", ";
          } else {
            file_out <<xtoffset[ij][jk]<<"};";
          }
          if ((jk+1)%8==0) file_out<<endl;
        }
      }
      file_out<<endl;
      
      file_out <<"double ytoffset[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
        file_out <<"// Y-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
          if (ij<nlayer-1 || jk<nstrip-1) {
            file_out <<ytoffset[ij][jk]<<", ";
          } else {
            file_out <<ytoffset[ij][jk]<<"};";
          }
          if ((jk+1)%8==0) file_out<<endl;
        }
      }
      
      file_out <<"double xtoffystr[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
        file_out <<"// X-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
          if (ij<nlayer-1 || jk<nstrip-1) {
            file_out <<xtoffystr[ij][jk]<<", ";
          } else {
            file_out <<xtoffystr[ij][jk]<<"};";
          }
          if ((jk+1)%8==0) file_out<<endl;
        }
      }
      file_out<<endl;
      
      file_out <<"double ytoffxstr[nlayer][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
        file_out <<"// Y-str "<<ij<<endl;
        for (int jk=0; jk<nstrip; jk++) {
          if (ij<nlayer-1 || jk<nstrip-1) {
            file_out <<ytoffxstr[ij][jk]<<", ";
          } else {
            file_out <<ytoffxstr[ij][jk]<<"};";
          }
          if ((jk+1)%8==0) file_out<<endl;
        }
      }
      file_out<<endl;
      
      file_out <<"double xt_slope_cor[nlayer][nstrip][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
        for (int jk=0; jk<nstrip; jk++) {
          for (int kl=0; kl<nstrip; kl++) {
            if (ij<nlayer-1 || jk<nstrip-1 || kl<nstrip-1) {
              file_out <<xt_slope_cor[ij][jk][kl]<<", ";
            } else {
              file_out <<xt_slope_cor[ij][jk][kl]<<"};";
            }
          }
          file_out <<"// X-str "<<ij<<" "<<jk<<endl;
        }
      }
      file_out<<endl;
      
      
      file_out <<"double yt_slope_cor[nlayer][nstrip][nstrip] = {"<<endl;
      for (int ij=0; ij<nlayer; ij++) {
        for (int jk=0; jk<nstrip; jk++) {
          for (int kl=0; kl<nstrip; kl++) {
            if (ij<nlayer-1 || jk<nstrip-1 || kl<nstrip-1) {
              file_out <<yt_slope_cor[ij][jk][kl]<<", ";
            } else {
              file_out <<yt_slope_cor[ij][jk][kl]<<"};";
            }
          }
          file_out <<"// Y-str "<<ij<<" "<<jk<<endl; 
        }
      }
      file_out<<endl;
    } // if (isTiming)
  } //  if (isalign>0)
  
  //  ps.Close();
  fileOut->cd();
  
#ifdef ISEFFICIENCY  
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=0; jk<2*nmxiter; jk++) {
      inefficiency_corx[ij][jk]->Divide(totalentry[ij][jk]);
      inefficiency_cory[ij][jk]->Divide(totalentry[ij][jk]);
      
      inefficiency_corx_5[ij][jk]->Divide(totalentry_5[ij][jk]);
      inefficiency_cory_5[ij][jk]->Divide(totalentry_5[ij][jk]);
          
      inefficiency_corx_6[ij][jk]->Divide(totalentry_6[ij][jk]);
      inefficiency_cory_6[ij][jk]->Divide(totalentry_6[ij][jk]);
      
      inefficiency_corx_7[ij][jk]->Divide(totalentry_7[ij][jk]);
      inefficiency_cory_7[ij][jk]->Divide(totalentry_7[ij][jk]);
      
      inefficiency_corx_8[ij][jk]->Divide(totalentry_8[ij][jk]);
      inefficiency_cory_8[ij][jk]->Divide(totalentry_8[ij][jk]);
      
      inefficiency_corx_9[ij][jk]->Divide(totalentry_9[ij][jk]);
      inefficiency_cory_9[ij][jk]->Divide(totalentry_9[ij][jk]);
      
      inefficiency_corx_10[ij][jk]->Divide(totalentry_10[ij][jk]);
      inefficiency_cory_10[ij][jk]->Divide(totalentry_10[ij][jk]);
      
      inefficiency_corx_11[ij][jk]->Divide(totalentry_11[ij][jk]);
      inefficiency_cory_11[ij][jk]->Divide(totalentry_11[ij][jk]);
      
      inefficiency_uncx_5[ij][jk]->Divide(totalentry_5[ij][jk]);
      inefficiency_uncy_5[ij][jk]->Divide(totalentry_5[ij][jk]);
      
      inefficiency_uncx_6[ij][jk]->Divide(totalentry_6[ij][jk]);
      inefficiency_uncy_6[ij][jk]->Divide(totalentry_6[ij][jk]);
      
      inefficiency_uncx_7[ij][jk]->Divide(totalentry_7[ij][jk]);
      inefficiency_uncy_7[ij][jk]->Divide(totalentry_7[ij][jk]);
      
      inefficiency_uncx_8[ij][jk]->Divide(totalentry_8[ij][jk]);
      inefficiency_uncy_8[ij][jk]->Divide(totalentry_8[ij][jk]);
      
      inefficiency_uncx_9[ij][jk]->Divide(totalentry_9[ij][jk]);
      inefficiency_uncy_9[ij][jk]->Divide(totalentry_9[ij][jk]);
      
      inefficiency_uncx_10[ij][jk]->Divide(totalentry_10[ij][jk]);
      inefficiency_uncy_10[ij][jk]->Divide(totalentry_10[ij][jk]);
      
      inefficiency_uncx_11[ij][jk]->Divide(totalentry_11[ij][jk]);
      inefficiency_uncy_11[ij][jk]->Divide(totalentry_11[ij][jk]);
      
      triggereffi_x_5[ij][jk]->Divide(totalentry_5[ij][jk]);
      triggereffi_y_5[ij][jk]->Divide(totalentry_5[ij][jk]);
      
      triggereffi_x_6[ij][jk]->Divide(totalentry_6[ij][jk]);
      triggereffi_y_6[ij][jk]->Divide(totalentry_6[ij][jk]);
      
      triggereffi_x_7[ij][jk]->Divide(totalentry_7[ij][jk]);
      triggereffi_y_7[ij][jk]->Divide(totalentry_7[ij][jk]);
      
      triggereffi_x_8[ij][jk]->Divide(totalentry_8[ij][jk]);
      triggereffi_y_8[ij][jk]->Divide(totalentry_8[ij][jk]);
      
      triggereffi_x_9[ij][jk]->Divide(totalentry_9[ij][jk]);
      triggereffi_y_9[ij][jk]->Divide(totalentry_9[ij][jk]);
      
      triggereffi_x_10[ij][jk]->Divide(totalentry_10[ij][jk]);
      triggereffi_y_10[ij][jk]->Divide(totalentry_10[ij][jk]);
      
      triggereffi_x_11[ij][jk]->Divide(totalentry_11[ij][jk]);
      triggereffi_y_11[ij][jk]->Divide(totalentry_11[ij][jk]);
      
      inefficiency_uncx[ij][jk]->Divide(totalentry[ij][jk]);
      inefficiency_uncy[ij][jk]->Divide(totalentry[ij][jk]);
      
      inefficiency_xpixel[ij][jk]->Divide(inefficiency_xallpixel[ij][jk]);
      inefficiency_ypixel[ij][jk]->Divide(inefficiency_yallpixel[ij][jk]);
      
      triggereffi_x[ij][jk]->Divide(totalentry[ij][jk]);
      triggereffi_y[ij][jk]->Divide(totalentry[ij][jk]);
      
      for (int kl=0; kl<totalentry[ij][jk]->GetNbinsX(); kl++) {
        for (int lm=0; lm <totalentry[ij][jk]->GetNbinsY(); lm++) {
          double antot = totalentry[ij][jk]->GetBinContent(kl+1, lm+1);
          if (antot <1) { antot = 1;};
          double xtot = triggereffi_xevt[ij][jk]->GetBinContent(kl+1, lm+1);
          double ytot = triggereffi_yevt[ij][jk]->GetBinContent(kl+1, lm+1);
          double effix = xtot/antot;
          if(effix>1.) { effix==1.;}
          if(effix<0.) { effix==0.;}
          double erreffix = pow(effix*(1-effix)/antot, 0.5);
          
          double effiy = ytot/antot;
          if(effiy>1.) { effiy==1.;}
          if(effiy<0.) { effiy==0.;}
          double erreffiy = pow(effiy*(1-effiy)/antot, 0.5);
          
          triggereffi_xevt[ij][jk]->SetBinContent(kl+1, lm+1, effix);
          triggereffi_xevt[ij][jk]->SetBinError(kl+1, lm+1, erreffix);
          
          triggereffi_yevt[ij][jk]->SetBinContent(kl+1, lm+1, effiy);
          triggereffi_yevt[ij][jk]->SetBinError(kl+1, lm+1, erreffiy);
        }
      }
      
      //   triggereffi_xevt[ij][jk]->Divide(totalentry[ij][jk]);
      //  triggereffi_yevt[ij][jk]->Divide(totalentry[ij][jk]);
      
      inefficiency_xt[ij][jk]->Divide(totalentry[ij][jk]); // (total_xt[ij][jk]);
      inefficiency_yt[ij][jk]->Divide(totalentry[ij][jk]); // (total_yt[ij][jk]);
      
      difefficiency_uncx[ij][jk]->Add(inefficiency_uncx[ij][jk], defefficiency_uncx[ij], 1., -1.);
      difefficiency_uncy[ij][jk]->Add(inefficiency_uncy[ij][jk], defefficiency_uncy[ij], 1., -1.);
      
      diftriggereffi_x[ij][jk]->Add(triggereffi_x[ij][jk], deftriggereffi_x[ij], 1., -1.);
      diftriggereffi_y[ij][jk]->Add(triggereffi_y[ij][jk], deftriggereffi_y[ij], 1., -1.);
      
      diftriggereffi_xevt[ij][jk]->Add(triggereffi_xevt[ij][jk], deftriggereffi_x[ij], 1., -1.);
      diftriggereffi_yevt[ij][jk]->Add(triggereffi_yevt[ij][jk], deftriggereffi_y[ij], 1., -1.);
      
      difefficiency_xt[ij][jk]->Add(inefficiency_uncx[ij][jk], inefficiency_xt[ij][jk], 1., -1.);
      difefficiency_yt[ij][jk]->Add(inefficiency_uncy[ij][jk], inefficiency_yt[ij][jk], 1., -1.);
    }
  }
  if (isTiming) { 
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<nmxiter+6; jk++) {
        xystr_xtdev[ij][jk]->Divide(nxystr_xtdev[ij][jk]);
        xystr_ytdev[ij][jk]->Divide(nxystr_ytdev[ij][jk]);
      }
    }
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 0NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 1NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX_5[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx_5[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx_5[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx_5[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY_5[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory_5[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory_5[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory_5[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 2NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX_6[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx_6[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx_6[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx_6[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY_6[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory_6[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory_6[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory_6[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 3NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX_7[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx_7[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx_7[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx_7[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY_7[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory_7[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory_7[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory_7[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 4NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX_8[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx_8[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx_8[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx_8[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY_8[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory_8[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory_8[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory_8[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 5NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX_9[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx_9[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx_9[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx_9[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY_9[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory_9[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory_9[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory_9[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 6NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX_10[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx_10[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx_10[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx_10[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY_10[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory_10[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory_10[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory_10[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// 8NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_corX_11[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_corx_11[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_corx_11[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_corx_11[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_corY_11[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_cory_11[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<" ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_cory_11[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_cory_11[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 0NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 1NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX_5[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx_5[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx_5[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx_5[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY_5[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy_5[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy_5[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy_5[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 2NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX_6[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx_6[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx_6[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx_6[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY_6[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy_6[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy_6[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy_6[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 4NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX_7[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx_7[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx_7[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx_7[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY_7[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy_7[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy_7[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy_7[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 4NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX_8[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx_8[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx_8[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx_8[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY_8[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy_8[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy_8[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy_8[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 5NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX_9[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx_9[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx_9[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx_9[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY_9[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy_9[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy_9[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy_9[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 6NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX_10[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx_10[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx_10[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx_10[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY_10[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy_10[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy_10[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy_10[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// Combined 7NMXITER========================"<<jk<<endl;
    file_out<<"  double ineffi_table_uncX_11[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncx_11[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer="<<ij<<" iter "<<jk<<"unc xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncx_11[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncx_11[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_uncY_11[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_uncy_11[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer="<<ij<<" iter "<<jk<<"unc ystrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_uncy_11[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<inefficiency_uncy_11[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 0Timing ===================="<<jk<<endl;
    file_out<<"  double ineffi_table_X_time[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_xt[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// XlayerTime ineffi="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_xt[ij][jk]->GetNbinsY(); iy++) {
          //Remove position ininefficincy to have extra inefficiency due to timing
          //Note that inefficiency_xt is an EFIICIENCY but inefficiency_cor/uncx is INEFFICIENCY
          file_out <<setprecision(4)<<
            inefficiency_xt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
          //	    -inefficiency_corx[ij][jk]->GetBinContent(ix+1, iy+1)
          //	    -inefficiency_uncx[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double ineffi_table_Y_time[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<inefficiency_yt[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// YlayerTime ineffi="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<inefficiency_yt[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<
            inefficiency_yt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
          //	    -inefficiency_cory[ij][jk]->GetBinContent(ix+1, iy+1)
          //	    -inefficiency_uncy[ij][jk]->GetBinContent(ix+1, iy+1) <<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"//NMXITER 1Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 2Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X_5[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x_5[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x_5[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x_5[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y_5[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y_5[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y_5[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y_5[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"//NMXITER 3Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X_6[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x_6[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x_6[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x_6[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y_6[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y_6[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y_6[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y_6[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 4Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X_7[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x_7[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x_7[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x_7[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y_7[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y_7[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y_7[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y_7[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 5Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X_8[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x_8[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x_8[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x_8[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y_8[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y_8[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y_8[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y_8[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 6Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X_9[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x_9[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x_9[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x_9[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y_9[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y_9[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y_9[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y_9[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 7Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X_10[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x_10[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x_10[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x_10[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y_10[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y_10[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y_10[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y_10[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 8Trigger========================"<<jk<<endl;
    file_out<<"  double effi_trig_X_11[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_x_11[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_x_11[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_x_11[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_Y_11[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_y_11[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer Trig="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_y_11[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_y_11[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
  
  
  for (int jk=0; jk<2*nmxiter; jk++) {
    file_out<<"// NMXITER 0TriggerEVT========================"<<jk<<endl;
    file_out<<"  double effi_trig_xevt[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_xevt[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Xlayer TrigEVT="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_xevt[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_xevt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
    
    file_out<<"  double effi_trig_yevt[nlayer][nstrip][nstrip]={"<<endl;
    for (int ij=0; ij<nlayer; ij++) {
      for (int ix=0; ix<triggereffi_yevt[ij][jk]->GetNbinsX(); ix++) {
        file_out<<"// Ylayer TrigEVT="<<ij<<" iter "<<jk<<" xstrip "<<ix<<endl; 
        for (int iy=0; iy<triggereffi_yevt[ij][jk]->GetNbinsY(); iy++) {
          file_out <<setprecision(4)<<triggereffi_yevt[ij][jk]->GetBinContent(ix+1, iy+1)<<",";
        }
        file_out<<endl;
      }
    }
    file_out<<"};"<<endl;
  }
#endif
  int nbinsx = costhe[0]->GetNbinsX();
  
  for (int ix=0; ix<5; ix++) {
    file_out <<"double ExpCount"<<ix<<"["<<nbinsx<<"]={";
    for (int ij=0; ij<nbinsx; ij++) {
      file_out <<setprecision(5)<<costhe[ix]->GetBinContent(ij+1)<<",";
    }
    file_out <<"};"<<endl;
  }
  
  int nbinsxx = costhe[5]->GetNbinsX();//10Nov
  file_out <<"double ExpCount"<<5<<"["<<nbinsxx<<"]={";
  for (int ij=0; ij<nbinsxx; ij++) {
    file_out <<setprecision(5)<<costhe[5]->GetBinContent(ij+1)<<",";
  }
  file_out <<"};"<<endl;
  
  if (isTiming) { 
    file_out<<"delta T covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatcov[ix][iy]<<", \t";
      }
      file_out<<endl;
    }
    file_out<<endl;
    file_out<<"delta T covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatCount[ix][iy]<<",";
      }
      file_out<<endl;
    }
    file_out<<endl;
    
    file_out<<"delta Tcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatcov[ix][iy]/TMath::Max(1,deltatCount[ix][iy])<<",";
        h_deltatcov->Fill(ix, iy, deltatcov[ix][iy]/TMath::Max(1,deltatCount[ix][iy]));
      }
      file_out<<endl;
    }
    file_out<<endl;
    
    file_out<<"delta T2 covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatcov2[ix][iy]<<", \t";
      }
      file_out<<endl;
    }
    file_out<<endl;
    file_out<<"delta T2 covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatCount2[ix][iy]<<",";
      }
      file_out<<endl;
    }
    file_out<<endl;
    
    file_out<<"delta T2covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatcov2[ix][iy]/TMath::Max(1,deltatCount2[ix][iy])<<",";
        h_deltatcov2->Fill(ix, iy, deltatcov2[ix][iy]/TMath::Max(1,deltatCount2[ix][iy]));
      }
      file_out<<endl;
    }
    file_out<<endl;
    
    //Y-side----------------------------------------------
    file_out<<"delta T covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatcovy[ix][iy]<<", \t";
      }
      file_out<<endl;
    }
    file_out<<endl;
    file_out<<"delta T covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatCounty[ix][iy]<<",";
      }
      file_out<<endl;
    }
    file_out<<endl;
    
    file_out<<"delta Tcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
    for(int ix=0;ix<nlayer;ix++){
      for(int iy=0;iy<nlayer;iy++){
        file_out<<deltatcovy[ix][iy]/TMath::Max(1,deltatCounty[ix][iy])<<",";
        h_deltatcovy->Fill(ix, iy, deltatcovy[ix][iy]/TMath::Max(1,deltatCounty[ix][iy]));
      }
      file_out<<endl;
    }
    file_out<<endl;
  }
  
  //Covariaance of X-pos
  file_out<<"delta Xpos covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposxcov[ix][iy]<<", \t";
    }
    file_out<<endl;
  }
  file_out<<endl;
  
  file_out<<"delta Xpos covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposxCount[ix][iy]<<",";
    }
    file_out<<endl;
  }
  file_out<<endl;
  
  file_out<<"delta Xposcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposxcov[ix][iy]/TMath::Max(1,deltaposxCount[ix][iy])<<",";
      h_deltaposxcov->Fill(ix, iy, deltaposxcov[ix][iy]/TMath::Max(1,deltaposxCount[ix][iy]));
    }
    file_out<<endl;
  }
  file_out<<endl;  
  
  //Covariaance of Y-pos
  file_out<<"delta Ypos covariance ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposycov[ix][iy]<<", \t";
    }
    file_out<<endl;
  }
  file_out<<endl;
  file_out<<"delta Ypos covariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposyCount[ix][iy]<<",";
    }
    file_out<<endl;
  }
  file_out<<endl;
  
  file_out<<"delta Yposcovariance count ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<deltaposycov[ix][iy]/TMath::Max(1,deltaposyCount[ix][iy])<<",";
      h_deltaposycov->Fill(ix, iy, deltaposycov[ix][iy]/TMath::Max(1,deltaposyCount[ix][iy]));
    }
    file_out<<endl;
  }
  file_out<<endl;  
  
  for (int ij=0; ij<8; ij++) {
    for (int jk=0; jk<nlayer; jk++) {
      for (int kl=0; kl<nstrip; kl++) {
        if (!posinuse[ij][jk][kl]) {file_out<<"posinuse["<<ij<<"]["<<jk<<"]["<<kl<<"] ="<<posinuse[ij][jk][kl]<<";"<<endl;}
      }
    }
  } 
  
  for (int ij=0; ij<8; ij++) {
    for (int jk=0; jk<nlayer; jk++) {
      for (int kl=0; kl<nstrip; kl++) {
        if (!timeinuse[ij][jk][kl]) {file_out<<"timeinuse["<<ij<<"]["<<jk<<"]["<<kl<<"] ="<<timeinuse[ij][jk][kl]<<";"<<endl;}
      }
    }
  } 
  if (isTiming) {
    for (int ij=0; ij<nlayer-1; ij++) {
      for (int jk=ij+1; jk<nlayer; jk++) {
        for (int kl=0; kl<=nmxiter; kl++) {
          double ent = timex_correl[ij][jk][kl]->GetEntries();
          if (ent>0) {
            double mean = timex_correl[ij][jk][kl]->GetMean();
            double rms = timex_correl[ij][jk][kl]->GetRMS();
            time_both_cormean[kl]->Fill(ij, jk, mean);
            time_both_corrms[kl]->Fill(ij, jk, rms);
          }
          ent = timey_correl[ij][jk][kl]->GetEntries();
          if (ent>0) {
            double mean = timey_correl[ij][jk][kl]->GetMean();
            double rms = timey_correl[ij][jk][kl]->GetRMS();
            time_both_cormean[kl]->Fill(jk+1, ij, mean);
            time_both_corrms[kl]->Fill(jk+1, ij, rms);
          }
        }
      }
    }
    
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<=nmxiter; jk++) {
        double ent = timex_shift[ij][jk]->GetEntries();
        if (ent>0) {
          double mean = timex_shift[ij][jk]->GetMean();
          double rms = timex_shift[ij][jk]->GetRMS();
          time_both_cormean[jk]->Fill(ij, ij, mean);
          time_both_corrms[jk]->Fill(ij, ij, rms);
        }
        
        ent = timey_shift[ij][jk]->GetEntries();
        if (ent>0) {
          double mean = timey_shift[ij][jk]->GetMean();
          double rms = timey_shift[ij][jk]->GetRMS();
          time_both_cormean[jk]->Fill(ij+1, ij, mean);
          time_both_corrms[jk]->Fill(ij+1, ij, rms);
        }
      }
    }
    
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=0; jk<npixel+1; jk++) {
        for (int kl=0; kl<=nmxiter; kl++) {
          double ent = timexy_correl[ij][jk][kl]->GetEntries();
          if (ent >10) {
            double mean = timexy_correl[ij][jk][kl]->GetMean();
            double rms = timexy_correl[ij][jk][kl]->GetRMS();
            timexy_cormean[jk]->Fill(ij, kl, mean);
            timexy_corrms[jk]->Fill(ij, kl, rms);
          }
        }
      }
    }
  } //if (isTiming)
  
  for (int ij=0; ij<nlayer; ij++) {
    for (int jk=1; jk<=nstrip; jk++) {
      double xx = max(1., h_raw_xcorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
        double yy = h_raw_xcorstrips[ij]->GetBinContent(jk, kl)/xx;
        h_raw_xcorstrips[ij]->SetBinContent(jk, kl, yy);
      }
      xx = max(1., h_raw_ycorstrips[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
        double yy = h_raw_ycorstrips[ij]->GetBinContent(jk, kl)/xx;
        h_raw_ycorstrips[ij]->SetBinContent(jk, kl, yy);
      }
      
      xx = max(1., h_raw_xstrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
        double yy = h_raw_xstrpnhits[ij]->GetBinContent(jk, kl)/xx;
        h_raw_xstrpnhits[ij]->SetBinContent(jk, kl, yy);
      }
      xx = max(1., h_raw_ystrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
        double yy = h_raw_ystrpnhits[ij]->GetBinContent(jk, kl)/xx;
        h_raw_ystrpnhits[ij]->SetBinContent(jk, kl, yy);
      }
      
      xx = max(1., h_raw_xystrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
        double yy = h_raw_xystrpnhits[ij]->GetBinContent(jk, kl)/xx;
        h_raw_xystrpnhits[ij]->SetBinContent(jk, kl, yy);
      }
      xx = max(1., h_raw_yxstrpnhits[ij]->GetBinContent(jk,0));
      for (int kl=1; kl<=nstrip; kl++) {
        double yy = h_raw_yxstrpnhits[ij]->GetBinContent(jk, kl)/xx;
        h_raw_yxstrpnhits[ij]->SetBinContent(jk, kl, yy);
      }
    }
  }
  
  if (!isOnlyCom) { 
    for (int ij=0; ij<nlayer; ij++) {
      for (int jk=1; jk<=nstrip; jk++) {
        double xx = max(1., h_xcorstrips[ij]->GetBinContent(jk,0));
        for (int kl=1; kl<=nstrip; kl++) {
          double yy = h_xcorstrips[ij]->GetBinContent(jk, kl)/xx;
          double yer = h_xcorstrips[ij]->GetBinError(jk, kl)/xx;
          h_xcorstrips[ij]->SetBinContent(jk, kl, yy);
          h_xcorstrips[ij]->SetBinError(jk, kl, yer);
        }
        
        xx = max(1., h_ycorstrips[ij]->GetBinContent(jk,0));
        for (int kl=1; kl<=nstrip; kl++) {
          double yy = h_ycorstrips[ij]->GetBinContent(jk, kl)/xx;
          double yer = h_ycorstrips[ij]->GetBinError(jk, kl)/xx;
          h_ycorstrips[ij]->SetBinContent(jk, kl, yy);
          h_ycorstrips[ij]->SetBinError(jk, kl, yer);
        }
        
        xx = max(1., h_xycorstrips[ij]->GetBinContent(jk,0));
        for (int kl=1; kl<=nstrip; kl++) {
          double yy = h_xycorstrips[ij]->GetBinContent(jk, kl)/xx;
          double yer = h_xycorstrips[ij]->GetBinError(jk, kl)/xx;
          h_xycorstrips[ij]->SetBinContent(jk, kl, yy);
          h_xycorstrips[ij]->SetBinError(jk, kl, yer);
        }
        
        xx = max(1., h_yxcorstrips[ij]->GetBinContent(jk,0));
        for (int kl=1; kl<=nstrip; kl++) {
          double yy = h_yxcorstrips[ij]->GetBinContent(jk, kl)/xx;
          double yer = h_yxcorstrips[ij]->GetBinError(jk, kl)/xx;
          h_yxcorstrips[ij]->SetBinContent(jk, kl, yy);
          h_yxcorstrips[ij]->SetBinError(jk, kl, yer);
        }
        
        if (isTiming) { 
          xx = max(1., h_xtcorstrips[ij]->GetBinContent(jk,0));
          for (int kl=1; kl<=nstrip; kl++) {
            double yy = h_xtcorstrips[ij]->GetBinContent(jk, kl)/xx;
            double yer = h_xtcorstrips[ij]->GetBinError(jk, kl)/xx;
            h_xtcorstrips[ij]->SetBinContent(jk, kl, yy);
            h_xtcorstrips[ij]->SetBinError(jk, kl, yer);
          }
          
          xx = max(1., h_ytcorstrips[ij]->GetBinContent(jk,0));
          for (int kl=1; kl<=nstrip; kl++) {
            double yy = h_ytcorstrips[ij]->GetBinContent(jk, kl)/xx;
            double yer = h_ytcorstrips[ij]->GetBinError(jk, kl)/xx;
            h_ytcorstrips[ij]->SetBinContent(jk, kl, yy);
            h_ytcorstrips[ij]->SetBinError(jk, kl, yer);
          }
          
          xx = max(1., h_xytcorstrips[ij]->GetBinContent(jk,0));
          for (int kl=1; kl<=nstrip; kl++) {
            double yy = h_xytcorstrips[ij]->GetBinContent(jk, kl)/xx;
            double yer = h_xytcorstrips[ij]->GetBinError(jk, kl)/xx;
            h_xytcorstrips[ij]->SetBinContent(jk, kl, yy);
            h_xytcorstrips[ij]->SetBinError(jk, kl, yer);
          }
          
          xx = max(1., h_yxtcorstrips[ij]->GetBinContent(jk,0));
          for (int kl=1; kl<=nstrip; kl++) {
            double yy = h_yxtcorstrips[ij]->GetBinContent(jk, kl)/xx;
            double yer = h_yxtcorstrips[ij]->GetBinError(jk, kl)/xx;
            h_yxtcorstrips[ij]->SetBinContent(jk, kl, yy);
            h_yxtcorstrips[ij]->SetBinError(jk, kl, yer);
          }
        } // if (isTiming)
      } //for (int jk=1; jk<=nstrip; jk++)
    } //for (int ij=0; ij<nlayer; ij++)
  } //   if (!isOnlyCom)
  
  ps.NewPage();
  gStyle->SetOptTitle(1);
  gStyle->SetOptLogy(1);
  gStyle->SetOptFit(101);
  gStyle->SetStatW(.36); //40);
  gStyle->SetStatH(.20); //30);
  gStyle->SetPadLeftMargin(0.07);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.07); //0.03);
  gStyle->SetPadRightMargin(0.01);
  TCanvas* c4c = new TCanvas("c4c", "c4c", 700, 900);
  c4c->Divide(3,4);
  gStyle->SetOptStat(1100);
  
  gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);
  
  if (nmxiter>0) { 
    ps.NewPage();
    int iitr = (isOnlyCom) ? 0 : 2*nmxiter-1; //which iteration should we use
    
    for (int ij=0; ij<nlayer; ij++) {
      c4c->cd(ij+1);
      xlayer_exterr[ij][iitr]->GetXaxis()->SetLabelSize(0.055);
      xlayer_exterr[ij][iitr]->SetLineColor(1); xlayer_exterr[ij][iitr]->Draw();
    }
    c4c->Update(); 
    gStyle->SetStatY(0.79); gStyle->SetStatTextColor(2);
    for (int ij=0; ij<nlayer; ij++) {
      c4c->cd(ij+1);
      ylayer_exterr[ij][iitr]->SetLineColor(2); ylayer_exterr[ij][iitr]->Draw("sames");
    }
    c4c->Update();
    
    if (isTiming) {
      ps.NewPage();
      gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);
      for (int ij=0; ij<nlayer; ij++) {
        c4c->cd(ij+1);
        xtime_exterr[ij][iitr]->GetXaxis()->SetLabelSize(0.055);
        xtime_exterr[ij][iitr]->SetLineColor(1); xtime_exterr[ij][iitr]->Draw();
      }
      c4c->Update(); 
      
      gStyle->SetStatY(0.79); gStyle->SetStatTextColor(2);
      for (int ij=0; ij<nlayer; ij++) {
        c4c->cd(ij+1);
        ytime_exterr[ij][iitr]->SetLineColor(2); ytime_exterr[ij][iitr]->Draw("sames");
      }
      c4c->Update();
    }
  }
  
  int iitr = (isOnlyCom) ? 0 : 2*nmxiter-1; //which iteration should we use
  ps.NewPage();
  gStyle->SetOptStat(1110);
  gStyle->SetStatW(.36);
  gStyle->SetStatH(.18);
  gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
  for (int ij=0; ij<nlayer; ij++) {
    c4c->cd(ij+1);
    gPad->SetLogy(1); 
    strp_xmul[ij][iitr]->GetXaxis()->SetLabelSize(0.07);
    strp_xmul[ij][iitr]->GetYaxis()->SetLabelSize(0.07);
    strp_xmul[ij][iitr]->SetLineColor(1);
    strp_xmul[ij][iitr]->ProjectionY()->Draw();
  } 
  c4c->Update();
  gStyle->SetStatY(.77); gStyle->SetStatTextColor(2);
  TH1F* tmphist[nlayer]={0};
  for (int ij=0; ij<nlayer; ij++) {
    c4c->cd(ij+1);
    gPad->SetLogy(1); 
    tmphist[ij] = (TH1F*)strp_ymul[ij][iitr]->ProjectionY();
    tmphist[ij]->SetLineColor(2);
    tmphist[ij]->Draw("sames");
  }
  
  c4c->Update();
  for (int ij=0; ij<nlayer; ij++) {
    if (tmphist[ij]) { delete tmphist[ij]; tmphist[ij]=0;}
  }
  
  ps.NewPage();
  gStyle->SetStatTextColor(1); gStyle->SetStatY(0.99);
  
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatW(.36); //40);
  gStyle->SetStatH(.28); //30);
  
  TCanvas* c4a = new TCanvas("c4a", "c4a", 700, 900);
  c4a->Divide(3,4);
  
  //  int iitr = 2*nmxiter-1; //which iteration should we use
  TH2F* histy[nlayer]={0};
  const int nmult=6;
  TGraph* gry[nlayer][nmult]={0};
  double content[100][nmult];
  double xvl[100];
  double yvl[100];
  double xerr[100]={0.0};
  double yerr[100];
  
  double total[100];
  const char* namex[nmult]={"null", "ZERO", "ONE", "TWO", "THREE", "More"};
  TLegend *tleg = new TLegend(.70, .62, .90, .92, "","brNDC");
  tleg->SetFillColor(10);
  tleg->SetTextSize(0.06);
  tleg->SetBorderSize(0);
  tleg->SetTextFont(42);
  
  int ishft = (isOnlyCom) ? 0 : nmxiter; //Without any iteration & proper efficiency
  for (int ixx=0; ixx<nmxiter; ixx++) {
    for (int jkl=0; jkl <2; jkl++) {
      if (ixx!=0 || jkl!=0) {ps.NewPage();}
      for (int ij=0; ij<nlayer; ij++) {
        c4a->cd(ij+1);
        switch(jkl) {
        case 0 : histy[ij] = (TH2F*)strp_xmul[ij][ishft+ixx]->Clone(); break;
        case 1 : histy[ij] = (TH2F*)strp_ymul[ij][ishft+ixx]->Clone(); break;
        default : histy[ij] = (TH2F*)strp_xmul[ij][ishft+ixx]->Clone(); break;
        }
        int nbinx=2+histy[ij]->GetNbinsX();
        int nbiny=2+histy[ij]->GetNbinsY();
        
        for (int jk=0; jk<nbinx; jk++) {
          total[jk]=0;
          xvl[jk] = histy[ij]->GetXaxis()->GetBinCenter(jk);
          for (int kl=0; kl<nbiny; kl++) {
            content[jk][kl] = histy[ij]->GetBinContent(jk, kl);
            total[jk] +=content[jk][kl];
          }
          if (total[jk]<1.0) total[jk]=1.0;
        }
        int icol=0;
        for (int kl=1; kl<nmult; kl++) {
          for (int jk=0; jk<nbinx; jk++) {
            yvl[jk] = content[jk][kl]/total[jk];
            yerr[jk] =sqrt(yvl[jk]*(1-yvl[jk])/total[jk]);
          }
          icol++; if (icol%5==0) icol++;
          
          gry[ij][kl]=new TGraphErrors(nbinx, xvl, yvl, xerr, yerr);
          gry[ij][kl]->SetTitle(histy[ij]->GetTitle());
          gry[ij][kl]->SetMarkerStyle(20);
          gry[ij][kl]->SetMarkerSize(0.5);
          gry[ij][kl]->SetMarkerColor(icol);
          gry[ij][kl]->GetXaxis()->SetLabelSize(.064);
          gry[ij][kl]->GetYaxis()->SetLabelSize(.054);
          gry[ij][kl]->GetXaxis()->SetLabelOffset(.001);
          gry[ij][kl]->GetYaxis()->SetLabelOffset(.01);
          gry[ij][kl]->GetXaxis()->SetTitle("Position in strip");
          gry[ij][kl]->GetXaxis()->CenterTitle();
          gry[ij][kl]->GetXaxis()->SetTitleSize(.08);
          gry[ij][kl]->GetXaxis()->SetTitleOffset(.8);
          gry[ij][kl]->GetXaxis()->SetRangeUser(-0.5, 0.5);
          gry[ij][kl]->SetMinimum(0);
          gry[ij][kl]->SetMaximum(0.9);
          //    gry[ij][kl]->SetNdivisions(506,"XY");
          if (kl==1) {
            gry[ij][kl]->Draw("AP");
          } else {
            gry[ij][kl]->Draw("P:same");
          }
          if (ixx==0 && ij==0 && jkl==0) {
            tleg->AddEntry(gry[ij][kl],namex[kl],"p");
          }
        }
        if (ij==0) tleg->Draw();
      }
      c4a->Update();
      for (int ij=0; ij<nlayer; ij++) {
        if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
        for (int jk=0; jk<6; jk++) {
          if (gry[ij][jk]) { delete gry[ij][jk]; gry[ij][jk]=0;}
        }
      }
    } //for (int jkl=0; jkl <2; jkl++)
  } //for (int ixx=0; ixx<nmxiter; ixx++)
  
  //  int iitr = 2*nmxiter-1; //which iteration should we use
  
  ps.NewPage(); 
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1F* histz[nlayer];
  if (nmxiter >0) { 
    latex.SetTextSize(0.08);
    for (int jkl=0; jkl <4; jkl++) { 
      if (isOnlyCom && jkl>=2) continue;
      for (int klm =0; klm<nmxusedhits; klm++) {
        ps.NewPage(); 
        for (int ij=0; ij<nlayer; ij++) {
          c4c->cd(ij+1);
          switch(jkl) {
          case 0 : histz[ij] = (TH1F*)xlayer_reso_mul[ij][2*nmxiter-2][klm]->Clone(); break;
          case 1 : histz[ij] = (TH1F*)ylayer_reso_mul[ij][2*nmxiter-2][klm]->Clone(); break;
          case 2 : histz[ij] = (TH1F*)xlayer_reso_mul[ij][2*nmxiter-1][klm]->Clone(); break;
          case 3 : histz[ij] = (TH1F*)ylayer_reso_mul[ij][2*nmxiter-1][klm]->Clone(); break;
          default : histz[ij] = (TH1F*)ylayer_reso_mul[ij][2*nmxiter-1][klm]->Clone(); break;
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
    
    if (isTiming) { 
      for (int jkl=0; jkl <4; jkl++) { 
        if (isOnlyCom && jkl>=2) continue;
        for (int klm =0; klm<nmxtimehit; klm++) {
          ps.NewPage(); 
          for (int ij=0; ij<nlayer; ij++) {
            c4c->cd(ij+1);
            switch(jkl) {
            case 0 : histz[ij] = (TH1F*)time_mulxreso[ij][2*nmxiter-2][klm]->Clone(); break;
            case 1 : histz[ij] = (TH1F*)time_mulyreso[ij][2*nmxiter-2][klm]->Clone(); break;
            case 2 : histz[ij] = (TH1F*)time_mulxreso[ij][2*nmxiter-1][klm]->Clone(); break;
            case 3 : histz[ij] = (TH1F*)time_mulyreso[ij][2*nmxiter-1][klm]->Clone(); break;
            default : histz[ij] = (TH1F*)time_mulxreso[ij][2*nmxiter-1][klm]->Clone(); break;
            }
            histz[ij]->GetXaxis()->SetLabelSize(.07);
            histz[ij]->GetYaxis()->SetLabelSize(.055);
            histz[ij]->GetXaxis()->SetTitle("#Delta t(ns)");
            histz[ij]->GetXaxis()->CenterTitle();
            histz[ij]->GetXaxis()->SetTitleOffset(0.9);
            histz[ij]->GetXaxis()->SetTitleSize(0.06); 
            
            TFitResultPtr ptr = histz[ij]->Fit("gaus","SQ");
            latex.DrawLatex(0.20, 0.76,Form("%g", int(1000*histz[ij]->GetMean())/1000.));
            latex.DrawLatex(0.20, 0.66,Form("%g", int(1000*histz[ij]->GetRMS())/1000.));
            int fitStatus = ptr;
            if (fitStatus==0 && histz[ij]->GetEntries()>3) { 
              latex.DrawLatex(0.70, 0.76,Form("%g", int(1000*ptr->Parameter(1))/1000.));
              latex.DrawLatex(0.70, 0.66,Form("%g", int(1000*ptr->Parameter(2))/1000.));
            }
          }
          c4c->Update();
          for (int ij=0; ij<nlayer; ij++) {
            if (histz[ij]) { delete histz[ij]; histz[ij]=0;}
          }
        }
      }
    }
    latex.SetTextSize(0.12);
  }
  
  ps.NewPage(); 
  gStyle->SetOptStat(1100); 
  gStyle->SetStatH(.30);
  ps.NewPage(); 
  c4c->cd(1); h_xprob->Draw();
  c4c->cd(2); h_reduchisqx->Draw();
  c4c->cd(3); h_xndf->Draw();
  
  c4c->cd(4); h_yprob->Draw();
  c4c->cd(5); h_reduchisqy->Draw();
  c4c->cd(6); h_yndf->Draw();
  
  if (isTiming) {
    c4c->cd(7); h_xtprob->Draw();
    c4c->cd(8); h_treduchisqx->Draw();
    c4c->cd(9); h_txndf->Draw();
    
    c4c->cd(10); h_ytprob->Draw();
    c4c->cd(11); h_treduchisqy->Draw();
    c4c->cd(12); h_tyndf->Draw();
  }
  c4c->Update();
  
  ps.NewPage(); 
  
  TH1F* hist1x[2*nprob];
  TH1F* hist1y[2*nprob];
  for (int jk=0; jk<nprob; jk++) {
    hist1x[jk] = (TH1F*)h_xnprob[jk]->Clone();
    hist1y[jk] = (TH1F*)h_ynprob[jk]->Clone();
    
    hist1x[jk+nprob] = (TH1F*)h_xtnprob[jk]->Clone();
    hist1y[jk+nprob] = (TH1F*)h_ytnprob[jk]->Clone();
  }
  gStyle->SetStatY(0.99); gStyle->SetStatTextColor(1);
  for (int jk=0; jk<2*nprob; jk++) {
    c4c->cd(jk+1);
    hist1x[jk]->SetLineColor(1); hist1x[jk]->Draw();
  }
  
  c4c->Update();
  gStyle->SetStatY(0.84); gStyle->SetStatTextColor(2);
  for ( int jk=0; jk<2*nprob; jk++) {
    c4c->cd(jk+1);
    hist1y[jk]->SetLineColor(2); hist1y[jk]->Draw("sames");
  }
  
  c4c->Update();
  
  for (int jk=0; jk<2*nprob; jk++) {
    if (hist1x[jk]) { delete hist1x[jk]; hist1x[jk]=0;}
    if (hist1y[jk]) { delete hist1y[jk]; hist1y[jk]=0;}
  }
  
  if (!isOnlyCom && plot_level>80) { 
    ps.NewPage(); 
    gStyle->SetStatTextColor(1); gStyle->SetStatY(0.99);
    gStyle->SetOptStat(0);
    for (int iitr=nmxiter-1; iitr<nmxiter; iitr++) {
      for (int jkl=0; jkl <4; jkl++) { 
        for (int klm =0; klm<nstr_posmx; klm++) {
          ps.NewPage(); 
          for (int ij=0; ij<nlayer; ij++) {
            c4a->cd(ij+1);
            switch(jkl) {
            case 0 : histy[ij] = (TH2F*)xstr_xdev[ij][iitr][klm]->Clone(); break;
            case 1 : histy[ij] = (TH2F*)ystr_ydev[ij][iitr][klm]->Clone(); break;
            case 2 : histy[ij] = (TH2F*)ystr_xdev[ij][iitr][klm]->Clone(); break;
            case 3 : histy[ij] = (TH2F*)xstr_ydev[ij][iitr][klm]->Clone(); break;
            default : histy[ij] = (TH2F*)xstr_xdev[ij][iitr][klm]->Clone(); break;
            }
            histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
            histy[ij]->GetXaxis()->SetTitle("Strip No");
            histy[ij]->GetXaxis()->CenterTitle();
            histy[ij]->GetXaxis()->SetTitleOffset(0.8);
            histy[ij]->GetXaxis()->SetTitleSize(0.07); 
            
            histy[ij]->GetYaxis()->SetTitle("#Delta X (Strip pitch)");
            histy[ij]->GetYaxis()->CenterTitle();
            histy[ij]->GetYaxis()->SetTitleOffset(0.8);
            histy[ij]->GetYaxis()->SetTitleSize(0.07); 
            
            histy[ij]->Draw("colz");
            histy[ij]->ProfileX()->Draw("sames");
          }
          c4a->Update();
          for (int ij=0; ij<nlayer; ij++) {
            if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
          }
        }
      }
    }
    
    if (isTiming) { 
      TProfile* str_tdev[nlayer]={0};
      for (int iitr=0; iitr<nmxiter; iitr++) {
        if (iitr>0 && iitr<nmxiter-2) continue;
        for (int jkl=0; jkl <4; jkl++) { 
          
          ps.NewPage(); 
          for (int ij=0; ij<nlayer; ij++) {
            c4a->cd(ij+1);
            switch(jkl) {
            case 0 : histy[ij] = (TH2F*)xstr_xtdev[ij][iitr]->Clone(); break;
            case 1 : histy[ij] = (TH2F*)ystr_ytdev[ij][iitr]->Clone(); break;
            case 2 : histy[ij] = (TH2F*)ystr_xtdev[ij][iitr]->Clone(); break;
            case 3 : histy[ij] = (TH2F*)xstr_ytdev[ij][iitr]->Clone(); break;
            default : histy[ij] = (TH2F*)xstr_xtdev[ij][iitr]->Clone(); break;
            }
            histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
            histy[ij]->Draw("colz");
            sprintf(name, "xxx_%i", ij);
            str_tdev[ij] = (TProfile*)histy[ij]->ProfileX(name);
            str_tdev[ij]->SetMarkerSize(0.36);
            str_tdev[ij]->SetMarkerStyle(24);
            str_tdev[ij]->Fit("pol4","Q","same");
          }
          c4a->Update();
          for (int ij=0; ij<nlayer; ij++) {
            if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
            if (str_tdev[ij]) { delete str_tdev[ij]; str_tdev[ij]=0;}
          }
        }
      }
    }
  }
      
  if (!isOnlyCom && plot_level>70) { 
    for (int lm=0; lm<nmxtimehit; lm++) { //multiplicity
      for (int jkl=0; jkl <12; jkl++) { 
        if (!isTiming && jkl>=4) continue;
        ps.NewPage(); 
        gStyle->SetOptStat(0);
        for (int ij=0; ij<nlayer; ij++) {
          c4a->cd(ij+1);
          switch(jkl) {
            
          case 0 : histy[ij] = (TH2F*)xpos_xdev[ij][lm]->Clone(); break;
          case 1 : histy[ij] = (TH2F*)xpos_ydev[ij][lm]->Clone(); break;
          case 2 : histy[ij] = (TH2F*)ypos_xdev[ij][lm]->Clone(); break;
          case 3 : histy[ij] = (TH2F*)ypos_ydev[ij][lm]->Clone(); break;     
            
          case 4 : histy[ij] = (TH2F*)xpos_xtdev_str[ij][lm]->Clone(); break;
          case 5 : histy[ij] = (TH2F*)ypos_ytdev_str[ij][lm]->Clone(); break;
            
          case 6 : histy[ij] = (TH2F*)xpos_xtdev[ij][lm]->Clone(); break;
          case 7 : histy[ij] = (TH2F*)ypos_ytdev[ij][lm]->Clone(); break;
            
          case 8 : histy[ij] = (TH2F*)xpos_xtdev_glb[ij][lm]->Clone(); break;
          case 9 : histy[ij] = (TH2F*)xpos_ytdev_glb[ij][lm]->Clone(); break;
          case 10 : histy[ij] = (TH2F*)ypos_xtdev_glb[ij][lm]->Clone(); break;	 
          case 11 : histy[ij] = (TH2F*)ypos_ytdev_glb[ij][lm]->Clone(); break;	 
            
            
          default : histy[ij] = (TH2F*)xpos_xdev[ij][lm]->Clone(); break;
          }
          if (jkl >=4 && jkl<=7 ) { 
            histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
          }
          histy[ij]->GetXaxis()->SetLabelSize(.07);
          histy[ij]->GetYaxis()->SetLabelSize(.07);
          
          histy[ij]->Draw("colz");
          histy[ij]->ProfileX()->Draw("sames");
        }
        c4a->Update();
        for (int ij=0; ij<nlayer; ij++) {
          if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
        }
      }
    }
  }
  
  if (isTiming) {
    ps.NewPage();
    for (int ij=0; ij<8; ij++) {
      c4a->cd(ij+1);
      nxypos_xytdev[ij]->GetXaxis()->SetRangeUser(-1., nusedstrip);
      nxypos_xytdev[ij]->GetYaxis()->SetRangeUser(-1., nusedstrip);
      nxypos_xytdev[ij]->Draw("colz");
    }
    for (int ij=0; ij<4; ij++) {
      c4a->cd(ij+9);
      xtdev_ytdev[ij]->Draw("colz");
    }
    c4a->Update();
  }
#ifdef ISEFFICIENCY
  int iitrx = (isOnlyCom) ? 0 : 2*nmxiter-1; //which iteration should we use
  if (plot_level>100) { 
    ps.NewPage();
    gStyle->SetPaintTextFormat("5.2f");
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadRightMargin(0.09);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadLeftMargin(0.08);
    gStyle->SetTitleFontSize(0.04);
    
    TCanvas *c45=new TCanvas ("c45","Strip occupancy",500,700);
    
    for (int ij=0; ij<nlayer; ij++) {
      ps.NewPage();
      c45->cd();
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetLabelSize(.035);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetLabelSize(.035);
      inefficiency_cory[ij][iitrx]->GetZaxis()->SetLabelSize(.035);
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetLabelOffset(-.001);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetLabelOffset(.01);
      
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetTitle("X-Strip No");
      inefficiency_cory[ij][iitrx]->GetXaxis()->CenterTitle();
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetTitleOffset(0.7);
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetTitleSize(0.04); 
      
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetTitle("Y-Strip No");
      inefficiency_cory[ij][iitrx]->GetYaxis()->CenterTitle();
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetTitleOffset(0.5);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetTitleSize(0.04); 
      inefficiency_cory[ij][iitrx]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_cory[ij][iitrx]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_cory[ij][iitrx]->SetMarkerSize(0.62);
      inefficiency_cory[ij][iitrx]->Draw("colz");
      c45->Update();
      ps.NewPage();
      inefficiency_cory[ij][iitrx]->Draw("text25");
      c45->Update();
      
      ps.NewPage();
      c45->cd();
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetLabelSize(.035);
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetLabelSize(.035);
      inefficiency_uncx[ij][iitrx]->GetZaxis()->SetLabelSize(.035);
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetTitle("X-Strip No");
      inefficiency_uncx[ij][iitrx]->GetXaxis()->CenterTitle();
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetTitleOffset(0.7);
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetTitleSize(0.04); 
      
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetTitle("Y-Strip No");
      inefficiency_uncx[ij][iitrx]->GetYaxis()->CenterTitle();
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetTitleOffset(0.5);
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetTitleSize(0.04); 
      
      inefficiency_uncx[ij][iitrx]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncx[ij][iitrx]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncx[ij][iitrx]->SetMarkerSize(0.62);
      inefficiency_uncx[ij][iitrx]->Draw("colz");
      c45->Update();
      ps.NewPage();
      inefficiency_uncx[ij][iitrx]->Draw("text25");
      c45->Update();
      
      ps.NewPage();
      c45->cd();
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetLabelSize(.035);
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetLabelSize(.035);
      inefficiency_uncy[ij][iitrx]->GetZaxis()->SetLabelSize(.035);
      
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetTitle("X-Strip No");
      inefficiency_uncy[ij][iitrx]->GetXaxis()->CenterTitle();
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetTitleOffset(0.7);
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetTitleSize(0.04); 
      
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetTitle("Y-Strip No");
      inefficiency_uncy[ij][iitrx]->GetYaxis()->CenterTitle();
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetTitleOffset(0.5);
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetTitleSize(0.04); 
      
      inefficiency_uncy[ij][iitrx]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncy[ij][iitrx]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
      inefficiency_uncy[ij][iitrx]->SetMarkerSize(0.62);
      inefficiency_uncy[ij][iitrx]->Draw("colz");
      c45->Update();
      ps.NewPage();
      inefficiency_uncy[ij][iitrx]->Draw("text25");
      c45->Update();
    }
    c45->Clear();
  }
  
  //  for (int jkl=0; jkl <26; jkl++) { //GMA
  for (int jkl=0; jkl <75/* 71 29*/; jkl++) {
    if (!isTiming && jkl>54/*50 8*/) continue;
    if (nmxiter==0 && (jkl<=47/*5*/ || jkl>=66/*62 20*/)) continue;
    if (nmxiter<2 && jkl>=70/*66 24*/) continue;
    if (isOnlyCom && jkl>55/* 51 9*/) continue;
    ps.NewPage();
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.11);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetTitleFontSize(0.07);
    double xmx=-1000000;
    double xmn=1000000;
    
    for (int ij=0; ij<nlayer; ij++) {
      switch(jkl) {
        
      case 0 : histy[ij] = (TH2F*)totalentry_5[ij][iitrx]->Clone(); break;
      case 1 : histy[ij] = (TH2F*)totalentry_6[ij][iitrx]->Clone(); break;
      case 2 : histy[ij] = (TH2F*)totalentry_7[ij][iitrx]->Clone(); break;
      case 3 : histy[ij] = (TH2F*)totalentry_8[ij][iitrx]->Clone(); break;
      case 4 : histy[ij] = (TH2F*)totalentry_9[ij][iitrx]->Clone(); break;
      case 5 : histy[ij] = (TH2F*)totalentry_10[ij][iitrx]->Clone(); break;
      case 6 : histy[ij] = (TH2F*)totalentry_11[ij][iitrx]->Clone(); break;
        
      case 7 : histy[ij] = (TH2F*)inefficiency_corx_5[ij][iitrx]->Clone(); break;
      case 8 : histy[ij] = (TH2F*)inefficiency_corx_6[ij][iitrx]->Clone(); break;
      case 9 : histy[ij] = (TH2F*)inefficiency_corx_7[ij][iitrx]->Clone(); break;
      case 10 : histy[ij] = (TH2F*)inefficiency_corx_8[ij][iitrx]->Clone(); break;
      case 11 : histy[ij] = (TH2F*)inefficiency_corx_9[ij][iitrx]->Clone(); break;
      case 12 : histy[ij] = (TH2F*)inefficiency_corx_10[ij][iitrx]->Clone(); break;
      case 13 : histy[ij] = (TH2F*)inefficiency_corx_11[ij][iitrx]->Clone(); break;
        
      case 14 : histy[ij] = (TH2F*)inefficiency_uncx_5[ij][iitrx]->Clone(); break;
      case 15 : histy[ij] = (TH2F*)inefficiency_uncx_6[ij][iitrx]->Clone(); break;
      case 16:  histy[ij] = (TH2F*)inefficiency_uncx_7[ij][iitrx]->Clone(); break;
      case 17 : histy[ij] = (TH2F*)inefficiency_uncx_8[ij][iitrx]->Clone(); break;
      case 18 : histy[ij] = (TH2F*)inefficiency_uncx_9[ij][iitrx]->Clone(); break;
      case 19 : histy[ij] = (TH2F*)inefficiency_uncx_10[ij][iitrx]->Clone(); break;
      case 20 : histy[ij] = (TH2F*)inefficiency_uncx_11[ij][iitrx]->Clone(); break;
        
      case 21 : histy[ij] = (TH2F*)inefficiency_uncy_5[ij][iitrx]->Clone(); break;
      case 22 : histy[ij] = (TH2F*)inefficiency_uncy_6[ij][iitrx]->Clone(); break;
      case 23 : histy[ij] = (TH2F*)inefficiency_uncy_7[ij][iitrx]->Clone(); break;
      case 24 : histy[ij] = (TH2F*)inefficiency_uncy_8[ij][iitrx]->Clone(); break;
      case 25 : histy[ij] = (TH2F*)inefficiency_uncy_9[ij][iitrx]->Clone(); break;
      case 26 : histy[ij] = (TH2F*)inefficiency_uncy_10[ij][iitrx]->Clone(); break;
      case 27 : histy[ij] = (TH2F*)inefficiency_uncy_11[ij][iitrx]->Clone(); break;
        
      case 28 : histy[ij] = (TH2F*)triggereffi_x_5[ij][iitrx]->Clone(); break;
      case 29 : histy[ij] = (TH2F*)triggereffi_x_6[ij][iitrx]->Clone(); break;
      case 30 : histy[ij] = (TH2F*)triggereffi_x_7[ij][iitrx]->Clone(); break;
      case 31 : histy[ij] = (TH2F*)triggereffi_x_8[ij][iitrx]->Clone(); break;
      case 32 : histy[ij] = (TH2F*)triggereffi_x_9[ij][iitrx]->Clone(); break;
      case 33 : histy[ij] = (TH2F*)triggereffi_x_10[ij][iitrx]->Clone(); break;
      case 34 : histy[ij] = (TH2F*)triggereffi_x_11[ij][iitrx]->Clone(); break;
        
      case 35 : histy[ij] = (TH2F*)triggereffi_y_5[ij][iitrx]->Clone(); break;
      case 36 : histy[ij] = (TH2F*)triggereffi_y_6[ij][iitrx]->Clone(); break;
      case 37 : histy[ij] = (TH2F*)triggereffi_y_7[ij][iitrx]->Clone(); break;
      case 38 : histy[ij] = (TH2F*)triggereffi_y_8[ij][iitrx]->Clone(); break;
      case 39 : histy[ij] = (TH2F*)triggereffi_y_9[ij][iitrx]->Clone(); break;
      case 40 : histy[ij] = (TH2F*)triggereffi_y_10[ij][iitrx]->Clone(); break;
      case 41 : histy[ij] = (TH2F*)triggereffi_y_11[ij][iitrx]->Clone(); break;
        
      case 42 : histy[ij] = (TH2F*)inefficiency_corx[ij][iitrx]->Clone(); break;
        
      case 43 : histy[ij] = (TH2F*)totalentry[ij][iitrx]->Clone(); break;
      case 44 : histy[ij] = (TH2F*)inefficiency_cory[ij][iitrx]->Clone(); break;
        
      case 45 : histy[ij] = (TH2F*)inefficiency_uncx[ij][iitrx]->Clone(); break;
      case 46 : histy[ij] = (TH2F*)inefficiency_uncy[ij][iitrx]->Clone(); break;
        
      case 47 : histy[ij] = (TH2F*)triggereffi_xevt[ij][iitrx]->Clone(); break;
      case 48 : histy[ij] = (TH2F*)triggereffi_yevt[ij][iitrx]->Clone(); break;
        
      case 49 : histy[ij] = (TH2F*)triggereffi_x[ij][iitrx]->Clone(); break;
      case 50 : histy[ij] = (TH2F*)triggereffi_y[ij][iitrx]->Clone(); break;
        
      case 51 : histy[ij] = (TH2F*)difefficiency_uncx[ij][iitrx]->Clone(); break;
      case 52 : histy[ij] = (TH2F*)difefficiency_uncy[ij][iitrx]->Clone(); break;
        
      case 53 : histy[ij] = (TH2F*)diftriggereffi_x[ij][iitrx]->Clone(); break;
      case 54 : histy[ij] = (TH2F*)diftriggereffi_y[ij][iitrx]->Clone(); break;
        
      case 55 : histy[ij] = (TH2F*)inefficiency_xt[ij][iitrx]->Clone(); break;
      case 56 : histy[ij] = (TH2F*)inefficiency_yt[ij][iitrx]->Clone(); break;
        
      case 57 : histy[ij] = (TH2F*)nxystr_xtdev[ij][nmxiter]->Clone(); break;
      case 58 : histy[ij] = (TH2F*)nxystr_ytdev[ij][nmxiter]->Clone(); break;
        
      case 59 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter]->Clone(); break;
      case 60 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter]->Clone(); break;
            
      case 61 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+1]->Clone(); break;
      case 62 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+1]->Clone(); break;
        
      case 63 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+2]->Clone(); break;
      case 64 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+2]->Clone(); break;
        
      case 65 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+3]->Clone(); break;
      case 66 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+3]->Clone(); break;
      case 67 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+4]->Clone(); break;
      case 68 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+4]->Clone(); break;
        
      case 69 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter+5]->Clone(); break;
      case 70 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter+5]->Clone(); break;
        
      case 71 : histy[ij] = (TH2F*)xystr_xtdev[ij][0]->Clone(); break;
      case 72 : histy[ij] = (TH2F*)xystr_ytdev[ij][0]->Clone(); break;
        
      case 73 : histy[ij] = (TH2F*)xystr_xtdev[ij][nmxiter-1]->Clone(); break;
      case 74 : histy[ij] = (TH2F*)xystr_ytdev[ij][nmxiter-1]->Clone(); break;
        
      default :histy[ij] = (TH2F*)inefficiency_uncx[ij][2*nmxiter-1]->Clone(); break; 
      }
      double yy1=histy[ij]->GetMaximum();
      if (yy1>xmx) xmx=yy1;
      
      double yy2=histy[ij]->GetMinimum();
      if (yy2<xmn) xmn=yy2;
    }
    if (jkl>=4 && xmn<0.1) { xmn=0.1;} //Avoid no hit area
    
    for (int ij=0; ij<nlayer; ij++) {
      c4a->cd(ij+1);
      //      histy[ij]->SetMinimum(xmn);
      //      histy[ij]->SetMaximum(xmx);
      histy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
      histy[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
      
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
    }
    c4a->Update();
    for (int ij=0; ij<nlayer; ij++) {
      if (histy[ij]) { delete histy[ij]; histy[ij]=0;}
    }
  }
#endif //define ISEFFICIENCY
  if (!isOnlyCom) { 
    TH2F* histyy[nlayer]={0};
    for (int ixx=0; ixx<8; ixx++) {
      if (!isTiming && ixx>=4) continue;
      ps.NewPage();
      
      for (int ij=0; ij<nlayer; ij++) {
        c4a->cd(ij+1);
        switch(ixx) {
        case 0 : histyy[ij] = (TH2F*)h_xcorstrips[ij]->Clone(); break;
        case 1 : histyy[ij] = (TH2F*)h_ycorstrips[ij]->Clone(); break;
        case 2 : histyy[ij] = (TH2F*)h_xycorstrips[ij]->Clone(); break;
        case 3 : histyy[ij] = (TH2F*)h_yxcorstrips[ij]->Clone(); break;
        case 4 : histyy[ij] = (TH2F*)h_xtcorstrips[ij]->Clone(); break;
        case 5 : histyy[ij] = (TH2F*)h_ytcorstrips[ij]->Clone(); break;
        case 6 : histyy[ij] = (TH2F*)h_xytcorstrips[ij]->Clone(); break;
        case 7 : histyy[ij] = (TH2F*)h_yxtcorstrips[ij]->Clone(); break;	     
        defalt : histyy[ij] = (TH2F*)h_xcorstrips[ij]->Clone(); break;
        }
        if (histyy[ij]->GetMaximum()>histyy[ij]->GetMinimum()+0.01) { 
          histyy[ij]->SetMaximum(min(1.0, histyy[ij]->GetMaximum()));
          histyy[ij]->SetMinimum(max(1.e-3, histyy[ij]->GetMinimum()));
        }
        histyy[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
        histyy[ij]->GetYaxis()->SetRangeUser(-0.5, nusedstrip);
        
        histyy[ij]->GetXaxis()->SetTitle("Extrapolated Strip");
        histyy[ij]->GetXaxis()->CenterTitle();
        histyy[ij]->GetXaxis()->SetTitleOffset(0.8);
        histyy[ij]->GetXaxis()->SetTitleSize(0.07); 
        
        histyy[ij]->GetYaxis()->SetTitle("Observed Strip");
        histyy[ij]->GetYaxis()->CenterTitle();
        histyy[ij]->GetYaxis()->SetTitleOffset(0.7);
        histyy[ij]->GetYaxis()->SetTitleSize(0.07); 
        
        histyy[ij]->GetZaxis()->SetLabelSize(.06);
        histyy[ij]->Draw("colz");
      }
      c4a->Update();
      for (int ij=0; ij<nlayer; ij++) {
        if (histyy[ij]) { delete histyy[ij]; histyy[ij]=0;}
      }
    }
  }
  
  ps.NewPage();
  gStyle->SetOptLogz(0);
  gStyle->SetPadTopMargin(0.09); //0.01);
  gStyle->SetOptStat(0);
  
  TCanvas* c6 = new TCanvas("c6", "c6", 700, 900);
  c6->Divide(2,4);    
  
  //  ps.NewPage();
  c6->cd(1); xlayer_alloccu->SetLineColor(1);    xlayer_alloccu->Draw(); 
  c6->cd(2); ylayer_alloccu->SetLineColor(1);    ylayer_alloccu->Draw(); 
  double entry = max(1.,h_xrawcorhits->GetBinContent(0,0));
  c6->cd(3);  h_xrawcorhits->Scale(1./entry);  h_xrawcorhits->Draw("colz");
  c6->cd(5);  h_yrawcorhits->Scale(1./entry);  h_yrawcorhits->Draw("colz");
  c6->cd(7);  h_xyrawcorhits->Scale(1./entry); h_xyrawcorhits->Draw("colz");
  
  if (isTiming) {
    c6->cd(4);  h_xtrawcorhits->Scale(1./entry);  h_xtrawcorhits->Draw("colz");
    c6->cd(6);  h_ytrawcorhits->Scale(1./entry);  h_ytrawcorhits->Draw("colz");
    c6->cd(8);  h_xytrawcorhits->Scale(1./entry); h_xytrawcorhits->Draw("colz");
  }
  
  c6->Update();
  
  ps.NewPage();
  c6->cd(1); xlayer_alloccu->SetLineColor(1);    xlayer_alloccu->Draw(); 
  xlayer_alloccusel->SetLineColor(2); xlayer_alloccusel->Draw("same");
  c6->cd(2); ylayer_alloccu->SetLineColor(1);    ylayer_alloccu->Draw(); 
  ylayer_alloccusel->SetLineColor(2); ylayer_alloccusel->Draw("same");
  c6->cd(3);  h_xcorhits->Scale(1./nTotalp);  h_xcorhits->Draw("colz");
  c6->cd(5);  h_ycorhits->Scale(1./nTotalp);  h_ycorhits->Draw("colz");
  c6->cd(7);  h_xycorhits->Scale(1./nTotalp); h_xycorhits->Draw("colz");
  
  if (isTiming) {
    c6->cd(4);  h_xtcorhits->Scale(1./nTotalt);  h_xtcorhits->Draw("colz");
    c6->cd(6);  h_ytcorhits->Scale(1./nTotalt);  h_ytcorhits->Draw("colz");
    c6->cd(8);  h_xytcorhits->Scale(1./nTotalt); h_xytcorhits->Draw("colz");
  }
  c6->Update();
  
  ps.NewPage();
  if (isalign>0) {
    c6->cd(1); shift_pos->Draw("colz");
    c6->cd(2); rms_pos->Draw("colz");
  }
  c6->cd(3); h_deltaposxcov->Draw("colz"); 
  c6->cd(4); h_deltaposycov->Draw("colz"); 
  
  c6->Update();
  
  if (isTiming) { 
    ps.NewPage();
    c6->cd(1); time_mean_reso->Draw("colz");
    c6->cd(2); time_rms_reso->Draw("colz");
    c6->cd(3); time_exterr_reso->Draw("colz");
    c6->cd(4); time_corrms_reso->Draw("colz");
    c6->cd(5); h_deltatcov->Draw("colz"); 
    c6->cd(6); h_deltatcov2->Draw("colz");
    c6->cd(7); h_deltatcovy->Draw("colz");
    
    c6->Update();
  }
  TH2F* hist_x_lay = new TH2F("hist_x_lay", "hist_x_lay", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  TH2F* hist_y_lay = new TH2F("hist_y_lay", "hist_y_lay", nlayer, -0.5, nlayer-0.5, nlayer, -0.5, nlayer-0.5);
  
  h_posxcoreff->Scale(1./nTotallx);
  h_posycoreff->Scale(1./nTotally);
  
  file_out<<"Xlay_Corr_eff ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<h_posxcoreff->GetBinContent(ix+1,iy+1)<<", \t";
      // hist_x_lay->Fill(ix,iy,h_posxcoreff->GetBinContent(ix,iy)/max(1.,posxcoreffCount));
    }
    file_out<<endl;
  }
  file_out<<"}"<<endl;
  
  file_out<<"Ylay_Corr_eff ["<<nlayer<<"]["<<nlayer<<"] ={"<<endl;
  for(int ix=0;ix<nlayer;ix++){
    for(int iy=0;iy<nlayer;iy++){
      file_out<<h_posycoreff->GetBinContent(ix+1,iy+1)<<", \t";
      // hist_y_lay->Fill(ix,iy,h_posycoreff->GetBinContent(ix,iy)/max(1.,posycoreffCount));
    }
    file_out<<endl;
  }
  file_out<<"}"<<endl;
  
  ps.NewPage();
  gStyle->SetOptLogz(0);
  gStyle->SetPadTopMargin(0.09); //0.01);
  gStyle->SetOptStat(0);
  TCanvas* c16 = new TCanvas("c16","c16",700,900);
  c16->Divide(1,2);
  c16->cd(1);
  h_posxcoreff->Draw("colz");
  h_posxcoreff->SetTitle("Cor_lay_eff_X");
  c16->Update();
  c16->cd(2);
  h_posycoreff->Draw("colz");
  h_posycoreff->SetTitle("Cor_lay_eff_Y");
  c16->Update();
  
  h_posxcoreff_dif->Add(h_posxcoreff,h_posxcoreff_def,1.,-1.);
  h_posycoreff_dif->Add(h_posycoreff,h_posycoreff_def,1.,-1.);
  ps.NewPage();
  gStyle->SetOptLogz(0);
  gStyle->SetPadTopMargin(0.09); //0.01);
  gStyle->SetOptStat(0);
  TCanvas* c17 = new TCanvas("c17","c17",700,900);
  c17->Divide(1,2);
  c17->cd(1);
  h_posxcoreff_dif->Draw("colz");
  h_posxcoreff_dif->SetTitle("Cor_lay_eff_X_dif");
  c17->Update();
  c17->cd(2);
  h_posycoreff_dif->Draw("colz");
  h_posycoreff_dif->SetTitle("Cor_lay_eff_Y_dif");
  c17->Update();
  
  if (isTiming) { 
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetPadTopMargin(0.10); 
    gStyle->SetPadBottomMargin(0.07);
    
    if (isalign>0) { 
      ps.NewPage();
      gStyle->SetPadLeftMargin(1.0);
      TCanvas* c5 = new TCanvas("c5", "c5", 700, 900);
      c5->Divide(3,4);
      
      for (int kl=0; kl<=nmxiter; kl++) {
        c5->cd(kl+1); 
        time_both_cormean[kl]->GetXaxis()->SetLabelSize(0.08);
        time_both_cormean[kl]->GetYaxis()->SetLabelSize(0.08);
        time_both_cormean[kl]->GetZaxis()->SetLabelSize(0.065);
        time_both_cormean[kl]->SetMaximum(0.5); //min(2.,time_both_cormean[kl]->GetMaximum()));
        time_both_cormean[kl]->SetMinimum(-0.5); //max(-2.,time_both_cormean[kl]->GetMinimum()));
        time_both_cormean[kl]->GetXaxis()->SetLabelOffset(0.001);
        time_both_cormean[kl]->Draw("colz");
      }
      c5->Update();
      
      ps.NewPage();
      for (int kl=0; kl<=nmxiter; kl++) {
        c5->cd(kl+1); 
        time_both_corrms[kl]->GetXaxis()->SetLabelSize(0.08);
        time_both_corrms[kl]->GetYaxis()->SetLabelSize(0.08);
        time_both_corrms[kl]->GetZaxis()->SetLabelSize(0.08);
        time_both_corrms[kl]->SetMaximum(2.5); //min(3.5,time_both_corrms[kl]->GetMaximum()));
        time_both_corrms[kl]->SetMinimum(0.5); //max(0.5,time_both_corrms[kl]->GetMinimum()));
        time_both_corrms[kl]->GetXaxis()->SetLabelOffset(0.001);
        time_both_corrms[kl]->Draw("colz");
      }
      c5->Update();
      
      ps.NewPage();
      gStyle->SetPadBottomMargin(0.07);
      TCanvas* c5b = new TCanvas("c5b", "c5b", 700, 900);
      c5b->Divide(4,5);
      for (int jk=0; jk<=npixel; jk++) { //150123
        c5b->cd(jk+1);
        timexy_cormean[jk]->GetXaxis()->SetLabelSize(0.08);
        timexy_cormean[jk]->GetYaxis()->SetLabelSize(0.08);
        
        timexy_cormean[jk]->GetZaxis()->SetLabelSize(0.08);
        if (timexy_cormean[jk]->GetMaximum()>timexy_cormean[jk]->GetMinimum()+0.1) { 
          timexy_cormean[jk]->SetMaximum(min(2.5,timexy_cormean[jk]->GetMaximum()));
          timexy_cormean[jk]->SetMinimum(max(-2.5,timexy_cormean[jk]->GetMinimum()));
        }
        timexy_cormean[jk]->Draw("colz");
      }
      
      c5b->Update();
      ps.NewPage();
      for (int jk=0; jk<=npixel; jk++) {
        c5b->cd(jk+1);
        timexy_corrms[jk]->GetXaxis()->SetLabelSize(0.08);
        timexy_corrms[jk]->GetYaxis()->SetLabelSize(0.08);
        timexy_corrms[jk]->GetZaxis()->SetLabelSize(0.08);
        timexy_corrms[jk]->SetMaximum(2.0); //min(2.4,timexy_corrms[jk]->GetMaximum()));
        timexy_corrms[jk]->SetMinimum(0.3); //max(0.4,timexy_corrms[jk]->GetMinimum()));
        timexy_corrms[jk]->Draw("colz");
      }
      
      c5b->Update();
      double amx[8] ={10.0, 0.6, 1.0, 3.0, 1.5, 6.0, 4.5, 3.5};
      double amn[8]={-10.0, -0.6, -1.0, 0.5, -1.0, -2.0, 0.0, 0.5};
      for (int ixx=0; ixx<8; ixx++) { 
        TH1F* tmphist[nmxiter];
        ps.NewPage();
        for (int jk=0; jk<nmxiter; jk++) {
          c5->cd(jk+1);
          switch(ixx) { 
          case 0 : tmphist[jk] = (TH1F*)time_offset[jk]->Clone(); break; 
          case 1 : tmphist[jk] = (TH1F*)shift_time_mnft[jk]->Clone(); break;
          case 2 : tmphist[jk] = (TH1F*)statmean_time[jk]->Clone(); break;
          case 3 : tmphist[jk] = (TH1F*)statrms_time[jk]->Clone(); break;
          case 4 : tmphist[jk] = (TH1F*)statskew_time[jk]->Clone(); break;
          case 5 : tmphist[jk] = (TH1F*)statkurt_time[jk]->Clone(); break;
          case 6 : tmphist[jk] = (TH1F*)rms_time[jk]->Clone(); break;
          case 7 : tmphist[jk] = (TH1F*)rms_timeused[jk]->Clone(); break;
          }
          
          if (jk>0 || tmphist[jk]->GetMaximum() >tmphist[jk]->GetMinimum()+0.1) { 
            tmphist[jk]->SetMaximum(min(amx[ixx], tmphist[jk]->GetMaximum()));
            tmphist[jk]->SetMinimum(max(amn[ixx], tmphist[jk]->GetMinimum()));
          }
          tmphist[jk]->GetXaxis()->SetLabelSize(0.08);
          tmphist[jk]->GetYaxis()->SetLabelSize(0.08);
          tmphist[jk]->GetXaxis()->SetLabelOffset(0.001);
          tmphist[jk]->Draw();
          if (ixx==6) {rms_timeused[jk]->Draw("same");}
        }
        c5->Update();
        for (int jk=0; jk<nmxiter; jk++) {
          if (tmphist[jk]) { delete tmphist[jk]; tmphist[jk]=0;}
        }
      }
      
      ps.NewPage();
      gStyle->SetOptLogz(0);
      gStyle->SetTitleFontSize(0.07);
      gStyle->SetPadLeftMargin(0.09);
      gStyle->SetPadBottomMargin(0.09);
      gStyle->SetPadTopMargin(0.09); //0.03);
      
      TCanvas* c4b = new TCanvas("c4b", "c4b", 700, 900);
      c4b->Divide(3,4);
      for (int jkl=0; jkl<4; jkl++) {
        if (jkl>0) { ps.NewPage();}
        TH2F* tmp2hist[nlayer];
        for (int ij=0; ij<nlayer; ij++) {
          switch(jkl) {
          case 0 : tmp2hist[ij] = (TH2F*)correction_xtime[ij]->Clone(); break;
          case 1 : tmp2hist[ij] = (TH2F*)fitted_rms_xtime[ij]->Clone(); break;
          case 2 : tmp2hist[ij] = (TH2F*)correction_ytime[ij]->Clone(); break;
          case 3 : tmp2hist[ij] = (TH2F*)fitted_rms_ytime[ij]->Clone(); break;
          }
          c4b->cd(ij+1);
          tmp2hist[ij]->GetXaxis()->SetLabelSize(0.07);
          tmp2hist[ij]->GetYaxis()->SetLabelSize(0.07);
          tmp2hist[ij]->GetZaxis()->SetLabelSize(0.06);
          if (jkl==0 || jkl==2) { 
            tmp2hist[ij]->SetMaximum(0.50); 
            tmp2hist[ij]->SetMinimum(-0.50);
          } else {
            tmp2hist[ij]->SetMaximum(4.00);
            tmp2hist[ij]->SetMinimum(0.50);
          }
          
          tmp2hist[ij]->GetXaxis()->SetTitle("Strip No");
          tmp2hist[ij]->GetXaxis()->CenterTitle();
          tmp2hist[ij]->GetXaxis()->SetTitleOffset(0.8);
          tmp2hist[ij]->GetXaxis()->SetTitleSize(0.08); 
          
          tmp2hist[ij]->GetYaxis()->SetTitle("# of iteration");
          tmp2hist[ij]->GetYaxis()->CenterTitle();
          tmp2hist[ij]->GetYaxis()->SetTitleOffset(0.5);
          tmp2hist[ij]->GetYaxis()->SetTitleSize(0.08); 
          
          tmp2hist[ij]->GetXaxis()->SetRangeUser(-0.5, nusedstrip);
          tmp2hist[ij]->Draw("colz");
        }
        c4b->Update();
        for (int ij=0; ij<nlayer; ij++) {
          if (tmp2hist[ij]) { delete tmp2hist[ij]; tmp2hist[ij]=0;}
        }
      }
      
      ps.NewPage();
      gStyle->SetOptLogy(0);
      gStyle->SetOptStat(1110);
      gStyle->SetStatW(0.36);
      gStyle->SetStatH(0.28); //30);
      gStyle->SetPadRightMargin(0.02);
      gStyle->SetTitleFontSize(0.07);
      gStyle->SetPadTopMargin(0.09);
      gStyle->SetPadBottomMargin(0.07);
      
      TCanvas* c5t = new TCanvas("c5t", "c5t", 700, 900);
      Time_diff_evt->Draw();
      Time_diff_evt->SetTitle("Delta_time_evt");
      Time_diff_evt->GetXaxis()->SetTitle("#Delta t");
      Time_diff_evt->GetXaxis()->CenterTitle();
      c5t->Update();
      
      ps.NewPage();
      gStyle->SetOptLogy(0);
      gStyle->SetOptStat(1110);
      gStyle->SetStatW(0.36);
      gStyle->SetStatH(0.28); //30);
      gStyle->SetPadRightMargin(0.02);
      gStyle->SetTitleFontSize(0.07);
      gStyle->SetPadTopMargin(0.09);
      gStyle->SetPadBottomMargin(0.07);
      
      TCanvas* c5a = new TCanvas("c5a", "c5a", 700, 900);
      c5a->Divide(3,4);
      
      TH1F* histx[nmxiter];
      TH1F* histy[nmxiter];
      
      for (int jkl=0; jkl<8; jkl++) {
        ps.NewPage();
        for (int jk=0; jk<nmxiter; jk++) {
          switch(jkl) {
          case 0 : histx[jk] = (TH1F*)time_offsetx[jk]->Clone();
            histy[jk] = (TH1F*)time_offsety[jk]->Clone(); break;
            
          case 1 : histx[jk] = (TH1F*)shift_time_mnftx[jk]->Clone();
            histy[jk] = (TH1F*)shift_time_mnfty[jk]->Clone(); break;
            
          case 2 : histx[jk] = (TH1F*)statmean_timex[jk]->Clone();
            histy[jk] = (TH1F*)statmean_timey[jk]->Clone(); break;
            
          case 3 : histx[jk] = (TH1F*)statrms_timex[jk]->Clone();
            histy[jk] = (TH1F*)statrms_timey[jk]->Clone(); break;
            
          case 4 : histx[jk] = (TH1F*)statskew_timex[jk]->Clone();
            histy[jk] = (TH1F*)statskew_timey[jk]->Clone(); break;
            
          case 5 : histx[jk] = (TH1F*)statkurt_timex[jk]->Clone();
            histy[jk] = (TH1F*)statkurt_timey[jk]->Clone(); break;
            
          case 6 : histx[jk] = (TH1F*)rms_timex[jk]->Clone();
            histy[jk] = (TH1F*)rms_timey[jk]->Clone(); break;
            
          case 7 : histx[jk] = (TH1F*)rms_timeusedx[jk]->Clone();
            histy[jk] = (TH1F*)rms_timeusedy[jk]->Clone(); break;
            
          default : histx[jk] = (TH1F*)statkurt_time[jk]->Clone();
            histy[jk] = (TH1F*)statkurt_timey[jk]->Clone(); break;
          }
        }
        gStyle->SetStatY(.99); gStyle->SetStatTextColor(1);
        for (int jk=0; jk<nmxiter; jk++) {
          c5a->cd(jk+1);
          histx[jk]->GetXaxis()->SetLabelSize(0.065);
          histx[jk]->GetYaxis()->SetLabelSize(0.065);
          histx[jk]->SetLineColor(1); histx[jk]->Draw();
        }
        c5a->Update(); 
        
        gStyle->SetStatY(0.77); gStyle->SetStatTextColor(2); 
        for (int jk=0; jk<nmxiter; jk++) {
          c5a->cd(jk+1);
          histy[jk]->SetLineColor(2); histy[jk]->Draw("sames");
        }
        
        c5a->Update();
        
        for (int jk=0; jk<nmxiter; jk++) {
          if (histx[jk]) { delete histx[jk]; histx[jk]=0;}
          if (histy[jk]) { delete histy[jk]; histy[jk]=0;}
        }
      }
      ps.NewPage();
    }  //if (isalign>0) 
  } // if (isTiming)
  
  ps.Close();
  
  //
  
  fileOut->Write();
  file_out.close();
  file_outstr.close();
  fileOut->Close();
  // filecorOut->Close();
  
  // filecorOut->cd();
  // T3->Write();
  // filecorOut->Close();
  return 0;
  //--------------------------------------------------------------
}
    

/*

*/
