#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;
#include "StraightLineFit.h"
const double cvalx = 29.979; // velocity of light in cm/ns
unsigned int mypow_2[32] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048}; //, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648};

StraightLineFit::StraightLineFit() {

  for (int ij=0; ij<nlayerx; ij++) {
    xval[ij] = yval[ij] = yerr[ij] = -100;
    ydev[ij] = yexp[ij] = yexperr[ij] =1000;
    isUsed[ij] = false;
  }
  isPosition=1;
  chiSquare = -1;
  isFailed=true;
  nused = occulyr = occulyr2 = layfirst = laylast = -1;
  allowedDev = intersect = slope = errlin = errcst = errcov = -1000.0;
  allLayerId = 0;
}


StraightLineFit::StraightLineFit(int type, double* xv, double* yv, double* ye, bool* used, int occu, int occu2, int first, int last, float mxdev) {

  isPosition = type;
  for (int ij=0; ij<nlayerx; ij++) {
    isUsed[ij] = used[ij];
    xval[ij] = xv[ij];
    yval[ij] = yv[ij]; 
    yerr[ij] = ye[ij]; 
    ydev[ij] = yexp[ij] = yexperr[ij] = 1000;
  }
  occulyr = occu;
  occulyr2 = occu2;
  layfirst = first;
  laylast = last;
  allowedDev = mxdev;
  allLayerId = 0;
  chiSquare = -1;
  isFailed=false;
  nused = -1;
  intersect = slope = errlin = errcst = errcov = -1000;
}

// xco = a zco + c
// szxy=0; // xco = a zco + c
// sz=0; //   z-layer
// sxy=0; //   x/y strip
// sn=0;  //Number
// sz2=0; //z-layer^2

void StraightLineFit::GetParameters(int& failed, double& cst, double& slp) {

  isFailed = false;
  bool layerP[nlayerx]; //Only layers, which are not used in fit, but want to calculate dev for it
  for (int ij=0; ij<nlayerx; ij++) { layerP[ij] = (yval[ij]>-50) ? true : false;}
  double tmpslope=0;
  for (int kl=0; kl<6; kl++) { //Four iterations to filter out noise hits
    nused = 0;
    chiSquare = 0.0;
    double szxy=0, sz=0, sxy=0, sn=0, sz2=0;
    tmpslope=0;
    for (int ix=0; ix<nlayerx; ix++) { ydev[ix] = 1000;}
    
    for(int ij=0;ij<nlayerx;ij++){
      if (ij>=layfirst && ij<=laylast && ij!=occulyr && ij!=occulyr2) {
        if (yval[ij]>-50 && isUsed[ij] && layerP[ij]) { 
          szxy +=xval[ij]*yval[ij]/yerr[ij];
          sz +=xval[ij]/yerr[ij];
          sz2 +=xval[ij]*xval[ij]/yerr[ij];
          sxy +=yval[ij]/yerr[ij];
          sn +=1/yerr[ij];
          nused++;	
          allLayerId +=mypow_2[ij];
          
          //  if (isPosition<=0) {
          //  cout<<"ij "<<ij<<" "<<xval[ij]<<" "<<yerr[ij]<<" "<<szxy<<" "<<sz<<" "<<sxy<<" "<<sn<< ""<<nused<<" "<<allLayerId<<endl;
          // }
        }
      }
    }
    
    //    cout<<"noused "<<nused<<" "<<layfirst<<" "<<laylast<<endl;
    if (sn >0 && sz2*sn - sz*sz !=0) { 
      slope = (szxy*sn - sz*sxy)/(sz2*sn - sz*sz);
      
      tmpslope = slope; 
      //      if (isPosition==0 && occulyr<nlayerx && occulyr2<nlayerx) { //time offset correction
      if (isPosition==0) { //time offset correction
        if (abs(cvalx*slope+1)<3.30) { 
          tmpslope = -1/cvalx;
        } else {
          nused = 0;
        }
      }	
      
      intersect = sxy/sn - tmpslope*sz/sn;  // <xy> = a + b <z>
      
      double determ = (sn*sz2 - sz*sz);
      errcst = sz2/determ;
      errcov = -sz/determ;
      errlin = sn/determ;
      slope2 = szxy/sz2;
    } else { 
      tmpslope = slope = slope2 = -999; intersect= -999;
      errcst =-999;
      errcov=-999; errlin=-999;
      isFailed = true;
    }

    ///////////////////////////////////////
    // If fit fails go out from the loop //
    ///////////////////////////////////////

    if (slope== -999 || nused<2) break;
      
    isFailed = false;
    for(int ij=0;ij<nlayerx;ij++){
      yexp[ij]=tmpslope*xval[ij]+intersect;
      
      //      cout <<int(isPosition)<<" ij "<< ij <<" "<<xval[ij]<<" "<<yval[ij]<<" "<<yexp[ij]<<" "<<tmpslope<<" "<<intersect<<endl;

      yexperr[ij] = errcst + 2*errcov*xval[ij] + errlin*xval[ij]*xval[ij];
      if (yexperr[ij]<0) {
	cout <<"Syexperr "<< ij<<" "<<nused<<" "<<sz2*sn - sz*sz<<" "<<sn<<" "<<yexperr[ij]<<" "<< errcst + 2*errcov*xval[ij] + errlin*xval[ij]*xval[ij]<<" "<<errcst<<" "<<errcov<<" "<<errlin<<endl;
	yexperr[ij] = sqrt(-yexperr[ij]);
      } else {
	yexperr[ij] = sqrt(yexperr[ij]);
      }

      if (layerP[ij]) { 
	ydev[ij] = (yval[ij]-yexp[ij]);
	if (yval[ij]>-50 && isUsed[ij]) { 

	  if (ij>=layfirst && ij<=laylast && ij!=occulyr && ij!=occulyr2) {
	    //	  if ((isPosition==1 && abs(ydev[ij])>allowedDev) ||  //Position in absolute scale
	    //	      (isPosition<=0 && abs(ydev[ij]/sqrt(yerr[ij]))>allowedDev)) { //Time interms of resolution
	    if (abs(ydev[ij])>allowedDev) { 
	      layerP[ij] = false; //Only for layer, whaich are not used in fit, but want to calculate dev
	      isFailed = true; // It is kept because if any event also rejected at least we will know and it takes care different filling accordingly
	      yval[ij] = -90.0;
	    }
	    
	    chiSquare +=  pow( ydev[ij], 2.)/yerr[ij]; 
	  }
	}
      }
    }
    if (!isFailed) break;
  } //  for (int kl=0; kl<2; kl++) 

  failed = isFailed;
  cst = intersect;
  slp = slope;

}

void StraightLineFit::GetError(double& lerr, double& slperr, double& cov) {
  lerr = errcst;
  slperr = errlin;
  cov = errcov;
}

void StraightLineFit::GetChisqure(int& ndof, double& chis) {
  ndof = nused;
  chis = chiSquare;
}

void StraightLineFit::GetFitValues(double* exp, double* dev, double* experr) {
  for (int ix=0; ix<nlayerx; ix++) {
    exp[ix] = yexp[ix];
    dev[ix] = ydev[ix];
    experr[ix] = yexperr[ix];
  }
}
