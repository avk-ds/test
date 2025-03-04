#ifndef STRAIGHTLINEFIT_H
#define STRAIGHTLINEFIT_H

//for time fit include usetime etc during the calcualtion of xtime etc.
// tmpxent & tmpyent : use the same criteria in main code, instead of carried it out from this class or use dev for this;

const int nlayerx=12;
const int nstripmx=32;
class StraightLineFit {
 public :
  StraightLineFit();
  StraightLineFit(int type, double* xv, double* yv, double* ye, bool* used, int occu, int occu2, int first, int last, float mxdev);
  void GetParameters(int& isfailed, double& a, double& b);
  void GetError(double& lerr, double& slperr, double& cov);
  void GetChisqure(int& ndof, double& chis);
  void GetFitValues(double* exp, double* dev, double* experr);
  int GetLayerIds() {return allLayerId ;}
  double GetSlope2() {return slope2 ;}
  ~StraightLineFit(){;};

 private :
  int occulyr;  //Layer which is under study
  int occulyr2;  //second layer for correlation study
  int layfirst; //Lowest layer
  int laylast;  //Highest layer
  double xval[nlayerx]; //x-valus
  double yval[nlayerx]; //Y-values
  double yerr[nlayerx]; //error square of Y
  bool   isUsed[nlayerx]; // Layer will be used in fit or not
  double yexp[nlayerx]; //Expected Y value after fit in each layer
  double ydev[nlayerx]; //Deviation from measurement
  double yexperr[nlayerx]; //Error on expected postion (from fit)
  
  float allowedDev;  //Allowed deviation to be included in fit
  double chiSquare;   // chi^2
  int nused;         // Total number of layer used in fit
  bool isFailed;     //Fit is failed or not
  int isPosition;   // Postion fit or time fit 
                    // 1 : position
                    // 0 : time fit and constraining slope to 1/c for correction
                    // -1 : No constraining
  double intersect;   // Constant fit para
  double   slope;     // Slope of the fit
  double   slope2;     // Slope of only time fit with the assumption that constant term is zero
  double errlin;     // Error in slope
  double errcst;     // Error in constant
  double errcov;     //Covariance term
  int allLayerId;   //Bitwise information of all layer used in fit
};
#endif
