//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep  3 09:43:07 2012 by ROOT version 5.32/03
// from TTree evetree/event tree
// found on file: /home/nitali/Programming/DATAFILE/NITALI_MUON_LT_20120505_004556.ire
//////////////////////////////////////////////////////////

#ifndef EveTree_h
#define EveTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TTimeStamp.h>
#include <TBits.h>


#define NL 12
#define NS 32
// Fixed size dimensions of array or collections stored in the TTree if any.
class EveTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TTimeStamp      *Evetime;
   Int_t           ENum;
   TBits           *xstriphitsL0;
   TBits           *ystriphitsL0;
   TBits           *xstriphitsL1;
   TBits           *ystriphitsL1;
   TBits           *xstriphitsL2;
   TBits           *ystriphitsL2;
   TBits           *xstriphitsL3;
   TBits           *ystriphitsL3;
   TBits           *xstriphitsL4;
   TBits           *ystriphitsL4;
   TBits           *xstriphitsL5;
   TBits           *ystriphitsL5;
   TBits           *xstriphitsL6;
   TBits           *ystriphitsL6;
   TBits           *xstriphitsL7;
   TBits           *ystriphitsL7;
   TBits           *xstriphitsL8;
   TBits           *ystriphitsL8;
   TBits           *xstriphitsL9;
   TBits           *ystriphitsL9;
   TBits           *xstriphitsL10;
   TBits           *ystriphitsL10;
   TBits           *xstriphitsL11;
   TBits           *ystriphitsL11;
   ULong64_t       tdcdata[32];
//   ULong64_t       tdcdata1[32];
//   ULong64_t       tdcdata2[32];
   Int_t           tdcmult[32];
   ULong64_t       tdcref;

   // List of branches
   TBranch        *b_Evetime;   //!
   TBranch        *b_ENum;   //!
   TBranch        *b_xstriphitsL0;   //!
   TBranch        *b_ystriphitsL0;   //!
   TBranch        *b_xstriphitsL1;   //!
   TBranch        *b_ystriphitsL1;   //!
   TBranch        *b_xstriphitsL2;   //!
   TBranch        *b_ystriphitsL2;   //!
   TBranch        *b_xstriphitsL3;   //!
   TBranch        *b_ystriphitsL3;   //!
   TBranch        *b_xstriphitsL4;   //!
   TBranch        *b_ystriphitsL4;   //!
   TBranch        *b_xstriphitsL5;   //!
   TBranch        *b_ystriphitsL5;   //!
   TBranch        *b_xstriphitsL6;   //!
   TBranch        *b_ystriphitsL6;   //!
   TBranch        *b_xstriphitsL7;   //!
   TBranch        *b_ystriphitsL7;   //!
   TBranch        *b_xstriphitsL8;   //!
   TBranch        *b_ystriphitsL8;   //!
   TBranch        *b_xstriphitsL9;   //!
   TBranch        *b_ystriphitsL9;   //!
   TBranch        *b_xstriphitsL10;   //!
   TBranch        *b_ystriphitsL10;   //!
   TBranch        *b_xstriphitsL11;   //!
   TBranch        *b_ystriphitsL11;   //!
   TBranch        *b_TDCdata;   //!
//   TBranch        *b_TDCdata1;   //!
//   TBranch        *b_TDCdata2;   //!
   TBranch        *b_TDCmult;   //!
   TBranch        *b_TDCref;   //!

   EveTree(TTree *tree=0);
   virtual ~EveTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TBits *xLayer[NL];
   TBits *yLayer[NL];
};

#endif

#ifdef EveTree_cxx
EveTree::EveTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/media/34A4E3FAA4E3BD0C/TESTRUN_4MERG_TOPMIDDLEBOTTOM.ire");
      if (!f || !f->IsOpen()) {
         f = new TFile("/media/34A4E3FAA4E3BD0C/TESTRUN_4MERG_TOPMIDDLEBOTTOM.ire");
      }
      f->GetObject("evetree",tree);

   }
   Init(tree);
}

//NITALI_MUON_LT_20120507_191134
//TESTRUN_BAKELITE_7MERG_TOPMIDDLEBOTTOM
//TESTRUN_4MERG_TOPMIDDLEBOTTOM
EveTree::~EveTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EveTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EveTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EveTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Evetime = 0;
   xstriphitsL0 = 0;
   ystriphitsL0 = 0;
   xstriphitsL1 = 0;
   ystriphitsL1 = 0;
   xstriphitsL2 = 0;
   ystriphitsL2 = 0;
   xstriphitsL3 = 0;
   ystriphitsL3 = 0;
   xstriphitsL4 = 0;
   ystriphitsL4 = 0;
   xstriphitsL5 = 0;
   ystriphitsL5 = 0;
   xstriphitsL6 = 0;
   ystriphitsL6 = 0;
   xstriphitsL7 = 0;
   ystriphitsL7 = 0;
   xstriphitsL8 = 0;
   ystriphitsL8 = 0;
   xstriphitsL9 = 0;
   ystriphitsL9 = 0;
   xstriphitsL10 = 0;
   ystriphitsL10 = 0;
   xstriphitsL11 = 0;
   ystriphitsL11 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Evetime", &Evetime, &b_Evetime);
   fChain->SetBranchAddress("ENum", &ENum, &b_ENum);
   fChain->SetBranchAddress("xstriphitsL0", &xstriphitsL0, &b_xstriphitsL0);
   fChain->SetBranchAddress("ystriphitsL0", &ystriphitsL0, &b_ystriphitsL0);
   fChain->SetBranchAddress("xstriphitsL1", &xstriphitsL1, &b_xstriphitsL1);
   fChain->SetBranchAddress("ystriphitsL1", &ystriphitsL1, &b_ystriphitsL1);
   fChain->SetBranchAddress("xstriphitsL2", &xstriphitsL2, &b_xstriphitsL2);
   fChain->SetBranchAddress("ystriphitsL2", &ystriphitsL2, &b_ystriphitsL2);
   fChain->SetBranchAddress("xstriphitsL3", &xstriphitsL3, &b_xstriphitsL3);
   fChain->SetBranchAddress("ystriphitsL3", &ystriphitsL3, &b_ystriphitsL3);
   fChain->SetBranchAddress("xstriphitsL4", &xstriphitsL4, &b_xstriphitsL4);
   fChain->SetBranchAddress("ystriphitsL4", &ystriphitsL4, &b_ystriphitsL4);
   fChain->SetBranchAddress("xstriphitsL5", &xstriphitsL5, &b_xstriphitsL5);
   fChain->SetBranchAddress("ystriphitsL5", &ystriphitsL5, &b_ystriphitsL5);
   fChain->SetBranchAddress("xstriphitsL6", &xstriphitsL6, &b_xstriphitsL6);
   fChain->SetBranchAddress("ystriphitsL6", &ystriphitsL6, &b_ystriphitsL6);
   fChain->SetBranchAddress("xstriphitsL7", &xstriphitsL7, &b_xstriphitsL7);
   fChain->SetBranchAddress("ystriphitsL7", &ystriphitsL7, &b_ystriphitsL7);
   fChain->SetBranchAddress("xstriphitsL8", &xstriphitsL8, &b_xstriphitsL8);
   fChain->SetBranchAddress("ystriphitsL8", &ystriphitsL8, &b_ystriphitsL8);
   fChain->SetBranchAddress("xstriphitsL9", &xstriphitsL9, &b_xstriphitsL9);
   fChain->SetBranchAddress("ystriphitsL9", &ystriphitsL9, &b_ystriphitsL9);
   fChain->SetBranchAddress("xstriphitsL10", &xstriphitsL10, &b_xstriphitsL10);
   fChain->SetBranchAddress("ystriphitsL10", &ystriphitsL10, &b_ystriphitsL10);
   fChain->SetBranchAddress("xstriphitsL11", &xstriphitsL11, &b_xstriphitsL11);
   fChain->SetBranchAddress("ystriphitsL11", &ystriphitsL11, &b_ystriphitsL11);
   fChain->SetBranchAddress("tdcdata", tdcdata, &b_TDCdata);
//   fChain->SetBranchAddress("tdcdata1", tdcdata1, &b_TDCdata1);
//   fChain->SetBranchAddress("tdcdata2", tdcdata2, &b_TDCdata2);
   fChain->SetBranchAddress("tdcmult", tdcmult, &b_TDCmult);
   fChain->SetBranchAddress("tdcref", &tdcref, &b_TDCref);
   Notify();
}

Bool_t EveTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EveTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EveTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EveTree_cxx
