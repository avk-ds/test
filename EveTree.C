#define EveTree_cxx
#include "EveTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2F.h"
#include "THStack.h"
#include <TApplication.h>


#include <iostream>
#include <cmath>
#include <stdlib.h>

using namespace std;

void EveTree::Loop()
{
    //   In a ROOT session, you can do:
    //      Root > .L EveTree.C
    //      Root > EveTree t
    //      Root > t.GetEntry(12); // Fill t data members with entry number 12
    //      Root > t.Show();       // Show values of entry 12
    //      Root > t.Show(16);     // Read and show values of entry 16
    //      Root > t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch
//    TApplication app("my_app",0,0);
//    TCanvas *c = new TCanvas("myCanvas","My Canvas");
//    TCanvas *c1 = new TCanvas("myCanvas1","My Canvas1");

    xLayer[0]=xstriphitsL0;
    yLayer[0]=ystriphitsL0;
    xLayer[1]=xstriphitsL1;
    yLayer[1]=ystriphitsL1;
    xLayer[2]=xstriphitsL2;
    yLayer[2]=ystriphitsL2;
    xLayer[3]=xstriphitsL3;
    yLayer[3]=ystriphitsL3;
    xLayer[4]=xstriphitsL4;
    yLayer[4]=ystriphitsL4;
    xLayer[5]=xstriphitsL5;
    yLayer[5]=ystriphitsL5;
    xLayer[6]=xstriphitsL6;
    yLayer[6]=ystriphitsL6;
    xLayer[7]=xstriphitsL7;
    yLayer[7]=ystriphitsL7;
    xLayer[8]=xstriphitsL8;
    yLayer[8]=ystriphitsL8;
    xLayer[9]=xstriphitsL9;
    yLayer[9]=ystriphitsL9;
    xLayer[10]=xstriphitsL10;
    yLayer[10]=ystriphitsL10;
    xLayer[11]=xstriphitsL11;
    yLayer[11]=ystriphitsL11;

    Long64_t nentries = fChain->GetEntriesFast();
    
    if (fChain == 0) return;
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
    }
}
