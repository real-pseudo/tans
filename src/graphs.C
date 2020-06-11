#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "geometry.h"
#include "particle.h"
#include "hit.h"
#include "Utility.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <algorithm>
#include <vector>
#include "TGraph.h"

void graphs2() {
  
  //grafica
  //leggo ntuple
  TFile filin("vtxreco.root");
  float sim,rec,diff_z,mult_rec,mult_sim,sim_z;
  TNtuple *pippo = (TNtuple*)filin.Get("nt");
      pippo->SetBranchAddress("zsim",&sim);
      pippo->SetBranchAddress("zrec",&rec);
      pippo->SetBranchAddress("diff",&diff_z);
      pippo->SetBranchAddress("mult",&mult_rec);
  TNtuple *nt1 = (TNtuple*)filin.Get("z_sim");
      nt1->SetBranchAddress("zsimtot",&sim_z);
      nt1->SetBranchAddress("multtot",&mult_sim);

  double sigma_z=5.3; //caratteristico del vtx generato
  int nsim=nt1->GetEntries();
  int nrec=pippo->GetEntries();
  
}