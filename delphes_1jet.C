#define delphes_1jet_cxx
#include "delphes_1jet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <vector>
#include <TTree.h>
#include <TFile.h>

void delphes_1jet::Loop()
{
    
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
    
    TFile* delphesFile = new TFile("delphes_n2n2_200_1jets.root"); //open the old root file from delphes
    TTree* delphesTree; //the old tree
    delphesFile->GetObject("Delphes",delphesTree); //access the old tree
    delphesTree->SetBranchStatus("*",0); //disable all branches
    delphesTree->SetBranchStatus("MissingET",1); //activate the needed branches
    delphesTree->SetBranchStatus("ScalarHT",1);
    delphesTree->SetBranchStatus("Electron",1);
    delphesTree->SetBranchStatus("Muon",1);
    
   TFile treefile("n2n2j_output.root","RECREATE"); //the output file
    
    TTree* n2n2_1jet = delphesTree->CloneTree(0); //copy the delphes tree to a new tree
    n2n2_1jet->CopyEntries(delphesTree); //get the entries
    
    TLorentzVector zeroTLV; //make a TLorentzVector (TLV) and set all entries to zero
    zeroTLV.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
    TLorentzVector n2A_MC; //make the MonteCarlo branches
    TLorentzVector n2B_MC;
    TLorentzVector n1A_MC;
    TLorentzVector n1B_MC;
    TLorentzVector lA_MC;
    TLorentzVector lB_MC;
    TLorentzVector lC_MC;
    TLorentzVector lD_MC;
    
    n2A_MC=zeroTLV; //set all MC branches to 0
    n2B_MC=zeroTLV;
    n1A_MC=zeroTLV;
    n1B_MC=zeroTLV;
    lA_MC=zeroTLV;
    lB_MC=zeroTLV;
    lC_MC=zeroTLV;
    lD_MC=zeroTLV;
    
    n2n2_1jet->Branch("n2A_MC","TLorentzVector",&n2A_MC); //add branches to the new tree
    n2n2_1jet->Branch("n2B_MC","TLorentzVector",&n2B_MC);
    n2n2_1jet->Branch("n1A_MC","TLorentzVector",&n1A_MC);
    n2n2_1jet->Branch("n1B_MC","TLorentzVector",&n1B_MC);
    n2n2_1jet->Branch("lA_MC","TLorentzVector",&lA_MC);
    n2n2_1jet->Branch("lB_MC","TLorentzVector",&lB_MC);
    n2n2_1jet->Branch("lC_MC","TLorentzVector",&lC_MC);
    n2n2_1jet->Branch("lD_MC","TLorentzVector",&lD_MC);
    
    
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //Event Loop
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       
       for(int i=0; i<Particle_size; i++) //Particle Loop for pythia/MC
       {
           if(Particle_PID[i]==1000023 && Particle_Status[i]==2) //n2
           {
               if(n2A_MC==zeroTLV) //fill the first n2
               {
                   n2A_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
               else //fill the second n2
               {
                   n2B_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
           }
           if(Particle_PID[i]==1000022 && Particle_Status[i]==1) //n1
           {
               if(n1A_MC==zeroTLV)
               {
                   n1A_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
               else
               {
                   n1B_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
           }
           if(Particle_PID[Particle_M1[i]]==23 && (abs(Particle_PID[i])==11 || abs(Particle_PID[i])==13)) //the leptons coming from the Zs
           {
               if(lA_MC==zeroTLV)
               {
                   lA_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
               else if(lB_MC==zeroTLV)
               {
                   lB_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
               else if(lC_MC==zeroTLV)
               {
                   lC_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
               else
               {
                   lD_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
           }
       }
       // in case a TLV is the same as another one
       if(n2A_MC==n2B_MC)
       {
           cout << "FAIL N2" << endl;
       }
       if(n1A_MC==n1B_MC)
       {
           cout << "FAIL N1" << endl;
       }
       if(lA_MC==lB_MC || lA_MC==lC_MC || lA_MC==lD_MC || lB_MC==lC_MC || lB_MC==lD_MC || lC_MC==lD_MC)
       {
           cout << "FAIL l" << endl;
       }
       n2n2_1jet->Fill(); //fill the tree
       n2A_MC=zeroTLV; //Reset the TLVs for the next event
       n2B_MC=zeroTLV;
       n1A_MC=zeroTLV;
       n1B_MC=zeroTLV;
       lA_MC=zeroTLV;
       lB_MC=zeroTLV;
       lC_MC=zeroTLV;
       lD_MC=zeroTLV;
   }
    //Write and Close the new tree/file
    n2n2_1jet->Write();
    treefile.Close();
}
