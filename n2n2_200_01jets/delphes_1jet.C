/* ----------------------------------------------------------------------------
 Author: Zach Flowers
 Description:  Makes a tree from the delphes output file
    Run before Analyze_n2n2j.C
 Date: April 23rd 2018
 * ---------------------------------------------------------------------------- */

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
    
   TFile treefile("n2n2j_output.root","RECREATE"); //create the output file
    
    TTree* n2n2_1jet = new TTree("n2n2j","Output"); //create the output tree
    
    TLorentzVector zeroTLV; //make a TLorentzVector (TLV) and set all entries to zero
    zeroTLV.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
    //make the MonteCarlo TLVs
    //second lightest neutralinos
    TLorentzVector n2A_MC;
    TLorentzVector n2B_MC;
    //lightest neutralinos
    TLorentzVector n1A_MC;
    TLorentzVector n1B_MC;
    //leptons
    TLorentzVector lA_MC;
    TLorentzVector lB_MC;
    TLorentzVector lC_MC;
    TLorentzVector lD_MC;
    
    //Store information about the jet(s) in an event
    vector<double> Jet_MC_PT;
    vector<double> Jet_MC_Eta;
    vector<double> Jet_MC_Phi;
    vector<double> Jet_MC_E;
    vector<double> Jet_MC_Mass;
    //number of entries or jets in the event
    Int_t Jet_MC_Size;
    
    //Delphes Branches
    
    //Missing Transverse Energy
    TVector3 MET;
    //Scalar HT
    Float_t HT=0.0;
    Int_t HT_Size;
    
    //Electron in Delphes
    vector<double> Electron_Delphes_PT;
    vector<double> Electron_Delphes_Eta;
    vector<double> Electron_Delphes_Phi;
    vector<double> Electron_Delphes_E;
    vector<double> Electron_Delphes_Charge;
    //number of electrons in an event
    Int_t Electron_Delphes_Size;
    
    //Muon in Delphes
    vector<double> Muon_Delphes_PT;
    vector<double> Muon_Delphes_Eta;
    vector<double> Muon_Delphes_Phi;
    vector<double> Muon_Delphes_E;
    vector<double> Muon_Delphes_Charge;
    //number of muons in an event
    Int_t Muon_Delphes_Size;
    
    //Jets in Delphes
    vector<double> Jet_Delphes_PT;
    vector<double> Jet_Delphes_Eta;
    vector<double> Jet_Delphes_Phi;
    vector<double> Jet_Delphes_E;
    vector<double> Jet_Delphes_Mass;
    //number of jets in an event
    Int_t Jet_Delphes_Size;
    
    //set all MC branches to 0
    n2A_MC=zeroTLV;
    n2B_MC=zeroTLV;
    n1A_MC=zeroTLV;
    n1B_MC=zeroTLV;
    lA_MC=zeroTLV;
    lB_MC=zeroTLV;
    lC_MC=zeroTLV;
    lD_MC=zeroTLV;
    
    //Make Branches in the tree for every variable declared above
    n2n2_1jet->Branch("MET","TVector3",&MET);
    n2n2_1jet->Branch("HT",&HT);
    n2n2_1jet->Branch("HT_Size",&HT_Size);
    n2n2_1jet->Branch("Electron_Delphes_PT","vector<double>",&Electron_Delphes_PT);
    n2n2_1jet->Branch("Electron_Delphes_Eta","vector<double>",&Electron_Delphes_Eta);
    n2n2_1jet->Branch("Electron_Delphes_Phi","vector<double>",&Electron_Delphes_Phi);
    n2n2_1jet->Branch("Electron_Delphes_E","vector<double>",&Electron_Delphes_E);
    n2n2_1jet->Branch("Electron_Delphes_Charge","vector<double>",&Electron_Delphes_Charge);
    n2n2_1jet->Branch("Electron_Delphes_Size",&Electron_Delphes_Size,"Electron_Delphes_Size/I");
    n2n2_1jet->Branch("Muon_Delphes_PT","vector<double>",&Muon_Delphes_PT);
    n2n2_1jet->Branch("Muon_Delphes_Eta","vector<double>",&Muon_Delphes_Eta);
    n2n2_1jet->Branch("Muon_Delphes_Phi","vector<double>",&Muon_Delphes_Phi);
    n2n2_1jet->Branch("Muon_Delphes_E","vector<double>",&Muon_Delphes_E);
    n2n2_1jet->Branch("Muon_Delphes_Charge","vector<double>",&Muon_Delphes_Charge);
    n2n2_1jet->Branch("Muon_Delphes_Size",&Muon_Delphes_Size,"Muon_Delphes_Size/I");
    n2n2_1jet->Branch("Jet_Delphes_PT","vector<double>",&Jet_Delphes_PT);
    n2n2_1jet->Branch("Jet_Delphes_Eta","vector<double>",&Jet_Delphes_Eta);
    n2n2_1jet->Branch("Jet_Delphes_Phi","vector<double>",&Jet_Delphes_Phi);
    n2n2_1jet->Branch("Jet_Delphes_E","vector<double>",&Jet_Delphes_E);
    n2n2_1jet->Branch("Jet_Delphes_Mass","vector<double>",&Jet_Delphes_Mass);
    n2n2_1jet->Branch("Jet_Delphes_Size",&Jet_Delphes_Size,"Jet_Delphes_Size/I");
    
    
    n2n2_1jet->Branch("n2A_MC","TLorentzVector",&n2A_MC);
    n2n2_1jet->Branch("n2B_MC","TLorentzVector",&n2B_MC);
    n2n2_1jet->Branch("n1A_MC","TLorentzVector",&n1A_MC);
    n2n2_1jet->Branch("n1B_MC","TLorentzVector",&n1B_MC);
    n2n2_1jet->Branch("lA_MC","TLorentzVector",&lA_MC);
    n2n2_1jet->Branch("lB_MC","TLorentzVector",&lB_MC);
    n2n2_1jet->Branch("lC_MC","TLorentzVector",&lC_MC);
    n2n2_1jet->Branch("lD_MC","TLorentzVector",&lD_MC);
    n2n2_1jet->Branch("Jet_MC_PT","vector<double>",&Jet_MC_PT);
    n2n2_1jet->Branch("Jet_MC_Eta","vector<double>",&Jet_MC_Eta);
    n2n2_1jet->Branch("Jet_MC_Phi","vector<double>",&Jet_MC_Phi);
    n2n2_1jet->Branch("Jet_MC_E","vector<double>",&Jet_MC_E);
    n2n2_1jet->Branch("Jet_MC_Mass","vector<double>",&Jet_MC_Mass);
    n2n2_1jet->Branch("Jet_MC_Size",&Jet_MC_Size,"Jet_MC_Size/I");
    
    HT_Size=nentries; //since there is only one value for HT for each event
    
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //Event Loop
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       
       //set the number of electrons, muons, and jets in the event from delphes
       Electron_Delphes_Size=Electron_size;
       Muon_Delphes_Size=Muon_size;
       Jet_Delphes_Size=Jet_size;
       
       //Looping over Delphes/RECO
       
       for(int i=0; i<MissingET_size; i++) //get MET
       {
           MET.SetPtEtaPhi(MissingET_MET[i],MissingET_Eta[i],MissingET_Phi[i]);
       }
       for(int i=0; i<ScalarHT_size; i++) //get HT
       {
           HT=ScalarHT_HT[i];
       }
       //get electrons from delphes
       for(int i=0; i<Electron_size; i++)
       {
           Electron_Delphes_PT.push_back(Electron_PT[i]);
           Electron_Delphes_Eta.push_back(Electron_Eta[i]);
           Electron_Delphes_Phi.push_back(Electron_Phi[i]);
           Electron_Delphes_E.push_back(Electron_T[i]);
           Electron_Delphes_Charge.push_back(Electron_Charge[i]);
       }
       //get muons from delphes
       for(int i=0; i<Muon_size; i++)
       {
           Muon_Delphes_PT.push_back(Muon_PT[i]);
           Muon_Delphes_Eta.push_back(Muon_Eta[i]);
           Muon_Delphes_Phi.push_back(Muon_Phi[i]);
           Muon_Delphes_E.push_back(Muon_T[i]);
           Muon_Delphes_Charge.push_back(Muon_Charge[i]);
       }
       //get jets from delphes
       for(int i=0; i<Jet_size; i++)
       {
           Jet_Delphes_PT.push_back(Jet_PT[i]);
           Jet_Delphes_Eta.push_back(Jet_Eta[i]);
           Jet_Delphes_Phi.push_back(Jet_Phi[i]);
           Jet_Delphes_E.push_back(Jet_T[i]);
           Jet_Delphes_Mass.push_back(Jet_Mass[i]);
       }
       
       //Done with Delphes
       
       //Next loop through Pythia/MC
       for(int i=0; i<Particle_size; i++) //Particle Loop for pythia/MC
       {
           if(Particle_PID[i]==1000023 && Particle_Status[i]==2) //n2
           {
               if(n2A_MC==zeroTLV) //fill the first n2
               {
                   n2A_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
               }
               else if(n2B_MC==zeroTLV) //fill the second n2
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
               else if(n1B_MC==zeroTLV)
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
           
           //pythia jet info
           if(((abs(Particle_PID[i])==21 || abs(Particle_PID[i])==1 || abs(Particle_PID[i])==2 || abs(Particle_PID[i])==3 || abs(Particle_PID[i])==4 || abs(Particle_PID[i])==5)) && Particle_Status[i]==2)
           {
               if(Particle_Status[i]==2 && Particle_Eta[i] < 800) //cut out a few gluons with unreasonable values
               {
                   Jet_MC_PT.push_back(Particle_PT[i]);
                   Jet_MC_Eta.push_back(Particle_Eta[i]);
                   Jet_MC_Phi.push_back(Particle_Phi[i]);
                   Jet_MC_E.push_back(Particle_E[i]);
                   Jet_MC_Mass.push_back(Particle_Mass[i]);
               }
           }
       }
       // in case a MC TLV is the same as another one
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
       Jet_MC_Size=Jet_MC_E.size();
       n2n2_1jet->Fill(); //fill the tree
       
       //reset all variables
       
       Electron_Delphes_PT.clear();
       Electron_Delphes_Eta.clear();
       Electron_Delphes_Phi.clear();
       Electron_Delphes_E.clear();
       Electron_Delphes_Charge.clear();
       
       Muon_Delphes_PT.clear();
       Muon_Delphes_Eta.clear();
       Muon_Delphes_Phi.clear();
       Muon_Delphes_E.clear();
       Muon_Delphes_Charge.clear();
       
       Jet_Delphes_PT.clear();
       Jet_Delphes_Eta.clear();
       Jet_Delphes_Phi.clear();
       Jet_Delphes_E.clear();
       Jet_Delphes_Mass.clear();
       
       Jet_MC_PT.clear();
       Jet_MC_Eta.clear();
       Jet_MC_Phi.clear();
       Jet_MC_E.clear();
       Jet_MC_Mass.clear();
       
       //Reset the MC TLVs for the next event
       n2A_MC=zeroTLV;
       n2B_MC=zeroTLV;
       n1A_MC=zeroTLV;
       n1B_MC=zeroTLV;
       lA_MC=zeroTLV;
       lB_MC=zeroTLV;
       lC_MC=zeroTLV;
       lD_MC=zeroTLV;
   }
    //Write and Close the tree and file
    n2n2_1jet->Write();
    treefile.Close();
}
