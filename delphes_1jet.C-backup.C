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
#include "Bonus.h"

void delphes_1jet::Loop()
{
    
   if (fChain == 0) return;
    
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
    
    TFile treefile("n2n2j_output.root","RECREATE"); //the output file
    
    TTree* n2n2_1jet = new TTree("n2n2j","Output");
    
    vector<double> n1_MC_PT;
    vector<double> n1_MC_Eta;
    vector<double> n1_MC_Phi;
    vector<double> n1_MC_E;
    Int_t n1_MC_Size=0;
    
    vector<double> n2_MC_PT;
    vector<double> n2_MC_Eta;
    vector<double> n2_MC_Phi;
    vector<double> n2_MC_E;
    Int_t n2_MC_Size=0;
    
    vector<double> GenJet_Delphes_PT;
    vector<double> GenJet_Delphes_Eta;
    vector<double> GenJet_Delphes_Phi;
    vector<double> GenJet_Delphes_E;
    vector<double> GenJet_Delphes_Mass;
    Int_t GenJet_Delphes_Size=0;
    
    vector<double> Jet_Delphes_PT;
    vector<double> Jet_Delphes_Eta;
    vector<double> Jet_Delphes_Phi;
    vector<double> Jet_Delphes_E;
    vector<double> Jet_Delphes_Mass;
    Int_t Jet_Delphes_Size=0;
    
    vector<double> JetPUPPI_Delphes_PT;
    vector<double> JetPUPPI_Delphes_Eta;
    vector<double> JetPUPPI_Delphes_Phi;
    vector<double> JetPUPPI_Delphes_E;
    vector<double> JetPUPPI_Delphes_Mass;
    Int_t JetPUPPI_Delphes_Size=0;
    
    vector<double> ElectronSignal_MC_PT;
    vector<double> ElectronSignal_MC_Eta;
    vector<double> ElectronSignal_MC_Phi;
    vector<double> ElectronSignal_MC_E;
    vector<double> ElectronSignal_MC_Charge;
    Int_t ElectronSignal_MC_Size=0;
    
    vector<double> MuonSignal_MC_PT;
    vector<double> MuonSignal_MC_Eta;
    vector<double> MuonSignal_MC_Phi;
    vector<double> MuonSignal_MC_E;
    vector<double> MuonSignal_MC_Charge;
    Int_t MuonSignal_MC_Size=0;
    
    vector<double> Electron_MC_PT;
    vector<double> Electron_MC_Eta;
    vector<double> Electron_MC_Phi;
    vector<double> Electron_MC_E;
    vector<double> Electron_MC_Charge;
    Int_t Electron_MC_Size=0;
    
    vector<double> Muon_MC_PT;
    vector<double> Muon_MC_Eta;
    vector<double> Muon_MC_Phi;
    vector<double> Muon_MC_E;
    vector<double> Muon_MC_Charge;
    Int_t Muon_MC_Size=0;
    
    Float_t MET_PT;
    Float_t MET_Eta;
    Float_t MET_Phi;
    
    Float_t PUPPI_MET_PT;
    Float_t PUPPI_MET_Eta;
    Float_t PUPPI_MET_Phi;
    
    Float_t HT=0.0;
    
    vector<double> Electron_Delphes_PT;
    vector<double> Electron_Delphes_Eta;
    vector<double> Electron_Delphes_Phi;
    vector<double> Electron_Delphes_E;
    vector<double> Electron_Delphes_Charge;
    Int_t Electron_Delphes_Size=0;
    
    vector<double> ElectronCHS_Delphes_PT;
    vector<double> ElectronCHS_Delphes_Eta;
    vector<double> ElectronCHS_Delphes_Phi;
    vector<double> ElectronCHS_Delphes_E;
    vector<double> ElectronCHS_Delphes_Charge;
    Int_t ElectronCHS_Delphes_Size=0;
    
    vector<double> MuonLoose_Delphes_PT;
    vector<double> MuonLoose_Delphes_Eta;
    vector<double> MuonLoose_Delphes_Phi;
    vector<double> MuonLoose_Delphes_E;
    vector<double> MuonLoose_Delphes_Charge;
    Int_t MuonLoose_Delphes_Size=0;
 
    vector<double> MuonLooseCHS_Delphes_PT;
    vector<double> MuonLooseCHS_Delphes_Eta;
    vector<double> MuonLooseCHS_Delphes_Phi;
    vector<double> MuonLooseCHS_Delphes_E;
    vector<double> MuonLooseCHS_Delphes_Charge;
    Int_t MuonLooseCHS_Delphes_Size=0;
    
    vector<double> MuonTight_Delphes_PT;
    vector<double> MuonTight_Delphes_Eta;
    vector<double> MuonTight_Delphes_Phi;
    vector<double> MuonTight_Delphes_E;
    vector<double> MuonTight_Delphes_Charge;
    Int_t MuonTight_Delphes_Size=0;
    
    vector<double> MuonTightCHS_Delphes_PT;
    vector<double> MuonTightCHS_Delphes_Eta;
    vector<double> MuonTightCHS_Delphes_Phi;
    vector<double> MuonTightCHS_Delphes_E;
    vector<double> MuonTightCHS_Delphes_Charge;
    Int_t MuonTightCHS_Delphes_Size=0;
    
    
    //Make Branches
    
    n2n2_1jet->Branch("n1_MC_PT","vector<double>",&n1_MC_PT);
    n2n2_1jet->Branch("n1_MC_Eta","vector<double>",&n1_MC_Eta);
    n2n2_1jet->Branch("n1_MC_Phi","vector<double>",&n1_MC_Phi);
    n2n2_1jet->Branch("n1_MC_E","vector<double>",&n1_MC_E);
    n2n2_1jet->Branch("n1_MC_Size",&n1_MC_Size,"n1_MC_Size/I");
    
    n2n2_1jet->Branch("n2_MC_PT","vector<double>",&n2_MC_PT);
    n2n2_1jet->Branch("n2_MC_Eta","vector<double>",&n2_MC_Eta);
    n2n2_1jet->Branch("n2_MC_Phi","vector<double>",&n2_MC_Phi);
    n2n2_1jet->Branch("n2_MC_E","vector<double>",&n2_MC_E);
    n2n2_1jet->Branch("n2_MC_Size",&n2_MC_Size,"n2_MC_Size/I");
    
    n2n2_1jet->Branch("GenJet_Delphes_PT","vector<double>",&GenJet_Delphes_PT);
    n2n2_1jet->Branch("GenJet_Delphes_Eta","vector<double>",&GenJet_Delphes_Eta);
    n2n2_1jet->Branch("GenJet_Delphes_Phi","vector<double>",&GenJet_Delphes_Phi);
    n2n2_1jet->Branch("GenJet_Delphes_E","vector<double>",&GenJet_Delphes_E);
    n2n2_1jet->Branch("GenJet_Delphes_Size",&GenJet_Delphes_Size,"GenJet_Delphes_Size/I");
    
    n2n2_1jet->Branch("Jet_Delphes_PT","vector<double>",&Jet_Delphes_PT);
    n2n2_1jet->Branch("Jet_Delphes_Eta","vector<double>",&Jet_Delphes_Eta);
    n2n2_1jet->Branch("Jet_Delphes_Phi","vector<double>",&Jet_Delphes_Phi);
    n2n2_1jet->Branch("Jet_Delphes_E","vector<double>",&Jet_Delphes_E);
    n2n2_1jet->Branch("Jet_Delphes_Size",&Jet_Delphes_Size,"Jet_Delphes_Size/I");
    
    n2n2_1jet->Branch("JetPUPPI_Delphes_PT","vector<double>",&JetPUPPI_Delphes_PT);
    n2n2_1jet->Branch("JetPUPPI_Delphes_Eta","vector<double>",&JetPUPPI_Delphes_Eta);
    n2n2_1jet->Branch("JetPUPPI_Delphes_Phi","vector<double>",&JetPUPPI_Delphes_Phi);
    n2n2_1jet->Branch("JetPUPPI_Delphes_E","vector<double>",&JetPUPPI_Delphes_E);
    n2n2_1jet->Branch("JetPUPPI_Delphes_Size",&JetPUPPI_Delphes_Size,"JetPUPPI_Delphes_Size/I");
    
    n2n2_1jet->Branch("ElectronSignal_MC_PT","vector<double>",&ElectronSignal_MC_PT);
    n2n2_1jet->Branch("ElectronSignal_MC_Eta","vector<double>",&ElectronSignal_MC_Eta);
    n2n2_1jet->Branch("ElectronSignal_MC_Phi","vector<double>",&ElectronSignal_MC_Phi);
    n2n2_1jet->Branch("ElectronSignal_MC_E","vector<double>",&ElectronSignal_MC_E);
    n2n2_1jet->Branch("ElectronSignal_MC_Size",&ElectronSignal_MC_Size,"ElectronSignal_MC_Size/I");
    
    n2n2_1jet->Branch("MuonSignal_MC_PT","vector<double>",&MuonSignal_MC_PT);
    n2n2_1jet->Branch("MuonSignal_MC_Eta","vector<double>",&MuonSignal_MC_Eta);
    n2n2_1jet->Branch("MuonSignal_MC_Phi","vector<double>",&MuonSignal_MC_Phi);
    n2n2_1jet->Branch("MuonSignal_MC_E","vector<double>",&MuonSignal_MC_E);
    n2n2_1jet->Branch("MuonSignal_MC_Size",&MuonSignal_MC_Size,"MuonSignal_MC_Size/I");
    
    n2n2_1jet->Branch("Electron_MC_PT","vector<double>",&Electron_MC_PT);
    n2n2_1jet->Branch("Electron_MC_Eta","vector<double>",&Electron_MC_Eta);
    n2n2_1jet->Branch("Electron_MC_Phi","vector<double>",&Electron_MC_Phi);
    n2n2_1jet->Branch("Electron_MC_E","vector<double>",&Electron_MC_E);
    n2n2_1jet->Branch("Electron_MC_Size",&Electron_MC_Size,"Electron_MC_Size/I");
    
    n2n2_1jet->Branch("Muon_MC_PT","vector<double>",&Muon_MC_PT);
    n2n2_1jet->Branch("Muon_MC_Eta","vector<double>",&Muon_MC_Eta);
    n2n2_1jet->Branch("Muon_MC_Phi","vector<double>",&Muon_MC_Phi);
    n2n2_1jet->Branch("Muon_MC_E","vector<double>",&Muon_MC_E);
    n2n2_1jet->Branch("Muon_MC_Size",&Muon_MC_Size,"Muon_MC_Size/I");
    
    n2n2_1jet->Branch("MET_PT",&MET_PT);
    n2n2_1jet->Branch("MET_Eta",&MET_Eta);
    n2n2_1jet->Branch("MET_Phi",&MET_Phi);
    
    n2n2_1jet->Branch("PUPPI_MET_PT",&PUPPI_MET_PT);
    n2n2_1jet->Branch("PUPPI_MET_Eta",&PUPPI_MET_Eta);
    n2n2_1jet->Branch("PUPPI_MET_Phi",&PUPPI_MET_Phi);
    
    n2n2_1jet->Branch("HT",&HT);
    
    n2n2_1jet->Branch("Electron_Delphes_PT","vector<double>",&Electron_Delphes_PT);
    n2n2_1jet->Branch("Electron_Delphes_Eta","vector<double>",&Electron_Delphes_Eta);
    n2n2_1jet->Branch("Electron_Delphes_Phi","vector<double>",&Electron_Delphes_Phi);
    n2n2_1jet->Branch("Electron_Delphes_E","vector<double>",&Electron_Delphes_E);
    n2n2_1jet->Branch("Electron_Delphes_Size",&Electron_Delphes_Size,"Electron_Delphes_Size/I");
    
    n2n2_1jet->Branch("ElectronCHS_Delphes_PT","vector<double>",&ElectronCHS_Delphes_PT);
    n2n2_1jet->Branch("ElectronCHS_Delphes_Eta","vector<double>",&ElectronCHS_Delphes_Eta);
    n2n2_1jet->Branch("ElectronCHS_Delphes_Phi","vector<double>",&ElectronCHS_Delphes_Phi);
    n2n2_1jet->Branch("ElectronCHS_Delphes_E","vector<double>",&ElectronCHS_Delphes_E);
    n2n2_1jet->Branch("ElectronCHS_Delphes_Size",&ElectronCHS_Delphes_Size,"ElectronCHS_Delphes_Size/I");
    
    n2n2_1jet->Branch("MuonLoose_Delphes_PT","vector<double>",&MuonLoose_Delphes_PT);
    n2n2_1jet->Branch("MuonLoose_Delphes_Eta","vector<double>",&MuonLoose_Delphes_Eta);
    n2n2_1jet->Branch("MuonLoose_Delphes_Phi","vector<double>",&MuonLoose_Delphes_Phi);
    n2n2_1jet->Branch("MuonLoose_Delphes_E","vector<double>",&MuonLoose_Delphes_E);
    n2n2_1jet->Branch("MuonLoose_Delphes_Size",&MuonLoose_Delphes_Size,"MuonLoose_Delphes_Size/I");
    
    n2n2_1jet->Branch("MuonLooseCHS_Delphes_PT","vector<double>",&MuonLooseCHS_Delphes_PT);
    n2n2_1jet->Branch("MuonLooseCHS_Delphes_Eta","vector<double>",&MuonLooseCHS_Delphes_Eta);
    n2n2_1jet->Branch("MuonLooseCHS_Delphes_Phi","vector<double>",&MuonLooseCHS_Delphes_Phi);
    n2n2_1jet->Branch("MuonLooseCHS_Delphes_E","vector<double>",&MuonLooseCHS_Delphes_E);
    n2n2_1jet->Branch("MuonLooseCHS_Delphes_Size",&MuonLooseCHS_Delphes_Size,"MuonLooseCHS_Delphes_Size/I");
    
    n2n2_1jet->Branch("MuonTight_Delphes_PT","vector<double>",&MuonTight_Delphes_PT);
    n2n2_1jet->Branch("MuonTight_Delphes_Eta","vector<double>",&MuonTight_Delphes_Eta);
    n2n2_1jet->Branch("MuonTight_Delphes_Phi","vector<double>",&MuonTight_Delphes_Phi);
    n2n2_1jet->Branch("MuonTight_Delphes_E","vector<double>",&MuonTight_Delphes_E);
    n2n2_1jet->Branch("MuonTight_Delphes_Size",&MuonTight_Delphes_Size,"MuonTight_Delphes_Size/I");
    
    n2n2_1jet->Branch("MuonTightCHS_Delphes_PT","vector<double>",&MuonTightCHS_Delphes_PT);
    n2n2_1jet->Branch("MuonTightCHS_Delphes_Eta","vector<double>",&MuonTightCHS_Delphes_Eta);
    n2n2_1jet->Branch("MuonTightCHS_Delphes_Phi","vector<double>",&MuonTightCHS_Delphes_Phi);
    n2n2_1jet->Branch("MuonTightCHS_Delphes_E","vector<double>",&MuonTightCHS_Delphes_E);
    n2n2_1jet->Branch("MuonTightCHS_Delphes_Size",&MuonTightCHS_Delphes_Size,"MuonTightCHS_Delphes_Size/I");
    
    //int dummycount=0;
    
    int jettag = 0;
    
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   { //Event Loop
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       
       n1_MC_PT.clear();
       n1_MC_Eta.clear();
       n1_MC_Phi.clear();
       n1_MC_E.clear();
       n1_MC_Size=0;
       
       n2_MC_PT.clear();
       n2_MC_Eta.clear();
       n2_MC_Phi.clear();
       n2_MC_E.clear();
       n2_MC_Size=0;
       
       GenJet_Delphes_PT.clear();
       GenJet_Delphes_Eta.clear();
       GenJet_Delphes_Phi.clear();
       GenJet_Delphes_E.clear();
       GenJet_Delphes_Mass.clear();
       GenJet_Delphes_Size=0;
       
       Jet_Delphes_PT.clear();
       Jet_Delphes_Eta.clear();
       Jet_Delphes_Phi.clear();
       Jet_Delphes_E.clear();
       Jet_Delphes_Mass.clear();
       Jet_Delphes_Size=0;
       
       JetPUPPI_Delphes_PT.clear();
       JetPUPPI_Delphes_Eta.clear();
       JetPUPPI_Delphes_Phi.clear();
       JetPUPPI_Delphes_E.clear();
       JetPUPPI_Delphes_Mass.clear();
       JetPUPPI_Delphes_Size=0;
       
       ElectronSignal_MC_PT.clear();
       ElectronSignal_MC_Eta.clear();
       ElectronSignal_MC_Phi.clear();
       ElectronSignal_MC_E.clear();
       ElectronSignal_MC_Charge.clear();
       ElectronSignal_MC_Size=0;
       
       MuonSignal_MC_PT.clear();
       MuonSignal_MC_Eta.clear();
       MuonSignal_MC_Phi.clear();
       MuonSignal_MC_E.clear();
       MuonSignal_MC_Charge.clear();
       MuonSignal_MC_Size=0;
       
       Electron_MC_PT.clear();
       Electron_MC_Eta.clear();
       Electron_MC_Phi.clear();
       Electron_MC_E.clear();
       Electron_MC_Charge.clear();
       Electron_MC_Size=0;
       
       Muon_MC_PT.clear();
       Muon_MC_Eta.clear();
       Muon_MC_Phi.clear();
       Muon_MC_E.clear();
       Muon_MC_Charge.clear();
       Muon_MC_Size=0;
        
       Electron_Delphes_PT.clear();
       Electron_Delphes_Eta.clear();
       Electron_Delphes_Phi.clear();
       Electron_Delphes_E.clear();
       Electron_Delphes_Charge.clear();
       Electron_Delphes_Size=Electron_size;
       
       ElectronCHS_Delphes_PT.clear();
       ElectronCHS_Delphes_Eta.clear();
       ElectronCHS_Delphes_Phi.clear();
       ElectronCHS_Delphes_E.clear();
       ElectronCHS_Delphes_Charge.clear();
       ElectronCHS_Delphes_Size=ElectronCHS_size;
       
       MuonLoose_Delphes_PT.clear();
       MuonLoose_Delphes_Eta.clear();
       MuonLoose_Delphes_Phi.clear();
       MuonLoose_Delphes_E.clear();
       MuonLoose_Delphes_Charge.clear();
       MuonLoose_Delphes_Size=MuonLoose_size;
       
       MuonLooseCHS_Delphes_PT.clear();
       MuonLooseCHS_Delphes_Eta.clear();
       MuonLooseCHS_Delphes_Phi.clear();
       MuonLooseCHS_Delphes_E.clear();
       MuonLooseCHS_Delphes_Charge.clear();
       MuonLooseCHS_Delphes_Size=MuonLooseCHS_size;
       
       MuonTight_Delphes_PT.clear();
       MuonTight_Delphes_Eta.clear();
       MuonTight_Delphes_Phi.clear();
       MuonTight_Delphes_E.clear();
       MuonTight_Delphes_Charge.clear();
       MuonTight_Delphes_Size=MuonTight_size;
       
       MuonTightCHS_Delphes_PT.clear();
       MuonTightCHS_Delphes_Eta.clear();
       MuonTightCHS_Delphes_Phi.clear();
       MuonTightCHS_Delphes_E.clear();
       MuonTightCHS_Delphes_Charge.clear();
       MuonTightCHS_Delphes_Size=MuonTightCHS_size;
       
       
       MET_PT=MissingET_MET[0];
       MET_Eta=MissingET_Eta[0];
       MET_Phi=MissingET_Phi[0];
       
       PUPPI_MET_PT=PuppiMissingET_MET[0];
       PUPPI_MET_Eta=PuppiMissingET_Eta[0];
       PUPPI_MET_Phi=PuppiMissingET_Phi[0];
       
       HT=ScalarHT_HT[0];
       for(int i = 0; i<Particle_size; i++)
       {
           if((abs(Particle_PID[i])==21 || abs(Particle_PID[i])==1 || abs(Particle_PID[i])==2 || abs(Particle_PID[i])==3 || abs(Particle_PID[i])==4) && Particle_Status[i]==2)
           {
               jettag++;
           }
       }
       if(jettag == 0) continue;
   
       for(int i=0; i<Particle_size; i++) //Particle Loop for pythia/MC
       {
           if(Particle_PID[i]==1000023 && Particle_Status[i]==2 && (Particle_PID[Particle_D1[i]] == 23 || Particle_PID[Particle_D1[i]] == 1000022)) //n2
           {
               n2_MC_PT.push_back(Particle_PT[i]);
               n2_MC_Eta.push_back(Particle_Eta[i]);
               n2_MC_Phi.push_back(Particle_Phi[i]);
               n2_MC_E.push_back(Particle_E[i]);
               n2_MC_Size++;
           }
           
           if(Particle_PID[i]==1000022 && Particle_Status[i]==1) //n1
           {
               n1_MC_PT.push_back(Particle_PT[i]);
               n1_MC_Eta.push_back(Particle_Eta[i]);
               n1_MC_Phi.push_back(Particle_Phi[i]);
               n1_MC_E.push_back(Particle_E[i]);
               n1_MC_Size++;
           }
           if((abs(Particle_PID[i])==11) && Particle_Status[i]==1)
           {
               ElectronSignal_MC_PT.push_back(Particle_PT[i]);
               ElectronSignal_MC_Eta.push_back(Particle_Eta[i]);
               ElectronSignal_MC_Phi.push_back(Particle_Phi[i]);
               ElectronSignal_MC_E.push_back(Particle_E[i]);
               ElectronSignal_MC_Size++;
           }
           if((abs(Particle_PID[i])==13) && Particle_Status[i]==1)
           {
               MuonSignal_MC_PT.push_back(Particle_PT[i]);
               MuonSignal_MC_Eta.push_back(Particle_Eta[i]);
               MuonSignal_MC_Phi.push_back(Particle_Phi[i]);
               MuonSignal_MC_E.push_back(Particle_E[i]);
               MuonSignal_MC_Size++;
           }
           if(abs(Particle_PID[i]) == 11)
           {
               Electron_MC_PT.push_back(Particle_PT[i]);
               Electron_MC_Eta.push_back(Particle_Eta[i]);
               Electron_MC_Phi.push_back(Particle_Phi[i]);
               Electron_MC_E.push_back(Particle_E[i]);
               Electron_MC_Size++;
           }
           if(abs(Particle_PID[i]) == 13)
           {
               Muon_MC_PT.push_back(Particle_PT[i]);
               Muon_MC_Eta.push_back(Particle_Eta[i]);
               Muon_MC_Phi.push_back(Particle_Phi[i]);
               Muon_MC_E.push_back(Particle_E[i]);
               Muon_MC_Size++;
           }
       }
       
       for(int i=0; i<GenJet_size; i++)
       {
           GenJet_Delphes_PT.push_back(GenJet_PT[i]);
           GenJet_Delphes_Eta.push_back(GenJet_Eta[i]);
           GenJet_Delphes_Phi.push_back(GenJet_Phi[i]);
           GenJet_Delphes_E.push_back(GenJet_T[i]);
           GenJet_Delphes_Mass.push_back(GenJet_Mass[i]);
       }
       GenJet_Delphes_Size=GenJet_size;
       
       for(int i=0; i<Jet_size; i++)
       {
           Jet_Delphes_PT.push_back(Jet_PT[i]);
           Jet_Delphes_Eta.push_back(Jet_Eta[i]);
           Jet_Delphes_Phi.push_back(Jet_Phi[i]);
           Jet_Delphes_E.push_back(Jet_T[i]);
           Jet_Delphes_Mass.push_back(Jet_Mass[i]);
       }
       Jet_Delphes_Size=Jet_size;
       
       for(int i=0; i<JetPUPPI_size; i++)
       {
           JetPUPPI_Delphes_PT.push_back(JetPUPPI_PT[i]);
           JetPUPPI_Delphes_Eta.push_back(JetPUPPI_Eta[i]);
           JetPUPPI_Delphes_Phi.push_back(JetPUPPI_Phi[i]);
           JetPUPPI_Delphes_E.push_back(JetPUPPI_T[i]);
           JetPUPPI_Delphes_Mass.push_back(JetPUPPI_Mass[i]);
       }
       JetPUPPI_Delphes_Size=JetPUPPI_size;
       
       for(int i=0; i<Electron_size; i++)
       {
           Electron_Delphes_PT.push_back(Electron_PT[i]);
           Electron_Delphes_Eta.push_back(Electron_Eta[i]);
           Electron_Delphes_Phi.push_back(Electron_Phi[i]);
           Electron_Delphes_E.push_back(Electron_T[i]);
           Electron_Delphes_Charge.push_back(Electron_Charge[i]);
       }
       Electron_Delphes_Size=Electron_size;
       
       for(int i=0; i<ElectronCHS_size; i++)
       {
           ElectronCHS_Delphes_PT.push_back(ElectronCHS_PT[i]);
           ElectronCHS_Delphes_Eta.push_back(ElectronCHS_Eta[i]);
           ElectronCHS_Delphes_Phi.push_back(ElectronCHS_Phi[i]);
           ElectronCHS_Delphes_E.push_back(ElectronCHS_T[i]);
           ElectronCHS_Delphes_Charge.push_back(ElectronCHS_Charge[i]);
       }
       ElectronCHS_Delphes_Size=ElectronCHS_size;
       
       for(int i=0; i<MuonLoose_size; i++)
       {
           MuonLoose_Delphes_PT.push_back(MuonLoose_PT[i]);
           MuonLoose_Delphes_Eta.push_back(MuonLoose_Eta[i]);
           MuonLoose_Delphes_Phi.push_back(MuonLoose_Phi[i]);
           MuonLoose_Delphes_E.push_back(MuonLoose_T[i]);
           MuonLoose_Delphes_Charge.push_back(MuonLoose_Charge[i]);
       }
       MuonLoose_Delphes_Size=MuonLoose_size;
       
       for(int i=0; i<MuonLooseCHS_size; i++)
       {
           MuonLooseCHS_Delphes_PT.push_back(MuonLooseCHS_PT[i]);
           MuonLooseCHS_Delphes_Eta.push_back(MuonLooseCHS_Eta[i]);
           MuonLooseCHS_Delphes_Phi.push_back(MuonLooseCHS_Phi[i]);
           MuonLooseCHS_Delphes_E.push_back(MuonLooseCHS_T[i]);
           MuonLooseCHS_Delphes_Charge.push_back(MuonLooseCHS_Charge[i]);
       }
       MuonLooseCHS_Delphes_Size=MuonLooseCHS_size;
       
       for(int i=0; i<MuonTight_size; i++)
       {
           MuonTight_Delphes_PT.push_back(MuonTight_PT[i]);
           MuonTight_Delphes_Eta.push_back(MuonTight_Eta[i]);
           MuonTight_Delphes_Phi.push_back(MuonTight_Phi[i]);
           MuonTight_Delphes_E.push_back(MuonTight_T[i]);
           MuonTight_Delphes_Charge.push_back(MuonTight_Charge[i]);
       }
       MuonTight_Delphes_Size=MuonTight_size;
       
       for(int i=0; i<MuonTightCHS_size; i++)
       {
           MuonTightCHS_Delphes_PT.push_back(MuonTightCHS_PT[i]);
           MuonTightCHS_Delphes_Eta.push_back(MuonTightCHS_Eta[i]);
           MuonTightCHS_Delphes_Phi.push_back(MuonTightCHS_Phi[i]);
           MuonTightCHS_Delphes_E.push_back(MuonTightCHS_T[i]);
           MuonTightCHS_Delphes_Charge.push_back(MuonTightCHS_Charge[i]);
       }
       MuonTightCHS_Delphes_Size=MuonTightCHS_size;
       
       n2n2_1jet->Fill(); //fill the tree
   }
    //cout << dummycount << endl;
    //Write and Close the new tree/file
    n2n2_1jet->Write("", TObject::kOverwrite);
    treefile.Close();
}
