#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <cmath>
#include <TVector2.h>

void Muon_MC_2000_2016_postVFP(const char* inFile,
                               const char* outFile = "Muon_MC_2000_2016_postVFP.root")
{
    TH1::SetDefaultSumw2();
    
    TH1D *h_Muon_pt = new TH1D("h_Muon_pt", "Muon p_{T}; p_{T} [GeV]; Events", 30, 65, 1500);
    TH1D *h_Muon_eta = new TH1D("h_Muon_eta", "Muon #eta; #eta; Events", 50, -3, 3);
    TH1D *h_Muon_phi = new TH1D("h_Muon_phi", "Muon #phi; #phi; Events", 50, -3.3, 3.3);
    TH1D *h_MET_pt = new TH1D("h_MET_pt", "MET p_{T}; MET p_{T} [GeV]; Events", 30, 85, 3500);
    TH1D *h_MET_phi = new TH1D("h_MET_phi", "MET #phi; #phi; Events", 50, -3.3, 3.3);
    TH1D *h_mT = new TH1D("h_mT", "m_{T}; m_{T} [GeV]; Events", 30, 200, 3500);

    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_iso(reader, "Muon_pfRelIso04_all");
    TTreeReaderArray<UChar_t> Muon_highPtId(reader, "Muon_highPtId");
    TTreeReaderArray<Bool_t> Muon_looseId(reader, "Muon_looseId");
    TTreeReaderValue<Bool_t> HLT_Mu50(reader, "HLT_Mu50");
    TTreeReaderValue<Float_t> MET_pt(reader, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(reader, "MET_phi");
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");
    TTreeReaderValue<Float_t> L1PreFiringWeight_Nom(reader, "L1PreFiringWeight_Nom");

    while(reader.Next())
    {   
        double w = *genWeight * (*L1PreFiringWeight_Nom);

        if(!(*HLT_Mu50)) continue;
        if(!(*MET_pt > 85.0)) continue;

        int selCount = 0;
        int selIdx = -1;

        bool ExtraMuon = false;

        const int n = Muon_pt.GetSize();

        for(int i = 0; i < n; ++i)
        {
            if(Muon_pt[i] > 65.0 && std::fabs(Muon_eta[i]) < 2.4 && Muon_iso[i] < 0.12 && Muon_highPtId[i] == 2)
            {
                ++selCount;
                selIdx = i;
            }
        }

        if(selCount == 1 && selIdx >= 0)
        {
            for(int i = 0; i < n; ++i)
            {
                if(i == selIdx) continue;

                if(Muon_pt[i] > 20.0 && Muon_looseId[i] && std::fabs(Muon_eta[i]) < 2.4)
                {
                    ExtraMuon = true;
                    break;
                }
            }
        }

        if(selCount == 1 && selIdx >= 0 && !ExtraMuon)
        {
            h_Muon_pt->Fill(Muon_pt[selIdx], w);
            h_Muon_eta->Fill(Muon_eta[selIdx], w);
            h_Muon_phi->Fill(Muon_phi[selIdx], w);
            h_MET_pt->Fill(*MET_pt, w);
            h_MET_phi->Fill(*MET_phi, w);

            double dphi = TVector2::Phi_mpi_pi(Muon_phi[selIdx] - *MET_phi);
            double mT = std::sqrt(2.0 * Muon_pt[selIdx] * (*MET_pt) * (1.0 - std::cos(dphi)));

            if(mT > 2000)
            h_mT->Fill(mT, w);
        }
    }

    TFile fout(outFile, "RECREATE");
    h_Muon_pt->Write();
    h_Muon_eta->Write();
    h_Muon_phi->Write();
    h_MET_pt->Write();
    h_MET_phi->Write();
    h_mT->Write();

    fout.Close();
}