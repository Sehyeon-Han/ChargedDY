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

void Muon_data_2016_preVFP(const char* inFile,
                           const char* outFile = "Muon_data_2016_preVFP.root")
{
    double mTBins[] = {200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1500, 2000, 3500};
    int nmTBins = sizeof(mTBins) / sizeof(double) - 1;

    TH1D *h_Muon_pt = new TH1D("h_Muon_pt", "Muon p_{T}; p_{T} [GeV]; Events", 30, 65, 1500);
    TH1D *h_Muon_eta = new TH1D("h_Muon_eta", "Muon #eta; #eta; Events", 50, -3, 3);
    TH1D *h_Muon_phi = new TH1D("h_Muon_phi", "Muon #phi; #phi; Events", 50, -3.3, 3.3);
    TH1D *h_MET_pt = new TH1D("h_MET_pt", "MET p_{T}; MET p_{T} [GeV]; Events", 30, 85, 3500);
    TH1D *h_MET_phi = new TH1D("h_MET_phi", "MET #phi; #phi; Events", 50, -3.3, 3.3);
    TH1D *h_mT = new TH1D("h_mT", "m_{T}; m_{T} [GeV]; Events", nmTBins, mTBins);

    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_iso(reader, "Muon_tkRelIso");
    TTreeReaderArray<UChar_t> Muon_highPtId(reader, "Muon_highPtId");
    TTreeReaderValue<Bool_t> HLT_Mu50(reader, "HLT_Mu50");
    TTreeReaderValue<Bool_t> HLT_TkMu50(reader, "HLT_TkMu50");
    TTreeReaderValue<Float_t> MET_pt(reader, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(reader, "MET_phi");

    while(reader.Next())
    {   
        if(!(*HLT_Mu50 || *HLT_TkMu50)) continue;
        if(!(*MET_pt > 85.0)) continue;

        int selCount = 0;
        int selIdx = -1;

        bool ExtraMuon = false;

        const int n = Muon_pt.GetSize();

        for(int i = 0; i < n; ++i)
        {
            if(Muon_pt[i] > 65.0 && std::fabs(Muon_eta[i]) < 2.4 && Muon_iso[i] < 0.10 && Muon_highPtId[i] == 2)
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

                if(Muon_pt[i] > 20.0)
                {
                    ExtraMuon = true;
                    break;
                }
            }
        }

        if(selCount == 1 && selIdx >= 0 && !ExtraMuon)
        {
            double dphi = TVector2::Phi_mpi_pi(Muon_phi[selIdx] - *MET_phi);
            double mT = std::sqrt(2.0 * Muon_pt[selIdx] * (*MET_pt) * (1.0 - std::cos(dphi)));

            if(mT > 200)
            {
                h_Muon_pt->Fill(Muon_pt[selIdx]);
                h_Muon_eta->Fill(Muon_eta[selIdx]);
                h_Muon_phi->Fill(Muon_phi[selIdx]);
                h_MET_pt->Fill(*MET_pt);
                h_MET_phi->Fill(*MET_phi);
                h_mT->Fill(mT);
            }

            if(mT > 3000)
            {
                std::cout << "mT" << mT << std::endl;
            }
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