#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <cmath>
#include <TVector2.h>
#include <algorithm>

void Muon_DYMC_2000_2016_preVFP(const char* inFile = "/pnfs/knu.ac.kr/data/cms/store/user/sungwon/DY_Run2_UL_NanoAOD/2016_preVFP/DYJetsToMuMu_M-2000toInf/*/*/*/*/*.root",
                                const char* outFile = "Muon_DYMC_2000_2016_preVFP.root")
{
    TH1::SetDefaultSumw2();

    double mTBins[] = {200, 250, 300, 350, 425, 500, 600, 750, 900, 1100, 1400, 2000, 5000};
    int nmTBins = sizeof(mTBins) / sizeof(double) - 1;
    
    TH1D *h_Muon_pt = new TH1D("h_Muon_pt", "Muon p_{T}; p_{T} [GeV]; Events", 30, 65, 1500);
    TH1D *h_Muon_eta = new TH1D("h_Muon_eta", "Muon #eta; #eta; Events", 50, -3, 3);
    TH1D *h_Muon_phi = new TH1D("h_Muon_phi", "Muon #phi; #phi; Events", 50, -3.3, 3.3);
    TH1D *h_MET_pt = new TH1D("h_MET_pt", "MET p_{T}; MET p_{T} [GeV]; Events", 30, 85, 3500);
    TH1D *h_MET_phi = new TH1D("h_MET_phi", "MET #phi; #phi; Events", 50, -3.3, 3.3);
    TH1D *h_mT = new TH1D("h_mT", "m_{T}; m_{T} [GeV]; Events", nmTBins, mTBins);

    TFile *f_SF_ID = TFile::Open("eff_mu_ID.root");
    TFile *f_SF_ISO = TFile::Open("eff_mu_ISO.root");
    TFile *f_SF_STRIG = TFile::Open("eff_mu_STRIG.root");
    TFile *f_PU = TFile::Open("PUreweight.root");

    TH2D* h2D_SF_ID = (TH2D*)f_SF_ID->Get("NUM_HighPtID_DEN_TrackerMuons_abseta_pt")->Clone();
    TH2D* h2D_SF_ISO = (TH2D*)f_SF_ISO->Get("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt")->Clone();
    TH2D* h2D_SF_STRIG = (TH2D*)f_SF_STRIG->Get("NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt")->Clone();
    TH1D* hPUw = (TH1D*)f_PU->Get("hPUw");

    auto GetSF = [](TH2D* h, double aeta, double pt)
    {
        const double xmax = h->GetYaxis()->GetXmax();
        const double pt_eval = std::min(pt, std::nextafter(xmax, 0.0));

        int binx = h->GetXaxis()->FindFixBin(aeta);
        int biny = h->GetYaxis()->FindFixBin(pt_eval);

        binx = std::max(1, std::min(binx, h->GetNbinsX()));
        biny = std::max(1, std::min(biny, h->GetNbinsY()));

        return h->GetBinContent(binx, biny);
    };

    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_iso(reader, "Muon_tkRelIso");
    TTreeReaderArray<UChar_t> Muon_highPtId(reader, "Muon_highPtId");
    TTreeReaderArray<Bool_t> Muon_looseId(reader, "Muon_looseId");
    TTreeReaderValue<Bool_t> HLT_TkMu50(reader, "HLT_TkMu50");
    TTreeReaderValue<Bool_t> HLT_Mu50(reader, "HLT_Mu50");
    TTreeReaderValue<Float_t> MET_pt(reader, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(reader, "MET_phi");
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");
    TTreeReaderValue<Float_t> L1PreFiringWeight_Nom(reader, "L1PreFiringWeight_Nom");
    TTreeReaderValue<Float_t> Pileup_nTrueInt(reader, "Pileup_nTrueInt");

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

                if(Muon_pt[i] > 20.0 && Muon_looseId[i] && std::fabs(Muon_eta[i]) < 2.4)
                {
                    ExtraMuon = true;
                    break;
                }
            }
        }

        if(selCount == 1 && selIdx >= 0 && !ExtraMuon)
        {
            double pt = Muon_pt[selIdx];
            double aeta = fabs(Muon_eta[selIdx]);

            double sf_id = GetSF(h2D_SF_ID, aeta, pt);
            double sf_iso = GetSF(h2D_SF_ISO, aeta, pt);
            double sf_trig = GetSF(h2D_SF_STRIG, aeta, pt);

            int bin_pu = hPUw->GetXaxis()->FindBin(*Pileup_nTrueInt);
            bin_pu = std::max(1, std::min(bin_pu, hPUw->GetNbinsX()));
            double wpu = hPUw->GetBinContent(bin_pu);

            double w = (double)(*genWeight) * sf_id * sf_iso * sf_trig * wpu * (*L1PreFiringWeight_Nom);
            
            double dphi = TVector2::Phi_mpi_pi(Muon_phi[selIdx] - *MET_phi);
            double mT = std::sqrt(2.0 * Muon_pt[selIdx] * (*MET_pt) * (1.0 - std::cos(dphi)));

            if(mT > 2000)
            {
                h_Muon_pt->Fill(Muon_pt[selIdx], w);
                h_Muon_eta->Fill(Muon_eta[selIdx], w);
                h_Muon_phi->Fill(Muon_phi[selIdx], w);
                h_MET_pt->Fill(*MET_pt, w);
                h_MET_phi->Fill(*MET_phi, w);
                h_mT->Fill(mT, w);
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
