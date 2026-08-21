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
#include <TStopwatch.h>
#include <TStyle.h>

void Muon_fake_rate_DYMC_1500_dependence(const char* inFile = "/pnfs/knu.ac.kr/data/cms/store/user/sungwon/DY_Run2_UL_NanoAOD/2016_preVFP/DYJetsToMuMu_M-1500to2000/*/*/*/*/*.root",
                                            const char* outFile = "fake_mu_DYMC_1500_merged_dependence.root")
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    double metBins[] = {0.0, 20.0, 40.0, 65.0};
    const int N_metBins = 3;

    TH1D *h_fake_den_MET = new TH1D("h_fake_den_MET", "Fake-rate denominator;MET [GeV];Events", N_metBins, metBins);
    TH1D *h_fake_num_MET = new TH1D("h_fake_num_MET", "Fake-rate numerator;MET [GeV];Events", N_metBins, metBins);

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
    TTreeReaderArray<UChar_t> Muon_highPtId(reader, "Muon_highPtId");
    TTreeReaderArray<Float_t> Muon_tkRelIso(reader, "Muon_tkRelIso");
    TTreeReaderArray<Float_t> Muon_dxy(reader, "Muon_dxy");
    TTreeReaderArray<Float_t> Muon_dxyErr(reader, "Muon_dxyErr");

    TTreeReaderValue<Bool_t> HLT_TkMu50(reader, "HLT_TkMu50");
    TTreeReaderValue<Bool_t> HLT_Mu50(reader, "HLT_Mu50");

    TTreeReaderArray<Float_t> Jet_pt(reader, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_phi(reader, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_eta(reader, "Jet_eta");
    TTreeReaderArray<Int_t> Jet_jetId(reader, "Jet_jetId");

    TTreeReaderValue<Float_t> MET_pt(reader, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(reader, "MET_phi");
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");
    TTreeReaderValue<Float_t> Pileup_nTrueInt(reader, "Pileup_nTrueInt");

    double N_den = 0.0;
    double N_num = 0.0;

    while(reader.Next())
    {
        if(!(*HLT_Mu50 || *HLT_TkMu50)) continue;
        if(*MET_pt >= 65.0) continue;

        const int nMuon = Muon_pt.GetSize();

        int looseCount = 0;
        int looseIdx = -1;

        for(int i = 0; i < nMuon; ++i)
        {
            bool passLoose = Muon_pt[i] > 65.0 && std::fabs(Muon_eta[i]) < 2.4 && Muon_highPtId[i] == 2 && Muon_tkRelIso[i] < 0.40;

            if(passLoose)
            {
                ++looseCount;
                looseIdx = i;
            }
        }

        if(looseCount != 1 || looseIdx < 0) continue;

        bool extraMuon = false;

        for(int i = 0; i < nMuon; ++i)
        {
            if(i == looseIdx) continue;

            if(Muon_pt[i] > 20.0)
            {
                extraMuon = true;
                break;
            }
        }

        if(extraMuon) continue;

        int i = looseIdx;

        double pt = Muon_pt[i];
        double eta = Muon_eta[i];
        double abseta = std::fabs(eta);
        double phi = Muon_phi[i];

        if(pt < 65.0 || pt >= 100.0) continue;
        if(abseta >= 0.9) continue;

        double sf_id = GetSF(h2D_SF_ID, abseta, pt);
        double sf_iso = GetSF(h2D_SF_ISO, abseta, pt);
        double sf_trig = GetSF(h2D_SF_STRIG, abseta, pt);

        int bin_pu = hPUw->GetXaxis()->FindBin(*Pileup_nTrueInt);
        bin_pu = std::max(1, std::min(bin_pu, hPUw->GetNbinsX()));
        double wpu = hPUw->GetBinContent(bin_pu);

        double w = (double)(*genWeight) * sf_id * sf_iso * sf_trig * wpu;

        double dphiMuonMET = std::fabs(TVector2::Phi_mpi_pi(phi - *MET_phi));

        if(dphiMuonMET >= M_PI / 6.0) continue;
        if(Muon_dxyErr[i] <= 0.0) continue;

        double d0Significance = std::fabs(Muon_dxy[i]) / Muon_dxyErr[i];

        if(d0Significance <= 1.5) continue;

        bool hasBackToBackJet = false;

        for(int j = 0; j < Jet_pt.GetSize(); ++j)
        {
            if(Jet_pt[j] <= 30.0) continue;
            if(std::fabs(Jet_eta[j]) >= 2.5) continue;
            if(!(Jet_jetId[j] & 2)) continue;

            double dphiJetMuon = std::fabs(TVector2::Phi_mpi_pi(Jet_phi[j] - phi));

            if(dphiJetMuon > 5.0 * M_PI / 6.0)
            {
                hasBackToBackJet = true;
                break;
            }
        }

        if(!hasBackToBackJet) continue;

        h_fake_den_MET->Fill(*MET_pt, w);
        N_den += w;

        bool passTight = Muon_highPtId[i] == 2 && Muon_tkRelIso[i] < 0.10;

        if(passTight)
        {
            h_fake_num_MET->Fill(*MET_pt, w);
            N_num += w;
        }
    }

    TFile *fout = TFile::Open(outFile, "RECREATE");

    h_fake_den_MET->Write();
    h_fake_num_MET->Write();

    fout->Close();

    delete fout;
}
