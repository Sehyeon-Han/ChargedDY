#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <cmath>
#include <algorithm>
#include <TVector2.h>

void Muon_QCD_CR_2016_preVFP(const char* inFile,
                             const char* outFile,
                             const char* realRateFile = "real_mu_rate.root",
                             const char* fakeRateFile = "fake_mu_normalized_results.root")
{
    TH1::SetDefaultSumw2();

    double ptBins[] = {65, 100, 150, 200, 300, 1500};
    const int N_ptBins = 5;

    TH1D *h_Muon_pt_CR = new TH1D("h_Muon_pt", "QCD;Muon p_{T} [GeV];Events", 30, 65, 1500);
    TH1D *h_Muon_eta_CR = new TH1D("h_Muon_eta", "QCD;Muon #eta;Events", 50, -3, 3);
    TH1D *h_Muon_phi_CR = new TH1D("h_Muon_phi", "QCD;Muon #phi;Events", 50, -3.5, 3.5);
    TH1D *h_MET_pt_CR = new TH1D("h_MET_pt", "QCD;MET p_{T} [GeV];Events", 40, 0, 500);
    TH1D *h_MET_phi_CR = new TH1D("h_MET_phi", "QCD;MET #phi;Events", 50, -3.5, 3.5);

    TFile *fReal = TFile::Open(realRateFile, "READ");
    TFile *fFake = TFile::Open(fakeRateFile, "READ");

    TH2D *h_real_rate = dynamic_cast<TH2D*>(fReal->Get("h_real_rate"));
    TH2D *h_fake_rate = dynamic_cast<TH2D*>(fFake->Get("h_fake_rate"));

    auto GetRate = [](TH2D *hist, double absEta, double pt)
    {
        double etaValue = std::clamp(absEta,
                                     hist->GetXaxis()->GetXmin() + 1.0e-6,
                                     hist->GetXaxis()->GetXmax() - 1.0e-6);

        double ptValue = std::clamp(pt,
                                    hist->GetYaxis()->GetXmin() + 1.0e-6,
                                    hist->GetYaxis()->GetXmax() - 1.0e-6);

        int etaBin = hist->GetXaxis()->FindFixBin(etaValue);
        int ptBin = hist->GetYaxis()->FindFixBin(ptValue);

        return hist->GetBinContent(etaBin, ptBin);
    };

    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<UChar_t> Muon_highPtId(reader, "Muon_highPtId");
    TTreeReaderArray<Float_t> Muon_tkRelIso(reader, "Muon_tkRelIso");
    TTreeReaderValue<Bool_t> HLT_TkMu50(reader, "HLT_TkMu50");
    TTreeReaderValue<Bool_t> HLT_Mu50(reader, "HLT_Mu50");
    TTreeReaderArray<Float_t> Muon_dxy(reader, "Muon_dxy");
    TTreeReaderArray<Float_t> Muon_dxyErr(reader, "Muon_dxyErr");
    TTreeReaderArray<Float_t> Jet_pt(reader, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_phi(reader, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_eta(reader, "Jet_eta");
    TTreeReaderArray<Int_t> Jet_jetId(reader, "Jet_jetId");
    TTreeReaderValue<Float_t> MET_pt(reader, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(reader, "MET_phi");

    Long64_t tightEvents = 0;
    Long64_t looseNotTightEvents = 0;

    double tightWeightSum = 0.0;
    double looseWeightSum = 0.0;
    while(reader.Next())
    {
        if(!(*HLT_Mu50 || *HLT_TkMu50)) continue;
        if(*MET_pt >= 65.0) continue;

        int looseCount = 0;
        int looseIdx = -1;

        const int nMuon = Muon_pt.GetSize();

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

        double muonPt = Muon_pt[looseIdx];
        double muonEta = Muon_eta[looseIdx];
        double muonPhi = Muon_phi[looseIdx];
        double muonIso = Muon_tkRelIso[looseIdx];

        bool passTight = muonIso < 0.10;
        bool passLooseNotTight = muonIso >= 0.10 && muonIso < 0.40;

        if(!passTight && !passLooseNotTight) continue;

        double epsilonR = GetRate(h_real_rate, std::fabs(muonEta), muonPt);
        double epsilonF = GetRate(h_fake_rate, std::fabs(muonEta), muonPt);

        double denominator = epsilonR - epsilonF;

        double weight = 0.0;

        if(passTight)
        {
            weight = -epsilonF * (1.0 - epsilonR) / denominator;
        }
        else
        {
            weight = epsilonF * epsilonR / denominator;
        }
        
        if(passTight)
        {
            ++tightEvents;
            tightWeightSum += weight;
        }
        else
        {
            ++looseNotTightEvents;
            looseWeightSum += weight;
        }

        h_Muon_pt_CR->Fill(muonPt, weight);
        h_Muon_eta_CR->Fill(muonEta, weight);
        h_Muon_phi_CR->Fill(muonPhi, weight);
        h_MET_pt_CR->Fill(*MET_pt, weight);
        h_MET_phi_CR->Fill(*MET_phi, weight);
    }

     TFile *fOut = TFile::Open(outFile, "RECREATE");

    h_Muon_pt_CR->Write();
    h_Muon_eta_CR->Write();
    h_Muon_phi_CR->Write();
    h_MET_pt_CR->Write();
    h_MET_phi_CR->Write();

    fOut->Close();
    fReal->Close();
    fFake->Close();

    delete fOut;
    delete fReal;
    delete fFake;
}