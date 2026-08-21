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

void Muon_fake_rate(const char* inFile,
                    const char* outFile)
{
    gStyle->SetOptStat(0);

    double ptBins[] = {65, 100, 150, 200, 300, 1500};
    const int N_ptBins = 5;

    double etaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4};
    const int N_etaBins = 4;

    TH2D *h_fake_den = new TH2D("h_fake_den", "Fake-rate denominator; |#eta|; p_{T} [GeV]", N_etaBins, etaBins, N_ptBins, ptBins);
    TH2D *h_fake_num = new TH2D("h_fake_num", "Fake-rate numerator; |#eta|; p_{T} [GeV]", N_etaBins, etaBins, N_ptBins, ptBins);

    TH1D *h_Muon_pt_den = new TH1D("h_Muon_pt_den", "Muon p_{T}; p_{T} [GeV]; Events", 30, 65, 1500);
    TH1D *h_Muon_eta_den = new TH1D("h_Muon_eta_den", "Muon #eta; #eta; Events", 50, -3, 3);
    TH1D *h_MET_pt_den = new TH1D("h_MET_pt_den", "MET p_{T}; MET p_{T} [GeV]; Events", 40, 0, 500);
    TH1D *h_Muon_phi_den = new TH1D("h_Muon_phi_den", "Muon #phi; #phi; Events", 50, -3.5, 3.5);
    TH1D *h_MET_phi_den = new TH1D("h_MET_phi_den", "MET #phi; MET #phi; Events", 50, -3.5, 3.5);
    TH1D *h_Muon_tkRelIso_den = new TH1D("h_Muon_tkRelIso_den", "Muon tkRelIso; iso; Events", 30, 0, 0.45);

    TH1D *h_Muon_pt_num = new TH1D("h_Muon_pt_num", "Muon p_{T}; p_{T} [GeV]; Events", 30, 65, 1500);
    TH1D *h_Muon_eta_num = new TH1D("h_Muon_eta_num", "Muon #eta; #eta; Events", 50, -3, 3);
    TH1D *h_MET_pt_num = new TH1D("h_MET_pt_num", "MET p_{T}; MET p_{T} [GeV]; Events", 40, 0, 500);
    TH1D *h_Muon_phi_num = new TH1D("h_Muon_phi_num", "Muon #phi; #phi; Events", 50, -3.5, 3.5);
    TH1D *h_MET_phi_num = new TH1D("h_MET_phi_num", "MET #phi; MET #phi; Events", 50, -3.5, 3.5);
    TH1D *h_Muon_tkRelIso_num = new TH1D("h_Muon_tkRelIso_num", "Muon tkRelIso; iso; Events", 30, 0, 0.45);

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

    double N_looseInclusive = 0.0;
    double N_tight = 0.0;

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

        bool passTight = Muon_highPtId[i] == 2 && Muon_tkRelIso[i] < 0.10;

        double ratePt = std::min(pt, ptBins[N_ptBins] - 1.0e-6);

        h_fake_den->Fill(abseta, ratePt);
        N_looseInclusive += 1.0;

        if(passTight)
        {
            h_fake_num->Fill(abseta, ratePt);
            N_tight += 1.0;

            h_Muon_pt_num->Fill(ratePt);
            h_Muon_eta_num->Fill(eta);
            h_Muon_phi_num->Fill(phi);
            h_MET_pt_num->Fill(*MET_pt);
            h_MET_phi_num->Fill(*MET_phi);
            h_Muon_tkRelIso_num->Fill(Muon_tkRelIso[i]);
        }
        
        h_Muon_pt_den->Fill(ratePt);
        h_Muon_eta_den->Fill(eta);
        h_Muon_phi_den->Fill(phi);
        h_MET_pt_den->Fill(*MET_pt);
        h_MET_phi_den->Fill(*MET_phi);
        h_Muon_tkRelIso_den->Fill(Muon_tkRelIso[i]);
    }

    TFile *fout = TFile::Open(outFile, "RECREATE");

    h_fake_den->Write();
    h_fake_num->Write();

    h_Muon_pt_den->Write();
    h_Muon_eta_den->Write();
    h_Muon_phi_den->Write();
    h_MET_pt_den->Write();
    h_MET_phi_den->Write();
    h_Muon_tkRelIso_den->Write();

    h_Muon_pt_num->Write();
    h_Muon_eta_num->Write();
    h_Muon_phi_num->Write();
    h_MET_pt_num->Write();
    h_MET_phi_num->Write();
    h_Muon_tkRelIso_num->Write();

    fout->Close();
}
