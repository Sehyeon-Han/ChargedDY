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

void Muon_real_rate(const char* inFile,
                    const char* outFile)
{
    TStopwatch timer;
    timer.Start();

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    double ptBins[] = {65, 100, 150, 200, 300, 1500};
    const int N_ptBins = 5;
    
    double etaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4};
    const int N_etaBins = 4;

    TH2D *h_real_den = new TH2D("h_real_den", "Real-rate denominator; |#eta|; p_{T} [GeV]", N_etaBins, etaBins, N_ptBins, ptBins);
    TH2D *h_real_num = new TH2D("h_real_num", "Real-rate numerator; |#eta|; p_{T} [GeV]", N_etaBins, etaBins, N_ptBins, ptBins);

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
    TTreeReaderArray<UChar_t> Muon_highPtId(reader, "Muon_highPtId");
    TTreeReaderArray<Int_t> Muon_genPartIdx(reader, "Muon_genPartIdx");
    TTreeReaderArray<Float_t> Muon_tkRelIso(reader, "Muon_tkRelIso");
    TTreeReaderValue<Bool_t> HLT_TkMu50(reader, "HLT_TkMu50");
    TTreeReaderValue<Bool_t> HLT_Mu50(reader, "HLT_Mu50");
    TTreeReaderArray<Int_t> GenPart_pdgId(reader, "GenPart_pdgId");
    TTreeReaderArray<Int_t> GenPart_status(reader, "GenPart_status");
    TTreeReaderArray<Int_t> GenPart_statusFlags(reader, "GenPart_statusFlags");
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");
    TTreeReaderValue<Float_t> L1PreFiringWeight_Nom(reader, "L1PreFiringWeight_Nom");
    TTreeReaderValue<Float_t> Pileup_nTrueInt(reader, "Pileup_nTrueInt");

    double N_L = 0.0;
    double N_T = 0.0;

    while(reader.Next())
    {
        for(int i = 0; i < Muon_pt.GetSize(); ++i)
        {
            if(!(*HLT_TkMu50 || *HLT_Mu50)) continue;

            double pt = Muon_pt[i];
            double aeta = std::fabs(Muon_eta[i]);

            double sf_id = GetSF(h2D_SF_ID, aeta, pt);
            double sf_iso = GetSF(h2D_SF_ISO, aeta, pt);
            double sf_trig = GetSF(h2D_SF_STRIG, aeta, pt);

            int bin_pu = hPUw->GetXaxis()->FindBin(*Pileup_nTrueInt);
            bin_pu = std::max(1, std::min(bin_pu, hPUw->GetNbinsX()));
            double wpu = hPUw->GetBinContent(bin_pu);

            double w = (double)(*genWeight) * sf_id * sf_iso * sf_trig * wpu * (*L1PreFiringWeight_Nom);

            if(pt < 65.0) continue;
            if(aeta >= 2.4) continue;

            int genIdx = Muon_genPartIdx[i];

            if(genIdx < 0) continue;

            if(!(std::abs(GenPart_pdgId[genIdx]) == 13)) continue;
            if(!(GenPart_status[genIdx] == 1)) continue;
            if(!(GenPart_statusFlags[genIdx] & 256)) continue;

            double fillPt = pt;

            if(fillPt >= 1500.0)
            fillPt = 1499.999;

            bool passLoose = (Muon_highPtId[i] == 2 && Muon_tkRelIso[i] < 0.40);
            bool passTight = (Muon_highPtId[i] == 2 && Muon_tkRelIso[i] < 0.10);

            if(!passLoose) continue;
            h_real_den->Fill(aeta, pt);
            N_L += w;

            if(passTight)
            {
                h_real_num->Fill(aeta, pt);
                N_T += w;
            }
        }
    }

    TFile *fout = TFile::Open(outFile, "RECREATE");

    h_real_den->Write();
    h_real_num->Write();

    fout->Close();

    timer.Stop();

    std::cout << "Real time: " << timer.RealTime() << " s" << std::endl;
    std::cout << "CPU time: " << timer.CpuTime() << " s" << std::endl;
}