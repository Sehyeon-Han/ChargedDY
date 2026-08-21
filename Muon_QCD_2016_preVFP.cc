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

void Muon_QCD_2016_preVFP(const char* inFile,
                          const char* outFile,
                          const char* realRateFile = "real_mu_rate.root",
                          const char* fakeRateFile = "fake_mu_normalized_results.root")
{
    TH1::SetDefaultSumw2();

    double mTBins[] = {200, 250, 300, 350, 400, 500, 600, 700, 800, 1000, 1500, 2000, 3500};
    int nmTBins = sizeof(mTBins) / sizeof(double) - 1;

    TH1D *h_Muon_pt = new TH1D("h_Muon_pt", "QCD;Muon p_{T} [GeV];Events", 30, 65, 1500);
    TH1D *h_Muon_eta = new TH1D("h_Muon_eta", "QCD;Muon #eta;Events", 50, -3, 3);
    TH1D *h_Muon_phi = new TH1D("h_Muon_phi", "QCD;Muon #phi;Events", 50, -3.3, 3.3);
    TH1D *h_MET_pt = new TH1D("h_MET_pt", "QCD;MET p_{T} [GeV];Events", 30, 85, 3500);
    TH1D *h_MET_phi = new TH1D("h_MET_phi", "QCD;MET #phi;Events", 50, -3.3, 3.3);
    TH1D *h_mT = new TH1D("h_mT", "QCD;m_{T} [GeV];Events", nmTBins, mTBins);
    TH1D *h_highMT_NT = new TH1D("h_high_mT_NT", "m_{T} >= 1500 GeV; Category; Events", 1, 0.0, 1.0);
    TH1D *h_highMT_NL = new TH1D("h_high_mT_NL", "m_{T} >= 1500 GeV; Category; Events", 1, 0.0, 1.0);

    TFile *fReal = TFile::Open(realRateFile, "READ");
    TFile *fFake = TFile::Open(fakeRateFile, "READ");

    if(!fReal || fReal->IsZombie())
    {
        std::cerr << "Cannot open real rate file: " << realRateFile << std::endl;
        return;
    }

    if(!fFake || fFake->IsZombie())
    {
        std::cerr << "Cannot open fake rate file: " << fakeRateFile << std::endl;
        fReal->Close();
        return;
    }

    TH2D *h_real_rate = dynamic_cast<TH2D*>(fReal->Get("h_real_rate"));
    TH2D *h_fake_rate = dynamic_cast<TH2D*>(fFake->Get("h_fake_rate"));

    if(!h_real_rate)
    {
        std::cerr << "Cannot find h_real_rate in " << realRateFile << std::endl;
        fReal->Close();
        fFake->Close();
        return;
    }

    if(!h_fake_rate)
    {
        std::cerr << "Cannot find h_fake_rate in " << fakeRateFile << std::endl;
        fReal->Close();
        fFake->Close();
        return;
    }

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

    if(chain.GetEntries() == 0)
    {
        std::cerr << "No entries found in input: " << inFile << std::endl;
        fReal->Close();
        fFake->Close();
        return;
    }

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

    Long64_t totalEvents = 0;
    Long64_t triggerEvents = 0;
    Long64_t selectedEvents = 0;
    Long64_t tightEvents = 0;
    Long64_t looseNotTightEvents = 0;
    Long64_t invalidRateEvents = 0;
    Long64_t finalEvents = 0;
    Long64_t highMT_TightEvents = 0;
    Long64_t highMT_LooseNotTightEvents = 0;

    double tightWeightSum = 0.0;
    double looseWeightSum = 0.0;
    double highMT_TightWeightSum = 0.0;
    double highMT_LooseWeightSum = 0.0;

    while(reader.Next())
    {
        ++totalEvents;

        if(!(*HLT_Mu50 || *HLT_TkMu50)) continue;
        ++triggerEvents;

        if(*MET_pt <= 85.0) continue;

        int looseCount = 0;
        int looseIdx = -1;

        const int nMuon = Muon_pt.GetSize();

        for(int i = 0; i < nMuon; ++i)
        {
            bool passLoose =
                Muon_pt[i] > 65.0 &&
                std::fabs(Muon_eta[i]) < 2.4 &&
                Muon_highPtId[i] == 2 &&
                Muon_iso[i] < 0.40;

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

        ++selectedEvents;

        double muonPt = Muon_pt[looseIdx];
        double muonEta = Muon_eta[looseIdx];
        double muonPhi = Muon_phi[looseIdx];
        double muonIso = Muon_iso[looseIdx];

        bool passTight = muonIso < 0.10;
        bool passLooseNotTight = muonIso >= 0.10 && muonIso < 0.40;

        if(!passTight && !passLooseNotTight) continue;

        double epsilonR = GetRate(h_real_rate, std::fabs(muonEta), muonPt);
        double epsilonF = GetRate(h_fake_rate, std::fabs(muonEta), muonPt);

        double denominator = epsilonR - epsilonF;

        if(epsilonR <= 0.0 || epsilonR > 1.0 ||
           epsilonF < 0.0 || epsilonF >= 1.0 ||
           std::fabs(denominator) < 1.0e-6)
        {
            ++invalidRateEvents;
            continue;
        }

        double weight = 0.0;

        if(passTight)
        {
            weight = -epsilonF * (1.0 - epsilonR) / denominator;
        }
        else
        {
            weight = epsilonF * epsilonR / denominator;
        }

        double dphi = TVector2::Phi_mpi_pi(muonPhi - *MET_phi);
        double mT = std::sqrt(2.0 * muonPt * (*MET_pt) * (1.0 - std::cos(dphi)));

        if(mT <= 200.0) continue;

        ++finalEvents;

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

        if(mT >= 1000.0)
        {
            if(passTight)
            {
                ++highMT_TightEvents;
                highMT_TightWeightSum += weight;
                h_highMT_NT->Fill(0.5);

            }
            else if(passLooseNotTight)
            {
                ++highMT_LooseNotTightEvents;
                highMT_LooseWeightSum += weight;
                h_highMT_NL->Fill(0.5);
            }
        }
            
        h_Muon_pt->Fill(muonPt, weight);
        h_Muon_eta->Fill(muonEta, weight);
        h_Muon_phi->Fill(muonPhi, weight);
        h_MET_pt->Fill(*MET_pt, weight);
        h_MET_phi->Fill(*MET_phi, weight);
        h_mT->Fill(mT, weight);
    }

    TFile *fOut = TFile::Open(outFile, "RECREATE");

    if(!fOut || fOut->IsZombie())
    {
        std::cerr << "Cannot create output file: " << outFile << std::endl;
        fReal->Close();
        fFake->Close();
        return;
    }

    h_Muon_pt->Write();
    h_Muon_eta->Write();
    h_Muon_phi->Write();
    h_MET_pt->Write();
    h_MET_phi->Write();
    h_mT->Write();
    h_highMT_NT->Write();
    h_highMT_NL->Write();

    fOut->Close();
    fReal->Close();
    fFake->Close();

    delete fOut;
    delete fReal;
    delete fFake;
}