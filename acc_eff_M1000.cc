#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TVector2.h>
#include <TH2D.h>

void acc_eff_M1000(const char* inFile,
                   const char* outFile)
{
    TH1::SetDefaultSumw2();

    double mTBins[] = {150, 200, 250, 300, 350, 425, 500, 600, 750, 900, 1100, 1400, 2000, 5000};
    const int N_mTBins = sizeof(mTBins) / sizeof(mTBins[0]) - 1;

    double mTBinsResp[] = {150, 200, 250, 300, 350, 425, 500, 600, 750, 900, 1100, 1400, 2000, 5000};
    const int N_mTBinsResp = 13;

    TH2D *h_response = new TH2D("h_response", "Response amtrix; GEN m_{T} [GeV]; Reco m_{T} [GeV]", N_mTBinsResp, mTBinsResp, N_mTBinsResp, mTBinsResp);

    TH1D *h_acc_den = new TH1D("h_acc_den", "acceptance denominator; m_{T}; Events", N_mTBins, mTBins);
    TH1D *h_acc_num_eff_den = new TH1D("h_acc_num", "acceptance numerator & efficiency denominator; m_{T}; Events", N_mTBins, mTBins);
    TH1D *h_eff_num = new TH1D("h_eff_num", "efficiency numerator; m_{T}; Events", N_mTBins, mTBins);

    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_tkRelIso(reader, "Muon_tkRelIso");
    TTreeReaderArray<UChar_t> Muon_highPtId(reader, "Muon_highPtId");
    TTreeReaderArray<Bool_t> Muon_looseId(reader, "Muon_looseId");
    TTreeReaderValue<Bool_t> HLT_Mu50(reader, "HLT_Mu50");
    TTreeReaderValue<Bool_t> HLT_TkMu50(reader, "HLT_TkMu50");
    TTreeReaderValue<Float_t> MET_pt(reader, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(reader, "MET_phi");
    TTreeReaderValue<Float_t> Pileup_nTrueInt(reader, "Pileup_nTrueInt");

    TTreeReaderArray<Float_t> GenPart_pt(reader, "GenPart_pt");
    TTreeReaderArray<Float_t> GenPart_phi(reader, "GenPart_phi");
    TTreeReaderArray<Float_t> GenPart_eta(reader, "GenPart_eta");
    TTreeReaderArray<Int_t> GenPart_pdgId(reader, "GenPart_pdgId");
    TTreeReaderValue<Float_t> GenMET_pt(reader, "GenMET_pt");
    TTreeReaderValue<Float_t> GenMET_phi(reader, "GenMET_phi");
    TTreeReaderArray<Int_t> GenPart_status(reader, "GenPart_status");
    TTreeReaderArray<Int_t> GenPart_statusFlags(reader, "GenPart_statusFlags");
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");
    TTreeReaderValue<Float_t> L1PreFiringWeight_Nom(reader, "L1PreFiringWeight_Nom");

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

    while(reader.Next())
    {
        bool pass_acc = false;

        int selIdx_gen = -1;
        double selected_gen_mT = -1.0;

        for(int i = 0; i < GenPart_pt.GetSize(); ++i)
        {
            double w_gen = double(*genWeight);

            if(abs(GenPart_pdgId[i]) != 13) continue;
            if(GenPart_status[i] != 1) continue;
            if(!(abs(GenPart_statusFlags[i]) & 256)) continue;

            selIdx_gen = i;
            
            double gen_dphi = TVector2::Phi_mpi_pi(GenPart_phi[i] - (*GenMET_phi));
            double gen_mT = std::sqrt(2.0 * GenPart_pt[i] * (*GenMET_pt) * (1.0 - std::cos(gen_dphi)));
            
            if(gen_mT >= 2000) continue;

            if(gen_mT < 150.0) continue;

            h_acc_den->Fill(gen_mT, w_gen);
            
            if(*GenMET_pt <= 85.0) continue;

            if(GenPart_pt[i] > 65.0 && std::fabs(GenPart_eta[i]) < 2.4)
            {
                selIdx_gen = i;
                selected_gen_mT = gen_mT;
                h_acc_num_eff_den->Fill(selected_gen_mT, w_gen);
                pass_acc = true;
                break;
            }
        }

        if(pass_acc)
        {
            if(!(*HLT_Mu50 || *HLT_TkMu50)) continue;
            if(*MET_pt <= 85.0) continue;

            int selCount_reco = 0;
            int selIdx_reco = -1;

            bool ExtraMuon = false;
            
            for(int i = 0; i < Muon_pt.GetSize(); ++i)
            {
                if(Muon_pt[i] > 65.0 && std::fabs(Muon_eta[i]) < 2.4 && Muon_tkRelIso[i] < 0.10 && Muon_highPtId[i] == 2)
                {
                    ++selCount_reco;
                    selIdx_reco = i;
                }
            }

            if(selCount_reco == 1 && selIdx_reco >= 0)
            {
                for(int i = 0; i < Muon_pt.GetSize(); ++i)
                {
                    if(i == selIdx_reco) continue;

                    if(Muon_pt[i] > 20.0 && Muon_looseId[i] && std::fabs(Muon_eta[i]) < 2.4)
                    {
                        ExtraMuon = true;
                        break;
                    }
                }
            }

            if(selCount_reco == 1 && selIdx_reco >= 0 && !ExtraMuon)
            {
                double pt = Muon_pt[selIdx_reco];
                double aeta = fabs(Muon_eta[selIdx_reco]);

                double sf_id = GetSF(h2D_SF_ID, aeta, pt);
                double sf_iso = GetSF(h2D_SF_ISO, aeta, pt);
                double sf_trig = GetSF(h2D_SF_STRIG, aeta, pt);

                int bin_pu = hPUw->GetXaxis()->FindBin(*Pileup_nTrueInt);
                bin_pu = std::max(1, std::min(bin_pu, hPUw->GetNbinsX()));
                double wpu = hPUw->GetBinContent(bin_pu);

                double w = double(*genWeight) * sf_id * sf_iso * sf_trig * wpu * (*L1PreFiringWeight_Nom);

                double reco_dphi = TVector2::Phi_mpi_pi(Muon_phi[selIdx_reco] - (*MET_phi));
                double reco_mT = std::sqrt(2.0 * Muon_pt[selIdx_reco] * (*MET_pt) * (1.0 - std::cos(reco_dphi)));

                if(reco_mT < 150.0) continue;

                h_eff_num->Fill(selected_gen_mT, w);
                h_response->Fill(selected_gen_mT, reco_mT, w);
            }
        }
    }

    TFile fout(outFile, "RECREATE");

    h_acc_den->Write();
    h_acc_num_eff_den->Write();
    h_eff_num->Write();
    h_response->Write();

    fout.Close();
}