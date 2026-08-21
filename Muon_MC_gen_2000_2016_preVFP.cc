#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TLorentzVector.h>

void Muon_MC_gen_2000_2016_preVFP(const char* inFile,
                                  const char* outFile = "Muon_MC_gen_2000_2016_preVFP.root")
{
    TH1::SetDefaultSumw2();

    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<Int_t> GenPart_pdgId(reader, "GenPart_pdgId");
    TTreeReaderArray<Int_t> GenPart_status(reader, "GenPart_status");
    TTreeReaderArray<Float_t> GenPart_pt(reader, "GenPart_pt");
    TTreeReaderArray<Float_t> GenPart_eta(reader, "GenPart_eta");
    TTreeReaderArray<Float_t> GenPart_phi(reader, "GenPart_phi");
    TTreeReaderArray<Float_t> GenPart_mass(reader, "GenPart_mass");
    TTreeReaderArray<Int_t> GenPart_statusFlags(reader, "GenPart_statusFlags");
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");

    TH1D *h_mW = new TH1D("h_mW", "m_{W}; m_{W} [GeV]; Events", 33, 200, 3500);

    while(reader.Next())
    {
        int muIdx = -1;
        int nuIdx = -1;

        int nMu = 0;
        int nNu = 0;

        for(int i = 0; i < GenPart_pdgId.GetSize(); ++i)
        {
            int absId = std::abs(GenPart_pdgId[i]);
            
            //if(!(GenPart_status[i] == 2)) continue;
            if(!(GenPart_statusFlags[i] & 128)) continue;

            if(absId == 15)
            {
                nMu++;
                muIdx = i;
            }

            if(absId == 16)
            {
                nNu++;
                nuIdx = i;
            }
        }

        if(!(nMu == 1 && nNu == 1)) continue;

        TLorentzVector mu, nu;

        mu.SetPtEtaPhiM(GenPart_pt[muIdx], GenPart_eta[muIdx], GenPart_phi[muIdx], GenPart_mass[muIdx]);

        nu.SetPtEtaPhiM(GenPart_pt[nuIdx], GenPart_eta[nuIdx], GenPart_phi[nuIdx], GenPart_mass[nuIdx]);

        double mW = (mu + nu).M();

        h_mW->Fill(mW, *genWeight);
    }

    TFile* fout = new TFile(outFile, "RECREATE");
    h_mW->Write();
    fout->Close();
}