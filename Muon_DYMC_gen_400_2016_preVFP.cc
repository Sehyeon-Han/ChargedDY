#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TLorentzVector.h>

void Muon_DYMC_gen_400_2016_preVFP(const char* inFile = "/pnfs/knu.ac.kr/data/cms/store/user/sungwon/DY_Run2_UL_NanoAOD/2016_preVFP/DYJetsToMuMu_M-400to500/*/*/*/*/*.root",
                                   const char* outFile = "Muon_DYMC_gen_400_2016_preVFP.root")
{
    TH1::SetDefaultSumw2();

    TChain chain("Events");
    chain.Add(inFile);

    std::cout << "Entries in chain = " << chain.GetEntries() << std::endl;

    TTreeReader reader(&chain);

    TTreeReaderArray<Int_t> LHEPart_pdgId(reader, "LHEPart_pdgId");
    TTreeReaderArray<Float_t> LHEPart_pt(reader, "LHEPart_pt");
    TTreeReaderArray<Float_t> LHEPart_eta(reader, "LHEPart_eta");
    TTreeReaderArray<Float_t> LHEPart_phi(reader, "LHEPart_phi");
    TTreeReaderArray<Float_t> LHEPart_mass(reader, "LHEPart_mass");
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");

    TH1D* h_mZ = new TH1D("h_mZ", "m_{Z}; m_{Z} [GeV]; Events", 33, 200, 3500);

    Long64_t nEvents = 0;
    Long64_t nFilled = 0;
    Long64_t nHasPair = 0;

    while(reader.Next())
    {
        nEvents++;

        TLorentzVector muMinus;
        TLorentzVector muPlus;

        bool hasMuMinus = false;
        bool hasMuPlus  = false;

        for(unsigned int i = 0; i < LHEPart_pdgId.GetSize(); ++i)
        {
            int pdgId = LHEPart_pdgId[i];

            if(pdgId == 13)
            {
                muMinus.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i], LHEPart_mass[i]);
                hasMuMinus = true;
            }

            if(pdgId == -13)
            {
                muPlus.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i], LHEPart_mass[i]);
                hasMuPlus = true;
            }
        }

        if(!(hasMuMinus && hasMuPlus)) continue;

        nHasPair++;

        TLorentzVector Z = muMinus + muPlus;
        double mZ = Z.M();

        h_mZ->Fill(mZ, *genWeight);
        nFilled++;
    }

    TFile* fout = new TFile(outFile, "RECREATE");

    h_mZ->Write();
    fout->Close();
}