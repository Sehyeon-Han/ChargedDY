#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TString.h>

void Muon_fake_rate_3D(const char* inFile,
                       const char* outFile)
{
    double ptBins[] = {65, 100, 150, 200, 300, 1500};
    const int N_ptBins = 5;

    double etaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4};
    const int N_etaBins = 4;

    double metBins[] = {0.0, 20.0, 40.0, 65.0};
    const int N_metBins = 3;
    
    TFile* fin = TFile::Open(inFile, "READ");

    if(!fin || fin->IsZombie())
    {
        std::cerr << "Cannot open input file: " << inFile << std::endl;

        return;
    }

    TH3D* h_fake_den_3D = new TH3D("h_fake_den_3D", "Fake-rate denominator;" "Muon p_{T} [GeV];" "|#eta|;" "MET [GeV]", N_ptBins, ptBins, N_etaBins, etaBins, N_metBins, metBins);
    TH3D* h_fake_num_3D = new TH3D("h_fake_num_3D", "Fake-rate numerator;" "Muon p_{T} [GeV];" "|#eta|;" "MET [GeV]", N_ptBins, ptBins, N_etaBins, etaBins, N_metBins, metBins);

    h_fake_den_3D->Sumw2();
    h_fake_num_3D->Sumw2();

    for(int ipt = 0; ipt < N_ptBins; ++ipt)
    {
        for(int ieta = 0; ieta < N_etaBins; ++ieta)
        {
            TString denName = Form("h_fake_den_MET_pt%d_eta%d", ipt + 1, ieta + 1);
            TString numName = Form("h_fake_num_MET_pt%d_eta%d", ipt + 1, ieta + 1);

            TH1D* hDen = dynamic_cast<TH1D*>(fin->Get(denName));
            TH1D* hNum = dynamic_cast<TH1D*>(fin->Get(numName));

            if(!hDen || !hNum)
            {
                std::cerr << "Missing histogram: " << denName << " or " << numName << std::endl;

                continue;
            }

            for(int imet = 0; imet < N_metBins; ++imet)
            {
                int bin1D = imet + 1;

                double den = hDen->GetBinContent(bin1D);
                double denErr = hDen->GetBinError(bin1D);
                double num = hNum->GetBinContent(bin1D);
                double numErr = hNum->GetBinError(bin1D);

                h_fake_den_3D->SetBinContent(ipt + 1, ieta + 1, imet + 1, den);
                h_fake_den_3D->SetBinError(ipt + 1, ieta + 1, imet + 1, denErr);
                h_fake_num_3D->SetBinContent(ipt + 1, ieta + 1, imet + 1, num);
                h_fake_num_3D->SetBinError(ipt + 1, ieta + 1, imet + 1, numErr);
            }
        }
    }

    TFile* fout = TFile::Open(outFile, "RECREATE");

    h_fake_den_3D->Write();
    h_fake_num_3D->Write();

    fout->Close();
    fin->Close();
}