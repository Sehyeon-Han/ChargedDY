#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

void make_Muon_fake_rate(const char* inFile = "fake_mu_merged.root",
                         const char* outFile = "fake_mu_rate.root")
{
    gStyle->SetOptStat(0);

    TFile *fin = TFile::Open(inFile, "READ");

    TH2D *h_fake_den = (TH2D*)fin->Get("h_fake_den");
    TH2D *h_fake_num = (TH2D*)fin->Get("h_fake_num");
    TH2D *h_fake_rate = (TH2D*)h_fake_num->Clone("h_fake_rate");

    h_fake_rate->SetDirectory(nullptr);
    h_fake_rate->Divide(h_fake_den);
    h_fake_rate->SetTitle("Muon fake rate; |#eta|; p_{T} [GeV]");

    std::cout << "\nMuon fake-rate bin contents\n";

    for(int ix = 1; ix <= h_fake_den->GetNbinsX(); ++ix)
    {
        for(int iy = 1; iy <= h_fake_den->GetNbinsY(); ++iy)
        {
            double den = h_fake_den->GetBinContent(ix, iy);
            double num = h_fake_num->GetBinContent(ix, iy);
            double rate = h_fake_rate->GetBinContent(ix, iy);
            double error = h_fake_rate->GetBinError(ix, iy);

            std::cout
                << "pt bin = " << iy
                << ", eta bin = " << ix
                << ", denominator = " << den
                << ", numerator = " << num
                << ", fake rate = " << rate
                << " +/- " << error
                << std::endl;
        }
    }

    TFile *fout = TFile::Open(outFile, "RECREATE");

    h_fake_den->Write();
    h_fake_num->Write();
    h_fake_rate->Write();

    TCanvas *c = new TCanvas("c", "Muon fake rate", 900, 800);

    h_fake_rate->SetMinimum(0.0);
    h_fake_rate->SetMaximum(1.0);
    h_fake_rate->Draw("COLZ TEXT");

    c->Write();
    c->SaveAs("fake_mu_rate.png");

    fout->Close();
    fin->Close();
}