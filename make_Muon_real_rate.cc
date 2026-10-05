#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

void make_Muon_real_rate(const char* inFile = "real_mu_merged.root",
                         const char* outFile = "real_mu_rate.root")
{
    gStyle->SetOptStat(0);

    TFile *fin = TFile::Open(inFile, "READ");

    TH2D *h_real_den = (TH2D*)fin->Get("h_real_den");
    TH2D *h_real_num = (TH2D*)fin->Get("h_real_num");
    TH2D *h_real_rate = (TH2D*)h_real_num->Clone("h_real_rate");

    h_real_rate->SetDirectory(nullptr);
    h_real_rate->Divide(h_real_den);
    h_real_rate->SetTitle("Muon real rate; |#eta|; p_{T} [GeV]");

    for(int ix = 1; ix <= h_real_den->GetNbinsX(); ++ix)
    {
        for(int iy = 1; iy <= h_real_den->GetNbinsY(); ++iy)
        {
            double den = h_real_den->GetBinContent(ix, iy);
            double num = h_real_num->GetBinContent(ix, iy);
            double rate = h_real_rate->GetBinContent(ix, iy);
            double error = h_real_rate->GetBinError(ix, iy);

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

    h_real_den->Write();
    h_real_num->Write();
    h_real_rate->Write();

    TCanvas *c = new TCanvas("c", "Muon real rate", 900, 800);

    h_real_rate->SetMinimum(0.94);
    h_real_rate->SetMaximum(1.0);
    h_real_rate->Draw("COLZ TEXT");

    c->Write();
    c->SaveAs("real_mu_rate.png");

    fout->Close();
    fin->Close();
}