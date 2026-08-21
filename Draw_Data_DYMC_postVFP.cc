#include <iostream>
#include <algorithm>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>

void Draw_Data_DYMC_postVFP(const double DYxsec50 = 6019.95,
                            const double lumi = 16810.0,
                            const double sumWDY50 = 3.3902455098858924e20)
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    TFile *fDYMC50 = TFile::Open("Muon_DYMC_50_2016_postVFP.root");

    TH1D *hDYMC50 = (TH1D*)fDYMC50->Get("h_mT");

    hDYMC50->SetDirectory(0);

    fDYMC50->Close();

    double kDY50 = (DYxsec50 * lumi) / sumWDY50;

    hDYMC50->Scale(kDY50);

    TCanvas *c = new TCanvas("c", "Transverse Mass", 900, 800);

    hDYMC50->Draw();

    c->cd();

    c->Update();

    c->SaveAs("Data_DYMC_2016_postVFP_mT.png");
}