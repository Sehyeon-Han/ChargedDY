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

void Draw_Data_DYMC_gen_preVFP(const double DYxsec200 = 2.78,
                               const double DYxsec400 = 0.15,
                               const double DYxsec500 = 0.084,
                               const double DYxsec700 = 0.013,
                               const double DYxsec800 = 0.011,
                               const double DYxsec1000 = 0.006,
                               const double DYxsec1500 = 0.00081,
                               const double DYxsec2000 = 0.0002,
                               const double lumi = 19520.0,
                               const double sumWDY200 = 2673760.516799,
                               const double sumWDY400 = 44344.277472,
                               const double sumWDY500 = 24718.574857,
                               const double sumWDY700 = 3951.644255,
                               const double sumWDY800 = 3345.560043,
                               const double sumWDY1000 = 1782.846123,
                               const double sumWDY1500 = 244.385548,
                               const double sumWDY2000 = 62.650483)
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    TFile *fDYMC200 = TFile::Open("Muon_DYMC_gen_200_2016_preVFP.root");
    TFile *fDYMC400 = TFile::Open("Muon_DYMC_gen_400_2016_preVFP.root");
    TFile *fDYMC500 = TFile::Open("Muon_DYMC_gen_500_2016_preVFP.root");
    TFile *fDYMC700 = TFile::Open("Muon_DYMC_gen_700_2016_preVFP.root");
    TFile *fDYMC800 = TFile::Open("Muon_DYMC_gen_800_2016_preVFP.root");
    TFile *fDYMC1000 = TFile::Open("Muon_DYMC_gen_1000_2016_preVFP.root");
    TFile *fDYMC1500 = TFile::Open("Muon_DYMC_gen_1500_2016_preVFP.root");
    TFile *fDYMC2000 = TFile::Open("Muon_DYMC_gen_2000_2016_preVFP.root");

    TH1D *hDYMC200 = (TH1D*)fDYMC200->Get("h_mZ");
    TH1D *hDYMC400 = (TH1D*)fDYMC400->Get("h_mZ");
    TH1D *hDYMC500 = (TH1D*)fDYMC500->Get("h_mZ");
    TH1D *hDYMC700 = (TH1D*)fDYMC700->Get("h_mZ");
    TH1D *hDYMC800 = (TH1D*)fDYMC800->Get("h_mZ");
    TH1D *hDYMC1000 = (TH1D*)fDYMC1000->Get("h_mZ");
    TH1D *hDYMC1500 = (TH1D*)fDYMC1500->Get("h_mZ");
    TH1D *hDYMC2000 = (TH1D*)fDYMC2000->Get("h_mZ");

    hDYMC200->SetDirectory(0);
    hDYMC400->SetDirectory(0);
    hDYMC500->SetDirectory(0);
    hDYMC700->SetDirectory(0);
    hDYMC800->SetDirectory(0);
    hDYMC1000->SetDirectory(0);
    hDYMC1500->SetDirectory(0);
    hDYMC2000->SetDirectory(0);

    fDYMC200->Close();
    fDYMC400->Close();
    fDYMC500->Close();
    fDYMC700->Close();
    fDYMC800->Close();
    fDYMC1000->Close();
    fDYMC1500->Close();
    fDYMC2000->Close();

    double kDY200 = (DYxsec200 * lumi) / sumWDY200;
    double kDY400 = (DYxsec400 * lumi) / sumWDY400;
    double kDY500 = (DYxsec500 * lumi) / sumWDY500;
    double kDY700 = (DYxsec700 * lumi) / sumWDY700;
    double kDY800 = (DYxsec800 * lumi) / sumWDY800;
    double kDY1000 = (DYxsec1000 * lumi) / sumWDY1000;
    double kDY1500 = (DYxsec1500 * lumi) / sumWDY1500;
    double kDY2000 = (DYxsec2000 * lumi) / sumWDY2000;

    hDYMC200->Scale(kDY200);
    hDYMC400->Scale(kDY400);
    hDYMC500->Scale(kDY500);
    hDYMC700->Scale(kDY700);
    hDYMC800->Scale(kDY800);
    hDYMC1000->Scale(kDY1000);
    hDYMC1500->Scale(kDY1500);
    hDYMC2000->Scale(kDY2000);

    TH1D *hMCsum = (TH1D*)hDYMC200->Clone("hMCsum");
    hMCsum->Add(hDYMC400);
    hMCsum->Add(hDYMC500);
    hMCsum->Add(hDYMC700);
    hMCsum->Add(hDYMC800);
    hMCsum->Add(hDYMC1000);
    hMCsum->Add(hDYMC1500);
    hMCsum->Add(hDYMC2000);

    TCanvas *c = new TCanvas("c", "m_{Z}", 900, 800);
    //c->SetLogx();
    c->SetLogy();

    THStack *stack = new THStack("hstack", "m_{Z} ;m_{Z} [GeV];Events");

    hDYMC200->SetLineColor(kGreen-2);
    hDYMC200->SetFillColor(kGreen-2);
    hDYMC400->SetLineColor(kGreen-1);
    hDYMC400->SetFillColor(kGreen-1);
    hDYMC500->SetLineColor(kGreen);
    hDYMC500->SetFillColor(kGreen);
    hDYMC700->SetLineColor(kGreen+1);
    hDYMC700->SetFillColor(kGreen+1);
    hDYMC800->SetLineColor(kGreen+2);
    hDYMC800->SetFillColor(kGreen+2);
    hDYMC1000->SetLineColor(kRed-1);
    hDYMC1000->SetFillColor(kRed-1);
    hDYMC1500->SetLineColor(kRed);
    hDYMC1500->SetFillColor(kRed);
    hDYMC2000->SetLineColor(kRed+1);
    hDYMC2000->SetFillColor(kRed+1);

    stack->Add(hDYMC200);
    stack->Add(hDYMC400);
    stack->Add(hDYMC500);
    stack->Add(hDYMC700);
    stack->Add(hDYMC800);
    stack->Add(hDYMC1000);
    stack->Add(hDYMC1500);
    stack->Add(hDYMC2000);
    stack->Draw("HIST");

    auto leg = new TLegend(0.65, 0.72, 0.88, 0.88);
    leg->AddEntry(hDYMC200, "Z#rightarrow#mu#mu", "f");

    c->Update();
    c->SaveAs("Data_DYMC_gen_2016_preVFP_mZ.png");
}