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

void Draw_Data_MC_gen_2016_preVFP(const double xsec200 = 8.78288653e00,
                                  const double xsec500 = 2.76542446e-01,
                                  const double xsec1000 = 1.58702653e-02,
                                  const double xsec2000 = 4.15066636e-04,
                                  const double lumi = 19520.0,
                                  const double sumW200 = 1250000,
                                  const double sumW500 = 750000,
                                  const double sumW1000 = 250000,
                                  const double sumW2000 = 50000)
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    TFile *fMC200 = TFile::Open("Muon_MC_gen_200_2016_preVFP.root");
    TFile *fMC500 = TFile::Open("Muon_MC_gen_500_2016_preVFP.root");
    TFile *fMC1000 = TFile::Open("Muon_MC_gen_1000_2016_preVFP.root");
    TFile *fMC2000 = TFile::Open("Muon_MC_gen_2000_2016_preVFP.root");

    TH1D *hMC200 = (TH1D*)fMC200->Get("h_mW");
    TH1D *hMC500 = (TH1D*)fMC500->Get("h_mW");
    TH1D *hMC1000 = (TH1D*)fMC1000->Get("h_mW");
    TH1D *hMC2000 = (TH1D*)fMC2000->Get("h_mW");

    hMC200->SetDirectory(0);
    hMC500->SetDirectory(0);
    hMC1000->SetDirectory(0);
    hMC2000->SetDirectory(0);

    fMC200->Close();
    fMC500->Close();
    fMC1000->Close();
    fMC2000->Close();

    double kW200 = (xsec200 * lumi) / sumW200;
    double kW500 = (xsec500 * lumi) / sumW500;
    double kW1000 = (xsec1000 * lumi) / sumW1000;
    double kW2000 = (xsec2000 * lumi) / sumW2000;

    hMC200->Scale(kW200);
    hMC500->Scale(kW500);
    hMC1000->Scale(kW1000);
    hMC2000->Scale(kW2000);

    TH1D *hMCsum = (TH1D*)hMC200->Clone("hMCsum");
    hMCsum->Add(hMC500);
    hMCsum->Add(hMC1000);
    hMCsum->Add(hMC2000);

    TCanvas *c = new TCanvas("c", "m_{W}", 900, 800);
    c->SetLogx();
    c->SetLogy();

    THStack *stack = new THStack("hstack", "m_{W} ;m_{W} [GeV]; Events");

    hMC200->SetLineColor(kGreen-2);
    hMC200->SetFillColor(kGreen-2);
    hMC500->SetLineColor(kGreen-1);
    hMC500->SetFillColor(kGreen-1);
    hMC1000->SetLineColor(kGreen);
    hMC1000->SetFillColor(kGreen);
    hMC2000->SetLineColor(kGreen+1);
    hMC2000->SetFillColor(kGreen+1);

    stack->Add(hMC200);
    stack->Add(hMC500);
    stack->Add(hMC1000);
    stack->Add(hMC2000);
    
    stack->Draw("HIST");

    auto leg = new TLegend(0.65, 0.72, 0.88, 0.88);
    leg->AddEntry(hMC200, "W#rightarrow#tau#nu", "f");

    c->Update();
    c->SaveAs("MC_gen_2016_preVFP_mW.png");
}