#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TString.h>

void Data_WMC(const double xsec200 = 8.78288653e00,
              const double xsec500 = 2.76542446e-01,
              const double xsec1000 = 1.58702653e-02,
              const double xsec2000 = 4.15066636e-04,
              const double lumi = 19520.0,
              const double sumW200 = 3250000,
              const double sumW500 = 750000,
              const double sumW1000 = 250000,
              const double sumW2000 = 50000)

{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    TFile *fData = TFile::Open("Muon_data_2016_preVFP.root");
    TFile *fMC200 = TFile::Open("Muon_MC_200_2016_preVFP.root");
    TFile *fMC500 = TFile::Open("Muon_MC_500_2016_preVFP.root");
    TFile *fMC1000 = TFile::Open("Muon_MC_1000_2016_preVFP.root");
    TFile *fMC2000 = TFile::Open("Muon_MC_2000_2016_preVFP.root");

    TH1D *hData = (TH1D*)fData->Get("h_mT");
    TH1D *hMC200 = (TH1D*)fMC200->Get("h_mT");
    TH1D *hMC500 = (TH1D*)fMC500->Get("h_mT");
    TH1D *hMC1000 = (TH1D*)fMC1000->Get("h_mT");
    TH1D *hMC2000 = (TH1D*)fMC2000->Get("h_mT");

    hData->SetDirectory(0);
    hMC200->SetDirectory(0);
    hMC500->SetDirectory(0);
    hMC1000->SetDirectory(0);
    hMC2000->SetDirectory(0);

    fData->Close();
    fMC200->Close();
    fMC500->Close();
    fMC1000->Close();
    fMC2000->Close();

    double kMC200 = (xsec200 * lumi) / sumW200;
    double kMC500 = (xsec500 * lumi) / sumW500;
    double kMC1000 = (xsec1000 * lumi) / sumW1000;
    double kMC2000 = (xsec2000 * lumi) / sumW2000;

    hMC200->Scale(kMC200);
    hMC500->Scale(kMC500);
    hMC1000->Scale(kMC1000);
    hMC2000->Scale(kMC2000);

    TH1D *hMCsum = (TH1D*)hMC200->Clone("hMCsum");
    hMCsum->Add(hMC500);
    hMCsum->Add(hMC1000);
    hMCsum->Add(hMC2000);

    TCanvas *c = new TCanvas("c", "m_{T}", 900, 600);

    THStack *stack = new THStack("hstack", "m_{T}; m_{T} [GeV]; Events");

    hMC200->SetFillColor(kGreen+1);
    hMC200->SetLineColor(kGreen+1);
    hMC500->SetFillColor(kGreen+1);
    hMC500->SetLineColor(kGreen+1);
    hMC1000->SetFillColor(kGreen+1);
    hMC1000->SetLineColor(kGreen+1);
    hMC2000->SetFillColor(kGreen+1);
    hMC2000->SetLineColor(kGreen+1);

    stack->Add(hMC200);
    stack->Add(hMC500);
    stack->Add(hMC1000);
    stack->Add(hMC2000);
    stack->Draw("HIST");

    hData->SetLineColor(kBlack);
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.0);
    hData->Draw("E1 SAME");

    c->Update();
}