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

void Draw_Data_MC_postVFP(const double xsec200 = 8.78288653e00,
                          const double xsec500 = 2.76542446e-01,
                          const double xsec1000 = 1.58702653e-02,
                          const double xsec2000 = 4.15066636e-04,
                          const double DYxsec200 = 2.78,
                          const double lumi = 16810.0,
                          const double sumW200 = 3250000,
                          const double sumW500 = 750000,
                          const double sumW1000 = 250000,
                          const double sumW2000 = 50000,
                          const double sumWDY200 = 2673760.516799)
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    TFile *fData = TFile::Open("Muon_data_2016_postVFP.root");
    TFile *fMC200 = TFile::Open("Muon_MC_200_2016_postVFP.root");
    TFile *fMC500 = TFile::Open("Muon_MC_500_2016_postVFP.root");
    TFile *fMC1000 = TFile::Open("Muon_MC_1000_2016_postVFP.root");
    TFile *fMC2000 = TFile::Open("Muon_MC_2000_2016_postVFP.root");
    TFile *fDYMC200 = TFile::Open("Muon_DYMC_200_2016_postVFP.root");

    TH1D *hData = (TH1D*)fData->Get("h_mT");
    TH1D *hMC200 = (TH1D*)fMC200->Get("h_mT");
    TH1D *hMC500 = (TH1D*)fMC500->Get("h_mT");
    TH1D *hMC1000 = (TH1D*)fMC1000->Get("h_mT");
    TH1D *hMC2000 = (TH1D*)fMC2000->Get("h_mT");
    TH1D *hDYMC200 = (TH1D*)fDYMC200->Get("h_mT");

    hData->SetDirectory(0);
    hMC200->SetDirectory(0);
    hMC500->SetDirectory(0);
    hMC1000->SetDirectory(0);
    hMC2000->SetDirectory(0);
    hDYMC200->SetDirectory(0);

    fData->Close();
    fMC200->Close();
    fMC500->Close();
    fMC1000->Close();
    fMC2000->Close();
    fDYMC200->Close();

    double kW200 = (xsec200 * lumi) / sumW200;
    double kW500 = (xsec500 * lumi) / sumW500;
    double kW1000 = (xsec1000 * lumi) / sumW1000;
    double kW2000 = (xsec2000 * lumi) / sumW2000;
    double kDY200 = (DYxsec200 * lumi) / sumWDY200;

    hMC200->Scale(kW200);
    hMC500->Scale(kW500);
    hMC1000->Scale(kW1000);
    hMC2000->Scale(kW2000);
    hDYMC200->Scale(kDY200);

    TH1D *hMCsum = (TH1D*)hMC200->Clone("hMCsum");
    hMCsum->Add(hMC500);
    hMCsum->Add(hMC1000);
    hMCsum->Add(hMC2000);
    hMCsum->Add(hDYMC200);

    const double ymax = std::max(hMCsum->GetMaximum(), hData->GetMaximum());

    TCanvas *c = new TCanvas("c", "Transverse Mass", 900, 800);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.18, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    THStack *stack = new THStack("hstack", "Transverse Mass;m_{T} [GeV];Events");

    hMC200->SetFillColor(kRed);
    hMC200->SetLineColor(kRed);
    hMC500->SetFillColor(kRed);
    hMC500->SetLineColor(kRed);
    hMC1000->SetLineColor(kRed);
    hMC1000->SetFillColor(kRed);
    hMC2000->SetLineColor(kRed);
    hMC2000->SetFillColor(kRed);
    hDYMC200->SetLineColor(kGreen);
    hDYMC200->SetFillColor(kGreen);

    stack->Add(hDYMC200);
    stack->Add(hMC200);
    stack->Add(hMC500);
    stack->Add(hMC1000);
    stack->Add(hMC2000);
    stack->Draw("HIST");

    stack->SetMinimum(1e-1);
    stack->SetMaximum(10.0 * ymax);

    hData->SetLineColor(kBlack);
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.0);
    hData->Draw("E1 SAME");

    stack->GetXaxis()->SetLabelSize(0);
    stack->GetXaxis()->SetTitleSize(0);
    stack->GetXaxis()->SetTickLength(0);

    auto leg = new TLegend(0.65, 0.72, 0.88, 0.88);
    leg->AddEntry(hData, "Data", "lep");
    leg->AddEntry(hMC200, "W#rightarrow#mu#nu", "f");
    leg->AddEntry(hDYMC200, "Z#rightarrow#mu#mu", "f");
    leg->Draw();

    pad1->Update();

    c->cd();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.18);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.38);
    pad2->SetLogx();
    pad2->Draw();
    pad2->cd();

    TH1D *hRatio = (TH1D*)hData->Clone("hRatio");
    hRatio->SetTitle("");
    hRatio->Divide(hMCsum);

    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.8);
    hRatio->SetLineColor(kBlack);

    hRatio->GetYaxis()->SetTitle("Data / MC");
    hRatio->GetYaxis()->SetNdivisions(303);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetTitleOffset(0.45);

    hRatio->GetXaxis()->SetTitle("m_{T} [GeV]");
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);

    hRatio->GetXaxis()->SetLabelSize(0);

    hRatio->SetMinimum(0.0);
    hRatio->SetMaximum(2.0);
    hRatio->Draw("E1");

    double y_min = hRatio->GetMinimum();
    double y_max = hRatio->GetMaximum();

    std::vector<double> marks = {500, 1000, 2000};

    TLatex latex;
    latex.SetTextSize(0.11);
    latex.SetTextFont(42);
    latex.SetTextAlign(23);
    double y_text = -0.35;

    for(double x : marks)
    {
        TLine *l = new TLine(x, y_min, x, y_max);
        l->SetLineStyle(3);
        l->SetLineColor(kGray);
        l->Draw("SAME");

        if (x == 500)
            latex.DrawLatex(x, y_text, "500");
        else if (x == 1000)
            latex.DrawLatex(x, y_text, "10^{3}");
        else if (x == 2000)
            latex.DrawLatex(x, y_text, "2#times10^{3}");
    }

    TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                            hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    pad2->Update();
    c->Update();

    c->SaveAs("Data_MC_2016_postVFP_mT.png");
}