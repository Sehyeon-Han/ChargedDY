#include <iostream>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <algorithm>

void Differential_cross_section()
{
    gStyle->SetOptStat(0);

    const int N = 12;

    double mTBins[] = {200, 250, 300, 350, 425, 500, 600, 750, 900, 1100, 1400, 2000, 5000};

    double x[N];

    for(int i = 0; i < N; ++i)
    {
        x[i] = std::sqrt(mTBins[i] * mTBins[i + 1]);
    }

    double ct18[N] = {0.0338, 0.013977, 0.0066322, 0.0030649, 0.0013568, 0.00061298, 0.00023586, 8.6252e-05, 3.2062e-05, 9.8179e-06, 1.8304e-06, 5.5968e-08};

    double result[N] = {0.0514, 0.016536, 0.0072009, 0.0030620, 0.0012422, 0.00047624, 0.00015065, 3.9952e-05, 2.9539e-05, 1.05825e-05, 1.8130e-06, 9.4555e-08};

    TH1D* hCT18   = new TH1D("hCT18", "", N, mTBins);
    TH1D* hResult = new TH1D("hResult", "", N, mTBins);

    for(int i = 0; i < N; ++i)
    {
        hCT18->SetBinContent(i + 1, ct18[i]);
        hResult->SetBinContent(i + 1, result[i]);

        hCT18->SetBinError(i + 1, 0.0);
        hResult->SetBinError(i + 1, 0.0);
    }

    hCT18->SetLineColor(kRed);
    hCT18->SetLineWidth(2);

    hResult->SetLineColor(kBlack);
    hResult->SetLineWidth(2);

    TCanvas* c = new TCanvas("c", "c", 900, 800);

    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.18, 1, 1);
    pad1->SetLeftMargin(0.13);
    pad1->SetBottomMargin(0.02);
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    hCT18->GetYaxis()->SetTitle("d#sigma/dm_{T}");
    hCT18->GetYaxis()->SetTitleSize(0.05);
    hCT18->GetYaxis()->SetLabelSize(0.04);
    hCT18->GetYaxis()->SetTitleOffset(1.1);

    hCT18->GetXaxis()->SetLabelSize(0);
    hCT18->GetXaxis()->SetTitleSize(0);

    hCT18->SetMinimum(1e-9);
    hCT18->SetMaximum(1e-1);

    hCT18->Draw("HIST");
    hResult->SetMarkerStyle(20);
    hResult->SetMarkerSize(1.0);
    hResult->SetMarkerColor(kBlack);
    hResult->SetLineColor(kBlack);

    hResult->Draw("P SAME");

    TLegend* leg = new TLegend(0.62, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hCT18, "CT18NNLO", "l");
    leg->AddEntry(hResult, "Data W#rightarrow#mu#nu", "P");
    leg->Draw();

    c->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.18);
    pad2->SetLeftMargin(0.13);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.38);
    pad2->SetLogx();
    pad2->Draw();
    pad2->cd();

    TH1D* hRatio = (TH1D*)hResult->Clone("hRatio");
    hRatio->Divide(hCT18);

    hRatio->GetYaxis()->SetTitle("Data / Theory");
    hRatio->GetYaxis()->SetNdivisions(303);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetTitleOffset(0.5);

    hRatio->GetXaxis()->SetTitle("m_{T} [GeV]");
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);

    hRatio->SetMinimum(0.75);
    hRatio->SetMaximum(1.25);

    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.8);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);

    hRatio->Draw("P");

    TLine* line = new TLine(200, 1.0, 5000, 1.0);
    line->SetLineStyle(2);
    line->Draw("SAME");

    c->SaveAs("CT18NNLO_data_ratio.png");
}