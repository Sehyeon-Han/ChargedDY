#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>
#include <cmath>

void Draw_FakeRate_pt_eta()
{
    gStyle->SetOptStat(0);

    const int Npt  = 5;
    const int Neta = 4;

    // pT bin edges
    double ptLow[Npt]  = {65, 100, 150, 200, 300};
    double ptHigh[Npt] = {100, 150, 200, 300, 1500};

    // x position: geometric center for log-x plot
    double x[Npt];
    double exl[Npt];
    double exh[Npt];

    for (int i = 0; i < Npt; i++) {
        x[i]   = std::sqrt(ptLow[i] * ptHigh[i]);
        exl[i] = x[i] - ptLow[i];
        exh[i] = ptHigh[i] - x[i];
    }

    // y errors = 0 (vertical error bars 안 그림)
    double eyl[Npt] = {0, 0, 0, 0, 0};
    double eyh[Npt] = {0, 0, 0, 0, 0};

    // =====================================================
    // 여기에 네 fake rate 값 넣으면 됨
    // fakeRate[eta bin][pt bin]
    // eta 0: 0.0 < |eta| < 0.9
    // eta 1: 0.9 < |eta| < 1.2
    // eta 2: 1.2 < |eta| < 2.1
    // eta 3: 2.1 < |eta| < 2.4
    // =====================================================

    double fakeRate[Neta][Npt] = {
        {0.992, 0.994, 0.993, 0.988, 0.942},  // eta bin 1
        {0.993, 0.994, 0.993, 0.990, 0.948},  // eta bin 2
        {0.994, 0.995, 0.994, 0.992, 0.969},  // eta bin 3
        {0.995, 0.996, 0.996, 0.994, 0.989}   // eta bin 4
    };

    TCanvas *c = new TCanvas("c", "Fake Rate", 800, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.08);
    c->SetLogx();

    TGraphAsymmErrors *gr[Neta];

    int colors[Neta]  = {kBlue+1, kOrange+7, kRed+1, kYellow+3};
    int markers[Neta] = {20, 21, 22, 23};

    for (int i = 0; i < Neta; i++) {
        gr[i] = new TGraphAsymmErrors(Npt, x, fakeRate[i], exl, exh, eyl, eyh);

        gr[i]->SetMarkerStyle(markers[i]);
        gr[i]->SetMarkerSize(1.2);
        gr[i]->SetMarkerColor(colors[i]);
        gr[i]->SetLineColor(colors[i]);
        gr[i]->SetLineWidth(2);
    }

    // 첫 번째 그래프로 축 생성
    gr[0]->SetTitle("");
    gr[0]->Draw("AP");

    gr[0]->GetXaxis()->SetTitle("Muon p_{T} (GeV)");
    gr[0]->GetYaxis()->SetTitle("Real rate");

    gr[0]->GetXaxis()->SetLimits(55, 2500);
    gr[0]->GetYaxis()->SetRangeUser(0.9, 1.05);

    gr[0]->GetXaxis()->SetTitleSize(0.055);
    gr[0]->GetYaxis()->SetTitleSize(0.055);
    gr[0]->GetXaxis()->SetLabelSize(0.045);
    gr[0]->GetYaxis()->SetLabelSize(0.045);
    gr[0]->GetYaxis()->SetTitleOffset(1.1);

    // 나머지 overlay
    for (int i = 1; i < Neta; i++) {
        gr[i]->Draw("P SAME");
    }

    // 범례
    TLegend *leg = new TLegend(0.58, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(gr[0], "0.0 < |#eta| < 0.9", "lp");
    leg->AddEntry(gr[1], "0.9 < |#eta| < 1.2", "lp");
    leg->AddEntry(gr[2], "1.2 < |#eta| < 2.1", "lp");
    leg->AddEntry(gr[3], "2.1 < |#eta| < 2.4", "lp");
    leg->Draw();

    c->SaveAs("FakeRate_pt_eta.png");
}