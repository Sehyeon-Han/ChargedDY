#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>

void Draw_PU()
{
    gStyle->SetOptStat(0);

    TFile *f = TFile::Open("PUreweight.root", "READ");

    if(!f || f->IsZombie())
    {
        std::cerr << "Cannot open PUreweight.root" << std::endl;
        return;
    }

    TH1D *hW = (TH1D*)f->Get("hPU_MC");
    TH1D *hData = (TH1D*)f->Get("hPU_Data");
    TH1D *hPUw = (TH1D*)f->Get("hPUw");

    if(!hW || !hData || !hPUw)
    {
        std::cerr << "Cannot find hW, hData, or hPUw" << std::endl;
        f->ls();
        return;
    }

    TH1D *hW_before = (TH1D*)hW->Clone("hW_before");
    TH1D *hW_after = (TH1D*)hW->Clone("hW_after");
    TH1D *hData_draw = (TH1D*)hData->Clone("hData_draw");

    hW_before->SetDirectory(0);
    hW_after->SetDirectory(0);
    hData_draw->SetDirectory(0);

    // PU-corrected MC distribution
    hW_after->Multiply(hPUw);

    // Shape comparison
    if(hW_before->Integral() > 0)
    {
        hW_before->Scale(1.0 / hW_before->Integral());
    }

    if(hW_after->Integral() > 0)
    {
        hW_after->Scale(1.0 / hW_after->Integral());
    }

    if(hData_draw->Integral() > 0)
    {
        hData_draw->Scale(1.0 / hData_draw->Integral());
    }

    hW_before->SetLineColor(kBlue);
    hW_before->SetLineWidth(2);

    hW_after->SetLineColor(kRed);
    hW_after->SetLineWidth(4);

    hData_draw->SetLineColor(kBlack);
    hData_draw->SetLineWidth(2);
    hData_draw->SetMarkerColor(kBlack);
    hData_draw->SetMarkerStyle(20);
    hData_draw->SetMarkerSize(0.8);

    double maximum = std::max({
        hW_before->GetMaximum(),
        hW_after->GetMaximum(),
        hData_draw->GetMaximum()
    });

    TCanvas *c = new TCanvas("c", "Pileup correction", 800, 700);

    hData_draw->SetTitle(
        "Pileup distribution;Number of true interactions;Normalized events"
    );

    hData_draw->SetMaximum(1.25 * maximum);
    hData_draw->SetMinimum(0.0);

    hData_draw->Draw("HIST");
    hW_before->Draw("HIST SAME");
    hW_after->Draw("HIST SAME");
    hData_draw->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.55, 0.68, 0.88, 0.88);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(hData_draw, "Data", "lep");
    leg->AddEntry(hW_before, "W MC before PU correction", "l");
    leg->AddEntry(hW_after, "W MC after PU correction", "l");

    leg->Draw();

    c->SaveAs("Pileup_correction.png");
    c->SaveAs("Pileup_correction.pdf");

    f->Close();
}