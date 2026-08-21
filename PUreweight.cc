#include <iostream>
#include <algorithm>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void PUreweight(const char* outFile = "PUreweight.root")
{
    gStyle->SetOptStat(0);

    const int N = 99;

    double MC2016_99[N] = {1.004024e-05, 5.764988e-05, 7.378914e-05, 1.109329e-04, 1.588577e-04, 3.686374e-04, 8.931141e-04, 1.897008e-03, 3.588802e-03,
                           6.360526e-03, 1.041740e-02, 1.581226e-02, 2.237857e-02, 2.991869e-02, 3.802759e-02, 4.543139e-02, 5.111811e-02, 5.474346e-02,
                           5.679062e-02, 5.771455e-02, 5.781769e-02, 5.712516e-02, 5.554565e-02, 5.313438e-02, 5.015190e-02, 4.668158e-02, 4.292446e-02,
                           3.895668e-02, 3.485072e-02, 3.073569e-02, 2.677121e-02, 2.297202e-02, 1.933887e-02, 1.596025e-02, 1.293105e-02, 1.028887e-02,
                           7.987828e-03, 6.066517e-03, 4.478209e-03, 3.215898e-03, 2.245042e-03, 1.514474e-03, 9.811837e-04, 6.096705e-04, 3.621934e-04,
                           2.115726e-04, 1.191524e-04, 6.491335e-05, 3.577958e-05, 1.990436e-05, 1.136393e-05, 6.496241e-06, 3.966262e-06, 2.379102e-06,
                           1.509974e-06, 1.098167e-06, 7.312985e-07, 6.103988e-07, 3.748458e-07, 2.651773e-07, 2.019235e-07, 1.393476e-07, 8.326001e-08,
                           6.049324e-08, 6.525366e-08, 5.905746e-08, 2.291625e-08, 1.972946e-08, 1.773110e-08, 3.575479e-09, 1.350398e-09, 8.500712e-09,
                           5.027919e-09, 4.937367e-10, 8.139197e-10, 5.627789e-09, 5.151406e-10, 8.216767e-10, 0.000000e-00, 1.491669e-09, 8.435180e-09,
                           0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00,
                           0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00, 0.000000e-00};

    double Data_2016_PreAPV[99] = {2.484602e-06, 2.224640e-05, 6.247184e-05, 1.084867e-04, 1.890218e-04, 4.346442e-04, 1.191636e-03, 2.467397e-03, 4.368343e-03,   
                                   7.481333e-03, 1.192522e-02, 1.742427e-02, 2.391670e-02, 3.155836e-02, 4.019179e-02, 4.873285e-02, 5.560076e-02, 6.024768e-02,   
                                   6.281161e-02, 6.354123e-02, 6.287627e-02, 6.122218e-02, 5.876006e-02, 5.550184e-02, 5.144212e-02, 4.667708e-02, 4.143100e-02,   
                                   3.601126e-02, 3.072492e-02, 2.580533e-02, 2.138478e-02, 1.750618e-02, 1.415330e-02, 1.128370e-02, 8.851128e-03, 6.814557e-03,   
                                   5.138221e-03, 3.787902e-03, 2.727272e-03, 1.916746e-03, 1.314788e-03, 8.804059e-04, 5.757299e-04, 3.678706e-04, 2.298147e-04,   
                                   1.404583e-04, 8.403983e-05, 4.925671e-05, 2.829831e-05, 1.594632e-05, 8.820875e-06, 4.794867e-06, 2.565141e-06, 1.353503e-06,   
                                   7.065553e-07, 3.664015e-07, 1.897361e-07, 9.870732e-08, 5.191258e-08, 2.775178e-08, 1.513368e-08, 8.425940e-09, 4.781222e-09,   
                                   2.754836e-09, 1.604379e-09, 9.402281e-10, 5.523667e-10, 3.243561e-10, 1.899794e-10, 1.108293e-10, 6.433391e-11, 3.713331e-11,   
                                   2.130094e-11, 1.213815e-11, 6.868220e-12, 3.857391e-12, 2.149397e-12, 1.187751e-12, 6.506233e-13, 3.531353e-13, 1.898346e-13,   
                                   1.010309e-13, 5.321180e-14, 2.772533e-14, 1.428600e-14, 7.277305e-15, 3.663773e-15, 1.822498e-15, 8.955267e-16, 4.345751e-16,   
                                   2.082285e-16, 9.849703e-17, 4.598754e-17, 2.118991e-17, 9.634240e-18, 4.322469e-18, 1.913089e-18, 8.349962e-19, 3.595710e-19};

    TH1D* hMC = new TH1D("hPU_MC", "PU distribution; Pileup nTrueInt; Probability [%]", N, -0.5, N - 0.5);
    TH1D* hData = new TH1D("hPU_Data", "PU distribution; Pileup nTrueInt; Probability [%]", N, -0.5, N - 0.5);

    for(int i = 0; i < N; ++i)
    {
       hMC->SetBinContent(i + 1, MC2016_99[i]);
       hData->SetBinContent(i + 1, Data_2016_PreAPV[i]);
    }

    TH1D* hW = (TH1D*)hData->Clone("hPUw");
    hW->SetTitle("PU weight (Data/MC); Pileup nTrueInt; Weight");
    hW->Divide(hMC);
    hMC->Scale(100.0);
    hData->Scale(100.0);

    TFile f(outFile, "RECREATE");
    hMC->Write();
    hData->Write();
    hW->Write();
    f.Close();

    TCanvas* c = new TCanvas("cPU", "PU Data vs MC", 900, 800);

    TPad* pad1 = new TPad("pad1","",0,0.30,1,1);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    hMC->GetXaxis()->SetLabelSize(0);
    hMC->GetXaxis()->SetTitleSize(0);
    hMC->GetXaxis()->SetTickLength(0);

    hData->GetXaxis()->SetLabelSize(0);
    hData->GetXaxis()->SetTitleSize(0);
    hData->GetXaxis()->SetTickLength(0);

    hMC->SetLineWidth(2);
    hData->SetLineWidth(2);
    hMC->SetLineColor(kRed);
    hData->SetLineColor(kBlack);

    double ymax = std::max(hMC->GetMaximum(), hData->GetMaximum());
    hMC->SetMaximum(1.25*ymax);

    hMC->Draw("HIST");
    hData->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.62,0.72,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hData,"Data PU","l");
    leg->AddEntry(hMC,"MC PU","l");
    leg->Draw();

    c->cd();

    TPad* pad2 = new TPad("pad2","",0,0,1,0.30);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.30);
    pad2->Draw();
    pad2->cd();

    hW->SetLineColor(kBlue);
    hW->SetLineWidth(2);
    hW->SetTitle("");
    hW->GetYaxis()->SetTitle("Data/MC");
    hW->GetYaxis()->SetNdivisions(505);
    hW->GetYaxis()->SetTitleSize(0.10);
    hW->GetYaxis()->SetLabelSize(0.09);
    hW->GetYaxis()->SetTitleOffset(0.45);

    hW->GetXaxis()->SetTitle("Pileup nTrueInt");
    hW->GetXaxis()->SetTitleSize(0.12);
    hW->GetXaxis()->SetLabelSize(0.10);

    hW->SetMinimum(0.0);
    hW->SetMaximum(3.0);  
    hW->Draw("HIST");

    TLine* line = new TLine(0, 1, N-1, 1);
    line->SetLineStyle(2);
    line->Draw("SAME");

    c->SaveAs("PU_Data_vs_MC_with_ratio.png");

    std::cout << "Saved: PU_Data_vs_MC_with_ratio.png\n";
}