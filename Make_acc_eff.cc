#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>

void Make_acc_eff()
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    const int N = 4;

    const char* files[N] = {"acc_eff_M200.root", "acc_eff_M500.root", "acc_eff_M1000.root", "acc_eff_M2000.root"};
    double xsec[N] = {8.78288653e00, 2.76542446e-01, 1.58702653e-02, 4.15066636e-04};
    double sumW[N] = {3250000, 750000, 250000, 50000};
    const double lumi = 19520.0;

    TH1D* h_acc_den_comb = nullptr;
    TH1D* h_acc_num_comb = nullptr;
    TH1D* h_eff_num_comb = nullptr;

    for(int i = 0; i < N; ++i)
    {
        TFile* f = TFile::Open(files[i]);

        if(!f || f->IsZombie())
        {
            std::cerr << "Cannot open " << files[i] << std::endl;
            continue;
        }

        TH1D* h_acc_den = (TH1D*)f->Get("h_acc_den");
        TH1D* h_acc_num = (TH1D*)f->Get("h_acc_num");
        TH1D* h_eff_num = (TH1D*)f->Get("h_eff_num");

        double norm = xsec[i] * lumi / sumW[i];

        TH1D* h_den_tmp = (TH1D*)h_acc_den->Clone(Form("h_acc_den_%d", i));
        TH1D* h_acc_tmp = (TH1D*)h_acc_num->Clone(Form("h_acc_num_%d", i));
        TH1D* h_eff_tmp = (TH1D*)h_eff_num->Clone(Form("h_eff_num_%d", i));

        h_den_tmp->SetDirectory(nullptr);
        h_acc_tmp->SetDirectory(nullptr);
        h_eff_tmp->SetDirectory(nullptr);

        h_den_tmp->Scale(norm);
        h_acc_tmp->Scale(norm);
        h_eff_tmp->Scale(norm);

        if(!h_acc_den_comb)
        {
            h_acc_den_comb = (TH1D*)h_den_tmp->Clone("h_acc_den_combined");
            h_acc_num_comb = (TH1D*)h_acc_tmp->Clone("h_acc_num_combined");
            h_eff_num_comb = (TH1D*)h_eff_tmp->Clone("h_eff_num_combined");

            h_acc_den_comb->SetDirectory(nullptr);
            h_acc_num_comb->SetDirectory(nullptr);
            h_eff_num_comb->SetDirectory(nullptr);
        }
        else
        {
            h_acc_den_comb->Add(h_den_tmp);
            h_acc_num_comb->Add(h_acc_tmp);
            h_eff_num_comb->Add(h_eff_tmp);
        }

        delete h_den_tmp;
        delete h_acc_tmp;
        delete h_eff_tmp;

        f->Close();
        delete f;
    }

    if(!h_acc_den_comb)
    {
        std::cerr << "No histograms loaded." << std::endl;
        return;
    }

    TH1D* h_acc = (TH1D*)h_acc_num_comb->Clone("h_acceptance");
    h_acc->SetTitle("Acceptance;m_{T}^{gen} [GeV];Acceptance");
    h_acc->Divide(h_acc_den_comb);

    TH1D* h_eff = (TH1D*)h_eff_num_comb->Clone("h_efficiency");
    h_eff->SetTitle("Efficiency;m_{T}^{gen} [GeV]; Efficiency");
    h_eff->Divide(h_acc_num_comb);

    TH1D* h_acc_eff = (TH1D*)h_eff_num_comb->Clone("h_acceptance_efficiency");
    h_acc_eff->SetTitle("Acceptance #times Efficiency;" "m_{T}^{gen} [GeV]; A#times#epsilon");
    h_acc_eff->Divide(h_acc_den_comb);

    std::cout << "\n====================================\n";
    std::cout << "       Acceptance / Efficiency\n";
    std::cout << "====================================\n";

    for(int i = 1; i <= h_acc->GetNbinsX(); ++i)
    {
        double low = h_acc->GetXaxis()->GetBinLowEdge(i);

        double high = h_acc->GetXaxis()->GetBinUpEdge(i);

        std::cout
            << "[" << low
            << ", " << high << "]"
            << "  A = "
            << h_acc->GetBinContent(i)
            << " +- "
            << h_acc->GetBinError(i)

            << "  eff = "
            << h_eff->GetBinContent(i)
            << " +- "
            << h_eff->GetBinError(i)

            << "  A*eff = "
            << h_acc_eff->GetBinContent(i)
            << " +- "
            << h_acc_eff->GetBinError(i)

            << std::endl;
    }

    TFile fout("acceptance_efficiency.root", "RECREATE");

    h_acc_den_comb->Write();
    h_acc_num_comb->Write();
    h_eff_num_comb->Write();

    h_acc->Write();
    h_eff->Write();
    h_acc_eff->Write();

    fout.Close();

    TCanvas* c1 = new TCanvas("c1", "", 800, 800);

    TPad* pad1_acc = new TPad("pad1_acc", "", 0.0, 0.30, 1.0, 1.0);
    TPad* pad2_acc = new TPad("pad2_acc", "", 0.0, 0.00, 1.0, 0.30);

    pad1_acc->SetBottomMargin(0.02);
    pad2_acc->SetTopMargin(0.03);
    pad2_acc->SetBottomMargin(0.30);

    pad1_acc->Draw();
    pad2_acc->Draw();

    pad1_acc->cd();
    pad1_acc->SetLogx();
    pad1_acc->SetLogy();

    h_acc_den_comb->SetTitle("");
    h_acc_den_comb->SetLineColor(kBlack);
    h_acc_den_comb->SetLineWidth(2);

    h_acc_num_comb->SetLineColor(kRed);
    h_acc_num_comb->SetLineWidth(2);

    h_acc_den_comb->GetYaxis()->SetTitle("Events");
    h_acc_den_comb->GetXaxis()->SetLabelSize(0);

    h_acc_den_comb->Draw("HIST E");
    h_acc_num_comb->Draw("HIST E SAME");

    TLegend* leg_acc = new TLegend(0.60, 0.72, 0.88, 0.88);
    leg_acc->SetBorderSize(0);
    leg_acc->AddEntry(h_acc_den_comb, "Acceptance denominator", "l");
    leg_acc->AddEntry(h_acc_num_comb, "Acceptance numerator", "l");
    leg_acc->Draw();

    pad2_acc->cd();
    pad2_acc->SetLogx();

    h_acc->SetTitle("");
    h_acc->SetMarkerStyle(20);
    h_acc->SetLineWidth(2);

    h_acc->GetYaxis()->SetTitle("Acceptance");
    h_acc->GetXaxis()->SetTitle("m_{T}^{gen} [GeV]");

    h_acc->GetYaxis()->SetRangeUser(0.0, 1.1);

    h_acc->GetYaxis()->SetTitleSize(0.10);
    h_acc->GetYaxis()->SetTitleOffset(0.45);
    h_acc->GetYaxis()->SetLabelSize(0.09);

    h_acc->GetXaxis()->SetTitleSize(0.11);
    h_acc->GetXaxis()->SetLabelSize(0.10);

    h_acc->Draw("E1");

    c1->SaveAs("acceptance.png");

    TCanvas* c2 = new TCanvas("c2", "", 800, 800);

    TPad* pad1_eff = new TPad("pad1_eff", "", 0.0, 0.30, 1.0, 1.0);
    TPad* pad2_eff = new TPad("pad2_eff", "", 0.0, 0.00, 1.0, 0.30);

    pad1_eff->SetBottomMargin(0.02);
    pad2_eff->SetTopMargin(0.03);
    pad2_eff->SetBottomMargin(0.30);

    pad1_eff->Draw();
    pad2_eff->Draw();

    pad1_eff->cd();
    pad1_eff->SetLogx();
    pad1_eff->SetLogy();

    h_acc_num_comb->SetTitle("");
    h_acc_num_comb->SetLineColor(kBlack);
    h_acc_num_comb->SetLineWidth(2);

    h_eff_num_comb->SetLineColor(kBlue);
    h_eff_num_comb->SetLineWidth(2);

    h_acc_num_comb->GetYaxis()->SetTitle("Events");
    h_acc_num_comb->GetXaxis()->SetLabelSize(0);

    h_acc_num_comb->Draw("HIST E");
    h_eff_num_comb->Draw("HIST E SAME");

    TLegend* leg_eff = new TLegend(0.60, 0.72, 0.88, 0.88);
    leg_eff->SetBorderSize(0);
    leg_eff->AddEntry(h_acc_num_comb, "Efficiency denominator", "l");
    leg_eff->AddEntry(h_eff_num_comb, "Efficiency numerator", "l");
    leg_eff->Draw();

    pad2_eff->cd();
    pad2_eff->SetLogx();

    h_eff->SetTitle("");
    h_eff->SetMarkerStyle(20);
    h_eff->SetLineWidth(2);

    h_eff->GetYaxis()->SetTitle("Efficiency");
    h_eff->GetXaxis()->SetTitle("m_{T}^{gen} [GeV]");

    h_eff->GetYaxis()->SetRangeUser(0.0, 1.1);

    h_eff->GetYaxis()->SetTitleSize(0.10);
    h_eff->GetYaxis()->SetTitleOffset(0.45);
    h_eff->GetYaxis()->SetLabelSize(0.09);

    h_eff->GetXaxis()->SetTitleSize(0.11);
    h_eff->GetXaxis()->SetLabelSize(0.10);

    h_eff->Draw("E1");

    c2->SaveAs("efficiency.png");

    TCanvas* c3 = new TCanvas("c3", "", 800, 700);

    h_acc_eff->SetMarkerStyle(20);
    h_acc_eff->SetLineWidth(2);
    h_acc_eff->GetYaxis()->SetRangeUser(0.0, 1.1);
    h_acc_eff->Draw("E1");
    c3->SetLogx();

    c3->SaveAs("acceptance_efficiency.png");
}