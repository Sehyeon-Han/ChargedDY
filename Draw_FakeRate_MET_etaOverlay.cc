#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <TString.h>

using namespace std;

void Draw_FakeRate_MET_etaOverlay()
{
    TFile *fData = TFile::Open("fake_mu_merged_dependence.root");

    struct Sample
    {
        const char* file;
        double xsec;
        double sumW;
    };

    const double lumi = 19520.0;

    vector<Sample> mcFiles =
    {
        {"fake_mu_merged_ST_s_channel_dependence.root", 10.32, 1.95962e+07},
        {"fake_mu_merged_ST_tW_antitop_dependence.root", 39.65, 7.47663e+07},
        {"fake_mu_merged_ST_tW_top_dependence.root", 39.65, 7.46247e+07},
        {"fake_mu_merged_ST_t_antitop_dependence.root", 80.0, 1.98386e+09},
        {"fake_mu_merged_ST_t_top_dependence.root", 134.2, 5.94814e+09},
        {"fake_mu_merged_TTTo2L2Nu_dependence.root", 89.31, 2.70453e+09},
        {"fake_mu_merged_TTToSemiLeptonic_dependence.root", 367.8, 3.97723e+10},
        {"fake_mu_merged_WToMuNu_M200_dependence.root", 8.78288653, 3250000},
        {"fake_mu_merged_WToMuNu_M500_dependence.root", 2.76542446e-1, 750000},
        {"fake_mu_merged_WToMuNu_M1000_dependence.root", 1.58702653e-2, 250000},
        {"fake_mu_merged_WToMuNu_M2000_dependence.root", 4.15066636e-4, 50000},
        {"fake_mu_merged_WToTauNu_M200_dependence.root", 8.78288653, 1250000},
        {"fake_mu_merged_WToTauNu_M500_dependence.root", 2.76542446e-1, 750000},
        {"fake_mu_merged_WToTauNu_M1000_dependence.root", 1.58702653e-2, 250000},
        {"fake_mu_merged_WToTauNu_M2000_dependence.root", 4.15066636e-4, 50000},
        {"fake_mu_merged_WWTo1L1Nu2Q_dependence.root", 51.65, 1.69101e+09},
        {"fake_mu_merged_WWTo2L2Nu_dependence.root", 11.09, 3.34565e+07},
        {"fake_mu_merged_WWTo4Q_dependence.root", 51.03, 1.67235e+09},
        {"fake_mu_merged_WZTo1L1Nu2Q_4f_dependence.root", 9.119, 5.58378e+07},
        {"fake_mu_merged_WZTo1L3Nu_4f_dependence.root", 3.414, 6.96859e+06},
        {"fake_mu_merged_WZTo2Q2Nu_4f_dependence.root", 6.331, 4.84591e+07},
        {"fake_mu_merged_WZTo3LNu_dependence.root", 5.213, 8.15165e+07},
        {"fake_mu_merged_ZZ_dependence.root", 12.17, 1282000},
        {"fake_mu_DYMC_200_merged_dependence.root", 2.78, 2673760.516799},
        {"fake_mu_DYMC_400_merged_dependence.root", 0.15, 44344.277472},
        {"fake_mu_DYMC_500_merged_dependence.root", 0.084, 24718.574857},
        {"fake_mu_DYMC_700_merged_dependence.root", 0.013, 3951.644255},
        {"fake_mu_DYMC_800_merged_dependence.root", 0.011, 3345.560043},
        {"fake_mu_DYMC_1000_merged_dependence.root", 0.006, 1782.846123},
        {"fake_mu_DYMC_1500_merged_dependence.root", 0.00081, 244.385548},
        {"fake_mu_DYMC_2000_merged_dependence.root", 0.0002, 62.650483}
    };

    int colors[4] = {kBlack, kRed, kBlue, kGreen + 2};
    int markers[4] = {20, 21, 22, 23};

    const char* etaLabel[4] = 
    {
        "0.0 < |#eta| < 0.9",
        "0.9 < |#eta| < 1.2",
        "1.2 < |#eta| < 2.1",
        "2.1 < |#eta| < 2.4" 
    };

    const char* ptLabel[5] = 
    {
        "65 < p_{T} < 100 GeV",
        "100 < p_{T} < 150 GeV",
        "150 < p_{T} < 200 GeV",
        "200 < p_{T} < 300 GeV",
        "300 < p_{T} < 1500 GeV"
    };

    for(int ipt = 1; ipt <= 5; ipt++)
    {
        TCanvas *c = new TCanvas(Form("c_pt%d", ipt), "", 800, 700);
        TLegend *leg = new TLegend(0.58, 0.65, 0.88, 0.88);
        leg->SetBorderSize(0);

        bool first = true;

        for(int ieta = 1; ieta <= 4; ieta++)
        {
            TString denName = Form("h_fake_den_MET_pt%d_eta%d", ipt, ieta);
            TString numName = Form("h_fake_num_MET_pt%d_eta%d", ipt, ieta);

            TH1D *hdataDen = (TH1D*)fData->Get(denName);
            TH1D *hdataNum = (TH1D*)fData->Get(numName);

            if(!hdataDen || !hdataNum)
            {
                std::cout << "Missing data hist: " << denName << " or " << numName << std::endl;
                continue;
            }

            TH1D *hMCDen = (TH1D*)hdataDen->Clone(Form("hMCDen_pt%d_eta%d", ipt, ieta));
            TH1D *hMCNum = (TH1D*)hdataNum->Clone(Form("hMCNum_pt%d_eta%d", ipt, ieta));

            hMCDen->Reset();
            hMCNum->Reset();

            for(auto sample : mcFiles)
            {
                TFile *fMC = TFile::Open(sample.file);

                TH1D *hDen = (TH1D*)fMC->Get(denName);
                TH1D *hNum = (TH1D*)fMC->Get(numName);

                double scale = lumi * sample.xsec / sample.sumW;

                if(hDen) hMCDen->Add(hDen, scale);
                if(hNum) hMCNum->Add(hNum, scale);

                fMC->Close();
            }

            TH1D *hFakeDen = (TH1D*)hdataDen->Clone(Form("hFakeDen_pt%d_eta%d", ipt, ieta));
            TH1D *hFakeNum = (TH1D*)hdataNum->Clone(Form("hFakeNum_pt%d_eta%d", ipt, ieta));

            hFakeDen->Add(hMCDen, -1.0);
            hFakeNum->Add(hMCNum, -1.0);

            TGraphErrors *gr = new TGraphErrors();
            int p = 0;

            for(int ibin = 1; ibin <= hFakeDen->GetNbinsX(); ibin++)
            {
                double den = hFakeDen->GetBinContent(ibin);
                double num = hFakeNum->GetBinContent(ibin);

                cout << "pt = " << ipt << ", eta = " << ieta << ", MET bin = " << ibin << ", den = " << den << ", num = " << num << endl;

                double denErr = hFakeDen->GetBinError(ibin);
                double numErr = hFakeNum->GetBinError(ibin);

                if(den <= 0) continue;
                if(num < 0) continue;

                double fr = num / den;
                double frErr = 0.0;

                if(num > 0)
                {
                    frErr = sqrt(pow(numErr / den, 2) + pow(num * denErr / (den * den), 2));
                }

                double x = hFakeDen->GetBinCenter(ibin);

                gr->SetPoint(p, x, fr);
                gr->SetPointError(p, 0, frErr);

                p++;
            }

            gr->SetMarkerStyle(markers[ieta - 1]);
            gr->SetMarkerColor(colors[ieta - 1]);
            gr->SetLineColor(colors[ieta - 1]);
            gr->SetLineWidth(2);

            if(first)
            {
                gr->SetTitle(Form("Fake rate vs MET (%s);MET [GeV];Fake rate", ptLabel[ipt-1]));
                gr->SetMinimum(0.0);
                gr->SetMaximum(1.0);
                gr->Draw("APL");

                first = false;
            }
            else
            {
                gr->Draw("PL SAME");
            }
            
            leg->AddEntry(gr, etaLabel[ieta - 1], "pl");
        }

        leg->Draw();
        c->SaveAs(Form("FakeRate_vs_MET_pt%d.png", ipt));
    }
    
    fData->Close();
}