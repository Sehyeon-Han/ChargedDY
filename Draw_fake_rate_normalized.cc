#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>

struct SampleInfo
{
    TString fileName;
    TString sampleName;
    TString processLabel;
    Color_t color;
    double xsec;
    double sumW;
};

double GetScale(const SampleInfo& sample, double lumi)
{
    if(sample.sumW == 0.0)
    {
        std::cerr << "[ERROR] sumW = 0: " << sample.sampleName << std::endl;
        return 0.0;
    }

    return sample.xsec * lumi / sample.sumW;
}

TString SafeName(TString name)
{
    name.ReplaceAll("#rightarrow", "To");
    name.ReplaceAll("#bar", "bar");
    name.ReplaceAll("#mu", "mu");
    name.ReplaceAll("#tau", "tau");
    name.ReplaceAll("#nu", "nu");
    name.ReplaceAll("#", "");
    name.ReplaceAll("{", "");
    name.ReplaceAll("}", "");
    name.ReplaceAll(" ", "_");
    name.ReplaceAll("-", "_");

    return name;
}

TH1D* LoadTH1D(const TString& fileName,
               const TString& histName,
               const TString& cloneName)
{
    TFile* file = TFile::Open(fileName, "READ");

    if(!file || file->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << fileName << std::endl;

        if(file)
            file->Close();

        return nullptr;
    }

    TH1D* source = dynamic_cast<TH1D*>(file->Get(histName));

    if(!source)
    {
        std::cerr << "[ERROR] Cannot find " << histName
                  << " in " << fileName << std::endl;

        file->Close();
        return nullptr;
    }

    TH1D* hist = dynamic_cast<TH1D*>(source->Clone(cloneName));
    hist->SetDirectory(nullptr);

    file->Close();

    return hist;
}

TH2D* LoadTH2D(const TString& fileName,
               const TString& histName,
               const TString& cloneName)
{
    TFile* file = TFile::Open(fileName, "READ");

    if(!file || file->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << fileName << std::endl;

        if(file)
            file->Close();

        return nullptr;
    }

    TH2D* source = dynamic_cast<TH2D*>(file->Get(histName));

    if(!source)
    {
        std::cerr << "[ERROR] Cannot find " << histName
                  << " in " << fileName << std::endl;

        file->Close();
        return nullptr;
    }

    TH2D* hist = dynamic_cast<TH2D*>(source->Clone(cloneName));
    hist->SetDirectory(nullptr);

    file->Close();

    return hist;
}

bool SameBinning1D(const TH1D* h1, const TH1D* h2)
{
    if(!h1 || !h2)
        return false;

    if(h1->GetNbinsX() != h2->GetNbinsX())
        return false;

    for(int i = 1; i <= h1->GetNbinsX(); ++i)
    {
        double edge1 = h1->GetXaxis()->GetBinLowEdge(i);
        double edge2 = h2->GetXaxis()->GetBinLowEdge(i);

        if(std::fabs(edge1 - edge2) > 1e-9)
            return false;
    }

    double last1 = h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX());
    double last2 = h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX());

    return std::fabs(last1 - last2) < 1e-9;
}

bool SameBinning2D(const TH2D* h1, const TH2D* h2)
{
    if(!h1 || !h2)
        return false;

    if(h1->GetNbinsX() != h2->GetNbinsX())
        return false;

    if(h1->GetNbinsY() != h2->GetNbinsY())
        return false;

    for(int ix = 1; ix <= h1->GetNbinsX(); ++ix)
    {
        double edge1 = h1->GetXaxis()->GetBinLowEdge(ix);
        double edge2 = h2->GetXaxis()->GetBinLowEdge(ix);

        if(std::fabs(edge1 - edge2) > 1e-9)
            return false;
    }

    for(int iy = 1; iy <= h1->GetNbinsY(); ++iy)
    {
        double edge1 = h1->GetYaxis()->GetBinLowEdge(iy);
        double edge2 = h2->GetYaxis()->GetBinLowEdge(iy);

        if(std::fabs(edge1 - edge2) > 1e-9)
            return false;
    }

    double xLast1 = h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX());
    double xLast2 = h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX());
    double yLast1 = h1->GetYaxis()->GetBinUpEdge(h1->GetNbinsY());
    double yLast2 = h2->GetYaxis()->GetBinUpEdge(h2->GetNbinsY());

    return std::fabs(xLast1 - xLast2) < 1e-9 &&
           std::fabs(yLast1 - yLast2) < 1e-9;
}

void DrawDataMCStack(const TString& dataFile,
                     const std::vector<SampleInfo>& samples,
                     const TString& histName,
                     const TString& xTitle,
                     const TString& outputName,
                     double lumi,
                     bool useLogX,
                     bool useLogY,
                     TFile* outputFile)
{
    TH1D* hData = LoadTH1D(dataFile,
                           histName,
                           Form("hData_%s", histName.Data()));

    if(!hData)
        return;

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.9);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);
    hData->SetLineWidth(1);

    std::map<std::string, TH1D*> processHists;
    std::map<std::string, Color_t> processColors;

    for(const SampleInfo& sample : samples)
    {
        TH1D* hist = LoadTH1D(sample.fileName,
                              histName,
                              Form("%s_%s",
                                   histName.Data(),
                                   sample.sampleName.Data()));

        if(!hist)
            continue;

        if(!SameBinning1D(hData, hist))
        {
            std::cerr << "[ERROR] Binning mismatch: "
                      << sample.fileName
                      << ", histogram = " << histName << std::endl;

            delete hist;
            continue;
        }

        double rawIntegral = hist->Integral(0, hist->GetNbinsX() + 1);
        double scale = GetScale(sample, lumi);

        if(scale == 0.0)
        {
            delete hist;
            continue;
        }

        hist->Scale(scale);

        double normalizedIntegral =
            hist->Integral(0, hist->GetNbinsX() + 1);

        std::cout << "[MC] " << sample.sampleName
                  << " | raw = " << rawIntegral
                  << " | scale = " << scale
                  << " | normalized = " << normalizedIntegral
                  << std::endl;

        std::string key = sample.processLabel.Data();

        if(processHists.find(key) == processHists.end())
        {
            TString safeName = SafeName(sample.processLabel);

            processHists[key] = dynamic_cast<TH1D*>(
                hist->Clone(Form("h_%s_%s",
                                 safeName.Data(),
                                 histName.Data()))
            );

            processHists[key]->SetDirectory(nullptr);
            processColors[key] = sample.color;
        }
        else
        {
            processHists[key]->Add(hist);
        }

        delete hist;
    }

    std::vector<std::string> stackOrder =
    {
        "W#rightarrow#mu#nu",
        "W#rightarrow#tau#nu",
        "Diboson",
        "Z#rightarrow#mu#mu",
        "single top",
        "t#bar{t}"
    };

    THStack* stack = new THStack(Form("stack_%s", histName.Data()), "");
    TH1D* hMCsum = nullptr;

    for(const std::string& process : stackOrder)
    {
        auto it = processHists.find(process);

        if(it == processHists.end())
            continue;

        TH1D* hist = it->second;

        hist->SetFillColor(processColors[process]);
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(1);

        stack->Add(hist);

        if(!hMCsum)
        {
            hMCsum = dynamic_cast<TH1D*>(
                hist->Clone(Form("hMCsum_%s", histName.Data()))
            );

            hMCsum->SetDirectory(nullptr);
        }
        else
        {
            hMCsum->Add(hist);
        }
    }

    if(!hMCsum)
    {
        std::cerr << "[ERROR] No MC histogram for "
                  << histName << std::endl;

        delete hData;
        delete stack;
        return;
    }

    TCanvas* canvas = new TCanvas(Form("c_%s", histName.Data()),
                                  histName,
                                  900,
                                  700);

    canvas->SetTopMargin(0.08);
    canvas->SetBottomMargin(0.13);
    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.05);

    if(useLogX)
        canvas->SetLogx();

    if(useLogY)
        canvas->SetLogy();

    canvas->cd();

    stack->Draw("HIST");

    stack->GetYaxis()->SetTitle("Events");
    stack->GetYaxis()->SetTitleSize(0.055);
    stack->GetYaxis()->SetLabelSize(0.045);
    stack->GetYaxis()->SetTitleOffset(1.1);

    stack->GetXaxis()->SetTitle(xTitle);
    stack->GetXaxis()->SetTitleSize(0.055);
    stack->GetXaxis()->SetLabelSize(0.045);
    stack->GetXaxis()->SetTitleOffset(1.05);

    double yMax = std::max(hData->GetMaximum(),
                           hMCsum->GetMaximum());

    if(useLogY)
    {
        stack->SetMinimum(1e-3);
        stack->SetMaximum(yMax > 0.0 ? 20.0 * yMax : 1.0);
    }
    else
    {
        stack->SetMinimum(0.0);
        stack->SetMaximum(yMax > 0.0 ? 1.45 * yMax : 1.0);
    }

    hData->Draw("E1 SAME");

    TLegend* legend = new TLegend(0.62, 0.55, 0.89, 0.88);

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.032);

    legend->AddEntry(hData, "Data", "lep");

    std::vector<std::string> legendOrder =
    {
        "t#bar{t}",
        "single top",
        "Z#rightarrow#mu#mu",
        "Diboson",
        "W#rightarrow#tau#nu",
        "W#rightarrow#mu#nu"
    };

    for(const std::string& process : legendOrder)
    {
        auto it = processHists.find(process);

        if(it != processHists.end())
            legend->AddEntry(it->second, process.c_str(), "f");
    }

    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.045);
    /*latex.DrawLatex(0.15, 0.93,
                    "#bf{CMS} #it{Work in progress}");*/

    latex.SetTextAlign(31);
    latex.SetTextSize(0.035);
    /*latex.DrawLatex(0.94, 0.93,
                    Form("%.2f fb^{-1} (13 TeV)",
                         lumi / 1000.0));*/

    canvas->SaveAs(outputName);

    if(outputFile)
    {
        outputFile->cd();

        canvas->Write();
        hData->Write();
        hMCsum->Write();

        for(const auto& item : processHists)
            item.second->Write();
    }

    std::cout << "========================================" << std::endl;
    std::cout << "Histogram     : " << histName << std::endl;
    std::cout << "Data integral : "
              << hData->Integral(0, hData->GetNbinsX() + 1)
              << std::endl;
    std::cout << "MC integral   : "
              << hMCsum->Integral(0, hMCsum->GetNbinsX() + 1)
              << std::endl;
    std::cout << "Output        : " << outputName << std::endl;
    std::cout << "========================================" << std::endl;
}

void MakeFakeRate(const TString& dataFile,
                  const std::vector<SampleInfo>& samples,
                  double lumi,
                  TFile* outputFile)
{
    TH2D* hDataDen = LoadTH2D(dataFile,
                              "h_fake_den",
                              "h_data_den");

    TH2D* hDataNum = LoadTH2D(dataFile,
                              "h_fake_num",
                              "h_data_num");

    if(!hDataDen || !hDataNum)
    {
        std::cerr << "[ERROR] Cannot load data fake-rate histograms."
                  << std::endl;

        delete hDataDen;
        delete hDataNum;
        return;
    }

    TH2D* hPromptDen = dynamic_cast<TH2D*>(
        hDataDen->Clone("h_prompt_den")
    );

    TH2D* hPromptNum = dynamic_cast<TH2D*>(
        hDataNum->Clone("h_prompt_num")
    );

    hPromptDen->Reset("ICES");
    hPromptNum->Reset("ICES");

    hPromptDen->SetDirectory(nullptr);
    hPromptNum->SetDirectory(nullptr);

    for(const SampleInfo& sample : samples)
    {
        TH2D* hDen = LoadTH2D(sample.fileName,
                              "h_fake_den",
                              Form("h_den_%s",
                                   sample.sampleName.Data()));

        TH2D* hNum = LoadTH2D(sample.fileName,
                              "h_fake_num",
                              Form("h_num_%s",
                                   sample.sampleName.Data()));

        if(!hDen || !hNum)
        {
            delete hDen;
            delete hNum;
            continue;
        }

        if(!SameBinning2D(hDataDen, hDen) ||
           !SameBinning2D(hDataNum, hNum))
        {
            std::cerr << "[ERROR] 2D binning mismatch: "
                      << sample.fileName << std::endl;

            delete hDen;
            delete hNum;
            continue;
        }

        double scale = GetScale(sample, lumi);

        if(scale == 0.0)
        {
            delete hDen;
            delete hNum;
            continue;
        }

        double rawDen = hDen->Integral(0, hDen->GetNbinsX() + 1,
                                      0, hDen->GetNbinsY() + 1);

        double rawNum = hNum->Integral(0, hNum->GetNbinsX() + 1,
                                      0, hNum->GetNbinsY() + 1);

        hDen->Scale(scale);
        hNum->Scale(scale);

        hPromptDen->Add(hDen);
        hPromptNum->Add(hNum);

        std::cout << "[FAKE RATE MC] " << sample.sampleName
                  << " | raw den = " << rawDen
                  << " | raw num = " << rawNum
                  << " | scale = " << scale
                  << " | norm den = "
                  << hDen->Integral(0, hDen->GetNbinsX() + 1,
                                    0, hDen->GetNbinsY() + 1)
                  << " | norm num = "
                  << hNum->Integral(0, hNum->GetNbinsX() + 1,
                                    0, hNum->GetNbinsY() + 1)
                  << std::endl;

        if(outputFile)
        {
            outputFile->cd();
            hDen->Write();
            hNum->Write();
        }

        delete hDen;
        delete hNum;
    }

    TH2D* hFakeDen = dynamic_cast<TH2D*>(
        hDataDen->Clone("h_fake_den_subtracted")
    );

    TH2D* hFakeNum = dynamic_cast<TH2D*>(
        hDataNum->Clone("h_fake_num_subtracted")
    );

    hFakeDen->SetDirectory(nullptr);
    hFakeNum->SetDirectory(nullptr);

    hFakeDen->Add(hPromptDen, -1.0);
    hFakeNum->Add(hPromptNum, -1.0);

    TH2D* hFakeRate = dynamic_cast<TH2D*>(
        hFakeNum->Clone("h_fake_rate")
    );

    hFakeRate->Reset("ICES");
    hFakeRate->SetDirectory(nullptr);
    hFakeRate->SetTitle("Muon fake rate;|#eta|;p_{T} [GeV]");

    for(int ix = 1; ix <= hFakeRate->GetNbinsX(); ++ix)
    {
        for(int iy = 1; iy <= hFakeRate->GetNbinsY(); ++iy)
        {
            double numerator = hFakeNum->GetBinContent(ix, iy);
            double denominator = hFakeDen->GetBinContent(ix, iy);
            double numeratorError = hFakeNum->GetBinError(ix, iy);
            double denominatorError = hFakeDen->GetBinError(ix, iy);

            if(denominator <= 0.0)
            {
                hFakeRate->SetBinContent(ix, iy, 0.0);
                hFakeRate->SetBinError(ix, iy, 0.0);

                std::cerr << "[WARNING] denominator <= 0"
                          << " | eta bin = " << ix
                          << " | pt bin = " << iy
                          << " | denominator = " << denominator
                          << std::endl;

                continue;
            }

            double rate = numerator / denominator;

            hFakeRate->SetBinContent(ix, iy, rate);

            std::cout << "[BIN] eta = " << ix
                      << ", pt = " << iy
                      << " | data den = "
                      << hDataDen->GetBinContent(ix, iy)
                      << " | prompt den = "
                      << hPromptDen->GetBinContent(ix, iy)
                      << " | fake den = " << denominator
                      << " | data num = "
                      << hDataNum->GetBinContent(ix, iy)
                      << " | prompt num = "
                      << hPromptNum->GetBinContent(ix, iy)
                      << " | fake num = " << numerator
                      << " | rate = " << rate
                      << std::endl;
        }
    }

    TCanvas* cRate = new TCanvas("c_fake_rate",
                                  "Fake rate",
                                  950,
                                  750);

    cRate->SetLeftMargin(0.12);
    cRate->SetRightMargin(0.18);
    cRate->SetBottomMargin(0.12);

    hFakeRate->SetMinimum(0.0);
    hFakeRate->SetMaximum(1.0);
    hFakeRate->SetMarkerSize(1.3);
    hFakeRate->Draw("COLZ TEXT");

    cRate->SaveAs("fake_mu_rate_subtracted.png");

    TCanvas* cDen = new TCanvas("c_fake_den_comparison",
                                 "Fake denominator comparison",
                                 1500,
                                 450);

    cDen->Divide(3, 1);

    cDen->cd(1);
    gPad->SetRightMargin(0.18);
    hDataDen->SetTitle("Data denominator;|#eta|;p_{T} [GeV]");
    hDataDen->Draw("COLZ TEXT");

    cDen->cd(2);
    gPad->SetRightMargin(0.18);
    hPromptDen->SetTitle("Prompt MC denominator;|#eta|;p_{T} [GeV]");
    hPromptDen->Draw("COLZ TEXT");

    cDen->cd(3);
    gPad->SetRightMargin(0.18);
    hFakeDen->SetTitle("Data - prompt MC denominator;|#eta|;p_{T} [GeV]");
    hFakeDen->Draw("COLZ TEXT");

    cDen->SaveAs("fake_mu_den_comparison.png");

    TCanvas* cNum = new TCanvas("c_fake_num_comparison",
                                 "Fake numerator comparison",
                                 1500,
                                 450);

    cNum->Divide(3, 1);

    cNum->cd(1);
    gPad->SetRightMargin(0.18);
    hDataNum->SetTitle("Data numerator;|#eta|;p_{T} [GeV]");
    hDataNum->Draw("COLZ TEXT");

    cNum->cd(2);
    gPad->SetRightMargin(0.18);
    hPromptNum->SetTitle("Prompt MC numerator;|#eta|;p_{T} [GeV]");
    hPromptNum->Draw("COLZ TEXT");

    cNum->cd(3);
    gPad->SetRightMargin(0.18);
    hFakeNum->SetTitle("Data - prompt MC numerator;|#eta|;p_{T} [GeV]");
    hFakeNum->Draw("COLZ TEXT");

    cNum->SaveAs("fake_mu_num_comparison.png");

    if(outputFile)
    {
        outputFile->cd();

        hDataDen->Write();
        hDataNum->Write();
        hPromptDen->Write();
        hPromptNum->Write();
        hFakeDen->Write();
        hFakeNum->Write();
        hFakeRate->Write();

        cRate->Write();
        cDen->Write();
        cNum->Write();
    }

    std::cout << "========================================" << std::endl;
    std::cout << "Data denominator   = "
              << hDataDen->Integral(0, hDataDen->GetNbinsX() + 1,
                                    0, hDataDen->GetNbinsY() + 1)
              << std::endl;
    std::cout << "Prompt denominator = "
              << hPromptDen->Integral(0, hPromptDen->GetNbinsX() + 1,
                                      0, hPromptDen->GetNbinsY() + 1)
              << std::endl;
    std::cout << "Fake denominator   = "
              << hFakeDen->Integral(0, hFakeDen->GetNbinsX() + 1,
                                    0, hFakeDen->GetNbinsY() + 1)
              << std::endl;
    std::cout << "Data numerator     = "
              << hDataNum->Integral(0, hDataNum->GetNbinsX() + 1,
                                    0, hDataNum->GetNbinsY() + 1)
              << std::endl;
    std::cout << "Prompt numerator   = "
              << hPromptNum->Integral(0, hPromptNum->GetNbinsX() + 1,
                                      0, hPromptNum->GetNbinsY() + 1)
              << std::endl;
    std::cout << "Fake numerator     = "
              << hFakeNum->Integral(0, hFakeNum->GetNbinsX() + 1,
                                    0, hFakeNum->GetNbinsY() + 1)
              << std::endl;
    std::cout << "========================================" << std::endl;
}

void Draw_fake_rate_normalized(const char* dataFile = "fake_mu_merged.root",
                    const char* outputRootFile = "fake_mu_normalized_results.root",
                    double lumi = 19520.0)
{
    TH1::SetDefaultSumw2();

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".3f");

    std::vector<SampleInfo> samples =
    {
        {"fake_mu_WToMuNu_M-200_merged.root",
         "WToMuNu_M200",
         "W#rightarrow#mu#nu",
         kGreen + 1,
         8.78288653e00,
         3250000.0},

        {"fake_mu_WToMuNu_M-500_merged.root",
         "WToMuNu_M500",
         "W#rightarrow#mu#nu",
         kGreen + 1,
         2.76542446e-01,
         750000.0},

        {"fake_mu_WToMuNu_M-1000_merged.root",
         "WToMuNu_M1000",
         "W#rightarrow#mu#nu",
         kGreen + 1,
         1.58702653e-02,
         250000.0},

        {"fake_mu_WToMuNu_M-2000_merged.root",
         "WToMuNu_M2000",
         "W#rightarrow#mu#nu",
         kGreen + 1,
         4.15066636e-04,
         50000.0},

        {"fake_mu_WToTauNu_M-200_merged.root",
         "WToTauNu_M200",
         "W#rightarrow#tau#nu",
         kGreen + 2,
         8.78288653e00,
         1250000.0},

        {"fake_mu_WToTauNu_M-500_merged.root",
         "WToTauNu_M500",
         "W#rightarrow#tau#nu",
         kGreen + 2,
         2.76542446e-01,
         750000.0},

        {"fake_mu_WToTauNu_M-1000_merged.root",
         "WToTauNu_M1000",
         "W#rightarrow#tau#nu",
         kGreen + 2,
         1.58702653e-02,
         250000.0},

        {"fake_mu_WToTauNu_M-2000_merged.root",
         "WToTauNu_M2000",
         "W#rightarrow#tau#nu",
         kGreen + 2,
         4.15066636e-04,
         50000.0},

        {"fake_mu_DYMC_200_merged.root",
         "DY200",
         "Z#rightarrow#mu#mu",
         kYellow,
         2.78,
         2673760.516799},

        {"fake_mu_DYMC_400_merged.root",
         "DY400",
         "Z#rightarrow#mu#mu",
         kYellow,
         0.15,
         44344.277472},

        {"fake_mu_DYMC_500_merged.root",
         "DY500",
         "Z#rightarrow#mu#mu",
         kYellow,
         0.084,
         24718.574857},

        {"fake_mu_DYMC_700_merged.root",
         "DY700",
         "Z#rightarrow#mu#mu",
         kYellow,
         0.013,
         3951.644255},

        {"fake_mu_DYMC_800_merged.root",
         "DY800",
         "Z#rightarrow#mu#mu",
         kYellow,
         0.011,
         3345.560043},

        {"fake_mu_DYMC_1000_merged.root",
         "DY1000",
         "Z#rightarrow#mu#mu",
         kYellow,
         0.006,
         1782.846123},

        {"fake_mu_DYMC_1500_merged.root",
         "DY1500",
         "Z#rightarrow#mu#mu",
         kYellow,
         0.00081,
         244.385548},

        {"fake_mu_DYMC_2000_merged.root",
         "DY2000",
         "Z#rightarrow#mu#mu",
         kYellow,
         0.0002,
         62.650483},

        {"fake_mu_TTTo2L2Nu_merged.root",
         "TTTo2L2Nu",
         "t#bar{t}",
         kRed,
         89.31,
         2.70453e09},

        {"fake_mu_TTToSemiLeptonic_merged.root",
         "TTToSemiLeptonic",
         "t#bar{t}",
         kRed,
         367.8,
         3.97723e10},

        {"fake_mu_ST_s-channel_4f_leptonDecays_merged.root",
         "ST_s_channel",
         "single top",
         kBlue,
         10.32,
         1.95962e07},

        {"fake_mu_ST_t-channel_anitop_4f_InclusiveDecays_merged.root",
         "ST_t_channel_antitop",
         "single top",
         kBlue,
         80.0,
         1.98386e09},

        {"fake_mu_ST_t-channel_top_4f_InclusiveDecays_merged.root",
         "ST_t_channel_top",
         "single top",
         kBlue,
         134.2,
         5.94814e09},

        {"fake_mu_ST_tW_antitop_5f_inclusiveDecays_merged.root",
         "ST_tW_antitop",
         "single top",
         kBlue,
         39.65,
         7.47663e07},

        {"fake_mu_ST_tW_top_5f_inclusiveDecays_merged.root",
         "ST_tW_top",
         "single top",
         kBlue,
         39.65,
         7.46247e07},

        {"fake_mu_WWTo1L1Nu2Q_merged.root",
         "WWTo1L1Nu2Q",
         "Diboson",
         kOrange,
         51.65,
         1.69101e09},

        {"fake_mu_WWTo2L2Nu_merged.root",
         "WWTo2L2Nu",
         "Diboson",
         kOrange,
         11.09,
         3.34565e07},

        {"fake_mu_WWTo4Q_4f_merged.root",
         "WWTo4Q",
         "Diboson",
         kOrange,
         51.03,
         1.67235e09},

        {"fake_mu_WZTo1L1Nu2Q_4f_merged.root",
         "WZTo1L1Nu2Q",
         "Diboson",
         kOrange,
         9.119,
         5.58378e07},

        {"fake_mu_WZTo1L3Nu_4f_merged.root",
         "WZTo1L3Nu",
         "Diboson",
         kOrange,
         3.414,
         6.96859e06},

        {"fake_mu_WZTo2Q2Nu_4f_merged.root",
         "WZTo2Q2Nu",
         "Diboson",
         kOrange,
         6.331,
         4.84591e07},

        {"fake_mu_WZTo3LNu_merged.root",
         "WZTo3LNu",
         "Diboson",
         kOrange,
         5.213,
         8.15165e07},

        {"fake_mu_ZZ_merged.root",
         "ZZ",
         "Diboson",
         kOrange,
         12.17,
         1282000.0}
    };

    TFile* outputFile = TFile::Open(outputRootFile, "RECREATE");

    if(!outputFile || outputFile->IsZombie())
    {
        std::cerr << "[ERROR] Cannot create output file: "
                  << outputRootFile << std::endl;

        return;
    }

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_pt_den",
                    "Muon p_{T} [GeV]",
                    "fake_CR_stack_Muon_pt_den.png",
                    lumi,
                    true,
                    true,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_pt_num",
                    "Muon p_{T} [GeV]",
                    "fake_CR_stack_Muon_pt_num.png",
                    lumi,
                    true,
                    true,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_eta_den",
                    "Muon #eta",
                    "fake_CR_stack_Muon_eta_den.png",
                    lumi,
                    false,
                    false,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_eta_num",
                    "Muon #eta",
                    "fake_CR_stack_Muon_eta_num.png",
                    lumi,
                    false,
                    false,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_phi_den",
                    "Muon #phi",
                    "fake_CR_stack_Muon_phi_den.png",
                    lumi,
                    false,
                    false,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_phi_num",
                    "Muon #phi",
                    "fake_CR_stack_Muon_phi_num.png",
                    lumi,
                    false,
                    false,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_MET_pt_den",
                    "MET p_{T} [GeV]",
                    "fake_CR_stack_MET_pt_den.png",
                    lumi,
                    false,
                    true,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_MET_pt_num",
                    "MET p_{T} [GeV]",
                    "fake_CR_stack_MET_pt_num.png",
                    lumi,
                    false,
                    true,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_MET_phi_den",
                    "MET #phi",
                    "fake_CR_stack_MET_phi_den.png",
                    lumi,
                    false,
                    false,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_MET_phi_num",
                    "MET #phi",
                    "fake_CR_stack_MET_phi_num.png",
                    lumi,
                    false,
                    false,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_tkRelIso_den",
                    "Muon tkRelIso",
                    "fake_CR_stack_Muon_tkRelIso_den.png",
                    lumi,
                    false,
                    true,
                    outputFile);

    DrawDataMCStack(dataFile,
                    samples,
                    "h_Muon_tkRelIso_num",
                    "Muon tkRelIso",
                    "fake_CR_stack_Muon_tkRelIso_num.png",
                    lumi,
                    false,
                    true,
                    outputFile);

    MakeFakeRate(dataFile,
                 samples,
                 lumi,
                 outputFile);

    outputFile->Close();

    std::cout << "All plots created." << std::endl;
    std::cout << "Output ROOT file: "
              << outputRootFile << std::endl;
}