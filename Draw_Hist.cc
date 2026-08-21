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
#include <cmath>

struct SampleInfo{
    TString fileName;
    TString label;
    Color_t color;
    double scale;
    bool isData;
};

TH1D* LoadHist(const SampleInfo& s, const char* histName, int idx)
{
    TFile *f = TFile::Open(s.fileName);
    if(!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open file: " << s.fileName << std::endl;
        return nullptr;
    }

    TH1D* h0 = (TH1D*)f->Get(histName);
    if(!h0)
    {
        std::cerr << "[ERROR] Cannot find histogram " << histName << " in file " << s.fileName << std::endl;
        f->Close();
        return nullptr;
    }

    TH1D* h = (TH1D*)h0->Clone(Form("%s_clone_%d", histName, idx));
    h->SetDirectory(0);

    if(TString(histName) == "h_mT")
    {
        int nBins = h->GetNbinsX();

        std::cout << s.fileName
                  << " | last bin = " << h->GetBinContent(nBins)
                  << " | overflow = " << h->GetBinContent(nBins + 1)
                  << std::endl;

        double content = h->GetBinContent(nBins) + h->GetBinContent(nBins + 1);
        double error = std::sqrt(std::pow(h->GetBinError(nBins), 2) + std::pow(h->GetBinError(nBins + 1), 2));

        h->SetBinContent(nBins, content);
        h->SetBinError(nBins, error);

        h->SetBinContent(nBins + 1, 0.0);
        h->SetBinError(nBins + 1, 0.0);
    }

    if(s.label == "QCD")
    {
        std::cout << "\n[QCD] " << histName << std::endl;

        for(int bin = 1; bin <= h->GetNbinsX(); ++bin)
        {
            std::cout << "bin " << bin
                      << " [" << h->GetXaxis()->GetBinLowEdge(bin)
                      << ", " << h->GetXaxis()->GetBinUpEdge(bin)
                      << "] = " << h->GetBinContent(bin)
                      << " +/- " << h->GetBinError(bin);

            if(h->GetBinContent(bin) <= 0.0)
                std::cout << "  <-- NON-POSITIVE";

            std::cout << std::endl;
        }
    }

    if(!s.isData)
    {
        h->Scale(s.scale);
        h->SetFillColor(s.color);
        h->SetLineColor(s.color);
    }
    else
    {
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.0);
    }

    f->Close();
    return h;
}

void DrawOneVariable(const std::vector<SampleInfo>& samples,
                     const char* histName,
                     const char* canvasTitle,
                     const char* xTitle,
                     const char* outName,
                     bool useLogx,
                     bool useLogy,
                     bool drawMassMarks)
{
    std::vector<TH1D*> hists;
    hists.reserve(samples.size());

    TH1D* hData = nullptr;
    TH1D* hMCsum = nullptr;

    THStack* stack = new THStack(Form("stack_%s", histName), Form("%s;%s;Events", canvasTitle, xTitle));

    TH1D* hLegend_W = nullptr;
    TH1D* hLegend_DY = nullptr;
    TH1D* hLegend_TT = nullptr;
    TH1D* hLegend_Dib = nullptr;
    TH1D* hLegend_QCD = nullptr;
    TH1D* hLegend_ST = nullptr;
    TH1D* hLegend_Tau = nullptr;

    for(size_t i = 0; i < samples.size(); ++i)
    {
        TH1D* h = LoadHist(samples[i], histName, i);
        if(!h) continue;

        hists.push_back(h);

        if(samples[i].isData)
        {
            hData = h;
            continue;
        }

        stack->Add(h);

        if(!hMCsum)
        {
            hMCsum = (TH1D*)h->Clone(Form("hMCsum_%s", histName));
            hMCsum->SetDirectory(0);
        }
        else
        {
            hMCsum->Add(h);
        }

        if(samples[i].label == "W#rightarrow#mu#nu" && !hLegend_W) hLegend_W = h;
        if(samples[i].label == "Z#rightarrow#mu#mu" && !hLegend_DY) hLegend_DY = h;
        if(samples[i].label == "t#bar{t}" && !hLegend_TT) hLegend_TT = h;
        if(samples[i].label == "Diboson" && !hLegend_Dib) hLegend_Dib = h;
        if(samples[i].label == "single top" && !hLegend_ST) hLegend_ST = h;
        if(samples[i].label == "W#rightarrow#tau#nu" && !hLegend_Tau) hLegend_Tau = h;
        if(samples[i].label == "QCD" && !hLegend_QCD) hLegend_QCD = h;
    }

    if(!hData || !hMCsum)
    {
        std::cerr << "[ERROR] Missing data or MC sum for " << histName << std::endl;
        return;
    }

    double ymax = std::max(hData->GetMaximum(), hMCsum->GetMaximum());

    TCanvas* c = new TCanvas(Form("c_%s", histName), canvasTitle, 900, 800);

    TPad* pad1 = new TPad(Form("pad1_%s", histName), "pad1", 0, 0.18, 1, 1);
    pad1->SetBottomMargin(0.02);
    if(useLogx) pad1->SetLogx();
    if(useLogy) pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    stack->Draw("HIST");

    if(useLogy)
    {
        stack->SetMinimum(1e-3);
        stack->SetMaximum(10.0 * ymax);
    }
    else
    {
        stack->SetMinimum(0.0);
        stack->SetMaximum(1.4 * ymax);
    }

    stack->GetXaxis()->SetLabelSize(0);
    stack->GetXaxis()->SetTitleSize(0);
    stack->GetXaxis()->SetTickLength(0);

    hData->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.65, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(hData, "Data", "lep");
    if(hLegend_W) leg->AddEntry(hLegend_W, "W#rightarrow#mu#nu", "f");
    if(hLegend_DY) leg->AddEntry(hLegend_DY, "Z#rightarrow#mu#mu", "f");
    if(hLegend_TT) leg->AddEntry(hLegend_TT, "t#bar{t}", "f");
    if(hLegend_Dib) leg->AddEntry(hLegend_Dib, "Diboson", "f");
    if(hLegend_ST) leg->AddEntry(hLegend_ST, "single top", "f");
    if(hLegend_Tau) leg->AddEntry(hLegend_Tau, "W#rightarrow#tau#nu", "f");
    if(hLegend_QCD) leg->AddEntry(hLegend_QCD, "QCD", "f");

    leg->Draw();

    pad1->Update();

    c->cd();

    TPad* pad2 = new TPad(Form("pad2_%s", histName), "pad2", 0, 0, 1, 0.18);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.38);
    if(useLogx) pad2->SetLogx();
    pad2->Draw();
    pad2->cd();

    TH1D* hRatio = (TH1D*)hData->Clone(Form("hRatio_%s", histName));
    hRatio->SetDirectory(0);
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

    hRatio->GetXaxis()->SetTitle(xTitle);
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);

    hRatio->SetMinimum(0.75);
    hRatio->SetMaximum(1.25);
    hRatio->Draw("E1");

    double y_min = hRatio->GetMinimum();
    double y_max = hRatio->GetMaximum();

    /*if(drawMassMarks)
    {
        std::vector<double> marks = {500, 1000, 2000};

        TLatex latex;
        latex.SetTextSize(0.11);
        latex.SetTextFont(42);
        latex.SetTextAlign(23);

        double y_text = -0.35;

        for(double x : marks)
        {
            TLine* l = new TLine(x, y_min, x, y_max);
            l->SetLineStyle(3);
            l->SetLineColor(kGray);
            l->Draw("SAME");

            if(x == 500)
            latex.DrawLatex(x, y_text, "500");
            else if(x == 1000)
            latex.DrawLatex(x, y_text, "10^{3}");
            else if(x == 2000)
            latex.DrawLatex(x, y_text, "2#times10^{3}");
        }
    }*/

    TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);

    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    pad2->Update();
    c->Update();

    c->SaveAs(outName);
}

bool SameMTBinning(const TH1D* h1, const TH1D* h2)
{
    if(!h1 || !h2)
        return false;

    if(h1->GetNbinsX() != h2->GetNbinsX())
        return false;

    for(int i = 1; i <= h1->GetNbinsX(); ++i)
    {
        double edge1 = h1->GetXaxis()->GetBinLowEdge(i);
        double edge2 = h2->GetXaxis()->GetBinLowEdge(i);

        if(std::abs(edge1 - edge2) > 1e-9)
            return false;
    }

    double lastEdge1 =
        h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX());

    double lastEdge2 =
        h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX());

    return std::abs(lastEdge1 - lastEdge2) < 1e-9;
}


void PrintMTBinning(const TH1D* h, const TString& fileName)
{
    if(!h)
        return;

    std::cout << "File  : " << fileName << std::endl;
    std::cout << "Nbins : " << h->GetNbinsX() << std::endl;
    std::cout << "Edges : ";

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        std::cout << h->GetXaxis()->GetBinLowEdge(i)
                  << ", ";
    }

    std::cout
        << h->GetXaxis()->GetBinUpEdge(h->GetNbinsX())
        << std::endl;
}

void CheckMTBinning(const std::vector<SampleInfo>& samples)
{
    TH1D* hReference = nullptr;
    TString referenceFile = "";

    std::cout << "\n====================================" << std::endl;
    std::cout << "Checking h_mT binning" << std::endl;
    std::cout << "====================================" << std::endl;

    for(size_t i = 0; i < samples.size(); ++i)
    {
        TFile* f = TFile::Open(samples[i].fileName);

        if(!f || f->IsZombie())
        {
            std::cerr << "[ERROR] Cannot open file: "
                      << samples[i].fileName << std::endl;

            if(f)
                f->Close();

            continue;
        }

        TH1D* h = dynamic_cast<TH1D*>(f->Get("h_mT"));

        if(!h)
        {
            std::cerr << "[ERROR] Cannot find h_mT in file: "
                      << samples[i].fileName << std::endl;

            f->Close();
            continue;
        }

        if(!hReference)
        {
            hReference =
                dynamic_cast<TH1D*>(h->Clone("h_mT_reference"));

            hReference->SetDirectory(0);
            referenceFile = samples[i].fileName;

            std::cout << "\n[REFERENCE]" << std::endl;
            PrintMTBinning(hReference, referenceFile);
        }
        else if(!SameMTBinning(hReference, h))
        {
            std::cerr << "\n[BINNING MISMATCH]" << std::endl;
            std::cerr << "Reference : "
                      << referenceFile << std::endl;
            std::cerr << "Different : "
                      << samples[i].fileName << std::endl;

            PrintMTBinning(h, samples[i].fileName);
        }

        f->Close();
    }

    delete hReference;

    std::cout << "\n====================================" << std::endl;
    std::cout << "h_mT binning check finished" << std::endl;
    std::cout << "====================================" << std::endl;
}

void Draw_Hist(const double xsec200 = 8.78288653e00,
               const double xsec500 = 2.76542446e-01,
               const double xsec1000 = 1.58702653e-02,
               const double xsec2000 = 4.15066636e-04,
               const double DYxsec200 = 2.78,
               const double DYxsec400 = 0.15,
               const double DYxsec500 = 0.084,
               const double DYxsec700 = 0.013,
               const double DYxsec800 = 0.011,
               const double DYxsec1000 = 0.006,
               const double DYxsec1500 = 0.00081,
               const double DYxsec2000 = 0.0002,
               const double TTTo2L2Nuxsec = 89.31,
               const double TTToSemiLeptonicxsec = 367.8,
               const double WWTo1L1Nu2Qxsec = 51.65,
               const double WWTo4Q_4fxsec = 51.03,
               const double WWTo2L2Nuxsec = 11.09,
               const double WZTo1L1Nu2Q_4fxsec = 9.119,
               const double WZTo3LNuxsec = 5.213,
               const double WZTo2Q2Nu_4fxsec = 6.331,
               const double WZTo1L3Nu_4fxsec = 3.414,
               const double ZZxsec = 12.17,
               const double ST_s_channel_4f_leptonDecays_xsec = 10.32,
               const double ST_t_channel_top_4f_InclusiveDecays_xsec = 134.2,
               const double ST_t_channel_anitop_4f_InclusiveDecays_xsec = 80.0,
               const double ST_tW_top_5f_inclusiveDecays_xsec = 39.65,
               const double ST_tW_antitop_5f_inclusiveDecays_xsec = 39.65,
               const double WToTauNu_M_200_xsec = 8.78288653e00,
               const double WToTauNu_M_500_xsec = 2.76542446e-01,
               const double WToTauNu_M_1000_xsec = 1.58702653e-02,
               const double WToTauNu_M_2000_xsec = 4.15066636e-04,
               const double lumi = 19520.0,
               const double sumW200 = 3250000,
               const double sumW500 = 750000,
               const double sumW1000 = 250000,
               const double sumW2000 = 50000,
               const double sumWDY200 = 2673760.516799,
               const double sumWDY400 = 44344.277472,
               const double sumWDY500 = 24718.574857,
               const double sumWDY700 = 3951.644255,
               const double sumWDY800 = 3345.560043,
               const double sumWDY1000 = 1782.846123,
               const double sumWDY1500 = 244.385548,
               const double sumWDY2000 = 62.650483,
               const double sumWTTTo2L2Nu = 2.70453e+09,
               const double sumWTTToSemiLeptonic = 3.97723e+10,
               const double sumWWWTo1L1Nu2Q = 1.69101e+09,
               const double sumWWWTo2L2Nu = 3.34565e+07,
               const double sumWWWTo4Q_4f = 1.67235e+09,
               const double sumWWZTo1L1Nu2Q_4f = 5.58378e+07,
               const double sumWWZTo1L3Nu_4f = 6.96859e+06,
               const double sumWWZTo2Q2Nu_4f = 4.84591e+07,
               const double sumWWZTo3LNu = 8.15165e+07,
               const double sumWZZ = 1282000,
               const double sumW_ST_s_channel_4f_leptonDecay = 1.95962e+07,
               const double sumW_ST_t_channel_anitop_4f_InclusiveDecays = 1.98386e+09,
               const double sumW_ST_t_channel_top_4f_InclusiveDecays = 5.94814e+09,
               const double sumW_ST_tW_top_5f_inclusiveDecays = 7.46247e+07,
               const double sumW_ST_tW_anitop_5f_inclusiveDecays = 7.47663e+07,
               const double sumW_WToTauNu_M_200 = 1250000,
               const double sumW_WToTauNu_M_500 = 750000,
               const double sumW_WToTauNu_M_1000 = 250000,
               const double sumW_WToTauNu_M_2000 = 50000)
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    double kW200 = (xsec200 * lumi) / sumW200;
    double kW500 = (xsec500 * lumi) / sumW500;
    double kW1000 = (xsec1000 * lumi) / sumW1000;
    double kW2000 = (xsec2000 * lumi) / sumW2000;

    double kDY200 = (DYxsec200 * lumi) / sumWDY200;
    double kDY400 = (DYxsec400 * lumi) / sumWDY400;
    double kDY500 = (DYxsec500 * lumi) / sumWDY500;
    double kDY700 = (DYxsec700 * lumi) / sumWDY700;
    double kDY800 = (DYxsec800 * lumi) / sumWDY800;
    double kDY1000 = (DYxsec1000 * lumi) / sumWDY1000;
    double kDY1500 = (DYxsec1500 * lumi) / sumWDY1500;
    double kDY2000 = (DYxsec2000 * lumi) / sumWDY2000;

    double kTTTo2L2Nu = (TTTo2L2Nuxsec * lumi) / sumWTTTo2L2Nu;
    double kTTToSemiLeptonic = (TTToSemiLeptonicxsec * lumi) / sumWTTToSemiLeptonic;

    double kWWTo1L1Nu2Q = (WWTo1L1Nu2Qxsec * lumi) / sumWWWTo1L1Nu2Q;
    double kWWTo2L2Nu = (WWTo2L2Nuxsec * lumi) / sumWWWTo2L2Nu;
    double kWWTo4Q_4f = (WWTo4Q_4fxsec * lumi) / sumWWWTo4Q_4f;
    double kWZTo1L1Nu2Q_4f = (WZTo1L1Nu2Q_4fxsec * lumi) / sumWWZTo1L1Nu2Q_4f;
    double kWZTo1L3Nu_4f = (WZTo1L3Nu_4fxsec * lumi) / sumWWZTo1L3Nu_4f;
    double kWZTo2Q2Nu_4f = (WZTo2Q2Nu_4fxsec * lumi) / sumWWZTo2Q2Nu_4f;
    double kWZTo3LNu = (WZTo3LNuxsec * lumi) / sumWWZTo3LNu;
    double kZZ = (ZZxsec * lumi) / sumWZZ;

    double kST_s_channel_4f_leptonDecays = (ST_s_channel_4f_leptonDecays_xsec * lumi) / sumW_ST_s_channel_4f_leptonDecay;
    double kST_t_channel_anitop_4f_InclusiveDecays = (ST_t_channel_anitop_4f_InclusiveDecays_xsec * lumi) / sumW_ST_t_channel_anitop_4f_InclusiveDecays;
    double kST_t_channel_top_4f_InclusiveDecays = (ST_t_channel_top_4f_InclusiveDecays_xsec * lumi) / sumW_ST_t_channel_top_4f_InclusiveDecays;
    double kST_tW_antitop_5f_inclusiveDecays = (ST_tW_antitop_5f_inclusiveDecays_xsec * lumi) / sumW_ST_tW_anitop_5f_inclusiveDecays;
    double kST_tW_top_5f_inclusiveDecays = (ST_tW_top_5f_inclusiveDecays_xsec * lumi) / sumW_ST_tW_top_5f_inclusiveDecays;
    
    double kWToTauNu_M_200 = (WToTauNu_M_200_xsec * lumi) / sumW_WToTauNu_M_200;
    double kWToTauNu_M_500 = (WToTauNu_M_500_xsec * lumi) / sumW_WToTauNu_M_500;
    double kWToTauNu_M_1000 = (WToTauNu_M_1000_xsec * lumi) / sumW_WToTauNu_M_1000;
    double kWToTauNu_M_2000 = (WToTauNu_M_2000_xsec * lumi) / sumW_WToTauNu_M_2000;

    std::vector<SampleInfo> samples = 
    {
        {"Muon_data_2016_preVFP.root", "Data", kBlack, 1.0, true},

        {"Muon_WToTauNu_M-200_2016_preVFP.root", "W#rightarrow#tau#nu", kGreen+2, kWToTauNu_M_200, false},
        {"Muon_WToTauNu_M-500_2016_preVFP.root", "W#rightarrow#tau#nu", kGreen+2, kWToTauNu_M_500, false},
        {"Muon_WToTauNu_M-1000_2016_preVFP.root", "W#rightarrow#tau#nu", kGreen+2, kWToTauNu_M_1000, false},
        {"Muon_WToTauNu_M-2000_2016_preVFP.root", "W#rightarrow#tau#nu", kGreen+2, kWToTauNu_M_2000, false},

        {"Muon_DYMC_200_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY200, false},
        {"Muon_DYMC_400_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY400, false},
        {"Muon_DYMC_500_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY500, false},
        {"Muon_DYMC_700_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY700, false},
        {"Muon_DYMC_800_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY800, false},
        {"Muon_DYMC_1000_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY1000, false},
        {"Muon_DYMC_1500_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY1500, false},
        {"Muon_DYMC_2000_2016_preVFP.root", "Z#rightarrow#mu#mu", kYellow, kDY2000, false},

        {"Muon_ST_s-channel_4f_leptonDecays_2016_preVFP.root", "single top", kBlue, kST_s_channel_4f_leptonDecays, false},
        {"Muon_ST_t-channel_antitop_4f_InclusiveDecays_2016_preVFP.root", "single top", kBlue, kST_t_channel_anitop_4f_InclusiveDecays, false},
        {"Muon_ST_t-channel_top_4f_InclusiveDecays_2016_preVFP.root", "single top", kBlue, kST_t_channel_top_4f_InclusiveDecays, false},
        {"Muon_ST_tW_antitop_5f_inclusiveDecays_2016_preVFP.root", "single top", kBlue, kST_tW_antitop_5f_inclusiveDecays, false},
        {"Muon_ST_tW_top_5f_inclusiveDecays_2016_preVFP.root", "single top", kBlue, kST_tW_top_5f_inclusiveDecays, false},

        {"Muon_QCD_2016_preVFP.root", "QCD", kMagenta, 1.0, false},

        {"Muon_WWTo1L1Nu2Q_2016_preVFP.root", "Diboson", kOrange, kWWTo1L1Nu2Q, false},
        {"Muon_WWTo2L2Nu_2016_preVFP.root", "Diboson", kOrange, kWWTo2L2Nu, false},
        {"Muon_WWTo4Q_4f_2016_preVFP.root", "Diboson", kOrange, kWWTo4Q_4f, false},
        {"Muon_WZTo1L1Nu2Q_4f_2016_preVFP.root", "Diboson", kOrange, kWZTo1L1Nu2Q_4f, false},
        {"Muon_WZTo1L3Nu_4f_2016_preVFP.root", "Diboson", kOrange, kWZTo1L3Nu_4f, false},
        {"Muon_WZTo2Q2Nu_4f_2016_preVFP.root", "Diboson", kOrange, kWZTo2Q2Nu_4f, false},
        {"Muon_WZTo3LNu_2016_preVFP.root", "Diboson", kOrange, kWZTo3LNu, false},
        {"Muon_ZZ_2016_preVFP.root", "Diboson", kOrange, kZZ, false},

        {"Muon_TTTo2L2Nu_2016_preVFP.root", "t#bar{t}", kRed, kTTTo2L2Nu, false},
        {"Muon_TTToSemiLeptonic_2016_preVFP.root", "t#bar{t}", kRed, kTTToSemiLeptonic, false},

        {"Muon_MC_200_2016_preVFP.root", "W#rightarrow#mu#nu", kGreen+1, kW200, false},
        {"Muon_MC_500_2016_preVFP.root", "W#rightarrow#mu#nu", kGreen+1, kW500, false},
        {"Muon_MC_1000_2016_preVFP.root", "W#rightarrow#mu#nu", kGreen+1, kW1000, false},
        {"Muon_MC_2000_2016_preVFP.root", "W#rightarrow#mu#nu", kGreen+1, kW2000, false}
    };

    CheckMTBinning(samples);

    DrawOneVariable(samples, "h_mT", "m_{T}", "m_{T} [GeV]", "Data_MC_2016_preVFP_mT.png", true, true, true);
    DrawOneVariable(samples, "h_Muon_pt", "Muon p_{T}", "p_{T} [GeV]", "Data_MC_2016_preVFP_pt.png", true, true, true);
    DrawOneVariable(samples, "h_Muon_eta", "Muon #eta", "#eta", "Data_MC_2016_preVFP_eta.png", false, false, false);
    DrawOneVariable(samples, "h_Muon_phi", "Muon #phi", "#phi", "Data_MC_2016_preVFP_phi.png", false, false, false);
    DrawOneVariable(samples, "h_MET_pt", "MET p_{T}", "MET p_{T}", "Data_MC_2016_preVFP_MET_pt.png", true, true, true);
    DrawOneVariable(samples, "h_MET_phi", "MET #phi", "#phi", "Data_MC_2016_preVFP_MET_phi.png", false, false, false);
}