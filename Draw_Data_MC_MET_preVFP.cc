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

void Draw_Data_MC_MET_preVFP(const double xsec200 = 8.78288653e00,
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
                             const double sumWTTToSemiLeptonic = 3.91287e+10,
                             const double sumWWWTo1L1Nu2Q = 1.69101e+09,
                             const double sumWWWTo2L2Nu = 3.34565e+07,
                             const double sumWWWTo4Q_4f = 1.58744e+09,
                             const double sumWWZTo1L1Nu2Q_4f = 5.58378e+07,
                             const double sumWWZTo1L3Nu_4f = 6.96859e+06,
                             const double sumWWZTo2Q2Nu_4f = 4.84591e+07,
                             const double sumWWZTo3LNu = 8.06709e+07,
                             const double sumWWZ = 7934000,
                             const double sumWZZ = 1282000,
                             const double sumW_ST_s_channel_4f_leptonDecay = 1.95962e+07,
                             const double sumW_ST_t_channel_anitop_4f_InclusiveDecays = 1.98386e+09,
                             const double sumW_ST_t_channel_top_4f_InclusiveDecays = 5.94814e+09,
                             const double sumW_ST_tW_top_5f_inclusiveDecays = 7.46247e+07,
                             const double sumW_ST_tW_anitop_5f_inclusiveDecays = 7.47663e+07)
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    TFile *fData = TFile::Open("Muon_data_2016_preVFP.root");
    TFile *fMC200 = TFile::Open("Muon_MC_200_2016_preVFP.root");
    TFile *fMC500 = TFile::Open("Muon_MC_500_2016_preVFP.root");
    TFile *fMC1000 = TFile::Open("Muon_MC_1000_2016_preVFP.root");
    TFile *fMC2000 = TFile::Open("Muon_MC_2000_2016_preVFP.root");
    TFile *fDYMC200 = TFile::Open("Muon_DYMC_200_2016_preVFP.root");
    TFile *fDYMC400 = TFile::Open("Muon_DYMC_400_2016_preVFP.root");
    TFile *fDYMC500 = TFile::Open("Muon_DYMC_500_2016_preVFP.root");
    TFile *fDYMC700 = TFile::Open("Muon_DYMC_700_2016_preVFP.root");
    TFile *fDYMC800 = TFile::Open("Muon_DYMC_800_2016_preVFP.root");
    TFile *fDYMC1000 = TFile::Open("Muon_DYMC_1000_2016_preVFP.root");
    TFile *fDYMC1500 = TFile::Open("Muon_DYMC_1500_2016_preVFP.root");
    TFile *fDYMC2000 = TFile::Open("Muon_DYMC_2000_2016_preVFP.root");
    TFile *fTTTo2L2Nu = TFile::Open("Muon_TTTo2L2Nu_2016_preVFP.root");
    TFile *fTTToSemiLeptonic = TFile::Open("Muon_TTToSemiLeptonic_2016_preVFP.root");
    TFile *fWWTo1L1Nu2Q = TFile::Open("Muon_WWTo1L1Nu2Q_2016_preVFP.root");
    TFile *fWWTo2L2Nu = TFile::Open("Muon_WWTo2L2Nu_2016_preVFP.root");
    TFile *fWWTo4Q_4f = TFile::Open("Muon_WWTo4Q_4f_2016_preVFP.root");
    TFile *fWZTo1L1Nu2Q_4f = TFile::Open("Muon_WZTo1L1Nu2Q_4f_2016_preVFP.root");
    TFile *fWZTo1L3Nu_4f = TFile::Open("Muon_WZTo1L3Nu_4f_2016_preVFP.root");
    TFile *fWZTo2Q2Nu_4f = TFile::Open("Muon_WZTo2Q2Nu_4f_2016_preVFP.root");
    TFile *fWZTo3LNu = TFile::Open("Muon_WZTo3LNu_2016_preVFP.root");
    TFile *fZZ = TFile::Open("Muon_ZZ_2016_preVFP.root");
    TFile *fST_s_channel_4f_leptonDecays = TFile::Open("Muon_ST_s-channel_4f_leptonDecays_2016_preVFP.root");
    TFile *fST_t_channel_anitop_4f_InclusiveDecays = TFile::Open("Muon_ST_t-channel_anitop_4f_InclusiveDecays_2016_preVFP.root");
    TFile *fST_t_channel_top_4f_InclusiveDecays = TFile::Open("Muon_ST_t-channel_top_4f_InclusiveDecays_2016_preVFP.root");
    TFile *fST_tW_antitop_5f_inclusiveDecays = TFile::Open("Muon_ST_tW_antitop_5f_inclusiveDecays_2016_preVFP.root");
    TFile *fST_tW_top_5f_inclusiveDecays = TFile::Open("Muon_ST_tW_top_5f_inclusiveDecays_2016_preVFP.root");

    TH1D *hData = (TH1D*)fData->Get("h_MET_pt");
    TH1D *hMC200 = (TH1D*)fMC200->Get("h_MET_pt");
    TH1D *hMC500 = (TH1D*)fMC500->Get("h_MET_pt");
    TH1D *hMC1000 = (TH1D*)fMC1000->Get("h_MET_pt");
    TH1D *hMC2000 = (TH1D*)fMC2000->Get("h_MET_pt");
    TH1D *hDYMC200 = (TH1D*)fDYMC200->Get("h_MET_pt");
    TH1D *hDYMC400 = (TH1D*)fDYMC400->Get("h_MET_pt");
    TH1D *hDYMC500 = (TH1D*)fDYMC500->Get("h_MET_pt");
    TH1D *hDYMC700 = (TH1D*)fDYMC700->Get("h_MET_pt");
    TH1D *hDYMC800 = (TH1D*)fDYMC800->Get("h_MET_pt");
    TH1D *hDYMC1000 = (TH1D*)fDYMC1000->Get("h_MET_pt");
    TH1D *hDYMC1500 = (TH1D*)fDYMC1500->Get("h_MET_pt");
    TH1D *hDYMC2000 = (TH1D*)fDYMC2000->Get("h_MET_pt");
    TH1D *hTTTo2L2Nu = (TH1D*)fTTTo2L2Nu->Get("h_MET_pt");
    TH1D *hTTToSemiLeptonic = (TH1D*)fTTToSemiLeptonic->Get("h_MET_pt");
    TH1D *hWWTo1L1Nu2Q = (TH1D*)fWWTo1L1Nu2Q->Get("h_MET_pt");
    TH1D *hWWTo2L2Nu = (TH1D*)fWWTo2L2Nu->Get("h_MET_pt");
    TH1D *hWWTo4Q_4f = (TH1D*)fWWTo4Q_4f->Get("h_MET_pt");
    TH1D *hWZTo1L1Nu2Q_4f = (TH1D*)fWZTo1L1Nu2Q_4f->Get("h_MET_pt");
    TH1D *hWZTo1L3Nu_4f = (TH1D*)fWZTo1L3Nu_4f->Get("h_MET_pt");
    TH1D *hWZTo2Q2Nu_4f = (TH1D*)fWZTo2Q2Nu_4f->Get("h_MET_pt");
    TH1D *hWZTo3LNu = (TH1D*)fWZTo3LNu->Get("h_MET_pt");
    TH1D *hZZ = (TH1D*)fZZ->Get("h_MET_pt");
    TH1D *hST_s_channel_4f_leptonDecays = (TH1D*)fST_s_channel_4f_leptonDecays->Get("h_MET_pt");
    TH1D *hST_t_channel_anitop_4f_InclusiveDecays = (TH1D*)fST_t_channel_anitop_4f_InclusiveDecays->Get("h_MET_pt");
    TH1D *hST_t_channel_top_4f_InclusiveDecays = (TH1D*)fST_t_channel_top_4f_InclusiveDecays->Get("h_MET_pt");
    TH1D *hST_tW_antitop_5f_inclusiveDecays = (TH1D*)fST_tW_antitop_5f_inclusiveDecays->Get("h_MET_pt");
    TH1D *hST_tW_top_5f_inclusiveDecays = (TH1D*)fST_tW_top_5f_inclusiveDecays->Get("h_MET_pt");

    hData->SetDirectory(0);
    hMC200->SetDirectory(0);
    hMC500->SetDirectory(0);
    hMC1000->SetDirectory(0);
    hMC2000->SetDirectory(0);
    hDYMC200->SetDirectory(0);
    hDYMC400->SetDirectory(0);
    hDYMC500->SetDirectory(0);
    hDYMC700->SetDirectory(0);
    hDYMC800->SetDirectory(0);
    hDYMC1000->SetDirectory(0);
    hDYMC1500->SetDirectory(0);
    hDYMC2000->SetDirectory(0);
    hTTTo2L2Nu->SetDirectory(0);
    hTTToSemiLeptonic->SetDirectory(0);
    hWWTo1L1Nu2Q->SetDirectory(0);
    hWWTo2L2Nu->SetDirectory(0);
    hWWTo4Q_4f->SetDirectory(0);
    hWZTo1L1Nu2Q_4f->SetDirectory(0);
    hWZTo1L3Nu_4f->SetDirectory(0);
    hWZTo2Q2Nu_4f->SetDirectory(0);
    hWZTo3LNu->SetDirectory(0);
    hZZ->SetDirectory(0);
    hST_s_channel_4f_leptonDecays->SetDirectory(0);
    hST_t_channel_anitop_4f_InclusiveDecays->SetDirectory(0);
    hST_t_channel_top_4f_InclusiveDecays->SetDirectory(0);
    hST_tW_antitop_5f_inclusiveDecays->SetDirectory(0);
    hST_tW_top_5f_inclusiveDecays->SetDirectory(0);

    fData->Close();
    fMC200->Close();
    fMC500->Close();
    fMC1000->Close();
    fMC2000->Close();
    fDYMC200->Close();
    fDYMC400->Close();
    fDYMC500->Close();
    fDYMC700->Close();
    fDYMC800->Close();
    fDYMC1000->Close();
    fDYMC1500->Close();
    fDYMC2000->Close();
    fTTTo2L2Nu->Close();
    fTTToSemiLeptonic->Close();
    fWWTo1L1Nu2Q->Close();
    fWWTo2L2Nu->Close();
    fWWTo4Q_4f->Close();
    fWZTo1L1Nu2Q_4f->Close();
    fWZTo1L3Nu_4f->Close();
    fWZTo2Q2Nu_4f->Close();
    fWZTo3LNu->Close();
    fZZ->Close();
    fST_s_channel_4f_leptonDecays->Close();
    fST_t_channel_anitop_4f_InclusiveDecays->Close();
    fST_t_channel_top_4f_InclusiveDecays->Close();
    fST_tW_antitop_5f_inclusiveDecays->Close();
    fST_tW_top_5f_inclusiveDecays->Close();

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

    hMC200->Scale(kW200);
    hMC500->Scale(kW500);
    hMC1000->Scale(kW1000);
    hMC2000->Scale(kW2000);
    hDYMC200->Scale(kDY200);
    hDYMC400->Scale(kDY400);
    hDYMC500->Scale(kDY500);
    hDYMC700->Scale(kDY700);
    hDYMC800->Scale(kDY800);
    hDYMC1000->Scale(kDY1000);
    hDYMC1500->Scale(kDY1500);
    hDYMC2000->Scale(kDY2000);
    hTTTo2L2Nu->Scale(kTTTo2L2Nu);
    hTTToSemiLeptonic->Scale(kTTToSemiLeptonic);
    hWWTo1L1Nu2Q->Scale(kWWTo1L1Nu2Q);
    hWWTo2L2Nu->Scale(kWWTo2L2Nu);
    hWWTo4Q_4f->Scale(kWWTo4Q_4f);
    hWZTo1L1Nu2Q_4f->Scale(kWZTo1L1Nu2Q_4f);
    hWZTo1L3Nu_4f->Scale(kWZTo1L3Nu_4f);
    hWZTo2Q2Nu_4f->Scale(kWZTo2Q2Nu_4f);
    hWZTo3LNu->Scale(kWZTo3LNu);
    hZZ->Scale(kZZ);
    hST_s_channel_4f_leptonDecays->Scale(kST_s_channel_4f_leptonDecays);
    hST_t_channel_anitop_4f_InclusiveDecays->Scale(kST_t_channel_anitop_4f_InclusiveDecays);
    hST_t_channel_top_4f_InclusiveDecays->Scale(kST_t_channel_top_4f_InclusiveDecays);
    hST_tW_antitop_5f_inclusiveDecays->Scale(kST_tW_antitop_5f_inclusiveDecays);
    hST_tW_top_5f_inclusiveDecays->Scale(kST_tW_top_5f_inclusiveDecays);

    TH1D *hMCsum = (TH1D*)hMC200->Clone("hMCsum");
    hMCsum->Add(hMC500);
    hMCsum->Add(hMC1000);
    hMCsum->Add(hMC2000);
    hMCsum->Add(hDYMC200);
    hMCsum->Add(hDYMC400);
    hMCsum->Add(hDYMC500);
    hMCsum->Add(hDYMC700);
    hMCsum->Add(hDYMC800);
    hMCsum->Add(hDYMC1000);
    hMCsum->Add(hDYMC1500);
    hMCsum->Add(hDYMC2000);
    hMCsum->Add(hTTTo2L2Nu);
    hMCsum->Add(hTTToSemiLeptonic);
    hMCsum->Add(hWWTo1L1Nu2Q);
    hMCsum->Add(hWWTo2L2Nu);
    hMCsum->Add(hWWTo4Q_4f);
    hMCsum->Add(hWZTo1L1Nu2Q_4f);
    hMCsum->Add(hWZTo1L3Nu_4f);
    hMCsum->Add(hWZTo2Q2Nu_4f);
    hMCsum->Add(hWZTo3LNu);
    hMCsum->Add(hZZ);
    hMCsum->Add(hST_s_channel_4f_leptonDecays);
    hMCsum->Add(hST_t_channel_anitop_4f_InclusiveDecays);
    hMCsum->Add(hST_t_channel_top_4f_InclusiveDecays);
    hMCsum->Add(hST_tW_antitop_5f_inclusiveDecays);
    hMCsum->Add(hST_tW_top_5f_inclusiveDecays);

    const double ymax = std::max(hMCsum->GetMaximum(), hData->GetMaximum());

    TCanvas *c = new TCanvas("c", "MET p_{T}", 900, 800);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.18, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    THStack *stack = new THStack("hstack", "MET ;MET p_{T} [GeV]; Events");

    hMC200->SetFillColor(kRed);
    hMC200->SetLineColor(kRed);
    hMC500->SetFillColor(kRed);
    hMC500->SetLineColor(kRed);
    hMC1000->SetLineColor(kRed);
    hMC1000->SetFillColor(kRed);
    hMC2000->SetLineColor(kRed);
    hMC2000->SetFillColor(kRed);
    hDYMC200->SetLineColor(kGreen+1);
    hDYMC200->SetFillColor(kGreen+1);
    hDYMC400->SetLineColor(kGreen+1);
    hDYMC400->SetFillColor(kGreen+1);
    hDYMC500->SetLineColor(kGreen+1);
    hDYMC500->SetFillColor(kGreen+1);
    hDYMC700->SetLineColor(kGreen+1);
    hDYMC700->SetFillColor(kGreen+1);
    hDYMC800->SetLineColor(kGreen+1);
    hDYMC800->SetFillColor(kGreen+1);
    hDYMC1000->SetLineColor(kGreen+1);
    hDYMC1000->SetFillColor(kGreen+1);
    hDYMC1500->SetLineColor(kGreen+1);
    hDYMC1500->SetFillColor(kGreen+1);
    hDYMC2000->SetLineColor(kGreen+1);
    hDYMC2000->SetFillColor(kGreen+1);
    hTTTo2L2Nu->SetLineColor(kOrange);
    hTTTo2L2Nu->SetFillColor(kOrange);
    hTTToSemiLeptonic->SetLineColor(kOrange);
    hTTToSemiLeptonic->SetFillColor(kOrange);
    hWWTo1L1Nu2Q->SetLineColor(kMagenta+1);
    hWWTo1L1Nu2Q->SetFillColor(kMagenta+1);
    hWWTo2L2Nu->SetLineColor(kMagenta+1);
    hWWTo2L2Nu->SetFillColor(kMagenta+1);
    hWWTo4Q_4f->SetLineColor(kMagenta+1);
    hWWTo4Q_4f->SetFillColor(kMagenta+1);
    hWZTo1L1Nu2Q_4f->SetLineColor(kMagenta+1);
    hWZTo1L1Nu2Q_4f->SetFillColor(kMagenta+1);
    hWZTo1L3Nu_4f->SetLineColor(kMagenta+1);
    hWZTo1L3Nu_4f->SetFillColor(kMagenta+1);
    hWZTo2Q2Nu_4f->SetLineColor(kMagenta+1);
    hWZTo2Q2Nu_4f->SetFillColor(kMagenta+1);
    hWZTo3LNu->SetLineColor(kMagenta+1);
    hWZTo3LNu->SetFillColor(kMagenta+1);
    hZZ->SetLineColor(kMagenta+1);
    hZZ->SetFillColor(kMagenta+1);
    hST_s_channel_4f_leptonDecays->SetLineColor(kBlue);
    hST_s_channel_4f_leptonDecays->SetFillColor(kBlue);
    hST_t_channel_anitop_4f_InclusiveDecays->SetLineColor(kBlue);
    hST_t_channel_anitop_4f_InclusiveDecays->SetFillColor(kBlue);
    hST_t_channel_top_4f_InclusiveDecays->SetLineColor(kBlue);
    hST_t_channel_top_4f_InclusiveDecays->SetFillColor(kBlue);
    hST_tW_antitop_5f_inclusiveDecays->SetLineColor(kBlue);
    hST_tW_antitop_5f_inclusiveDecays->SetFillColor(kBlue);
    hST_tW_top_5f_inclusiveDecays->SetLineColor(kBlue);
    hST_tW_top_5f_inclusiveDecays->SetFillColor(kBlue);

    stack->Add(hDYMC200);
    stack->Add(hDYMC400);
    stack->Add(hDYMC500);
    stack->Add(hDYMC700);
    stack->Add(hDYMC800);
    stack->Add(hDYMC1000);
    stack->Add(hDYMC1500);
    stack->Add(hDYMC2000);
    stack->Add(hST_s_channel_4f_leptonDecays);
    stack->Add(hST_t_channel_anitop_4f_InclusiveDecays);
    stack->Add(hST_t_channel_top_4f_InclusiveDecays);
    stack->Add(hST_tW_antitop_5f_inclusiveDecays);
    stack->Add(hST_tW_top_5f_inclusiveDecays);
    stack->Add(hWWTo1L1Nu2Q);
    stack->Add(hWWTo2L2Nu);
    stack->Add(hWWTo4Q_4f);
    stack->Add(hWZTo1L1Nu2Q_4f);
    stack->Add(hWZTo1L3Nu_4f);
    stack->Add(hWZTo2Q2Nu_4f);
    stack->Add(hWZTo3LNu);
    stack->Add(hZZ);
    stack->Add(hTTTo2L2Nu);
    stack->Add(hTTToSemiLeptonic);
    stack->Add(hMC200);
    stack->Add(hMC500);
    stack->Add(hMC1000);
    stack->Add(hMC2000);
    stack->Draw("HIST");

    stack->SetMinimum(1e-3);
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
    leg->AddEntry(hTTTo2L2Nu, "t#bar{t}", "f");
    leg->AddEntry(hZZ, "Diboson", "f");
    leg->AddEntry(hST_s_channel_4f_leptonDecays, "single top", "f");
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

    hRatio->GetXaxis()->SetTitle("MET p_{T} [GeV]");
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

    c->SaveAs("Data_MC_2016_preVFP_MET_pt.png");
}