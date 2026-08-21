#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TColor.h>


// ============================================================
// Sample information
// ============================================================

struct SampleInfo
{
    std::string fileName;
    double xsec;   // pb
    double sumW;
};


// ============================================================
// Normalize each MC sample and sum
//
// Assumption:
// histogram itself was filled with MC event weights/genWeight etc.
// Final sample normalization = lumi * xsec / sumW
// ============================================================

TH1D* SumNormalizedHistograms(
    const std::vector<SampleInfo>& samples,
    const std::string& histName,
    const std::string& newName,
    double lumi)
{
    TH1D* hSum = nullptr;

    for(const auto& sample : samples)
    {
        TFile* f = TFile::Open(sample.fileName.c_str(), "READ");

        if(!f || f->IsZombie())
        {
            std::cerr
                << "[ERROR] Cannot open "
                << sample.fileName
                << std::endl;

            if(f)
            {
                f->Close();
                delete f;
            }

            continue;
        }


        TH1D* h =
            dynamic_cast<TH1D*>(f->Get(histName.c_str()));

        if(!h)
        {
            std::cerr
                << "[ERROR] Cannot find "
                << histName
                << " in "
                << sample.fileName
                << std::endl;

            f->Close();
            delete f;
            continue;
        }


        if(sample.sumW == 0.0)
        {
            std::cerr
                << "[ERROR] sumW = 0 for "
                << sample.fileName
                << std::endl;

            f->Close();
            delete f;
            continue;
        }


        // ----------------------------------------------------
        // Clone original histogram
        // ----------------------------------------------------

        std::string tempName =
            newName + "_temp_" +
            std::to_string(reinterpret_cast<std::uintptr_t>(h));

        TH1D* hTemp =
            dynamic_cast<TH1D*>(h->Clone(tempName.c_str()));

        hTemp->SetDirectory(nullptr);


        // ----------------------------------------------------
        // MC normalization
        // ----------------------------------------------------

        double scale =
            lumi * sample.xsec / sample.sumW;

        hTemp->Scale(scale);


        std::cout
            << sample.fileName
            << "\n"
            << "    xsec   = " << sample.xsec << " pb\n"
            << "    sumW   = " << sample.sumW << "\n"
            << "    scale  = " << scale << "\n"
            << "    yield  = "
            << hTemp->Integral(0, hTemp->GetNbinsX() + 1)
            << std::endl;


        // ----------------------------------------------------
        // Add to process histogram
        // ----------------------------------------------------

        if(!hSum)
        {
            hSum =
                dynamic_cast<TH1D*>(
                    hTemp->Clone(newName.c_str())
                );

            hSum->SetDirectory(nullptr);
        }
        else
        {
            // basic binning consistency check
            if(hSum->GetNbinsX() != hTemp->GetNbinsX())
            {
                std::cerr
                    << "[ERROR] Different binning in "
                    << sample.fileName
                    << std::endl;

                delete hTemp;

                f->Close();
                delete f;

                continue;
            }

            hSum->Add(hTemp);
        }


        delete hTemp;

        f->Close();
        delete f;
    }

    return hSum;
}


// ============================================================
// Main plotting macro
// ============================================================

void Draw_QCD_CR_Validation()
{
    gStyle->SetOptStat(0);

    const double lumi = 19520.0; // pb^-1


    // ========================================================
    // Data
    //
    // 실제 Data fake-rate output filename으로 바꿔줘.
    // ========================================================

    TFile* fData =
        TFile::Open("fake_mu_merged.root", "READ");


    // ========================================================
    // Data-driven QCD prediction
    //
    // 실제 QCD CR output filename으로 바꿔줘.
    // ========================================================

    TFile* fQCD =
        TFile::Open("Muon_QCD_CR_2016_preVFP.root", "READ");


    if(!fData || fData->IsZombie())
    {
        std::cerr
            << "[ERROR] Cannot open Data file."
            << std::endl;

        return;
    }

    if(!fQCD || fQCD->IsZombie())
    {
        std::cerr
            << "[ERROR] Cannot open QCD file."
            << std::endl;

        return;
    }


    // ========================================================
    // W -> mu nu
    // ========================================================

    std::vector<SampleInfo> WMuSamples = {

        {
            "fake_mu_WToMuNu_M-200_merged.root",
            8.78288653e+00,
            3250000
        },

        {
            "fake_mu_WToMuNu_M-500_merged.root",
            2.76542446e-01,
            750000
        },

        {
            "fake_mu_WToMuNu_M-1000_merged.root",
            1.58702653e-02,
            250000
        },

        {
            "fake_mu_WToMuNu_M-2000_merged.root",
            4.15066636e-04,
            50000
        }
    };


    // ========================================================
    // W -> tau nu
    // ========================================================

    std::vector<SampleInfo> WTauSamples = {

        {
            "fake_mu_WToTauNu_M-200_merged.root",
            8.78288653e+00,
            1250000
        },

        {
            "fake_mu_WToTauNu_M-500_merged.root",
            2.76542446e-01,
            750000
        },

        {
            "fake_mu_WToTauNu_M-1000_merged.root",
            1.58702653e-02,
            250000
        },

        {
            "fake_mu_WToTauNu_M-2000_merged.root",
            4.15066636e-04,
            50000
        }
    };


    // ========================================================
    // DY
    // ========================================================

    std::vector<SampleInfo> DYSamples = {

        {
            "fake_mu_DYMC_200_merged.root",
            2.78,
            2673760.516799
        },

        {
            "fake_mu_DYMC_400_merged.root",
            0.15,
            44344.277472
        },

        {
            "fake_mu_DYMC_500_merged.root",
            0.084,
            24718.574857
        },

        {
            "fake_mu_DYMC_700_merged.root",
            0.013,
            3951.644255
        },

        {
            "fake_mu_DYMC_800_merged.root",
            0.011,
            3345.560043
        },

        {
            "fake_mu_DYMC_1000_merged.root",
            0.006,
            1782.846123
        },

        {
            "fake_mu_DYMC_1500_merged.root",
            0.00081,
            244.385548
        },

        {
            "fake_mu_DYMC_2000_merged.root",
            0.0002,
            62.650483
        }
    };


    // ========================================================
    // ttbar
    // ========================================================

    std::vector<SampleInfo> TTSamples = {

        {
            "fake_mu_TTTo2L2Nu_merged.root",
            89.31,
            2.70453e+09
        },

        {
            "fake_mu_TTToSemiLeptonic_merged.root",
            367.8,
            3.97723e+10
        }
    };


    // ========================================================
    // Single top
    // ========================================================

    std::vector<SampleInfo> STSamples = {

        {
            "fake_mu_ST_s-channel_4f_leptonDecays_merged.root",
            10.32,
            1.95962e+07
        },

        // 파일 이름이 실제로 anitop으로 되어 있으므로 그대로 사용
        {
            "fake_mu_ST_t-channel_anitop_4f_InclusiveDecays_merged.root",
            80.0,
            1.98386e+09
        },

        {
            "fake_mu_ST_t-channel_top_4f_InclusiveDecays_merged.root",
            134.2,
            5.94814e+09
        },

        {
            "fake_mu_ST_tW_antitop_5f_inclusiveDecays_merged.root",
            39.65,
            7.47663e+07
        },

        {
            "fake_mu_ST_tW_top_5f_inclusiveDecays_merged.root",
            39.65,
            7.46247e+07
        }
    };


    // ========================================================
    // Diboson
    // ========================================================

    std::vector<SampleInfo> VVSamples = {

        // WW
        {
            "fake_mu_WWTo1L1Nu2Q_merged.root",
            51.65,
            1.69101e+09
        },

        {
            "fake_mu_WWTo2L2Nu_merged.root",
            11.09,
            3.34565e+07
        },

        {
            "fake_mu_WWTo4Q_4f_merged.root",
            51.03,
            1.67235e+09
        },


        // WZ
        {
            "fake_mu_WZTo1L1Nu2Q_4f_merged.root",
            9.119,
            5.58378e+07
        },

        {
            "fake_mu_WZTo1L3Nu_4f_merged.root",
            3.414,
            6.96859e+06
        },

        {
            "fake_mu_WZTo2Q2Nu_4f_merged.root",
            6.331,
            4.84591e+07
        },

        {
            "fake_mu_WZTo3LNu_merged.root",
            5.213,
            8.15165e+07
        },


        // ZZ
        {
            "fake_mu_ZZ_merged.root",
            12.17,
            1282000
        }
    };


    // ========================================================
    // Plot variables
    // ========================================================

    struct PlotInfo
    {
        std::string numHist;
        std::string qcdHist;
        std::string xTitle;
        std::string output;
    };


    std::vector<PlotInfo> plots = {

        {
            "h_Muon_pt_num",
            "h_Muon_pt",
            "Muon p_{T} [GeV]",
            "CR_validation_Muon_pt.png"
        },

        {
            "h_Muon_eta_num",
            "h_Muon_eta",
            "Muon #eta",
            "CR_validation_Muon_eta.png"
        },

        {
            "h_Muon_phi_num",
            "h_Muon_phi",
            "Muon #phi",
            "CR_validation_Muon_phi.png"
        },

        {
            "h_MET_pt_num",
            "h_MET_pt",
            "MET p_{T} [GeV]",
            "CR_validation_MET_pt.png"
        },

        {
            "h_MET_phi_num",
            "h_MET_phi",
            "MET #phi",
            "CR_validation_MET_phi.png"
        }
    };


    // ========================================================
    // Plot loop
    // ========================================================

    for(const auto& p : plots)
    {
        const bool useLogX = (p.numHist == "h_Muon_pt_num");

        const bool useLogY =
            (p.numHist == "h_Muon_pt_num" ||
             p.numHist == "h_MET_pt_num");

        std::cout
            << "\n\n========================================"
            << std::endl;

        std::cout
            << "Variable: "
            << p.xTitle
            << std::endl;

        std::cout
            << "========================================"
            << std::endl;


        // ====================================================
        // Data numerator = observed tight CR
        // ====================================================

        TH1D* hDataOriginal =
            dynamic_cast<TH1D*>(
                fData->Get(p.numHist.c_str())
            );

        if(!hDataOriginal)
        {
            std::cerr
                << "[ERROR] Missing Data histogram: "
                << p.numHist
                << std::endl;

            continue;
        }


        TH1D* hData =
            dynamic_cast<TH1D*>(
                hDataOriginal->Clone("hData")
            );

        hData->SetDirectory(nullptr);


        // ====================================================
        // QCD prediction
        //
        // IMPORTANT:
        // Do NOT luminosity-normalize this histogram.
        // It is data-driven.
        // ====================================================

        TH1D* hQCDOriginal =
            dynamic_cast<TH1D*>(
                fQCD->Get(p.qcdHist.c_str())
            );

        if(!hQCDOriginal)
        {
            std::cerr
                << "[ERROR] Missing QCD histogram: "
                << p.qcdHist
                << std::endl;

            delete hData;
            continue;
        }


        TH1D* hQCD =
            dynamic_cast<TH1D*>(
                hQCDOriginal->Clone("hQCD")
            );

        hQCD->SetDirectory(nullptr);


        // ====================================================
        // Normalize and sum MC samples
        // ====================================================

        TH1D* hWMu =
            SumNormalizedHistograms(
                WMuSamples,
                p.numHist,
                "hWMu",
                lumi
            );


        TH1D* hWTau =
            SumNormalizedHistograms(
                WTauSamples,
                p.numHist,
                "hWTau",
                lumi
            );


        TH1D* hDY =
            SumNormalizedHistograms(
                DYSamples,
                p.numHist,
                "hDY",
                lumi
            );


        TH1D* hTT =
            SumNormalizedHistograms(
                TTSamples,
                p.numHist,
                "hTT",
                lumi
            );


        TH1D* hST =
            SumNormalizedHistograms(
                STSamples,
                p.numHist,
                "hST",
                lumi
            );


        TH1D* hVV =
            SumNormalizedHistograms(
                VVSamples,
                p.numHist,
                "hVV",
                lumi
            );


        if(!hWMu ||
           !hWTau ||
           !hDY ||
           !hTT ||
           !hST ||
           !hVV)
        {
            std::cerr
                << "[ERROR] One or more MC groups are missing."
                << std::endl;

            delete hData;
            delete hQCD;

            if(hWMu) delete hWMu;
            if(hWTau) delete hWTau;
            if(hDY) delete hDY;
            if(hTT) delete hTT;
            if(hST) delete hST;
            if(hVV) delete hVV;

            continue;
        }


        // ====================================================
        // Colors
        // ====================================================

        hQCD->SetFillColor(kMagenta);

        hWMu->SetFillColor(kGreen + 1);
        hWTau->SetFillColor(kGreen + 2);

        hDY->SetFillColor(kYellow);

        hTT->SetFillColor(kRed);

        hST->SetFillColor(kBlue);

        hVV->SetFillColor(kOrange);


        for(TH1D* h :
            {hQCD, hWMu, hWTau, hDY, hTT, hST, hVV})
        {
            h->SetLineColor(kBlack);
            h->SetLineWidth(1);
        }


        hData->SetMarkerStyle(20);
        hData->SetMarkerSize(1.0);
        hData->SetLineColor(kBlack);


        // ====================================================
        // Stack
        // ====================================================

        THStack* stack =
            new THStack("stack", "");


        stack->Add(hWMu);
        stack->Add(hWTau);
        stack->Add(hDY);
        stack->Add(hST);
        stack->Add(hVV);
        stack->Add(hTT);
        stack->Add(hQCD);


        // ====================================================
        // Total prediction
        // ====================================================

        TH1D* hTotal =
            dynamic_cast<TH1D*>(
                hQCD->Clone("hTotal")
            );


        hTotal->Add(hVV);
        hTotal->Add(hST);
        hTotal->Add(hTT);
        hTotal->Add(hDY);
        hTotal->Add(hWTau);
        hTotal->Add(hWMu);


        // ====================================================
        // Print yields
        // ====================================================

        auto IntegralAll =
            [](TH1D* h)
        {
            return h->Integral(
                0,
                h->GetNbinsX() + 1
            );
        };


        std::cout << "\n--------- CR yields ---------\n";

        std::cout
            << "Data             = "
            << IntegralAll(hData)
            << std::endl;

        std::cout
            << "W -> mu nu       = "
            << IntegralAll(hWMu)
            << std::endl;

        std::cout
            << "W -> tau nu      = "
            << IntegralAll(hWTau)
            << std::endl;

        std::cout
            << "DY               = "
            << IntegralAll(hDY)
            << std::endl;

        std::cout
            << "ttbar            = "
            << IntegralAll(hTT)
            << std::endl;

        std::cout
            << "Single top       = "
            << IntegralAll(hST)
            << std::endl;

        std::cout
            << "Diboson          = "
            << IntegralAll(hVV)
            << std::endl;

        std::cout
            << "Multijet         = "
            << IntegralAll(hQCD)
            << std::endl;

        std::cout
            << "------------------------------\n";

        std::cout
            << "Total prediction = "
            << IntegralAll(hTotal)
            << std::endl;


        double totalPrediction =
            IntegralAll(hTotal);

        if(totalPrediction != 0.0)
        {
            std::cout
                << "Data / Prediction= "
                << IntegralAll(hData) /
                   totalPrediction
                << std::endl;
        }


        // ====================================================
        // Canvas
        // ====================================================

        TCanvas* c =
            new TCanvas(
                ("c_" + p.numHist).c_str(),
                "c",
                900,
                800
            );


        TPad* pad1 =
            new TPad(
                ("pad1_" + p.numHist).c_str(),
                "pad1",
                0,
                0.18,
                1,
                1
            );


        pad1->SetLeftMargin(0.13);
        pad1->SetBottomMargin(0.02);

        if(useLogX)
            pad1->SetLogx();

        if(useLogY)
            pad1->SetLogy();

        pad1->Draw();
        pad1->cd();


        stack->Draw("HIST");

        stack->GetYaxis()->SetTitle(
            "Events"
        );


        stack->GetYaxis()->SetTitleSize(0.045);

        stack->GetYaxis()->SetLabelSize(0.04);

        stack->GetXaxis()->SetTitleSize(0);
        stack->GetXaxis()->SetLabelSize(0);
        stack->GetXaxis()->SetTickLength(0);


        double maxY =
            std::max(
                stack->GetMaximum(),
                hData->GetMaximum()
            );


        if(useLogY)
        {
            stack->SetMinimum(1e-3);

            if(maxY > 0)
                stack->SetMaximum(10.0 * maxY);
        }
        else if(maxY > 0)
        {
            stack->SetMinimum(0.0);
            stack->SetMaximum(1.45 * maxY);
        }


        // Data over stack
        hData->Draw("E1 SAME");


        // ====================================================
        // Legend
        // ====================================================

        TLegend* leg =
            new TLegend(
                0.61,
                0.52,
                0.88,
                0.88
            );


        leg->SetBorderSize(0);
        leg->SetFillStyle(0);


        leg->AddEntry(
            hData,
            "Data",
            "lep"
        );


        leg->AddEntry(
            hWMu,
            "W #rightarrow #mu#nu",
            "f"
        );


        leg->AddEntry(
            hWTau,
            "W #rightarrow #tau#nu",
            "f"
        );


        leg->AddEntry(
            hDY,
            "DY",
            "f"
        );


        leg->AddEntry(
            hTT,
            "t#bar{t}",
            "f"
        );


        leg->AddEntry(
            hST,
            "Single top",
            "f"
        );


        leg->AddEntry(
            hVV,
            "Diboson",
            "f"
        );


        leg->AddEntry(
            hQCD,
            "Multijet",
            "f"
        );


        leg->Draw();


        // ====================================================
        // Data / MC ratio
        // ====================================================

        pad1->Update();
        c->cd();

        TPad* pad2 =
            new TPad(
                ("pad2_" + p.numHist).c_str(),
                "pad2",
                0,
                0,
                1,
                0.18
            );


        pad2->SetLeftMargin(0.13);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.38);

        if(useLogX)
            pad2->SetLogx();

        pad2->Draw();
        pad2->cd();


        TH1D* hRatio =
            dynamic_cast<TH1D*>(
                hData->Clone(
                    ("hRatio_" + p.numHist).c_str()
                )
            );


        hRatio->SetDirectory(nullptr);
        hRatio->SetTitle("");
        hRatio->Divide(hTotal);

        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerSize(0.8);
        hRatio->SetLineColor(kBlack);

        hRatio->GetYaxis()->SetTitle("Data / MC");
        hRatio->GetYaxis()->SetNdivisions(303);
        hRatio->GetYaxis()->SetTitleSize(0.10);
        hRatio->GetYaxis()->SetLabelSize(0.08);
        hRatio->GetYaxis()->SetTitleOffset(0.45);

        hRatio->GetXaxis()->SetTitle(p.xTitle.c_str());
        hRatio->GetXaxis()->SetTitleSize(0.12);
        hRatio->GetXaxis()->SetLabelSize(0.10);
        hRatio->GetXaxis()->SetTitleOffset(1.0);

        hRatio->SetMinimum(0.75);
        hRatio->SetMaximum(1.25);
        hRatio->Draw("E1");


        TLine* ratioLine =
            new TLine(
                hRatio->GetXaxis()->GetXmin(),
                1.0,
                hRatio->GetXaxis()->GetXmax(),
                1.0
            );


        ratioLine->SetLineStyle(2);
        ratioLine->SetLineWidth(1);
        ratioLine->Draw("SAME");

        pad2->Update();
        c->Update();


        // ====================================================
        // Save
        // ====================================================

        c->SaveAs(
            p.output.c_str()
        );


        // ====================================================
        // Cleanup
        // ====================================================
    }


    // ========================================================
    // Close files
    // ========================================================

    fData->Close();
    fQCD->Close();
}
