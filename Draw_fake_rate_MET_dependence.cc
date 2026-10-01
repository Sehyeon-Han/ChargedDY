#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>


// ================================================================
// Sample information
// ================================================================
struct SampleInfo
{
    TString fileName;
    TString sampleName;
    double xsec;
    double sumW;
};


// ================================================================
// MC normalization
//
// scale = xsec * luminosity / sumW
// ================================================================
double GetScale(const SampleInfo& sample, double lumi)
{
    if(sample.sumW == 0.0)
    {
        std::cerr
            << "[ERROR] sumW = 0 : "
            << sample.sampleName
            << std::endl;

        return 0.0;
    }

    return sample.xsec * lumi / sample.sumW;
}


// ================================================================
// Load TH1D
// ================================================================
TH1D* LoadTH1D(const TString& fileName,
               const TString& histName,
               const TString& cloneName)
{
    TFile* file = TFile::Open(fileName, "READ");

    if(!file || file->IsZombie())
    {
        std::cerr
            << "[ERROR] Cannot open file: "
            << fileName
            << std::endl;

        if(file)
            file->Close();

        return nullptr;
    }

    TH1D* source =
        dynamic_cast<TH1D*>(file->Get(histName));

    if(!source)
    {
        std::cerr
            << "[ERROR] Cannot find histogram "
            << histName
            << " in "
            << fileName
            << std::endl;

        file->Close();
        delete file;

        return nullptr;
    }

    TH1D* hist =
        dynamic_cast<TH1D*>(
            source->Clone(cloneName)
        );

    hist->SetDirectory(nullptr);

    file->Close();
    delete file;

    return hist;
}


// ================================================================
// Check 1D histogram binning
// ================================================================
bool SameBinning1D(const TH1D* h1,
                   const TH1D* h2)
{
    if(!h1 || !h2)
        return false;

    if(h1->GetNbinsX() != h2->GetNbinsX())
        return false;

    for(int i = 1; i <= h1->GetNbinsX(); ++i)
    {
        double edge1 =
            h1->GetXaxis()->GetBinLowEdge(i);

        double edge2 =
            h2->GetXaxis()->GetBinLowEdge(i);

        if(std::fabs(edge1 - edge2) > 1.0e-9)
            return false;
    }

    double last1 =
        h1->GetXaxis()->GetBinUpEdge(
            h1->GetNbinsX()
        );

    double last2 =
        h2->GetXaxis()->GetBinUpEdge(
            h2->GetNbinsX()
        );

    return std::fabs(last1 - last2) < 1.0e-9;
}


// ================================================================
// Make MET-dependent fake rates for all 20 (pT, |eta|) bins
//
// fake rate =
// (Data numerator - Prompt numerator)
// -----------------------------------
// (Data denominator - Prompt denominator)
//
// Input histogram names:
//   h_fake_den_MET_pt1_eta1 ... h_fake_den_MET_pt5_eta4
//   h_fake_num_MET_pt1_eta1 ... h_fake_num_MET_pt5_eta4
//
// MET bins:
//   0 - 20 GeV
//   20 - 40 GeV
//   40 - 65 GeV
// ================================================================
void MakeFakeRateMET(const TString& dataFile,
                     const std::vector<SampleInfo>& samples,
                     double lumi,
                     TFile* outputFile)
{
    const int N_ptBins = 5;
    const int N_etaBins = 4;

    double ptBins[]  = {65, 100, 150, 200, 300, 1500};
    double etaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4};

    for(int ipt = 0; ipt < N_ptBins; ++ipt)
    {
        for(int ieta = 0; ieta < N_etaBins; ++ieta)
        {
            TString denHistName =
                Form(
                    "h_fake_den_MET_pt%d_eta%d",
                    ipt + 1,
                    ieta + 1
                );

            TString numHistName =
                Form(
                    "h_fake_num_MET_pt%d_eta%d",
                    ipt + 1,
                    ieta + 1
                );


            // ====================================================
            // Load Data
            // ====================================================
            TH1D* hDataDen =
                LoadTH1D(
                    dataFile,
                    denHistName,
                    Form(
                        "h_data_den_MET_pt%d_eta%d",
                        ipt + 1,
                        ieta + 1
                    )
                );

            TH1D* hDataNum =
                LoadTH1D(
                    dataFile,
                    numHistName,
                    Form(
                        "h_data_num_MET_pt%d_eta%d",
                        ipt + 1,
                        ieta + 1
                    )
                );

            if(!hDataDen || !hDataNum)
            {
                std::cerr
                    << "[ERROR] Cannot load Data histograms for "
                    << "pt" << ipt + 1
                    << " eta" << ieta + 1
                    << std::endl;

                delete hDataDen;
                delete hDataNum;

                continue;
            }


            // ====================================================
            // Prompt MC total histograms
            // ====================================================
            TH1D* hPromptDen =
                dynamic_cast<TH1D*>(
                    hDataDen->Clone(
                        Form(
                            "h_prompt_den_MET_pt%d_eta%d",
                            ipt + 1,
                            ieta + 1
                        )
                    )
                );

            TH1D* hPromptNum =
                dynamic_cast<TH1D*>(
                    hDataNum->Clone(
                        Form(
                            "h_prompt_num_MET_pt%d_eta%d",
                            ipt + 1,
                            ieta + 1
                        )
                    )
                );

            hPromptDen->Reset("ICES");
            hPromptNum->Reset("ICES");

            hPromptDen->SetDirectory(nullptr);
            hPromptNum->SetDirectory(nullptr);


            // ====================================================
            // Loop over prompt MC samples
            // ====================================================
            for(const SampleInfo& sample : samples)
            {
                TH1D* hDen =
                    LoadTH1D(
                        sample.fileName,
                        denHistName,
                        Form(
                            "h_den_MET_%s_pt%d_eta%d",
                            sample.sampleName.Data(),
                            ipt + 1,
                            ieta + 1
                        )
                    );

                TH1D* hNum =
                    LoadTH1D(
                        sample.fileName,
                        numHistName,
                        Form(
                            "h_num_MET_%s_pt%d_eta%d",
                            sample.sampleName.Data(),
                            ipt + 1,
                            ieta + 1
                        )
                    );

                if(!hDen || !hNum)
                {
                    delete hDen;
                    delete hNum;
                    continue;
                }


                // -----------------------------------------------
                // Check binning
                // -----------------------------------------------
                if(!SameBinning1D(hDataDen, hDen) ||
                   !SameBinning1D(hDataNum, hNum))
                {
                    std::cerr
                        << "[ERROR] MET binning mismatch: "
                        << sample.fileName
                        << " | pt" << ipt + 1
                        << " eta" << ieta + 1
                        << std::endl;

                    delete hDen;
                    delete hNum;

                    continue;
                }


                // -----------------------------------------------
                // MC normalization
                // scale = xsec * lumi / sumW
                // -----------------------------------------------
                double scale =
                    GetScale(sample, lumi);

                if(scale == 0.0)
                {
                    delete hDen;
                    delete hNum;
                    continue;
                }


                double rawDen =
                    hDen->Integral(
                        0,
                        hDen->GetNbinsX() + 1
                    );

                double rawNum =
                    hNum->Integral(
                        0,
                        hNum->GetNbinsX() + 1
                    );


                hDen->Scale(scale);
                hNum->Scale(scale);


                // -----------------------------------------------
                // Add to total prompt MC
                // -----------------------------------------------
                hPromptDen->Add(hDen);
                hPromptNum->Add(hNum);


                std::cout
                    << "[MC] "
                    << sample.sampleName
                    << " | pt" << ipt + 1
                    << " eta" << ieta + 1
                    << " | raw den = "
                    << rawDen
                    << " | raw num = "
                    << rawNum
                    << " | scale = "
                    << scale
                    << " | norm den = "
                    << hDen->Integral(
                           0,
                           hDen->GetNbinsX() + 1
                       )
                    << " | norm num = "
                    << hNum->Integral(
                           0,
                           hNum->GetNbinsX() + 1
                       )
                    << std::endl;


                delete hDen;
                delete hNum;
            }


            // ====================================================
            // QCD = Data - Prompt MC
            // ====================================================
            TH1D* hFakeDen =
                dynamic_cast<TH1D*>(
                    hDataDen->Clone(
                        Form(
                            "h_fake_den_MET_subtracted_pt%d_eta%d",
                            ipt + 1,
                            ieta + 1
                        )
                    )
                );

            TH1D* hFakeNum =
                dynamic_cast<TH1D*>(
                    hDataNum->Clone(
                        Form(
                            "h_fake_num_MET_subtracted_pt%d_eta%d",
                            ipt + 1,
                            ieta + 1
                        )
                    )
                );

            hFakeDen->SetDirectory(nullptr);
            hFakeNum->SetDirectory(nullptr);

            hFakeDen->Add(hPromptDen, -1.0);
            hFakeNum->Add(hPromptNum, -1.0);


            // ====================================================
            // Fake-rate histogram
            // ====================================================
            TH1D* hFakeRate =
                dynamic_cast<TH1D*>(
                    hFakeNum->Clone(
                        Form(
                            "h_fake_rate_MET_pt%d_eta%d",
                            ipt + 1,
                            ieta + 1
                        )
                    )
                );

            hFakeRate->Reset("ICES");
            hFakeRate->SetDirectory(nullptr);


            TString ptLabel;

            if(ipt == N_ptBins - 1)
            {
                ptLabel =
                    Form(
                        "p_{T} #geq %.0f GeV",
                        ptBins[ipt]
                    );
            }
            else
            {
                ptLabel =
                    Form(
                        "%.0f #leq p_{T} < %.0f GeV",
                        ptBins[ipt],
                        ptBins[ipt + 1]
                    );
            }

            TString etaLabel =
                Form(
                    "%.1f #leq |#eta| < %.1f",
                    etaBins[ieta],
                    etaBins[ieta + 1]
                );

            hFakeRate->SetTitle(
                Form(
                    "%s, %s;MET [GeV];Fake rate",
                    ptLabel.Data(),
                    etaLabel.Data()
                )
            );


            // ====================================================
            // Calculate fake rate
            // ====================================================
            std::cout << std::endl;
            std::cout
                << "=============================================="
                << std::endl;

            std::cout
                << " MET DEPENDENT FAKE RATE"
                << " | pt" << ipt + 1
                << " eta" << ieta + 1
                << std::endl;

            std::cout
                << " "
                << ptLabel
                << ", "
                << etaLabel
                << std::endl;

            std::cout
                << "=============================================="
                << std::endl;


            for(int bin = 1;
                bin <= hFakeRate->GetNbinsX();
                ++bin)
            {
                double low =
                    hFakeRate->GetXaxis()
                             ->GetBinLowEdge(bin);

                double high =
                    hFakeRate->GetXaxis()
                             ->GetBinUpEdge(bin);


                double dataDen =
                    hDataDen->GetBinContent(bin);

                double dataNum =
                    hDataNum->GetBinContent(bin);

                double promptDen =
                    hPromptDen->GetBinContent(bin);

                double promptNum =
                    hPromptNum->GetBinContent(bin);

                double denominator =
                    hFakeDen->GetBinContent(bin);

                double numerator =
                    hFakeNum->GetBinContent(bin);


                double denominatorError =
                    hFakeDen->GetBinError(bin);

                double numeratorError =
                    hFakeNum->GetBinError(bin);


                // -----------------------------------------------
                // Check denominator
                // -----------------------------------------------
                if(denominator <= 0.0)
                {
                    hFakeRate->SetBinContent(bin, 0.0);
                    hFakeRate->SetBinError(bin, 0.0);

                    std::cerr
                        << "[WARNING] denominator <= 0"
                        << " | pt" << ipt + 1
                        << " eta" << ieta + 1
                        << " | MET = "
                        << low
                        << " - "
                        << high
                        << " GeV"
                        << " | denominator = "
                        << denominator
                        << std::endl;

                    continue;
                }


                // -----------------------------------------------
                // Fake rate
                // -----------------------------------------------
                double rate =
                    numerator / denominator;


                // -----------------------------------------------
                // Error propagation
                //
                // Same convention as the original macro:
                // numerator and denominator uncertainties are
                // treated independently after prompt subtraction.
                // -----------------------------------------------
                double rateError = 0.0;

                if(numerator != 0.0)
                {
                    rateError =
                        std::fabs(rate) *
                        std::sqrt(
                            std::pow(
                                numeratorError / numerator,
                                2
                            )
                            +
                            std::pow(
                                denominatorError / denominator,
                                2
                            )
                        );
                }
                else
                {
                    rateError =
                        numeratorError / denominator;
                }


                hFakeRate->SetBinContent(
                    bin,
                    rate
                );

                hFakeRate->SetBinError(
                    bin,
                    rateError
                );


                // -----------------------------------------------
                // Print bin information
                // -----------------------------------------------
                std::cout << std::endl;

                std::cout
                    << "MET "
                    << low
                    << " - "
                    << high
                    << " GeV"
                    << std::endl;

                std::cout
                    << "  Data denominator   = "
                    << dataDen
                    << std::endl;

                std::cout
                    << "  Prompt denominator = "
                    << promptDen
                    << std::endl;

                std::cout
                    << "  Fake denominator   = "
                    << denominator
                    << std::endl;

                std::cout
                    << "  Data numerator     = "
                    << dataNum
                    << std::endl;

                std::cout
                    << "  Prompt numerator   = "
                    << promptNum
                    << std::endl;

                std::cout
                    << "  Fake numerator     = "
                    << numerator
                    << std::endl;

                std::cout
                    << "  Fake rate          = "
                    << rate
                    << " +/- "
                    << rateError
                    << std::endl;
            }


            std::cout << std::endl;
            std::cout
                << "=============================================="
                << std::endl;


            // ====================================================
            // Canvas
            // ====================================================
            TCanvas* cFakeRate =
                new TCanvas(
                    Form(
                        "c_fake_rate_MET_pt%d_eta%d",
                        ipt + 1,
                        ieta + 1
                    ),
                    "Fake rate vs MET",
                    900,
                    700
                );

            cFakeRate->SetTopMargin(0.08);
            cFakeRate->SetBottomMargin(0.13);
            cFakeRate->SetLeftMargin(0.13);
            cFakeRate->SetRightMargin(0.05);


            hFakeRate->SetMarkerStyle(20);
            hFakeRate->SetMarkerSize(1.2);
            hFakeRate->SetLineWidth(2);

            hFakeRate->SetMinimum(0.0);
            hFakeRate->SetMaximum(1.0);


            hFakeRate->GetXaxis()
                     ->SetTitle("MET [GeV]");

            hFakeRate->GetYaxis()
                     ->SetTitle("Fake rate");


            hFakeRate->GetXaxis()
                     ->SetTitleSize(0.055);

            hFakeRate->GetYaxis()
                     ->SetTitleSize(0.055);


            hFakeRate->GetXaxis()
                     ->SetLabelSize(0.045);

            hFakeRate->GetYaxis()
                     ->SetLabelSize(0.045);


            hFakeRate->GetXaxis()
                     ->SetTitleOffset(1.05);

            hFakeRate->GetYaxis()
                     ->SetTitleOffset(1.10);


            hFakeRate->Draw("E1");


            // ====================================================
            // Save PNG
            // ====================================================
            cFakeRate->SaveAs(
                Form(
                    "fake_mu_rate_MET_dependence_pt%d_eta%d.png",
                    ipt + 1,
                    ieta + 1
                )
            );


            // ====================================================
            // Save ROOT objects
            // ====================================================
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

                cFakeRate->Write();
            }


            // ====================================================
            // Cleanup
            // ====================================================
            delete cFakeRate;

            delete hDataDen;
            delete hDataNum;

            delete hPromptDen;
            delete hPromptNum;

            delete hFakeDen;
            delete hFakeNum;

            delete hFakeRate;
        }
    }
}


// ================================================================
// Main function
// ================================================================
void Draw_fake_rate_MET_dependence(
    const char* dataFile =
        "fake_mu_merged_dependence.root",

    const char* outputRootFile =
        "fake_mu_MET_dependence_results.root",

    double lumi = 19520.0
)
{
    TH1::SetDefaultSumw2();

    gStyle->SetOptStat(0);


    // ============================================================
    // Prompt MC samples
    // ============================================================
    std::vector<SampleInfo> samples =
    {
        // ========================================================
        // W -> mu nu
        // ========================================================
        {
            "fake_mu_merged_WToMuNu_M200_dependence.root",
            "WToMuNu_M200",
            8.78288653e00,
            3250000.0
        },

        {
            "fake_mu_merged_WToMuNu_M500_dependence.root",
            "WToMuNu_M500",
            2.76542446e-01,
            750000.0
        },

        {
            "fake_mu_merged_WToMuNu_M1000_dependence.root",
            "WToMuNu_M1000",
            1.58702653e-02,
            250000.0
        },

        {
            "fake_mu_merged_WToMuNu_M2000_dependence.root",
            "WToMuNu_M2000",
            4.15066636e-04,
            50000.0
        },


        // ========================================================
        // W -> tau nu
        // ========================================================
        {
            "fake_mu_merged_WToTauNu_M200_dependence.root",
            "WToTauNu_M200",
            8.78288653e00,
            1250000.0
        },

        {
            "fake_mu_merged_WToTauNu_M500_dependence.root",
            "WToTauNu_M500",
            2.76542446e-01,
            750000.0
        },

        {
            "fake_mu_merged_WToTauNu_M1000_dependence.root",
            "WToTauNu_M1000",
            1.58702653e-02,
            250000.0
        },

        {
            "fake_mu_merged_WToTauNu_M2000_dependence.root",
            "WToTauNu_M2000",
            4.15066636e-04,
            50000.0
        },


        // ========================================================
        // DY
        // ========================================================
        {
            "fake_mu_DYMC_200_merged_dependence.root",
            "DY200",
            2.78,
            2673760.516799
        },

        {
            "fake_mu_DYMC_400_merged_dependence.root",
            "DY400",
            0.15,
            44344.277472
        },

        {
            "fake_mu_DYMC_500_merged_dependence.root",
            "DY500",
            0.084,
            24718.574857
        },

        {
            "fake_mu_DYMC_700_merged_dependence.root",
            "DY700",
            0.013,
            3951.644255
        },

        {
            "fake_mu_DYMC_800_merged_dependence.root",
            "DY800",
            0.011,
            3345.560043
        },

        {
            "fake_mu_DYMC_1000_merged_dependence.root",
            "DY1000",
            0.006,
            1782.846123
        },

        {
            "fake_mu_DYMC_1500_merged_dependence.root",
            "DY1500",
            0.00081,
            244.385548
        },

        {
            "fake_mu_DYMC_2000_merged_dependence.root",
            "DY2000",
            0.0002,
            62.650483
        },


        // ========================================================
        // ttbar
        // ========================================================
        {
            "fake_mu_merged_TTTo2L2Nu_dependence.root",
            "TTTo2L2Nu",
            89.31,
            2.70453e09
        },

        {
            "fake_mu_merged_TTToSemiLeptonic_dependence.root",
            "TTToSemiLeptonic",
            367.8,
            3.97723e10
        },


        // ========================================================
        // Single top
        // ========================================================
        {
            "fake_mu_merged_ST_s_channel_dependence.root",
            "ST_s_channel",
            10.32,
            1.95962e07
        },

        {
            "fake_mu_merged_ST_t_antitop_dependence.root",
            "ST_t_channel_antitop",
            80.0,
            1.98386e09
        },

        {
            "fake_mu_merged_ST_t_top_dependence.root",
            "ST_t_channel_top",
            134.2,
            5.94814e09
        },

        {
            "fake_mu_merged_ST_tW_antitop_dependence.root",
            "ST_tW_antitop",
            39.65,
            7.47663e07
        },

        {
            "fake_mu_merged_ST_tW_top_dependence.root",
            "ST_tW_top",
            39.65,
            7.46247e07
        },


        // ========================================================
        // Diboson
        // ========================================================
        {
            "fake_mu_merged_WWTo1L1Nu2Q_dependence.root",
            "WWTo1L1Nu2Q",
            51.65,
            1.69101e09
        },

        {
            "fake_mu_merged_WWTo2L2Nu_dependence.root",
            "WWTo2L2Nu",
            11.09,
            3.34565e07
        },

        {
            "fake_mu_merged_WWTo4Q_dependence.root",
            "WWTo4Q",
            51.03,
            1.67235e09
        },

        {
            "fake_mu_merged_WZTo1L1Nu2Q_4f_dependence.root",
            "WZTo1L1Nu2Q",
            9.119,
            5.58378e07
        },

        {
            "fake_mu_merged_WZTo1L3Nu_4f_dependence.root",
            "WZTo1L3Nu",
            3.414,
            6.96859e06
        },

        {
            "fake_mu_merged_WZTo2Q2Nu_4f_dependence.root",
            "WZTo2Q2Nu",
            6.331,
            4.84591e07
        },

        {
            "fake_mu_merged_WZTo3LNu_dependence.root",
            "WZTo3LNu",
            5.213,
            8.15165e07
        },

        {
            "fake_mu_merged_ZZ_dependence.root",
            "ZZ",
            12.17,
            1282000.0
        }
    };


    // ============================================================
    // Output ROOT file
    // ============================================================
    TFile* outputFile =
        TFile::Open(
            outputRootFile,
            "RECREATE"
        );

    if(!outputFile ||
       outputFile->IsZombie())
    {
        std::cerr
            << "[ERROR] Cannot create output file: "
            << outputRootFile
            << std::endl;

        return;
    }


    // ============================================================
    // Make MET-dependent fake rate
    // ============================================================
    MakeFakeRateMET(
        dataFile,
        samples,
        lumi,
        outputFile
    );


    outputFile->Close();
    delete outputFile;


    std::cout << std::endl;

    std::cout
        << "=============================================="
        << std::endl;

    std::cout
        << "MET fake-rate dependence completed."
        << std::endl;

    std::cout
        << "Output ROOT file : "
        << outputRootFile
        << std::endl;

    std::cout
        << "Output PNGs      : "
        << "fake_mu_rate_MET_dependence_pt1_eta1.png"
        << " ... "
        << "fake_mu_rate_MET_dependence_pt5_eta4.png"
        << std::endl;

    std::cout
        << "=============================================="
        << std::endl;
}
