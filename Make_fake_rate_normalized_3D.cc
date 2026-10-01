#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH3D.h>
#include <TString.h>
#include <cmath>
#include <string>

struct SampleInfo
{
    std::string fileName;
    double xsec;
    double sumW;
};

void Make_fake_rate_normalized_3D(const char* dataFile = "fake_mu_merged_3D.root",
                                  const char* outFile = "fake_rate_3D.root")
{
    const double xsec200 = 8.78288653e00;
    const double xsec500 = 2.76542446e-01;
    const double xsec1000 = 1.58702653e-02;
    const double xsec2000 = 4.15066636e-04;
    const double DYxsec200 = 2.78;
    const double DYxsec400 = 0.15;
    const double DYxsec500 = 0.084;
    const double DYxsec700 = 0.013;
    const double DYxsec800 = 0.011;
    const double DYxsec1000 = 0.006;
    const double DYxsec1500 = 0.00081;
    const double DYxsec2000 = 0.0002;
    const double TTTo2L2Nuxsec = 89.31;
    const double TTToSemiLeptonicxsec = 367.8;
    const double WWTo1L1Nu2Qxsec = 51.65;
    const double WWTo4Q_4fxsec = 51.03;
    const double WWTo2L2Nuxsec = 11.09;
    const double WZTo1L1Nu2Q_4fxsec = 9.119;
    const double WZTo3LNuxsec = 5.213;
    const double WZTo2Q2Nu_4fxsec = 6.331;
    const double WZTo1L3Nu_4fxsec = 3.414;
    const double ZZxsec = 12.17;
    const double ST_s_channel_4f_leptonDecays_xsec = 10.32;
    const double ST_t_channel_top_4f_InclusiveDecays_xsec = 134.2;
    const double ST_t_channel_anitop_4f_InclusiveDecays_xsec = 80.0;
    const double ST_tW_top_5f_inclusiveDecays_xsec = 39.65;
    const double ST_tW_antitop_5f_inclusiveDecays_xsec = 39.65;
    const double WToTauNu_M_200_xsec = 8.78288653e00;
    const double WToTauNu_M_500_xsec = 2.76542446e-01;
    const double WToTauNu_M_1000_xsec = 1.58702653e-02;
    const double WToTauNu_M_2000_xsec = 4.15066636e-04;
    const double lumi = 19520.0;
    const double sumW200 = 3250000;
    const double sumW500 = 750000;
    const double sumW1000 = 250000;
    const double sumW2000 = 50000;
    const double sumWDY200 = 2673760.516799;
    const double sumWDY400 = 44344.277472;
    const double sumWDY500 = 24718.574857;
    const double sumWDY700 = 3951.644255;
    const double sumWDY800 = 3345.560043;
    const double sumWDY1000 = 1782.846123;
    const double sumWDY1500 = 244.385548;
    const double sumWDY2000 = 62.650483;
    const double sumWTTTo2L2Nu = 2.70453e+09;
    const double sumWTTToSemiLeptonic = 3.97723e+10;
    const double sumWWWTo1L1Nu2Q = 1.69101e+09;
    const double sumWWWTo2L2Nu = 3.34565e+07;
    const double sumWWWTo4Q_4f = 1.67235e+09;
    const double sumWWZTo1L1Nu2Q_4f = 5.58378e+07;
    const double sumWWZTo1L3Nu_4f = 6.96859e+06;
    const double sumWWZTo2Q2Nu_4f = 4.84591e+07;
    const double sumWWZTo3LNu = 8.15165e+07;
    const double sumWZZ = 1282000;
    const double sumW_ST_s_channel_4f_leptonDecay = 1.95962e+07;
    const double sumW_ST_t_channel_anitop_4f_InclusiveDecays = 1.98386e+09;
    const double sumW_ST_t_channel_top_4f_InclusiveDecays = 5.94814e+09;
    const double sumW_ST_tW_top_5f_inclusiveDecays = 7.46247e+07;
    const double sumW_ST_tW_anitop_5f_inclusiveDecays = 7.47663e+07;
    const double sumW_WToTauNu_M_200 = 1250000;
    const double sumW_WToTauNu_M_500 = 750000;
    const double sumW_WToTauNu_M_1000 = 250000;
    const double sumW_WToTauNu_M_2000 = 50000;

    std::vector<SampleInfo> samples = 
    {
        {"fake_mu_merged_WToMuNu_M200_3D.root", xsec200, sumW200},
        {"fake_mu_merged_WToMuNu_M500_3D.root", xsec500, sumW500},
        {"fake_mu_merged_WToMuNu_M1000_3D.root", xsec1000, sumW1000},
        {"fake_mu_merged_WToMuNu_M2000_3D.root", xsec2000, sumW2000},
        {"fake_mu_merged_WToTauNu_M200_3D.root", xsec200, sumW_WToTauNu_M_200},
        {"fake_mu_merged_WToTauNu_M500_3D.root", xsec500, sumW_WToTauNu_M_500},
        {"fake_mu_merged_WToTauNu_M1000_3D.root", xsec1000,sumW_WToTauNu_M_1000},
        {"fake_mu_merged_WToTauNu_M2000_3D.root", xsec2000, sumW_WToTauNu_M_2000},
        {"fake_mu_merged_DYMC_200_3D.root", DYxsec200, sumWDY200},
        {"fake_mu_merged_DYMC_400_3D.root", DYxsec400, sumWDY400},
        {"fake_mu_merged_DYMC_500_3D.root", DYxsec500, sumWDY500},
        {"fake_mu_merged_DYMC_700_3D.root", DYxsec700, sumWDY700},
        {"fake_mu_merged_DYMC_800_3D.root", DYxsec800, sumWDY800},
        {"fake_mu_merged_DYMC_1000_3D.root", DYxsec1000, sumWDY1000},
        {"fake_mu_merged_DYMC_1500_3D.root", DYxsec1500, sumWDY1500},
        {"fake_mu_merged_DYMC_2000_3D.root", DYxsec2000, sumWDY2000},
        {"fake_mu_merged_TTTo2L2Nu_3D.root", TTTo2L2Nuxsec,sumWTTTo2L2Nu},
        {"fake_mu_merged_TTToSemiLeptonic_3D.root", TTToSemiLeptonicxsec, sumWTTToSemiLeptonic},
        {"fake_mu_merged_ST_s_channel_3D.root", ST_s_channel_4f_leptonDecays_xsec, sumW_ST_s_channel_4f_leptonDecay},
        {"fake_mu_merged_ST_t_top_3D.root", ST_t_channel_top_4f_InclusiveDecays_xsec, sumW_ST_t_channel_top_4f_InclusiveDecays},
        {"fake_mu_merged_ST_t_antitop_3D.root", ST_t_channel_anitop_4f_InclusiveDecays_xsec, sumW_ST_t_channel_anitop_4f_InclusiveDecays},
        {"fake_mu_merged_ST_tW_top_3D.root", ST_tW_top_5f_inclusiveDecays_xsec, sumW_ST_tW_top_5f_inclusiveDecays},
        {"fake_mu_merged_ST_tW_antitop_3D.root", ST_tW_antitop_5f_inclusiveDecays_xsec, sumW_ST_tW_anitop_5f_inclusiveDecays},
        {"fake_mu_merged_WWTo1L1Nu2Q_3D.root", WWTo1L1Nu2Qxsec, sumWWWTo1L1Nu2Q},
        {"fake_mu_merged_WWTo2L2Nu_3D.root", WWTo2L2Nuxsec, sumWWWTo2L2Nu},
        {"fake_mu_merged_WWTo4Q_3D.root", WWTo4Q_4fxsec, sumWWWTo4Q_4f},
        {"fake_mu_merged_WZTo1L1Nu2Q_4f_3D.root", WZTo1L1Nu2Q_4fxsec, sumWWZTo1L1Nu2Q_4f},
        {"fake_mu_merged_WZTo1L3Nu_4f_3D.root", WZTo1L3Nu_4fxsec, sumWWZTo1L3Nu_4f},
        {"fake_mu_merged_WZTo2Q2Nu_4f_3D.root", WZTo2Q2Nu_4fxsec, sumWWZTo2Q2Nu_4f},
        {"fake_mu_merged_WZTo3LNu_3D.root", WZTo3LNuxsec, sumWWZTo3LNu},
        {"fake_mu_merged_ZZ_3D.root", ZZxsec, sumWZZ}
    };

    TFile* fData = TFile::Open(dataFile, "READ");

    if(!fData || fData->IsZombie())
    {
        std::cerr << "Cannot open DATA file: " << dataFile << std::endl;

        return;
    }

    TH3D* h_data_den_in = dynamic_cast<TH3D*>(fData->Get("h_fake_den_3D"));
    TH3D* h_data_num_in = dynamic_cast<TH3D*>(fData->Get("h_fake_num_3D"));

    TH3D* h_data_den = dynamic_cast<TH3D*>(h_data_den_in->Clone("h_fake_den_3D_data"));
    TH3D* h_data_num = dynamic_cast<TH3D*>(h_data_num_in->Clone("h_fake_num_3D_data"));

    if(!h_data_den_in || !h_data_num_in)
    {
        std::cerr << "[ERROR] Missing DATA histograms." << std::endl;
        fData->Close();
        return;
    }

    h_data_den->SetDirectory(nullptr);
    h_data_num->SetDirectory(nullptr);

    fData->Close();

    TH3D* h_prompt_den = dynamic_cast<TH3D*>(h_data_den->Clone("h_fake_den_3D_prompt"));
    TH3D* h_prompt_num = dynamic_cast<TH3D*>(h_data_num->Clone("h_fake_num_3D_prompt"));

    h_prompt_den->Reset("ICES");
    h_prompt_num->Reset("ICES");

    h_prompt_den->Sumw2();
    h_prompt_num->Sumw2();

    for(const auto& sample : samples)
    {
        TFile* fMC = TFile::Open(sample.fileName.c_str(), "READ");

        if(!fMC || fMC->IsZombie())
        {
            std::cerr << "[WARNING] Cannot open: " << sample.fileName << std::endl;

            continue;
        }

        TH3D* h_den_in = dynamic_cast<TH3D*>(fMC->Get("h_fake_den_3D"));
        TH3D* h_num_in = dynamic_cast<TH3D*>(fMC->Get("h_fake_num_3D"));

        if(!h_den_in || !h_num_in)
        {
            std::cerr << "[WARNING] Missing histogram in " << sample.fileName << std::endl;

            fMC->Close();
            continue;
        }

        TH3D* h_den = dynamic_cast<TH3D*>(h_den_in->Clone());
        TH3D* h_num = dynamic_cast<TH3D*>(h_num_in->Clone());

        h_den->SetDirectory(nullptr);
        h_num->SetDirectory(nullptr);

        fMC->Close();

        const double norm = lumi * sample.xsec / sample.sumW;

        h_den->Scale(norm);
        h_num->Scale(norm);

        h_prompt_den->Add(h_den);
        h_prompt_num->Add(h_num);

        std::cout << sample.fileName << " norm = " << norm << " den = " << h_den->Integral() << " num = " << h_num->Integral() << std::endl;

        delete h_den;
        delete h_num;
    }

    TH3D* h_fake_den = dynamic_cast<TH3D*>(h_data_den->Clone("h_fake_den_3D_subtracted"));
    TH3D* h_fake_num = dynamic_cast<TH3D*>(h_data_num->Clone("h_fake_num_3D_subtracted"));

    h_fake_den->Add(h_prompt_den, -1.0);
    h_fake_num->Add(h_prompt_num, -1.0);

    TH3D* h_fake_rate = dynamic_cast<TH3D*>(h_fake_num->Clone("h_fake_rate_3D"));

    h_fake_rate->Reset("ICES");

       h_fake_rate->SetTitle("Prompt-subtracted fake rate;" "Muon p_{T} [GeV];" "|#eta|;" "MET [GeV]");

    for(int ix = 1; ix <= h_fake_den->GetNbinsX(); ++ix)
    {
        for(int iy = 1; iy <= h_fake_den->GetNbinsY(); ++iy)
        {
            for(int iz = 1; iz <= h_fake_den->GetNbinsZ(); ++iz)
            {
                const double den = h_fake_den->GetBinContent(ix, iy, iz);
                const double num = h_fake_num->GetBinContent(ix, iy, iz);
                    
                if(den <= 0.0)
                {
                    std::cerr << "[WARNING] non-positive denominator" << " pT bin = " << ix << " eta bin = " << iy << " MET bin = " << iz << " den = " << den << " num = " << num << std::endl;

                    h_fake_rate->SetBinContent(ix, iy, iz, 0.0);
                    h_fake_rate->SetBinError(ix, iy, iz, 0.0);

                    continue;
                }

                const double rate = num / den;

                h_fake_rate->SetBinContent(ix, iy, iz, rate);

                double numErr = h_fake_num->GetBinError(ix, iy, iz);
                double denErr = h_fake_den->GetBinError(ix, iy, iz);

                double rateErr = 0.0;

                if(num != 0.0)
                {
                    rateErr = std::fabs(rate) * std::sqrt((numErr * numErr) / (num * num) + (denErr * denErr) / (den * den));
                }

                h_fake_rate->SetBinError(ix, iy, iz, rateErr);

                if(rate < 0.0 || rate > 1.0)
                {
                    std::cerr << "[WARNING] unusual fake rate" << " pT bin = " << ix << "eta bin = " << iy << " MET bin = " << iz << " rate = " << rate << " data den = " << h_data_den->GetBinContent(ix, iy, iz) << " prompt dem = " << h_prompt_den->GetBinContent(ix, iy, iz) << std::endl;
                }
            }
        }
    }

    std::cout << "\n===== FINAL FAKE RATE =====\n";

    for(int ix = 1; ix <= h_fake_rate->GetNbinsX(); ++ix)
    {
        for(int iy = 1; iy <= h_fake_rate->GetNbinsY(); ++iy)
        {
            for(int iz = 1; iz <= h_fake_rate->GetNbinsZ(); ++iz)
            {
                std::cout
                    << "pT=" << ix
                    << " eta=" << iy
                    << " MET=" << iz
                    << " | data den="
                    << h_data_den->GetBinContent(ix, iy, iz)
                    << " prompt den="
                    << h_prompt_den->GetBinContent(ix, iy, iz)
                    << " fake den="
                    << h_fake_den->GetBinContent(ix, iy, iz)
                    << " | data num="
                    << h_data_num->GetBinContent(ix, iy, iz)
                    << " prompt num="
                    << h_prompt_num->GetBinContent(ix, iy, iz)
                    << " fake num="
                    << h_fake_num->GetBinContent(ix, iy, iz)
                    << " | rate="
                    << h_fake_rate->GetBinContent(ix, iy, iz)
                    << " +/- "
                    << h_fake_rate->GetBinError(ix, iy, iz)
                    << std::endl;
            }
        }
    }

    TFile* fout = TFile::Open(outFile, "RECREATE");

    h_data_den->Write();
    h_data_num->Write();
    h_prompt_den->Write();
    h_prompt_num->Write();
    h_fake_den->Write();
    h_fake_num->Write();
    h_fake_rate->Write();
    fout->Close();
}