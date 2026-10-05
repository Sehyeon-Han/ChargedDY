#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <iostream>

void Make_QCD_BDT_Hist()
{
    TFile *fin = TFile::Open("Data_QCD_enriched.root");

    if(!fin || fin->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open Data_QCD_enriched.root" << std::endl;
        return;
    }

    TTree *tree = (TTree*)fin->Get("BDTTree");

    if(!tree)
    {
        std::cerr << "[ERROR] Cannot find BDTTree" << std::endl;
        fin->Close();
        return;
    }

    TH1::SetDefaultSumw2();

    // -----------------------------------------
    // BDT score cut
    // -----------------------------------------
    const char* bdtCut = "bdt_score < -0.20 && mt >= 150";

    // -----------------------------------------
    // mT binning
    // -----------------------------------------
    double mTBins[] = {150, 200, 250, 300, 350, 425, 500, 600, 750, 900, 1100, 1400, 2000, 5000};
    const int N_mTBins = sizeof(mTBins) / sizeof(mTBins[0]) - 1;

    // -----------------------------------------
    // Histograms
    // -----------------------------------------

    TH1D *h_mT = new TH1D("h_mT", "QCD;m_{T} [GeV];Events", N_mTBins, mTBins);
    TH1D *h_Muon_pt = new TH1D( "h_Muon_pt", "QCD;Muon p_{T} [GeV];Events", 30, 65, 1500);
    TH1D *h_Muon_eta = new TH1D("h_Muon_eta", "QCD;Muon #eta;Events", 50, -3, 3);
    TH1D *h_Muon_phi = new TH1D("h_Muon_phi", "QCD;Muon #phi;Events", 50, -3.3, 3.3);
    TH1D *h_MET_pt = new TH1D("h_MET_pt", "QCD;MET p_{T} [GeV];Events", 30, 85, 3500);
    TH1D *h_MET_phi = new TH1D("h_MET_phi", "QCD;MET #phi;Events", 50, -3.3, 3.3);

    // -----------------------------------------
    // Fill histograms with BDT cut
    // -----------------------------------------

    tree->Draw("mt>>h_mT", bdtCut, "goff");
    tree->Draw("mu_pt>>h_Muon_pt", bdtCut, "goff");
    tree->Draw("mu_eta>>h_Muon_eta", bdtCut, "goff");
    tree->Draw("mu_phi>>h_Muon_phi", bdtCut, "goff");
    tree->Draw("met>>h_MET_pt", bdtCut, "goff");
    tree->Draw("MET_phi>>h_MET_phi", bdtCut, "goff");

    // -----------------------------------------
    // Print information
    // -----------------------------------------

    Long64_t nTotal = tree->GetEntries();
    Long64_t nPass  = tree->GetEntries(bdtCut);

    std::cout << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << " QCD BDT selection" << std::endl;
    std::cout << "======================================" << std::endl;

    std::cout << "BDT cut       : "
              << bdtCut << std::endl;

    std::cout << "Total entries : "
              << nTotal << std::endl;

    std::cout << "Passed entries: "
              << nPass << std::endl;

    if(nTotal > 0)
    {
        std::cout << "Efficiency    : "
                  << static_cast<double>(nPass) /
                     static_cast<double>(nTotal)
                  << std::endl;
    }

    std::cout << std::endl;

    std::cout << "h_mT integral       = "
              << h_mT->Integral() << std::endl;

    std::cout << "h_Muon_pt integral  = "
              << h_Muon_pt->Integral() << std::endl;

    std::cout << "h_Muon_eta integral = "
              << h_Muon_eta->Integral() << std::endl;

    std::cout << "h_Muon_phi integral = "
              << h_Muon_phi->Integral() << std::endl;

    std::cout << "h_MET_pt integral   = "
              << h_MET_pt->Integral() << std::endl;

    std::cout << "h_MET_phi integral  = "
              << h_MET_phi->Integral() << std::endl;

    // -----------------------------------------
    // Save ROOT file
    // -----------------------------------------

    TFile *fout = new TFile("Muon_QCD_BDT_2016_preVFP.root", "RECREATE");

    h_mT->Write();
    h_Muon_pt->Write();
    h_Muon_eta->Write();
    h_Muon_phi->Write();
    h_MET_pt->Write();
    h_MET_phi->Write();

    fout->Close();
    fin->Close();

    std::cout << std::endl;
    std::cout
        << "Saved: Muon_QCD_BDT_2016_preVFP.root"
        << std::endl;
}