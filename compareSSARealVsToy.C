#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TStyle.h>
#include <TLatex.h>

#include <vector>
#include <map>
#include <set>
#include <iostream>

//void compareSSARealVsToy(const char* realFile = "/Users/amandeep/Desktop/OuterTracker/full_results/10k_events/1sigma/Results.root",
void compareSSARealVsToy(const char* realFile = "/Users/amandeep/Desktop/OuterTracker/full_results/10k_events/1sigma/Results_2sigma.root",
                          const char* treeName = "noiseTree_Hyb1_SSA",
                          //const char* toyFile  = "ToyHits_Hyb0_2sigma.root", 
                          //const char* toyFile  = "ToyHits_Hyb0_2sigma_flat_f_corr.root", 
                          const char* toyFile  = "ToyHits_Hyb1_2sigma_subChip_f_corr.root", 
                          const char* outDir = "SSA_comparison")
{
    gStyle->SetOptStat(0);

    // Open the ROOT file with real data
    TFile* file = TFile::Open(realFile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open " << realFile << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "Error: cannot find tree " << treeName << std::endl;
        file->Close();
        return;
    }

        // --- Open toy model file ---
        TFile* fToy = TFile::Open(toyFile, "READ");
        if (!fToy || fToy->IsZombie()) {
            std::cerr << "Error: cannot open toy file " << toyFile << std::endl;
            return;
        }

    int eventID;
    std::vector<short>* hitRows = nullptr;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("hitRows", &hitRows);

    const int nChips = 8;
    const int nRowsPerChip = 120;

    // Loop over chips
    for (int chip = 0; chip < nChips; ++chip) {
        int rowStart = chip * nRowsPerChip;
        int rowEnd   = rowStart + nRowsPerChip - 1;

        // Map: eventID -> total hits in this chip
        std::map<int,int> hitsPerEvent;

        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; ++i) {
            tree->GetEntry(i);
            if (!hitRows) continue;

            int countThisEntry = 0;
            for (size_t j = 0; j < hitRows->size(); ++j) {
                short row = hitRows->at(j);
                if (row >= rowStart && row <= rowEnd) {
                    countThisEntry++;
                }
            }

            // Accumulate hits per event
            hitsPerEvent[eventID] += countThisEntry;
        }

        // Make histogram for this chip
        TString histName = TString::Format("hHits_chip%d", chip);
        TH1F* hHits = new TH1F(histName, TString::Format("SSA Chip %d Hits;Hits per Event;Events", chip),
                               120, 0, 120);

        for (const auto& evt : hitsPerEvent) {
            hHits->Fill(evt.second);
        }

        // --- Load toy histogram from file ---
        TString hNameToy = TString::Format("hToy_chip%d", chip);
        TH1F* hToy = (TH1F*)fToy->Get(hNameToy);
        if (!hToy) {
            std::cerr << "Warning: toy histogram not found for chip " << chip << std::endl;
            continue;
        }

        // Normalize both histograms to unity for shape comparison
        hHits->Scale(1.0 / hHits->Integral());
        hToy->Scale(1.0 / hToy->Integral());

        hHits->SetLineColor(kBlue);
        hHits->SetLineWidth(2);
        hToy->SetLineColor(kRed);
        hToy->SetLineStyle(2);
        hToy->SetLineWidth(2);

        // --- Draw overlay ---
        TCanvas* c1 = new TCanvas(TString::Format("c_chip%d", chip), TString::Format("SSA Chip %d", chip), 800, 600);
        c1->SetLogy();

        hHits->Draw("HIST");
        hToy->Draw("HIST SAME");

        // Legend
        TLegend* leg = new TLegend(0.55, 0.75, 0.85, 0.9);
        leg->AddEntry(hHits, "Real Data", "l");
        leg->AddEntry(hToy, "Toy Model", "l");
        leg->Draw();

        // Stats text
        double meanReal = hHits->GetMean();
        double sigmaReal = hHits->GetStdDev();
        double meanToy = hToy->GetMean();
        double sigmaToy = hToy->GetStdDev();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.15, 0.85, TString::Format("Real:  #mu=%.2f, #sigma=%.2f", meanReal, sigmaReal));
        latex.DrawLatex(0.15, 0.80, TString::Format("Toy:   #mu=%.2f, #sigma=%.2f", meanToy, sigmaToy));

        //TString pdfName = TString::Format("%s/Overlay_chip%d_Hyb0_2sigma.pdf", outDir, chip);
        //TString pdfName = TString::Format("%s/Overlay_chip%d_Hyb0_2sigma_flat_f_corr.pdf", outDir, chip);
        TString pdfName = TString::Format("%s/Overlay_chip%d_Hyb1_2sigma_subChip_f_corr.pdf", outDir, chip);
        c1->SaveAs(pdfName);

        delete c1;
        delete hHits;
        // hToy belongs to toy file, no delete
    }

    file->Close();
    fToy->Close();
    std::cout << "Overlay plots saved in " << outDir << std::endl;
//        double mean = hHits->GetMean();
//        double sigma = hHits->GetStdDev();
//        std::cout << "Chip " << chip << ": mean hits = " << mean << ", sigma = " << sigma << std::endl;
//
//        // Draw and save
//        TCanvas* c1 = new TCanvas(TString::Format("c_chip%d", chip), TString::Format("SSA Chip %d", chip), 800, 600);
//        c1->SetLogy();
//        hHits->Draw();
//
//        TLatex latex;
//        latex.SetNDC();
//        latex.SetTextSize(0.04);
//        latex.DrawLatex(0.15, 0.85, TString::Format("Mean = %.2f, #sigma = %.2f", mean, sigma));
//
//        TString pdfName = TString::Format("SSA_chip_%d_hits.pdf", chip);
//        //TString pdfName = TString::Format("%s/SSA_chip_%d_hits.pdf", outDir, chip);
//        c1->SaveAs(pdfName);
//
//        delete c1;
//        delete hHits;
//    }
//
//    file->Close();
//    std::cout << "Done! Histograms saved in " << outDir << std::endl;
}

