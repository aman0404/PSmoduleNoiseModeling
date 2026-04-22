#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

void compute_fCorr_perChip(
    const char* realFile = "Results_2sigma.root",
    const char* treeName = "noiseTree_Hyb1_SSA",
    const char* outfile  = "fCorr_perChip_Hyb1.txt"
) {
    const int nChips       = 8;
    const int nRowsPerChip = 120;

    // ============ OPEN REAL FILE ============
    TFile *file = TFile::Open(realFile);
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open " << realFile << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "ERROR: cannot find tree " << treeName << std::endl;
        file->Close();
        return;
    }

    int eventID;
    std::vector<short>* hitRows = nullptr;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("hitRows", &hitRows);

    // ============ BUILD hitsPerEvent PER CHIP ============
    std::map<int, std::map<int,int>> hitsPerEvent;  // chip -> (event -> hits)

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        if (!hitRows) continue;

        for (int chip = 0; chip < nChips; chip++) {
            int rowStart = chip * nRowsPerChip;
            int rowEnd   = rowStart + nRowsPerChip - 1;

            int hitsThisEntry = 0;

            for (size_t j = 0; j < hitRows->size(); j++) {
                short row = hitRows->at(j);

                if (row >= rowStart && row <= rowEnd)
                    hitsThisEntry++;
            }

            hitsPerEvent[chip][eventID] += hitsThisEntry;
        }
    }

    // ============ COMPUTE f_corr PER CHIP ============
    std::map<int,double> fCorr;

    std::cout << "\n=== Computing f_corr per chip ===\n" << std::endl;

    for (int chip = 0; chip < nChips; chip++) {
        
        // Build histogram for this chip
        TString hName = TString::Format("hHits_chip%d", chip);
        TH1F* h = new TH1F(hName, "tmp", 200, 0, 200);

        for (auto &kv : hitsPerEvent[chip])
            h->Fill(kv.second);

        double meanReal  = h->GetMean();
        double sigmaReal = h->GetStdDev();

        // Estimate random noise expectation
        double p = meanReal / nRowsPerChip;  // occupancy per channel

        double sigmaRand = sqrt(nRowsPerChip * p * (1 - p));

        double val = 0.0;
        if (sigmaReal > 0)
            val = 1.0 - (sigmaRand*sigmaRand)/(sigmaReal*sigmaReal);

        if (val < 0) val = 0.0;
        if (val > 1) val = 1.0;

        fCorr[chip] = val;

        std::cout << "Chip " << chip 
                  << "  mean=" << meanReal
                  << "  sigma=" << sigmaReal
                  << "  sigmaRand=" << sigmaRand
                  << "  f_corr=" << val
                  << std::endl;

        delete h;
    }

    // ============ WRITE OUTPUT TXT FILE ============
    std::ofstream out(outfile);
    out << "{\n";

    for (int chip = 0; chip < nChips; chip++) {
        TString key = TString::Format("\"chip%d\"", chip);
        out << "  " << key << ": " << fCorr[chip];

        if (chip < nChips - 1) out << ",";
        out << "\n";
    }

    out << "}\n";
    out.close();

    std::cout << "\nSaved f_corr values to " << outfile << std::endl;
    std::cout << "==========================================\n";

    file->Close();
}

