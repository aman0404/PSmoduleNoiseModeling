#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TMath.h>
#include <TLatex.h>

#include <vector>
#include <set>
#include <map>
#include <iostream>

// Reads SSA hits from a single hybrid tree (e.g. "noiseTree_Hyb0_SSA" or "noiseTree_Hyb1_SSA")
// Builds per-chip, per-channel occupancy, converts to thresholds, and plots:
//    X = channel (0..119), Y = Q_thr = sqrt(2)*Erf^{-1}(1-2*occupancy)
// One canvas with 8 pads (2x4), one per SSA chip.
void plotThresholds_perChip(const char* filename = "Results_3sigma.root",
                                    const char* treename = "noiseTree_Hyb0_SSA",
                                    const char* outbase  = "ThresholdsFromOcc_Hyb0_3sigma")
{
    gStyle->SetOptStat(0);

    // Open file and tree
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get(treename);
    if (!tree) {
        std::cerr << "Error: could not find tree " << treename << " in file " << filename << std::endl;
        file->Close();
        return;
    }

    // Branches we need
    int eventID = -1;
    std::vector<short>* hitRows = nullptr;   // SSA has 1 column, 120 rows per chip; row encodes chip*120 + channel
    std::vector<short>* hitCols = nullptr;   // not used for SSA (1 col), but we bind to be safe if present

    tree->SetBranchAddress("eventID",  &eventID);
    tree->SetBranchAddress("hitRows",  &hitRows);
    // If the branch doesn't exist, ROOT will keep hitCols as nullptr (fine)
    tree->SetBranchAddress("hitCols",  &hitCols);

    const int nChips = 8;
    const int nChan  = 120;

    // Collect unique events and, for each chip, map: eventID -> set of channels that fired
    std::set<int> uniqueEvents;
    std::map<int, std::map<int, std::set<int>>> chansPerEventPerChip; // [chip][eventID] -> {channels}

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if (!hitRows) continue;

        uniqueEvents.insert(eventID);

        // For each hit row, decode chip & channel and record into the set (deduplicates within an event)
        for (size_t j = 0; j < hitRows->size(); ++j) {
            int row = static_cast<int>(hitRows->at(j));       // 0..959
            if (row < 0) continue;
            int chip    = row / nChan;                        // 0..7
            int channel = row % nChan;                        // 0..119
            if (chip < 0 || chip >= nChips) continue;
            if (channel < 0 || channel >= nChan) continue;

            chansPerEventPerChip[chip][eventID].insert(channel);
        }
    }

    const int nEvents = static_cast<int>(uniqueEvents.size());
    if (nEvents == 0) {
        std::cerr << "No events found. Aborting.\n";
        file->Close();
        return;
    }

    // Compute occupancy per chip/channel: p = (# events where channel fired) / nEvents
    // Then convert to threshold: Q_thr = sqrt(2)*Erf^{-1}(1 - 2*p)
    // Protect against p=0 or p=1 with a small epsilon clamp.
    const double eps = 1e-12;

    // thresholds[chip][channel]
    std::vector<std::vector<double>> thresholds(nChips, std::vector<double>(nChan, 0.0));

    for (int chip = 0; chip < nChips; ++chip) {
        // Count events per channel
        std::vector<int> firedCount(nChan, 0);

        for (const auto& evt_pair : chansPerEventPerChip[chip]) {
            const std::set<int>& chset = evt_pair.second;
            for (int ch : chset) {
                if (ch >= 0 && ch < nChan) firedCount[ch]++;
            }
        }

        // Convert occupancy to threshold
        for (int ch = 0; ch < nChan; ++ch) {
            double p = (nEvents > 0) ? static_cast<double>(firedCount[ch]) / nEvents : 0.0;
            // clamp to (eps, 1-eps) to keep ErfInverse finite
            if (p < eps) p = eps;
            if (p > 1.0 - eps) p = 1.0 - eps;
            thresholds[chip][ch] = std::sqrt(2.0) * TMath::ErfInverse(1.0 - 2.0 * p);
        }
    }
    // ---------- Save thresholds as Python dictionary ----------
    std::ofstream fout(Form("%s.txt", outbase));
    fout << "{\n";
    for (int chip = 0; chip < nChips; ++chip) {
        fout << "  " << chip << ": [";
        for (int ch = 0; ch < nChan; ++ch) {
            fout << thresholds[chip][ch];
            if (ch < nChan - 1) fout << ", ";
        }
        fout << "]";
        if (chip < nChips - 1) fout << ",";
        fout << "\n";
    }
    fout << "}\n";
    fout.close();
    std::cout << "Thresholds saved to " << outbase << ".txt\n";

    // ---------- Plot: one pad per chip (2x4) ----------
    TCanvas* c = new TCanvas("c", "SSA Threshold vs Channel (from occupancy)", 1400, 900);
    c->Divide(2, 4);

    for (int chip = 0; chip < nChips; ++chip) {
        c->cd(chip + 1);

        // Build graph (channel -> threshold)
        TGraph* gr = new TGraph(nChan);
        for (int ch = 0; ch < nChan; ++ch) {
            gr->SetPoint(ch, ch, thresholds[chip][ch]);
        }

        gr->SetTitle(Form("SSA Chip %d;Channel;Threshold (units of #sigma)", chip));
        gr->SetLineColor(kBlue + (chip % 4));
        gr->SetLineWidth(2);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);
        gr->SetMarkerColor(kBlue + (chip % 4));

        gr->Draw("APL");

        // Right-side Y axis ticks/labels
        gPad->Update();
        TGaxis* rightAxis = new TGaxis(gPad->GetUxmax(),
                                       gPad->GetUymin(),
                                       gPad->GetUxmax(),
                                       gPad->GetUymax(),
                                       gPad->GetUymin(),
                                       gPad->GetUymax(),
                                       510, "+L");
        rightAxis->SetTitle("Threshold (units of #sigma)");
        rightAxis->SetTitleOffset(1.2);
        rightAxis->SetLabelSize(0.035);
        rightAxis->Draw();

        // Small note with the formula used
        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.04);
        lat.DrawLatex(0.12, 0.88, "#it{Q_{thr} = #sqrt{2} Erf^{-1}(1-2p)}");
    }

    // Save
    TString pdfName = TString::Format("%s.pdf", outbase);
    TString pngName = TString::Format("%s.png", outbase);
    c->SaveAs(pdfName);
    c->SaveAs(pngName);

    file->Close();
}

