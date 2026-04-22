void plotSSA_vs_SSA_Correlations_MultiPage(){
    TFile* file = TFile::Open("Results_2sigma.root");
    TTree* treeSSA = (TTree*)file->Get("noiseTree_Hyb1_SSA");

    int eventID;
    std::vector<short>* hitCols_SSA = nullptr;
    std::vector<short>* hitRows_SSA = nullptr;

    treeSSA->SetBranchAddress("eventID", &eventID);
    treeSSA->SetBranchAddress("hitCols", &hitCols_SSA);
    treeSSA->SetBranchAddress("hitRows", &hitRows_SSA);

    // Map: SSA chip -> map of eventID -> hit count
    std::map<int, std::map<int, int>> ssaHitsPerEvent;

    Long64_t nentries = treeSSA->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        treeSSA->GetEntry(i);

        for (int chip = 0; chip <= 7; ++chip) {
            int rowStart = chip * 120;
            int rowEnd = rowStart + 119;

            int hitsThisChip = 0;
            for (size_t j = 0; j < hitRows_SSA->size(); ++j) {
                short row = hitRows_SSA->at(j);
                if (row >= rowStart && row <= rowEnd)
                    hitsThisChip++;
            }
            ssaHitsPerEvent[chip][eventID] += hitsThisChip;
        }
    }

    TString pdfFile = "ssaHits/SSA_vs_SSA_Correlations_Hyb1_2sigma.pdf";
    bool firstPage = true;

    int nChips = 8;
    // Open file once outside the loops
    std::ofstream corrFile("ChipCorrelationCoefficients_Hyb1_2sigma.txt");
    corrFile << "{\n";


    for (int chip1 = 0; chip1 <= 6; ++chip1) {
        for (int chip2 = chip1 + 1; chip2 <= 7; ++chip2) {
            TString histName = Form("hCorr_SSA%d_SSA%d", chip1, chip2);
            TH2F* hCorr = new TH2F(histName,
                Form("SSA%d hits vs SSA%d hits;SSA%d hits;SSA%d hits",
                     chip2, chip1, chip1, chip2),
                100, 0, 100, 100, 0, 100);

            double sum_x = 0, sum_y = 0, sum_xx = 0, sum_yy = 0, sum_xy = 0;
            int n = 0;

            for (const auto& evtHits : ssaHitsPerEvent[chip1]) {
                int evt = evtHits.first;
                int hits1 = evtHits.second;
                int hits2 = ssaHitsPerEvent[chip2][evt];

                hCorr->Fill(hits1, hits2);

                sum_x += hits1;
                sum_y += hits2;
                sum_xx += hits1 * hits1;
                sum_yy += hits2 * hits2;
                sum_xy += hits1 * hits2;
                n++;
            }

            double corr = 0;
            if (n > 1) {
                double num = n * sum_xy - sum_x * sum_y;
                double den_x = n * sum_xx - sum_x * sum_x;
                double den_y = n * sum_yy - sum_y * sum_y;
                if (den_x > 0 && den_y > 0) {
                    corr = num / (sqrt(den_x) * sqrt(den_y));
                }
            }

            TCanvas* c = new TCanvas(Form("c_SSA%d_SSA%d", chip1, chip2),
                                     Form("SSA%d vs SSA%d Correlation", chip1, chip2),
                                     800, 600);
            hCorr->SetStats(0);
            hCorr->Draw("COLZ");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.DrawLatex(0.15, 0.85, Form("Pearson #rho = %.3f", corr));

            // Save coefficient into file (Python dict style)
            corrFile << "  \"" << chip1 << "-" << chip2 << "\": " << corr;
            if (!(chip1 == nChips-1 && chip2 == nChips-1)) corrFile << ",";
            corrFile << "\n";
            
            if (firstPage) {
                c->SaveAs(pdfFile + "(");
                firstPage = false;
            } else if (chip1 == 6 && chip2 == 7) {
                c->SaveAs(pdfFile + ")");
            } else {
                c->SaveAs(pdfFile);
            }

            delete hCorr;
            delete c;
        }
    }
            corrFile << "}\n";
            corrFile.close();


    file->Close();
}

