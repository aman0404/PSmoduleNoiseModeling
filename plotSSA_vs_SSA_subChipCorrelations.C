void plotSSA_vs_SSA_subChipCorrelations(){
    TFile* file = TFile::Open("Results_2sigma.root");
    TTree* treeSSA = (TTree*)file->Get("noiseTree_Hyb1_SSA");

    int eventID;
    std::vector<short>* hitCols_SSA = nullptr;
    std::vector<short>* hitRows_SSA = nullptr;

    treeSSA->SetBranchAddress("eventID", &eventID);
    treeSSA->SetBranchAddress("hitCols", &hitCols_SSA);
    treeSSA->SetBranchAddress("hitRows", &hitRows_SSA);

// subchip = chip*2 + region
// region 0 = channels [0..59] of this chip
// region 1 = channels [60..119] of this chip

std::map<int, std::map<int,int>> subchipHitsPerEvent; // subchip -> (event -> hits)

Long64_t nentries = treeSSA->GetEntries();
for (Long64_t i = 0; i < nentries; ++i) {
    treeSSA->GetEntry(i);

    for (int chip = 0; chip < 8; chip++) {
        for (int region = 0; region < 2; region++) {

            int sub = chip * 2 + region;   // subchip index 0..15

            int rowStart = chip * 120 + (region == 0 ? 0 : 60);
            int rowEnd   = rowStart + 59;

            int hitsThisSubchip = 0;

            for (size_t j = 0; j < hitRows_SSA->size(); ++j) {
                short row = hitRows_SSA->at(j);
                if (row >= rowStart && row <= rowEnd)
                    hitsThisSubchip++;
            }

            subchipHitsPerEvent[sub][eventID] += hitsThisSubchip;
        }
    }
}

/////////////////////////////////////////////////////////
// 16×16 CORRELATION MATRIX
/////////////////////////////////////////////////////////

TString pdfFile = "ssaHits/SSA_subchip16_Correlations_Hyb1.pdf";
bool firstPage = true;

const int nSubchips = 16;

ofstream corrFile("SubchipCorrelationCoefficients_16x16_Hyb1.txt");
corrFile << "{\n";

for (int s1 = 0; s1 < nSubchips; s1++) {
    for (int s2 = s1+1; s2 < nSubchips; s2++) {

        TString histName = Form("hCorr_sub%d_sub%d", s1, s2);
        TH2F* hCorr = new TH2F(histName,
            Form("Subchip %d hits vs Subchip %d hits;Subchip %d hits;Subchip %d hits",
                 s1, s2, s1, s2),
            50, 0, 50, 50, 0, 50);

        double sum_x = 0, sum_y = 0, sum_xx = 0, sum_yy = 0, sum_xy = 0;
        int n = 0;

        // Loop over events for subchip s1
        for (const auto& evtHits : subchipHitsPerEvent[s1]) {
            int evt = evtHits.first;
            int hits1 = evtHits.second;
            int hits2 = subchipHitsPerEvent[s2][evt];

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
            double num   = n * sum_xy - sum_x * sum_y;
            double den_x = n * sum_xx - sum_x * sum_x;
            double den_y = n * sum_yy - sum_y * sum_y;
        
            if (den_x > 0 && den_y > 0) {
                corr = num / (sqrt(den_x * den_y));
            }
        }

        TCanvas* c = new TCanvas(Form("c_sub%d_sub%d", s1, s2),
                                 Form("Subchip %d vs %d", s1, s2), 800, 600);
        hCorr->SetStats(0);
        hCorr->Draw("COLZ");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.15, 0.85, Form("Pearson #rho = %.3f", corr));

        corrFile << Form("  \"%d-%d\": %.4f,\n", s1, s2, corr);

        if (firstPage) {
            c->SaveAs(pdfFile + "(");
            firstPage = false;
        } else if (s1 == nSubchips-2 && s2 == nSubchips-1) {
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
}
