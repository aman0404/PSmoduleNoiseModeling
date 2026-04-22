void plotSSAoccupancy_perChip() {
    TFile* file = TFile::Open("Results.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("noiseTree_Hyb0_SSA");
    if (!tree) {
        std::cerr << "Could not retrieve tree!" << std::endl;
        return;
    }

    int eventID;
    std::vector<short>* hitRows = nullptr;
    std::vector<short>* hitCols = nullptr;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("hitRows", &hitRows);
    tree->SetBranchAddress("hitCols", &hitCols);

    // loop over 8 SSA chips in this hybrid
    for (int currentChip = 0; currentChip < 8; ++currentChip) {

        int rowStart = currentChip * 120;
        int rowEnd   = rowStart + 119;

        // 120 channels per SSA chip
        TH1F* hOcc = new TH1F(Form("hOcc_chip%d", currentChip),
                              Form("SSA Occupancy - Chip %d;Channel;Occupancy", currentChip),
                              120, 0, 120);

        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; ++i) {
            tree->GetEntry(i);

            if (!hitRows) continue;

            for (size_t j = 0; j < hitRows->size(); ++j) {
                short row = hitRows->at(j);

                if (row >= rowStart && row <= rowEnd) {
                    int channel = row - rowStart; // 0..119 within chip
                    hOcc->Fill(channel);
                }
            }
        }

        // draw
        TCanvas* c = new TCanvas(Form("c_chip%d", currentChip),
                                 Form("SSA Chip %d Occupancy", currentChip), 1000, 600);
        hOcc->SetStats(0);
        hOcc->Draw("HIST");

        TString pdfName = Form("SSA_chip_%d_occupancy.pdf", currentChip);
        c->SaveAs(pdfName);

        delete c;
        delete hOcc;
    }

    file->Close();
}

