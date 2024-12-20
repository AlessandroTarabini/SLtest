// File: bootstrapping.C

void bootstrapping(int n_bootstraps = 500) {
    const std::vector<std::string> folders = {"htt_tt_0_8TeV", "htt_tt_1_8TeV", "htt_tt_2_8TeV"};
    
    // Open the original ROOT file in read mode
    TFile* original_file = TFile::Open("inputs/htt_tt.input.root", "READ");
    if (!original_file || original_file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Canvas to plot the histograms
    TCanvas* canvas = new TCanvas("canvas", "canvas", 1000, 1500);
    canvas->Divide(1, 3);

    TRandom3 randGen;

    for (size_t folder_idx = 0; folder_idx < folders.size(); ++folder_idx) {
        const std::string& folder_name = folders[folder_idx];
        std::cout << "Processing folder " << folder_name << std::endl;

        // Access the folder
        TDirectory* folder = (TDirectory*)original_file->Get(folder_name.c_str());
        if (!folder) {
            std::cerr << "Folder " << folder_name << " not found!" << std::endl;
            continue;
        }

        // Access the histogram
        TH1F* original_hist = (TH1F*)folder->Get("data_obs");
        if (!original_hist) {
            std::cerr << "Histogram data_obs not found in folder " << folder_name << std::endl;
            continue;
        }

        // Plot the inclusive dataset
        canvas->cd(folder_idx + 1);
        original_hist->Draw();

        // Create arrays to store x values and weights
        std::vector<double> x_values;
        std::vector<std::vector<double>> poiss_weights;

        // Fill x_values and generate Poisson weights
        for (int i = 1; i <= original_hist->GetNbinsX(); ++i) {
            int bin_content = static_cast<int>(original_hist->GetBinContent(i));
            if (bin_content == 0) continue; // If we do not have data, skip it, we do not generate data where there were not (TBC)
            double bin_center = original_hist->GetBinCenter(i);
            
            for (int j = 0; j < bin_content; ++j) {
                x_values.push_back(bin_center);
                std::vector<double> weights;
                for (int k = 0; k < n_bootstraps; ++k) {
                    weights.push_back(randGen.Poisson(1));
                }
                poiss_weights.push_back(weights);
            }
        }

        // For each bootstrap sample
        for (int bootstrap_idx = 0; bootstrap_idx < n_bootstraps; ++bootstrap_idx) {
            TFile* new_file = TFile::Open(Form("inputs/htt_tt.input_%d.root", bootstrap_idx), "UPDATE");
            new_file->mkdir(folder_name.c_str());
            new_file->cd(folder_name.c_str());

            TH1F* hist_new = (TH1F*)original_hist->Clone("data_obs");
            hist_new->Reset();

            // Fill the histogram with the current bootstrap sample
            std::vector<double> current_weights;
            for (size_t i = 0; i < x_values.size(); ++i) {
                current_weights.push_back(poiss_weights[i][bootstrap_idx]);
            }
            
            hist_new->FillN(x_values.size(), &x_values[0], &current_weights[0]);
            hist_new->Write();

            // Copy other TH1F histograms
            TIter next(folder->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)next())) {
                TObject* obj = key->ReadObj();
                if (obj->InheritsFrom(TH1F::Class()) && std::string(obj->GetName()) != "data_obs") {
                    obj->Write();
                }
                delete obj;
            }
            
            delete hist_new;
            new_file->Close();
            delete new_file;
        }
    }

    canvas->SaveAs("data_obs.pdf");
    std::cout << "Done!" << std::endl;
    original_file->Close();
    delete original_file;
    delete canvas;
    std::cout << "Original file closed" << std::endl;
}