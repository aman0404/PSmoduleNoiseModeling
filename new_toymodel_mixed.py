import numpy as np
import ast
import matplotlib.pyplot as plt
import ROOT

def load_thresholds(filename):
    with open(filename) as f:
        thresholds = ast.literal_eval(f.read())
    return thresholds  # dict: chip -> [thresholds per channel]

####adding correlations
def load_correlations(filename, n_chips=8):
    with open(filename) as f:
        corr_dict = ast.literal_eval(f.read())

    corr_matrix = np.eye(n_chips)
    
    for key, val in corr_dict.items():
        c1, c2 = map(int, key.split("-"))
        corr_matrix[c1, c2] = val
        corr_matrix[c2, c1] = val
    
    #print("corr_matrix ", corr_matrix)
    return corr_matrix

def simulate_noise_hits(thresholds, corr_matrix, n_events=10000, seed=123, f_corr=0.5):
    rng = np.random.default_rng(seed)
    n_chips = len(thresholds)
    print("n_chips ", n_chips)
    n_channels = len(next(iter(thresholds.values())))
    chips = sorted(thresholds.keys())

    hits_per_chip = {chip: np.zeros(n_events, dtype=int) for chip in chips}

    # Cholesky decomposition for correlated Gaussian
    L = np.linalg.cholesky(corr_matrix)


    for evt in range(n_events):
        # chip-level correlated Gaussian offsets
        z = rng.normal(size=n_chips)
        correlated_offsets = L @ z
        for i, chip in enumerate(chips):
            thr_array = np.array(thresholds[chip])

            # Uncorrelated channel noise
            noise_uncorr = rng.normal(loc=0, scale=1, size=n_channels)
            # Blend correlated and uncorrelated parts
            #noise = (1 - f_corr) * noise_uncorr + f_corr * (noise_uncorr + correlated_offsets[i])
            noise = np.sqrt(1.0 - f_corr) * noise_uncorr + np.sqrt(f_corr) * correlated_offsets[i]

            hits = np.sum(noise > thr_array)
            hits_per_chip[chip][evt] = hits

    return hits_per_chip

#######



if __name__ == "__main__":
    # Example usage
    thresholds = load_thresholds("ThresholdsFromOcc_Hyb0_2sigma.txt")
    #hits = simulate_noise_hits(thresholds, n_events=1000)
    n_chips = len(thresholds)
    corr_matrix = load_correlations("ChipCorrelationCoefficients_Hyb0_2sigma.txt", n_chips)
    hits_per_chip = simulate_noise_hits(thresholds, corr_matrix, n_events=10000, f_corr=0.35) # try 30% correlated, 70% uncorrelated
    #hits_per_chip = simulate_noise_hits(thresholds, corr_matrix, n_events=10000, f_corr=0.30) # try 30% correlated, 70% uncorrelated
    #hits_per_chip = simulate_noise_hits(thresholds, corr_matrix, n_events=10000)  ##default 10k, 10000
    #hits_per_chip = simulate_noise_hits(thresholds, n_events=10000)

    # Print summary
    for chip in hits_per_chip:
    #for chip in hits:
        mean_hits = np.mean(hits_per_chip[chip])
        std_hits  = np.std(hits_per_chip[chip])
        print(f"Chip {chip}: mean hits = {mean_hits:.2f}, sigma = {std_hits:.2f}")

    for chip, hits in hits_per_chip.items():
        plt.figure(figsize=(8,6))
        max_hits = hits.max()
        bins = np.arange(0, max_hits + 2) - 0.5  # integer bins
        plt.hist(hits, bins=bins, histtype='stepfilled', color='skyblue', edgecolor='blue', linewidth=1.5)
        plt.xlabel("Number of Hits per Event")
        plt.ylabel("Number of Events")
        plt.title(f"Chip {chip} - Toy Noise Hits")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.yscale('log')  # optional: log-y scale like ROOT
        png_name = f"NoiseHits_Chip{chip}.png"
        #plt.savefig(png_name, dpi=300)
        plt.close()
        print(f"Saved {png_name}")


    root_file = ROOT.TFile("ToyHits_Hyb0_2sigma.root", "RECREATE")
    for chip, hits in hits_per_chip.items():
        max_hits = int(hits.max())   # convert to int
        hName = f"hToy_chip{chip}"
        h = ROOT.TH1F(hName,
                  f"Toy Noise Hits Chip {chip};Hits per Event;Events",
                  max_hits + 1, float(-0.5), float(max_hits + 0.5))
        for hit in hits:
            h.Fill(int(hit))  # fill expects a number, ensure it's int
        h.Write()

    tree = ROOT.TTree("toyTree", "Toy model hits per chip")

    n_chips = len(hits_per_chip)
    # Create numpy arrays as buffers for each branch
    hit_buffers = {}
    for chip in sorted(hits_per_chip.keys()):
        hit_buffers[chip] = np.zeros(1, dtype=np.int32)
        tree.Branch(f"hits_chip{chip}", hit_buffers[chip], f"hits_chip{chip}/I")

    n_events = len(next(iter(hits_per_chip.values())))
    for evt in range(n_events):
        for chip in sorted(hits_per_chip.keys()):
            hit_buffers[chip][0] = int(hits_per_chip[chip][evt])
        tree.Fill()

    tree.Write()
    root_file.Close()

#    print("Saved ToyHits_Hyb1_2sigma.root with histograms and toyTree")

    # Plot histogram of hits per event per chip
#    plt.figure(figsize=(12,6))
#    colors = plt.cm.tab10.colors  # color for each chip
#
#    max_hits = 0
#    for chip, hits in hits_per_chip.items():
#        max_hits = max(max_hits, hits.max())
#
#    bins = np.arange(0, max_hits+2) - 0.5  # integer bins
#
#    for chip, hits in hits_per_chip.items():
#        plt.hist(hits, bins=bins, histtype='step', linewidth=2, color=colors[chip % 10], label=f'Chip {chip}')
#
#    plt.xlabel("Number of Hits per Event")
#    plt.ylabel("Number of Events")
#    plt.title("Toy-model Noise Hits per Chip")
#    plt.legend()
#    plt.grid(True)
#    plt.tight_layout()
#    plt.show()

