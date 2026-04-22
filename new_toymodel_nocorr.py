import numpy as np
import ast
import matplotlib.pyplot as plt
import ROOT

def load_thresholds(filename):
    """Load thresholds dictionary saved from ROOT macro."""
    with open(filename) as f:
        thresholds = ast.literal_eval(f.read())
    return thresholds  # dict: chip -> [thresholds per channel]

####adding correlations
def load_correlations(filename, n_chips=8):
    """Load correlation coefficients into a symmetric matrix."""
    with open(filename) as f:
        corr_dict = ast.literal_eval(f.read())

    corr_matrix = np.eye(n_chips)
    for key, val in corr_dict.items():
        c1, c2 = map(int, key.split("-"))
        corr_matrix[c1, c2] = val
        corr_matrix[c2, c1] = val
    return corr_matrix

#def simulate_noise_hits(thresholds, corr_matrix, n_events=10000, seed=123):
def simulate_noise_hits(thresholds, n_events=10000, seed=123):
    """
    Simulate noise hits per chip using thresholds.
    
    Args:
      thresholds: dict[chip] -> list of thresholds
      n_events: number of toy events to simulate
      seed: random seed
      
    Returns:
      hits_per_chip: dict[chip] -> array of length n_events with number of hits
    """
    rng = np.random.default_rng(seed)
    n_chips = len(thresholds)
    n_channels = len(next(iter(thresholds.values())))
#    hits_per_chip = {chip: np.zeros(n_events, dtype=int) for chip in thresholds}
#    chips = sorted(thresholds.keys())  # for consistent indexing

####no correlation
    hits_per_chip = {chip: np.zeros(n_events, dtype=int) for chip in thresholds}

    for chip, thr_list in thresholds.items():
        thr_array = np.array(thr_list)
        # For each event, generate Gaussian noise per channel
        for evt in range(n_events):
            noise = rng.normal(loc=0, scale=1, size=n_channels)
            hits = np.sum(noise > thr_array)
            hits_per_chip[chip][evt] = hits

    return hits_per_chip
#######

###with correlation
#    hits_per_chip = {chip: np.zeros(n_events, dtype=int) for chip in chips}
#
#    # Cholesky decomposition for correlated Gaussian
#    L = np.linalg.cholesky(corr_matrix)
#
#    for evt in range(n_events):
#        # event-level correlated Gaussian offsets (one per chip)
#        z = rng.normal(size=n_chips)
#        correlated_offsets = L @ z
#
#        for i, chip in enumerate(chips):
#            thr_array = np.array(thresholds[chip])
#            # add correlated offset to all channels of this chip
#            noise = rng.normal(loc=0, scale=1, size=n_channels) + correlated_offsets[i]
#            hits = np.sum(noise > thr_array)
#            hits_per_chip[chip][evt] = hits
#
#    return hits_per_chip
#######



if __name__ == "__main__":
    # Example usage
    #thresholds = load_thresholds("ThresholdsFromOcc_Hyb0.txt")
    thresholds = load_thresholds("ThresholdsFromOcc_Hyb1_2sigma.txt")
    #n_chips = len(thresholds)
    #corr_matrix = load_correlations("ChipCorrelationCoefficients_Hyb1_2sigma.txt", n_chips)
    #hits_per_chip = simulate_noise_hits(thresholds, corr_matrix, n_events=10000)
    hits_per_chip = simulate_noise_hits(thresholds, n_events=10000)

    # Print summary
    for chip in hits_per_chip:
    #for chip in hits:
        mean_hits = np.mean(hits_per_chip[chip])
        std_hits  = np.std(hits_per_chip[chip])
        print(f"Chip {chip}: mean hits = {mean_hits:.2f}, sigma = {std_hits:.2f}")

    #for chip, hits in hits_per_chip.items():
    #    plt.figure(figsize=(8,6))
    #    max_hits = hits.max()
    #    bins = np.arange(0, max_hits + 2) - 0.5  # integer bins
    #    plt.hist(hits, bins=bins, histtype='stepfilled', color='skyblue', edgecolor='blue', linewidth=1.5)
    #    plt.xlabel("Number of Hits per Event")
    #    plt.ylabel("Number of Events")
    #    plt.title(f"Chip {chip} - Toy Noise Hits")
    #    plt.grid(True, linestyle='--', alpha=0.6)
    #    plt.tight_layout()
    #    plt.yscale('log')  # optional: log-y scale like ROOT
    #    png_name = f"NoiseHits_Chip{chip}.png"
    #    plt.savefig(png_name, dpi=300)
    #    plt.close()
    #    print(f"Saved {png_name}")


    #root_file = ROOT.TFile("ToyHits_Hyb0_noCorr.root", "RECREATE")
    root_file = ROOT.TFile("ToyHits_Hyb1_2sigma.root", "RECREATE")
    for chip, hits in hits_per_chip.items():
        max_hits = int(hits.max())   # convert to int
        hName = f"hToy_chip{chip}"
        h = ROOT.TH1F(hName,
                  f"Toy Noise Hits Chip {chip};Hits per Event;Events",
                  max_hits + 1, float(-0.5), float(max_hits + 0.5))
        for hit in hits:
            h.Fill(int(hit))  # fill expects a number, ensure it's int
        h.Write()

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

