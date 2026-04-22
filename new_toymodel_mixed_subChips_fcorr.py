import numpy as np
import ast
import matplotlib.pyplot as plt
import ROOT

import numpy as np
import ast
import matplotlib.pyplot as plt

def load_fcorr(filename):
    import ast
    with open(filename, "r") as f:
        fcorr_dict = ast.literal_eval(f.read())

    f_corr = np.zeros(8)
    for i in range(8):
        f_corr[i] = fcorr_dict[f"chip{i}"]

    return f_corr


def load_thresholds(filename):
    with open(filename) as f:
        thresholds = ast.literal_eval(f.read())
    return thresholds  # dict: chip -> [thresholds per channel]

############################################################
# LOAD CORRELATION MATRIX (16×16)
############################################################

def load_correlations(filename, n_units=16):
    """
    Load a correlation dictionary of form:
        { "0-1": rho01, "0-2": rho02, ... }
    and return a 16×16 symmetric correlation matrix.
    """
    with open(filename, "r") as f:
        corr_dict = ast.literal_eval(f.read())

    corr_matrix = np.eye(n_units)

    for key, val in corr_dict.items():
        u1, u2 = map(int, key.split("-"))
        corr_matrix[u1, u2] = val
        corr_matrix[u2, u1] = val

    return corr_matrix

f_corr_chip = load_fcorr("/Users/amandeep/Desktop/OuterTracker/full_results/10k_events/pedenoise/Hybrid1_fcorr.txt")
#f_corr_chip = load_fcorr("/Users/amandeep/Desktop/OuterTracker/full_results/10k_events/1sigma/fCorr_perChip_Hyb1.txt")

############################################################
# SIMULATE NOISE HITS WITH 16×16 SUBCHIP CORRELATIONS
############################################################

def simulate_noise_hits(thresholds, corr_matrix, f_corr_chip, n_events=10000, seed=123):
#def simulate_noise_hits(thresholds, corr_matrix, n_events=10000, seed=123, f_corr=0.3):
    """
    thresholds: dict {chip → list of 120 threshold values}
    corr_matrix : 16×16 correlation matrix for 2×8 subchips
    f_corr: CM fraction to share between correlated & random noise
    """

    rng = np.random.default_rng(seed)
    
    chips = sorted(thresholds.keys())      # [0..7]
    n_chips = len(chips)                   # 8
    n_channels = len(next(iter(thresholds.values())))   # 120
    n_units = n_chips * 2                  # 16 subchips (two per chip)
    
    # Hits recorded for each chip and event
    hits_per_chip = {chip: np.zeros(n_events, dtype=int) for chip in chips}

    # Cholesky factor for 16×16 correlation matrix
    L = np.linalg.cholesky(corr_matrix)

    ########################################################
    # EVENT LOOP
    ########################################################

    for evt in range(n_events):

        # Generate correlated offsets for the 16 subchips
        z = rng.normal(size=n_units)
        correlated_offsets = L @ z

        for i, chip in enumerate(chips):

            thr_array = np.array(thresholds[chip])   # 120 thresholds
            
            # noise vector for all 120 channels
            noise_uncorr = rng.normal(size=n_channels)

            # Identify subchips for this chip
            sub0 = 2 * i       # region 0
            sub1 = 2 * i + 1   # region 1

            offset_region0 = correlated_offsets[sub0]
            offset_region1 = correlated_offsets[sub1]

            # Construct full-channel noise
            noise = np.zeros(n_channels)

            f_corr = f_corr_chip[chip]  
            # Region 0 → channels 0–59
            noise[0:60] = (
                np.sqrt(1 - f_corr) * noise_uncorr[0:60] +
                np.sqrt(f_corr) * offset_region0
            )

            # Region 1 → channels 60–119
            noise[60:120] = (
                np.sqrt(1 - f_corr) * noise_uncorr[60:120] +
                np.sqrt(f_corr) * offset_region1
            )

            # Count hits
            hits = np.sum(noise > thr_array)
            hits_per_chip[chip][evt] = hits

    return hits_per_chip


############################################################
# PLOT RESULTS
############################################################

def plot_results(hits_per_chip, outname="toy_hits.png"):
    """
    Plot histograms of hits per chip.
    """

    plt.figure(figsize=(10, 6))

    for chip, hits in hits_per_chip.items():
        plt.hist(hits, bins=60, alpha=0.5, label=f"Chip {chip}")

    plt.xlabel("Hits per event")
    plt.ylabel("Events")
    plt.legend()
    plt.grid(True)
    plt.savefig(outname, dpi=200)
    plt.close()


############################################################
# MAIN EXAMPLE (you may modify)
############################################################

if __name__ == "__main__":
    
    thresholds = load_thresholds("ThresholdsFromOcc_Hyb1_2sigma.txt")

    # Load 16×16 correlation matrix
    corr_matrix = load_correlations("/Users/amandeep/Desktop/OuterTracker/full_results/10k_events/1sigma/SubchipCorrelationCoefficients_16x16_Hyb1.txt")

    # Simulate hits
    hits_per_chip = simulate_noise_hits(
        thresholds,
        corr_matrix,
        f_corr_chip,
        n_events=50000,
        seed=42
    )


############################################################
# SUMMARY PRINT
############################################################

print("\n======== Toy Model Summary ========\n")
for chip in sorted(hits_per_chip.keys()):
    mean_hits = np.mean(hits_per_chip[chip])
    std_hits  = np.std(hits_per_chip[chip])
    print(f"Chip {chip}: mean hits = {mean_hits:.2f}, sigma = {std_hits:.2f}")

############################################################
# PER-CHIP HISTOGRAM PLOTS (PNG)
############################################################

for chip, hits in hits_per_chip.items():
    plt.figure(figsize=(8,6))
    max_hits = hits.max()
    bins = np.arange(0, max_hits + 2) - 0.5  # integer bins
    
    plt.hist(hits, bins=bins, histtype='stepfilled',
             color='skyblue', edgecolor='blue', linewidth=1.5)
    
    plt.xlabel("Number of Hits per Event")
    plt.ylabel("Number of Events")
    plt.title(f"Chip {chip} - Toy Noise Hits")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.yscale('log')  # optional
    
    png_name = f"NoiseHits_Chip{chip}.png"
    #plt.savefig(png_name, dpi=300)
    plt.close()

    print(f"Saved {png_name}")

############################################################
# SAVE ROOT FILE WITH HISTOGRAMS + TTREE
############################################################

root_file = ROOT.TFile("ToyHits_Hyb1_2sigma_subChip_f_corr.root", "RECREATE")

# 1. Save TH1F histograms
for chip, hits in hits_per_chip.items():
    max_hits = int(np.max(hits))
    hName = f"hToy_chip{chip}"
    
    h = ROOT.TH1F(hName,
                  f"Toy Noise Hits Chip {chip};Hits per Event;Events",
                  max_hits + 1, -0.5, max_hits + 0.5)
    
    for hit in hits:
        h.Fill(int(hit))
    
    h.Write()

# 2. Create TTree storing per-chip hits
tree = ROOT.TTree("toyTree", "Toy model hits per chip")

n_chips = len(hits_per_chip)
hit_buffers = {}

# numpy arrays as buffers for ROOT branches
for chip in sorted(hits_per_chip.keys()):
    hit_buffers[chip] = np.zeros(1, dtype=np.int32)
    tree.Branch(f"hits_chip{chip}", hit_buffers[chip],
                f"hits_chip{chip}/I")

n_events = len(next(iter(hits_per_chip.values())))
for evt in range(n_events):
    for chip in sorted(hits_per_chip.keys()):
        hit_buffers[chip][0] = int(hits_per_chip[chip][evt])
    tree.Fill()

tree.Write()
root_file.Close()

print("\nSaved ToyHits_Hyb1_2sigma.root with histograms and toyTree")
print("==============================================\n")


