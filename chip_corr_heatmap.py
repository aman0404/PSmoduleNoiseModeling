#!/usr/bin/env python3
import argparse, os, sys, ast
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import ROOT
    have_root = True
except Exception:
    have_root = False
    print("[warn] PyROOT not found; ROOT saving will be skipped.", file=sys.stderr)

# ---------------------------
# Helpers: load from ROOT + optional occupancy
# ---------------------------
def load_Kmatrix_from_root(rootfile, treename="toyTree"):
    """
    Reads toyTree where branches are hits_chip0 ... hits_chipN.
    Returns:
      K : (n_events, n_chips) array
      chips : sorted list of chip ids
    """
    f = ROOT.TFile.Open(rootfile, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {rootfile}")
    t = f.Get(treename)
    if not t:
        raise RuntimeError(f"Tree '{treename}' not found in {rootfile}")

    # discover branches
    brs = [br.GetName() for br in t.GetListOfBranches()]
    chip_branches = sorted([b for b in brs if b.startswith("hits_chip")],
                           key=lambda s: int(s.replace("hits_chip","")))
    if not chip_branches:
        raise RuntimeError("No branches named hits_chipX found.")

    chips = [int(b.replace("hits_chip","")) for b in chip_branches]
    n_chips = len(chips)
    n_events = t.GetEntries()

    # read into matrix
    K = np.zeros((n_events, n_chips), dtype=float)
    for ievt in range(n_events):
        t.GetEntry(ievt)
        for j, chip in enumerate(chips):
            K[ievt, j] = getattr(t, f"hits_chip{chip}")

    f.Close()
    return K, chips

def load_occupancy(path, expect_chips=None, expect_channels=120):
    """
    Occupancy dict file saved by ROOT macro: {chip: [p_ch0..p_ch119]}
    """
    with open(path) as f:
        d = ast.literal_eval(f.read())
    occ = {int(k): np.array(v, dtype=float) for k, v in d.items()}
    if expect_chips is not None:
        for c in expect_chips:
            if c not in occ:
                raise ValueError(f"Occupancy file missing chip {c}")
            if occ[c].shape[0] != expect_channels:
                raise ValueError(f"Occupancy chip {c} length {occ[c].shape[0]} != {expect_channels}")
    return occ

# ---------------------------
# Correlations
# ---------------------------
def corr_from_K(K_matrix):
    """Pearson correlation across chips (columns of K)."""
    return np.corrcoef(K_matrix, rowvar=False)

def corr_residualized(K_matrix, occupancy_dict, chip_ids):
    """
    Residualized Pearson after removing Poisson-binomial expectation:
    Z_e^(i) = (K - sum p) / sqrt(sum p(1-p))
    """
    Z = np.zeros_like(K_matrix, dtype=float)
    for col, chip in enumerate(chip_ids):
        p = occupancy_dict[chip]
        mu = float(np.sum(p))                  # E[K]
        var0 = float(np.sum(p * (1.0 - p)))    # Var[K]
        sd0 = np.sqrt(var0) if var0 > 0 else 1.0
        Z[:, col] = (K_matrix[:, col] - mu) / sd0
    return np.corrcoef(Z, rowvar=False)

# ---------------------------
# Plotting
# ---------------------------
def plot_corr_heatmap(ax, R, chips, title, annotate=True):
    im = ax.imshow(R, vmin=-1.0, vmax=+1.0, origin='lower', cmap="coolwarm")
    ax.set_title(title)
    ax.set_xlabel("Chip")
    ax.set_ylabel("Chip")
    ax.set_xticks(range(len(chips)))
    ax.set_yticks(range(len(chips)))
    ax.set_xticklabels(chips)
    ax.set_yticklabels(chips)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r")
    if annotate:
        n = R.shape[0]
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{R[i,j]:.2f}", ha="center", va="center", fontsize=8)

# ---------------------------
# Save matrices to ROOT (optional)
# ---------------------------
def write_corr_to_root(R_mats, names, chips, outroot):
    if not have_root:
        print("[root] skipped (PyROOT not available)")
        return
    mode = "RECREATE" if not os.path.exists(outroot) else "UPDATE"
    f = ROOT.TFile(outroot, mode)
    n = len(chips)
    for R, nm in zip(R_mats, names):
        h2 = ROOT.TH2D(nm, f"{nm};Chip;Chip", n, -0.5, n-0.5, n, -0.5, n-0.5)
        for i in range(n):
            for j in range(n):
                h2.SetBinContent(i+1, j+1, float(R[i,j]))
        for i in range(n):
            h2.GetXaxis().SetBinLabel(i+1, f"{chips[i]}")
            h2.GetYaxis().SetBinLabel(i+1, f"{chips[i]}")
        h2.Write()
    f.Close()
    print(f"[root] wrote {', '.join(names)} to {outroot}")

# ---------------------------
# CLI
# ---------------------------
def main():
    ap = argparse.ArgumentParser(description="Build chip-by-chip correlation heatmaps from ToyHits ROOT file")
    ap.add_argument("--root", required=True, help="Toy ROOT file (e.g. ToyHits_Hyb1_2sigma.root)")
    ap.add_argument("--tree", default="toyTree", help="Tree name (default: toyTree)")
    ap.add_argument("--occupancy", help="Optional Occupancy dict (for residualized corr)")
    ap.add_argument("--pdf", required=True, help="Output multi-page PDF (heatmaps)")
    ap.add_argument("--rootout", help="Optional ROOT output to store TH2 correlation matrices")
    ap.add_argument("--no-annotate", action="store_true", help="Disable cell value annotation")
    args = ap.parse_args()

    # Load K (events × chips)
    if not have_root:
        raise RuntimeError("PyROOT is required to read the toy ROOT file.")
    K, chips = load_Kmatrix_from_root(args.root, treename=args.tree)
    print(f"[info] loaded K with shape {K.shape} for chips {chips}")

    # Raw correlation
    R_raw = corr_from_K(K)

    # Residualized correlation (optional)
    R_resid = None
    if args.occupancy:
        occ = load_occupancy(args.occupancy, expect_chips=chips)
        R_resid = corr_residualized(K, occ, chip_ids=chips)

    # Write multipage PDF
    with PdfPages(args.pdf) as pdf:
        # Page 1: raw
        fig, ax = plt.subplots(figsize=(6.5,5.5))
        plot_corr_heatmap(ax, R_raw, chips, "Chip-by-chip correlation (raw)",
                          annotate=(not args.no_annotate))
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2: residualized (if provided)
        if R_resid is not None:
            fig, ax = plt.subplots(figsize=(6.5,5.5))
            plot_corr_heatmap(ax, R_resid, chips, "Chip-by-chip correlation (residualized)",
                              annotate=(not args.no_annotate))
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    print(f"[pdf] wrote {args.pdf}")

    # Optional: save to ROOT as TH2
    if args.rootout:
        mats = [R_raw]
        names = ["R_raw"]
        if R_resid is not None:
            mats.append(R_resid); names.append("R_resid")
        write_corr_to_root(mats, names, chips, args.rootout)

if __name__ == "__main__":
    main()

