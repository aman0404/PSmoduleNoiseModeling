#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import ast

def load_corr_dict(txt_file):
    with open(txt_file, "r") as f:
        return ast.literal_eval(f.read())

corr_dict = load_corr_dict("/Users/amandeep/Desktop/OuterTracker/full_results/10k_events/1sigma/SubchipCorrelationCoefficients_16x16_Hyb1.txt")

  

n_chips = 16
R = np.eye(n_chips)
for k, v in corr_dict.items():
    i, j = map(int, k.split("-"))
    R[i, j] = v
    R[j, i] = v

chips = list(range(n_chips))

# --- Mask diagonal ---
R_masked = R.copy()
np.fill_diagonal(R_masked, np.nan)  # diagonal becomes NaN

fig, ax = plt.subplots(figsize=(6,5))
im = ax.imshow(R_masked, vmin=0, vmax=+1, cmap="coolwarm", origin="lower")
ax.set_xticks(chips)
ax.set_yticks(chips)
ax.set_xticklabels(chips)
ax.set_yticklabels(chips)
ax.set_xlabel("Chip i")
ax.set_ylabel("Chip j")
ax.set_title("Chip-by-chip Correlation")

# Annotate values
for i in range(n_chips):
    for j in range(n_chips):
        if i != j: 
           ax.text(j, i, f"{R[i,j]:.2f}", ha="center", va="center", fontsize=8)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Pearson r")
plt.tight_layout()
plt.savefig("subChip_corr_map__Hyb1_2sigma.png", dpi=200)
plt.close()
print("Saved chip_corr_map.png")

