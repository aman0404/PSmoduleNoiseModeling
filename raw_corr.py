#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import ast

def load_corr_dict(txt_file):
    with open(txt_file, "r") as f:
        return ast.literal_eval(f.read())

corr_dict = load_corr_dict("ChipCorrelationCoefficients_Hyb1_2sigma.txt")

#corr_dict = {
#  "0-1": 0.847571,
#  "0-2": 0.733431,
#  "0-3": 0.659392,
#  "0-4": 0.584195,
#  "0-5": 0.587293,
#  "0-6": 0.54212,
#  "0-7": 0.569727,
#  "1-2": 0.783561,
#  "1-3": 0.701717,
#  "1-4": 0.563087,
#  "1-5": 0.580365,
#  "1-6": 0.555861,
#  "1-7": 0.639052,
#  "2-3": 0.826283,
#  "2-4": 0.59398,
#  "2-5": 0.542605,
#  "2-6": 0.564796,
#  "2-7": 0.625263,
#  "3-4": 0.651323,
#  "3-5": 0.56787,
#  "3-6": 0.598003,
#  "3-7": 0.619405,
#  "4-5": 0.756444,
#  "4-6": 0.737042,
#  "4-7": 0.617137,
#  "5-6": 0.764437,
#  "5-7": 0.695462,
#  "6-7": 0.718911,
#}
 
  

n_chips = 8
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
plt.savefig("chip_corr_map__Hyb1_2sigma.png", dpi=200)
plt.close()
print("Saved chip_corr_map.png")

