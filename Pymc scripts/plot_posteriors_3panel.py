import arviz as az
import matplotlib.pyplot as plt
import numpy as np

# Set style for publication-quality figures with larger fonts
plt.style.use('default')
plt.rcParams['font.size'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.titlesize'] = 18

# 1. load your trace
idata = az.from_netcdf("/blue/juannanzhou/fahimeh.rahimi/O2_data/Results/trace_t2000_d7000_20250518_223910.nc")

# 2. pull out posterior samples (flatten chains & draws)
mu1  = idata.posterior["mu_g1"].values.flatten()
mu2  = idata.posterior["mu_g2"].values.flatten()
sig1 = idata.posterior["sigma_g1"].values.flatten()
sig2 = idata.posterior["sigma_g2"].values.flatten()
sig_e = idata.posterior["sigma_e"].values.flatten()

# 3. Calculate ratios for third panel
sigma_g1_2 = sig1 ** 2
sigma_g2_2 = sig2 ** 2
sigma_e2 = sig_e ** 2
ratio_g1_e = sigma_g1_2 / sigma_e2
ratio_g2_e = sigma_g2_2 / sigma_e2

# 4. set up a 1×3 grid for plotting
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Define colors for better contrast
color1 = '#1f77b4'  # Blue
color2 = '#ff7f0e'  # Orange

# --- first plot: both mu's together ---
axes[0].hist(mu1, bins=60, density=True, alpha=0.7, color=color1, label="O1", edgecolor='black', linewidth=0.5)
axes[0].hist(mu2, bins=60, density=True, alpha=0.7, color=color2, label="O2", edgecolor='black', linewidth=0.5)
axes[0].set_xlabel("μ")
axes[0].set_ylabel("Density")
axes[0].text(0.05, 0.95, 'a', transform=axes[0].transAxes, fontsize=16, fontweight='bold')
axes[0].legend(loc='upper right', frameon=True, fancybox=True, shadow=True)

# --- second plot: both sigmas together ---
axes[1].hist(sig1, bins=60, density=True, alpha=0.7, color=color1, label="O1", edgecolor='black', linewidth=0.5)
axes[1].hist(sig2, bins=60, density=True, alpha=0.7, color=color2, label="O2", edgecolor='black', linewidth=0.5)
axes[1].set_xlabel("σ")
axes[1].text(0.05, 0.95, 'b', transform=axes[1].transAxes, fontsize=16, fontweight='bold')

# --- third plot: both ratio distributions with medians ---
axes[2].hist(ratio_g1_e, bins=60, density=True, alpha=0.7, color=color1, label="O1", edgecolor='black', linewidth=0.5)
axes[2].hist(ratio_g2_e, bins=60, density=True, alpha=0.7, color=color2, label="O2", edgecolor='black', linewidth=0.5)

axes[2].set_xlabel("σ²_g / σ²_e")
axes[2].text(0.05, 0.95, 'c', transform=axes[2].transAxes, fontsize=16, fontweight='bold')
# Make x-axis ticks less dense
axes[2].xaxis.set_major_locator(plt.MaxNLocator(5))

plt.tight_layout()
plt.show() 