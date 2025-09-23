#!/usr/bin/env python

import os, time, datetime
import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as tt
import arviz as az
import matplotlib.pyplot as plt

def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)

def resample_mutations(
    observed_mutations, 
    false_positive_lambda, 
    mutation_rate, 
    failure_to_recall_rate, 
    genome_length=1e9, 
    n_generations=150
):
    # Replace NaN with 0, and clip negatives to 0
    n_obs = np.array(observed_mutations)
    n_obs = np.nan_to_num(n_obs, nan=0.0)
    n_obs = np.clip(n_obs, 0, None)
    n_obs = n_obs.astype(int)
    
    # False positives: Poisson distribution with expectation of false_positive_lambda
    # But limit to not exceed observed count to avoid negative results
    n_false_positives = np.random.poisson(false_positive_lambda, size=n_obs.shape)
    n_false_positives = np.minimum(n_false_positives, n_obs)  # Don't subtract more than observed
    
    # False negatives (missed mutations): Poisson distribution
    # mutation_rate * n_generations * genome_length * failure_to_recall_rate
    n_missed = np.random.poisson(
        mutation_rate * n_generations * genome_length * failure_to_recall_rate, 
        size=n_obs.shape
    )
    
    # Final corrected count: observed - false_positives + missed_mutations
    n_corrected = n_obs - n_false_positives + n_missed
    
    # Ensure non-negative values
    n_corrected = np.clip(n_corrected, 0, None)
    
    return n_corrected

def build_G_matrices(df):
    """Build G1 and G2 matrices from a dataframe"""
    # G1: only rows with G150 or G300
    mask = df['Genassay'].isin(["G150","G300"])
    g1 = df.loc[mask].reset_index()
    bases = g1['Base'].unique()
    n1 = len(bases)
    
    G1_small = np.zeros((n1, n1))
    for i, base in enumerate(bases):
        G1_small[i, i] = g1.loc[g1['Base']==base, 'O1_totalmutsx'].iloc[0]
    
    base_to_idx = {b:i for i,b in enumerate(bases)}
    idx_G1 = np.array([
        base_to_idx[row['Base']] if row['Genassay'] in ["G150","G300"] else -1
        for _, row in df.iterrows()
    ], dtype=int)
    
    # G2: only rows with G300
    g2 = df[df['Genassay']=="G300"].reset_index()
    G2_full = np.diag(g2['O2_totalmutsx'].values)
    keep = ~g2.duplicated(['Base','Rep'])
    sel_small = np.where(keep)[0]
    G2_small = G2_full[np.ix_(sel_small,sel_small)]
    ids_small = g2.loc[sel_small,'Sublinex'].astype(str).values
    lookup = {s:i for i,s in enumerate(ids_small)}
    idx_G2 = np.array([lookup.get(s,-1) for s in df['Sublinex'].astype(str)], dtype=int)
    
    return G1_small, G2_small, idx_G1, idx_G2

# ========= PARAMETERS ==========
DATA_PATH  = "/blue/juannanzhou/fahimeh.rahimi/O2_data/Data/clean_df_with_G0.csv"
OUTPUT_DIR = "/blue/juannanzhou/fahimeh.rahimi/O2_data/Results"
TUNE, DRAWS = 2000, 20000
N_BOOTSTRAPS = 8000  # Increased to compensate for filtering


FALSE_POSITIVE_LAMBDA = 3  
MUTATION_RATE_O1 = 2.126878e-09  
MUTATION_RATE_O2 = 2.474332e-09   
FAILURE_TO_RECALL_RATE = 0.025  
GENOME_LENGTH = 1e8 
N_GENERATIONS = 150  

INTERCEPT_SD, BLOCK_EFFECT_SD = 10, 5
SIGMA_G_SD, MU_G_SD, SIGMA_E_SD = 5, 10, 5
TARGET_ACCEPT = 0.95

# ========= SETUP ==========
os.makedirs(OUTPUT_DIR, exist_ok=True)
t0 = time.time()

log(f"Loading {DATA_PATH}")
df_orig = pd.read_csv(DATA_PATH).reset_index(drop=True)
# parse Base & Rep
split = df_orig['Sublinex'].astype(str).str.split('.', expand=True)
df_orig['Base'], df_orig['Rep'] = split[0].astype(int), split[1].astype(int)
log(f"Loaded {len(df_orig):,} rows")

# ========= PRE-GENERATE ALL RESAMPLED G MATRICES ==========
log(f"Pre-generating {N_BOOTSTRAPS} resampled G matrices...")

all_G1_matrices = []
all_G2_matrices = []
all_idx_G1 = []
all_idx_G2 = []

for bootstrap_num in range(N_BOOTSTRAPS):
    if bootstrap_num % 500 == 0:
        log(f"Generated {bootstrap_num}/{N_BOOTSTRAPS} G matrices")
    
    # Resample O1 and O2 mutations for this bootstrap
    df = df_orig.copy()
    df['O1_totalmutsx'] = resample_mutations(
        df['O1_totalmutsx'],
        0,  # FALSE_POSITIVE_LAMBDA set to zero for O1
        MUTATION_RATE_O1,
        FAILURE_TO_RECALL_RATE,
        GENOME_LENGTH,
        N_GENERATIONS
    )
    df['O2_totalmutsx'] = resample_mutations(
        df['O2_totalmutsx'],
        FALSE_POSITIVE_LAMBDA,
        MUTATION_RATE_O2,
        FAILURE_TO_RECALL_RATE,
        GENOME_LENGTH,
        N_GENERATIONS
    )
    
    # Build G matrices from resampled counts
    G1_small, G2_small, idx_G1, idx_G2 = build_G_matrices(df)
    
    # Ensure G matrices have non-negative values (required for covariance matrices)
    G1_small = np.maximum(G1_small, 0)
    G2_small = np.maximum(G2_small, 0)
    
    # Skip matrices with zero diagonals (would cause singular covariance)
    if np.any(np.diag(G1_small) == 0) or np.any(np.diag(G2_small) == 0):
        continue
    
    all_G1_matrices.append(G1_small)
    all_G2_matrices.append(G2_small)
    all_idx_G1.append(idx_G1)
    all_idx_G2.append(idx_G2)

log(f"Pre-generated {len(all_G1_matrices)} G matrices")

# ========= SINGLE PYMCM MODEL WITH BOOTSTRAP LOOP ==========
log("Starting single PyMC model with bootstrap loop...")

with pm.Model() as model:
    # intercept + block
    c   = pm.Normal("c", mu=0, sigma=INTERCEPT_SD)
    b   = pm.Normal("b", mu=0, sigma=BLOCK_EFFECT_SD, shape=df_orig['Block'].nunique())
    blk = b[df_orig['Block'].values - 1]

    # genetic effects parameters
    sigma_g1 = pm.HalfNormal("sigma_g1", sigma=SIGMA_G_SD)
    mu_g1    = pm.Normal("mu_g1", mu=0, sigma=MU_G_SD)
    sigma_g2 = pm.HalfNormal("sigma_g2", sigma=SIGMA_G_SD)
    mu_g2    = pm.Normal("mu_g2", mu=0, sigma=MU_G_SD)

    # residual
    sigma_e = pm.HalfNormal("sigma_e", sigma=SIGMA_E_SD)

    # Random index to select which G matrix to use at each step
    actual_n_matrices = len(all_G1_matrices)
    bootstrap_idx = pm.DiscreteUniform("bootstrap_idx", lower=0, upper=actual_n_matrices-1)
    
    # Convert lists to numpy arrays and then to PyTensor constants
    G1_array = np.array(all_G1_matrices)
    G2_array = np.array(all_G2_matrices)
    idx_G1_array = np.array(all_idx_G1)
    idx_G2_array = np.array(all_idx_G2)
    
    # Convert to PyTensor constants
    G1_tensor = tt.constant(G1_array)
    G2_tensor = tt.constant(G2_array)
    idx_G1_tensor = tt.constant(idx_G1_array)
    idx_G2_tensor = tt.constant(idx_G2_array)
    
    # Select G matrices based on bootstrap index using PyTensor indexing
    G1_selected = G1_tensor[bootstrap_idx]
    G2_selected = G2_tensor[bootstrap_idx]
    idx_G1_selected = idx_G1_tensor[bootstrap_idx]
    idx_G2_selected = idx_G2_tensor[bootstrap_idx]
    
    # Build genetic effects using selected G matrices
    G1_c = G1_selected  # Already a tensor
    mu1 = tt.dot(G1_c, tt.ones(G1_c.shape[0]) * mu_g1)
    g1_small = pm.MvNormal(
        "g1_small",
        mu=mu1,
        cov=G1_c * sigma_g1**2,
        shape=G1_c.shape[0]
    )

    G2_c = G2_selected  # Already a tensor
    mu2 = tt.dot(G2_c, tt.ones(G2_c.shape[0]) * mu_g2)
    g2_small = pm.MvNormal(
        "g2_small",
        mu=mu2,
        cov=G2_c * sigma_g2**2,
        shape=G2_c.shape[0]
    )

    # Expand to full observations
    n_obs = len(df_orig)
    g1_full = tt.zeros(n_obs)
    g1_full = tt.set_subtensor(g1_full[idx_G1_selected >= 0], g1_small[idx_G1_selected[idx_G1_selected >= 0]])
    g2_full = tt.zeros(n_obs)
    g2_full = tt.set_subtensor(g2_full[idx_G2_selected >= 0], g2_small[idx_G2_selected[idx_G2_selected >= 0]])

    # likelihood
    mu = c + blk + g1_full + g2_full
    pm.Normal("likelihood", mu=mu, sigma=sigma_e, observed=df_orig['lnCI'].values)

    log(f"Sampling {DRAWS} draws ({TUNE} tune; target_accept={TARGET_ACCEPT})")
    trace = pm.sample(
        draws         = DRAWS,
        tune          = TUNE,
        target_accept = TARGET_ACCEPT,
        return_inferencedata=True
    )
    log("Sampling finished")

# ========= SAVE RESULTS ==========
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
trace.to_netcdf(f"{OUTPUT_DIR}/bootstrap_trace_{timestamp}.nc")
summary = az.summary(trace)
summary.to_csv(f"{OUTPUT_DIR}/bootstrap_summary_{timestamp}.csv", index=False)

log(f"All done in {(time.time()-t0)/60:.1f} minutes")
log(f"Completed single model with {N_BOOTSTRAPS} pre-generated G matrices") 