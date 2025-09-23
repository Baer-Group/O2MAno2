# O2 Mutation Accumulation Analysis

This repository contains the Bayesian mixed-effects modeling analysis for quantifying how mutations accumulated across sequential mutation accumulation (MA) experiments contribute to competitive fitness in *C. elegans*.

## Project Overview

We compared one-stage and two-stage Bayesian models to understand the fitness effects of first-order (O1) and second-order (O2) mutations, incorporating a bootstrapping approach to account for measurement errors in mutation counts.

## Repository Structure

```
O2_data/
├── Data/
│   ├── clean_df_with_G0.csv          # Main dataset with fitness measurements
│   ├── G_matrix.csv                  # Genetic covariance matrix
│   └── G_matrix_indices.csv          # Matrix indexing information
├── Results/                          # Model outputs and figures
├── scripts/                          # Analysis scripts
├── bootstrap_pymc_model.py           # Main bootstrapped two-stage model
├── plot_posteriors_3panel.py        # Generates publication figure
├── generate_slurm_jobs.py           # SLURM job generator
└── methods_results.tex              # LaTeX manuscript
```

## Analysis Workflow

### 1. Data Preparation
- **Input**: Raw mutation accumulation data with fitness measurements (lnCI)
- **Processing**: Data cleaning, G0 handling, block effect assignment
- **Output**: `clean_df_with_G0.csv`

### 2. Model Development

#### One-Stage Model
- **File**: Single G matrix approach (referenced in manuscript)
- **Approach**: Unified genetic covariance matrix combining O1 and O2 mutations
- **Parameters**: Single set of genetic parameters (μ_g, σ_g)

#### Two-Stage Model  
- **File**: `bootstrap_pymc_model.py`
- **Approach**: Separate G1 and G2 matrices for O1 and O2 mutations
- **Parameters**: Independent genetic parameters (μ_g1, σ_g1, μ_g2, σ_g2)

### 3. Bootstrapping Approach

The two-stage model incorporates measurement error through bootstrapping:

**False Positives (O2 only):**
- Distribution: Poisson(λ = 3)
- Applied to: O2 mutations only
- Process: Subtract from observed counts

**False Negatives (O1 and O2):**
- Distribution: Poisson(λ_missed)
- λ_missed = μ_mutation × n_generations × genome_length × failure_rate
- Mutation rates: 2.126878×10⁻⁹ (O1), 2.474332×10⁻⁹ (O2)
- Generations: 150, Genome length: 10⁸ bp, Failure rate: 0.025
- Process: Add to observed counts

**Implementation:**
- Pre-generates 8,000 resampled G matrices
- Filters out matrices with zero diagonals
- Uses PyTensor indexing for dynamic matrix selection during MCMC

### 4. Model Fitting

**Framework**: PyMC with NUTS sampler
**Priors**: 
- Intercepts and genetic means: Normal(0, 10)
- Standard deviations: HalfNormal(5)
**Sampling**: 2,000 tuning + 20,000 posterior draws
**Target acceptance**: 0.95

### 5. Model Selection

**Metrics**: WAIC and LOO cross-validation
**Results**:
- Two-stage model: WAIC = -1378.04, LOO = -1378.93
- One-stage model: WAIC = -1377.09, LOO = -1377.63
- Conclusion: Two-stage model preferred (lower is better)

## Key Scripts

### `bootstrap_pymc_model.py`
Main analysis script implementing the two-stage model with bootstrapping.

**Key Functions:**
- `resample_mutations()`: Handles false positive/negative simulation
- `build_G_matrices()`: Constructs G1 and G2 matrices from resampled data
- PyMC model with dynamic G matrix selection

**Usage:**
```bash
python bootstrap_pymc_model.py
```

### `plot_posteriors_3panel.py`
Generates the main publication figure showing posterior distributions.

**Outputs:**
- Panel A: Mean genetic effects (μ) for O1 vs O2
- Panel B: Standard deviations (σ) for O1 vs O2  
- Panel C: Genetic variance ratios (σ²_g / σ²_e) for O1 vs O2

**Usage:**
```bash
python plot_posteriors_3panel.py
```

### `generate_slurm_jobs.py`
Creates SLURM job submission scripts for HPC execution.

**Configuration:**
- 4 CPU cores, 16GB memory, 2-day time limit
- Conda environment: pymc_env
- Email notifications included

**Usage:**
```bash
python generate_slurm_jobs.py
sbatch generated_slurms/run_tune2000_draws20000.slurm
```

## Results

### Model Comparison
The two-stage model showed superior performance across both WAIC and LOO metrics, with small but consistent advantages over the one-stage approach.

### Posterior Distributions
- **O1 mutations**: Negative mean effects (-0.02), lower variability
- **O2 mutations**: Near-zero mean effects (0.00), higher variability
- **Genetic variance ratios**: O1 concentrated near zero, O2 with broader spread

### Key Findings
- O1 and O2 mutations show distinct fitness effect patterns
- Two-stage modeling captures biological relevance of separating mutation classes
- Bootstrapping approach successfully propagates measurement uncertainty

## Dependencies

**Python packages:**
- PyMC (≥5.0)
- ArviZ
- NumPy
- Pandas
- Matplotlib
- PyTensor

**System requirements:**
- Python 3.8+
- SLURM (for HPC execution)
- Conda environment management

## File Outputs

**Model traces:**
- `Results/bootstrap_trace_YYYYMMDD_HHMMSS.nc`
- `Results/bootstrap_summary_YYYYMMDD_HHMMSS.csv`

**Figures:**
- Posterior distribution plots (3-panel publication figure)
- Trace plots for diagnostics

**Manuscripts:**
- `methods_results.tex`: LaTeX source for Methods and Results sections

## Citation

If you use this analysis approach, please cite:
[Your publication details here]

## Contact

[Your contact information]

---

*Analysis conducted using Bayesian mixed-effects modeling with PyMC for mutation accumulation fitness effects in C. elegans.* 