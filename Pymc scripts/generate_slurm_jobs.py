#!/usr/bin/env python

import os

# Configuration
model_script = "bootstrap_pymc_model.py"
slurm_dir = "./generated_slurms"

os.makedirs(slurm_dir, exist_ok=True)

# Your desired tuning/sample settings
tunes = [2000]
draws = [20000]  # Using your current setting

# Loop over each combination of tune & draw
for tune in tunes:
    for draw in draws:
        tag = f"tune{tune}_draws{draw}"
        
        # Filename for the SLURM submission script
        slurm_filename = f"run_{tag}.slurm"

        # Write out the corresponding SLURM submission script
        with open(os.path.join(slurm_dir, slurm_filename), "w") as f_slurm:
            f_slurm.write(f"""#!/bin/bash
#SBATCH --job-name=pmc_{tag}
#SBATCH --output=pmc_{tag}_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fahimeh.rahimi@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=2-00:00:00
#SBATCH --partition=hpg-default

source ~/.bashrc
conda activate pymc_env

# Run the bootstrap model script
python {model_script}
""")

print("✅ Generated all SLURM files.")
print(f"Generated {len(tunes) * len(draws)} SLURM scripts in {slurm_dir}/")
print("To submit jobs, use: sbatch generated_slurms/run_*.slurm") 