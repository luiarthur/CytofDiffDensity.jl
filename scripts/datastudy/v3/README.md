# README

This is just like v2. But the priors for tau, omega, psi are looser.

# SLURM files
- `run-single.jl` for running one experiment.
- `concat-results.jl` for concatenating results.
- `slurm-spawn-jobs.sh` SLURM job for doing all simulations.
- `slurm-concat-results.sh` SLURM job for `concat-results.jl`
- `slurm-submit.sh` for submitting to SLURM via `source slurm-submit.sh`. Runs
  `slurm-spawn-jobs.sh` and `slurm-concat-results.sh`
