# README

This is just like v2. But the priors for tau, omega, psi are looser.

# SLURM files
- `run-single.jl` for running one experiment.
- `concat-results.jl` for concatenating results.
- `spawn-jobs.sh` SLURM job for doing all simulations.
- `concat-results.sh` SLURM job for `concat-results.jl`
- `-jobs.sh` SLURM job for doing all simulations.
- `submit.sh` for submitting to SLURM via `source submit.sh`. Runs
  `spawn-jobs.sh` and `concat-results.sh`
