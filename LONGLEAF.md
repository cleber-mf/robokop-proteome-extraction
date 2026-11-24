# Running the all-protein ROBOKOP export on Longleaf

## What this job requests

`scripts/submit_robokop_all.sbatch` is intentionally CPU-only.  It requests
one exclusive node and `--mem=0`, which asks Slurm for all memory available on
that node.  It tries `hov,general` in that order, so Slurm can place the job in
either partition.  If your account does not have access to `hov`, change the
directive to `#SBATCH --partition=general`.

This reservation may queue for a long time.  The extractor is mostly limited
by the remote ROBOKOP API and disk writes, rather than CPU or GPU speed.  A GPU
is neither requested nor used.

Before submitting, inspect the live partition policy on Longleaf:

```bash
scontrol show partition hov
scontrol show partition general
```

Adjust the partition, time limit, or memory request if Longleaf reports an
account, QOS, or resource-limit error.

## 1. Put the repository and output in suitable storage

Clone the repository into a project or work filesystem.  Set `ROBO_RUN_DIR` to
another high-capacity directory on that filesystem; do not use your home
directory for the all-protein output.  The directory holds the SQLite state,
raw CSV, split CSVs, and temporary SQLite files.

```bash
git clone https://github.com/cleber-mf/robokop-proteome-extraction.git
cd robokop-proteome-extraction
export ROBO_RUN_DIR=/work/your_group/your_user/robokop_all
```

## 2. Create the Conda environment once

The environment has no CUDA or pip dependencies.  It includes Python 3.11,
SQLite, TLS certificates, and OpenSSL; the extraction code otherwise uses the
standard library.

```bash
bash scripts/create_conda_env.sh
```

If neither `miniforge` nor `anaconda` is a current Longleaf module, run
`module avail` and load the supported Conda module before rerunning the command.
You may place the environment outside the repository by setting
`ROBO_CONDA_ENV=/path/to/robokop-proteome`.

## 3. Submit and monitor

```bash
ROBO_RUN_DIR=/work/your_group/your_user/robokop_all \
  sbatch scripts/submit_robokop_all.sbatch

squeue -u "$USER"
```

If the job ends early or is requeued, submit the same command again.  The
state database prevents completed protein batches from being queried again.
After a completed job, inspect resource use before asking for another full
node:

```bash
seff <job_id>
```
