# Microbiome HPC (SLURM + Snakemake)

This repository demonstrates a robust pattern for HPC automation with SLURM and Snakemake. It emphasizes reproducibility, graceful failure handling, and efficient scheduling.

## Quickstart
1. **Edit paths** in `config/snake_config.yaml` and the scripts for your environment.
2. **Create environment**: `source env_setup.sh`
3. **Run checks**: `bash preflight_checks.sh`
4. **Launch workflow**: `bash snakemake_workflow.sh`
5. (Optional) **Manual array** demo: `sbatch --array=1-$(wc -l < samples.txt) job_array.sh samples.txt`
6. **Sync outputs**: `bash rsync.sh`

## Input Mapping Methods
You have **three** ways to bind `samples.txt` IDs to raw reads:


1. **Convention (default)**
- Paths inferred from `config.yaml`:
- `raw_dir` + `<ID>` + `r1_suffix` / `r2_suffix`
- Example: `/.../raw/S01_R1.fastq.gz`, `/.../raw/S01_R2.fastq.gz`


2. **Manifest (recommended for heterogeneous inputs)**
- Create `samples.tsv` with columns: `sample`, `r1_path`, `r2_path`.
- Set `use_manifest: true` in `config.yaml`.


3. **Auto‑discovery of lanes (L001/L002/…)**
- Enable `autodiscover: true` and set `lane_glob`.
- The pipeline will **glob** per‑lane files and **merge** to per‑sample `.fastq.gz`.


## First Step
- `fastp` runs **directly on .fastq.gz**, producing trimmed `.fastq.gz` and reports.
- Swap in `kneaddata` or other tools easily—most read `.gz` natively.


## Notes
- Replace the placeholder "compute" rule with your actual tools (fastp, kneaddata, MetaPhlAn, HUMAnN, MaAsLin3, etc.). Map tool threads to `resources.cpus`.
- Keep modules and package versions pinned for reproducibility.
- Consider adding CI (e.g., `snakemake --lint`) to catch DAG issues early.