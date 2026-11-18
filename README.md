# Normative Modeling of Structural and Diffusion MRI for Epilepsy Diagnosis and Classification

Analysis code for Zhou et al., medRxiv 2026, [doi:10.64898/2026.11.14.26368080](https://doi.org/10.64898/2026.11.14.26368080). The study tests whether normative deviations from clinical T1-weighted and diffusion MRI distinguish epilepsy from functional seizures, focal from generalized epilepsy, and drug-resistant from drug-responsive epilepsy.

## Set up the environment

Use Python 3.12 and install the pinned packages:

```
pip install -r requirements.txt
```

For the exact package set that produced the results, install `environment/requirements-lock.txt` instead.

## Get the data

The derived data (CentileBrain z-scores, tract profiles and node-level z-scores) will be released on Pennsieve upon publication. Raw MRI is not shared. Place the input files under `data/` as listed in [data/README.md](data/README.md).

## Run the analysis

```
bash run_analysis.sh
```

This runs the five notebooks in order and writes executed copies, tables and figures to `outputs/`. It takes about five minutes on a recent laptop. The Figure 5 cell downloads atlas meshes on first use, so it needs network access once.

To run the notebooks yourself in Jupyter, open them from `analysis/` and run `supp_table_discrimination` last, since it reads tables that `composite_figures` and `abnormality_burden_summary` write.

## Check the results

Tables and figures are written to `outputs/`. Compare them with the outputs stored in the notebooks under `analysis/` and with the tables in the paper.

## Preprocessing

`preprocessing/` holds the FreeSurfer, CentileBrain input and tractometry scripts used on the Penn cluster. They and the `supplementary/` scripts read site-specific paths from `.env`: copy `.env.example` to `.env`, fill it in, and run `source env.sh` first. The notebooks do not use `.env`. [preprocessing/freesurfer/README.md](preprocessing/freesurfer/README.md) gives the FreeSurfer command and container.

## Citation

Zhou DJ, Jaskir M, Hutchinson R, et al. Normative Modeling of Structural and Diffusion MRI for Epilepsy Diagnosis and Classification. medRxiv. 2026. doi:10.64898/2026.11.14.26368080

## License

MIT, see `LICENSE`.
