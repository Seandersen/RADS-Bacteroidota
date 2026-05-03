# RADS — Bacteroidota Defense System Analysis

Analysis of defense systems flanking EFB_0058 (recombinase) across **2,035 Bacteroidota genomes** using the [RADS pipeline](https://github.com/Seandersen/RADS).

RADS was developed by Shelby E Andersen in collaboration with Joshua M Kirsch, Jay R Hesselberth, and Breck A Duerkop.

## Dataset Summary

| Metric | Value |
|--------|-------|
| Phylum | Bacteroidota (Bacteroidetes) |
| Query | EFB_0058 recombinase |
| Total genomes | 2,035 |
| BLAST hits | 290 |
| Hits per Mb | 0.035 |
| Contigs analyzed | 278 |
| Defense genes identified | 525 |
| Co-transcribed ORFs | 77 |

Results are in `results/efb0058_bacteroidota_withbinom_slurm/`.

---

## Viewing the Interactive Dashboard

### Option A — Pixi (recommended)

```bash
git clone https://github.com/Seandersen/RADS-Bacteroidota.git
cd RADS-Bacteroidota
curl -fsSL https://pixi.sh/install.sh | bash   # skip if Pixi already installed
pixi install
pixi run dashboard
```

Open **http://localhost:8000** in your browser.

### Option B — Conda/Mamba

```bash
git clone https://github.com/Seandersen/RADS-Bacteroidota.git
cd RADS-Bacteroidota
conda env create -f environment.yaml
conda activate rads
python -m shiny run dashboard/app.py --host 0.0.0.0 --port 8000
```

Open **http://localhost:8000** in your browser.

### HPC / Remote server access

If running on an HPC node or remote server, use an SSH tunnel:

```bash
# On your laptop (replace <user>, <server>, <node>, <port>)
ssh -N -L 8000:<node>:8000 <user>@<server>
```

Then open **http://localhost:8000** in your browser.

Alternatively, if your cluster provides Open OnDemand, start the dashboard with `pixi run dashboard-hpc` and navigate to:

```
https://<your-ood-url>/rnode/<hostname>/8000/
```

---

## Generating the Shareable HTML Report

The self-contained HTML report requires no server and can be shared as a single file.

### Pixi

```bash
pixi run python workflow/scripts/generate_report.py \
    --results results/efb0058_bacteroidota_withbinom_slurm \
    --output  results/efb0058_bacteroidota_withbinom_slurm/report.html
```

### Conda

```bash
conda activate rads
python workflow/scripts/generate_report.py \
    --results results/efb0058_bacteroidota_withbinom_slurm \
    --output  results/efb0058_bacteroidota_withbinom_slurm/report.html
```

Open `report.html` in any browser — no server required.  
Add `--include-locus-viewer` to embed gene-arrow diagrams (produces a larger file).

---

## Publication Figures

Defense system figures are in `results/efb0058_bacteroidota_withbinom_slurm/figures/`:

| File | Content |
|------|---------|
| `1_AB_categories_heatmap.pdf` | A: Category bar chart; B: Family × category heatmap |
| `2_top_two_categories.pdf` | Top two categories broken down by system type |
| `4_cotrans_scores.pdf` | Defense score distribution for co-transcribed ORFs |

To regenerate figures:

```bash
# Pixi
pixi run python results/efb0058_bacteroidota_withbinom_slurm/defense_plots.py

# Conda
conda activate rads
python results/efb0058_bacteroidota_withbinom_slurm/defense_plots.py
```

---

## Results Structure

```
results/efb0058_bacteroidota_withbinom_slurm/
├── blast_results/master_blast.txt         # BLAST hits (query vs all genomes)
├── defensefinder/
│   ├── defense_finder_genes.tsv           # Per-gene defense annotations
│   └── defense_finder_systems.tsv         # Per-system annotations
├── cotranscription/                        # Co-transcribed ORF analysis
├── interproscan_results.tsv               # Domain annotations (InterPro/Pfam)
├── BinomialAnalysis.csv                   # Enriched domains (binomial test)
├── defense_scores.tsv                     # Defense context scores
├── metrics/pipeline_metrics.json          # Run statistics
├── figures/                               # Publication PDFs
└── defense_plots.py                       # Figure generation script
```

---

## Pipeline

This dataset was generated with the [RADS pipeline](https://github.com/Seandersen/RADS). See the pipeline repo for full documentation, including how to run your own analysis.
