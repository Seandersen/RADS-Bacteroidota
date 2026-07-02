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
| Defense genes identified | 234 |
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
python -m shiny run dashboard/app.py --port 8000
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
