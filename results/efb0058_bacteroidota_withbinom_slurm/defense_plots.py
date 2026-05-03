#!/usr/bin/env python3
"""
Defense system visualization for Bacteroidota RADS dataset.
Publication-quality figures matching the reference aesthetic:
  - White background, single dark-teal bars, no top/right spines
  - Bold A/B panel labels
  - PDF output

Figures produced:
  1_AB_categories_heatmap.pdf  — A: category bar chart; B: family × category heatmap
  2_top_two_categories.pdf     — A/B: top two categories broken down by system type
  4_cotrans_scores.pdf         — Histogram of defense scores for co-transcribed ORFs
"""

import re
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RESULTS     = Path("/Volumes/Duerkop-lab/Shelby/017_RADSBacteroidotaFinal/RADS/results/efb0058_bacteroidota_withbinom_slurm")
GENOMES_DIR = Path("/Volumes/Duerkop-lab/Shelby/017_RADSBacteroidotaFinal/RADS/RADS/results/efb0058_bacteroidota_withbinom_slurm/genomes")
DEFENSE_TSV = RESULTS / "defensefinder" / "defense_finder_genes.tsv"
SCORES_TSV  = RESULTS / "defense_scores.tsv"
OUT_DIR     = RESULTS / "figures"
OUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Publication aesthetic
# ---------------------------------------------------------------------------
BAR_COLOR  = "#3d7a7a"          # single teal for all bars
TEAL_DARK  = "#2d4a4a"          # darkest teal (heatmap max)
TEAL_CMAP  = LinearSegmentedColormap.from_list("teal_pub", ["#ffffff", TEAL_DARK])

plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica Neue", "DejaVu Sans"],
    "font.size":          8,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "xtick.major.size":   3,
    "ytick.major.size":   3,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
    "pdf.fonttype":       42,   # embed fonts as TrueType
    "ps.fonttype":        42,
})

# ---------------------------------------------------------------------------
# Defense system category mapping
# ---------------------------------------------------------------------------
_CAT_EXACT = {
    "Cas": "CRISPR-Cas", "CRISPR-Cas": "CRISPR-Cas",
    "CBASS": "CBASS", "Thoeris": "CBASS", "Pycsar": "CBASS",
    "RM": "Nucleic Acid Restriction", "BREX": "Nucleic Acid Restriction",
    "DISARM": "Nucleic Acid Restriction", "Dnd": "Nucleic Acid Restriction",
    "Dpd": "Nucleic Acid Restriction", "Wadjet": "Nucleic Acid Restriction",
    "Zorya": "Nucleic Acid Restriction", "Shedu": "Nucleic Acid Restriction",
    "NixI": "Nucleic Acid Restriction", "Nhi": "Nucleic Acid Restriction",
    "pAgo": "Nucleic Acid Restriction", "RADAR": "Nucleic Acid Restriction",
    "SspBCDE": "Nucleic Acid Restriction",
    "PrrC": "tRNA Degradation", "RloC": "tRNA Degradation", "CapRel": "tRNA Degradation",
    "Retron": "Retron",
    "DRT": "Toxin/Antitoxin", "MazEF": "Toxin/Antitoxin", "RexAB": "Toxin/Antitoxin",
    "RnlAB": "Toxin/Antitoxin", "SoFIC": "Toxin/Antitoxin",
    "RosmerTA": "Toxin/Antitoxin", "ShosTA": "Toxin/Antitoxin",
    "SanaTA": "Toxin/Antitoxin",
    "Gabija": "Abortive Infection", "Druantia": "Abortive Infection",
    "Hachiman": "Abortive Infection", "Lamassu-Fam": "Abortive Infection",
    "Lamassu": "Abortive Infection", "PARIS": "Abortive Infection",
    "Paris": "Abortive Infection", "Avs": "Abortive Infection",
    "BstA": "Abortive Infection", "Kiwa": "Abortive Infection",
    "Lit": "Abortive Infection", "Shango": "Abortive Infection",
    "JukAB": "Abortive Infection", "Septu": "Abortive Infection",
    "SEFIR": "Abortive Infection", "GasderMIN": "Abortive Infection",
    "SpbK": "Abortive Infection", "Stk2": "Abortive Infection",
    "Pif": "Abortive Infection", "DdmDE": "Abortive Infection",
    "Dsr": "Abortive Infection", "Viperin": "Abortive Infection",
    "DarTG": "Abortive Infection", "Detocs": "Abortive Infection",
    "Borvo": "Abortive Infection", "Abi": "Abortive Infection",
    "AbiD": "Abortive Infection", "AbiH": "Abortive Infection",
    "AbiZ": "Abortive Infection", "AbiAlpha": "Abortive Infection",
    "Menshen": "Unknown Mechanism", "Mokosh": "Unknown Mechanism",
    "Aditi": "Unknown Mechanism", "Dazbog": "Unknown Mechanism",
    "Tiamat": "Unknown Mechanism", "Dodola": "Unknown Mechanism",
    "Eleos": "Unknown Mechanism", "NLR": "Unknown Mechanism",
    "Azaca": "Unknown Mechanism", "Bunzi": "Unknown Mechanism",
    "Uzume": "Unknown Mechanism", "ISG15-like": "Unknown Mechanism",
    "MADS": "Unknown Mechanism", "Veles": "Unknown Mechanism",
    "DS-1": "Unknown Mechanism",
    "Anti_RM": "Anti-defense", "ardc": "Anti-defense",
    "arda_ardu": "Anti-defense", "ADF": "Anti-defense",
}
_CAT_PREFIX = [
    ("Cas_Type_",  "CRISPR-Cas"),
    ("CBASS_Type_","CBASS"),
    ("RM_Type_",   "Nucleic Acid Restriction"),
    ("Retron_",    "Retron"),
    ("Abi",        "Abortive Infection"),
    ("Gao_",       "Unknown Mechanism"),
    ("PD-Lambda-", "Unknown Mechanism"),
    ("PD-T7-",     "Unknown Mechanism"),
    ("PD-T4-",     "Abortive Infection"),
    ("FS_",        "Nucleic Acid Restriction"),
    ("Lamassu",    "Abortive Infection"),
]

# Desired display order for categories (matches reference figure)
CAT_ORDER = [
    "Toxin/Antitoxin",
    "Nucleic Acid Restriction",
    "Abortive Infection",
    "CBASS",
    "Unknown Mechanism",
    "Retron",
    "tRNA Degradation",
    "CRISPR-Cas",
    "Anti-defense",
]


def category(system_type: str, subtype: str = "") -> str:
    for field in (system_type, subtype):
        if not field or field == "Unknown":
            continue
        if field in _CAT_EXACT:
            return _CAT_EXACT[field]
        for prefix, cat in _CAT_PREFIX:
            if field.startswith(prefix):
                return cat
    return "Unknown Mechanism"


# ---------------------------------------------------------------------------
# Bacteroidota genus → (family, order)
# ---------------------------------------------------------------------------
GENUS_TO_FAMILY = {
    # Bacteroidales
    "Bacteroides":    ("Bacteroidaceae",    "Bacteroidales"),
    "Phocaeicola":    ("Bacteroidaceae",    "Bacteroidales"),
    "Prevotella":     ("Prevotellaceae",    "Bacteroidales"),
    "Hallella":       ("Prevotellaceae",    "Bacteroidales"),
    "Hoylesella":     ("Prevotellaceae",    "Bacteroidales"),
    "Segatella":      ("Prevotellaceae",    "Bacteroidales"),
    "Porphyromonas":  ("Porphyromonadaceae","Bacteroidales"),
    "Parabacteroides":("Tannerellaceae",    "Bacteroidales"),
    "Tannerella":     ("Tannerellaceae",    "Bacteroidales"),
    "Macellibacteroides":("Tannerellaceae", "Bacteroidales"),
    "Alistipes":      ("Rikenellaceae",     "Bacteroidales"),
    "Rikenella":      ("Rikenellaceae",     "Bacteroidales"),
    "Dysgonomonas":   ("Dysgonomonadaceae", "Bacteroidales"),
    "Indibacter":     ("Dysgonomonadaceae", "Bacteroidales"),
    "Barnesiella":    ("Barnesiellaceae",   "Bacteroidales"),
    "Coprobacter":    ("Barnesiellaceae",   "Bacteroidales"),
    "Odoribacter":    ("Odoribacteraceae",  "Bacteroidales"),
    "Butyricimonas":  ("Odoribacteraceae",  "Bacteroidales"),
    "Blattabacterium":("Blattabacteriaceae","Bacteroidales"),
    "Muribaculum":    ("Muribaculaceae",    "Bacteroidales"),
    "Duncaniella":    ("Muribaculaceae",    "Bacteroidales"),
    "Marinifilum":    ("Marinifilaceae",    "Bacteroidales"),
    "Prolixibacter":  ("Marinifilaceae",    "Bacteroidales"),
    "Alkaliflexus":   ("Marinifilaceae",    "Bacteroidales"),
    "Candidatus":     ("Candidatus",        "Bacteroidales"),
    # Flavobacteriales
    "Flavobacterium": ("Flavobacteriaceae", "Flavobacteriales"),
    "Tenacibaculum":  ("Flavobacteriaceae", "Flavobacteriales"),
    "Polaribacter":   ("Flavobacteriaceae", "Flavobacteriales"),
    "Maribacter":     ("Flavobacteriaceae", "Flavobacteriales"),
    "Leeuwenhoekiella":("Flavobacteriaceae","Flavobacteriales"),
    "Capnocytophaga": ("Flavobacteriaceae", "Flavobacteriales"),
    "Zobellia":       ("Flavobacteriaceae", "Flavobacteriales"),
    "Cellulophaga":   ("Flavobacteriaceae", "Flavobacteriales"),
    "Formosa":        ("Flavobacteriaceae", "Flavobacteriales"),
    "Psychroflexus":  ("Flavobacteriaceae", "Flavobacteriales"),
    "Gillisia":       ("Flavobacteriaceae", "Flavobacteriales"),
    "Lacinutrix":     ("Flavobacteriaceae", "Flavobacteriales"),
    "Muricauda":      ("Flavobacteriaceae", "Flavobacteriales"),
    "Gelidibacter":   ("Flavobacteriaceae", "Flavobacteriales"),
    "Wautersiella":   ("Flavobacteriaceae", "Flavobacteriales"),
    "Psychroserpens": ("Flavobacteriaceae", "Flavobacteriales"),
    "Croceibacter":   ("Flavobacteriaceae", "Flavobacteriales"),
    "Ulvibacter":     ("Flavobacteriaceae", "Flavobacteriales"),
    "Daejeonia":      ("Flavobacteriaceae", "Flavobacteriales"),
    "Myroides":       ("Flavobacteriaceae", "Flavobacteriales"),
    "Elizabethkingia":("Weeksellaceae",     "Flavobacteriales"),
    "Chryseobacterium":("Weeksellaceae",    "Flavobacteriales"),
    "Cloacibacterium":("Weeksellaceae",     "Flavobacteriales"),
    "Bergeyella":     ("Weeksellaceae",     "Flavobacteriales"),
    "Ornithobacterium":("Weeksellaceae",    "Flavobacteriales"),
    "Coenonia":       ("Weeksellaceae",     "Flavobacteriales"),
    "Empedobacter":   ("Weeksellaceae",     "Flavobacteriales"),
    "Weeksella":      ("Weeksellaceae",     "Flavobacteriales"),
    "Vaginella":      ("Weeksellaceae",     "Flavobacteriales"),
    "Riemerella":     ("Weeksellaceae",      "Flavobacteriales"),
    # Sphingobacteriales
    "Sphingobacterium":("Sphingobacteriaceae","Sphingobacteriales"),
    "Pedobacter":     ("Sphingobacteriaceae","Sphingobacteriales"),
    "Olivibacter":    ("Sphingobacteriaceae","Sphingobacteriales"),
    # Chitinophagales
    "Chitinophaga":   ("Chitinophagaceae",  "Chitinophagales"),
    "Niastella":      ("Chitinophagaceae",  "Chitinophagales"),
    "Terrimonas":     ("Chitinophagaceae",  "Chitinophagales"),
    "Niabella":       ("Chitinophagaceae",  "Chitinophagales"),
    "Filimonas":      ("Chitinophagaceae",  "Chitinophagales"),
    "Flavisolibacter":("Chitinophagaceae",  "Chitinophagales"),
    "Lacibacter":     ("Chitinophagaceae",  "Chitinophagales"),
    "Ferruginibacter":("Chitinophagaceae",  "Chitinophagales"),
    "Arachidicoccus": ("Chitinophagaceae",  "Chitinophagales"),
    "Sediminibacterium":("Chitinophagaceae","Chitinophagales"),
    "Gracilimonas":   ("Chitinophagaceae",  "Chitinophagales"),
    "Mucilaginibacter":("Chitinophagaceae", "Chitinophagales"),
    "Rhizosphaericola":("Chitinophagaceae", "Chitinophagales"),
    "Haliscomenobacter":("Saprospiraceae",  "Saprospirales"),
    "Saprospira":     ("Saprospiraceae",    "Saprospirales"),
    "Lewinella":      ("Saprospiraceae",    "Saprospirales"),
    # Cytophagales
    "Catalinimonas":  ("Catalimonadaceae",  "Cytophagales"),
    "Cytophaga":      ("Cytophagaceae",     "Cytophagales"),
    "Hymenobacter":   ("Hymenobacteraceae", "Cytophagales"),
    "Pontibacter":    ("Cytophagaceae",     "Cytophagales"),
    "Runella":        ("Cytophagaceae",     "Cytophagales"),
    "Flectobacillus": ("Cytophagaceae",     "Cytophagales"),
    "Larkinella":     ("Cytophagaceae",     "Cytophagales"),
    "Spirosoma":      ("Cytophagaceae",     "Cytophagales"),
    "Cyclobacterium": ("Cyclobacteriaceae", "Cytophagales"),
    "Algoriphagus":   ("Cyclobacteriaceae", "Cytophagales"),
    "Aquiflexum":     ("Cyclobacteriaceae", "Cytophagales"),
    "Porifericola":   ("Cyclobacteriaceae", "Cytophagales"),
    # MAGs
    "MAG":     ("Unclassified", "Unclassified"),
    "Unknown": ("Unclassified", "Unclassified"),
}

# Family row order for heatmap: Unclassified first, then by order
FAMILY_ORDER = [
    # Unclassified
    "Unclassified",
    # Sphingobacteriales
    "Sphingobacteriaceae",
    # Chitinophagales
    "Chitinophagaceae",
    # Saprospirales
    "Saprospiraceae",
    # Flavobacteriales
    "Weeksellaceae",
    "Flavobacteriaceae",

    # Cytophagales
    "Cytophagaceae",
    "Cyclobacteriaceae",
    "Hymenobacteraceae",
    "Catalimonadaceae",
    # Bacteroidales
    "Bacteroidaceae",
    "Prevotellaceae",
    "Porphyromonadaceae",
    "Tannerellaceae",
    "Rikenellaceae",
    "Odoribacteraceae",
    "Muribaculaceae",
    "Dysgonomonadaceae",
    "Barnesiellaceae",
    "Blattabacteriaceae",
    "Marinifilaceae",
    "Candidatus",
]

# Order → display label and which families it spans
ORDER_GROUPS = [
    ("Sphingo-\nbacteriales", ["Sphingobacteriaceae"]),
    ("Chitino-\nphagales",    ["Chitinophagaceae"]),
    ("Saprospirales",         ["Saprospiraceae"]),
    ("Flavobac-\nteriales",   ["Weeksellaceae", "Flavobacteriaceae"]),
    ("Cytophagales",          ["Cytophagaceae", "Cyclobacteriaceae", "Hymenobacteraceae",
                                "Catalimonadaceae"]),
    ("Bacteroidales",         ["Bacteroidaceae", "Prevotellaceae", "Porphyromonadaceae",
                                "Tannerellaceae", "Rikenellaceae", "Odoribacteraceae",
                                "Muribaculaceae", "Dysgonomonadaceae", "Barnesiellaceae",
                                "Blattabacteriaceae", "Marinifilaceae", "Candidatus"]),
]


def genus_to_family_order(genus: str) -> tuple[str, str]:
    return GENUS_TO_FAMILY.get(genus, ("Other Bacteroidota", "Unclassified"))


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_defense_systems() -> pd.DataFrame:
    """Load unique defense genes from defense_finder_genes.tsv.

    The genes TSV can have repeated header lines (one per genome processed),
    so we skip any line where hit_id == 'hit_id'.  We deduplicate by hit_id
    so each unique gene is counted once.
    """
    rows = []
    seen_hits = set()
    with open(DEFENSE_TSV) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("replicon\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            hit_id = parts[1] if len(parts) > 1 else ""
            if hit_id == "hit_id":          # repeated header line
                continue
            if hit_id in seen_hits:
                continue
            seen_hits.add(hit_id)
            model_fqn = parts[4] if len(parts) > 4 else ""
            fqn_parts = model_fqn.split("/")
            subtype     = fqn_parts[-1] if fqn_parts else "Unknown"
            system_type = fqn_parts[-2] if len(fqn_parts) >= 2 else subtype
            sys_id = parts[5] if len(parts) > 5 else ""
            m = re.match(r"([A-Z]{2}_?[A-Z0-9]+\.\d+)_", hit_id)
            chrom_acc = m.group(1) if m else ""
            rows.append({
                "hit_id":      hit_id,
                "sys_id":      sys_id,
                "system_type": system_type,
                "subtype":     subtype,
                "category":    category(system_type, subtype),
                "chrom_acc":   chrom_acc,
            })
    return pd.DataFrame(rows)


def load_interpro_typeIV(existing_hit_ids: set) -> pd.DataFrame:
    """Return InterPro-derived putative Type IV restriction enzyme ORFs.

    Includes:
      - DUF262 (PF03235 / IPR004919)  — putative Type IV restriction enzyme
      - Type IV Mrr (IPR007560)        — Type IV restriction endonuclease

    ORFs already present in the DefenseFinder genes TSV are excluded to
    prevent double-counting.
    """
    TYPE4_PFAM = {"PF03235"}          # DUF262
    TYPE4_IPR  = {"IPR007560"}        # Restriction endonuclease type IV, Mrr

    rows = []
    seen = set()
    with open(RESULTS / "interproscan_results.tsv") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            orf_id   = parts[0]
            pfam_acc = parts[4]
            ipr_acc  = parts[11]
            match_ok = parts[9] == "T"
            if not match_ok:
                continue
            if pfam_acc not in TYPE4_PFAM and ipr_acc not in TYPE4_IPR:
                continue
            if orf_id in existing_hit_ids or orf_id in seen:
                continue
            seen.add(orf_id)
            m = re.match(r"([A-Z]{2}_?[A-Z0-9]+\.\d+)_", orf_id)
            chrom_acc = m.group(1) if m else ""
            if pfam_acc in TYPE4_PFAM:
                label = "DUF262 (putative Type IV)"
            else:
                label = "Type IV Mrr"
            rows.append({
                "hit_id":      orf_id,
                "sys_id":      f"interpro_{label.replace(' ', '_')}",
                "system_type": "RM",
                "subtype":     label,
                "category":    "Nucleic Acid Restriction",
                "chrom_acc":   chrom_acc,
            })
    return pd.DataFrame(rows)


def load_acc_to_genus(cache: Path = Path("/tmp/genus_cache.json")) -> dict[str, str]:
    if cache.exists():
        with open(cache) as f:
            return json.load(f)
    acc_to_genus: dict[str, str] = {}
    for fna in GENOMES_DIR.glob("*.fna"):
        try:
            with open(fna, encoding="utf-8", errors="replace") as fh:
                first = fh.readline()
            if not first.startswith(">"):
                continue
            parts = first[1:].split(None, 1)
            acc = parts[0]
            genus = "Unknown"
            if len(parts) >= 2:
                words = parts[1].strip().split()
                genus = re.sub(r"[^A-Za-z]", "", words[0]) if words else "Unknown"
            acc_to_genus[acc] = genus
        except Exception:
            pass
    with open(cache, "w") as f:
        json.dump(acc_to_genus, f)
    return acc_to_genus


def load_defense_scores() -> pd.DataFrame:
    return pd.read_csv(SCORES_TSV, sep="\t")


# ---------------------------------------------------------------------------
# Helper: add bold panel label (A, B, …)
# ---------------------------------------------------------------------------
def panel_label(ax, letter: str, fontsize: int = 12, x: float = -0.12, y: float = 1.05):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------------------
# Helper: draw order-group bracket annotations on heatmap y-axis
# ---------------------------------------------------------------------------
def add_order_brackets(ax, families_present: list[str], x_offset: float = -0.38):
    """Draw vertical bracket + rotated order label to the left of the heatmap y-axis."""
    n = len(families_present)
    for label, fam_list in ORDER_GROUPS:
        # Which rows in families_present are in this order group?
        idxs = [i for i, f in enumerate(families_present) if f in fam_list]
        if not idxs:
            continue
        # In heatmap coordinates, row 0 is top → y positions are reversed
        # seaborn heatmap: row i has center at i + 0.5 in data coords
        y_top = n - min(idxs) - 0.5
        y_bot = n - max(idxs) - 0.5
        y_mid = (y_top + y_bot) / 2

        # Bracket line
        ax.annotate("", xy=(x_offset + 0.04, y_bot), xytext=(x_offset + 0.04, y_top),
                    xycoords=("axes fraction", "data"),
                    textcoords=("axes fraction", "data"),
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color="black", lw=0.8,
                                    connectionstyle="arc3,rad=0"))
        # Short horizontal ticks at top and bottom of bracket
        for y_tick in (y_top, y_bot):
            ax.annotate("", xy=(x_offset + 0.04, y_tick),
                        xytext=(x_offset + 0.08, y_tick),
                        xycoords=("axes fraction", "data"),
                        textcoords=("axes fraction", "data"),
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", color="black", lw=0.8))
        # Rotated label
        ax.text(x_offset, y_mid, label,
                transform=ax.get_yaxis_transform(),
                rotation=90, ha="center", va="center",
                fontsize=7, fontweight="bold")


# ===========================================================================
# FIGURE 1 — A/B: Categories bar chart + Family heatmap
# ===========================================================================
def make_figure_AB(df_sys: pd.DataFrame, acc_to_genus: dict,
                   out_path: Path):

    # ---- Data preparation ----
    # Category counts sorted most → least abundant
    cat_counts = df_sys["category"].value_counts()
    cats_present = list(cat_counts.index)   # already most→least
    counts = [cat_counts[c] for c in cats_present]

    # Heatmap: family × category
    df = df_sys.copy()
    df["genus"]  = df["chrom_acc"].map(acc_to_genus).fillna("Unknown")
    df["family"] = df["genus"].map(lambda g: genus_to_family_order(g)[0])
    pivot = df.pivot_table(index="family", columns="category",
                           values="hit_id", aggfunc="count", fill_value=0)
    # Reorder rows to match FAMILY_ORDER (only families present in data)
    fam_order = [f for f in FAMILY_ORDER if f in pivot.index]
    fam_order += [f for f in pivot.index if f not in fam_order]
    pivot = pivot.reindex(fam_order).fillna(0)
    # Reorder columns most → least abundant (matching bar chart)
    col_order = [c for c in cats_present if c in pivot.columns]
    col_order += [c for c in pivot.columns if c not in col_order]
    pivot = pivot[col_order]
    # Drop zero rows/cols
    pivot = pivot.loc[pivot.sum(axis=1) > 0, pivot.sum(axis=0) > 0]
    families_present = list(pivot.index)

    # ---- Layout ----
    n_fam = len(families_present)
    n_cat = len(pivot.columns)
    fig_w = 4.5 + n_cat * 0.55
    fig_h = max(3.5, n_fam * 0.45 + 1.5)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = gridspec.GridSpec(1, 2,
                            width_ratios=[1.0, max(1.5, n_cat * 0.38)],
                            wspace=0.5)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])

    # ---- Panel A: vertical bar chart ----
    x_pos = np.arange(len(cats_present))
    ax_a.bar(x_pos, counts, color=BAR_COLOR, width=0.65, zorder=3)
    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels(cats_present, rotation=45, ha="right", fontsize=7)
    ax_a.set_ylabel("Number of Instances", fontsize=8)
    ax_a.set_ylim(0, max(counts) * 1.15)
    ax_a.yaxis.set_tick_params(labelsize=8)
    ax_a.spines["left"].set_linewidth(0.8)
    ax_a.spines["bottom"].set_linewidth(0.8)
    panel_label(ax_a, "A")

    # ---- Panel B: heatmap ----
    # Font size for annotations (counts)
    annot_kws = {"size": 6.5}
    # Replace zeros with empty string for cleaner display
    annot_data = pivot.astype(int).astype(str).replace("0", "")

    sns.heatmap(
        pivot,
        ax=ax_b,
        cmap=TEAL_CMAP,
        linewidths=0.4,
        linecolor="#e0e0e0",
        annot=annot_data,
        fmt="",
        annot_kws=annot_kws,
        cbar_kws={"label": "", "shrink": 0.5, "pad": 0.02},
        yticklabels=True,
        xticklabels=True,
    )

    # Style heatmap axes
    ax_b.set_xlabel("", fontsize=0)
    ax_b.set_ylabel("", fontsize=0)
    ax_b.set_xticklabels(ax_b.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    ax_b.set_yticklabels(ax_b.get_yticklabels(), rotation=0, fontsize=7)
    ax_b.tick_params(left=False, bottom=False)

    # Colorbar label
    cbar = ax_b.collections[0].colorbar
    cbar.set_label("Count", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # Order group brackets
    add_order_brackets(ax_b, families_present)

    panel_label(ax_b, "B", x=-0.55)

    plt.savefig(out_path, dpi=300, bbox_inches="tight", format="pdf")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ===========================================================================
# FIGURE 2 — A/B: Top two categories broken down by system type
# ===========================================================================
def make_figure_top_two(df_sys: pd.DataFrame, out_path: Path):

    cat_counts = df_sys["category"].value_counts()
    top2 = cat_counts.index[:2].tolist()

    fig = plt.figure(figsize=(8, 4))
    gs  = gridspec.GridSpec(1, 2, wspace=0.45)
    letters = ["A", "B"]

    for idx, (cat, letter) in enumerate(zip(top2, letters)):
        ax = fig.add_subplot(gs[idx])
        subset = df_sys[df_sys["category"] == cat].copy()
        type_counts = subset["subtype"].value_counts()

        x_pos = np.arange(len(type_counts))
        ax.bar(x_pos, type_counts.values, color=BAR_COLOR, width=0.65, zorder=3)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(type_counts.index, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Number of Instances", fontsize=8)
        ax.set_title(cat, fontsize=9, fontweight="bold", pad=6)
        ax.set_ylim(0, type_counts.max() * 1.18)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.spines["left"].set_linewidth(0.8)
        ax.spines["bottom"].set_linewidth(0.8)
        panel_label(ax, letter)

    plt.savefig(out_path, dpi=300, bbox_inches="tight", format="pdf")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ===========================================================================
# FIGURE 4 — Histogram: defense scores for co-transcribed ORFs
# ===========================================================================
def make_figure_scores(df_scores: pd.DataFrame, out_path: Path):

    n    = len(df_scores)
    high = df_scores[df_scores["defense_score"] >= 0.5]
    low  = df_scores[df_scores["defense_score"] <  0.5]
    bins = np.linspace(0, 1, 21)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hist(low["defense_score"],  bins=bins, color=BAR_COLOR,
            edgecolor="white", linewidth=0.5, alpha=0.5,
            label=f"Score < 0.5  (n={len(low)})")
    ax.hist(high["defense_score"], bins=bins, color=TEAL_DARK,
            edgecolor="white", linewidth=0.5, alpha=0.9,
            label=f"Score ≥ 0.5  (n={len(high)})")

    med = df_scores["defense_score"].median()
    ax.axvline(med, color="#555555", linestyle="--", linewidth=1.0)
    ax.text(med + 0.02, ax.get_ylim()[1] * 0.9,
            f"median = {med:.2f}", fontsize=7, color="#555555")

    ax.set_xlabel("Defense Score", fontsize=8)
    ax.set_ylabel("Count", fontsize=8)
    ax.set_title(f"Defense Scores — Co-transcribed ORFs (n={n})",
                 fontsize=9, fontweight="bold", pad=6)
    ax.set_xlim(-0.02, 1.02)
    ax.yaxis.set_tick_params(labelsize=8)
    ax.xaxis.set_tick_params(labelsize=8)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.legend(fontsize=7, frameon=False)
    panel_label(ax, "A")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight", format="pdf")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ===========================================================================
# Main
# ===========================================================================
if __name__ == "__main__":
    print("Loading defense systems...")
    df_sys = load_defense_systems()
    df_typeIV = load_interpro_typeIV(existing_hit_ids=set(df_sys["hit_id"]))
    if not df_typeIV.empty:
        print(f"  Adding {len(df_typeIV)} InterPro Type IV hits "
              f"({df_typeIV['subtype'].value_counts().to_dict()})")
        df_sys = pd.concat([df_sys, df_typeIV], ignore_index=True)
    print(f"  {len(df_sys)} unique genes | "
          f"categories: {df_sys['category'].value_counts().to_dict()}")

    print("Loading genus map (from cache)...")
    acc_to_genus = load_acc_to_genus()
    print(f"  {sum(df_sys['chrom_acc'].map(acc_to_genus).notna())} / "
          f"{len(df_sys)} systems with genus")

    print("Loading defense scores...")
    df_scores = load_defense_scores()
    print(f"  {len(df_scores)} co-transcribed ORFs")

    print("\nGenerating PDFs...")
    make_figure_AB(df_sys, acc_to_genus,
                   OUT_DIR / "1_AB_categories_heatmap.pdf")
    make_figure_top_two(df_sys,
                        OUT_DIR / "2_top_two_categories.pdf")
    make_figure_scores(df_scores,
                       OUT_DIR / "4_cotrans_scores.pdf")
    print(f"\nAll figures saved to: {OUT_DIR}")
