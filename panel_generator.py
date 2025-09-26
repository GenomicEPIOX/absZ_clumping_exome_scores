#!/usr/bin/env python3
"""
make_panel.py

Build a per-protein variant panel from Dhindsa et al. (Supp. Table 2):
- Filters to a protein and variant type (cis/trans/combined/trans_cis_pos)
- Prunes variants in LD (within 10Mb and r² >= threshold) by keeping the SNP with the largest |Z|
- Writes:
  - <output>/<PROTEIN>/<PROTEIN>_<variant_type>_<r2>.panel  (Genotype, A1, weight)
  - <output>/<PROTEIN>/<PROTEIN>_<variant_type>_<r2>.log    (concise pruning log)

Requires: pandas, numpy, PyTables (for HDF5).
"""

import argparse
import sys
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd


# --------------------- CLI ---------------------

def check_arg(args=None):
    parser = argparse.ArgumentParser(
        description=(
            'Creates exome-score panels, selecting the variant with the largest |Z| '
            'among variants in LD (r² >= --r2) within 10Mb.'
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-p', '--protein_name', help='Protein name', required=True)
    parser.add_argument('-o', '--output',
                        help='Directory to write results (a subfolder named after the protein will be created)',
                        default='./results/')
    parser.add_argument('-ld', '--ld_matrix',
                        help='plink LD matrix in HDF5 format with keys chr1..chr23 containing columns ID_A, ID_B, PHASED_R2',
                        required=True)
    parser.add_argument('-r2', '--r2', type=float, default=0.05,
                        help='R-square threshold for pruning (default: 0.05)')
    parser.add_argument('-ex', '--excel_file',
                        help='Excel file from Dhindsa et al., Supplementary Table 2',
                        required=True)
    parser.add_argument('-mhc', '--mhc_filter', type=int, choices=[0, 1], default=0,
                        help='Set 1 to exclude variants in the classical MHC region (chr6:28,510,120–33,480,577)')
    parser.add_argument('-v', '--variant_type',
                        help='Variant set to use: cis, trans, combined, or trans_cis_pos',
                        choices=['cis', 'trans', 'combined', 'trans_cis_pos'], required=True)
    parser.add_argument('-m', '--missing',
                        help='File with variants to remove (e.g., missing from genotype set). One variant ID per line; no header.',
                        required=True)
    results = parser.parse_args(args)
    return (results.protein_name, results.output, results.ld_matrix, results.r2,
            results.excel_file, results.mhc_filter, results.variant_type, results.missing)


# --------------------- Utilities ---------------------

class UnionFind:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def filter_MHC(df):
    """Exclude classical MHC region: chr6:28,510,120–33,480,577."""
    in_mhc = (df["CHR"] == 6) & (df["BP"].between(28_510_120, 33_480_577))
    return df.loc[~in_mhc].copy()


def make_score(df):
    """Return 3-column score file: Genotype, A1, weight (weights >= 0)."""
    small = df[["Genotype", "REF", "ALT", "beta"]].copy()

    small_bad = small.loc[small["beta"] < 0].copy()
    small_good = small.loc[small["beta"] > 0].copy()

    small_bad["weight"] = -small_bad["beta"]
    small_bad["A1"] = small_bad["REF"]
    small_bad = small_bad[["Genotype", "A1", "weight"]]

    small_good["weight"] = small_good["beta"]
    small_good["A1"] = small_good["ALT"]
    small_good = small_good[["Genotype", "A1", "weight"]]

    score = pd.concat([small_good, small_bad], axis=0)
    score["Genotype"] = score["Genotype"].str.replace("-", ":", regex=False)
    return score


# --------------------- Data prep ---------------------

def make_temp_file(Dhindsa_excel_file, protein, variant_type, missing_variants) :
    missing = pd.read_csv(missing_variants, sep="\t", header=None, names=["Genotype"])
    raw = pd.read_excel(Dhindsa_excel_file)

    small = raw[raw.get("model", "genotypic").eq("genotypic")] if "model" in raw.columns else raw
    small = small[["Protein", "Genotype", "beta", "se", "Match(?Cis)", "cis_trans_position_1mb_coding_region"]].copy()
    small = small[small["Protein"] == protein].copy()

    if small.empty:
        print(f"Unable to find {protein}")
        return None

    # Normalize IDs FIRST so the overlap with 'missing' is on the same format
    small["Genotype"] = (
        small["Genotype"].astype(str)
        .str.replace("-", ":", regex=False)
        .str.replace("X:", "23:", regex=False)
    )

    # --- count overlap with --missing BEFORE we filter them out
    input_before_missing = small["Genotype"].nunique()
    excluded_by_missing = small["Genotype"].isin(missing["Genotype"]).sum()
    print(f"[MISSING] {protein}: excluded {excluded_by_missing} SNPs (overlap with --missing)")

    # Now actually drop them
    small = small[~small["Genotype"].isin(missing["Genotype"])].copy()
    if small.empty:
        print(f"No genetic variants at all for protein : {protein} after --missing exclusion")
        return None

    parts = small["Genotype"].str.split(":", n=3, expand=True)
    small["CHR"] = pd.to_numeric(parts[0], errors="coerce").astype("Int64")
    small["BP"]  = pd.to_numeric(parts[1], errors="coerce")
    small["REF"] = parts[2]
    small["ALT"] = parts[3]
    small = small.dropna(subset=["CHR", "BP"]).copy()
    small["CHR"] = small["CHR"].astype(int)
    small["Z"] = small["beta"] / small["se"]
    small["abs_Z"] = small["Z"].abs()

    # Select by variant_type
    if variant_type == "cis":
        out = small[small["Match(?Cis)"] == "Yes"].copy()
        if out.empty: print(f"No cis variants for protein : {protein}")
    elif variant_type == "trans":
        out = small[small["Match(?Cis)"] == "No"].copy()
        if out.empty: print(f"No trans variants for protein : {protein}")
    elif variant_type == "combined":
        out = small
    elif variant_type == "trans_cis_pos":
        out = small[small["cis_trans_position_1mb_coding_region"] == "cis-position, trans-gene"].copy()
        if out.empty:
            print(f"No cis-position, trans-gene variants for protein : {protein}")
            return None
    else:
        print("You have not selected cis, trans, combined or trans_cis_pos - check your input for variant type")
        return None

    # carry the counts forward (for process_protein to use)
    out.attrs["excluded_by_missing"] = int(excluded_by_missing)
    out.attrs["input_before_missing"] = int(input_before_missing)
    return out


# --------------------- Fast LD pruning ---------------------

def prune_ld_fast(chr_df, ld_df, r2_threshold):
    """
    Given a single-chromosome dataframe (Genotype, BP, abs_Z) and its LD table (ID_A, ID_B, PHASED_R2),
    build components for edges with DISTANCE <= 10Mb and PHASED_R2 >= r2_threshold.
    Keep the SNP with the largest |Z| per component.
    """
    ids = chr_df["Genotype"].unique()
    if len(ids) == 0:
        empty = pd.DataFrame(columns=["SNP1", "SNP2", "DISTANCE", "R2", "NOTE", "FLAG"])
        return empty, []

    id_set = set(ids)
    ld = ld_df.loc[
        ld_df["ID_A"].isin(id_set) & ld_df["ID_B"].isin(id_set),
        ["ID_A", "ID_B", "PHASED_R2"]
    ].copy()

    pos = chr_df.set_index("Genotype")["BP"]
    absZ = chr_df.set_index("Genotype")["abs_Z"]

    # Distance and r² filters
    if not ld.empty:
        ld["DISTANCE"] = (pos.reindex(ld["ID_A"]).to_numpy() - pos.reindex(ld["ID_B"]).to_numpy())
        ld["DISTANCE"] = np.abs(ld["DISTANCE"])
        ld = ld.loc[ld["DISTANCE"] <= 10_000_000]
        ld = ld.loc[ld["PHASED_R2"] >= r2_threshold]

    if ld.empty:
        # Everything independent
        empty = pd.DataFrame(columns=["SNP1", "SNP2", "DISTANCE", "R2", "NOTE", "FLAG"])
        return empty, ids.tolist()

    # Build components
    uf = UnionFind()
    for a, b in ld[["ID_A", "ID_B"]].itertuples(index=False):
        uf.union(a, b)

    comp = pd.Series({g: uf.find(g) for g in ids}, name="COMP")
    nodes = chr_df[["Genotype", "abs_Z"]].copy()
    nodes["COMP"] = nodes["Genotype"].map(comp)

    # Winner = max |Z|
    idx = nodes.groupby("COMP")["abs_Z"].idxmax()
    winners_df = nodes.loc[idx, ["COMP", "Genotype", "abs_Z"]].rename(
        columns={"Genotype": "WINNER", "abs_Z": "WINNER_absZ"}
    )
    nodes = nodes.merge(winners_df, on="COMP", how="left")

    keep = winners_df["WINNER"].tolist()
    losers = nodes.loc[~nodes["Genotype"].isin(keep)].copy()

    # Build concise log: loser -> winner
    ld_idx = ld_df.set_index(["ID_A", "ID_B"])["PHASED_R2"]

    def get_r2(a, b):
        try:
            return float(ld_idx.loc[(a, b)])
        except KeyError:
            try:
                return float(ld_idx.loc[(b, a)])
            except KeyError:
                return np.nan

    logs = []
    for g, w in zip(losers["Genotype"], losers["WINNER"]):
        dist = abs(int(pos[g]) - int(pos[w]))
        r2_val = get_r2(g, w)
        note = f"{w} wins with abs z {float(absZ[w])} compared to {g} with {float(absZ[g])}"
        logs.append({
            "SNP1": w,
            "SNP2": g,
            "DISTANCE": dist,
            "R2": r2_val,
            "NOTE": note,
            "FLAG": "OK"
        })

    log_df = pd.DataFrame(logs, columns=["SNP1", "SNP2", "DISTANCE", "R2", "NOTE", "FLAG"])
    return log_df, keep


def process_chromosome_fast(chromosome, ld_file, df, r2):
    chr_df = df.loc[df["CHR"] == chromosome].copy()
    if chr_df.empty:
        empty = pd.DataFrame(columns=["SNP1", "SNP2", "DISTANCE", "R2", "NOTE", "FLAG"])
        return empty, [], []

    chr_key = f"chr{chromosome}"
    try:
        ld_matrix_chr = pd.read_hdf(ld_file, key=chr_key)
    except (ValueError, KeyError, FileNotFoundError) as e:
        raise RuntimeError(f"Failed to read LD matrix for {chr_key} from {ld_file}: {e}") from e

    needed = {"ID_A", "ID_B", "PHASED_R2"}
    if not needed.issubset(ld_matrix_chr.columns):
        raise ValueError(f"LD matrix for {chr_key} is missing columns: {needed - set(ld_matrix_chr.columns)}")

    log_df, keep = prune_ld_fast(chr_df, ld_matrix_chr, r2)
    to_remove = [g for g in chr_df["Genotype"].unique() if g not in set(keep)]
    return log_df, to_remove, keep


# --------------------- Orchestration ---------------------

def process_protein(protein_name, r2, ld_file, variant_type, mhc_flag, out_dir, PROTEIN_DF):
    if PROTEIN_DF is None or PROTEIN_DF.empty:
        print(f"{protein_name} contains no {variant_type} hits")
        return

    df = PROTEIN_DF.copy()
    df["CHR"] = df["CHR"].astype(int)
    df = df.sort_values(["CHR", "BP"])

    # MHC audit (purely informational; does not drop anything)
    mhc_mask = (df["CHR"] == 6) & df["BP"].between(28_510_120, 33_480_577)
    total_mhc = df.loc[mhc_mask, "Genotype"].nunique()
    print(f"[MHC] input: {total_mhc} SNPs in chr6:28,510,120–33,480,577")

    # —— per-chromosome pruning (use your fast or original function) ——
    chromosomes = sorted(df["CHR"].unique().tolist())
    logs = []
    results = {}

    for chrom in chromosomes:
        log_df, dropped, keep = process_chromosome_fast(chrom, ld_file, df, r2)  # or process_chromosome
        if not log_df.empty:
            logs.append(log_df)
        results[chrom] = keep

        if chrom == 6:
            kept_chr6 = len(set(keep) & set(df.loc[mhc_mask, "Genotype"]))
            print(f"[MHC][chr6] kept after pruning: {kept_chr6}/{total_mhc} at r2={r2}")

    main_log = (pd.concat(logs, ignore_index=True)
                if logs else pd.DataFrame(columns=["SNP1","SNP2","DISTANCE","R2","NOTE","FLAG"]))
    flattened_keep = [g for keep in results.values() for g in keep]
    keep_df = df.loc[df["Genotype"].isin(flattened_keep)].copy()

    # —— summaries ——
    excluded_by_missing = int(PROTEIN_DF.attrs.get("excluded_by_missing", 0))
    input_before_missing = int(PROTEIN_DF.attrs.get("input_before_missing", 0))
    n_input_after_missing = df["Genotype"].nunique()
    n_kept = keep_df["Genotype"].nunique()
    removed_by_pruning = n_input_after_missing - n_kept

    kept_mhc = keep_df.loc[mhc_mask, "Genotype"].nunique()
    mhc_pct = (kept_mhc / total_mhc * 100) if total_mhc else 0.0

    print(f"[MISSING] overlap excluded: {excluded_by_missing} (from {input_before_missing} pre-missing candidates)")
    print(f"[PRUNING] removed by LD pruning: {removed_by_pruning}  (kept {n_kept}/{n_input_after_missing})")
    print(f"[MHC] kept after pruning: {kept_mhc}/{total_mhc} ({mhc_pct:.1f}%)")

    # Build outputs
    score_file = make_score(keep_df)
    folder_name = os.path.join(out_dir, protein_name)
    os.makedirs(folder_name, exist_ok=True)
    out_panel = os.path.join(folder_name, f"{protein_name}_{variant_type}_{r2}.panel")
    out_log   = os.path.join(folder_name, f"{protein_name}_{variant_type}_{r2}.log")
    out_sum   = os.path.join(folder_name, f"{protein_name}_{variant_type}_{r2}.summary.txt")

    main_log.to_csv(out_log, sep="\t", index=False, header=False)
    if not score_file.empty:
        score_file.to_csv(out_panel, sep="\t", index=False, header=False)

    # Write a concise summary file
    with open(out_sum, "w") as fh:
        fh.write(f"protein\t{protein_name}\n")
        fh.write(f"variant_type\t{variant_type}\n")
        fh.write(f"r2\t{r2}\n")
        fh.write(f"input_before_missing\t{input_before_missing}\n")
        fh.write(f"excluded_by_missing\t{excluded_by_missing}\n")
        fh.write(f"input_after_missing\t{n_input_after_missing}\n")
        fh.write(f"kept_after_pruning\t{n_kept}\n")
        fh.write(f"removed_by_pruning\t{removed_by_pruning}\n")
        fh.write(f"mhc_input\t{total_mhc}\n")
        fh.write(f"mhc_kept\t{kept_mhc}\n")
        fh.write(f"mhc_kept_pct\t{mhc_pct:.2f}\n")

    print(f"[DONE] panel: {out_panel}")
    print(f"[DONE] log:   {out_log}")
    print(f"[DONE] summary: {out_sum}")


def main(protein_name, output, ld_matrix, r2, excel_file, mhc_filter, variant_type, missing_variants):
    # Basic validations
    if not (0.0 <= float(r2) <= 1.0):
        raise ValueError(f"--r2 must be between 0 and 1 (got {r2})")
    if not os.path.exists(ld_matrix):
        raise FileNotFoundError(f"LD matrix file not found: {ld_matrix}")
    if not os.path.exists(excel_file):
        raise FileNotFoundError(f"Excel file not found: {excel_file}")
    if not os.path.exists(missing_variants):
        raise FileNotFoundError(f"Missing-variants file not found: {missing_variants}")

    r2 = float(r2)

    var_df = make_temp_file(excel_file, protein_name, variant_type, missing_variants)
    if var_df is None or var_df.empty:
        print(f"{protein_name} contains no {variant_type} hits after filtering.")
        return

    process_protein(protein_name, r2, ld_matrix, variant_type, mhc_filter, output, var_df)


# --------------------- Entrypoint ---------------------

if __name__ == '__main__':
    try:
        protein_name, output, ld_matrix, r2, excel_file, mhc_filter, variant_type, missing_variants = check_arg(sys.argv[1:])
        main(protein_name, output, ld_matrix, r2, excel_file, mhc_filter, variant_type, missing_variants)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
