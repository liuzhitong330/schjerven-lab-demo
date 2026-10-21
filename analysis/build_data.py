#!/usr/bin/env python3
"""Rebuild a descriptive, model-stratified IK1-response assay shortlist.

No alignment, DESeq2, significance testing or full ssvQC execution is claimed.
Requires Python >=3.10 and numpy >=1.24. See README.md for source restrictions.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys
import urllib.request

import numpy as np

SOURCE_FILES = {
    "counts": {
        "file": "GSE310958_MXP5_PDX2_raw_counts.tsv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE310nnn/GSE310958/suppl/GSE310958_MXP5_PDX2_raw_counts.tsv.gz",
        "sha256": "67552dfdb59d2f2d2d65518b99e14052e712bc986ed3e2b8cb334077a9d61d5c",
    },
    "gencode": {
        "file": "wgEncodeGencodeAttrsV36.txt.gz",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeAttrsV36.txt.gz",
        "sha256": "2098b0e7cd1ae412b1edfe93d267405175eadb63f63b9e8315b6bddda0a7c6be",
    },
    "gencode_schema": {
        "file": "wgEncodeGencodeAttrsV36.sql",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeAttrsV36.sql",
        "sha256": "5c238ca6074985eabebe572b0f182382c965221439108c727c9668c3ca3dac22",
    },
    "metadata": {
        "file": "GSE310958_family.soft.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE310nnn/GSE310958/soft/GSE310958_family.soft.gz",
        "sha256": "e7d03a48c73671930c0b1c0f2a9847e5bb5367beb20d573a8dee83ce6593be88",
    },
    "legacy_metadata": {
        "file": "GSE90659_family.soft.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90659/soft/GSE90659_family.soft.gz",
        "sha256": "603968f4615d468d040d7f9f802570d070fc5c0a4eb4de81c6c512b0335d1793",
    },
    "ssvqc_source": {
        "file": "ssvQC_correlation.R",
        "url": "https://raw.githubusercontent.com/FrietzeLabUVM/ssvQC/90a38bf72a97db43cdcdbb4e68bea9a4ceb389c5/R/ssvQC_correlation.R",
        "sha256": "62c0fcf160938cd7e7009a0d386e1b6cc9088f7c254a9c27049bf4f1985328be",
    },
}

SAMPLES = [
    {"id": "GSM9313472", "column": "MXP5_GFP_rep1", "label": "MXP5 GFP 1", "model": "MXP5", "condition": "GFP", "replicate": 1, "origin": "GSE310958"},
    {"id": "GSM9313473", "column": "MXP5_GFP_rep2", "label": "MXP5 GFP 2", "model": "MXP5", "condition": "GFP", "replicate": 2, "origin": "GSE310958"},
    {"id": "GSM9313474", "column": "MXP5_IK1_rep1", "label": "MXP5 IK1 1", "model": "MXP5", "condition": "IK1", "replicate": 1, "origin": "GSE310958"},
    {"id": "GSM9313475", "column": "MXP5_IK1_rep2", "label": "MXP5 IK1 2", "model": "MXP5", "condition": "IK1", "replicate": 2, "origin": "GSE310958"},
    {"id": "GSM2410132", "column": "GSM2410132_PDX2_GFP_rep1", "label": "PDX2 GFP 1", "model": "PDX2", "condition": "GFP", "replicate": 1, "origin": "GSE90659, reprocessed in GSE310958"},
    {"id": "GSM2410133", "column": "GSM2410133_PDX2_GFP_rep2", "label": "PDX2 GFP 2", "model": "PDX2", "condition": "GFP", "replicate": 2, "origin": "GSE90659, reprocessed in GSE310958"},
    {"id": "GSM2410140", "column": "GSM2410140_PDX2_IK_rep1", "label": "PDX2 IK1 1", "model": "PDX2", "condition": "IK1", "replicate": 1, "origin": "GSE90659, reprocessed in GSE310958"},
    {"id": "GSM2410141", "column": "GSM2410141_PDX2_IK_rep2", "label": "PDX2 IK1 2", "model": "PDX2", "condition": "IK1", "replicate": 2, "origin": "GSE90659, reprocessed in GSE310958"},
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def acquire(folder: Path, allow_download: bool) -> dict[str, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, spec in SOURCE_FILES.items():
        path = folder / spec["file"]
        if not path.exists():
            if not allow_download:
                raise FileNotFoundError(f"Missing {path}; use --download or supply all inputs")
            request = urllib.request.Request(spec["url"], headers={"User-Agent": "Research-reproducibility/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as target:
                while chunk := response.read(1024 * 1024):
                    target.write(chunk)
        actual = sha256(path)
        if actual != spec["sha256"]:
            raise ValueError(f"Source hash mismatch for {path.name}: {actual}. Review changed source explicitly.")
        paths[key] = path
    return paths


def load_counts(path: Path):
    counts_by_id = {}
    multiplicity = {}
    original_rows = 0
    with gzip.open(path, "rt") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        assert header == ["GeneID"] + [s["column"] for s in SAMPLES], header
        for row in reader:
            assert len(row) == 9 and re.fullmatch(r"ENSG\d{11}", row[0])
            values = [int(x) for x in row[1:]]
            assert all(x >= 0 for x in values)
            original_rows += 1
            if row[0] not in counts_by_id:
                counts_by_id[row[0]] = values
            else:
                counts_by_id[row[0]] = [a + b for a, b in zip(counts_by_id[row[0]], values)]
            multiplicity[row[0]] = multiplicity.get(row[0], 0) + 1
    ids = list(counts_by_id)
    matrix = np.array([counts_by_id[x] for x in ids], dtype=np.int64)
    return ids, matrix, {"original_rows": original_rows, "unique_ids": len(ids), "duplicated_ids": [key for key, n in multiplicity.items() if n > 1]}


def load_annotation(path: Path):
    annotations = {}
    transcript_rows = 0
    with gzip.open(path, "rt") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 14, fields
            transcript_rows += 1
            stable_id = fields[0].split(".")[0]
            assert re.fullmatch(r"ENSG\d{11}", stable_id), stable_id
            annotations.setdefault(stable_id, set()).add((fields[1], fields[2]))
    unambiguous = {key: next(iter(value)) for key, value in annotations.items() if len(value) == 1}
    ambiguous = {key: sorted(value) for key, value in annotations.items() if len(value) > 1}
    return unambiguous, {"source": "Official UCSC hg38 wgEncodeGencodeAttrsV36 transcript attributes collapsed to unambiguous geneId/geneName/geneType tuples", "transcript_records": transcript_rows, "unique_ids": len(annotations), "ambiguous_ids": ambiguous}


def median_ratio_factors(matrix):
    """Classical positive-in-all-samples median ratio, independently per model.

    This is only a size-factor calculation, not DESeq2 dispersion/model fitting.
    Factor geometric mean is rescaled to one within each four-library model.
    """
    output = np.empty(matrix.shape[1], dtype=float)
    reference_counts = {}
    for model, offset in [("MXP5", 0), ("PDX2", 4)]:
        part = matrix[:, offset:offset + 4]
        keep = np.all(part > 0, axis=1)
        reference_counts[model] = int(np.sum(keep))
        assert reference_counts[model] > 0
        reference = np.exp(np.mean(np.log(part[keep]), axis=1))
        factors = np.median(part[keep] / reference[:, None], axis=0)
        factors = factors / np.exp(np.mean(np.log(factors)))
        assert np.all(factors > 0)
        output[offset:offset + 4] = factors
    return output, reference_counts


def pearson_from_ssvqc(matrix):
    """Source-derived adaptation of ssvQC's feature-by-sample cor(mat).

    Original: R/ssvQC_correlation.R, commit pinned in SOURCE_FILES. Original
    features are peak read counts. Here rows are genes, never called peak QC.
    Pearson's correlation is unchanged; neither full R nor ssvQC is executed.
    """
    centered = matrix.astype(float) - np.mean(matrix, axis=0)
    covariance_numerator = centered.T @ centered
    norms = np.sqrt(np.diag(covariance_numerator))
    if np.any(norms == 0):
        raise ValueError("Cannot correlate a constant sample column")
    return covariance_numerator / np.outer(norms, norms)


def scalar_pearson(x, y):
    """Independent scalar formulation used for checks (Python math.fsum)."""
    mx, my = math.fsum(x) / len(x), math.fsum(y) / len(y)
    xc, yc = [v - mx for v in x], [v - my for v in y]
    numerator = math.fsum(a * b for a, b in zip(xc, yc))
    return numerator / math.sqrt(math.fsum(v * v for v in xc) * math.fsum(v * v for v in yc))


def safe_log_ratio(numerator, denominator):
    result = np.full_like(numerator, np.nan, dtype=float)
    valid = (numerator > 0) & (denominator > 0)
    result[valid] = np.log2(numerator[valid] / denominator[valid])
    return result


def metrics(matrix, divisor, cpm):
    normalized = matrix / divisor
    mean_effects, pair_effects, expression = [], [], []
    for offset in (0, 4):
        ctrl = normalized[:, offset:offset + 2]
        ik1 = normalized[:, offset + 2:offset + 4]
        mean_effects.append(safe_log_ratio(np.mean(ik1, axis=1), np.mean(ctrl, axis=1)))
        pair_effects.append(np.stack([safe_log_ratio(ik1[:, j], ctrl[:, k]) for j, k in itertools.product(range(2), repeat=2)], axis=1))
        expression.append(np.minimum(np.mean(cpm[:, offset:offset + 2], axis=1), np.mean(cpm[:, offset + 2:offset + 4], axis=1)))
    return {"effect": np.stack(mean_effects, axis=1), "pairs": np.stack(pair_effects, axis=1), "expression": np.stack(expression, axis=1)}


def selected(met, effect=0.5, min_cpm=1.0, strict=True):
    means, pairs, expression = met["effect"], met["pairs"], met["expression"]
    eligible = np.all(expression >= min_cpm, axis=1)
    up = np.all(means >= effect, axis=1)
    down = np.all(means <= -effect, axis=1)
    if strict:
        up &= np.all(pairs >= effect, axis=(1, 2))
        down &= np.all(pairs <= -effect, axis=(1, 2))
    return eligible & (up | down)


def finite_list(values):
    return [round(float(x), 9) if np.isfinite(x) else None for x in values]


def run_checks(matrix, factors, raw_corr, log_corr, log_cpm, met):
    fixture = np.array([[1., 3., 10.], [2., 5., 8.], [3., 7., 6.], [4., 9., 4.]])
    expected = np.array([[1., 1., -1.], [1., 1., -1.], [-1., -1., 1.]])
    assert np.allclose(pearson_from_ssvqc(fixture), expected, atol=1e-14)
    raw_err = float(np.max(np.abs(raw_corr - np.corrcoef(matrix.T))))
    log_err = float(np.max(np.abs(log_corr - np.corrcoef(log_cpm.T))))
    scalar_err = 0.0
    for i in range(8):
        for j in range(i, 8):
            scalar_err = max(scalar_err, abs(raw_corr[i, j] - scalar_pearson(matrix[:, i].tolist(), matrix[:, j].tolist())))
    assert raw_err < 1e-12 and log_err < 1e-12 and scalar_err < 1e-12
    assert np.allclose(np.diag(raw_corr), 1)
    assert np.allclose(raw_corr, raw_corr.T)
    assert abs(float(np.prod(factors[:4])) - 1) < 1e-12
    assert abs(float(np.prod(factors[4:])) - 1) < 1e-12
    assert np.isnan(safe_log_ratio(np.array([0., 1.]), np.array([1., 0.]))).all()
    assert np.allclose(safe_log_ratio(np.array([4., 1.]), np.array([1., 4.])), [2., -2.])
    normalization_fixture = np.array([[2, 4, 8, 16, 2, 4, 8, 16], [6, 12, 24, 48, 6, 12, 24, 48]])
    fixture_factors, _ = median_ratio_factors(normalization_fixture)
    expected_factors = np.tile(np.array([1, 2, 4, 8]) / np.sqrt(8), 2)
    assert np.allclose(fixture_factors, expected_factors)
    independent_factors = []
    for offset in (0, 4):
        ratios = [[] for _ in range(4)]
        for row in matrix[:, offset:offset + 4].tolist():
            if min(row) > 0:
                geo = math.exp(math.fsum(math.log(v) for v in row) / 4)
                for j in range(4):
                    ratios[j].append(row[j] / geo)
        import statistics
        medians = [statistics.median(column) for column in ratios]
        scaling = math.exp(math.fsum(math.log(v) for v in medians) / 4)
        independent_factors.extend(v / scaling for v in medians)
    size_factor_err = float(np.max(np.abs(factors - independent_factors)))
    assert size_factor_err < 1e-12
    for norm in met.values():
        assert np.all(~selected(norm, strict=True) | selected(norm, strict=False))
        assert np.all(~selected(norm, effect=1) | selected(norm, effect=0.5))
    return {"status": "PASS", "tests": ["Pearson known linear and reversed fixture", "Pearson NumPy corrcoef equivalence", "Pearson independent scalar equivalence on every sample pair", "Correlation symmetry and diagonal", "Size-factor geometric means are one within model", "Known proportional-count size-factor fixture", "Independent scalar median-ratio calculation on full data", "Zero ratios produce missing values", "Exact ratio sanity check", "Strict shortlist subset of mean-only shortlist", "Increasing effect floor never increases shortlist"], "max_abs_error_numpy_raw": raw_err, "max_abs_error_numpy_log_cpm": log_err, "max_abs_error_scalar_raw": float(scalar_err), "max_abs_error_scalar_size_factors": size_factor_err, "R_or_full_ssvQC_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path("inputs"), help="Folder containing pinned public inputs")
    parser.add_argument("--download", action="store_true", help="Download missing public inputs; require exact SHA256")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1], help="Site folder receiving data.js and analysis outputs")
    args = parser.parse_args()
    paths = acquire(args.inputs, args.download)
    assert "cor_mat = cor(mat)" in paths["ssvqc_source"].read_text()
    ids, matrix, count_audit = load_counts(paths["counts"])
    annotation, annotation_audit = load_annotation(paths["gencode"])
    denominators = np.sum(matrix, axis=0)
    cpm = matrix / denominators * 1e6
    factors, reference_counts = median_ratio_factors(matrix)
    raw_corr = pearson_from_ssvqc(matrix)
    log_cpm = np.log2(cpm + 1)
    log_corr = pearson_from_ssvqc(log_cpm)
    met = {"median_ratio": metrics(matrix, factors, cpm), "cpm": metrics(matrix, denominators / 1e6, cpm)}
    checks = run_checks(matrix, factors, raw_corr, log_corr, log_cpm, met)
    genes = [{"id": gene, "symbol": annotation.get(gene, (gene, "unmapped"))[0], "type": annotation.get(gene, (gene, "unmapped"))[1], "counts": counts.tolist()} for gene, counts in zip(ids, matrix)]
    samples = [{**s, "induction": "1 microgram/mL doxycycline", "time_hours": 24, "paired": False, "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + s["id"]} for s in SAMPLES]
    provenance = {
        "schema_version": "1.0", "checked_date": "2026-10-21", "study_title": "IKAROS Gene Regulatory Network Reveal ERG as a Vulnerability in B-cell Acute Lymphoblastic Leukemia", "doi": "10.64898/2025.12.29.696402", "study_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12776271/", "study_status": "Preprint version 1; not peer reviewed at source inspection", "geo": "GSE310958", "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310958", "sources": SOURCE_FILES, "count_audit": count_audit, "annotation_audit": annotation_audit, "unmapped_ids": [gene for gene in ids if gene not in annotation], "genes_excluded": 0,
        "correlation_source": {"repository": "https://github.com/FrietzeLabUVM/ssvQC", "commit": "90a38bf72a97db43cdcdbb4e68bea9a4ceb389c5", "file": "R/ssvQC_correlation.R", "relationship": "Schjerven-coauthored collaborator software, verified through the 2021 primary publication", "adaptation": "Source-defined Pearson feature-by-sample count matrix primitive adapted from peak counts to gene counts; not full ssvQC or R execution. Added log2(CPM+1) sensitivity view.", "license": "MIT, Copyright 2021 Joseph R Boyd"},
        "limitations": ["Eight libraries: n=2 biological replicates per condition and model; four cross-condition ratios are sensitivity comparisons, not four replicates or confidence intervals.", "PDX2 data were reused from the 2017 GSE90659 deposit and reprocessed in GSE310958. Model and historical preparation/instrument are confounded; compare within-model effects, not absolute cross-model expression.", "All effects are descriptive; no p-values, differential-expression significance, causal target discovery or clinical inference.", "RNA counts do not establish protein activity or surface abundance; gene-level IKZF1 cannot separate endogenous mutant and induced IK1 transcripts.", "Selection criteria are exploratory assay-prioritization choices, not validated biological/QC cutoffs. Independent assay validation and appropriate controls are required.", "Median-ratio normalization assumes sufficiently many comparable genes; no spike-ins or absolute RNA measurement. Normalization sensitivity is shown rather than assumed away."],
        "formulas": {"mean_effect": "log2(mean(normalized IK1 counts) / mean(normalized GFP counts)); missing if either mean is zero; no pseudocount", "pair_range": "Range of four log2(normalized IK1_j / normalized GFP_k) combinations within each model; missing if either count is zero", "median_ratio": "Per model: include genes positive in all four samples; divide each count by its gene geometric mean; take each column median; rescale four factors to geometric mean one", "expression_floor": "Both GFP and IK1 mean CPM >= selected CPM floor in each model (all four condition means must pass)", "strict_support": "Same mean direction in both models and all four finite cross-condition ratios per model meet the selected effect floor in that direction", "log_cpm_QC": "log2(assigned-count CPM + 1), used only for contextual sample correlation and not effect ratios"}
    }
    data = {"samples": samples, "genes": genes, "normalization": {"cpm_denominators": denominators.tolist(), "median_ratio_factors": factors.tolist(), "median_ratio_reference_gene_counts": reference_counts}, "qc": {"raw_correlation": raw_corr.tolist(), "log_cpm_correlation": log_corr.tolist(), "nonzero_genes": int(np.sum(np.any(matrix > 0, axis=1))), "checks": checks}, "provenance": provenance}
    out = args.output
    (out / "analysis").mkdir(parents=True, exist_ok=True)
    (out / "data.js").write_text("window.SCHJERVEN_DATA=" + json.dumps(data, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + ";\n")
    (out / "analysis" / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    (out / "analysis" / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    with (out / "analysis" / "sample-manifest.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
        writer.writeheader(); writer.writerows(samples)
    summary = {"row_count": len(ids), "annotation_mapped": len(ids) - len(provenance["unmapped_ids"]), "assigned_count_totals": denominators.tolist(), "median_ratio_factors": factors.tolist(), "median_ratio_reference_gene_counts": reference_counts, "defaults": {"effect_floor_log2": 0.5, "expression_floor_CPM": 1.0, "strict_all_comparisons": True}, "normalizations": {}, "selected_gene_examples": {}}
    masks = {}
    for norm, measurements in met.items():
        strict = selected(measurements, strict=True); permissive = selected(measurements, strict=False)
        masks[norm] = strict
        scores = np.min(np.abs(measurements["pairs"]), axis=(1, 2))
        order = sorted(np.flatnonzero(strict), key=lambda i: (-scores[i], genes[i]["id"]))
        summary["normalizations"][norm] = {"strict_count": int(np.sum(strict)), "mean_only_count": int(np.sum(permissive)), "removed_by_pair_rule": int(np.sum(permissive & ~strict)), "top20": [{"id": ids[i], "symbol": genes[i]["symbol"], "type": genes[i]["type"], "mean_log2_effects": finite_list(measurements["effect"][i]), "pair_ranges": [finite_list([np.min(row), np.max(row)]) for row in measurements["pairs"][i]], "min_abs_pair_effect": float(scores[i])} for i in order[:20]]}
    summary["strict_both_normalizations"] = int(np.sum(masks["median_ratio"] & masks["cpm"]))
    summary["strict_median_ratio_only"] = int(np.sum(masks["median_ratio"] & ~masks["cpm"]))
    summary["strict_cpm_only"] = int(np.sum(masks["cpm"] & ~masks["median_ratio"]))
    wanted = {"ERG", "IKZF1", "CCND2", "CCDC26", "MYC", "PAX5", "BCL6", "KLF2", "CD19", "CD79A", "CD79B", "MS4A1", "VPREB1", "IGLL1", "DNTT", "RAG1", "RAG2", "DUSP6", "GAPDH", "ACTB", "RPLP0", "HPRT1", "TBP"}
    for i, gene in enumerate(genes):
        if gene["symbol"] in wanted:
            summary["selected_gene_examples"][gene["symbol"]] = {"id": gene["id"], "counts": gene["counts"], "CPM": finite_list(cpm[i]), **{norm: {"mean_log2_effects": finite_list(measurements["effect"][i]), "pair_ranges": [finite_list([np.min(row), np.max(row)]) for row in measurements["pairs"][i]], "strict_selected": bool(masks[norm][i])} for norm, measurements in met.items()}}
    (out / "analysis" / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"gene_count": len(genes), "mapped": summary["annotation_mapped"], "strict_median_ratio": summary["normalizations"]["median_ratio"]["strict_count"], "strict_cpm": summary["normalizations"]["cpm"]["strict_count"], "checks": checks["status"], "data_bytes": (out / "data.js").stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
