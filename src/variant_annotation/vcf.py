from __future__ import annotations

import csv
import shutil
import subprocess
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .clinvar import ClinVarStore


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "outputs"
PATHOGENIC_LABELS = {"pathogenic", "likely pathogenic"}


@dataclass(frozen=True)
class VcfVariant:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: str


def resolve_project_path(path: str | Path) -> Path:
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    return candidate.resolve()


def parse_vcf(vcf_path: str | Path) -> list[VcfVariant]:
    variants: list[VcfVariant] = []
    with resolve_project_path(vcf_path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Expected at least 8 VCF columns, got {len(fields)}: {line.rstrip()}")
            chrom, pos, variant_id, ref, alt, qual, filter_value, info = fields[:8]
            variants.append(
                VcfVariant(
                    chrom=chrom,
                    pos=int(pos),
                    variant_id=variant_id,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    filter=filter_value,
                    info=info,
                )
            )
    return variants


def _python_normalize_vcf(input_path: Path, output_path: Path) -> int:
    """Split multiallelic rows into one row per ALT allele."""
    variant_count = 0
    with input_path.open(encoding="utf-8") as source, output_path.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith("#"):
                target.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Expected at least 8 VCF columns, got {len(fields)}: {line.rstrip()}")
            for alt in fields[4].split(","):
                normalized = fields.copy()
                normalized[4] = alt
                target.write("\t".join(normalized) + "\n")
                variant_count += 1
    return variant_count


def normalize_vcf(input_vcf: str, output_vcf: str | None = None) -> dict[str, Any]:
    input_path = resolve_project_path(input_vcf)
    if output_vcf is None:
        OUTPUT_DIR.mkdir(exist_ok=True)
        output_path = OUTPUT_DIR / f"{input_path.stem}.normalized.vcf"
    else:
        output_path = resolve_project_path(output_vcf)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    bcftools = shutil.which("bcftools")
    method = "custom python function"
    warning = None
    if bcftools:
        command = [bcftools, "norm", "-m", "-any", str(input_path), "-o", str(output_path)]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode == 0:
            method = "bcftools norm -m -any"
            variant_count = len(parse_vcf(output_path))
        else:
            warning = (
                "bcftools was found but normalization failed; proceeding using a "
                "custom Python function that only splits multiallelic ALT rows."
            )
            warnings.warn(warning, RuntimeWarning, stacklevel=2)
            variant_count = _python_normalize_vcf(input_path, output_path)
    else:
        warning = (
            "bcftools was not found; normalization will proceed using a custom "
            "Python function that only splits multiallelic ALT rows."
        )
        warnings.warn(warning, RuntimeWarning, stacklevel=2)
        variant_count = _python_normalize_vcf(input_path, output_path)

    result = {
        "input_vcf": str(input_path),
        "output_vcf": str(output_path),
        "method": method,
        "variant_count": variant_count,
    }
    if warning:
        result["warning"] = warning
    return result


def annotate_variants(variants: list[VcfVariant], store: ClinVarStore) -> list[dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    for variant in variants:
        result = store.query_variant(
            chrom=variant.chrom,
            pos=variant.pos,
            ref=variant.ref,
            alt=variant.alt,
        )
        annotations.append(
            {
                "chrom": variant.chrom,
                "pos": variant.pos,
                "id": variant.variant_id,
                "ref": variant.ref,
                "alt": variant.alt,
                "found": result["found"],
                "variation_id": result.get("variation_id", ""),
                "clinical_significance": result.get("clinical_significance", ""),
                "gene": result.get("gene", ""),
                "condition": result.get("condition", ""),
            }
        )
    return annotations


def write_annotation_tsv(rows: list[dict[str, Any]], output_tsv: Path) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "found",
        "variation_id",
        "clinical_significance",
        "gene",
        "condition",
    ]
    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_pathogenic_report(rows: list[dict[str, Any]], report_path: Path) -> list[dict[str, Any]]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    candidates = [
        row
        for row in rows
        if str(row["clinical_significance"]).lower() in PATHOGENIC_LABELS
    ]
    lines = [
        "# Candidate Pathogenic Variants",
        "",
        f"Total annotated variants: {len(rows)}",
        f"Candidate pathogenic/likely pathogenic variants: {len(candidates)}",
        "",
    ]
    if candidates:
        lines.append("| Variant | Clinical significance | Gene | Condition | Variation ID |")
        lines.append("| --- | --- | --- | --- | --- |")
        for row in candidates:
            variant = f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}"
            lines.append(
                f"| {variant} | {row['clinical_significance']} | {row['gene']} | "
                f"{row['condition']} | {row['variation_id']} |"
            )
    else:
        lines.append("No pathogenic or likely pathogenic variants were found.")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return candidates


def summarize_annotations(
    input_vcf: str,
    output_tsv: str | None = None,
    report_path: str | None = None,
    normalize: bool = True,
    store: ClinVarStore | None = None,
) -> dict[str, Any]:
    store = store or ClinVarStore()
    input_path = resolve_project_path(input_vcf)
    OUTPUT_DIR.mkdir(exist_ok=True)

    if normalize:
        normalized = normalize_vcf(str(input_path))
        vcf_for_annotation = Path(normalized["output_vcf"])
        normalization = normalized
    else:
        vcf_for_annotation = input_path
        normalization = None

    tsv_path = resolve_project_path(output_tsv) if output_tsv else OUTPUT_DIR / f"{input_path.stem}.annotations.tsv"
    report_output = (
        resolve_project_path(report_path)
        if report_path
        else OUTPUT_DIR / f"{input_path.stem}.pathogenic_report.md"
    )

    rows = annotate_variants(parse_vcf(vcf_for_annotation), store)
    write_annotation_tsv(rows, tsv_path)
    candidates = write_pathogenic_report(rows, report_output)

    return {
        "input_vcf": str(input_path),
        "normalized_vcf": str(vcf_for_annotation) if normalize else None,
        "normalization": normalization,
        "annotation_tsv": str(tsv_path),
        "pathogenic_report": str(report_output),
        "variant_count": len(rows),
        "found_count": sum(1 for row in rows if row["found"]),
        "missing_count": sum(1 for row in rows if not row["found"]),
        "candidate_count": len(candidates),
        "candidates": candidates,
    }
