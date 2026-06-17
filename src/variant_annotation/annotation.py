from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .clinvar import ClinVarStore
from .vcf import OUTPUT_DIR, PATHOGENIC_LABELS, VcfVariant, normalize_vcf, parse_vcf, resolve_project_path


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
DEFAULT_GNOMAD_PATH = DATA_DIR / "synthetic_gnomad.parquet"
DEFAULT_OMIM_PATH = DATA_DIR / "synthetic_omim.parquet"
WAREHOUSE_DIR = ROOT / "warehouse"


@dataclass(frozen=True)
class ParquetVariantStore:
    parquet_path: Path
    not_found_message: str

    def query_variant(self, chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
        frame = pd.read_parquet(self.parquet_path)
        chrom_value = str(chrom).removeprefix("chr")
        matches = frame[
            (frame["chrom"].astype(str).str.removeprefix("chr") == chrom_value)
            & (frame["pos"].astype(int) == int(pos))
            & (frame["ref"].astype(str).str.upper() == ref.upper())
            & (frame["alt"].astype(str).str.upper() == alt.upper())
        ]
        if matches.empty:
            return {
                "found": False,
                "query": {"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt},
                "message": self.not_found_message,
            }

        record = matches.iloc[0].to_dict()
        record["pos"] = int(record["pos"])
        record["found"] = True
        return record


@dataclass(frozen=True)
class OmimStore:
    parquet_path: Path = DEFAULT_OMIM_PATH

    def query_gene(self, gene: str) -> dict[str, Any]:
        frame = pd.read_parquet(self.parquet_path)
        matches = frame[frame["gene"].astype(str).str.upper() == gene.upper()]
        if matches.empty:
            return {
                "found": False,
                "query": {"gene": gene},
                "message": "Gene is not present in the synthetic OMIM parquet.",
            }

        record = matches.iloc[0].to_dict()
        record["found"] = True
        return record


@dataclass(frozen=True)
class AnnotationStores:
    clinvar: ClinVarStore
    gnomad: ParquetVariantStore
    omim: OmimStore

    @classmethod
    def defaults(cls) -> "AnnotationStores":
        return cls(
            clinvar=ClinVarStore(),
            gnomad=ParquetVariantStore(
                parquet_path=DEFAULT_GNOMAD_PATH,
                not_found_message="Variant is not present in the synthetic gnomAD parquet.",
            ),
            omim=OmimStore(),
        )


def _fallback_consequence(variant: VcfVariant) -> str:
    if len(variant.ref) != len(variant.alt):
        return "frameshift_variant" if abs(len(variant.ref) - len(variant.alt)) % 3 else "inframe_indel"
    if len(variant.ref) == 1 and len(variant.alt) == 1:
        return "missense_variant"
    return "sequence_variant"


def run_vep(
    input_vcf: str,
    output_tsv: str | None = None,
    cache_dir: str | None = None,
    offline: bool = True,
) -> dict[str, Any]:
    """Run VEP when available; otherwise create deterministic VEP-like annotations."""
    input_path = resolve_project_path(input_vcf)
    OUTPUT_DIR.mkdir(exist_ok=True)
    output_path = resolve_project_path(output_tsv) if output_tsv else OUTPUT_DIR / f"{input_path.stem}.vep.tsv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    vep = shutil.which("vep")
    if vep:
        command = [
            vep,
            "--input_file",
            str(input_path),
            "--output_file",
            str(output_path),
            "--tab",
            "--force_overwrite",
            "--no_stats",
            "--symbol",
            "--canonical",
        ]
        if offline:
            command.append("--offline")
        if cache_dir:
            command.extend(["--dir_cache", str(resolve_project_path(cache_dir))])
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        return {
            "input_vcf": str(input_path),
            "output_tsv": str(output_path),
            "method": "vep",
            "command": command,
            "returncode": completed.returncode,
            "stderr": completed.stderr,
            "stdout": completed.stdout,
        }

    rows = []
    for variant in parse_vcf(input_path):
        rows.append(
            {
                "chrom": variant.chrom,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "most_severe_consequence": _fallback_consequence(variant),
                "impact": "HIGH" if len(variant.ref) != len(variant.alt) else "MODERATE",
                "transcript_id": f"ENST_SYN_{str(variant.pos)[-6:]}",
                "vep_source": "fallback",
            }
        )
    _write_tsv(rows, output_path)
    return {
        "input_vcf": str(input_path),
        "output_tsv": str(output_path),
        "method": "deterministic fallback",
        "warning": "VEP executable was not found; wrote VEP-like annotations for local testing.",
        "variant_count": len(rows),
    }


def _write_tsv(rows: list[dict[str, Any]], output_path: Path) -> None:
    if not rows:
        output_path.write_text("", encoding="utf-8")
        return
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _annotation_rows(variants: list[VcfVariant], stores: AnnotationStores) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for variant in variants:
        clinvar = stores.clinvar.query_variant(variant.chrom, variant.pos, variant.ref, variant.alt)
        gnomad = stores.gnomad.query_variant(variant.chrom, variant.pos, variant.ref, variant.alt)
        gene = str(clinvar.get("gene", ""))
        omim = stores.omim.query_gene(gene) if gene else {"found": False}
        rows.append(
            {
                "chrom": variant.chrom,
                "pos": variant.pos,
                "id": variant.variant_id,
                "ref": variant.ref,
                "alt": variant.alt,
                "qual": variant.qual,
                "filter": variant.filter,
                "vep_consequence": _fallback_consequence(variant),
                "clinvar_found": clinvar["found"],
                "clinvar_variation_id": clinvar.get("variation_id", ""),
                "clinical_significance": clinvar.get("clinical_significance", ""),
                "gene": gene,
                "condition": clinvar.get("condition", ""),
                "gnomad_found": gnomad["found"],
                "gnomad_af": gnomad.get("allele_frequency", ""),
                "gnomad_ac": gnomad.get("allele_count", ""),
                "gnomad_an": gnomad.get("allele_number", ""),
                "gnomad_popmax": gnomad.get("popmax", ""),
                "omim_found": omim["found"],
                "omim_mim_number": omim.get("mim_number", ""),
                "omim_phenotype": omim.get("phenotype", ""),
                "omim_inheritance": omim.get("inheritance", ""),
            }
        )
    return rows


def _write_iceberg_like_table(rows: list[dict[str, Any]], table_path: Path) -> dict[str, Any]:
    table_path.mkdir(parents=True, exist_ok=True)
    data_path = table_path / "data.parquet"
    metadata_path = table_path / "metadata.json"
    frame = pd.DataFrame(rows)
    frame.to_parquet(data_path, index=False)
    metadata = {
        "format": "iceberg-like-local-parquet",
        "note": "Local stand-in for an Iceberg table. Use submit_spark_job for a real catalog write.",
        "row_count": len(frame),
        "columns": list(frame.columns),
        "data_files": [str(data_path)],
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return {"table_path": str(table_path), "data_path": str(data_path), "metadata_path": str(metadata_path)}


def _write_pathogenic_report(rows: list[dict[str, Any]], report_path: Path) -> list[dict[str, Any]]:
    candidates = [
        row
        for row in rows
        if str(row["clinical_significance"]).lower() in PATHOGENIC_LABELS
    ]
    lines = [
        "# Pathogenic/Likely Pathogenic Variant Report",
        "",
        f"Total variants: {len(rows)}",
        f"Pathogenic or likely pathogenic variants: {len(candidates)}",
        "",
    ]
    if candidates:
        lines.append("| Variant | Gene | Significance | gnomAD AF | OMIM phenotype |")
        lines.append("| --- | --- | --- | --- | --- |")
        for row in candidates:
            variant = f"{row['chrom']}:{row['pos']} {row['ref']}>{row['alt']}"
            lines.append(
                f"| {variant} | {row['gene']} | {row['clinical_significance']} | "
                f"{row['gnomad_af']} | {row['omim_phenotype']} |"
            )
    else:
        lines.append("No pathogenic or likely pathogenic variants were found.")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return candidates


def _write_qc_summary(rows: list[dict[str, Any]], qc_path: Path, normalization: dict[str, Any] | None) -> dict[str, Any]:
    variant_count = len(rows)
    qc = {
        "variant_count": variant_count,
        "pass_filter_count": sum(1 for row in rows if row["filter"] in {"PASS", "."}),
        "clinvar_hit_count": sum(1 for row in rows if row["clinvar_found"]),
        "gnomad_hit_count": sum(1 for row in rows if row["gnomad_found"]),
        "omim_hit_count": sum(1 for row in rows if row["omim_found"]),
        "pathogenic_or_likely_pathogenic_count": sum(
            1 for row in rows if str(row["clinical_significance"]).lower() in PATHOGENIC_LABELS
        ),
        "normalization": normalization,
    }
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    return qc


def annotate_vcf_realistic(
    input_vcf: str,
    output_table: str | None = None,
    pathogenic_report: str | None = None,
    qc_summary: str | None = None,
    normalize: bool = True,
    run_vep_step: bool = True,
    stores: AnnotationStores | None = None,
) -> dict[str, Any]:
    stores = stores or AnnotationStores.defaults()
    input_path = resolve_project_path(input_vcf)
    OUTPUT_DIR.mkdir(exist_ok=True)
    WAREHOUSE_DIR.mkdir(exist_ok=True)

    if normalize:
        normalization = normalize_vcf(str(input_path))
        vcf_for_annotation = Path(normalization["output_vcf"])
    else:
        normalization = None
        vcf_for_annotation = input_path

    vep_result = run_vep(str(vcf_for_annotation)) if run_vep_step else None
    rows = _annotation_rows(parse_vcf(vcf_for_annotation), stores)

    default_name = input_path.stem.replace(".", "_")
    table_path = resolve_project_path(output_table) if output_table else WAREHOUSE_DIR / f"{default_name}_annotations"
    report_path = (
        resolve_project_path(pathogenic_report)
        if pathogenic_report
        else OUTPUT_DIR / f"{input_path.stem}.pathogenic_likely_pathogenic.md"
    )
    qc_path = resolve_project_path(qc_summary) if qc_summary else OUTPUT_DIR / f"{input_path.stem}.qc_summary.json"

    table = _write_iceberg_like_table(rows, table_path)
    candidates = _write_pathogenic_report(rows, report_path)
    qc = _write_qc_summary(rows, qc_path, normalization)

    return {
        "input_vcf": str(input_path),
        "normalized_vcf": str(vcf_for_annotation) if normalize else None,
        "vep": vep_result,
        "annotated_table": table,
        "pathogenic_report": str(report_path),
        "qc_summary": str(qc_path),
        "variant_count": len(rows),
        "candidate_count": len(candidates),
        "qc": qc,
    }


def submit_spark_job(
    app_path: str,
    args: list[str] | None = None,
    master: str = "local[*]",
    dry_run: bool = True,
) -> dict[str, Any]:
    spark_submit = shutil.which("spark-submit")
    command = ["spark-submit", "--master", master, str(resolve_project_path(app_path)), *(args or [])]
    if dry_run or not spark_submit:
        return {
            "submitted": False,
            "dry_run": dry_run,
            "command": command,
            "reason": "spark-submit executable was not found." if not spark_submit else "dry_run was true.",
        }
    command[0] = spark_submit
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    return {
        "submitted": completed.returncode == 0,
        "dry_run": False,
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def run_nextflow_workflow(
    workflow: str,
    params: dict[str, str] | None = None,
    profile: str | None = None,
    work_dir: str | None = None,
    dry_run: bool = True,
) -> dict[str, Any]:
    nextflow = shutil.which("nextflow")
    command = ["nextflow", "run", workflow]
    for key, value in (params or {}).items():
        command.extend([f"--{key}", value])
    if profile:
        command.extend(["-profile", profile])
    if work_dir:
        command.extend(["-work-dir", str(resolve_project_path(work_dir))])
    if dry_run or not nextflow:
        return {
            "submitted": False,
            "dry_run": dry_run,
            "command": command,
            "reason": "nextflow executable was not found." if not nextflow else "dry_run was true.",
        }
    command[0] = nextflow
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    return {
        "submitted": completed.returncode == 0,
        "dry_run": False,
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }
