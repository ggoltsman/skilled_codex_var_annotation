from __future__ import annotations

from pathlib import Path
from typing import Any

from mcp.server.fastmcp import FastMCP

from .annotation import (
    DEFAULT_GNOMAD_PATH,
    GnomadStore,
    OmimStore,
    annotate_vcf_realistic,
    run_nextflow_workflow as run_nextflow_workflow_job,
    run_vep as run_vep_annotation_step,
    submit_spark_job as submit_spark_job_step,
)
from .clinvar import DEFAULT_CLINVAR_PATH, ClinVarStore
from .vcf import normalize_vcf as normalize_vcf_file
from .vcf import summarize_annotations as summarize_vcf_annotations


mcp = FastMCP("variant-annotation")
store = ClinVarStore(Path(DEFAULT_CLINVAR_PATH))
gnomad_store = GnomadStore(
    parquet_path=Path(DEFAULT_GNOMAD_PATH),
    not_found_message="Variant is not present in the synthetic gnomAD parquet.",
)
omim_store = OmimStore()


@mcp.tool()
def query_variant(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
    """Return a synthetic ClinVar-style annotation for one variant."""
    return store.query_variant(chrom=chrom, pos=pos, ref=ref, alt=alt)


@mcp.tool()
def query_clinvar(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
    """Return a synthetic ClinVar-style annotation for one variant."""
    return store.query_variant(chrom=chrom, pos=pos, ref=ref, alt=alt)


@mcp.tool()
def query_gnomad(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any]:
    """Return a synthetic gnomAD population-frequency annotation for one variant."""
    return gnomad_store.query_variant(chrom=chrom, pos=pos, ref=ref, alt=alt)


@mcp.tool()
def query_omim(gene: str) -> dict[str, Any]:
    """Return a synthetic OMIM disease annotation for one gene symbol."""
    return omim_store.query_gene(gene=gene)


@mcp.tool()
def normalize_vcf(input_vcf: str, output_vcf: str | None = None) -> dict[str, Any]:
    """Normalize a VCF by splitting multiallelic rows into one row per ALT allele."""
    return normalize_vcf_file(input_vcf=input_vcf, output_vcf=output_vcf)


@mcp.tool()
def summarize_annotations(
    input_vcf: str,
    output_tsv: str | None = None,
    report_path: str | None = None,
    normalize: bool = True,
) -> dict[str, Any]:
    """Annotate a VCF and write an annotation TSV plus pathogenic candidate report."""
    return summarize_vcf_annotations(
        input_vcf=input_vcf,
        output_tsv=output_tsv,
        report_path=report_path,
        normalize=normalize,
        store=store,
    )


@mcp.tool()
def run_vep_annotation(
    input_vcf: str,
    output_tsv: str | None = None,
    cache_dir: str | None = None,
    offline: bool = True,
) -> dict[str, Any]:
    """Run VEP when installed; otherwise emit deterministic VEP-like annotations."""
    return run_vep_annotation_step(
        input_vcf=input_vcf,
        output_tsv=output_tsv,
        cache_dir=cache_dir,
        offline=offline,
    )


@mcp.tool()
def annotate_vcf(
    input_vcf: str,
    output_table: str | None = None,
    pathogenic_report: str | None = None,
    qc_summary: str | None = None,
    llm_summary_report: str | None = None,
    normalize: bool = True,
    run_vep_step: bool = True,
    generate_llm_summary: bool = True,
    llm_model: str | None = None,
) -> dict[str, Any]:
    """Run the realistic VCF annotation workflow and write table/report/QC/LLM outputs."""
    return annotate_vcf_realistic(
        input_vcf=input_vcf,
        output_table=output_table,
        pathogenic_report=pathogenic_report,
        qc_summary=qc_summary,
        llm_summary_report=llm_summary_report,
        normalize=normalize,
        run_vep_step=run_vep_step,
        generate_llm_summary=generate_llm_summary,
        llm_model=llm_model,
    )


@mcp.tool()
def submit_spark_job(
    app_path: str,
    args: list[str] | None = None,
    master: str = "local[*]",
    dry_run: bool = True,
) -> dict[str, Any]:
    """Submit a Spark job, or return the spark-submit command in dry-run mode."""
    return submit_spark_job_step(app_path=app_path, args=args, master=master, dry_run=dry_run)


@mcp.tool()
def run_nextflow_workflow(
    workflow: str,
    params: dict[str, str] | None = None,
    profile: str | None = None,
    work_dir: str | None = None,
    dry_run: bool = True,
) -> dict[str, Any]:
    """Run a Nextflow workflow, or return the nextflow command in dry-run mode."""
    return run_nextflow_workflow_job(
        workflow=workflow,
        params=params,
        profile=profile,
        work_dir=work_dir,
        dry_run=dry_run,
    )


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
