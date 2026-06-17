---
name: variant-annotation
description: Annotate VCFs with a realistic local VEP-style workflow that joins synthetic ClinVar, gnomAD, and OMIM parquet data, writes an Iceberg-like annotated table, and produces pathogenic reports plus QC summaries through the variant-annotation MCP server.
---

# Variant Annotation

Use this skill when the user asks to annotate a VCF, inspect variant evidence, query local ClinVar/gnomAD/OMIM data, submit a Spark annotation job, run a Nextflow workflow, or explain how the local variant annotation MCP server works.

## Data Sources

The project contains local synthetic parquet files at:

- `data/synthetic_clinvar.parquet`
- `data/synthetic_gnomad.parquet`
- `data/synthetic_omim.parquet`

ClinVar is keyed by variant and includes:

- `chrom`
- `pos`
- `ref`
- `alt`
- `variation_id`
- `clinical_significance`
- `gene`
- `condition`

gnomAD is keyed by variant and includes:

- `chrom`
- `pos`
- `ref`
- `alt`
- `allele_frequency`
- `allele_count`
- `allele_number`
- `popmax`

OMIM is keyed by gene and includes:

- `gene`
- `mim_number`
- `phenotype`
- `inheritance`

## Output Table Schema

The annotated Iceberg-like table schema is defined outside the workflow code at:

`src/variant_annotation/schemas/annotated_variant_table.schema.json`

`annotate_vcf()` loads this schema when writing `warehouse/*_annotations/data.parquet`. The schema controls column order, logical types, nullability, and field descriptions. The generated `metadata.json` records the schema name, version, source path, and SHA-256 checksum for traceability.

## MCP Tool

Prefer the MCP server tool when it is available:

```python
normalize_vcf(input_vcf: str, output_vcf: str | None = None)
query_variant(chrom: str, pos: int, ref: str, alt: str)  # backwards-compatible ClinVar lookup
query_clinvar(chrom: str, pos: int, ref: str, alt: str)
query_gnomad(chrom: str, pos: int, ref: str, alt: str)
query_omim(gene: str)
run_vep_annotation(
    input_vcf: str,
    output_tsv: str | None = None,
    cache_dir: str | None = None,
    offline: bool = True,
)
annotate_vcf(
    input_vcf: str,
    output_table: str | None = None,
    pathogenic_report: str | None = None,
    qc_summary: str | None = None,
    normalize: bool = True,
    run_vep_step: bool = True,
)
summarize_annotations(
    input_vcf: str,
    output_tsv: str | None = None,
    report_path: str | None = None,
    normalize: bool = True,
)
submit_spark_job(
    app_path: str,
    args: list[str] | None = None,
    master: str = "local[*]",
    dry_run: bool = True,
)
run_nextflow_workflow(
    workflow: str,
    params: dict[str, str] | None = None,
    profile: str | None = None,
    work_dir: str | None = None,
    dry_run: bool = True,
)
```

Call it with VCF coordinates exactly as they appear in the VCF row. Chromosomes may be passed with or without a `chr` prefix.

`annotate_vcf()` is the preferred end-to-end workflow. It normalizes the VCF, runs VEP if the `vep` executable is available, falls back to deterministic VEP-like consequence labels when VEP is unavailable, joins local ClinVar/gnomAD/OMIM parquet data, writes an annotated Iceberg-like parquet table under `warehouse/`, writes a pathogenic/likely pathogenic Markdown report, and writes a QC JSON summary.

`summarize_annotations()` is retained as the minimal legacy ClinVar-only workflow.

## Workflow

1. Inspect the input VCF.
2. Normalize variant representation with `normalize_vcf`. It defaults to `bcftools norm -m -any` first. If `bcftools` is not found, warn that normalization will proceed using a custom Python function for multiallelic splitting.
3. Use MCP tools when available:
   - normalize_vcf
   - run_vep_annotation
   - query_clinvar
   - query_gnomad
   - query_omim
   - annotate_vcf
   - submit_spark_job
   - run_nextflow_workflow
4. Prefer local reference data over public web queries.
5. Produce:
   - annotated Iceberg-like parquet table
   - candidate pathogenic/likely pathogenic variant report
   - QC summary
6. For missing variants, state which local synthetic parquet source has no matching record.

For production-shaped requests, use `submit_spark_job()` to launch an existing Spark app that writes to a real Iceberg catalog. Use `dry_run=True` unless the user explicitly wants to execute the job. Use `run_nextflow_workflow()` when the user asks to run a workflow wrapper around normalization, VEP, annotation, Spark table creation, and report generation.

## Example

VCF row:

```text
1	11008	.	C	G	.	PASS	.
```

Tool call:

```python
query_variant(chrom="1", pos=11008, ref="C", alt="G")
```

Expected result:
1       11008	C	G	a synthetic `Pathogenic` annotation for `OR4F5`.

End-to-end local annotation:

```python
annotate_vcf(input_vcf="examples/sample.vcf")
```

Expected outputs:

- `warehouse/sample_annotations/data.parquet`
- `warehouse/sample_annotations/metadata.json`
- `outputs/sample.pathogenic_likely_pathogenic.md`
- `outputs/sample.qc_summary.json`
