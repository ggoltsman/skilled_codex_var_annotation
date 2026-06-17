# Skilled MCP Variant Annotation

This is a Codex bioinformatics project that wires together:

- a Codex skill named `variant-annotation`
- a Python MCP server for local VCF annotation workflows
- synthetic ClinVar, gnomAD, and OMIM parquet files
- a sample VCF containing variants that hit and miss the synthetic database
- local outputs shaped like a production workflow: an annotated Iceberg-like parquet table, a pathogenic/likely pathogenic report, and a QC summary

## Layout

```text
.codex/skills/variant-annotation/SKILL.md  Codex skill instructions
src/variant_annotation/                   Python package and MCP server
scripts/generate_synthetic_data.py         Creates synthetic parquet data and sample VCF
scripts/smoke_test_mcp.py                  Verifies the MCP tool over stdio
data/synthetic_clinvar.parquet             Synthetic ClinVar records keyed by variant
data/synthetic_gnomad.parquet              Synthetic gnomAD population frequencies keyed by variant
data/synthetic_omim.parquet                Synthetic OMIM disease records keyed by gene
examples/sample.vcf                        Small VCF for testing
outputs/                                  Generated normalized VCFs, VEP TSVs, reports, and QC summaries
warehouse/                                Generated Iceberg-like local parquet annotation tables
requirements.txt                           Runtime dependencies
```

## Setup

Dependencies are installed into the project-local virtual environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/generate_synthetic_data.py
```

`normalize_vcf()` defaults to `bcftools norm -m -any` first. If `bcftools` is not found, it warns and proceeds with a custom Python function that splits comma-separated ALT alleles without reference-based left alignment.

Smoke-test the MCP tool:

```bash
.venv/bin/python scripts/smoke_test_mcp.py
```

## MCP Server

Run the server over stdio:

```bash
.venv/bin/python -m variant_annotation.server
```

Example Codex MCP configuration:

```json
{
  "mcpServers": {
    "variant-annotation": {
      "command": "/Users/gee/Training/Codex/skilled_mcp_var_annotation/.venv/bin/python",
      "args": ["-m", "variant_annotation.server"],
      "cwd": "/Users/gee/Training/Codex/skilled_mcp_var_annotation"
    }
  }
}
```

The exposed tools are:

```python
normalize_vcf(input_vcf: str, output_vcf: str | None = None) -> dict
query_variant(chrom: str, pos: int, ref: str, alt: str) -> dict
query_clinvar(chrom: str, pos: int, ref: str, alt: str) -> dict
query_gnomad(chrom: str, pos: int, ref: str, alt: str) -> dict
query_omim(gene: str) -> dict
run_vep_annotation(
    input_vcf: str,
    output_tsv: str | None = None,
    cache_dir: str | None = None,
    offline: bool = True,
) -> dict
annotate_vcf(
    input_vcf: str,
    output_table: str | None = None,
    pathogenic_report: str | None = None,
    qc_summary: str | None = None,
    normalize: bool = True,
    run_vep_step: bool = True,
) -> dict
summarize_annotations(
    input_vcf: str,
    output_tsv: str | None = None,
    report_path: str | None = None,
    normalize: bool = True,
) -> dict
submit_spark_job(
    app_path: str,
    args: list[str] | None = None,
    master: str = "local[*]",
    dry_run: bool = True,
) -> dict
run_nextflow_workflow(
    workflow: str,
    params: dict[str, str] | None = None,
    profile: str | None = None,
    work_dir: str | None = None,
    dry_run: bool = True,
) -> dict
```

`query_variant()` is a backwards-compatible ClinVar lookup. Prefer `query_clinvar()`, `query_gnomad()`, and `query_omim()` for explicit source-specific lookups.

`normalize_vcf()` writes a normalized VCF under `outputs/` by default.

`run_vep_annotation()` runs the `vep` executable when it is installed. If VEP is unavailable, it writes deterministic VEP-like consequence annotations so the local workflow remains testable.

`annotate_vcf()` is the preferred end-to-end workflow. It normalizes a VCF, runs the VEP step, joins synthetic ClinVar, gnomAD, and OMIM parquet data, writes an Iceberg-like local parquet table, writes a Markdown pathogenic/likely pathogenic report, and writes a JSON QC summary.

`summarize_annotations()` is retained as the minimal legacy ClinVar-only workflow.

`submit_spark_job()` and `run_nextflow_workflow()` expose production-shaped orchestration hooks. Both default to `dry_run=True` and return the command they would run unless execution is explicitly requested and the executable is available.

## Example Workflow

Run the realistic local annotation workflow through MCP or directly in Python:

```python
from variant_annotation.annotation import annotate_vcf_realistic

annotate_vcf_realistic("examples/sample.vcf")
```

Expected outputs:

```text
outputs/sample.normalized.vcf
outputs/sample.normalized.vep.tsv
outputs/sample.pathogenic_likely_pathogenic.md
outputs/sample.qc_summary.json
warehouse/sample_annotations/data.parquet
warehouse/sample_annotations/metadata.json
```

## How The Pieces Fit

The skill teaches Codex when and how to use the local variant annotation MCP tools. The MCP server loads synthetic parquet files, exposes source-specific lookups, and provides a batch workflow for normalization, VEP-style consequence annotation, ClinVar/gnomAD/OMIM joins, table creation, pathogenic reporting, and QC. The sample VCF provides realistic input rows that can be queried one at a time or processed as a batch through `annotate_vcf()`.

For example, this VCF row:

```text
1	11008	.	C	G	.	PASS	.
```

maps to:

```python
query_variant(chrom="1", pos=11008, ref="C", alt="G")
```

which returns a synthetic pathogenic `OR4F5` annotation.
