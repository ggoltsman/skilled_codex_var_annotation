# Skilled MCP Variant Annotation

This is a minimal Codex bioinformatics project that wires together:

- a Codex skill named `variant-annotation`
- a Python MCP server exposing `query_variant()`, `normalize_vcf()`, and `summarize_annotations()`
- a synthetic ClinVar-style parquet file with 20 variants
- a sample VCF containing variants that hit and miss the synthetic database

## Layout

```text
.codex/skills/variant-annotation/SKILL.md  Codex skill instructions
src/variant_annotation/                   Python package and MCP server
scripts/generate_synthetic_data.py         Creates the parquet and sample VCF
scripts/smoke_test_mcp.py                  Verifies the MCP tool over stdio
data/synthetic_clinvar.parquet             Synthetic ClinVar records
examples/sample.vcf                        Small VCF for testing
outputs/                                  Generated normalized VCFs, TSVs, and reports
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
query_variant(chrom: str, pos: int, ref: str, alt: str) -> dict
normalize_vcf(input_vcf: str, output_vcf: str | None = None) -> dict
summarize_annotations(
    input_vcf: str,
    output_tsv: str | None = None,
    report_path: str | None = None,
    normalize: bool = True,
) -> dict
```

`query_variant()` returns `found: true` plus ClinVar-style fields when the variant exists in `data/synthetic_clinvar.parquet`; otherwise it returns `found: false` and the original query.

`normalize_vcf()` writes a normalized VCF under `outputs/` by default.

`summarize_annotations()` normalizes a VCF, annotates every row against the synthetic ClinVar parquet, writes an annotated TSV, and writes a Markdown report containing pathogenic or likely pathogenic candidates.

## How The Pieces Fit

The skill teaches Codex when and how to use the local variant annotation MCP tools. The MCP server loads the synthetic parquet file, exposes a single-variant lookup, and provides a small batch workflow for VCF normalization plus annotation summaries. The sample VCF provides realistic input rows that can be queried one at a time through `query_variant()` or processed as a batch through `summarize_annotations()`.

For example, this VCF row:

```text
1	11008	.	C	G	.	PASS	.
```

maps to:

```python
query_variant(chrom="1", pos=11008, ref="C", alt="G")
```

which returns a synthetic pathogenic `OR4F5` annotation.
