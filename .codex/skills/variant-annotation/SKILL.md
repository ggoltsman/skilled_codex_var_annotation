---
name: variant-annotation
description: Use the local synthetic ClinVar parquet through the variant-annotation MCP server to annotate VCF variants or individual chrom/pos/ref/alt variant queries.
---

# Variant Annotation

Use this skill when the user asks to annotate a genomic variant, inspect the sample VCF, or explain how the local variant annotation MCP server works.

## Data Source

The project contains a synthetic ClinVar-style parquet file at:

`data/synthetic_clinvar.parquet`

It has 20 records with these fields:

- `chrom`
- `pos`
- `ref`
- `alt`
- `variation_id`
- `clinical_significance`
- `gene`
- `condition`

## MCP Tool

Prefer the MCP server tool when it is available:

```python
normalize_vcf(input_vcf: str, output_vcf: str | None = None)
query_variant(chrom: str, pos: int, ref: str, alt: str)
summarize_annotations(
    input_vcf: str,
    output_tsv: str | None = None,
    report_path: str | None = None,
    normalize: bool = True,
)
```

Call it with VCF coordinates exactly as they appear in the VCF row. Chromosomes may be passed with or without a `chr` prefix.

## Workflow

1. Inspect the input VCF.
2. Normalize variant representation with `normalize_vcf`. It defaults to `bcftools norm -m -any` first. If `bcftools` is not found, warn that normalization will proceed using a custom Python function for multiallelic splitting.
3. Use MCP tools when available:
   - normalize_vcf
   - query_variant
   - summarize_annotations
4. Prefer local reference data over public web queries.
5. Produce:
   - annotated TSV
   - candidate pathogenic/likely pathogenic variant report
6. For missing variants, state that the synthetic ClinVar parquet has no matching record.

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
