---
name: validate-vcf-to-parquet
description: Parse and validate VCF files before converting them to parquet. Use when Codex needs to inspect a .vcf or .vcf.gz file, validate VCF headers, required field names, INFO/FORMAT/FILTER declarations, fixed-column values, INFO and FORMAT field types, or produce a parquet table plus validation report from VCF records.
---

# Validate VCF To Parquet

## Quick Start

Use the bundled parser for deterministic validation and conversion:

```bash
python .codex/skills/validate-vcf-to-parquet/scripts/vcf_to_parquet.py input.vcf output.parquet --report validation.json
```

For gzipped VCFs, pass the `.vcf.gz` path directly.

## Workflow

1. Inspect the user's requested VCF path and choose an output parquet path.
2. Run `scripts/vcf_to_parquet.py`.
3. Read the JSON report before declaring success.
4. If `valid` is false, report the blocking errors with line numbers and do not treat the parquet as analysis-ready.
5. If `valid` is true but warnings exist, summarize the warnings and provide the parquet path.

## What The Script Validates

The script validates:

- exactly one `#CHROM` header line
- required fixed columns in this exact order: `#CHROM POS ID REF ALT QUAL FILTER INFO`
- `FORMAT` presence when sample columns exist
- duplicate field names in the header
- `##INFO`, `##FORMAT`, and `##FILTER` declaration syntax
- `POS` as positive integer
- `QUAL` as float or `.`
- `INFO` IDs against header declarations, including basic `Number` and `Type` checks
- `FORMAT` IDs and sample values against header declarations, including basic `Number` and `Type` checks
- `FILTER` values against declared filters, warning on undeclared filters by default

See `references/vcf-validation.md` before changing validation behavior or interpreting edge cases.

## Output Table

The parquet output contains one row per VCF record. Multi-allelic records remain a single row with comma-separated `alt` so the table preserves the input record shape.

Columns:

- `chrom`, `pos`, `id`, `ref`, `alt`, `qual`, `filter`
- `info_json`: parsed INFO fields as JSON
- `format`: FORMAT column or null
- `samples_json`: sample values keyed by sample name as JSON
- `source_line`: 1-based input line number

## Common Options

- Use `--allow-undeclared-info` for permissive files with INFO fields not declared in the header.
- Use `--allow-undeclared-format` for permissive files with FORMAT fields not declared in the header.
- Use `--strict-filter` to make undeclared FILTER values errors instead of warnings.

Prefer strict validation unless the user explicitly wants a best-effort conversion.
