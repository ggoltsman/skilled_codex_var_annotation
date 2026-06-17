# VCF Validation Contract

Use this reference when interpreting or modifying `scripts/vcf_to_parquet.py`.

## Header Rules

- Metadata lines must start with `##`.
- The column header must start with the required fixed VCF columns in this exact order:
  `#CHROM POS ID REF ALT QUAL FILTER INFO`.
- Sample columns are allowed only after `FORMAT`.
- If sample columns exist, `FORMAT` is required as column 9.
- Header column names must be unique.

## Metadata Declarations

The parser recognizes declarations shaped like:

```text
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FILTER=<ID=q10,Description="Quality below 10">
```

Supported INFO and FORMAT types:

- `Integer`
- `Float`
- `Flag`
- `Character`
- `String`

Supported Number checks:

- `0`: no value, used for `Flag`
- `1`: exactly one value
- `A`: one value per ALT allele
- `R`: one value for REF plus each ALT allele
- `.`: any count

The script intentionally treats `G` as permissive because ploidy and allele ordering require deeper genotype modeling.

## Output Shape

Parquet rows preserve VCF records rather than splitting ALT alleles. Store structured fields as JSON strings so the parquet table is easy to inspect with pandas, DuckDB, Spark, or PyArrow without nested-schema surprises.

## Error vs Warning

Use errors for malformed structure, invalid fixed-column values, declared type mismatches, and undeclared INFO/FORMAT fields under strict defaults.

Use warnings for undeclared FILTER values by default because many small VCFs omit FILTER declarations. Enable `--strict-filter` when a pipeline requires complete FILTER metadata.
