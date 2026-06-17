from __future__ import annotations

import argparse
import gzip
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


REQUIRED_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
SUPPORTED_TYPES = {"Integer", "Float", "Flag", "Character", "String"}
META_RE = re.compile(r"^##(?P<section>INFO|FORMAT|FILTER)=<(?P<body>.*)>$")


@dataclass
class FieldDef:
    id: str
    number: str | None
    type: str | None
    description: str | None


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def split_meta_body(body: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    in_quotes = False
    escaped = False
    for char in body:
        if escaped:
            current.append(char)
            escaped = False
            continue
        if char == "\\":
            current.append(char)
            escaped = True
            continue
        if char == '"':
            in_quotes = not in_quotes
            current.append(char)
            continue
        if char == "," and not in_quotes:
            parts.append("".join(current))
            current = []
            continue
        current.append(char)
    parts.append("".join(current))
    return parts


def parse_meta_definition(line: str) -> tuple[str, FieldDef] | None:
    match = META_RE.match(line.rstrip("\n"))
    if not match:
        return None
    values: dict[str, str] = {}
    for part in split_meta_body(match.group("body")):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        values[key] = value.strip().strip('"')
    field_id = values.get("ID")
    if not field_id:
        return None
    return (
        match.group("section"),
        FieldDef(
            id=field_id,
            number=values.get("Number"),
            type=values.get("Type"),
            description=values.get("Description"),
        ),
    )


def validate_type(value: str, field_type: str) -> bool:
    if value == ".":
        return True
    if field_type == "Integer":
        try:
            int(value)
            return True
        except ValueError:
            return False
    if field_type == "Float":
        try:
            float(value)
            return True
        except ValueError:
            return False
    if field_type == "Character":
        return len(value) == 1
    return True


def expected_count(number: str | None, alt_count: int) -> int | None:
    if number is None or number in {".", "G"}:
        return None
    if number == "A":
        return alt_count
    if number == "R":
        return alt_count + 1
    try:
        return int(number)
    except ValueError:
        return None


def validate_values(
    value: str | None,
    field: FieldDef,
    alt_count: int,
    line_no: int,
    label: str,
    errors: list[str],
) -> None:
    if field.type not in SUPPORTED_TYPES:
        errors.append(f"line {line_no}: {label} {field.id} has unsupported declared Type={field.type}")
        return
    if field.type == "Flag":
        if value is not None:
            errors.append(f"line {line_no}: Flag field {label} {field.id} must not have a value")
        return
    if value is None:
        errors.append(f"line {line_no}: {label} {field.id} requires a value")
        return

    values = [] if value == "" else value.split(",")
    count = expected_count(field.number, alt_count)
    if count is not None and len(values) != count:
        errors.append(
            f"line {line_no}: {label} {field.id} expected {count} value(s) from Number={field.number}, got {len(values)}"
        )
    for item in values:
        if not validate_type(item, field.type or "String"):
            errors.append(f"line {line_no}: {label} {field.id} value {item!r} does not match Type={field.type}")


def parse_info(
    value: str,
    info_defs: dict[str, FieldDef],
    alt_count: int,
    line_no: int,
    allow_undeclared: bool,
    errors: list[str],
) -> dict[str, str | bool]:
    if value == ".":
        return {}
    parsed: dict[str, str | bool] = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" in item:
            key, raw_value = item.split("=", 1)
            parsed[key] = raw_value
        else:
            key, raw_value = item, None
            parsed[key] = True
        field = info_defs.get(key)
        if field is None:
            if not allow_undeclared:
                errors.append(f"line {line_no}: INFO field {key} is not declared in the header")
            continue
        validate_values(raw_value, field, alt_count, line_no, "INFO", errors)
    return parsed


def validate_format_and_samples(
    format_value: str | None,
    sample_values: list[str],
    sample_names: list[str],
    format_defs: dict[str, FieldDef],
    alt_count: int,
    line_no: int,
    allow_undeclared: bool,
    errors: list[str],
) -> dict[str, dict[str, str | None]]:
    if not sample_names:
        return {}
    if not format_value or format_value == ".":
        errors.append(f"line {line_no}: sample columns are present but FORMAT is missing")
        return {}

    keys = format_value.split(":")
    for key in keys:
        field = format_defs.get(key)
        if field is None:
            if key != "GT" and not allow_undeclared:
                errors.append(f"line {line_no}: FORMAT field {key} is not declared in the header")
            continue

    samples: dict[str, dict[str, str | None]] = {}
    for sample_name, raw_sample in zip(sample_names, sample_values):
        values = raw_sample.split(":")
        sample = {key: values[index] if index < len(values) else None for index, key in enumerate(keys)}
        samples[sample_name] = sample
        for key, raw_value in sample.items():
            field = format_defs.get(key)
            if field is None or key == "GT":
                continue
            validate_values(raw_value, field, alt_count, line_no, f"FORMAT {sample_name}", errors)
    return samples


def parse_vcf_to_rows(
    input_vcf: Path,
    allow_undeclared_info: bool,
    allow_undeclared_format: bool,
    strict_filter: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    errors: list[str] = []
    warnings: list[str] = []
    info_defs: dict[str, FieldDef] = {}
    format_defs: dict[str, FieldDef] = {}
    filter_defs: dict[str, FieldDef] = {}
    header: list[str] | None = None
    rows: list[dict[str, Any]] = []

    with open_text(input_vcf) as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                parsed_meta = parse_meta_definition(line)
                if parsed_meta:
                    section, field = parsed_meta
                    if section == "INFO":
                        info_defs[field.id] = field
                    elif section == "FORMAT":
                        format_defs[field.id] = field
                    elif section == "FILTER":
                        filter_defs[field.id] = field
                continue
            if line.startswith("#"):
                if header is not None:
                    errors.append(f"line {line_no}: duplicate #CHROM header")
                    continue
                header = line.split("\t")
                if header[: len(REQUIRED_COLUMNS)] != REQUIRED_COLUMNS:
                    errors.append(f"line {line_no}: header must start with {' '.join(REQUIRED_COLUMNS)}")
                if len(header) != len(set(header)):
                    errors.append(f"line {line_no}: header contains duplicate column names")
                if len(header) > 8 and header[8] != "FORMAT":
                    errors.append(f"line {line_no}: FORMAT column is required before sample columns")
                continue

            if header is None:
                errors.append(f"line {line_no}: data row appears before #CHROM header")
                continue

            fields = line.split("\t")
            if len(fields) != len(header):
                errors.append(f"line {line_no}: expected {len(header)} columns, got {len(fields)}")
                continue

            chrom, pos, variant_id, ref, alt, qual, filter_value, info = fields[:8]
            alt_values = alt.split(",") if alt not in {"", "."} else []
            alt_count = len(alt_values)
            try:
                pos_value = int(pos)
                if pos_value <= 0:
                    errors.append(f"line {line_no}: POS must be a positive integer")
            except ValueError:
                pos_value = pd.NA
                errors.append(f"line {line_no}: POS {pos!r} is not an integer")
            if qual != ".":
                try:
                    float(qual)
                except ValueError:
                    errors.append(f"line {line_no}: QUAL {qual!r} is not a float or '.'")
            if filter_value not in {".", "PASS"}:
                for filter_id in filter_value.split(";"):
                    if filter_id not in filter_defs:
                        message = f"line {line_no}: FILTER value {filter_id} is not declared in the header"
                        if strict_filter:
                            errors.append(message)
                        else:
                            warnings.append(message)

            info_json = parse_info(info, info_defs, alt_count, line_no, allow_undeclared_info, errors)
            sample_names = header[9:] if len(header) > 9 else []
            format_value = fields[8] if len(fields) > 8 else None
            sample_values = fields[9:] if len(fields) > 9 else []
            samples_json = validate_format_and_samples(
                format_value,
                sample_values,
                sample_names,
                format_defs,
                alt_count,
                line_no,
                allow_undeclared_format,
                errors,
            )
            rows.append(
                {
                    "chrom": chrom,
                    "pos": pos_value,
                    "id": variant_id,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "filter": filter_value,
                    "info_json": json.dumps(info_json, sort_keys=True),
                    "format": format_value,
                    "samples_json": json.dumps(samples_json, sort_keys=True),
                    "source_line": line_no,
                }
            )

    if header is None:
        errors.append("missing #CHROM header")

    report = {
        "input_vcf": str(input_vcf),
        "valid": not errors,
        "error_count": len(errors),
        "warning_count": len(warnings),
        "errors": errors,
        "warnings": warnings,
        "record_count": len(rows),
        "sample_count": max(len(header or []) - 9, 0),
        "samples": (header or [])[9:],
        "declared_info_fields": sorted(info_defs),
        "declared_format_fields": sorted(format_defs),
        "declared_filter_fields": sorted(filter_defs),
    }
    return rows, report


def write_parquet(rows: list[dict[str, Any]], output_parquet: Path) -> None:
    output_parquet.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    columns = [
        "chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "qual",
        "filter",
        "info_json",
        "format",
        "samples_json",
        "source_line",
    ]
    for column in columns:
        if column not in frame.columns:
            frame[column] = pd.NA
    frame = frame[columns]
    frame["chrom"] = frame["chrom"].astype("string")
    frame["pos"] = frame["pos"].astype("Int64")
    frame["id"] = frame["id"].astype("string")
    frame["ref"] = frame["ref"].astype("string")
    frame["alt"] = frame["alt"].astype("string")
    frame["qual"] = frame["qual"].astype("string")
    frame["filter"] = frame["filter"].astype("string")
    frame["info_json"] = frame["info_json"].astype("string")
    frame["format"] = frame["format"].astype("string")
    frame["samples_json"] = frame["samples_json"].astype("string")
    frame["source_line"] = frame["source_line"].astype("Int64")
    frame.to_parquet(output_parquet, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a VCF and convert records to parquet.")
    parser.add_argument("input_vcf", type=Path)
    parser.add_argument("output_parquet", type=Path)
    parser.add_argument("--report", type=Path, help="Path for JSON validation report.")
    parser.add_argument("--allow-undeclared-info", action="store_true")
    parser.add_argument("--allow-undeclared-format", action="store_true")
    parser.add_argument("--strict-filter", action="store_true")
    args = parser.parse_args()

    rows, report = parse_vcf_to_rows(
        input_vcf=args.input_vcf,
        allow_undeclared_info=args.allow_undeclared_info,
        allow_undeclared_format=args.allow_undeclared_format,
        strict_filter=args.strict_filter,
    )
    write_parquet(rows, args.output_parquet)
    report["output_parquet"] = str(args.output_parquet)

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    else:
        print(json.dumps(report, indent=2))

    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
