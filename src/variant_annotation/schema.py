from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from importlib import resources
from typing import Any

import pandas as pd


ANNOTATED_TABLE_SCHEMA_RESOURCE = "schemas/annotated_variant_table.schema.json"
PANDAS_DTYPES = {
    "boolean": "boolean",
    "double": "Float64",
    "int64": "Int64",
    "string": "string",
}
MAX_COERCION_ERROR_VALUES = 5


@dataclass(frozen=True)
class TableField:
    name: str
    type: str
    nullable: bool
    description: str


@dataclass(frozen=True)
class TableSchema:
    name: str
    version: str
    description: str
    fields: list[TableField]
    source: str
    sha256: str

    @property
    def column_names(self) -> list[str]:
        return [field.name for field in self.fields]

    def to_metadata(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "source": self.source,
            "sha256": self.sha256,
            "fields": [
                {
                    "name": field.name,
                    "type": field.type,
                    "nullable": field.nullable,
                    "description": field.description,
                }
                for field in self.fields
            ],
        }


def load_annotated_table_schema() -> TableSchema:
    package = resources.files("variant_annotation")
    schema_file = package.joinpath(ANNOTATED_TABLE_SCHEMA_RESOURCE)
    raw = schema_file.read_bytes()
    payload = json.loads(raw.decode("utf-8"))
    return TableSchema(
        name=payload["name"],
        version=payload["version"],
        description=payload["description"],
        fields=[TableField(**field) for field in payload["fields"]],
        source=f"variant_annotation/{ANNOTATED_TABLE_SCHEMA_RESOURCE}",
        sha256=hashlib.sha256(raw).hexdigest(),
    )


def _format_sample_values(values: pd.Series) -> str:
    sample = values.drop_duplicates().head(MAX_COERCION_ERROR_VALUES).tolist()
    return ", ".join(repr(value) for value in sample)


def apply_table_schema(rows: list[dict[str, Any]], schema: TableSchema) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    for field in schema.fields:
        if field.name not in frame.columns:
            frame[field.name] = pd.NA

    frame = frame[schema.column_names]
    for field in schema.fields:
        try:
            dtype = PANDAS_DTYPES[field.type]
        except KeyError as exc:
            raise ValueError(
                f"Annotated table schema field {field.name!r} has unsupported type "
                f"{field.type!r}."
            ) from exc

        try:
            frame[field.name] = frame[field.name].astype(dtype)
        except (TypeError, ValueError) as exc:
            values = _format_sample_values(frame[field.name].dropna())
            sample_message = f" Sample values: {values}." if values else ""
            raise ValueError(
                f"Annotated table column {field.name!r} cannot be coerced from "
                f"{frame[field.name].dtype} to schema type {field.type!r} "
                f"(pandas dtype {dtype!r}).{sample_message}"
            ) from exc

    if not frame.empty:
        required_nulls = [
            field.name
            for field in schema.fields
            if not field.nullable and frame[field.name].isna().any()
        ]
        if required_nulls:
            columns = ", ".join(required_nulls)
            raise ValueError(f"Annotated table has null values in non-nullable columns: {columns}")
    return frame
