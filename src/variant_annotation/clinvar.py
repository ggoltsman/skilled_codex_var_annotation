from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_CLINVAR_PATH = Path(__file__).resolve().parents[2] / "data" / "synthetic_clinvar.parquet"


@dataclass(frozen=True)
class ClinVarStore:
    """Tiny parquet-backed lookup store keyed by chrom, pos, ref, and alt."""

    parquet_path: Path = DEFAULT_CLINVAR_PATH

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
                "query": {
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                },
                "message": "Variant is not present in the synthetic ClinVar parquet.",
            }

        record = matches.iloc[0].to_dict()
        record["pos"] = int(record["pos"])
        record["found"] = True
        return record
