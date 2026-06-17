from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
EXAMPLES_DIR = ROOT / "examples"


VARIANTS = [
    ("1", 10177, "A", "AC", "RCV000000001", "Benign", "ACVRL1", "not provided"),
    ("1", 10352, "T", "TA", "RCV000000002", "Likely benign", "WASH7P", "not specified"),
    ("1", 10616, "CCGCCGTTGCAAAGGCGCGCCG", "C", "RCV000000003", "Uncertain significance", "DDX11L1", "developmental disorder"),
    ("1", 11008, "C", "G", "RCV000000004", "Pathogenic", "OR4F5", "congenital anomaly"),
    ("1", 13110, "G", "A", "RCV000000005", "Likely pathogenic", "OR4F5", "cardiomyopathy"),
    ("2", 20001, "A", "G", "RCV000000006", "Benign", "BRCA2", "hereditary cancer"),
    ("2", 20045, "G", "T", "RCV000000007", "Pathogenic", "BRCA2", "breast-ovarian cancer"),
    ("3", 30017, "T", "C", "RCV000000008", "Likely benign", "MLH1", "Lynch syndrome"),
    ("3", 30122, "C", "CT", "RCV000000009", "Uncertain significance", "MLH1", "colorectal cancer"),
    ("4", 40009, "G", "A", "RCV000000010", "Pathogenic", "HTT", "Huntington disease"),
    ("5", 50033, "A", "T", "RCV000000011", "Benign", "APC", "familial adenomatous polyposis"),
    ("6", 60044, "C", "A", "RCV000000012", "Likely pathogenic", "HFE", "hemochromatosis"),
    ("7", 70055, "G", "C", "RCV000000013", "Pathogenic", "CFTR", "cystic fibrosis"),
    ("8", 80066, "T", "G", "RCV000000014", "Benign", "WRN", "Werner syndrome"),
    ("9", 90077, "A", "C", "RCV000000015", "Uncertain significance", "TSC1", "tuberous sclerosis"),
    ("10", 100088, "G", "GA", "RCV000000016", "Likely pathogenic", "RET", "multiple endocrine neoplasia"),
    ("11", 110099, "C", "T", "RCV000000017", "Pathogenic", "HBB", "sickle cell disease"),
    ("12", 120111, "T", "A", "RCV000000018", "Likely benign", "PAH", "phenylketonuria"),
    ("13", 130123, "G", "T", "RCV000000019", "Benign", "RB1", "retinoblastoma"),
    ("14", 140135, "A", "G", "RCV000000020", "Pathogenic", "SERPINA1", "alpha-1 antitrypsin deficiency"),
]

GNOMAD_VARIANTS = [
    ("1", 10177, "A", "AC", 0.021, 1280, 61000, "NFE:0.018;AFR:0.030"),
    ("1", 10352, "T", "TA", 0.015, 940, 61500, "NFE:0.010;AMR:0.020"),
    ("1", 10616, "CCGCCGTTGCAAAGGCGCGCCG", "C", 0.0004, 24, 60012, "SAS:0.001"),
    ("1", 11008, "C", "G", 0.0, 0, 60200, ""),
    ("1", 13110, "G", "A", 0.00002, 1, 60150, "NFE:0.00003"),
    ("2", 20001, "A", "G", 0.12, 7200, 60000, "NFE:0.130;AFR:0.100"),
    ("2", 20045, "G", "T", 0.00001, 1, 61200, "EAS:0.00002"),
    ("3", 30017, "T", "C", 0.008, 480, 60220, "AMR:0.011"),
    ("3", 30122, "C", "CT", 0.0007, 44, 60310, "NFE:0.0008"),
    ("4", 40009, "G", "A", 0.0, 0, 62000, ""),
    ("5", 50033, "A", "T", 0.045, 2700, 60040, "AFR:0.060"),
    ("6", 60044, "C", "A", 0.00003, 2, 60100, "NFE:0.00005"),
    ("7", 70055, "G", "C", 0.00001, 1, 60090, "NFE:0.00002"),
    ("8", 80066, "T", "G", 0.032, 1920, 60000, "NFE:0.029"),
    ("9", 90077, "A", "C", 0.0012, 72, 60200, "SAS:0.002"),
    ("10", 100088, "G", "GA", 0.0, 0, 60050, ""),
    ("11", 110099, "C", "T", 0.0002, 12, 60400, "AFR:0.001"),
    ("12", 120111, "T", "A", 0.009, 545, 60200, "NFE:0.008"),
    ("13", 130123, "G", "T", 0.018, 1090, 60500, "AMR:0.022"),
    ("14", 140135, "A", "G", 0.00004, 3, 61000, "NFE:0.00006"),
]

OMIM_GENES = [
    ("ACVRL1", "600376", "Hereditary hemorrhagic telangiectasia", "AD"),
    ("WASH7P", "", "No curated OMIM phenotype in synthetic data", ""),
    ("DDX11L1", "", "No curated OMIM phenotype in synthetic data", ""),
    ("OR4F5", "608565", "Synthetic congenital anomaly susceptibility", "AD"),
    ("BRCA2", "600185", "Breast-ovarian cancer, familial, 2", "AD"),
    ("MLH1", "120436", "Lynch syndrome 2", "AD"),
    ("HTT", "143100", "Huntington disease", "AD"),
    ("APC", "611731", "Familial adenomatous polyposis 1", "AD"),
    ("HFE", "613609", "Hemochromatosis type 1", "AR"),
    ("CFTR", "602421", "Cystic fibrosis", "AR"),
    ("WRN", "604611", "Werner syndrome", "AR"),
    ("TSC1", "605284", "Tuberous sclerosis 1", "AD"),
    ("RET", "164761", "Multiple endocrine neoplasia IIA", "AD"),
    ("HBB", "141900", "Sickle cell disease", "AR"),
    ("PAH", "612349", "Phenylketonuria", "AR"),
    ("RB1", "180200", "Retinoblastoma", "AD"),
    ("SERPINA1", "613490", "Alpha-1 antitrypsin deficiency", "AR"),
]


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)
    EXAMPLES_DIR.mkdir(exist_ok=True)

    frame = pd.DataFrame(
        VARIANTS,
        columns=[
            "chrom",
            "pos",
            "ref",
            "alt",
            "variation_id",
            "clinical_significance",
            "gene",
            "condition",
        ],
    )
    frame.to_parquet(DATA_DIR / "synthetic_clinvar.parquet", index=False)

    gnomad = pd.DataFrame(
        GNOMAD_VARIANTS,
        columns=[
            "chrom",
            "pos",
            "ref",
            "alt",
            "allele_frequency",
            "allele_count",
            "allele_number",
            "popmax",
        ],
    )
    gnomad.to_parquet(DATA_DIR / "synthetic_gnomad.parquet", index=False)

    omim = pd.DataFrame(
        OMIM_GENES,
        columns=["gene", "mim_number", "phenotype", "inheritance"],
    )
    omim.to_parquet(DATA_DIR / "synthetic_omim.parquet", index=False)

    sample_rows = VARIANTS[:7] + [
        ("3", 30017, "T", "C,G", ".", ".", ".", "."),
        ("15", 150001, "C", "G", ".", ".", ".", "."),
        ("16", 160002, "G", "A", ".", ".", ".", "."),
    ]
    vcf_lines = [
        "##fileformat=VCFv4.2",
        "##source=skilled-mcp-var-annotation",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for chrom, pos, ref, alt, *_ in sample_rows:
        vcf_lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.")

    (EXAMPLES_DIR / "sample.vcf").write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
