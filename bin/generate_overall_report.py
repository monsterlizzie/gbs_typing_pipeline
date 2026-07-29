#!/usr/bin/env python3

import sys
import os
import re
import glob
import pandas as pd

from itertools import chain
from typing import List


# ---------- Columns to include ----------
COLUMNS_BY_CATEGORY = {
    'IDENTIFICATION': ['Sample_ID'],
    'QC': [
        'Read_QC',
        'Assembly_QC',
        'Mapping_QC',
        'Taxonomy_QC',
        'Overall_QC'
    ],
    'READ': [
        'Bases'
    ],
    'ASSEMBLY': [
        'Contigs#',
        'Assembly_Length',
        'Seq_Depth'
    ],
    'MAPPING': [
        'Ref_Cov_%',
        'Het-SNP#'
    ],
    'TAXONOMY': [
        'S.agalactiae_%',
        'Top_Non-agalactiae_Species',
        'Top_Non-agalactiae_Species_%'
    ],
}

QC_FIXED = list(
    chain.from_iterable(COLUMNS_BY_CATEGORY.values())
)

PBP_COLUMNS = [
    'pbp1A',
    'pbp2B',
    'pbp2X'
]

VALID_ID = re.compile(r'^[A-Za-z0-9_.:-]+$')


def infer_sample_id_from_path(p: str) -> str:
    parent = os.path.basename(
        os.path.dirname(p)
    )

    if not VALID_ID.match(parent):
        stem = os.path.splitext(
            os.path.basename(p)
        )[0]

        parent = re.sub(
            r'_report$',
            '',
            stem
        )

    return parent


def normalise_id_to_sample_id(
    df: pd.DataFrame,
    want: str = 'Sample_ID'
) -> pd.DataFrame:

    if want in df.columns:
        return df

    for cand in [
        'sample_id',
        'Sample_id',
        'sample',
        'Sample',
        'isolate',
        'Isolate',
        'id',
        'ID'
    ]:
        if cand in df.columns:
            return df.rename(
                columns={cand: want}
            )

    return df


def read_qc_stack(qc_glob: str) -> pd.DataFrame:
    paths = sorted(
        glob.glob(qc_glob)
    )

    if not paths:
        sys.exit(
            f"[combine] No QC files matched: {qc_glob}"
        )

    dfs = []

    for p in paths:
        df = pd.read_csv(
            p,
            dtype=str
        )

        df = normalise_id_to_sample_id(
            df,
            'Sample_ID'
        )

        if 'Sample_ID' not in df.columns:
            df.insert(
                0,
                'Sample_ID',
                infer_sample_id_from_path(p)
            )

        df.replace(
            "",
            pd.NA,
            inplace=True
        )

        dfs.append(df)

    return pd.concat(
        dfs,
        ignore_index=True
    )


def read_typer(path_or_none: str) -> pd.DataFrame:
    if path_or_none == 'NONE':
        return pd.DataFrame(
            columns=['Sample_ID']
        )

    try:
        df = pd.read_csv(
            path_or_none,
            sep='\t',
            dtype=str
        )
    except Exception:
        df = pd.read_csv(
            path_or_none,
            sep=',',
            dtype=str
        )

    df = normalise_id_to_sample_id(
        df,
        'Sample_ID'
    )

    if 'Sample_ID' not in df.columns:
        sys.exit(
            "[combine] Typer table missing a sample ID column "
            "(for example, Sample_ID)"
        )

    df.replace(
        "",
        pd.NA,
        inplace=True
    )

    return df


def ensure_columns(
    df: pd.DataFrame,
    cols: List[str]
) -> None:

    for c in cols:
        if c not in df.columns:
            df[c] = pd.NA


def read_pbp(path_or_none: str) -> pd.DataFrame:
    """
    Read the combined PBP call table and convert it from long format:

        Sample_ID    PBP_gene    PBP_call
        SRR9088888   pbp1A       3
        SRR9088888   pbp2B       1
        SRR9088888   pbp2X       4

    into wide format:

        Sample_ID    pbp1A    pbp2B    pbp2X
        SRR9088888   3        1        4
    """

    output_columns = [
        'Sample_ID',
        'pbp1A',
        'pbp2B',
        'pbp2X'
    ]

    if path_or_none == 'NONE':
        return pd.DataFrame(
            columns=output_columns
        )

    if not os.path.exists(path_or_none):
        sys.exit(
            f"[combine] PBP results file does not exist: "
            f"{path_or_none}"
        )

    try:
        df = pd.read_csv(
            path_or_none,
            sep='\t',
            dtype=str
        )
    except Exception:
        df = pd.read_csv(
            path_or_none,
            sep=',',
            dtype=str
        )

    df = normalise_id_to_sample_id(
        df,
        'Sample_ID'
    )

    df.replace(
        "",
        pd.NA,
        inplace=True
    )

    if 'Sample_ID' not in df.columns:
        sys.exit(
            "[combine] PBP table missing a sample ID column "
            "(for example, Sample_ID or Sample_id)"
        )

    required_columns = {
    'PBP_allele'
    }

    missing_columns = required_columns.difference(df.columns)

    if missing_columns:
        sys.exit(
            "[combine] PBP table missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    # Extract gene name from entries such as:
    # 3||GBS_1A
    # 1||GBS_2B
    # 4||GBS_2X

    df["PBP_gene"] = (
        df["PBP_allele"]
        .str.extract(r"GBS_(1A|2B|2X)")
        .iloc[:, 0]
        .map({
            "1A": "pbp1A",
            "2B": "pbp2B",
            "2X": "pbp2X"
        })
    )

    # Extract allele number (everything before ||)

    df["PBP_call"] = (
        df["PBP_allele"]
        .str.split("||", n=1, regex=False)
        .str[0]
    )

    cleaned = df.dropna(
        subset=[
            "Sample_ID",
            "PBP_gene",
            "PBP_call"
        ]
    ).copy()

    duplicate_calls = (
        cleaned
        .groupby(
            ["Sample_ID", "PBP_gene"]
        )["PBP_call"]
        .nunique()
    )

    conflicting = duplicate_calls[
        duplicate_calls > 1
    ]

    if not conflicting.empty:
        conflicts = [
            f"{sample}:{gene}"
            for sample, gene in conflicting.index
        ]

        sys.exit(
            "[combine] Conflicting PBP calls detected for: "
            + ", ".join(conflicts)
        )

    cleaned = cleaned.drop_duplicates(
        subset=[
            "Sample_ID",
            "PBP_gene",
            "PBP_call"
        ]
    )

    pbp_wide = (
        cleaned
        .pivot(
            index="Sample_ID",
            columns="PBP_gene",
            values="PBP_call"
        )
        .reset_index()
    )

    pbp_wide.columns.name = None

    ensure_columns(
        pbp_wide,
        PBP_COLUMNS
    )

    return pbp_wide[output_columns]


def main():
    if len(sys.argv) != 5:
        sys.exit(
            "Usage: generate_overall_report.py "
            "<qc_glob> "
            "<typer_path_or_NONE> "
            "<pbp_path_or_NONE> "
            "<output_csv>"
        )

    qc_glob = sys.argv[1]
    typer_path = sys.argv[2]
    pbp_path = sys.argv[3]
    out_csv = sys.argv[4]

    # 1) Load data
    qc_all = read_qc_stack(
        qc_glob
    )

    typer_df = read_typer(
        typer_path
    )

    pbp_df = read_pbp(
        pbp_path
    )

    # 2) Preserve typer columns in their original order
    typer_cols_in_order = [
        c
        for c in list(typer_df.columns)
        if c != 'Sample_ID'
    ]

    # PBP columns are added separately
    typer_cols_in_order = [
        c
        for c in typer_cols_in_order
        if c not in PBP_COLUMNS
    ]

    # 3) Outer-join QC and typer results on Sample_ID
    merged = pd.merge(
        qc_all,
        typer_df,
        on='Sample_ID',
        how='outer'
    )

    # Add PBP results by matching Sample_ID
    # PBP results do not create new rows in the summary
    merged = pd.merge(
        merged,
        pbp_df,
        on='Sample_ID',
        how='left'
    )

    # 4) Ensure fixed QC and PBP columns exist
    ensure_columns(
        merged,
        QC_FIXED
    )

    ensure_columns(
        merged,
        PBP_COLUMNS
    )

    # Insert pbp1A, pbp2B and pbp2X immediately before
    # typer_pipeline_version
    if 'typer_pipeline_version' in typer_cols_in_order:
        version_index = typer_cols_in_order.index(
            'typer_pipeline_version'
        )

        ordered_typer_cols = (
            typer_cols_in_order[:version_index]
            + PBP_COLUMNS
            + typer_cols_in_order[version_index:]
        )
    else:
        ordered_typer_cols = (
            typer_cols_in_order
            + PBP_COLUMNS
        )

    # Keep any remaining columns that were not included above
    extras = [
        c
        for c in merged.columns
        if c not in set(
            ['Sample_ID']
            + QC_FIXED
            + ordered_typer_cols
        )
    ]

    final_cols = (
        ['Sample_ID']
        + QC_FIXED[1:]
        + ordered_typer_cols
        + extras
    )

    # Guarantee unique columns while preserving order
    seen = set()
    ordered = []

    for c in final_cols:
        if c in merged.columns and c not in seen:
            seen.add(c)
            ordered.append(c)

    out = (
        merged
        .reindex(columns=ordered)
        .sort_values('Sample_ID')
    )

    cols_to_drop = [
        c
        for c in (
            'Top_Non-Agalactiae_Species',
        )
        if c in out.columns
    ]

    if cols_to_drop:
        out.drop(
            columns=cols_to_drop,
            inplace=True
        )

    # 5) For PASS rows, missing in-silico fields become
    # MODULE FAILURE
    if 'Overall_QC' in out.columns:
        mask = (
            out['Overall_QC'] == 'PASS'
        )

        # These columns are handled separately or allowed to be empty
        skip_cols = [
            '23S1_SNP',
            '23S3_SNP',
            'gyrA_SNP',
            'parC_SNP',
            'typer_pipeline_version',
            'cps_type',
            'pbp1A',
            'pbp2B',
            'pbp2X'
        ]

        fill_cols = [
            c
            for c in out.columns
            if c not in skip_cols
        ]

        out.loc[
            mask,
            fill_cols
        ] = (
            out.loc[
                mask,
                fill_cols
            ]
            .fillna('MODULE FAILURE')
        )

        # When PBP typing ran, every PASS sample should have
        # one call for each PBP gene.
        #
        # Valid values should already be:
        # allele number, "new", or "neg".
        #
        # A missing PBP row therefore indicates module failure.
        if pbp_path != 'NONE':
            out.loc[
                mask,
                PBP_COLUMNS
            ] = (
                out.loc[
                    mask,
                    PBP_COLUMNS
                ]
                .fillna('MODULE FAILURE')
            )

    # All remaining empty values become NA
    out = out.fillna('NA')

    # 6) Write final summary
    out.to_csv(
        out_csv,
        index=False
    )

    print(
        f"[combine] Wrote {out_csv} with "
        f"{len(out)} rows × "
        f"{len(out.columns)} cols"
    )


if __name__ == "__main__":
    main()