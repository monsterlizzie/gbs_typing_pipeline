#!/usr/bin/env python3
"""Combine staged QC reports without inventing unperformed typing results."""
import sys
from pathlib import Path
from generate_overall_report import read_qc_stack

QC_COLUMNS = ['Sample_ID', 'Read_QC', 'Assembly_QC', 'Mapping_QC',
              'Taxonomy_QC', 'Overall_QC', 'Bases', 'Contigs#',
              'Assembly_Length', 'Seq_Depth', 'Ref_Cov_%', 'Het-SNP#',
              'S.agalactiae_%', 'Top_Non-agalactiae_Species',
              'Top_Non-agalactiae_Species_%']

def write_report(qc_glob, version_file, output):
    frame = read_qc_stack(qc_glob)
    if frame['Sample_ID'].isna().any() or frame['Sample_ID'].duplicated().any():
        raise ValueError('QC sample identifiers must be present and unique')
    missing = set(QC_COLUMNS[:6]) - set(frame.columns)
    if missing:
        raise ValueError(f'Missing QC status columns: {sorted(missing)}')
    version = Path(version_file).read_text().splitlines()
    if len(version) != 2 or version[0] != 'version' or not version[1].strip():
        raise ValueError('Expected version header and one nonempty version value')
    frame = frame.reindex(columns=QC_COLUMNS).sort_values('Sample_ID')
    frame['pipeline_version'] = version[1].strip()
    frame.fillna('NA').to_csv(output, index=False)

if __name__ == '__main__':
    write_report(*sys.argv[1:])
