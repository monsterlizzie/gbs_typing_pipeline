#!/usr/bin/env python3
from collections import defaultdict
import argparse, sys
from lib.seq_data import SeqData, BlastData

IDENTITY_THRESHOLD = 50
FRAGMENT_LENGTH_THRESHOLD = 0.5


def write_content(content, output_file):
    with open(output_file, 'w') as out:
        for item in content:
            out.write(item)

def get_pbp_call_from_output_prefix(output_prefix):
    if 'GBS1A-1' in output_prefix:
        return 'NEW||GBS_1A'

    if 'GBS2B-1' in output_prefix:
        return 'NEW||GBS_2B'

    if 'GBS2X-1' in output_prefix:
        return 'NEW||GBS_2X'

    raise ValueError(
        f'Cannot determine PBP type from output prefix: '
        f'{output_prefix}'
    )

def get_imperfect_allele(best_blast_hits, query_seq_data):
    seq_lengths = query_seq_data.calculate_seq_length()
    return ['>' + hit + '\n' + query_seq_data.get_data()[hit] + '\n' for hit in best_blast_hits if (100 > float(best_blast_hits[hit][1]) >= IDENTITY_THRESHOLD) and (seq_lengths[hit]/int(best_blast_hits[hit][2]) >= FRAGMENT_LENGTH_THRESHOLD)]


def get_identical_allele(best_blast_hits):
    return ['Contig\tPBP_allele\n' + hit + '\t' + best_blast_hits[hit][0] + '\n' for hit in best_blast_hits if float(best_blast_hits[hit][1]) == 100]


def get_arguments():
    parser = argparse.ArgumentParser(description='Output the start and end positions of b-lactam genes from contigs in a BED file.')
    parser.add_argument('--blast_out_file', '-b', dest='blast_out', required=True,
                        help='Input BLAST results file.', type = str)
    parser.add_argument('--query_fasta', '-f', dest='fasta_qu', required=True,
                        help='Fasta file query.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                        help='Output prefix.', type = str)
    return parser


def main():
    args = get_arguments().parse_args()

    # Load blast data
    blast_data = BlastData(args.blast_out)

    # Get blast best hit
    best_blast_hits = blast_data.get_best_hit()

    # Load sequence data
    seq_data = SeqData(args.fasta_qu)

    # Get identical allele
    identical_allele = get_identical_allele(best_blast_hits)

    # Get imperfect allele
    imperfect_allele = get_imperfect_allele(best_blast_hits, seq_data)

    # Write data
    
    if identical_allele:
        write_content(
            identical_allele,
            args.output + '_existing_allele.txt'
        )

    elif imperfect_allele:
        # Preserve the novel amino-acid sequence
        write_content(
            imperfect_allele,
            args.output + '_new_allele.faa'
        )

        # Also create a machine-readable PBP result so the
        # downstream report records this as NEW.
        pbp_call = get_pbp_call_from_output_prefix(
            args.output
        )

        contig = next(iter(best_blast_hits))

        new_result = [
            'Contig\tPBP_allele\n'
            + contig
            + '\t'
            + pbp_call
            + '\n'
        ]

        write_content(
            new_result,
            args.output + '_existing_allele.txt'
        )

    else:
        print(
            'Warning: No qualifying PBP allele hit found.',
            file=sys.stderr
        )


if __name__ == "__main__":
    sys.exit(main())
