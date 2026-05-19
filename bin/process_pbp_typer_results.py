#!/usr/bin/env python3

import argparse
import glob
import csv
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq


def parse_args():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--srst2_prefix", required=True)
    parser.add_argument("--pbp_ref", required=True)
    parser.add_argument("--output_prefix", required=True)

    return parser.parse_args()


def clean_sample_id(sample):
    """
    Remove SRST2 preprocessing prefix.
    """

    return sample.replace("processed-", "")


def allele_to_number(allele):
    """
    Convert allele name to compact allele number.

    Example:
    GBS1A-1 -> 1
    GBS2B-17 -> 17
    """

    if "-" in allele:
        return allele.split("-")[-1].replace("*", "")

    return allele.replace("*", "")


def classify_from_diffs(allele, diffs):
    """
    Determine final allele classification.

    Exact match:
        1

    Novel SNPs:
        new
    """

    diffs = str(diffs).strip()

    # Exact allele match
    if diffs in ["", "nan", "0", "0snp", "None"]:
        return allele_to_number(allele)

    # Novel SNP-containing allele
    return "new"


def parse_results(results_file):
    """
    Parse SRST2 fullgenes results table.
    """

    calls = {
        "pbp1A": "neg",
        "pbp2B": "neg",
        "pbp2X": "neg"
    }

    novel_genes = []
    sample_id = None

    with open(results_file) as f:

        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}

        for line in f:

            row = line.rstrip("\n").split("\t")

            sample = clean_sample_id(row[idx["Sample"]])
            gene = row[idx["gene"]]
            allele = row[idx["allele"]]
            diffs = row[idx["diffs"]]

            sample_id = sample

            call = classify_from_diffs(allele, diffs)

            # PBP1A
            if "GBS1A" in gene:

                calls["pbp1A"] = call

                if "new" in call:
                    novel_genes.append(gene)

            # PBP2B
            elif "GBS2B" in gene:

                calls["pbp2B"] = call

                if "new" in call:
                    novel_genes.append(gene)

            # PBP2X
            elif "GBS2X" in gene:

                calls["pbp2X"] = call

                if "new" in call:
                    novel_genes.append(gene)

    return sample_id, calls, novel_genes


def parse_mpileup_bases(bases, ref_base):
    """
    Parse mpileup base string into usable base calls.
    """

    calls = []
    i = 0

    while i < len(bases):

        c = bases[i]

        # Read-start marker
        if c == "^":
            i += 2
            continue

        # Read-end marker
        if c == "$":
            i += 1
            continue

        # Indels
        if c in "+-":

            i += 1
            num = ""

            while i < len(bases) and bases[i].isdigit():
                num += bases[i]
                i += 1

            if num:
                i += int(num)

            continue

        # Reference base
        if c in ".,":

            calls.append(ref_base.upper())

        # SNP base
        elif c.upper() in ["A", "C", "G", "T"]:

            calls.append(c.upper())

        i += 1

    return calls


def load_reference(ref_fasta):
    """
    Load reference nucleotide sequences.
    """

    return {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(ref_fasta, "fasta")
    }


def build_consensus_from_pileup(pileup_file, ref_sequences):
    """
    Reconstruct consensus nucleotide sequence from pileup.
    """

    consensus = {
        gene: list(seq)
        for gene, seq in ref_sequences.items()
    }

    with open(pileup_file) as f:

        for line in f:

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 5:
                continue

            gene = parts[0]
            pos = int(parts[1]) - 1
            ref_base = parts[2].upper()
            depth = int(parts[3])
            bases = parts[4]

            if gene not in consensus:
                continue

            if depth == 0:
                continue

            base_calls = parse_mpileup_bases(bases, ref_base)

            if not base_calls:
                continue

            majority_base = Counter(base_calls).most_common(1)[0][0]

            if pos < len(consensus[gene]):
                consensus[gene][pos] = majority_base

    return {
        gene: "".join(seq)
        for gene, seq in consensus.items()
    }


def translate_nt(seq):
    """
    Translate nucleotide sequence to amino acid sequence.
    """

    seq = seq.upper().replace("-", "")

    return str(Seq(seq).translate(to_stop=True))


def write_summary(output_prefix, sample_id, calls):
    """
    Write compact allele summary table.
    """

    out_file = f"{output_prefix}_pbp_alleles.txt"

    with open(out_file, "w") as out:

        out.write("Sample_ID\tpbp1A\tpbp2B\tpbp2X\n")

        out.write(
            f"{sample_id}\t"
            f"{calls['pbp1A']}\t"
            f"{calls['pbp2B']}\t"
            f"{calls['pbp2X']}\n"
        )


def write_novel_outputs(output_prefix, sample_id, novel_genes, consensus_nt):
    """
    Write:
    1. FASTA of novel amino acid sequences
    2. TSV summary table
    """

    if not novel_genes:
        return

    faa_file = f"{output_prefix}_PBP_new_allele.faa"
    report_file = f"{output_prefix}_PBP_new_alleles_report.tsv"

    with open(faa_file, "w") as faa, open(report_file, "w", newline="") as report:

        writer = csv.writer(report, delimiter="\t")

        writer.writerow([
            "source_file",
            "Contig",
            "sequence"
        ])

        for gene in novel_genes:

            if gene not in consensus_nt:
                continue

            aa_seq = translate_nt(consensus_nt[gene])

            source_file = f"{sample_id}_{gene}_PBP_new_allele.faa"

            contig = gene

            faa.write(f">{sample_id}_{gene}_PBP_new_allele\n")
            faa.write(f"{aa_seq}\n")

            writer.writerow([
                source_file,
                contig,
                aa_seq
            ])


def main():

    args = parse_args()

    # Locate SRST2 outputs
    results_matches = glob.glob(
        f"{args.srst2_prefix}*fullgenes*results.txt"
    )

    pileup_matches = glob.glob(
        f"{args.srst2_prefix}*.pileup"
    )

    # No results found
    if not results_matches:

        sample_id = args.output_prefix

        calls = {
            "pbp1A": "neg",
            "pbp2B": "neg",
            "pbp2X": "neg"
        }

        write_summary(args.output_prefix, sample_id, calls)

        return

    results_file = results_matches[0]

    sample_id, calls, novel_genes = parse_results(results_file)

    # Write summary table
    write_summary(args.output_prefix, sample_id, calls)

    # Generate novel allele outputs
    if novel_genes and pileup_matches:

        ref_sequences = load_reference(args.pbp_ref)

        consensus_nt = build_consensus_from_pileup(
            pileup_matches[0],
            ref_sequences
        )

        write_novel_outputs(
            args.output_prefix,
            sample_id,
            novel_genes,
            consensus_nt
        )


if __name__ == "__main__":
    main()