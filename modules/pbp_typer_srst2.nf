process pbp_typer_srst2 {

    label 'srst2'
    label 'farm_mid'

    publishDir "${params.output}/typer", mode: 'copy'

    input:
    tuple val(pair_id), path(read1), path(read2), path(unpaired)

    path(gbs_blactam_ref)
    path(gbs_blactam_1A_db)
    path(gbs_blactam_2B_db)
    path(gbs_blactam_2X_db)

    val(min_coverage)
    val(max_divergence)

    output:
    tuple val(pair_id), file("${pair_id}_pbp_alleles.txt"), emit: pbp_out
    tuple val(pair_id), file("${pair_id}_PBP_new_allele.faa"), optional: true, emit: new_pbp
    tuple val(pair_id), file("${pair_id}_PBP_new_alleles_report.tsv"), optional: true, emit: novel_report

    script:
    """
    set +e

    srst2 \\
        --samtools_args '\\-A' \\
        --input_pe ${read1} ${read2} \\
        --output ${pair_id}_PBP \\
        --log \\
        --save_scores \\
        --min_coverage ${min_coverage} \\
        --max_divergence ${max_divergence} \\
        --gene_db ${gbs_blactam_ref}

    python3 ${projectDir}/bin/process_pbp_typer_results.py \\
        --srst2_prefix ${pair_id}_PBP \\
        --pbp_ref ${gbs_blactam_ref} \\
        --output_prefix ${pair_id}

    touch ${pair_id}_pbp_alleles.txt
    """
}