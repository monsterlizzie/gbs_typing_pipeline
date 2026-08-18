process GENERATE_SAMPLE_REPORT {
    label 'python_container'
    label 'farm_low'

    tag "$sample_id"

    publishDir "${params.output}/qc_reports", mode: 'copy'

    input:
    tuple val(sample_id), path("${sample_id}_process_report_?.csv")

    output:
    path "${sample_id}_qc.csv", emit: report

    script:
    """
    generate_sample_report.py "$sample_id" "${sample_id}_qc.csv" ${sample_id}_process_report_*.csv
    """
}


process GENERATE_OVERALL_REPORT {
    label 'python_container'
    label 'farm_low'
  
    publishDir "${params.output}", mode: 'copy'

    input:
    val qc_glob
    val typer_path
    val pbp_path

    // NEW
    val pbp_target_status_path
    val pbp_allele_status_path

    output:
    path "summary.csv"

    script:
    """
    set -e

    echo "QC glob: ${qc_glob}"
    echo "Typer: ${typer_path}"
    echo "PBP alleles: ${pbp_path}"
    echo "PBP target status: ${pbp_target_status_path}"
    echo "PBP allele status: ${pbp_allele_status_path}"

    python3 ${projectDir}/bin/generate_overall_report.py \
        "${qc_glob}" \
        "${typer_path}" \
        "${pbp_path}" \
        "${pbp_target_status_path}" \
        "${pbp_allele_status_path}" \
        summary.csv
    """
}