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

  output:
  path "summary.csv"

  script:
  """
  set -e
  echo "QC glob: ${qc_glob}"
  echo "Typer: ${typer_path}"
  python3 ${projectDir}/bin/generate_overall_report.py \
      "${qc_glob}" \
      "${typer_path}" \
      summary.csv
  """
}
