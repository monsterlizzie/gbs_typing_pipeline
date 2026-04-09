nextflow.enable.dsl=2
// import modules for process
include { FILE_VALIDATION; PREPROCESS; READ_QC } from "$projectDir/modules/preprocess"
include { ASSEMBLY_UNICYCLER; ASSEMBLY_SHOVILL; ASSEMBLY_ASSESS; ASSEMBLY_QC; ASSEMBLY_QC_FALLBACK } from "$projectDir/modules/assembly"
include { GET_REF_GENOME_BWA_DB; MAPPING; SAM_TO_SORTED_BAM; SNP_CALL; HET_SNP_COUNT; MAPPING_QC; MAPPING_QC_FALLBACK} from "$projectDir/modules/mapping"
include { GET_KRAKEN2_DB; TAXONOMY; BRACKEN; TAXONOMY_QC; TAXONOMY_QC_FALLBACK } from "$projectDir/modules/taxonomy"
include { OVERALL_QC } from "$projectDir/modules/overall_qc"
include { GENERATE_SAMPLE_REPORT; GENERATE_OVERALL_REPORT } from "$projectDir/modules/output"
include { serotyping } from './modules/serotyping.nf'
include { srst2_for_res_typing; split_target_RES_seq_from_sam_file; split_target_RES_sequences; freebayes } from './modules/res_alignments.nf'
include { res_typer } from './modules/res_typer.nf'
include { surface_typer } from './modules/surface_typer.nf'
include { getmlst_for_srst2; srst2_for_mlst; get_mlst_allele_and_pileup} from './modules/mlst.nf'
include { get_pbp_genes; get_pbp_alleles } from './modules/pbp_typer.nf'
include { finalise_sero_res_results; finalise_surface_typer_results; finalise_pbp_existing_allele_results; combine_results } from './modules/combine.nf'
include { get_version } from './modules/version.nf'

// this utility process to ensure 'databases/' exists
process INIT_DB_DIR {
    label 'bash_container'
    label 'farm_local'
    publishDir "${params.db}"

    output:
    val "${params.db}", emit: db_dir
    path "do_not_modify", emit: dummy

    script:
    """
    touch do_not_modify
    """
}

// Main pipeline workflow
workflow {

    main:

    INIT_DB_DIR()
    db_dir_ch = INIT_DB_DIR.out.db_dir

    // Get path and prefix of Reference Genome BWA Database, generate from assembly if necessary
    GET_REF_GENOME_BWA_DB(params.ref_genome, db_dir_ch)

    // Get path to Kraken2 Database, download if necessary
    GET_KRAKEN2_DB(params.kraken2_db_remote, db_dir_ch)


    // ----------------------------------------------------------
    // READS: build canonical (deduplicated) read-pairs channel
    // ----------------------------------------------------------
    raw_read_pairs_ch = Channel.fromFilePairs(
        "$params.reads/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}",
        checkIfExists: true
    )

    RAW_READS_ONE_ch = raw_read_pairs_ch
        .map { id, reads -> tuple(id as String, (reads as List).collect { it.toString() }) }
        .groupTuple()
        .map { id, lists ->
            def pick = lists.find { pair -> pair.every { it.endsWith('.gz') } } ?: lists[0]
            tuple(id, pick.collect { file(it) })
        }

    // ----------------------------------------------------------
    // QC SECTION — can be skipped using --skip_qc true
    // ----------------------------------------------------------
    if (!params.skip_qc) {
    // Basic input files validation
    // Output into Channel FILE_VALIDATION.out.result
    FILE_VALIDATION(RAW_READS_ONE_ch)

    // From RAW_READS_ONE_ch, only output valid reads of samples based on FILE_VALIDATION
    VALID_READS_ch = FILE_VALIDATION.out.result
                        .join(RAW_READS_ONE_ch, failOnDuplicate: true)   // (id, status, [R1,R2])
                        .filter { it[1] == 'PASS' }
                        .map { id, status, reads -> tuple(id, reads) }   // keep (id, [R1,R2])

    // Preprocess valid read pairs
    // Output into Channels PREPROCESS.out.processed_reads & PREPROCESS.out.json
    PREPROCESS(VALID_READS_ch)

    // From PREPROCESS.out.json, provide Read QC status
    // Output into Channels READ_QC.out.bases, READ_QC.out.result, READ_QC.out.report
    READ_QC(PREPROCESS.out.json, params.length_low, params.depth)

    // From PREPROCESS.out.processed_reads, only output reads of samples passed Read QC
    READ_QC_PASSED_READS_ch = READ_QC.out.result
                        .join(PREPROCESS.out.processed_reads, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    // Assembler
    switch (params.assembler) {
        case 'shovill':
            ASSEMBLY_ch = ASSEMBLY_SHOVILL(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break
        case 'unicycler':
            ASSEMBLY_ch = ASSEMBLY_UNICYCLER(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break
    }

    ASSEMBLY_ASSESS(ASSEMBLY_ch)

    ASSEMBLY_QC(
        ASSEMBLY_ASSESS.out.report
        .join(READ_QC.out.bases, failOnDuplicate: true),
        params.contigs,
        params.length_low,
        params.length_high,
        params.depth
    )

    // Fallback Assembly QC for samples that failed READ_QC
    READ_QC_FAIL_SAMPLE_ID_ch = READ_QC.out.result.filter { it[1] == "FAIL" }.map { it[0] }
    ASSEMBLY_QC_FALLBACK(READ_QC_FAIL_SAMPLE_ID_ch)

    assembly_qc_all = ASSEMBLY_QC.out.result.mix(ASSEMBLY_QC_FALLBACK.out.result)
    assembly_qc_report_all = ASSEMBLY_QC.out.report.mix(ASSEMBLY_QC_FALLBACK.out.report)

    // Mapping & QC
    MAPPING(GET_REF_GENOME_BWA_DB.out.path, GET_REF_GENOME_BWA_DB.out.prefix, READ_QC_PASSED_READS_ch)
    SAM_TO_SORTED_BAM(MAPPING.out.sam, params.lite)
    SNP_CALL(params.ref_genome, SAM_TO_SORTED_BAM.out.sorted_bam, params.lite)
    HET_SNP_COUNT(SNP_CALL.out.vcf)

    MAPPING_QC(
        SAM_TO_SORTED_BAM.out.ref_coverage
            .join(HET_SNP_COUNT.out.result, failOnDuplicate: true),
        params.ref_coverage,
        params.het_snp_site
    )

    READ_QC_FAIL_ch = READ_QC.out.result.filter { it[1] == "FAIL" }
    MAPPING_QC_FALLBACK(READ_QC_FAIL_ch, params.ref_coverage, params.het_snp_site)

    mapping_qc_all = MAPPING_QC.out.result.mix(MAPPING_QC_FALLBACK.out.result)
    mapping_qc_report_all = MAPPING_QC.out.report.mix(MAPPING_QC_FALLBACK.out.report)

    // Taxonomy & QC
    TAXONOMY(GET_KRAKEN2_DB.out.path, params.kraken2_memory_mapping, READ_QC_PASSED_READS_ch)
    BRACKEN(
        GET_KRAKEN2_DB.out.path,
        TAXONOMY.out.report,
        params.read_len,
        params.classification_level,
        params.threshold
    )
    TAXONOMY_QC(BRACKEN.out.bracken_report, params.sagalactiae_percentage, params.top_non_agalactiae_species_percentage)

    TAXONOMY_QC_FALLBACK(READ_QC_FAIL_SAMPLE_ID_ch)
    taxonomy_qc_all = TAXONOMY_QC.out.result.mix(TAXONOMY_QC_FALLBACK.out.result)
    taxonomy_qc_report_all = TAXONOMY_QC.out.report.mix(TAXONOMY_QC_FALLBACK.out.report)

    // OVERALL_QC (leftmost seed must be unique IDs!)
    OVERALL_QC(
        RAW_READS_ONE_ch.map{ it[0] }     // << use deduped IDs
        .join(FILE_VALIDATION.out.result, failOnDuplicate: true, remainder: true)
        .join(READ_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(assembly_qc_all, failOnDuplicate: true, remainder: true)
        .join(mapping_qc_all, failOnDuplicate: true, remainder: true)
        .join(taxonomy_qc_all, failOnDuplicate: true, remainder: true)
    )

    // Reads that passed overall QC
    OVERALL_QC_PASSED_READS_ch = OVERALL_QC.out.result
                        .join(READ_QC_PASSED_READS_ch, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    // Extract only paired reads (R1, R2)
    OVERALL_QC_PASSED_PAIRED_READS_ch = OVERALL_QC_PASSED_READS_ch.map { id, r1, r2, unpaired -> [id, [r1, r2]] }

    // Assemblies for passed samples
    OVERALL_QC_PASSED_ASSEMBLIES_ch = OVERALL_QC.out.result
                            .join(ASSEMBLY_ch, failOnDuplicate: true)
                            .filter { it[1] == 'PASS' }
                            .map { it[0, 2..-1] }
                            

    // Sample reports
    GENERATE_SAMPLE_REPORT(
        RAW_READS_ONE_ch.map{ it[0] }       // << use deduped IDs
        .join(READ_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(assembly_qc_report_all, failOnDuplicate: true, remainder: true)
        .join(mapping_qc_report_all, failOnDuplicate: true, remainder: true)
        .join(taxonomy_qc_report_all, failOnDuplicate: true, remainder: true)
        .join(OVERALL_QC.out.report, failOnDuplicate: true, failOnMismatch: true)
        .map { row ->
            def sample_id = row[0]
            def report_paths = row[1..-1].findAll { it != null && it != "NA" }
            [sample_id, report_paths]
        }
    )
} // <--- end skip_qc block

// ----------------------------------------------------------
// Fallback if QC is skipped — use reads directly
// ----------------------------------------------------------
if (params.skip_qc) {
    log.info "Skipping QC - using raw reads directly for typer modules."
    OVERALL_QC_PASSED_PAIRED_READS_ch = Channel.fromFilePairs(
        "$params.reads/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}",
        checkIfExists: true
    ).map { id, reads -> tuple(id, reads) }
}

 // Guard clause — allow QC-only runs but prevent 'nothing to do'
if (!params.run_sero_res && !params.run_surfacetyper && !params.run_mlst && !params.run_pbptyper) {
    if (params.skip_qc) {
        println(" Error: all typer modules disabled and QC skipped — nothing to run.")
        System.exit(1)
    } else {
        log.info " Running QC-only mode (no typing modules enabled)."
    }
    Channel.fromFilePairs( params.reads, checkIfExists: true )
            .set { read_pairs_ch }
}

   

    // Check outputs dir
    if (params.output == ""){
        println("Please specify the results directory with --params.output.")
        println("Print help with nextflow main.nf --help")
        System.exit(1)
    }

    // Param range checks 
    if (params.gbs_res_min_coverage < 0 | params.gbs_res_min_coverage > 100){
        println("--gbs_res_min_coverage value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
    if (params.gbs_res_max_divergence < 0 | params.gbs_res_max_divergence > 100){
        println("--gbs_res_max_divergence value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
    other_res_min_coverage_list = params.other_res_min_coverage.toString().tokenize(' ')
    for (other_res_min_coverage in other_res_min_coverage_list){
        if (other_res_min_coverage.toDouble() < 0 | other_res_min_coverage.toDouble() > 100){
            println("--other_res_min_coverage value(s) not in range. Please specify a value between 0 and 100.")
            System.exit(1)
        }
    }
    other_res_max_divergence_list = params.other_res_max_divergence.toString().tokenize(' ')
    for (other_res_max_divergence in other_res_max_divergence_list){
        if (other_res_max_divergence.toDouble() < 0 | other_res_max_divergence.toDouble() > 100){
            println("--other_res_max_divergence value(s) not in range. Please specify a value between 0 and 100.")
            System.exit(1)
        }
    }
    if (params.restyper_min_read_depth < 0){
        println("--restyper_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }
    if (params.serotyper_min_read_depth < 0){
        println("--serotyper_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }
    if (params.mlst_min_coverage < 0 | params.mlst_min_coverage > 100){
        println("--mlst_min_coverage value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
    if (params.mlst_min_read_depth < 0){
        println("--mlst_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }
    if (params.surfacetyper_min_coverage < 0 | params.surfacetyper_min_coverage > 100){
        println("--surfacetyper_min_coverage value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
    if (params.surfacetyper_max_divergence < 0 | params.surfacetyper_max_divergence > 100){
        println("--surfacetyper_max_divergence value not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
    if (params.surfacetyper_min_read_depth < 0){
        println("--surfacetyper_min_read_depth value not in range. Please specify a value of 0 or above.")
        System.exit(1)
    }

    // -------- Workflows --------

    workflow GBS_RES {
        take: reads
        main:
            gbs_res_typer_db = file(params.gbs_res_typer_db, checkIfExists: true)
            gbs_res_targets_db = file(params.gbs_res_targets_db, checkIfExists: true)
            split_target_RES_sequences(gbs_res_typer_db, gbs_res_targets_db)
            srst2_for_res_typing(reads, gbs_res_typer_db, params.gbs_res_min_coverage, params.gbs_res_max_divergence)
            fullgenes = srst2_for_res_typing.out.fullgenes
            split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, gbs_res_targets_db)
            freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out)
            consensus = freebayes.out.consensus
        emit:
            fullgenes
            consensus
    }

    workflow OTHER_RES {
        take: reads
        main:
            other_res_db = file(params.other_res_db, checkIfExists: true)
            srst2_for_res_typing(reads, other_res_db, params.other_res_min_coverage, params.other_res_max_divergence)
            fullgenes = srst2_for_res_typing.out.fullgenes
        emit:
            fullgenes
    }

    workflow MLST {
        take: reads
        main:
            getmlst_for_srst2()
            srst2_for_mlst(getmlst_for_srst2.out.getmlst_results, reads, params.mlst_min_coverage)
            get_mlst_allele_and_pileup(srst2_for_mlst.out.bam_and_srst2_results, params.mlst_min_read_depth)
            new_alleles = get_mlst_allele_and_pileup.out.new_alleles
            pileup = get_mlst_allele_and_pileup.out.pileup
            existing_alleles = get_mlst_allele_and_pileup.out.existing_alleles
            status = get_mlst_allele_and_pileup.out.new_alleles_status
            srst2_results = srst2_for_mlst.out.srst2_results
        emit:
            new_alleles
            pileup
            existing_alleles
            status
            srst2_results
    }

    workflow PBP1A {
        take: pbp_typer_output
        main:
            get_pbp_alleles(pbp_typer_output, 'GBS1A-1', file(params.gbs_blactam_1A_db, checkIfExists: true))
            get_pbp_alleles.out.new_pbp.subscribe { it -> it.copyTo(file("${params.output}")) }
            finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))
        emit:
            finalise_pbp_existing_allele_results.out
    }

    workflow PBP2B {
        take: pbp_typer_output
        main:
            get_pbp_alleles(pbp_typer_output, 'GBS2B-1', file(params.gbs_blactam_2B_db, checkIfExists: true))
            get_pbp_alleles.out.new_pbp.subscribe { it -> it.copyTo(file("${params.output}")) }
            finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))
        emit:
            finalise_pbp_existing_allele_results.out
    }

    workflow PBP2X {
        take: pbp_typer_output
        main:
            get_pbp_alleles(pbp_typer_output, 'GBS2X-1', file(params.gbs_blactam_2X_db, checkIfExists: true))
            get_pbp_alleles.out.new_pbp.subscribe { it -> it.copyTo(file("${params.output}")) }
            finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))
        emit:
            finalise_pbp_existing_allele_results.out
    }

    if (params.run_sero_res){
        serotyping(OVERALL_QC_PASSED_PAIRED_READS_ch, file(params.sero_gene_db, checkIfExists: true), params.serotyper_min_read_depth)

        GBS_RES(OVERALL_QC_PASSED_PAIRED_READS_ch)
        OTHER_RES(OVERALL_QC_PASSED_PAIRED_READS_ch)

        GBS_RES.out.fullgenes
            .join(GBS_RES.out.consensus)
            .join(OTHER_RES.out.fullgenes)
            .set { res_files_ch }

        res_typer(res_files_ch, params.restyper_min_read_depth, file(params.config, checkIfExists: true))

        sero_res_ch = serotyping.out.join(res_typer.out.res_out)
        finalise_sero_res_results(sero_res_ch, file(params.config, checkIfExists: true))

        finalise_sero_res_results.out.sero_res_incidence
            .collectFile(name: file("${params.output}/${params.sero_res_incidence_out}"), keepHeader: true)
        finalise_sero_res_results.out.res_alleles_variants
            .collectFile(name: file("${params.output}/${params.alleles_variants_out}"), keepHeader: true)
        finalise_sero_res_results.out.res_variants
            .collectFile(name: file("${params.output}/${params.variants_out}"), keepHeader: true)

        res_typer.out.res_accessions
            .collectFile(name: file("${params.output}/${params.res_accessions_out}"))
    }

    if (params.run_mlst){
        MLST(OVERALL_QC_PASSED_PAIRED_READS_ch)
        MLST.out.new_alleles.subscribe { it -> it.copyTo(file("${params.output}")) }
        MLST.out.pileup.subscribe { it -> it.copyTo(file("${params.output}")) }
        MLST.out.existing_alleles
            .collectFile(name: file("${params.output}/${params.existing_mlst_alleles_out}"), keepHeader: true, sort: true)
        MLST.out.status
            .collectFile(name: file("${params.output}/${params.new_mlst_alleles_status}"), keepHeader: false, sort: true)
    }

    if (params.run_surfacetyper){
        surface_typer(OVERALL_QC_PASSED_PAIRED_READS_ch, file(params.gbs_surface_typer_db, checkIfExists: true),
            params.surfacetyper_min_read_depth, params.surfacetyper_min_coverage, params.surfacetyper_max_divergence)

        finalise_surface_typer_results(surface_typer.out, file(params.config, checkIfExists: true))

        finalise_surface_typer_results.out.surface_protein_incidence
            .collectFile(name: file("${params.output}/${params.surface_protein_incidence_out}"), keepHeader: true)
        finalise_surface_typer_results.out.surface_protein_variants
            .collectFile(name: file("${params.output}/${params.surface_protein_variants_out}"), keepHeader: true)
    }

    if (params.run_pbptyper) {
        if (!params.pbp_contig) {
        println("Please specify contigs with --pbp_contig.")
        println("Print help with --pbp_contig")
        System.exit(1)
    }

    contig_paths = Channel
        .fromPath(params.pbp_contig, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

    get_pbp_genes(
        contig_paths,
        file(params.gbs_blactam_db, checkIfExists: true),
        params.pbp_frac_align_threshold,
        params.pbp_frac_identity_threshold
    )

    PBP1A(get_pbp_genes.out)
    PBP2B(get_pbp_genes.out)
    PBP2X(get_pbp_genes.out)

    PBP1A.out
        .concat(PBP2B.out, PBP2X.out)
        .set { PBP_all }

    PBP_all
        .collectFile(name: file("${params.output}/${params.existing_pbp_alleles_out}"), keepHeader: true, sort: true)
    }

    if (params.run_sero_res & params.run_surfacetyper & params.run_mlst){
        get_version()
        version_ch = get_version.out

        combined_ch = serotyping.out
            .join(res_typer.out.res_out)
            .join(surface_typer.out)
            .join(MLST.out.srst2_results)

        combine_results(combined_ch, file(params.config, checkIfExists: true), version_ch)

        combine_results.out
            .collectFile(name: file("${params.output}/${params.gbs_typer_report}"), keepHeader: true, sort: true)
    }

    // Typer CSV channel from the collectFile 
    typer_csv_ch = Channel.empty()
    if (params.run_sero_res & params.run_surfacetyper & params.run_mlst) {
        typer_csv_ch = combine_results.out.collectFile(
            name: file("${params.output}/${params.gbs_typer_report}"),
            keepHeader: true,
            sort: true
        )
    }

    // Barrier for overall report
    if (!params.skip_qc) {
        done_ch    = GENERATE_SAMPLE_REPORT.out.collect()
        qc_glob_ch = done_ch.map { "${params.output}/qc_reports/*_qc.csv" }
    } else {
        // If QC was skipped, create a dummy channel so downstream logic still works
        qc_glob_ch = Channel.value('NONE')
    }
   

    // Fallback for typer channel (QC-only run)
    if (!params.run_sero_res && !params.run_surfacetyper && !params.run_mlst && !params.run_pbptyper) {
    typer_csv_ch = Channel.value('NONE')
    }

    // Typer path (or NONE if typer wasn’t run)
    typer_path_ch = typer_csv_ch.ifEmpty { Channel.value('NONE') }.map { it.toString() }

    // Fire the overall report (still runs, even if typer_path_ch == 'NONE')
    GENERATE_OVERALL_REPORT(qc_glob_ch, typer_path_ch)



}
