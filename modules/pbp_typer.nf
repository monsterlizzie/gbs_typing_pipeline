process get_pbp_genes {
    label 'farm_mid'

    input:
    tuple val(pair_id), file(contigs)
    file(blactam_ref)
    val(frac_align_len_threshold)
    val(frac_identity_len_threshold)

    output:
    tuple val(pair_id), file("${pair_id}_*bed"), file(contigs),
        optional: true,
        emit: pbp_beds

    tuple val(pair_id),
        file("${pair_id}_PBP_target_status.txt"),
        emit: pbp_target_status

    """
    # Build a blast reference database from the assemblies
    makeblastdb \
        -in ${contigs} \
        -dbtype nucl \
        -out ${pair_id}_contig_blast_db

    # Blast the b-lactam database against the assembly
    blastn \
        -db ${pair_id}_contig_blast_db \
        -query ${blactam_ref} \
        -outfmt 6 \
        -word_size 7 \
        -out ${pair_id}_blast_blactam.out

    # Get BED file of PBP fragments
    get_pbp_genes_from_contigs.py \
        --blast_out_file ${pair_id}_blast_blactam.out \
        --query_fasta ${blactam_ref} \
        --frac_align_len_threshold ${frac_align_len_threshold} \
        --frac_identity_threshold ${frac_identity_len_threshold} \
        --output_prefix ${pair_id}_

    # Record which expected PBP targets were detected
    STATUS="${pair_id}_PBP_target_status.txt"
    : > "\$STATUS"

    for PBP in GBS1A-1 GBS2B-1 GBS2X-1
    do
        if [ -s "${pair_id}_\${PBP}.bed" ]; then
            echo "\${PBP}=TARGET_FOUND" >> "\$STATUS"
        else
            echo "\${PBP}=PBP_NOT_DETECTED" >> "\$STATUS"
        fi
    done

    unlink ${blactam_ref}
    """
}


process get_pbp_alleles {
    label 'farm_mid'

    input:
    tuple val(pair_id), file(bed_file), file(contigs)
    val(pbp_type)
    file(gbs_blactam_db)

    output:
    path "${pair_id}_${pbp_type}_PBP_new_allele.faa",
         optional: true,
         emit: new_pbp

    tuple val(pair_id),
          file("${pair_id}_${pbp_type}_PBP_existing_allele.txt"),
          optional: true,
          emit: existing_pbp

    // NEW: one small status file per PBP
    tuple val(pair_id),
          val(pbp_type),
          file("${pair_id}_${pbp_type}_PBP_status.txt"),
          emit: pbp_status

    """
    STATUS="${pair_id}_${pbp_type}_PBP_status.txt"
    BED="${pair_id}_${pbp_type}.bed"
    FRAG="${pair_id}_${pbp_type}_contig_fragments.fasta"
    PROT="${pair_id}_${pbp_type}.faa"
    BLAST="${pair_id}_blast_${pbp_type}.out"

    : > "\$STATUS"

    # ---------------------------------------------------------
    # PBP target missing
    # ---------------------------------------------------------
    if [ ! -s "\$BED" ]; then
        echo "PBP_NOT_DETECTED" > "\$STATUS"
        exit 0
    fi

    # ---------------------------------------------------------
    # Read BED coordinates
    # ---------------------------------------------------------
    chrom=\$(awk -F '\\t' 'NR==1 {print \$1}' "\$BED")
    start=\$(awk -F '\\t' 'NR==1 {print \$2}' "\$BED")
    end=\$(awk -F '\\t' 'NR==1 {print \$3}' "\$BED")

    # Invalid/unreadable BED coordinates
    if ! [[ "\$start" =~ ^-?[0-9]+\$ && "\$end" =~ ^-?[0-9]+\$ ]]; then
        echo "MODULE_FAILURE" > "\$STATUS"
        echo "Invalid BED coordinates for ${pair_id} ${pbp_type}: start=\$start end=\$end" >&2
        exit 0
    fi

    # ---------------------------------------------------------
    # Partial PBP at start of contig
    # ---------------------------------------------------------
    if [ "\$start" -lt 0 ]; then
        missing=\$((0 - start))

        echo "PARTIAL_PBP_GENE" > "\$STATUS"
        echo \
            "${pair_id} ${pbp_type}: predicted PBP extends \$missing bp before start of \$chrom." \
            >&2

        exit 0
    fi

    # ---------------------------------------------------------
    # Get actual contig length
    # ---------------------------------------------------------
    if ! samtools faidx ${contigs}; then
        echo "MODULE_FAILURE" > "\$STATUS"
        echo "samtools faidx failed for ${pair_id} ${pbp_type}" >&2
        exit 0
    fi

    contig_len=\$(awk -v c="\$chrom" '\$1==c {print \$2}' ${contigs}.fai)

    if [ -z "\$contig_len" ]; then
        echo "MODULE_FAILURE" > "\$STATUS"
        echo "Contig \$chrom not found in assembly for ${pair_id} ${pbp_type}" >&2
        exit 0
    fi

    # ---------------------------------------------------------
    # Partial PBP at end of contig
    # ---------------------------------------------------------
    if [ "\$end" -gt "\$contig_len" ]; then
        missing=\$((end - contig_len))

        echo "PARTIAL_PBP_GENE" > "\$STATUS"
        echo \
            "${pair_id} ${pbp_type}: predicted PBP extends \$missing bp beyond end of \$chrom." \
            >&2

        exit 0
    fi

    # ---------------------------------------------------------
    # Original PBP allele workflow
    # ---------------------------------------------------------

    if ! bedtools getfasta \
        -s \
        -fi ${contigs} \
        -bed "\$BED" \
        -fo "\$FRAG"
    then
        echo "MODULE_FAILURE" > "\$STATUS"
        exit 0
    fi

    if [ ! -s "\$FRAG" ]; then
        echo "PARTIAL_PBP_GENE" > "\$STATUS"
        exit 0
    fi

    if ! translate_pbp_genes.py \
        --blactam_fasta "\$FRAG" \
        --output_file "\$PROT"
    then
        echo "MODULE_FAILURE" > "\$STATUS"
        exit 0
    fi

    if [ ! -s "\$PROT" ]; then
        echo "NO_ALLELE_MATCH" > "\$STATUS"
        exit 0
    fi

    if ! makeblastdb \
        -in ${gbs_blactam_db} \
        -dbtype prot \
        -out ${gbs_blactam_db}
    then
        echo "MODULE_FAILURE" > "\$STATUS"
        exit 0
    fi

    if ! blastp \
        -db ${gbs_blactam_db} \
        -query "\$PROT" \
        -outfmt 6 \
        -out "\$BLAST"
    then
        echo "MODULE_FAILURE" > "\$STATUS"
        exit 0
    fi

    # BLAST ran successfully but returned zero hits
    if [ ! -s "\$BLAST" ]; then
        echo "NO_BLAST_HIT" > "\$STATUS"
        exit 0
    fi

    if ! get_pbp_alleles.py \
        --blast_out_file "\$BLAST" \
        --query_fasta "\$PROT" \
        --output_prefix ${pair_id}_${pbp_type}_PBP
    then
        echo "MODULE_FAILURE" > "\$STATUS"
        exit 0
    fi

    # ---------------------------------------------------------
    # Interpret original outputs
    # ---------------------------------------------------------

    if [ -s ${pair_id}_${pbp_type}_PBP_existing_allele.txt ]; then

        echo "EXISTING_ALLELE" > "\$STATUS"

    elif [ -s ${pair_id}_${pbp_type}_PBP_new_allele.faa ]; then

        echo "NEW_ALLELE" > "\$STATUS"

    else

        echo "NO_ALLELE_MATCH" > "\$STATUS"
    fi

    # Cleanup
    unlink "\$BED" 2>/dev/null || true
    unlink ${contigs} 2>/dev/null || true
    unlink ${gbs_blactam_db} 2>/dev/null || true
    """
}