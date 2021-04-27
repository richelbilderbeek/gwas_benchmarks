process regenie_step_1 {
    label "regenie"
    input:
        tuple path(hardcalls_bed), path(hardcalls_bim), path(hardcalls_fam), path(pheno)
    output:
        path("regenie_step_1_pred.list"), emit: regenie_predictions
    script:
        """
        regenie \
            --step 1 \
            --bed merged_hardcalls \
            --phenoFile ${pheno} \
            --covarFile ${pheno} \
            --phenoCol standing_height \
            --covarCol genetic_sex \
            --bsize 1000 \
            --lowmem \
            --lowmem-prefix \$SNIC_TMP \
            --threads ${task.cpus} \
            --out regenie_step_1
        """
}


process regenie_step_2 {
    label "regenie"
    publishDir "results/regenie/", mode: 'copy'
    input:
        tuple val(prefix), path(genotypes), path(pheno), path(pred)
    output:
        path("*.regenie"), emit: results
    script:
        """
        regenie \
            --step 2 \
            --bgen ${bgen} \
            --ref-first \
            --sample ${sample} \
            --phenoFile ${pheno} \
            --covarFile ${pheno} \
            --phenoCol standing_height \
            --covarCol genetic_sex \
            --firth 0.01 --approx \
            --pred ${step_1_pred} \
            --bsize 400 \
            --split \
            --threads ${task.cpus} \
            --out ${prefix}
        """
}