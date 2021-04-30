process bolt_lmm {
    label "bolt"
    publishDir "results/bolt/", mode: 'copy'
    tag "all"
    input:
        path(hardcalls)
        path(imputed)
        path(pheno)
    output:
        path("bolt_all_chr.*"), emit: bolt_results
    script:
        """
        for bgen in *bgen
            do
                prefix=\$(basename \$bgen .bgen)
                echo \$bgen \$prefix.sample >> list_bgen.txt
            done
        bolt \
            --bed chr{1:22}.bed \
            --bim chr{1:22}.bim \
            --fam chr1.fam \
            --phenoFile "${pheno}" \
            --phenoCol standing_height \
            --covarFile "${pheno}" \
            --covarCol genetic_sex \
            --qCovarCol PC{1:15} \
            --LDscoresFile \$BOLT_LMM_ROOT/tables/LDSCORE.1000G_EUR.tab.gz \
            --geneticMapFile \$BOLT_LMM_ROOT/tables/genetic_map_hg19_withX.txt.gz \
            --lmmForceNonInf \
            --numThreads ${task.cpus} \
            --statsFile bolt_all_chr.hardcalls.txt \
            --bgenSampleFileList list_bgen.txt \
            --statsFileBgenSnps bolt_all_chr.bgen.txt \
            --verboseStats
        """
}